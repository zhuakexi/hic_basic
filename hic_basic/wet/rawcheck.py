import concurrent.futures
import hashlib
import re
import shutil
from pathlib import Path

def calculate_md5_for_file(file_path):
    md5 = hashlib.md5()
    with open(file_path, "rb") as f:
        for chunk in iter(lambda: f.read(4096), b""):
            md5.update(chunk)
    return file_path, md5.hexdigest()
def ls_md5_files(download_dir):
    """
    Find all md5.txt files in download_dir.
    Input:
        download_dir: directory containing downloaded files
    Output:
        md5_lines: list of tuples, (parent_dir, md5, filename)
    """
    md5_files = Path(download_dir).glob('**/md5.txt')
    md5_lines = [(md5_file.parent, line.split()) for md5_file in md5_files for line in open(md5_file).readlines()]
    return md5_lines
def verify_md5(download_dir, result_file):
    """
    Verify md5 checksums of downloaded files.
    Input:
        download_dir: directory containing downloaded files
        result_file: file to write md5 checksum results
    Output:
        True if all md5 checksums are OK, False otherwise
    """
    # --- create result file --- #
    result_file = Path(result_file)
    result_file.parent.mkdir(parents=True, exist_ok=True)
    if result_file.exists():
        result_file.unlink()
    result_file.touch()

    # --- find all md5.txt files in download_dir --- #
    md5_lines = ls_md5_files(download_dir)
    # --- verify md5 checksums --- #
    nfails = 0
    with concurrent.futures.ProcessPoolExecutor(8) as executor:
        future_to_md5_line = {executor.submit(calculate_md5_for_file, parent_dir / filename): (parent_dir, filename, expected_md5) 
                              for parent_dir, (expected_md5, filename) in md5_lines}
        for future in concurrent.futures.as_completed(future_to_md5_line):
            parent_dir, filename, expected_md5 = future_to_md5_line[future]
            try:
                file, md5 = future.result()
                md5_results = filename + (' OK' if md5 == expected_md5 else ' FAILED')
                with open(result_file, 'a') as f:
                    f.write(md5_results + '\n')
                if md5 != expected_md5:
                    print(f'ERROR: md5 checksum of {file} does not match.')
                    nfails += 1
            except Exception as exc:
                print('%r generated an exception: %s' % (file, exc))
    if nfails == 0:
        print('All md5 checksums are OK.')
        return True
    else:
        print(f'{nfails} md5 checksums are not OK.')
        return False
def ls_samples(download_dir, show=True):
    """
    Find all samples in download_dir. Samples are defined as the directories containing fastq files.
    Input:
        download_dir: directory containing downloaded files
        show: print samples if True
    Output:
        samples: list of sample names
    """
    # --- find all fastq files in download_dir --- #
    # fastq_files ends with .fastq.gz or fq.gz
    fastq_files = [path for path in Path(download_dir).rglob('*')
                if re.search(r'\.(fastq|fq)\.gz$', str(path), re.IGNORECASE)]
    # --- find all samples in download_dir --- #
    samples = {}
    for fastq_file in fastq_files:
        sample_name = fastq_file.parent.name
        sample = fastq_file.parent # each sample directory may contain multiple fastq files
        if sample_name not in samples:
            samples[sample_name] = set()
        samples[sample_name].add(sample)
    # --- sort sample by key name --- #
    if show:
        print(f'Found {len(samples)} samples:')
        for sample_name, sample in sorted(samples.items(), key=lambda x: x[0]):
            print(f'{sample_name}')
    return samples
def split_download_dir(download_dir, task_dirps, sub="Rawdata"):
    """
    Split download_dir into multiple directories, each containing multiple samples and 1 md5 file for each sample.
    Input:
        download_dir: directory containing downloaded files
        task_dirps: dictionary, {task_dirp: [sample_names]}
        sub: only consider samples that have this substring in their path
    Output:
        None
    """
    # --- find all samples and md5 lines in download_dir --- #
    real_samples = ls_samples(download_dir)
    md5_lines = ls_md5_files(download_dir)
    # --- split samples into multiple task directories --- #
    for task_dirp, sample_names in task_dirps.items():
        task_dirp = Path(task_dirp)
        task_dirp.mkdir(parents=True, exist_ok=True)
        new_md5_file = task_dirp.parent / "md5.txt" # md5 file is outside of task directory
        if new_md5_file.exists():
            new_md5_file.unlink()
        new_md5_file.touch()
        for sample_name in sample_names:
            destination_dirp = task_dirp / sample_name
            source_dirps = [] # all possible source dirp, but only one is valid in the end
            # --- find all fastq files of this sample --- #
            if sub:
                source_dirps = [source_dirp for source_dirp in real_samples[sample_name] if sub in str(source_dirp)]
            else:
                source_dirps = list(real_samples[sample_name])
            if len(source_dirps) == 0:
                print(f'WARNING: no sample directory contains {sub} in {sample_name}')
                continue
            elif len(source_dirps) > 1:
                raise ValueError(f'More than one sample directory for {sample_name}: {source_dirps}')
            else: # only this case is valid
                source_dirp = source_dirps[0]
            # --- add new line to new md5 file --- #
            valid_md5_lines = 0
            for parent_dir, (expected_md5, filename) in md5_lines:
                filename = Path(filename)
                if filename.parent.name == sample_name:
                    if sub and sub not in str(filename):
                        continue
                    relative_filename = Path(filename.parent.parent.name) /  Path(filename.parent.name) / Path(filename.name)
                    with open(new_md5_file, 'a') as f:
                        f.write(f'{expected_md5} {relative_filename}\n')
                    valid_md5_lines += 1
            if valid_md5_lines == 0:
                print(f'WARNING: no md5 line for {sample_name}')
            # --- copy sample directory to task directory --- #
            if destination_dirp.exists():
                shutil.rmtree(destination_dirp)
            shutil.copytree(source_dirp, destination_dirp)