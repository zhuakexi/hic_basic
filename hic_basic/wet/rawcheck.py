import concurrent.futures
import hashlib
import re
import shutil
from pathlib import Path

def calculate_md5_for_file(file_path):
    """
    Calculate md5 checksum for a file.
    This give exactly the same result as md5sum command.
    Input:
        file_path: path to the file
    Output:
        file_path: path to the file
        md5: md5 checksum of the file
    """
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
    subs = ['md5.txt', 'MD5.txt']
    md5_files = [path for path in Path(download_dir).rglob('*') if path.name in subs]
    #md5_lines = [(md5_file.parent, line.split()) for md5_file in md5_files for line in open(md5_file).readlines()]
    md5_lines = []
    for md5_file in md5_files:
        md5_file = Path(md5_file)
        parent_dir = md5_file.parent
        with open(md5_file) as f:
            for line in f:
                line = line.split()
                if len(line) == 2:
                    md5_lines.append((parent_dir, line))
                else:
                    print(f'WARNING: invalid md5 line in {md5_file}: {line}')
    return md5_lines
def calc_md5(input_files, file_labels=None, result_file=None, force=False, threads=8):
    """
    Calculate md5 checksums of input files.
    Input:
        input_files: list of files to calculate md5 checksums
        file_labels: rename files in md5 checksum results, for example: /dir/file -> file
        result_file: file to write md5 checksum results
    Output:
        True if all md5 checksums are OK, False otherwise
    """
    # --- create result file --- #
    if result_file is not None:
        result_file = Path(result_file)
        if result_file.exists() and not force:
            print(f'{result_file} already exists. Use force=True to overwrite.')
            return True
        result_file.parent.mkdir(parents=True, exist_ok=True)
        if result_file.exists():
            # remove and create new file
            result_file.unlink()
            result_file.touch()
        else:
            # create new file
            result_file.touch()
        
    if file_labels is None:
        # use input_files as file_labels
        file_labels = input_files
    # --- verify md5 checksums --- #
    nfails = 0
    with concurrent.futures.ProcessPoolExecutor(threads) as executor:
        future_to_md5_line = {executor.submit(calculate_md5_for_file, input_file): (input_file, file_label)
                              for input_file, file_label in zip(input_files, file_labels)}
        for future in concurrent.futures.as_completed(future_to_md5_line):
            input_file, file_label = future_to_md5_line[future]
            try:
                file, md5 = future.result() # file is the same as input_file
                with open(result_file, 'a') as f:
                    f.write(f'{md5}\t{file_label}\n')
            except Exception as exc:
                print('%r generated an exception: %s' % (input_file, exc))
                nfails += 1
    if nfails == 0:
        print('All md5 checksums are successfully calculated.')
        return True
    else:
        print(f'{nfails} md5 checksums are not successfully calculated.')
        return False
def verify_md5(download_dir, result_file, threads=8):
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
    with concurrent.futures.ProcessPoolExecutor(threads) as executor:
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
def ls_samples(download_dir, show=False):
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
def gen_sample_table(task_dirp, R1_file="_R1.fq.gz", R2_file="_R2.fq.gz"):
    """
    Create sample table for a raw directory.
    Input:
        task_dirp: directory containing samples
        R1_file: R1 file suffix
        R2_file: R2 file suffix
    Output:
        sample_table will be written to task_dirp/../sample_table.tsv
    """
    task_dirp = Path(task_dirp)
    sample_table = task_dirp.parent / "sample_table.csv"
    if sample_table.exists():
        sample_table.unlink()
    sample_table.touch()
    for sample_dirp in task_dirp.glob('*'):
        if not sample_dirp.is_dir():
            continue
        sample_name = sample_dirp.name
        R1_files = list(sample_dirp.glob(f'*{R1_file}'))
        R2_files = list(sample_dirp.glob(f'*{R2_file}'))
        if len(R1_files) == 0:
            print(f'WARNING: no R1 file for {sample_name}')
            continue
        elif len(R1_files) > 1:
            raise ValueError(f'More than one R1 file for {sample_name}: {R1_files}')
        else:
            R1_file = R1_files[0]
        if len(R2_files) == 0:
            print(f'WARNING: no R2 file for {sample_name}')
            continue
        elif len(R2_files) > 1:
            raise ValueError(f'More than one R2 file for {sample_name}: {R2_files}')
        else:
            R2_file = R2_files[0]
        with open(sample_table, 'a') as f:
            f.write(f'{sample_name},{R1_file},{R2_file}\n')
    print(f'Sample table written to {sample_table}')