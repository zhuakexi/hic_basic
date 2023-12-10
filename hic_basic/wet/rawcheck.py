import concurrent.futures
import hashlib
from pathlib import Path

def calculate_md5_for_file(file_path):
    md5 = hashlib.md5()
    with open(file_path, "rb") as f:
        for chunk in iter(lambda: f.read(4096), b""):
            md5.update(chunk)
    return file_path, md5.hexdigest()

def verify_md5(download_dir, result_file):
    """
    Verify md5 checksums of downloaded files.
    Input:
        download_dir: directory containing downloaded files
        result_file: file to write md5 checksum results
    Output:
        True if all md5 checksums are OK, False otherwise
    """
    # --- find all md5.txt files in download_dir --- #
    md5_files = Path(download_dir).glob('**/md5.txt')
    md5_lines = [(md5_file.parent, line.split()) for md5_file in md5_files for line in open(md5_file).readlines()]
    
    # --- verify md5 checksums --- #
    nfails = 0
    with concurrent.futures.ProcessPoolExecutor(8) as executor:
        future_to_md5_line = {executor.submit(calculate_md5_for_file, parent_dir / filename): (parent_dir, filename, expected_md5) 
                              for parent_dir, (expected_md5, filename) in md5_lines}
        for future in concurrent.futures.as_completed(future_to_md5_line):
            parent_dir, filename, expected_md5 = future_to_md5_line[future]
            try:
                file, md5 = future.result()
                md5_results = filename + (' OK' if md5 == expected_md5 else ' ERROR')
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