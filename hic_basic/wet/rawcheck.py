import hashlib
from pathlib import Path

def calculate_md5(file_path):
    md5 = hashlib.md5()
    with open(file_path, "rb") as f:
        # 以块的方式更新 MD5 值，适用于大文件
        for chunk in iter(lambda: f.read(4096), b""):
            md5.update(chunk)
    return md5.hexdigest()
def verify_md5(download_dir, result_file):
    """
    Verify md5 checksums of downloaded files.
    """
    # --- find all md5.txt files in download_dir --- #
    md5_files = Path(download_dir).glob('**/md5.txt')
    parent_dirs = [md5_file.parent for md5_file in md5_files]
    
    # --- verify md5 checksums --- #
    nfails = 0
    for parent_dir in parent_dirs:
        md5_file = parent_dir / 'md5.txt'
        with open(md5_file, 'r') as f:
            md5_lines = f.readlines()
        md5_results = []
        for md5_line in md5_lines:
            md5, filename = md5_line.split()
            filepath = parent_dir / filename
            if md5 != calculate_md5(filepath):
                print(f'ERROR: md5 checksum of {filepath} does not match.')
                nfails += 1
                md5_results.append(filename + ' ERROR')
            else:
                md5_results.append(filename + ' OK')
        md5_results = '\n'.join(md5_results)
        with open(result_file, 'a') as f:
            f.write(md5_results)
            f.write('\n')
    if nfails == 0:
        print('All md5 checksums are OK.')
        return True
    else:
        print(f'{nfails} md5 checksums are not OK.')
        return False