import sys
sys.path.insert(0, "/share/home/ychi/dev/hic_basic")
import unittest
from pathlib import Path
from hic_basic.wet.rawcheck import verify_md5

class TestVerifyMd5(unittest.TestCase):
    def test_verify_md5(self):
        download_dir = "/sharec/ychi/downloads/1208"
        result_file = Path(__file__).parent / "output/wet/md5.result"

        # Run the verify_md5 function
        verify_md5(download_dir, str(result_file))

        # Count the number of lines in the result file
        with open(result_file, 'r') as f:
            result_lines = sum(1 for line in f)

        # Count the total number of lines in all md5.txt files
        md5_files = Path(download_dir).glob('**/md5.txt')
        md5_lines = sum(1 for md5_file in md5_files for line in open(md5_file))

        # Assert that the number of lines in the result file equals the total number of lines in all md5.txt files
        self.assertEqual(result_lines, md5_lines)

if __name__ == '__main__':
    unittest.main()