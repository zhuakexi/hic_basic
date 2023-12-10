import sys
sys.path.insert(0, "/share/home/ychi/dev/hic_basic")
import unittest
from pathlib import Path
from hic_basic.wet.rawcheck import verify_md5, ls_samples, split_download_dir

class TestVerifyMd5(unittest.TestCase):
    def setUp(self) -> None:
        #self.download_dir = "/sharec/ychi/downloads/1208"
        self.download_dir = "/sharec/ychi/downloads/1208/PM-XS01KF2022110123-01/ANNO_XS01KF2022110123_PM-XS01KF2022110123-01_2023-12-07_21-29-50_222LFTLT4"
    def test_verify_md5(self):
        download_dir = self.download_dir
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
    def test_ls_samples(self):
        download_dir = self.download_dir
        samples = ls_samples(download_dir, show=True)
        self.assertGreater(len(samples), 0)
    def test_split_download_dir(self):
        download_dir = self.download_dir
        split_download_dir(
            download_dir,
            {
                "/sharec/ychi/raw/sperm44_mESC/Rawdata" : [f"ESCO100{i}" for i in [2,4,5,7,8]],
                "/sharec/ychi/raw/sperm45_bulk_mESC/Rawdata" : ["ESCOBulk"],
                "/sharec/ychi/raw/sperm46_bulk_GM/Rawdata" : ["GMOBulk"]
            }
            )
        self.assertTrue(True)
if __name__ == '__main__':
    unittest.main()