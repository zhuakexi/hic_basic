import json
import unittest
from pathlib import Path

import numpy as np
from unittest.mock import patch, MagicMock

from hic_basic.atac import (
    calculate_fragment_lengths,
    plot_fragment_length_distribution,
    get_tss_regions,
    calculate_coverage,
    normalize_coverage
)

class TestATACFunctions(unittest.TestCase):
    def setUp(self):
        # Setup code if needed
        self.chrom_size = "/share/home/ychi/dev/hic_basic/hic_basic/ref/hg38.chrom.sizes.txt"
        self.bam_path = "/shared/ychi/repo/hsv_15/mapping/E250083200_L01_DYJ-kd3h.filtered_markeddup.bam"
        self.bam_path1 = "/shared/ychi/repo/hsv_15/mapping/E250083200_L01_DYJ-wt3h.filtered_markeddup.bam"
        self.gtf = "/share/Data/ychi/genome/GRCh38/gencode.v33.primary_assembly.annotation.gtf"
        self.out_dir = Path(__file__).parent / "output" / "atac"
        if not self.out_dir.exists():
            self.out_dir.mkdir(parents=True)
    # def test_calculate_fragment_lengths(self):
    #     fragment_lengths = calculate_fragment_lengths(
    #         self.bam_path,
    #         sampling = 10000
    #     )
    #     with open(self.out_dir / "fragment_lengths.json", "w") as f:
    #         json.dump(fragment_lengths.tolist(), f)
    #     print(f"Fragment lengths: {fragment_lengths}")
    def test_plot_fragment_length_distribution(self):
        fragment_lengths = np.random.randint(50, 5000, size=1000)
        fig = plot_fragment_length_distribution(fragment_lengths)
        fig.write_image(self.out_dir / "fragment_length_distribution.png")
        print("Fragment length distribution plot saved.")
    def test_plot_fragment_length_distribution_real(self):
        with open(self.out_dir / "fragment_lengths.json", "r") as f:
            fragment_lengths = np.array(json.load(f))
        fig = plot_fragment_length_distribution(
            fragment_lengths,
            range_x = [0, 800]
            )
        fig.write_image(self.out_dir / "fragment_length_distribution_real.png")
        print("Real fragment length distribution plot saved.")
    def test_calculate_tss_enrichment(self):
        tss_ext = get_tss_regions(
            "/share/home/ychi/dev/hic_basic/hic_basic/ref/cache/GRCh38.TSS.bed",
            self.chrom_size, window_size=2000)
        coverage_array, _ = calculate_coverage(
            self.bam_path, tss_ext, 150, sampling=100
            )
        normalized_array = normalize_coverage(
            coverage_array,
            read_length,
            greenleaf_norm=True)
if __name__ == "__main__":
    unittest.main()