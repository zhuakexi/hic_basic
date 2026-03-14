#!/usr/bin/env python
"""
Create a synthetic BAM file for testing phasing functions.

This script generates a small BAM file with ~2000 paired-end reads
distributed across chromosomes with specific ratios:
- chr1: ~32.2% (967 reads)
- chr2: ~21.5% (645 reads)
- chr3: ~10.7% (322 reads)
- chrX: 3% (60 reads)
- chrY: 0.3% (6 reads)

Total: 2000 reads (1000 pairs)

Usage:
    python create_synthetic_bam.py

Output:
    tests/data/synthetic_pe.bam
    tests/data/synthetic_pe.bam.bai
    tests/data/synthetic_bam_README.md
"""

import pysam
import random
import os
from pathlib import Path

# Set seed for reproducibility
random.seed(42)

# Configuration
OUTPUT_DIR = Path(__file__).parent
OUTPUT_BAM = OUTPUT_DIR / "synthetic_pe.bam"
TOTAL_READS = 2000  # Total individual reads (1000 pairs)

# Chromosome distribution
# X: 3%, Y: 0.3%, rest split as chr1:chr2:chr3 = 3:2:1
CHR_X_COUNT = int(TOTAL_READS * 0.03)  # 60
CHR_Y_COUNT = int(TOTAL_READS * 0.003)  # 6
CHR_REMAINING = TOTAL_READS - CHR_X_COUNT - CHR_Y_COUNT  # 1934

# Split remaining as 3:2:1 ratio (total 6 parts)
CHR1_COUNT = int(CHR_REMAINING * 3 / 6)  # 967
CHR2_COUNT = int(CHR_REMAINING * 2 / 6)  # 645
CHR3_COUNT = CHR_REMAINING - CHR1_COUNT - CHR2_COUNT  # 322

# Verify counts add up
assert CHR1_COUNT + CHR2_COUNT + CHR3_COUNT + CHR_X_COUNT + CHR_Y_COUNT == TOTAL_READS

# Chromosome sizes (small for testing)
CHROM_SIZES = {
    "chr1": 100000,
    "chr2": 80000,
    "chr3": 60000,
    "chrX": 50000,
    "chrY": 40000,
}

# Distribution config
CHROM_DIST = [
    ("chr1", CHR1_COUNT),
    ("chr2", CHR2_COUNT),
    ("chr3", CHR3_COUNT),
    ("chrX", CHR_X_COUNT),
    ("chrY", CHR_Y_COUNT),
]

# Read configuration
READ_LENGTH = 100
FRAGMENT_SIZE_MEAN = 300
FRAGMENT_SIZE_STD = 50


def generate_read_sequence(length: int) -> str:
    """Generate a random DNA sequence."""
    return ''.join(random.choices('ACGT', k=length))


def generate_quality_string(length: int, mean_qual: int = 35) -> str:
    """Generate quality scores as Phred+33 encoded string."""
    quals = [max(0, min(93, random.gauss(mean_qual, 5))) for _ in range(length)]
    return ''.join(chr(int(q) + 33) for q in quals)


def create_read_pair(read_num: int, chrom: str, pos: int, fragment_size: int):
    """
    Create a pair of reads (read1 and read2) at the given position.

    Args:
        read_num: Unique read pair number
        chrom: Chromosome name
        pos: Start position (0-based)
        fragment_size: Fragment size (distance between read starts)

    Returns:
        Two pysam.AlignedSegment objects (read1, read2)
    """
    # Create reference header info
    header = {
        'HD': {'VN': '1.0'},
        'SQ': [{'SN': chrom, 'LN': CHROM_SIZES[chrom]}]
    }

    # We'll create the reads and add them to a file with proper header
    read_name = f"read{read_num:05d}"
    seq1 = generate_read_sequence(READ_LENGTH)
    seq2 = generate_read_sequence(READ_LENGTH)
    qual1 = generate_quality_string(READ_LENGTH)
    qual2 = generate_quality_string(READ_LENGTH)

    # Read 1: forward strand
    read1 = pysam.AlignedSegment()
    read1.query_name = read_name
    read1.query_sequence = seq1
    read1.query_qualities = [ord(c) - 33 for c in qual1]
    read1.flag = 0x41  # paired, read1, unmapped (will be updated)
    read1.reference_id = 0  # Will be set when added to file
    read1.reference_start = pos
    read1.mapping_quality = 60
    read1.cigarstring = f"{READ_LENGTH}M"
    read1.next_reference_id = 0
    read1.next_reference_start = pos + fragment_size
    read1.template_length = fragment_size

    # Read 2: reverse strand (mate position)
    read2 = pysam.AlignedSegment()
    read2.query_name = read_name
    read2.query_sequence = seq2
    read2.query_qualities = [ord(c) - 33 for c in qual2]
    read2.flag = 0x81  # paired, read2, unmapped (will be updated)
    read2.reference_id = 0
    read2.reference_start = pos + fragment_size
    read2.mapping_quality = 60
    read2.cigarstring = f"{READ_LENGTH}M"
    read2.next_reference_id = 0
    read2.next_reference_start = pos
    read2.template_length = -fragment_size

    return read1, read2


def main():
    print(f"Creating synthetic BAM file with {TOTAL_READS} reads ({TOTAL_READS//2} pairs)")
    print(f"Chromosome distribution:")
    for chrom, count in CHROM_DIST:
        pct = count / TOTAL_READS * 100
        print(f"  {chrom}: {count} reads ({pct:.1f}%)")

    # Collect all reads with their chromosome info
    all_reads = []
    for chrom, count in CHROM_DIST:
        for i in range(count):
            # Random position ensuring reads fit on chromosome
            max_pos = CHROM_SIZES[chrom] - FRAGMENT_SIZE_MEAN - READ_LENGTH
            pos = random.randint(0, max_pos)
            fragment_size = int(random.gauss(FRAGMENT_SIZE_MEAN, FRAGMENT_SIZE_STD))
            fragment_size = max(100, min(1000, fragment_size))  # Clamp to reasonable range
            all_reads.append((chrom, pos, fragment_size))

    # Shuffle to mix chromosomes
    random.shuffle(all_reads)

    # Create BAM file with full header
    header = {
        'HD': {'VN': '1.0', 'SO': 'unsorted'},
        'SQ': [{'SN': chrom, 'LN': size} for chrom, size in CHROM_SIZES.items()]
    }

    read_num = 0
    # Build chromosome to ref_id mapping
    chrom_to_ref_id = {chrom: i for i, chrom in enumerate(CHROM_SIZES.keys())}

    with pysam.AlignmentFile(OUTPUT_BAM, 'wb', header=header) as bam_out:
        for chrom, pos, fragment_size in all_reads:
            read1, read2 = create_read_pair(read_num, chrom, pos, fragment_size)

            # Get reference ID for this chromosome
            ref_id = chrom_to_ref_id[chrom]
            read1.reference_id = ref_id
            read2.reference_id = ref_id
            read1.next_reference_id = ref_id
            read2.next_reference_id = ref_id

            # Update flags for proper pairing
            # Read 1: paired (0x1), proper pair (0x2), read1 (0x40)
            read1.flag = 0x41 | 0x2
            # Read 2: paired (0x1), proper pair (0x2), read2 (0x80), reverse (0x10)
            read2.flag = 0x81 | 0x2 | 0x10

            bam_out.write(read1)
            bam_out.write(read2)
            read_num += 1

    # Sort and index the BAM file
    print(f"\nSorting and indexing BAM file...")
    sorted_bam = OUTPUT_BAM.with_suffix('.sorted.bam')
    pysam.sort("-o", str(sorted_bam), str(OUTPUT_BAM))
    os.rename(sorted_bam, OUTPUT_BAM)  # Replace unsorted with sorted

    pysam.index(str(OUTPUT_BAM))

    # Verify the file
    print(f"\nVerifying BAM file...")
    with pysam.AlignmentFile(OUTPUT_BAM, 'rb') as bam_in:
        total = 0
        chrom_counts = {}
        for read in bam_in:
            if not read.is_unmapped:
                total += 1
                chrom = bam_in.get_reference_name(read.reference_id)
                chrom_counts[chrom] = chrom_counts.get(chrom, 0) + 1

        print(f"Total reads in file: {total}")
        print(f"Chromosome distribution:")
        for chrom in sorted(chrom_counts.keys()):
            count = chrom_counts[chrom]
            pct = count / total * 100
            print(f"  {chrom}: {count} ({pct:.2f}%)")

    # Create README (use escaped braces for literal {{}} in markdown)
    readme_content = f"""# Synthetic BAM Test File

## Overview
This BAM file (`synthetic_pe.bam`) contains synthetic paired-end reads generated for testing phasing and chromosome analysis functions.

## Generation Parameters
- **Total reads**: {TOTAL_READS} ({TOTAL_READS//2} read pairs)
- **Read length**: {READ_LENGTH} bp
- **Fragment size**: ~{FRAGMENT_SIZE_MEAN} bp (±{FRAGMENT_SIZE_STD} bp)
- **Random seed**: 42 (reproducible)

## Chromosome Distribution
| Chromosome | Read Count | Percentage |
|------------|------------|------------|
| chr1       | {CHR1_COUNT:>6}   | {CHR1_COUNT/TOTAL_READS*100:.2f}% |
| chr2       | {CHR2_COUNT:>6}   | {CHR2_COUNT/TOTAL_READS*100:.2f}% |
| chr3       | {CHR3_COUNT:>6}   | {CHR3_COUNT/TOTAL_READS*100:.2f}% |
| chrX       | {CHR_X_COUNT:>6}   | {CHR_X_COUNT/TOTAL_READS*100:.2f}% |
| chrY       | {CHR_Y_COUNT:>6}   | {CHR_Y_COUNT/TOTAL_READS*100:.2f}% |

## Expected Ratios
- X chromosome: 3% of total reads
- Y chromosome: 0.3% of total reads
- Autosomal (chr1:chr2:chr3): 3:2:1 ratio for remaining reads

## File Information
- **Format**: BAM (binary SAM)
- **Index**: synthetic_pe.bam.bai
- **Sort order**: Coordinate sorted
- **Reference**: Custom 5-chromosome reference (chr1, chr2, chr3, chrX, chrY)

## Usage Example
```python
import pysam

# Open BAM file
bam = pysam.AlignmentFile("synthetic_pe.bam", "rb")

# Count reads per chromosome
chrom_counts = {{}}
for read in bam.fetch():
    chrom = bam.get_reference_name(read.reference_id)
    chrom_counts[chrom] = chrom_counts.get(chrom, 0) + 1

print(chrom_counts)
```

## Generated By
```bash
python create_synthetic_bam.py
```

## Notes
- This is a synthetic file for testing purposes only
- Read sequences are random (not biologically meaningful)
- All reads map to a single chromosome (no chimeric reads)
- Reads are properly paired with realistic fragment sizes
"""

    readme_path = OUTPUT_DIR / "synthetic_bam_README.md"
    with open(readme_path, 'w') as f:
        f.write(readme_content)

    print(f"\nFiles created:")
    print(f"  - {OUTPUT_BAM}")
    print(f"  - {OUTPUT_BAM}.bai")
    print(f"  - {readme_path}")
    print("\nDone!")


if __name__ == "__main__":
    main()
