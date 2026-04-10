"""Built-in file format registry for bioinformatics file types.

Each entry contains:
    - format: Human-readable name
    - extension: Common file extensions
    - description: Short description
    - reader: hic_basic function (module.function) for reading
    - writer: hic_basic function (module.function) for writing

Usage:
    >>> from hic_basic.data.file_formats import FILE_FORMATS, get_format
    >>> get_format("cool")
    {'format': 'Cooler', 'extension': '.cool/.mcool', ...}
"""

FILE_FORMATS: dict[str, dict] = {
    # Hi-C / chromatin contact matrices
    "cool": {
        "format": "Cooler",
        "extension": ".cool/.mcool",
        "description": "Sparse contact matrix stored in HDF5 format (cooler ecosystem). "
                       ".cool is single-resolution; .mcool stores multiple resolutions.",
        "reader": "cooler.Cooler",
        "writer": "cooler.create_cooler; hic_basic.coolstuff.pairs2cool",
    },
    "scool": {
        "format": "sCool",
        "extension": ".scool",
        "description": "Single-cell Cooler format. Stores multiple single-cell "
                       "cool objects within one HDF5 file (one cell = one group).",
        "reader": "cooler.fileops.open_scool; cooler.ScoolReader",
        "writer": "cooler.create_scool; hic_basic.coolstuff.pairs2scool; hic_basic.coolstuff.cools2scool",
    },
    "pairs": {
        "format": "Pairs",
        "extension": ".pairs/.pairsam",
        "description": "Contact pairs in TSV format (4DN/Pairix standard). "
                       "Header lines start with '#' containing chromosome sizes and column names.",
        "reader": "hic_basic.hicio.parse_pairs",
        "writer": "N/A (typically produced by alignment pipelines)",
    },
    "bedpe": {
        "format": "BEDPE",
        "extension": ".bedpe",
        "description": "Paired-end BED format for chromatin interactions. "
                       "Two genomic intervals per row (chrom1,start1,end1, chrom2,start2,end2).",
        "reader": "pandas.read_csv (sep='\\t'); cooler.core.read_table",
        "writer": "pandas.DataFrame.to_csv (sep='\\t')",
    },
    # Alignment
    "bam": {
        "format": "BAM",
        "extension": ".bam",
        "description": "Binary alignment format, sorted by coordinate, indexed (.bai).",
        "reader": "pysam.AlignmentFile; hic_basic.sam (module)",
        "writer": "N/A (produced by aligners such as bwa, bowtie2)",
    },
    "sam": {
        "format": "SAM",
        "extension": ".sam",
        "description": "Text alignment format (Sequence Alignment/Map).",
        "reader": "hic_basic.hicio.parse_sam_line",
        "writer": "N/A (produced by aligners)",
    },
    "cram": {
        "format": "CRAM",
        "extension": ".cram",
        "description": "Compressed reference-based alignment format.",
        "reader": "pysam.AlignmentFile",
        "writer": "N/A (produced by samtools view -C)",
    },
    # Genomic intervals
    "bed": {
        "format": "BED",
        "extension": ".bed",
        "description": "Browser Extensible Data format for genomic intervals.",
        "reader": "hic_basic.hicio.parse_bed",
        "writer": "pandas.DataFrame.to_csv (sep='\\t', header=False, index=False)",
    },
    "bedgraph": {
        "format": "BEDGraph",
        "extension": ".bedgraph/.bedGraph",
        "description": "BED format with score column for continuous signal.",
        "reader": "pandas.read_csv (sep='\\s+')",
        "writer": "pandas.DataFrame.to_csv (sep='\\t', header=False, index=False)",
    },
    "bigwig": {
        "format": "bigWig",
        "extension": ".bw/.bigwig",
        "description": "Binary indexed signal track format (UCSC).",
        "reader": "pybigwig.open; hic_basic.plot.track.BigWigTrack",
        "writer": "pyBigWig (module); bedGraphToBigWig (UCSC tool)",
    },
    "gff": {
        "format": "GFF/GTF",
        "extension": ".gff/.gff3/.gtf",
        "description": "General Feature Format / Gene Transfer Format for annotations.",
        "reader": "hic_basic.hicio.parse_gtf; hic_basic.hicio.parse_gff",
        "writer": "N/A (typically from annotation databases)",
    },
    # Variant / genotype
    "vcf": {
        "format": "VCF",
        "extension": ".vcf/.vcf.gz",
        "description": "Variant Call Format for SNPs/indels/structural variants.",
        "reader": "hic_basic.hicio.parse_vcf; cyvcf2.VCF",
        "writer": "N/A (produced by variant callers)",
    },
    # Sequence
    "fasta": {
        "format": "FASTA",
        "extension": ".fa/.fasta/.fna",
        "description": "Sequence text format (nucleotide or protein).",
        "reader": "Bio.SeqIO.read (BioPython); pysam.FastaFile",
        "writer": "Bio.SeqIO.write (BioPython)",
    },
    "fastq": {
        "format": "FASTQ",
        "extension": ".fq/.fastq",
        "description": "Sequencing reads with quality scores.",
        "reader": "hic_basic.sequence.read_fastq_ids; Bio.SeqIO.parse (BioPython)",
        "writer": "Bio.SeqIO.write (BioPython)",
    },
    # Matrix / table
    "tsv": {
        "format": "TSV",
        "extension": ".tsv/.tab",
        "description": "Tab-separated values table.",
        "reader": "pandas.read_csv (sep='\\t')",
        "writer": "pandas.DataFrame.to_csv (sep='\\t', index=False)",
    },
    "csv": {
        "format": "CSV",
        "extension": ".csv",
        "description": "Comma-separated values table.",
        "reader": "pandas.read_csv",
        "writer": "pandas.DataFrame.to_csv (index=False)",
    },
}


def get_format(key: str) -> dict | None:
    """Look up a file format definition by key.

    Parameters
    ----------
    key : str
        Format key (e.g., ``"cool"``, ``"bam"``).

    Returns
    -------
    dict or None
        Format definition dict, or None if not found.
    """
    return FILE_FORMATS.get(key)


def list_formats() -> list[str]:
    """Return sorted list of all registered format keys."""
    return sorted(FILE_FORMATS)
