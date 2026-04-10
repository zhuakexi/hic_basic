"""Built-in file format registry for bioinformatics file types.

Usage:
    >>> from hic_basic.data.file_formats import FILE_FORMATS, get_format
    >>> get_format("cool")
    {'format': 'Cooler', 'extension': '.cool/.mcool', ...}
"""

FILE_FORMATS: dict[str, dict] = {
    # Hi-C / chromatin
    "cool": {
        "format": "Cooler",
        "extension": ".cool/.mcool",
        "description": "Sparse contact matrix stored in HDF5 format (cooler ecosystem)",
    },
    "scool": {
        "format": "sCool",
        "extension": ".scool",
        "description": "Multi-resolution cooler format (HDF5)",
    },
    "pairs": {
        "format": "Pairs",
        "extension": ".pairs/.pairsam",
        "description": "Contact pairs in TSV format (4DN/Pairix standard)",
    },
    "bedpe": {
        "format": "BEDPE",
        "extension": ".bedpe",
        "description": "Paired-end BED format for chromatin interactions",
    },
    # Alignment
    "bam": {
        "format": "BAM",
        "extension": ".bam",
        "description": "Binary alignment format, sorted by coordinate, indexed (.bai)",
    },
    "sam": {
        "format": "SAM",
        "extension": ".sam",
        "description": "Text alignment format (Sequence Alignment/Map)",
    },
    "cram": {
        "format": "CRAM",
        "extension": ".cram",
        "description": "Compressed reference-based alignment format",
    },
    # Genomic intervals
    "bed": {
        "format": "BED",
        "extension": ".bed",
        "description": "Browser Extensible Data format for genomic intervals",
    },
    "bedgraph": {
        "format": "BEDGraph",
        "extension": ".bedgraph/.bedGraph",
        "description": "BED format with score column for continuous signal",
    },
    "bigwig": {
        "format": "bigWig",
        "extension": ".bw/.bigwig",
        "description": "Binary indexed signal track format (UCSC)",
    },
    "gff": {
        "format": "GFF/GTF",
        "extension": ".gff/.gff3/.gtf",
        "description": "General Feature Format / Gene Transfer Format for annotations",
    },
    # Variant / genotype
    "vcf": {
        "format": "VCF",
        "extension": ".vcf/.vcf.gz",
        "description": "Variant Call Format for SNPs/indels/structural variants",
    },
    # Sequence
    "fasta": {
        "format": "FASTA",
        "extension": ".fa/.fasta/.fna",
        "description": "Sequence text format (nucleotide or protein)",
    },
    "fastq": {
        "format": "FASTQ",
        "extension": ".fq/.fastq",
        "description": "Sequencing reads with quality scores",
    },
    # Matrix / table
    "tsv": {
        "format": "TSV",
        "extension": ".tsv/.tab",
        "description": "Tab-separated values table",
    },
    "csv": {
        "format": "CSV",
        "extension": ".csv",
        "description": "Comma-separated values table",
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
