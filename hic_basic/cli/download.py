name = "download"
description = "Download data from GEO or SRA."

import os

from ..data import download_geo
from ..data import download_sra


def _parse_accessions(accession_arg: str):
    """Parse an accession argument.

    Supports:
      - comma-separated list: "SRR1,SRR2"
      - a plain accession: "SRR123"
      - a text file path: each non-empty, non-comment line is an accession
    """
    if accession_arg is None:
        return []
    accession_arg = str(accession_arg).strip()
    if accession_arg == "":
        return []

    # If it looks like a file, read accessions from it
    if os.path.isfile(accession_arg):
        accessions = []
        with open(accession_arg, "rt") as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith("#"):
                    continue
                accessions.append(line)
        return accessions

    # Otherwise treat as comma-separated values
    return [a.strip() for a in accession_arg.split(",") if a.strip()]

def add_arguments(subparser):
    subparser.add_argument(
        "-m", "--mode",
        choices=["geo", "sra"],
        default="geo",
        help="Download mode: 'geo' for GEO series/samples; 'sra' for SRA run accessions (SRR...)."
    )
    subparser.add_argument(
        "-a", "--accession",
        type=str,
        required=True,
        help=(
            "Accession: GEO ID (e.g., GSE123456) for geo mode; for sra mode either "
            "(1) a run accession (SRR...), (2) a comma-separated list, or (3) a text file path "
            "with one accession per line."
        )
    )
    subparser.add_argument(
        "-o", "--output",
        type=str,
        required=True,
        help="Output directory to save the downloaded data."
    )
    subparser.add_argument(
        "--method",
        type=str,
        default=None,
        help=(
            "Download method. For geo: 'ftp' or 'https' (default ftp). "
            "For sra: 'ena' (default) or 'sra-tools' to use prefetch/fasterq-dump."
        )
    )

def run(args):
    if args.mode == "geo":
        method = args.method or "ftp"
        download_geo(
            geo_id=args.accession,
            outdir=args.output,
            method=method
        )
    elif args.mode == "sra":
        method = args.method or "ena"
        accessions = _parse_accessions(args.accession)
        if len(accessions) == 0:
            raise ValueError("No SRA accessions found. Provide SRR..., a comma-separated list, or a txt file.")
        download_sra(
            accessions=accessions,
            outdir=args.output,
            method=method
        )
    else:
        raise ValueError("mode must be 'geo' or 'sra'")