name = "download"
description = "Download data from usual sources like GEO."

from ..data import download_geo

def add_arguments(subparser):
    subparser.add_argument(
        "-a", "--accession",
        type=str,
        required=True,
        help="GEO accession ID (e.g., GSE123456)."
    )
    subparser.add_argument(
        "-o", "--output",
        type=str,
        required=True,
        help="Output directory to save the downloaded data."
    )
    subparser.add_argument(
        "--method",
        choices=["ftp", "https"],
        default="ftp",
        help="Method to download data (default: ftp)."
    )

def run(args):
    download_geo(
        geo_id = args.accession,
        outdir = args.output,
        method = args.method
    )