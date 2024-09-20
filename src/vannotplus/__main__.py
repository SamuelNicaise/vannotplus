"""
@Goal: Perform VANNOT's steps who depend on ped9 : exomiser and familybarcode annotations
@Author: Samuel Nicaise
@Date: September 2024
"""

import argparse
import logging as log

import vannotplus
from vannotplus.family.barcode import main_barcode


def set_log_level(verbosity):
    verbosity = verbosity.lower()
    configs = {
        "debug": log.DEBUG,
        "info": log.INFO,
        "warning": log.WARNING,
        "error": log.ERROR,
        "critical": log.CRITICAL,
    }
    if verbosity not in configs.keys():
        raise ValueError(
            f"Unknown verbosity level: {verbosity}\nPlease use any in: {configs.keys()}"
        )
    log.basicConfig(
        format="%(asctime)s [%(levelname)s] %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
        level=configs[verbosity],
    )


def main():
    parser = argparse.ArgumentParser(prog="vannotplus")
    parser.add_argument(
        "--version",
        action="version",
        version=f"{parser.prog} {vannotplus.__version__}",
    )

    subparsers = parser.add_subparsers(help="sub-command help")

    barcode_parser = subparsers.add_parser(
        "barcode",
        help="",
        formatter_class=argparse.MetavarTypeHelpFormatter,
    )
    barcode_parser.set_defaults(subparser="barcode")

    barcode_parser.add_argument(
        "-i",
        "--input",
        type=str,
        required=True,
        help="Input VCF containing all samples of interest",
    )
    barcode_parser.add_argument(
        "-o",
        "--output",
        type=str,
        required=True,
        help="Output VCF containing all samples of interest",
    )
    barcode_parser.add_argument(
        "-p",
        "--ped-dir",
        type=str,
        required=True,
        help="Ped_raw's data directory, containing all ped files",
    )
    barcode_parser.add_argument(
        "-a",
        "--app",
        type=str,
        required=True,
        help="STARK application",
    )

    for subparser in (barcode_parser,):
        subparser.add_argument(
            "-v",
            "--verbosity",
            type=str,
            default="info",
            help="Verbosity level  [info]",
        )

    args = parser.parse_args()
    if not hasattr(args, "subparser"):
        parser.print_help()
    else:
        set_log_level(args.verbosity)
        log.debug(f"Args: {str(args)}")
        if args.subparser == "barcode":
            main_barcode(args.input, args.output, args.ped_dir, args.app)


if __name__ == "__main__":
    main()
