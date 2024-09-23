"""
@Goal: Perform VANNOT's steps who depend on ped9 : exomiser and familybarcode annotations
@Author: Samuel Nicaise
@Date: September 2024
"""

import argparse
import logging as log

import vannotplus
from vannotplus.commons import set_log_level, load_config
from vannotplus.family.barcode import main_barcode
from vannotplus.exomiser.exomiser import main_exomiser


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

    exomiser_parser = subparsers.add_parser(
        "exomiser",
        help="",
        formatter_class=argparse.MetavarTypeHelpFormatter,
    )
    exomiser_parser.set_defaults(subparser="exomiser")

    for subparser in (barcode_parser, exomiser_parser):
        subparser.add_argument(
            "-i",
            "--input",
            type=str,
            required=True,
            help="Input VCF containing all samples of interest",
        )
        subparser.add_argument(
            "-o",
            "--output",
            type=str,
            required=True,
            help="Output VCF containing all samples of interest",
        )
        subparser.add_argument(
            "-a",
            "--app",
            type=str,
            required=True,
            help="STARK application",
        )
        subparser.add_argument(
            "-c",
            "--config",
            type=str,
            required=True,
            help="YAML config file",
        )
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
        config = load_config(args.config)
        log.debug(f"config: {config}")
        if args.subparser == "barcode":
            main_barcode(args.input, args.output, args.app, config)
        elif args.subparser == "exomiser":
            main_exomiser(args.input, args.output, args.app, config)


if __name__ == "__main__":
    main()
