"""
@Goal: Expand VANNOT's capabilities with ped9 (family barcode, exomiser), GMC and scoring
@Author: Samuel Nicaise
@Date: September 2024
"""

import argparse
import logging as log
from os.path import dirname, join as osj
import shutil

import vannotplus
from vannotplus.commons import set_log_level, load_config
from vannotplus.family.barcode import main_barcode, main_barcode_fast
from vannotplus.exomiser.exomiser import main_exomiser
from vannotplus.annot.score import main_annot


def main_config(output):
    default_config = osj(dirname(vannotplus.__file__), "config.yml")
    shutil.copy(default_config, output)


def main():
    parser = argparse.ArgumentParser(prog="vannotplus")
    parser.add_argument(
        "--version",
        action="version",
        version=f"{parser.prog} {vannotplus.__version__}",
    )

    subparsers = parser.add_subparsers(help="sub-command help")

    config_parser = subparsers.add_parser(
        "config",
        help="Create a default config file",
        formatter_class=argparse.MetavarTypeHelpFormatter,
    )
    config_parser.set_defaults(subparser="config")
    config_parser.add_argument(
        "-c", "--config", default="./config.yml", help="Config file path [./config.yml]"
    )

    barcode_parser = subparsers.add_parser(
        "barcode",
        help="Compute family barcode based on Ped file",
        formatter_class=argparse.MetavarTypeHelpFormatter,
    )
    barcode_parser.set_defaults(subparser="barcode")

    exomiser_parser = subparsers.add_parser(
        "exomiser",
        help="Add exomiser scores for each sample in input VCF based on HPOs in Ped file",
        formatter_class=argparse.MetavarTypeHelpFormatter,
    )
    exomiser_parser.set_defaults(subparser="exomiser")

    score_parser = subparsers.add_parser(
        "annot",
        help="Add Gene Variations Count (GVC) for each sample and vannotscore for each variant to input VCF",
        formatter_class=argparse.MetavarTypeHelpFormatter,
    )
    score_parser.set_defaults(subparser="score")

    for subparser in (barcode_parser, exomiser_parser, score_parser):
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

    for subparser in (barcode_parser, exomiser_parser):
        subparser.add_argument(
            "-a",
            "--app",
            type=str,
            required=True,
            help="STARK application",
        )

    # separate loops to keep args in order
    default_config = osj(dirname(vannotplus.__file__), "config.yml")
    for subparser in (barcode_parser, exomiser_parser, score_parser):
        subparser.add_argument(
            "-c",
            "--config",
            type=str,
            default=default_config,
            help=f"YAML config file [{default_config}]",
        )

    score_parser.add_argument(
        "-vs",
        "--vannotscore",
        action="store_true",
        help="Compute vannotscore if set to True. If False, only GMC is computed [False]",
    )

    for subparser in (barcode_parser, exomiser_parser, score_parser, config_parser):
        subparser.add_argument(
            "-v",
            "--verbosity",
            type=str,
            default="info",
            help="Verbosity level [info]",
        )

    args = parser.parse_args()
    if not hasattr(args, "subparser"):
        parser.print_help()
    else:
        set_log_level(args.verbosity)
        log.debug(f"Args: {str(args)}")
        if args.subparser == "config":
            main_config(args.config)
        config = load_config(args.config)
        log.debug(f"config: {config}")
        if args.subparser == "barcode":
            main_barcode_fast(args.input, args.output, args.app, config)
        elif args.subparser == "exomiser":
            main_exomiser(args.input, args.output, args.app, config)
        elif args.subparser == "score":
            main_annot(args.input, args.output, config, do_vannotscore=args.vannotscore)


if __name__ == "__main__":
    main()
