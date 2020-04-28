# -*- coding: utf-8 -*-

# This file is part of the magpurify2 package, available at:
# https://github.com/apcamargo/magpurify2
#
# Magpurify2 is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program. If not, see <https://www.gnu.org/licenses/>.
#
# Contact: antoniop.camargo@gmail.com

import argparse
import multiprocessing
import sys

from biolib.logger import logger_setup

import magpurify2
from magpurify2.main import (
    composition_module,
    coverage_module,
    taxonomy_module,
    filter_module,
)


def cli():
    parser = argparse.ArgumentParser(
        description="Identify and remove contaminants from metagenome-assembled genomes.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "--version", action="version", version=f"{parser.prog} v{magpurify2.__version__}",
    )
    subparsers = parser.add_subparsers()
    composition_subcommand = subparsers.add_parser(
        "composition",
        help="Identify putative contaminants using tetranucleotide frequences.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    composition_parser(composition_subcommand)
    coverage_subcommand = subparsers.add_parser(
        "coverage",
        help="Identify putative contaminants using coverage profiles.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    coverage_parser(coverage_subcommand)

    taxonomy_subcommand = subparsers.add_parser(
        "taxonomy",
        help="Identify putative contaminants through taxonomy assignment.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    taxonomy_parser(taxonomy_subcommand)

    filter_subcommand = subparsers.add_parser(
        "filter",
        help="Remove identified contaminants from input MAGs.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    filter_parser(filter_subcommand)
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(0)
    elif len(sys.argv) == 2:
        if sys.argv[1] == "composition":
            composition_subcommand.print_help()
            sys.exit(0)
        elif sys.argv[1] == "coverage":
            coverage_subcommand.print_help()
            sys.exit(0)
        elif sys.argv[1] == "coding_density":
            taxonomy_subcommand.print_help()
            sys.exit(0)
        elif sys.argv[1] == "filter":
            filter_subcommand.print_help()
            sys.exit(0)
    args = parser.parse_args()
    logger_setup(
        args.output_directory,
        f"{sys.argv[1]}.log",
        "MAGpurify2",
        magpurify2.__version__,
        args.quiet,
    )
    args.func(args)


def composition_parser(parser):
    parser.set_defaults(func=composition_module)
    parser.add_argument("genomes", nargs="+", help="Input genomes in the FASTA format.")
    parser.add_argument(
        "output_directory", help="Directory to write the output files to.",
    )
    parser.add_argument(
        "-s",
        "--strictness",
        help=(
            "Strictness of the contaminant detection algorithm. "
            "Must be a number between 0 (less strict) and 1 (more strict)."
        ),
        default=0.5,
        type=float,
    )
    parser.add_argument(
        "-t",
        "--threads",
        default=multiprocessing.cpu_count(),
        type=int,
        help="Number of threads to use. All by default.",
    )
    parser.add_argument("--quiet", help="Suppress the logger output", action="store_true")


def coverage_parser(parser):
    parser.set_defaults(func=coverage_module)
    parser.add_argument("genomes", nargs="+", help="Input genomes in the FASTA format.")
    parser.add_argument(
        "output_directory", help="Directory to write the output files to.",
    )
    parser.add_argument(
        "--bam_files", required=True, nargs="+", help="Input sorted BAM files."
    )
    parser.add_argument(
        "-s",
        "--strictness",
        help=(
            "Strictness of the contaminant detection algorithm. "
            "Must be a number between 0 (less strict) and 1 (more strict)."
        ),
        default=0.5,
        type=float,
    )
    parser.add_argument(
        "-t",
        "--threads",
        default=multiprocessing.cpu_count(),
        type=int,
        help="Number of threads to use. All by default.",
    )
    parser.add_argument("--quiet", help="Suppress the logger output", action="store_true")


def taxonomy_parser(parser):
    parser.set_defaults(func=taxonomy_module)
    parser.add_argument("genomes", nargs="+", help="Input genomes in the FASTA format.")
    parser.add_argument(
        "output_directory", help="Directory to write the output files to.",
    )
    parser.add_argument("database", help="Path to MAGpurify2's database directory.")
    parser.add_argument(
        "-s",
        "--strictness",
        help=(
            "Strictness of the contaminant detection algorithm. "
            "Must be a number between 0 (less strict) and 1 (more strict)."
        ),
        default=0.5,
        type=float,
    )
    parser.add_argument(
        "-t",
        "--threads",
        default=multiprocessing.cpu_count(),
        type=int,
        help="Number of threads to use. All by default.",
    )
    parser.add_argument("--quiet", help="Suppress the logger output", action="store_true")


def filter_parser(parser):
    parser.set_defaults(func=filter_module)
    parser.add_argument("genomes", nargs="+", help="Input genomes in the FASTA format.")
    parser.add_argument(
        "output_directory",
        help="Directory with the outputs of the contaminant identification modules.",
    )
    parser.add_argument(
        "filtered_output_directory",
        help="Directory where the filtered MAGs will be written to.",
    )
    parser.add_argument(
        "-f",
        "--filtering",
        default="any",
        type=str,
        choices=["any", "all"],
        help="Remove contigs that were flagged as contaminants in 'any' or in 'all' modules",
    )
    parser.add_argument(
        "-t",
        "--threads",
        default=multiprocessing.cpu_count(),
        type=int,
        help="Number of threads to use. All by default.",
    )
    parser.add_argument("--quiet", help="Suppress the logger output", action="store_true")