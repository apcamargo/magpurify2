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
from pathlib import Path

from biolib.logger import logger_setup

import magpurify2


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
    parser.set_defaults(func=magpurify2.composition.main)
    parser.add_argument(
        "genomes", nargs="+", help="Input genomes in the FASTA format.", type=Path,
    )
    parser.add_argument(
        "output_directory", help="Directory to write the output files to.", type=Path,
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
    parser.set_defaults(func=magpurify2.coverage.main)
    parser.add_argument(
        "genomes", nargs="+", help="Input genomes in the FASTA format.", type=Path,
    )
    parser.add_argument(
        "output_directory", help="Directory to write the output files to.", type=Path,
    )
    parser.add_argument(
        "--bam_files",
        required=True,
        nargs="+",
        help="Input sorted BAM files.",
        type=Path,
    )
    parser.add_argument(
        "--min_identity",
        default=0.97,
        help="Exclude reads by overall identity to the reference sequences.",
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
    parser.set_defaults(func=magpurify2.taxonomy.main)
    parser.add_argument(
        "genomes", nargs="+", help="Input genomes in the FASTA format.", type=Path,
    )
    parser.add_argument(
        "output_directory", help="Directory to write the output files to.", type=Path,
    )
    parser.add_argument(
        "database", help="Path to MAGpurify2's database directory.", type=Path,
    )
    parser.add_argument(
        "--contig_min_fraction",
        default=0.75,
        type=float,
        help="The contig-level taxonomy must agree with the taxonomy of at least "
        "`contig_min_fraction` of its genes. This value must be equal to or greater than "
        "0.5 and less than 1.",
    )
    parser.add_argument(
        "--genome_min_fraction",
        default=0.75,
        type=float,
        help="The genome-level taxonomy must agree with the taxonomy of at least "
        "`genome_min_fraction` of its contigs (weighted by length). This value must be "
        "equal to or greater than 0.5 and less than 1.",
    )
    parser.add_argument(
        "--allow_genus",
        help="Allows genus-level taxonomic assignment.",
        action="store_true",
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
    parser.set_defaults(func=magpurify2.filter.main)
    parser.add_argument(
        "genomes", nargs="+", help="Input genomes in the FASTA format.", type=Path,
    )
    parser.add_argument(
        "output_directory",
        help="Directory with the outputs of the contaminant identification modules.",
        type=Path,
    )
    parser.add_argument(
        "filtered_output_directory",
        help="Directory where the filtered MAGs will be written to.",
        type=Path,
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
