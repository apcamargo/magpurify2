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
        add_help=False,
    )
    options = parser.add_argument_group("Options")
    options.add_argument(
        "--version", action="version", version=f"{parser.prog} v{magpurify2.__version__}",
    )
    options.add_argument(
        "-h", "--help", action="help", help="Show this help message and exit"
    )
    subparsers = parser.add_subparsers(title="Modules")
    composition_subcommand = subparsers.add_parser(
        "composition",
        help="Identify putative contaminants using tetranucleotide frequences.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        add_help=False,
    )
    composition_parser(composition_subcommand)
    coverage_subcommand = subparsers.add_parser(
        "coverage",
        help="Identify putative contaminants using coverage profiles.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        add_help=False,
    )
    coverage_parser(coverage_subcommand)
    taxonomy_subcommand = subparsers.add_parser(
        "taxonomy",
        help="Identify putative contaminants through taxonomy assignment.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        add_help=False,
    )
    taxonomy_parser(taxonomy_subcommand)
    filter_subcommand = subparsers.add_parser(
        "filter",
        help="Remove identified contaminants from input MAGs.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        add_help=False,
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
        elif sys.argv[1] == "taxonomy":
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
    required = parser.add_argument_group("Required arguments")
    options = parser.add_argument_group("Data processing options")
    other = parser.add_argument_group("Other options")
    required.add_argument(
        "genomes", nargs="+", help="Input genomes in the FASTA format.", type=Path,
    )
    required.add_argument(
        "output_directory", help="Directory to write the output files to.", type=Path,
    )
    options.add_argument(
        "--n_iterations",
        default=4,
        help="Number of iterations of data embedding and cluster selection with different seeds.",
        type=int,
    )
    options.add_argument(
        "--n_components",
        default=3,
        help="The dimension of the space to embed the data into (UMAP parameter).",
        type=int,
    )
    options.add_argument(
        "--min_dist",
        default=0.1,
        help="The effective minimum distance between embedded points (UMAP parameter).",
        type=float,
    )
    options.add_argument(
        "--n_neighbors",
        default=15,
        help="The size of local neighborhood used for manifold approximation (UMAP parameter).",
        type=int,
    )
    options.add_argument(
        "--set_op_mix_ratio",
        default=1.0,
        help="Interpolate between the union (1.0) and intersection (0.0) to combine the local simplicial sets (UMAP parameter).",
        type=float,
    )
    other.add_argument(
        "-t",
        "--threads",
        default=multiprocessing.cpu_count(),
        type=int,
        help="Number of threads to use. All by default.",
    )
    other.add_argument(
        "-q", "--quiet", help="Suppress the logger output", action="store_true"
    )
    other.add_argument(
        "-h", "--help", action="help", help="Show this help message and exit"
    )


def coverage_parser(parser):
    parser.set_defaults(func=magpurify2.coverage.main)
    required = parser.add_argument_group("Required arguments")
    options = parser.add_argument_group("Data processing options")
    other = parser.add_argument_group("Other options")
    required.add_argument(
        "genomes", nargs="+", help="Input genomes in the FASTA format.", type=Path,
    )
    required.add_argument(
        "output_directory", help="Directory to write the output files to.", type=Path,
    )
    required.add_argument(
        "--bam_files",
        required=True,
        nargs="+",
        help="Input sorted BAM files.",
        type=Path,
    )
    options.add_argument(
        "--min_identity",
        default=0.97,
        help="Exclude reads by overall identity to the reference sequences.",
        type=float,
    )
    options.add_argument(
        "--n_iterations",
        default=4,
        help="Number of iterations of data embedding and cluster selection with different seeds.",
        type=int,
    )
    options.add_argument(
        "--n_components",
        default=3,
        help="The dimension of the space to embed the data into (UMAP parameter).",
        type=int,
    )
    options.add_argument(
        "--min_dist",
        default=0.15,
        help="The effective minimum distance between embedded points (UMAP parameter).",
        type=float,
    )
    options.add_argument(
        "--n_neighbors",
        default=15,
        help="The size of local neighborhood used for manifold approximation (UMAP parameter).",
        type=int,
    )
    options.add_argument(
        "--set_op_mix_ratio",
        default=0.3,
        help="Interpolate between the union (1.0) and intersection (0.0) to combine the local simplicial sets (UMAP parameter).",
        type=float,
    )
    other.add_argument(
        "-t",
        "--threads",
        default=multiprocessing.cpu_count(),
        type=int,
        help="Number of threads to use. All by default.",
    )
    other.add_argument(
        "-q", "--quiet", help="Suppress the logger output", action="store_true"
    )
    other.add_argument(
        "-h", "--help", action="help", help="Show this help message and exit"
    )


def taxonomy_parser(parser):
    parser.set_defaults(func=magpurify2.taxonomy.main)
    required = parser.add_argument_group("Required arguments")
    options = parser.add_argument_group("Data processing options")
    other = parser.add_argument_group("Other options")
    required.add_argument(
        "genomes", nargs="+", help="Input genomes in the FASTA format.", type=Path,
    )
    required.add_argument(
        "output_directory", help="Directory to write the output files to.", type=Path,
    )
    required.add_argument(
        "database", help="Path to MAGpurify2's database directory.", type=Path,
    )
    options.add_argument(
        "--contig_min_fraction",
        default=0.75,
        type=float,
        help="The contig-level taxonomy must agree with the taxonomy of at least "
        "`contig_min_fraction` of its genes. This value must be equal to or greater than "
        "0.5 and less than 1.",
    )
    options.add_argument(
        "--genome_min_fraction",
        default=0.75,
        type=float,
        help="The genome-level taxonomy must agree with the taxonomy of at least "
        "`genome_min_fraction` of its contigs (weighted by length). This value must be "
        "equal to or greater than 0.5 and less than 1.",
    )
    options.add_argument(
        "--allow_genus",
        help="Allows genus-level taxonomic assignment.",
        action="store_true",
    )
    other.add_argument(
        "-t",
        "--threads",
        default=multiprocessing.cpu_count(),
        type=int,
        help="Number of threads to use. All by default.",
    )
    other.add_argument(
        "-q", "--quiet", help="Suppress the logger output", action="store_true"
    )
    other.add_argument(
        "-h", "--help", action="help", help="Show this help message and exit"
    )


def filter_parser(parser):
    parser.set_defaults(func=magpurify2.filter.main)
    required = parser.add_argument_group("Required arguments")
    options = parser.add_argument_group("Data processing options")
    other = parser.add_argument_group("Other options")
    required.add_argument(
        "genomes", nargs="+", help="Input genomes in the FASTA format.", type=Path,
    )
    required.add_argument(
        "output_directory",
        help="Directory with the outputs of the contaminant identification modules.",
        type=Path,
    )
    required.add_argument(
        "filtered_output_directory",
        help="Directory where the filtered MAGs will be written to.",
        type=Path,
    )
    options.add_argument(
        "--composition_threshold",
        default=0.15,
        type=float,
        help="Minimum score for a contig not to be flagged as a contaminant by the 'composition' module.",
    )
    options.add_argument(
        "--coverage_threshold",
        default=0.15,
        type=float,
        help="Minimum score for a contig not to be flagged as a contaminant by the 'coverage' module.",
    )
    options.add_argument(
        "--taxonomy_threshold",
        default=0.25,
        type=float,
        help="Minimum score for a contig not to be flagged as a contaminant by the 'taxonomy' module.",
    )
    options.add_argument(
        "--mode",
        default="any",
        type=str,
        choices=["any", "all"],
        help="Remove contigs that were flagged as contaminants by 'any' or by 'all' modules",
    )
    other.add_argument(
        "-t",
        "--threads",
        default=multiprocessing.cpu_count(),
        type=int,
        help="Number of threads to use. All by default.",
    )
    other.add_argument(
        "-q", "--quiet", help="Suppress the logger output", action="store_true"
    )
    other.add_argument(
        "-h", "--help", action="help", help="Show this help message and exit"
    )
