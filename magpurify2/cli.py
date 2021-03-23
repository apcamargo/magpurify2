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

default_values = {
    "composition": {
        "n_iterations": 4,
        "n_components": 3,
        "min_dist": 0.1,
        "n_neighbors": 15,
        "set_op_mix_ratio": 1.0,
    },
    "coverage": {
        "min_identity": 0.97,
        "trim_lower": 0.05,
        "trim_upper": 0.05,
        "min_average_coverage": 1.0,
        "n_iterations": 4,
        "n_components": 3,
        "min_dist": 0.15,
        "n_neighbors": 15,
        "set_op_mix_ratio": 0.6,
    },
    "codon_usage": {
        "n_iterations": 4,
        "n_components": 3,
        "min_dist": 0.15,
        "n_neighbors": 15,
        "set_op_mix_ratio": 0.6,
    },
    "taxonomy": {
        "contig_min_fraction": 0.5,
        "genome_min_fraction": 0.5,
        "min_genus_identity": 0.83,
    },
    "filter": {"probability_threshold": 0.5, "suffix": "filtered"},
}


def cli():
    parser = argparse.ArgumentParser(
        description="Identify and remove contaminants from metagenome-assembled genomes.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        add_help=False,
    )
    subparsers = parser.add_subparsers(title="Modules")
    end_to_end_subcommand = subparsers.add_parser(
        "end_to_end",
        help="Execute the complete MAGpurify2 pipeline to identify and remove putative "
        "contaminants from MAGs.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        add_help=False,
    )
    end_to_end_parser(end_to_end_subcommand)
    composition_subcommand = subparsers.add_parser(
        "composition",
        help="Identify putative contaminants using tetranucleotide frequencies.",
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
    codon_usage_subcommand = subparsers.add_parser(
        "codon_usage",
        help="Identify putative contaminants using gene codon usage profiles.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        add_help=False,
    )
    codon_usage_parser(codon_usage_subcommand)
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
    options = parser.add_argument_group("Options")
    options.add_argument(
        "--version", action="version", version=f"{parser.prog} v{magpurify2.__version__}",
    )
    options.add_argument(
        "-h", "--help", action="help", help="Show this help message and exit"
    )

    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(0)
    elif len(sys.argv) == 2:
        if sys.argv[1] == "end_to_end":
            end_to_end_subcommand.print_help()
            sys.exit(0)
        elif sys.argv[1] == "composition":
            composition_subcommand.print_help()
            sys.exit(0)
        elif sys.argv[1] == "coverage":
            coverage_subcommand.print_help()
            sys.exit(0)
        elif sys.argv[1] == "codon_usage":
            codon_usage_subcommand.print_help()
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


def end_to_end_parser(parser):
    parser.set_defaults(func=magpurify2.end_to_end.main)
    required = parser.add_argument_group("Required arguments")
    coverage_input = parser.add_argument_group(
        "Coverage data input (required, mutually exclusive arguments)",
    )
    coverage_input_mex = coverage_input.add_mutually_exclusive_group(required=True)
    database = parser.add_argument_group(
        "Database input (required when not using fast mode)"
    )
    options = parser.add_argument_group("Data processing options")
    other = parser.add_argument_group("Other options")
    required.add_argument(
        "genomes",
        nargs="+",
        help="Input genomes in the FASTA format. '.bz2', '.gz' and '.xz' compressed "
        "files are also supported.",
        type=Path,
    )
    required.add_argument(
        "output_directory",
        help="Directory where the logs, temporary and supporting data, and "
        "final score files will be written to.",
        type=Path,
    )
    required.add_argument(
        "filtered_output_directory",
        help="Directory where the filtered MAGs will be written to.",
        type=Path,
    )
    coverage_input_mex.add_argument(
        "--bam_files", nargs="+", help="Input sorted BAM files.", type=Path,
    )
    coverage_input_mex.add_argument(
        "--coverage_file",
        help="Input a tabular file containing contig coverage data.",
        type=Path,
    )
    database.add_argument(
        "--taxonomy_database",
        help="Path to the MAGpurify2 taxonomy database directory.",
        type=Path,
    )
    options.add_argument(
        "--checkm_file",
        type=Path,
        help="CheckM tabular file (i.e. the output of `checkm qa --tab_table …`) used "
        "to define thresholds as a function of MAG quality (MAGs with higher "
        "contamination levels will undergo less conservative filtering). "
        "Overrides `--probability_threshold`.",
    )
    options.add_argument(
        "--probability_threshold",
        default=default_values["filter"]["probability_threshold"],
        type=float,
        help="Contigs whose estimated probability of being a contaminant is above "
        "`probability_threshold` will be filtered out. The higher this value, the more "
        "conservative will be the contig filtering.",
    )
    options.add_argument(
        "--fast_mode",
        help="Identify contaminants using only the composition and coverage modules. The "
        "codon_usage and taxonomy modules will not be executed and the `--taxonomy_database` "
        "parameter is not required.",
        action="store_true",
    )
    options.add_argument(
        "--suffix",
        default=default_values["filter"]["suffix"],
        type=str,
        help="Suffix to be added to the filenames of the filtered genomes, before "
        "the file extension (e.g.: 'genome.suffix.fna').",
    )
    other.add_argument(
        "-t",
        "--threads",
        default=multiprocessing.cpu_count(),
        type=int,
        help="Number of threads to use. All by default.",
    )
    other.add_argument(
        "-q", "--quiet", help="Suppress the logger output", action="store_true",
    )
    other.add_argument(
        "-h", "--help", action="help", help="Show this help message and exit",
    )


def composition_parser(parser):
    parser.set_defaults(func=magpurify2.composition.main)
    required = parser.add_argument_group("Required arguments")
    options = parser.add_argument_group("Data processing options")
    other = parser.add_argument_group("Other options")
    required.add_argument(
        "genomes",
        nargs="+",
        help="Input genomes in the FASTA format. '.bz2', '.gz' and '.xz' compressed "
        "files are also supported.",
        type=Path,
    )
    required.add_argument(
        "output_directory",
        help="Directory where the logs, temporary and supporting data, and "
        "final score files will be written to.",
        type=Path,
    )
    options.add_argument(
        "--n_iterations",
        default=default_values["composition"]["n_iterations"],
        help="Number of iterations of data embedding and cluster selection with "
        "different seeds.",
        type=int,
    )
    options.add_argument(
        "--n_components",
        default=default_values["composition"]["n_components"],
        help="The dimension of the space to embed the data into (UMAP parameter).",
        type=int,
    )
    options.add_argument(
        "--min_dist",
        default=default_values["composition"]["min_dist"],
        help="The effective minimum distance between embedded points (UMAP parameter).",
        type=float,
    )
    options.add_argument(
        "--n_neighbors",
        default=default_values["composition"]["n_neighbors"],
        help="The size of local neighborhood used for manifold approximation "
        "(UMAP parameter).",
        type=int,
    )
    options.add_argument(
        "--set_op_mix_ratio",
        default=default_values["composition"]["set_op_mix_ratio"],
        help="Interpolate between the union (1.0) and intersection (0.0) to combine the "
        "local simplicial sets (UMAP parameter).",
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
        "-q", "--quiet", help="Suppress the logger output", action="store_true",
    )
    other.add_argument(
        "-h", "--help", action="help", help="Show this help message and exit",
    )


def coverage_parser(parser):
    parser.set_defaults(func=magpurify2.coverage.main)
    required = parser.add_argument_group("Required arguments")
    coverage_input = parser.add_argument_group(
        "Coverage data input (required, mutually exclusive arguments)",
    )
    coverage_input_mex = coverage_input.add_mutually_exclusive_group(required=True)
    options = parser.add_argument_group("Data processing options")
    other = parser.add_argument_group("Other options")
    required.add_argument(
        "genomes",
        nargs="+",
        help="Input genomes in the FASTA format. '.bz2', '.gz' and '.xz' compressed "
        "files are also supported.",
        type=Path,
    )
    required.add_argument(
        "output_directory",
        help="Directory where the logs, temporary and supporting data, and "
        "final score files will be written to.",
        type=Path,
    )
    coverage_input_mex.add_argument(
        "--bam_files", nargs="+", help="Input sorted BAM files.", type=Path,
    )
    coverage_input_mex.add_argument(
        "--coverage_file",
        help="Input a tabular file containing contig coverage data.",
        type=Path,
    )
    options.add_argument(
        "--min_identity",
        default=default_values["coverage"]["min_identity"],
        help="Exclude reads by overall identity to the reference sequences.",
        type=float,
    )
    options.add_argument(
        "--trim_lower",
        default=default_values["coverage"]["trim_lower"],
        help="Fraction to trim from the lower tail of the coverage distribution to "
        "compute contig mean coverages.",
        type=float,
    )
    options.add_argument(
        "--trim_upper",
        default=default_values["coverage"]["trim_upper"],
        help="Fraction to trim from the upper tail of the coverage distribution to "
        "compute contig mean coverages.",
        type=float,
    )
    options.add_argument(
        "--min_average_coverage",
        default=default_values["coverage"]["min_average_coverage"],
        help="Do not compute scores using coverage data from samples where the average "
        "genome coverage is less than `min_average_coverage`.",
        type=float,
    )
    options.add_argument(
        "--n_iterations",
        default=default_values["coverage"]["n_iterations"],
        help="Number of iterations of data embedding and cluster selection with "
        "different seeds.",
        type=int,
    )
    options.add_argument(
        "--n_components",
        default=default_values["coverage"]["n_components"],
        help="The dimension of the space to embed the data into (UMAP parameter).",
        type=int,
    )
    options.add_argument(
        "--min_dist",
        default=default_values["coverage"]["min_dist"],
        help="The effective minimum distance between embedded points (UMAP parameter).",
        type=float,
    )
    options.add_argument(
        "--n_neighbors",
        default=default_values["coverage"]["n_neighbors"],
        help="The size of local neighborhood used for manifold approximation "
        "(UMAP parameter).",
        type=int,
    )
    options.add_argument(
        "--set_op_mix_ratio",
        default=default_values["coverage"]["set_op_mix_ratio"],
        help="Interpolate between the union (1.0) and intersection (0.0) to combine the "
        "local simplicial sets (UMAP parameter).",
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
        "-q", "--quiet", help="Suppress the logger output", action="store_true",
    )
    other.add_argument(
        "-h", "--help", action="help", help="Show this help message and exit",
    )


def codon_usage_parser(parser):
    parser.set_defaults(func=magpurify2.codon_usage.main)
    required = parser.add_argument_group("Required arguments")
    options = parser.add_argument_group("Data processing options")
    other = parser.add_argument_group("Other options")
    required.add_argument(
        "genomes",
        nargs="+",
        help="Input genomes in the FASTA format. '.bz2', '.gz' and '.xz' compressed "
        "files are also supported.",
        type=Path,
    )
    required.add_argument(
        "output_directory",
        help="Directory where the logs, temporary and supporting data, and "
        "final score files will be written to.",
        type=Path,
    )
    options.add_argument(
        "--n_iterations",
        default=default_values["codon_usage"]["n_iterations"],
        help="Number of iterations of data embedding and cluster selection with "
        "different seeds.",
        type=int,
    )
    options.add_argument(
        "--n_components",
        default=default_values["codon_usage"]["n_components"],
        help="The dimension of the space to embed the data into (UMAP parameter).",
        type=int,
    )
    options.add_argument(
        "--min_dist",
        default=default_values["codon_usage"]["min_dist"],
        help="The effective minimum distance between embedded points (UMAP parameter).",
        type=float,
    )
    options.add_argument(
        "--n_neighbors",
        default=default_values["codon_usage"]["n_neighbors"],
        help="The size of local neighborhood used for manifold approximation "
        "(UMAP parameter).",
        type=int,
    )
    options.add_argument(
        "--set_op_mix_ratio",
        default=default_values["codon_usage"]["set_op_mix_ratio"],
        help="Interpolate between the union (1.0) and intersection (0.0) to combine the "
        "local simplicial sets (UMAP parameter).",
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
        "-q", "--quiet", help="Suppress the logger output", action="store_true",
    )
    other.add_argument(
        "-h", "--help", action="help", help="Show this help message and exit",
    )


def taxonomy_parser(parser):
    parser.set_defaults(func=magpurify2.taxonomy.main)
    required = parser.add_argument_group("Required arguments")
    options = parser.add_argument_group("Data processing options")
    other = parser.add_argument_group("Other options")
    required.add_argument(
        "genomes",
        nargs="+",
        help="Input genomes in the FASTA format. '.bz2', '.gz' and '.xz' compressed "
        "files are also supported.",
        type=Path,
    )
    required.add_argument(
        "output_directory",
        help="Directory where the logs, temporary and supporting data, and "
        "final score files will be written to.",
        type=Path,
    )
    required.add_argument(
        "taxonomy_database",
        help="Path to the MAGpurify2 taxonomy database directory.",
        type=Path,
    )
    options.add_argument(
        "--contig_min_fraction",
        default=default_values["taxonomy"]["contig_min_fraction"],
        type=float,
        help="The contig-level taxonomy must agree with the taxonomy of at least "
        "`contig_min_fraction` of its genes (weighted by the bitscore of the alignment). "
        "The value must be equal to or greater than 0.5 and less than 1.",
    )
    options.add_argument(
        "--genome_min_fraction",
        default=default_values["taxonomy"]["genome_min_fraction"],
        type=float,
        help="The genome-level taxonomy must agree with the taxonomy of at least "
        "`genome_min_fraction` of its genes (weighted by the bitscore of the alignment). "
        "The value must be equal to or greater than 0.5 and less than 1.",
    )
    options.add_argument(
        "--min_genus_identity",
        default=default_values["taxonomy"]["min_genus_identity"],
        type=float,
        help="A given gene will only be assigned up to the genus rank if the value "
        "obtained from 'alignment identity * target coverage' is equal to or greater "
        "than `min_genus_identity`. In case this condition is not satisfied, the gene "
        "will be assigned up to the family rank.",
    )
    other.add_argument(
        "-t",
        "--threads",
        default=multiprocessing.cpu_count(),
        type=int,
        help="Number of threads to use. All by default.",
    )
    other.add_argument(
        "-q", "--quiet", help="Suppress the logger output", action="store_true",
    )
    other.add_argument(
        "-h", "--help", action="help", help="Show this help message and exit",
    )


def filter_parser(parser):
    parser.set_defaults(func=magpurify2.filter.main)
    required = parser.add_argument_group("Required arguments")
    options = parser.add_argument_group("Data processing options")
    other = parser.add_argument_group("Other options")
    required.add_argument(
        "genomes",
        nargs="+",
        help="Input genomes in the FASTA format. '.bz2', '.gz' and '.xz' compressed "
        "files are also supported.",
        type=Path,
    )
    required.add_argument(
        "output_directory",
        help="Directory where the logs, temporary and supporting data, and "
        "final score files will be written to.",
        type=Path,
    )
    required.add_argument(
        "filtered_output_directory",
        help="Directory where the filtered MAGs will be written to.",
        type=Path,
    )
    options.add_argument(
        "--checkm_file",
        type=Path,
        help="CheckM tabular file (i.e. the output of `checkm qa --tab_table …`) used "
        "to define thresholds as a function of MAG quality (MAGs with higher "
        "contamination levels will undergo less conservative filtering). "
        "Overrides `--probability_threshold`.",
    )
    options.add_argument(
        "--probability_threshold",
        default=default_values["filter"]["probability_threshold"],
        type=float,
        help="Contigs whose estimated probability of being a contaminant is above "
        "`probability_threshold` will be filtered out. The higher this value, the more "
        "conservative will be the contig filtering.",
    )
    options.add_argument(
        "--fast_mode",
        help="Use only the outputs of the composition and coverage modules. The "
        "outputs of the codon_usage and taxonomy modules are not required.",
        action="store_true",
    )
    options.add_argument(
        "--suffix",
        default=default_values["filter"]["suffix"],
        type=str,
        help="Suffix to be added to the filenames of the filtered genomes, before "
        "the file extension (e.g.: 'genome.suffix.fna').",
    )
    other.add_argument(
        "-t",
        "--threads",
        default=multiprocessing.cpu_count(),
        type=int,
        help="Number of threads to use. All by default.",
    )
    other.add_argument(
        "-q", "--quiet", help="Suppress the logger output", action="store_true",
    )
    other.add_argument(
        "-h", "--help", action="help", help="Show this help message and exit",
    )
