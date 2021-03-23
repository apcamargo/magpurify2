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

import imp
import logging
import sys
from pathlib import Path

import numpy as np

from magpurify2 import tools
from magpurify2.core import ContigClassifier, Mag


def main(args):
    logger = logging.getLogger("timestamp")
    args.contig_min_fraction = tools.validade_input(
        args.probability_threshold, "probability_threshold", [0.0, 1.0]
    )
    if not args.output_directory.is_dir():
        logger.error(
            f"The directory '{args.output_directory}' was not found. Please sure you have provided "
            "the correct path and that MAGpurify2's modules have been executed."
        )
        sys.exit(1)
    scores_directory = args.output_directory.joinpath("scores")
    if not scores_directory.is_dir():
        logger.error(
            f"The directory '{scores_directory}' was not found. Please sure you have executed "
            "MAGpurify2's modules."
        )
        sys.exit(1)
    if not args.filtered_output_directory.is_dir():
        logger.info(
            f"'{args.filtered_output_directory}' does not exist. Creating it now."
        )
        args.filtered_output_directory.mkdir()
    contig_prob_file = scores_directory.joinpath("contamination_probabilities.tsv")

    logger.info(f"Reading {len(args.genomes)} genomes.")
    mag_list = [Mag(genome) for genome in args.genomes]

    # Determine the required modules and check if their outputs exist
    package_path = imp.find_module("magpurify2")[1]
    if args.fast_mode:
        module_list = ["composition", "coverage"]
        model_file = Path(package_path).joinpath("models", "fast_model.json")
    else:
        module_list = ["composition", "coverage", "codon_usage", "taxonomy"]
        model_file = Path(package_path).joinpath("models", "full_model.json")
    logger.info(f"Reading contig scores from the {', '.join(module_list)} modules.")
    module_path_list = [
        scores_directory.joinpath(f"{module}_scores.tsv") for module in module_list
    ]
    missing_modules = [
        str(module_path)
        for module_path in module_path_list
        if not module_path.exists()
    ]

    if missing_modules:
        logger.error(
            f"The following files were not found: '{', '.join(missing_modules)}'."
        )
        sys.exit(1)

    # Load data from the modules outputs
    genome_contig_matrix = np.genfromtxt(
        module_path_list[0],
        dtype=str,
        skip_header=1,
        comments=None,
        delimiter="\t",
        usecols=[0, 1],
    )
    feature_matrix = []
    module_ncol_dict = {"composition": 5, "coverage": 5, "codon_usage": 7, "taxonomy": 6}
    for module, module_path in zip(module_list, module_path_list):
        if not (
            np.genfromtxt(
                module_path,
                dtype=str,
                skip_header=1,
                comments=None,
                delimiter="\t",
                usecols=[0, 1],
            )
            == genome_contig_matrix
        ).all():
            logger.error(
                "It seems that MAGpurify2's modules were executed with different inputs. "
                "Please make sure to execute all the modules with the exact same inputs in the "
                "same order."
            )
            sys.exit(1)
        module_ncol = module_ncol_dict[module]
        module_feature_matrix = np.genfromtxt(
            module_path,
            dtype=float,
            skip_header=1,
            comments=None,
            delimiter="\t",
            usecols=range(2, module_ncol),
        )
        feature_matrix.append(module_feature_matrix)
    feature_matrix = np.column_stack(feature_matrix)
    # Transform 'contig_length'
    feature_matrix[:, 2] = np.tanh(np.log10(feature_matrix[:, 2] + 1) * 0.2)
    # Transform 'n_samples'
    feature_matrix[:, 5] = np.clip((feature_matrix[:, 5] / 10), 0, 1)
    if not args.fast_mode:
        # Transform 'n_genes'
        feature_matrix[:, 8] = np.tanh(np.log10(feature_matrix[:, 8] + 1))
        # Transform 'total_cds_length'
        feature_matrix[:, 9] = np.tanh(np.log10(feature_matrix[:, 9] + 1) * 0.2)
        # Transform 'genome_rank'
        feature_matrix[:, 12] = feature_matrix[:, 12] / 6
        # Transform 'contig_rank'
        feature_matrix[:, 13] = feature_matrix[:, 13] / 6

    logger.info(
        f"Estimating contamination probabilities and identifying contaminant contigs."
    )
    contig_classification = ContigClassifier(
        genome_contig_matrix,
        feature_matrix,
        model_file,
        args.fast_mode,
        args.probability_threshold,
        args.checkm_file,
        args.threads,
    )

    logger.info(f"Writing contig contamination probabilities to: '{contig_prob_file}'.")
    tools.write_module_output([contig_classification], contig_prob_file)

    logger.info(f"Writing filtered genomes to '{args.filtered_output_directory}'.")
    input_genomes = {mag.genome for mag in mag_list}
    not_in_mags_contaminants = input_genomes.difference(
        contig_classification.mags_contaminants_dict.keys()
    )
    if not_in_mags_contaminants:
        logger.warning(
            "The following genomes were not previously analysed by any of MAGpurify2's "
            f"modules and will be skipped: {', '.join(not_in_mags_contaminants)}."
        )
    for mag in mag_list:
        tools.write_filtered_genome(
            mag,
            contig_classification.mags_contaminants_dict,
            args.filtered_output_directory,
            args.suffix,
        )
