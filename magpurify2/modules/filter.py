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

import logging
import sys
from collections import defaultdict

import numpy as np
from joblib import Parallel, delayed

from magpurify2 import tools
from magpurify2.core import ContigClassifier, Mag


def main(args):
    logger = logging.getLogger("timestamp")
    args.contig_min_fraction = tools.validade_input(
        args.probability_threshold, "probability_threshold", [0.0, 1.0]
    )
    if not args.output_directory.exists():
        logger.error(
            f"The directory '{args.output_directory}' was not found. Please sure you have provided "
            "the correct path and that MAGpurify2's modules have been executed."
        )
        sys.exit(1)
    scores_directory = args.output_directory.joinpath("scores")
    if not scores_directory.exists():
        logger.error(
            f"The directory '{scores_directory}' was not found. Please sure you have executed "
            "MAGpurify2's modules."
        )
        sys.exit(1)
    contig_prob_file = scores_directory.joinpath("contamination_probabilities.tsv")

    logger.info(f"Reading {len(args.genomes)} genomes.")
    mag_list = [Mag(genome) for genome in args.genomes]

    # Determine the required modules and check if their outputs exist
    if args.fast_mode:
        module_list = ["composition", "coverage"]
    else:
        module_list = ["composition", "coverage", "codon_usage", "taxonomy"]
    logger.info(f"Reading contig scores from the {', '.join(module_list)} modules.")
    module_path_list = [
        scores_directory.joinpath(f"{module}_scores.tsv") for module in module_list
    ]
    missing_modules = []
    for module_path in module_path_list:
        if not module_path.exists():
            missing_modules.append(str(module_path))
    if missing_modules:
        logger.error(f"The following files were not found: '{', '.join(missing_modules)}'.")
        sys.exit(1)

    # Load data from the modules outputs
    genome_contig_matrix = np.loadtxt(
        module_path_list[0], dtype=str, skiprows=1, usecols=[0, 1]
    )
    feature_matrix = []
    module_ncol_dict = {"composition": 5, "coverage": 5, "codon_usage": 6, "taxonomy": 3}
    for module, module_path in zip(module_list, module_path_list):
        if (
            np.loadtxt(module_path, dtype=str, skiprows=1, usecols=[0, 1])
            != genome_contig_matrix
        ).all():
            logger.error(
                "It seems that MAGpurify2's modules were executed with different inputs. "
                "Please make sure to execute all the modules with the exact same inputs in the "
                "same order."
            )
            sys.exit(1)
        module_ncol = module_ncol_dict[module]
        module_feature_matrix = np.loadtxt(
            module_path, dtype=float, skiprows=1, usecols=range(2, module_ncol)
        )
        feature_matrix.append(module_feature_matrix)
    feature_matrix = np.column_stack(feature_matrix)
    # Transform 'contig_length'
    feature_matrix[:, 2] = np.tanh(np.log10(feature_matrix[:, 2] + 1) * 0.2)
    # Transform 'n_samples'
    feature_matrix[:, 5] = np.clip((feature_matrix[:, 5] / 10), 0, 1)
    if not args.fast_mode:
        # Transform 'n_genes'
        feature_matrix[:, 7] = np.tanh(np.log10(feature_matrix.loc[:, "n_genes"] + 1))
        # Transform 'total_cds_length'
        feature_matrix[:, 8] = np.tanh(np.log10(feature_matrix[:, 8] + 1) * 0.2)

    logger.info(
        f"Estimating contamination probabilities and identifying contaminant contigs."
    )
    contig_classification = ContigClassifier(
        genome_contig_matrix,
        feature_matrix,
        args.model_file,
        args.probability_threshold,
        args.threads,
    )

    logger.info(f"Writing contig contamination probabilities to: '{contig_prob_file}'.")
    tools.write_module_output([contig_classification], contig_prob_file)

    logger.info(f"Writing filtered genomes to '{args.filtered_output_directory}'.")
    input_genomes = {mag.genome for mag in mag_list}
    not_in_mags_contaminants = input_genomes.difference(mags_contaminants.keys())
    if not_in_mags_contaminants:
        logger.warning(
            "The following genomes were not previously analysed by any of MAGpurify2's "
            f"modules and will be skipped: {', '.join(not_in_mags_contaminants)}."
        )
    Parallel(n_jobs=args.threads)(
        delayed(tools.write_filtered_genome)(
            mag, contig_classification, args.filtered_output_directory
        )
        for mag in mag_list
    )
