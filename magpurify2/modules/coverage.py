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

import gzip
import logging
import pickle
from functools import reduce

import numpy as np
from joblib import Parallel, delayed

from magpurify2 import tools
from magpurify2.core import Coverage, Mag


def main(args):
    logger = logging.getLogger("timestamp")
    logger.info("Executing MAGpurify2 coverage module.")
    args.contig_min_fraction = tools.validade_input(
        args.min_identity, "min_identity", [0.0, 1.0]
    )
    args.trim_lower = tools.validade_input(args.trim_lower, "trim_lower", [0.0, 1.0])
    args.trim_upper = tools.validade_input(args.trim_upper, "trim_upper", [0.0, 1.0])
    args.min_average_coverage = tools.validade_input(
        args.min_average_coverage, "min_average_coverage", [0, 999]
    )
    args.n_iterations = tools.validade_input(args.n_iterations, "n_iterations", [1, 999])
    args.n_components = tools.validade_input(args.n_components, "n_components", [1, 999])
    args.n_neighbors = tools.validade_input(args.n_neighbors, "n_neighbors", [1, 999])
    args.set_op_mix_ratio = tools.validade_input(
        args.set_op_mix_ratio, "set_op_mix_ratio", [0.0, 1.0]
    )
    tools.check_output_directory(args.output_directory)
    scores_directory = args.output_directory.joinpath("scores")
    scores_directory.mkdir(exist_ok=True)
    coverage_score_file = scores_directory.joinpath("coverage_scores.tsv")

    # Check if coverages were provided as BAMs or a tabular file
    if args.bam_files:
        tools.check_bam_files(args.bam_files)
        input_coverage = args.bam_files
    else:
        input_coverage = [args.coverage_file]

    # Read input genomes
    logger.info(f"Reading {len(args.genomes)} genomes.")
    mag_list = [Mag(genome, store_sequences=False) for genome in args.genomes]
    coverage_data_file = args.output_directory.joinpath("coverages.pickle.gz")
    signatures = [tools.get_file_signature(filepath) for filepath in input_coverage]

    # Check if the coverage computation needs to be executed again. To do that, we check
    # if a pickle dump exists. If it does, we check if the signatures of all input BAM
    # files are in it.
    skip_coverage = False
    if coverage_data_file.exists():
        with gzip.open(coverage_data_file) as fin:
            loaded_signatures, contig_names, coverage_matrix = pickle.load(fin)
            if set(signatures).issubset(loaded_signatures):
                logger.info("Skipping contig coverages computation.")
                skip_coverage = True
                # If the number of input BAM files is less than the number of input BAM
                # files, select the coverage data of the input BAMs.
                if len(loaded_signatures) > len(signatures):
                    selected_input_files = [
                        signature in signatures for signature in loaded_signatures
                    ]
                    coverage_matrix = coverage_matrix[:, selected_input_files]
            else:
                coverage_data_file.unlink()
    if not skip_coverage:
        if args.bam_files:
            logger.info(
                f"Computing contig coverages from {len(input_coverage)} BAM files."
            )
            contig_names, coverage_matrix = tools.get_bam_coverages(
                [str(filepath) for filepath in input_coverage],
                min_identity=args.min_identity,
                trim_lower=args.trim_lower,
                trim_upper=args.trim_upper,
                contig_set=reduce(lambda x, y: x.union(y.contigs), mag_list, set()),
                threads=args.threads,
            )
            contig_names = np.array(contig_names)
        elif args.coverage_file:
            logger.info(f"Reading contig coverages from '{input_coverage[0]}'.")
            contig_names, coverage_matrix = tools.get_tsv_coverages(
                input_coverage[0],
                contig_set=reduce(lambda x, y: x.union(y.contigs), mag_list, set()),
            )
        logger.info(f"Saving contig coverages data to '{coverage_data_file}'.")
        with gzip.open(coverage_data_file, "wb") as fout:
            pickle.dump((signatures, contig_names, coverage_matrix), fout)

    # Build a dictionary where the keys are genome names and the values are numpy matrices
    # of the coverage values
    coverage_dict = Parallel(n_jobs=args.threads)(
        delayed(lambda x, y, z: z[[np.where(y == i)[0][0] for i in x.contigs]])(
            mag, contig_names, coverage_matrix,
        )
        for mag in mag_list
    )
    coverage_dict = dict(zip([mag.genome for mag in mag_list], coverage_dict))

    logger.info("Computing contig scores.")

    mag_coverage_list = Parallel(n_jobs=args.threads)(
        delayed(Coverage)(
            mag,
            coverage_dict[mag.genome],
            args.min_average_coverage,
            args.n_iterations,
            args.n_components,
            args.min_dist,
            args.n_neighbors,
            args.set_op_mix_ratio,
        )
        for mag in mag_list
    )
    logger.info(f"Writing output to: '{coverage_score_file}'.")
    tools.write_module_output(mag_coverage_list, coverage_score_file)
