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
from itertools import compress

from joblib import Parallel, delayed

from magpurify2 import tools
from magpurify2.core import Coverage, Mag


def main(args):
    logger = logging.getLogger("timestamp")
    args.contig_min_fraction = tools.validade_input(
        args.min_identity, "min_identity", [0.0, 1.0], logger
    )
    tools.check_bam_files(args.bam_files, logger)
    tools.check_output_directory(args.output_directory, logger)
    logger.info(f"Reading {len(args.genomes)} genomes.")
    mag_list = [Mag(genome, store_sequences=False) for genome in args.genomes]
    coverage_data_file = args.output_directory.joinpath("coverages.pickle.gz")
    bam_signatures = [tools.get_file_signature(filepath) for filepath in args.bam_files]

    # Check if the coverage computation needs to be executed again. To do that, we check
    # if a pickle dump exists. If it does, we check if the signatures of all input BAM
    # files are in it.
    skip_coverage = False
    if coverage_data_file.exists():
        with gzip.open(coverage_data_file) as fin:
            loaded_bam_signatures, coverage_dict = pickle.load(fin)
            if set(bam_signatures).issubset(loaded_bam_signatures):
                logger.info("Skipping contig coverages computation.")
                skip_coverage = True
                # If the number of input BAM files is less than the number of input BAM
                # files, select the coverage data of the input BAMs.
                if loaded_bam_signatures > bam_signatures:
                    selected_bams = [
                        signature in bam_signatures for signature in loaded_bam_signatures
                    ]
                    coverage_dict = {
                        contig: list(compress(coverages, selected_bams))
                        for contig, coverages in coverage_dict.items()
                    }
            else:
                coverage_data_file.unlink()
    if not skip_coverage:
        logger.info(f"Computing contig coverages from {len(args.bam_files)} BAM files.")
        coverage_dict = tools.get_coverages(
            [str(filepath) for filepath in args.bam_files],
            min_identity=args.min_identity,
            threads=args.threads,
        )
        logger.info(f"Saving contig coverages data to '{coverage_data_file}'.")
        with gzip.open(coverage_data_file, "wb") as fout:
            pickle.dump((bam_signatures, coverage_dict), fout)

    logger.info("Computing contig scores.")
    mag_coverage_list = Parallel(n_jobs=args.threads)(
        delayed(Coverage)(mag, coverage_dict) for mag in mag_list
    )
    coverage_score_file = args.output_directory.joinpath("scores_coverage.tsv")
    logger.info(f"Writing output to: '{coverage_score_file}'.")
    tools.write_contig_score_output(mag_coverage_list, coverage_score_file)
