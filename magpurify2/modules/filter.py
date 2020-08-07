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
from collections import defaultdict

from joblib import Parallel, delayed

from magpurify2 import tools
from magpurify2.core import Mag


def main(args):
    logger = logging.getLogger("timestamp")
    logger.info(f"Reading {len(args.genomes)} genomes.")
    mag_list = [Mag(genome) for genome in args.genomes]
    input_genomes = {mag.genome for mag in mag_list}
    scores_directory = args.output_directory.joinpath("scores")
    scores_file_list = scores_directory.glob("*_scores.tsv")
    mags_contaminants = defaultdict(lambda: defaultdict(list))
    module_thresholds = {
        "codon_usage": args.codon_usage,
        "composition": args.composition_threshold,
        "coverage": args.coverage_threshold,
        "coverage_clust": args.coverage_clustering_threshold,
        "taxonomy": args.taxonomy_threshold,
    }
    if not args.filtered_output_directory.is_dir():
        logger.warning(
            f"'{args.filtered_output_directory}' does not exist. Creating it now."
        )
        args.filtered_output_directory.mkdir()
    for score_file in scores_file_list:
        logger.info(f"Reading '{score_file}'.")
        current_module = score_file.stem.replace("_scores", "")
        with open(score_file) as fin:
            for line in fin:
                genome, contig, score = line.split()
                if genome in input_genomes:
                    if float(score) < module_thresholds[current_module]:
                        mags_contaminants[genome][contig].append(False)
                    else:
                        mags_contaminants[genome][contig].append(True)
    logger.info("Writing filtered genomes.")
    not_in_mags_contaminants = set(input_genomes).difference(mags_contaminants.keys())
    if not_in_mags_contaminants:
        logger.warning(
            "The following genomes were not previously analysed by any of MAGpurify2's "
            f"modules and will be skipped: {', '.join(not_in_mags_contaminants)}."
        )
    Parallel(n_jobs=args.threads)(
        delayed(tools.write_filtered_genome)(
            mag, mags_contaminants, args.mode, args.filtered_output_directory
        )
        for mag in mag_list
    )
