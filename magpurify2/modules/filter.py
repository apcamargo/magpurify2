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
import shutil
import sys
from collections import defaultdict

import taxopy
from joblib import Parallel, delayed

from magpurify2 import external, tools
from magpurify2.core import Composition, Coverage, Mag, Taxonomy


def main(args):
    logger = logging.getLogger("timestamp")
    logger.info(f"Reading {len(args.genomes)} genomes.")
    mag_list = [Mag(genome) for genome in args.genomes]
    input_genomes = {mag.genome for mag in mag_list}
    contamination_file_list = args.output_directory.glob("contaminants*.tsv")
    mags_contaminants = defaultdict(lambda: defaultdict(list))
    if not args.filtered_output_directory.is_dir():
        logger.warning(
            "Output directory for the filtered genomes does not exist. Creating it now."
        )
        args.filtered_output_directory.mkdir()
    for contamination_file in contamination_file_list:
        logger.info(f"Reading '{contamination_file}'.")
        with open(contamination_file) as fin:
            for line in fin:
                genome, contig, contamination = line.split()
                if genome in input_genomes:
                    if contamination == "Yes":
                        mags_contaminants[genome][contig].append(False)
                    elif contamination == "No":
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
            mag, mags_contaminants, args.filtering, args.filtered_output_directory
        )
        for mag in mag_list
    )
