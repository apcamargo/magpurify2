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
    args.strictness = tools.check_strictness(args.strictness, logger)
    tools.check_output_directory(args.output_directory, logger)
    logger.info(f"Reading {len(args.genomes)} genomes.")
    mag_list = [Mag(genome) for genome in args.genomes]
    composition_data_file = args.output_directory.joinpath(
        "tnfs.pickle.gz"
    )

    # Check if the TNF computation needs to be executed again. To do that, we check if
    # a pickle dump exists. If it does, we check if all the input genomes are in it.
    skip_tnf = False
    if composition_data_file.exists():
        with gzip.open(composition_data_file) as fin:
            composition_dict = pickle.load(fin)
            if {mag.genome for mag in mag_list}.issubset(composition_dict.keys()):
                logger.info("Skipping tetranucleotide frequencies computation.")
                skip_tnf = True
            else:
                composition_data_file.unlink()
    if not skip_tnf:
        logger.info("Computing tetranucleotide frequencies.")
        composition_dict = {
            mag.genome: tools.get_tnf(mag.sequences, threads=args.threads)
            for mag in mag_list
        }
        logger.info(
            f"Saving tetranucleotide frequencies data to '{composition_data_file}'."
        )
        with gzip.open(composition_data_file, "wb") as fout:
            pickle.dump(composition_dict, fout)

    logger.info("Identifying putative contaminants.")
    mag_composition_list = Parallel(n_jobs=args.threads)(
        delayed(Composition)(mag, composition_dict, args.strictness) for mag in mag_list
    )
    for mag_composition in mag_composition_list:
        logger.info(
            f"{mag_composition.genome}: {sum(mag_composition.contaminants)}/"
            f"{len(mag_composition)} contigs flagged as contaminants."
        )
    composition_output_file = args.output_directory.joinpath(
        "contaminants_composition.tsv"
    )
    logger.info(f"Writing output to: '{composition_output_file}'.")
    tools.write_contamination_output(mag_composition_list, composition_output_file)