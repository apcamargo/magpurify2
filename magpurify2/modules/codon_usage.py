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
import shutil
import sys

import taxopy

from magpurify2 import external, tools
from magpurify2.core import Mag, CodonUsage


def main(args):
    logger = logging.getLogger("timestamp")
    args.min_genes = tools.validade_input(args.min_genes, "min_genes", [1, 999], logger)
    # Check if Prodigal is an executable in the user PATH.
    missing_executables = [
        executable
        for executable in ["prodigal"]
        if not tools.check_executable(executable)
    ]
    if missing_executables:
        logger.error(
            "The following dependencies are missing and must be installed: "
            f"{', '.join(missing_executables)}."
        )
        sys.exit(1)
    tools.check_output_directory(args.output_directory, logger)
    scores_directory = args.output_directory.joinpath("scores")
    scores_directory.mkdir(exist_ok=True)
    codon_usage_score_file = scores_directory.joinpath("codon_usage_scores.tsv")

    # Read input genomes
    logger.info(f"Reading {len(args.genomes)} genomes.")
    mag_list = [Mag(genome, store_sequences=False) for genome in args.genomes]

    external.prodigal(args.genomes, args.output_directory, logger, args.threads)
    prodigal_output_directory = args.output_directory.joinpath("prodigal")

    logger.info("Computing contig scores.")
    mag_codon_usage_list = [
        CodonUsage(
            mag,
            args.min_genes,
            prodigal_output_directory.joinpath(mag.genome + "_genes.fna"),
        )
        for mag in mag_list
    ]
    # Write contig score file
    logger.info(f"Writing output to: '{codon_usage_score_file}'.")
    tools.write_contig_score_output(mag_codon_usage_list, codon_usage_score_file)
