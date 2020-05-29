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
from magpurify2.core import Mag, Taxonomy


def main(args):
    logger = logging.getLogger("timestamp")
    args.contig_min_fraction = tools.validade_input(
        args.contig_min_fraction, "contig_min_fraction", [0.5, 1.0], logger
    )
    args.genome_min_fraction = tools.validade_input(
        args.genome_min_fraction, "genome_min_fraction", [0.5, 1.0], logger
    )
    # Check if Prodigal and MMSeqs2 are executables in the user PATH.
    missing_executables = [
        executable
        for executable in ["prodigal", "mmseqs"]
        if not tools.check_executable(executable)
    ]
    if missing_executables:
        logger.error(
            "The following dependencies are missing and must be installed: "
            f"{', '.join(missing_executables)}."
        )
        sys.exit(1)
    database = external.Database(args.database, logger)
    logger.info(f"Using database '{database.name}' version {database.version}.")
    tools.check_output_directory(args.output_directory, logger)
    logger.info(f"Reading {len(args.genomes)} genomes.")
    mag_list = [Mag(genome, store_sequences=False) for genome in args.genomes]

    mmseqs2_output_directory = args.output_directory.joinpath("mmseqs2")
    mmseqs2_input_file = mmseqs2_output_directory.joinpath("mmseqs2_input.faa")
    mmseqs2_output_file = mmseqs2_output_directory.joinpath("mmseqs2_output.tsv")

    # Check if MMSeqs2 needs to be executed again. To do that, we check if a MMSeqs2
    # output file exists. If it does, we check if all the input genomes are in it.
    skip_mmseqs = False
    if mmseqs2_output_directory.is_dir():
        if mmseqs2_output_file.exists():
            with open(mmseqs2_output_file) as fin:
                searched_genomes = {line.split("~")[0] for line in fin}
            if {mag.genome for mag in mag_list}.issubset(searched_genomes):
                logger.info("Skipping MMSeqs2 search.")
                skip_mmseqs = True
            else:
                shutil.rmtree(mmseqs2_output_directory)
                mmseqs2_output_directory.mkdir()
        else:
            shutil.rmtree(mmseqs2_output_directory)
            mmseqs2_output_directory.mkdir()
    else:
        mmseqs2_output_directory.mkdir()

    if not skip_mmseqs:
        external.prodigal(args.genomes, args.output_directory, logger, args.threads)
        logger.info(f"Writing MMSeqs2 input file to: '{mmseqs2_input_file}'.")
        tools.write_mmseqs2_input(args.output_directory)
        logger.info("Running MMSeqs2 for gene taxonomic assignment.")
        external.mmseqs2(args.output_directory, database, logger, args.threads)

    # Create a TaxDb object containing GTDB's taxonomic data.
    taxdb = taxopy.TaxDb(
        nodes_dmp=database.nodes_dmp, names_dmp=database.names_dmp, keep_files=True,
    )
    logger.info(f"Reading MMSeqs2 output file.")
    taxonomy_dict = tools.get_taxonomy_dict(mmseqs2_output_file, taxdb)
    logger.info("Computing contig scores.")
    mag_taxonomy_list = [
        Taxonomy(
            mag,
            taxonomy_dict,
            args.contig_min_fraction,
            args.genome_min_fraction,
            args.allow_genus,
            taxdb,
        )
        for mag in mag_list
    ]
    taxonomy_score_file = args.output_directory.joinpath("scores_taxonomy.tsv")
    contig_taxonomy_file = args.output_directory.joinpath("taxonomic_assignment.tsv")
    logger.info(
        f"Writing output to: '{taxonomy_score_file} and '{contig_taxonomy_file}'."
    )
    tools.write_contig_score_output(mag_taxonomy_list, taxonomy_score_file)
    tools.write_contig_taxonomy_output(mag_taxonomy_list, contig_taxonomy_file)
