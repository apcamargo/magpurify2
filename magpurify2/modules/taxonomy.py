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
    logger.info("Executing MAGpurify2 taxonomy module.")
    args.contig_min_fraction = tools.validade_input(
        args.contig_min_fraction, "contig_min_fraction", [0.5, 1.0]
    )
    args.genome_min_fraction = tools.validade_input(
        args.genome_min_fraction, "genome_min_fraction", [0.01, 0.99]
    )
    args.min_genus_identity = tools.validade_input(
        args.min_genus_identity, "min_genus_identity", [0.0, 1.0]
    )
    # Check if Prodigal and MMseqs2 are executables in the user PATH.
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
    database = external.Database(args.taxonomy_database)
    logger.info(f"Using database '{database.name}' version {database.version}.")
    tools.check_output_directory(args.output_directory)
    scores_directory = args.output_directory.joinpath("scores")
    scores_directory.mkdir(exist_ok=True)
    taxonomy_score_file = scores_directory.joinpath("taxonomy_scores.tsv")
    contig_taxonomy_file = args.output_directory.joinpath("taxonomic_assignment.tsv")

    # Read input genomes
    logger.info(f"Reading {len(args.genomes)} genomes.")
    mag_list = [Mag(genome, store_sequences=False) for genome in args.genomes]

    mmseqs2_output_directory = args.output_directory.joinpath("mmseqs2")
    mmseqs2_input_file = mmseqs2_output_directory.joinpath("mmseqs2_input.faa")
    mmseqs2_taxonomy_file = mmseqs2_output_directory.joinpath("mmseqs2_taxonomy.tsv")
    mmseqs2_alignment_file = mmseqs2_output_directory.joinpath("mmseqs2_alignment.tsv")

    # Check if MMseqs2 needs to be executed again. To do that, we check if a MMseqs2
    # output file exists. If it does, we check if all the input genomes are in it.
    skip_mmseqs = False
    if mmseqs2_output_directory.is_dir():
        if mmseqs2_taxonomy_file.exists() and mmseqs2_alignment_file.exists():
            with open(mmseqs2_taxonomy_file) as fin:
                searched_genomes_taxonomy = {line.split("~")[0] for line in fin}
            with open(mmseqs2_alignment_file) as fin:
                searched_genomes_alignment = {line.split("~")[0] for line in fin}
            if {mag.genome for mag in mag_list}.issubset(searched_genomes_taxonomy):
                logger.info("Skipping MMseqs2 search.")
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
        external.prodigal(args.genomes, args.output_directory, args.threads)
        logger.info(f"Writing MMseqs2 input file to: '{mmseqs2_input_file}'.")
        tools.write_mmseqs2_input(args.output_directory)
        logger.info("Running MMseqs2 for gene taxonomic assignment.")
        external.mmseqs2(args.output_directory, database, args.threads)

    # Create a TaxDb object containing GTDB's taxonomic data.
    taxdb = taxopy.TaxDb(
        nodes_dmp=database.nodes_dmp, names_dmp=database.names_dmp, keep_files=True,
    )
    logger.info(f"Reading MMseqs2 output files.")
    mmseqs2_dict = tools.get_mmseqs2(str(mmseqs2_taxonomy_file), str(mmseqs2_alignment_file))
    logger.info("Computing contig scores.")
    mag_taxonomy_list = [
        Taxonomy(
            mag,
            mmseqs2_dict.get(mag.genome, None),
            args.contig_min_fraction,
            args.genome_min_fraction,
            args.min_genus_identity,
            taxdb,
        )
        for mag in mag_list
    ]
    # Write contig score file
    logger.info(
        f"Writing output to: '{taxonomy_score_file}' and '{contig_taxonomy_file}'."
    )
    tools.write_module_output(mag_taxonomy_list, taxonomy_score_file)
    tools.write_contig_taxonomy_output(mag_taxonomy_list, contig_taxonomy_file)
