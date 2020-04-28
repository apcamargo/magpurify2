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
from collections import defaultdict

import taxopy
from joblib import Parallel, delayed

from magpurify2 import external, tools
from magpurify2.core import Composition, Coverage, Mag, Taxonomy


def composition_module(args):
    logger = logging.getLogger("timestamp")
    args.strictness = tools.check_strictness(args.strictness, logger)
    tools.check_output_directory(args.output_directory, logger)
    logger.info(f"Reading {len(args.genomes)} genomes.")
    mag_list = [Mag(genome) for genome in args.genomes]
    logger.info("Computing tetranucleotide frequencies.")
    composition_dict = {
        mag.genome: tools.get_tnf(mag.sequences, threads=args.threads) for mag in mag_list
    }
    logger.info("Identifying putative contaminants.")
    mag_composition_list = Parallel(n_jobs=args.threads)(
        delayed(Composition)(mag, composition_dict, args.strictness, args.threads)
        for mag in mag_list
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


def coverage_module(args):
    logger = logging.getLogger("timestamp")
    args.strictness = tools.check_strictness(args.strictness, logger)
    tools.check_bam_files(args.bam_files, logger)
    tools.check_output_directory(args.output_directory, logger)
    logger.info(f"Reading {len(args.genomes)} genomes.")
    mag_list = [Mag(genome, store_sequences=False) for genome in args.genomes]
    logger.info(f"Computing contig coverages from {len(args.bam_files)} BAM files.")
    coverage_dict = tools.get_coverages(
        [str(filepath) for filepath in args.bam_files], threads=args.threads
    )
    logger.info("Identifying putative contaminants.")
    mag_coverage_list = Parallel(n_jobs=args.threads)(
        delayed(Coverage)(mag, coverage_dict, args.strictness, args.threads)
        for mag in mag_list
    )
    for mag_coverage in mag_coverage_list:
        logger.info(
            f"{mag_coverage.genome}: {sum(mag_coverage.contaminants)}/"
            f"{len(mag_coverage)} contigs flagged as contaminants."
        )
    coverage_output_file = args.output_directory.joinpath("contaminants_coverage.tsv")
    logger.info(f"Writing output to: '{coverage_output_file}'.")
    tools.write_contamination_output(mag_coverage_list, coverage_output_file)


def taxonomy_module(args):
    logger = logging.getLogger("timestamp")
    args.strictness = tools.check_strictness(args.strictness, logger)
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
    taxonomy_dict = tools.get_taxonomy_dict(mmseqs2_output_file, taxdb, 0.75)
    logger.info("Identifying putative contaminants.")
    mag_taxonomy_list = Parallel(n_jobs=args.threads)(
        delayed(Taxonomy)(mag, taxonomy_dict, args.strictness, taxdb) for mag in mag_list
    )
    for mag_taxonomy in mag_taxonomy_list:
        logger.info(
            f"{mag_taxonomy.genome}: {sum(mag_taxonomy.contaminants)}/"
            f"{len(mag_taxonomy)} contigs flagged as contaminants."
        )
    taxonomy_output_file = args.output_directory.joinpath("contaminants_taxonomy.tsv")
    logger.info(f"Writing output to: '{taxonomy_output_file}'.")
    tools.write_contamination_output(mag_taxonomy_list, taxonomy_output_file)


def filter_module(args):
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
