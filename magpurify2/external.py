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

import bz2
import glob
import gzip
import lzma
import os
import sys
import shutil
import subprocess

from biolib.external.prodigal import Prodigal

from magpurify2 import tools


class Database:
    def __init__(self, path, logger):
        if not os.path.isdir(path):
            logger.error("The provided database path does not point to a directory")
            sys.exit(1)
        else:
            expected_files = {
                "magpurify2DB",
                "magpurify2DB_h",
                "magpurify2DB.dbtype",
                "magpurify2DB.index",
                "magpurify2DB.lookup",
                "magpurify2DB.source",
                "magpurify2DB_h.dbtype",
                "magpurify2DB_h.index",
                "magpurify2DB_mapping",
                "magpurify2DB_merged.dmp",
                "magpurify2DB_names.dmp",
                "magpurify2DB_nodes.dmp",
                "magpurify2DB_name",
                "magpurify2DB_version",
            }
            found_files = set(os.listdir(path))
            if expected_files.issubset(found_files):
                with open(os.path.join(path, "magpurify2DB_version")) as fin:
                    self.version = fin.read().strip()
                with open(os.path.join(path, "magpurify2DB_name")) as fin:
                    self.name = fin.read().strip()
                self.mmseqs_db = os.path.join(path, self.name)
                self.nodes_dmp = os.path.join(path, "magpurify2DB_nodes.dmp")
                self.names_dmp = os.path.join(path, "magpurify2DB_names.dmp")
            else:
                missing_files = expected_files - found_files
                logger.error(
                    "The following files are missing from the database: "
                    f"{', '.join(missing_files)}."
                )
                sys.exit(1)


def prodigal(genomes, output_directory, logger, threads):
    prodigal_output_directory = os.path.join(output_directory, "prodigal")
    # Identify which of the input genomes already have Prodigal predictions in the
    # output directory. Those genomes won't be analyzed by Prodigal again.
    genomes_without_prediction, genomes_with_prediction = tools.check_prediction(
        genomes, output_directory
    )
    if genomes_with_prediction:
        logger.info(
            f"Skipping gene prediction for {len(genomes_with_prediction)} genomes."
        )
    if genomes_without_prediction:
        # Prodigal doesn't support compressed inputs, so we need to identify and
        # decompress them first.
        n_compressed_genomes = sum(
            tools.is_compressed(filepath) != tools.Compression.noncompressed
            for filepath in genomes_without_prediction
        )
        extracted_genomes_directory = os.path.join(output_directory, "extracted_genomes")
        if n_compressed_genomes:
            # If there is any compressed input to Prodigal, a new directory will be
            # created and compressed files will be decompressed into it. Non-compressed
            # inputs will simply be copied into this new directory.
            logger.info(f"Extracting {n_compressed_genomes} compressed genomes")
            if os.path.isdir(extracted_genomes_directory):
                shutil.rmtree(extracted_genomes_directory)
            os.mkdir(extracted_genomes_directory)
            for filepath in genomes_without_prediction:
                if tools.is_compressed(filepath) != tools.Compression.noncompressed:
                    filepath_compression = tools.is_compressed(filepath)
                    genome = tools.get_genome_name(filepath)
                    extraction_output = os.path.join(
                        extracted_genomes_directory, f"{genome}.fna"
                    )
                    if filepath_compression == tools.Compression.gzip:
                        fin = gzip.open(filepath, "rt")
                    elif filepath_compression == tools.Compression.bzip2:
                        fin = bz2.open(filepath, "rt")
                    elif filepath_compression == tools.Compression.xz:
                        fin = lzma.open(filepath, "rt")
                    with open(extraction_output, "w") as fout:
                        shutil.copyfileobj(fin, fout)
                    fin.close()
                else:
                    genome = tools.get_genome_name(filepath)
                    copy_output = os.path.join(
                        extracted_genomes_directory, f"{genome}.fna"
                    )
                    shutil.copyfile(filepath, copy_output)
            # Change the path to Prodigal inputs to the new directory containing the
            # decompressed FASTA files.
            genomes_without_prediction = glob.glob(
                os.path.join(extracted_genomes_directory, "*.fna")
            )
        # Execute Prodigal for the genomes without predictions.
        logger.info(
            f"Predicting proteins from {len(genomes_without_prediction)} genomes using Prodigal."
        )
        prodigal_runner = Prodigal(threads, verbose=False)
        prodigal_runner.run(genomes_without_prediction, prodigal_output_directory)
        # If it exists, delete the directory containing decompressed FASTA files.
        if os.path.isdir(extracted_genomes_directory):
            shutil.rmtree(extracted_genomes_directory)
    # Delete the ".fna" and ".gff" generated by Prodigal.
    fna_gff_files = glob.glob(
        os.path.join(prodigal_output_directory, "*.fna")
    ) + glob.glob(os.path.join(prodigal_output_directory, "*.gff"))
    for filepath in fna_gff_files:
        os.remove(filepath)


def mmseqs2(output_directory, database, logger, threads):
    log_file = os.path.join(output_directory, "mmseqs2.log")
    mmseqs2_output_directory = os.path.join(output_directory, "mmseqs2")
    mmseqs2_input_file = os.path.join(mmseqs2_output_directory, "mmseqs2_input.faa")
    mmseqs2_output_file = os.path.join(mmseqs2_output_directory, "mmseqs2_output.tsv")
    querydb_directory = os.path.join(mmseqs2_output_directory, "querydb")
    querydb_files = os.path.join(querydb_directory, "querydb")
    taxonomydb_directory = os.path.join(mmseqs2_output_directory, "taxonomydb")
    taxonomydb_files = os.path.join(taxonomydb_directory, "taxonomydb")
    tmp_directory = os.path.join(mmseqs2_output_directory, "tmp")
    # Create intermediate directories.
    os.mkdir(querydb_directory)
    os.mkdir(taxonomydb_directory)
    os.mkdir(tmp_directory)
    # Define the MMSeqs2 commands:
    first_command = ["mmseqs", "createdb", mmseqs2_input_file, querydb_files]
    second_command = [
        "mmseqs",
        "taxonomy",
        querydb_files,
        database.mmseqs_db,
        taxonomydb_files,
        tmp_directory,
        "-s",
        "3.5",
        "--lca-mode",
        "3",
        "--threads",
        str(threads),
    ]
    third_command = [
        "mmseqs",
        "createtsv",
        querydb_files,
        taxonomydb_files,
        mmseqs2_output_file,
    ]
    with open(log_file, "w") as fout:
        for command in [first_command, second_command, third_command]:
            try:
                subprocess.run(command, stdout=fout, stderr=fout, check=True)
            except subprocess.CalledProcessError:
                logger.error(f"'{' '.join(command)}' failed.")
                sys.exit(1)
    # Remove intermediate directories.
    shutil.rmtree(querydb_directory)
    shutil.rmtree(taxonomydb_directory)
    shutil.rmtree(tmp_directory)
