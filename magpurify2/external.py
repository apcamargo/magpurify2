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
import gzip
import logging
import lzma
import shutil
import subprocess
import sys

from biolib.external.prodigal import Prodigal

from magpurify2 import tools

logger = logging.getLogger("timestamp")


class Database:
    def __init__(self, path):
        if not path.is_dir():
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
            found_files = {filepath.name for filepath in path.iterdir()}
            if expected_files.issubset(found_files):
                with open(path.joinpath("magpurify2DB_version")) as fin:
                    self.version = fin.read().strip()
                with open(path.joinpath("magpurify2DB_name")) as fin:
                    self.name = fin.read().strip()
                self.mmseqs_db = path.joinpath(self.name)
                self.nodes_dmp = path.joinpath("magpurify2DB_nodes.dmp")
                self.names_dmp = path.joinpath("magpurify2DB_names.dmp")
            else:
                missing_files = expected_files - found_files
                logger.error(
                    "The following files are missing from the database: "
                    f"{', '.join(missing_files)}."
                )
                sys.exit(1)


def prodigal(genomes, output_directory, threads):
    prodigal_output_directory = output_directory.joinpath("prodigal")
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
        extracted_genomes_directory = output_directory.joinpath("extracted_genomes")
        if n_compressed_genomes:
            # If there is any compressed input to Prodigal, a new directory will be
            # created and compressed files will be decompressed into it. Non-compressed
            # inputs will simply be copied into this new directory.
            logger.info(f"Extracting {n_compressed_genomes} compressed genomes")
            if extracted_genomes_directory.is_dir():
                shutil.rmtree(extracted_genomes_directory)
            extracted_genomes_directory.mkdir()
            for filepath in genomes_without_prediction:
                genome = filepath.stem
                if tools.is_compressed(filepath) != tools.Compression.noncompressed:
                    genome = genome.rsplit(".", 1)[0]
                    filepath_compression = tools.is_compressed(filepath)
                    extraction_output = extracted_genomes_directory.joinpath(
                        f"{genome}.fna"
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
                    copy_output = extracted_genomes_directory.joinpath(f"{genome}.fna")
                    shutil.copyfile(filepath, copy_output)
            # Change the path to Prodigal inputs to the new directory containing the
            # decompressed FASTA files.
            genomes_without_prediction = list(extracted_genomes_directory.glob("*.fna"))
        # Execute Prodigal for the genomes without predictions.
        logger.info(
            f"Predicting proteins from {len(genomes_without_prediction)} genomes using Prodigal."
        )
        prodigal_runner = Prodigal(threads, verbose=False)
        # Convert Path objects to strings in order to use biolib's Prodigal wrapper.
        genomes_without_prediction = [
            str(filepath) for filepath in genomes_without_prediction
        ]
        prodigal_runner.run(genomes_without_prediction, prodigal_output_directory)
        # If it exists, delete the directory containing decompressed FASTA files.
        if extracted_genomes_directory.is_dir():
            shutil.rmtree(extracted_genomes_directory)
    # Delete the "*.gff" files generated by Prodigal.
    for filepath in prodigal_output_directory.glob("*.gff"):
        filepath.unlink()


def mmseqs2(output_directory, database, threads):
    log_file = output_directory.joinpath("mmseqs2.log")
    mmseqs2_output_directory = output_directory.joinpath("mmseqs2")
    mmseqs2_input_file = mmseqs2_output_directory.joinpath("mmseqs2_input.faa")
    mmseqs2_taxonomy_file = mmseqs2_output_directory.joinpath("mmseqs2_taxonomy.tsv")
    mmseqs2_alignment_file = mmseqs2_output_directory.joinpath("mmseqs2_alignment.tsv")
    querydb_directory = mmseqs2_output_directory.joinpath("querydb")
    querydb_prefix = querydb_directory.joinpath("querydb")
    taxonomydb_directory = mmseqs2_output_directory.joinpath("taxonomydb")
    taxonomydb_prefix = taxonomydb_directory.joinpath("taxonomydb")
    taxonomydb_aln_prefix = taxonomydb_directory.joinpath("taxonomydb_aln")
    besthitdb_directory = mmseqs2_output_directory.joinpath("besthitdb")
    besthitdb_prefix = besthitdb_directory.joinpath("besthitdb")
    tmp_directory = mmseqs2_output_directory.joinpath("tmp")
    # Create intermediate directories.
    querydb_directory.mkdir()
    taxonomydb_directory.mkdir()
    besthitdb_directory.mkdir()
    tmp_directory.mkdir()
    # Define the MMseqs2 commands:
    createdb_command = ["mmseqs", "createdb", mmseqs2_input_file, querydb_prefix]
    taxonomy_command = [
        "mmseqs",
        "taxonomy",
        querydb_prefix,
        database.mmseqs_db,
        taxonomydb_prefix,
        tmp_directory,
        "-s",
        "3.0",
        "--blacklist",
        "''",
        "--lca-mode",
        "3",
        "--tax-output-mode",
        "2",
        "--threads",
        str(threads),
    ]
    createtsv_command = [
        "mmseqs",
        "createtsv",
        querydb_prefix,
        taxonomydb_prefix,
        mmseqs2_taxonomy_file,
    ]
    best_hit_command = [
        "mmseqs",
        "filterdb",
        taxonomydb_aln_prefix,
        besthitdb_prefix,
        "--extract-lines",
        "1",
        "--threads",
        str(threads),
    ]
    convertalis_command = [
        "mmseqs",
        "convertalis",
        querydb_prefix,
        database.mmseqs_db,
        besthitdb_prefix,
        mmseqs2_alignment_file,
        "--format-output",
        "query,fident,tcov,bits",
        "--threads",
        str(threads),
    ]
    with open(log_file, "w") as fout:
        for command in [
            createdb_command,
            taxonomy_command,
            createtsv_command,
            best_hit_command,
            convertalis_command,
        ]:
            try:
                subprocess.run(command, stdout=fout, stderr=fout, check=True)
            except subprocess.CalledProcessError:
                command_str = " ".join([str(i) for i in command])
                logger.error(f"'{command_str}' failed.")
                sys.exit(1)
    # Remove intermediate directories.
    shutil.rmtree(querydb_directory)
    shutil.rmtree(taxonomydb_directory)
    shutil.rmtree(besthitdb_directory)
    shutil.rmtree(tmp_directory)
