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
import hashlib
import lzma
import os
import sys
import textwrap
from collections import defaultdict
from enum import Enum, auto
from pathlib import Path

import taxopy
import umap
from Bio import SeqIO, bgzf
from sklearn.cluster import MeanShift

from magpurify2._coverage import get_coverages
from magpurify2._tnf import get_tnf


class Compression(Enum):
    gzip = auto()
    bzip2 = auto()
    xz = auto()
    noncompressed = auto()


def check_strictness(strictness, logger):
    """
    Checks the value of the `strictness` argument. If it is between 0 and 1, the
    same value is returned. If the value is less than 0 or greater than 1, the
    logger will give a warning and the returned value will be 0 (if the input
    was less than 0) or 1 (is the input was greater than 1)

    Parameters
    ----------
    strictness : float
        The value of the `strictness` argument.
    logger : Logger
        The logger of the program.

    Returns
    -------
    float
        A float between 0 and 1.
    """
    if strictness > 1 or strictness < 0:
        strictness = min(max(strictness, 0.0), 1.0)
        logger.warning(
            "The value of the `strictness` parameter must be between 0 and 1. "
            f"Setting it to {strictness}."
        )
    return strictness


def check_output_directory(output_directory, logger):
    """
    Checks if the directory exists. If it doesn't, the logger will give a
    warning and the directory will be created.

    Parameters
    ----------
    output_directory : Path
        Path object pointing to the output directory of the program.
    logger : Logger
        The logger of the program.
    """
    if not output_directory.is_dir():
        logger.warning("Output directory does not exist. Creating it now.")
        output_directory.mkdir()


def check_executable(executable):
    """
    Checks if a executable is available in the PATH.

    Parameters
    ----------
    executable : str
        Name of the executable.

    Returns
    -------
    bool
        Returns `True` if the executable is in the PATH and `False` if it isn't.
    """
    found = False
    for path in os.environ["PATH"].split(os.pathsep):
        executable_path = os.path.join(path.strip('"'), executable)
        if os.path.isfile(executable_path) and os.access(executable_path, os.X_OK):
            found = True
            break
    return found


def is_compressed(filepath):
    """
    Checks if a file is compressed (gzip, bzip2 or xz).

    Parameters
    ----------
    filepath : Path
        Path object pointing to a file.

    Returns
    -------
    Compression
        Returns a `Compression` enum with one of the following attributes:
        `gzip`, `bzip2`, `xz`, or `noncompressed`
    """
    with open(filepath, "rb") as fin:
        signature = fin.peek(8)[:8]
        if tuple(signature[:2]) == (0x1F, 0x8B):
            return Compression.gzip
        elif tuple(signature[:3]) == (0x42, 0x5A, 0x68):
            return Compression.bzip2
        elif tuple(signature[:7]) == (0xFD, 0x37, 0x7A, 0x58, 0x5A, 0x00, 0x00):
            return Compression.xz
        else:
            return Compression.noncompressed


def get_file_signature(filepath):
    """
    Returns the SHA-256 hash of the last 128 kilobytes of the input file.

    Parameters
    ----------
    filepath : Path
        Path object pointing to a file.

    Returns
    -------
    bytes
        SHA-256 hash of the last 128 kilobytes of the input file.
    """
    with open(filepath, "rb") as fin:
        fin.seek(-128000, os.SEEK_END)
        return hashlib.sha256(fin.read(128000)).hexdigest()


def is_sorted_bam(filepath):
    """
    Checks if a BAM file is sorted by coordinate.

    Parameters
    ----------
    filepath : Path
        Path object pointing to BAM file.

    Returns
    -------
    bool
        Returns `True` if the BAM file is sorted by coordinate and `False`
        otherwise.
    """
    with bgzf.BgzfReader(filepath, "rb") as fin:
        bam_header = fin.readline().strip()
        return b"SO:coordinate" in bam_header


def check_bam_files(bam_files, logger):
    """
    Checks if the input BAM files exist and if they are sorted. Exits with an
    error otherwise.

    Parameters
    ----------
    bam_files : list
        List with the Path objects pointing to the input BAM files.
    logger : Logger
        The logger of the program.
    """
    if not all(filepath.exists() for filepath in bam_files):
        logger.error("At least one of the supplied BAM files does not exist.")
        sys.exit(1)
    if not all(is_sorted_bam(filepath) for filepath in bam_files):
        logger.error(
            "At least one of the supplied BAM files is not sorted by coordinate "
            "or does not contain a valid header. You can sort it using `samtools sort`."
        )
        sys.exit(1)


def read_fasta(filepath):
    """
    Read FASTA file and yields the records names and sequences. Automatically
    detects if the file is compressed (gzip, bzip2 or xz) and opens it properly.

    Parameters
    ----------
    filepath : Path
        Path object pointing to a FASTA file.

    Yields
    -------
    str
        Name of a record.
    str
        Description of a record.
    str
        Sequence of a record.
    """
    filepath_compression = is_compressed(filepath)
    if filepath_compression == Compression.gzip:
        fin = gzip.open(filepath, "rt")
    elif filepath_compression == Compression.bzip2:
        fin = bz2.open(filepath, "rt")
    elif filepath_compression == Compression.xz:
        fin = lzma.open(filepath, "rt")
    else:
        fin = open(filepath, "r")
    for record in SeqIO.parse(fin, "fasta"):
        yield record.id.replace("~", "_"), record.description, str(record.seq)
    fin.close()


def check_prediction(genome_list, output_directory):
    """
    Checks which genomes have been processed by Prodigal and which ones still
    need to go through gene prediction.

    Parameters
    ----------
    genome_list : list
        List containing the Path objects pointing to the input genomes.
    output_directory : Path
        Path object pointing to the output directory of the program.

    Returns
    -------
    list
        List of the genomes for which gene prediction has not been performed.
    list
        List of the genomes for which gene prediction was already performed.
    """
    prodigal_output_directory = output_directory.joinpath("prodigal")
    genomes_with_prediction = []
    genomes_without_prediction = []
    for filepath in genome_list:
        genome_name = filepath.stem
        genome_faa = prodigal_output_directory.joinpath(f"{genome_name}_genes.faa")
        if not genome_faa.exists():
            genomes_without_prediction.append(filepath)
        else:
            genomes_with_prediction.append(filepath)
    return genomes_without_prediction, genomes_with_prediction


def write_mmseqs2_input(output_directory):
    """
    Prepare the FASTA input for MMSeqs2 by the concatenation of Prodigal
    outputs. The header of each record in the FASTA is formatted as
    ">genome~contig~gene" so that the results for each contig in each genome can
    be identified in downstream steps.

    Parameters
    ----------
    output_directory : Path
        Path object pointing to the output directory of the program.
    """
    mmseqs2_output_directory = output_directory.joinpath("mmseqs2")
    mmseqs2_input_file = mmseqs2_output_directory.joinpath("mmseqs2_input.faa")
    with open(mmseqs2_input_file, "w") as fout:
        prodigal_output_directory = output_directory.joinpath("prodigal")
        for filepath in prodigal_output_directory.glob("*.faa"):
            genome = Path(filepath).stem.replace("_genes", "")
            for name, _, sequence in read_fasta(filepath):
                contig, gene_number = name.rsplit("_", 1)
                contig.replace("~", "_")
                fout.write(f">{genome}~{contig}~{gene_number}\n")
                fout.write(f"{textwrap.fill(sequence, 70)}\n")


def get_taxonomy_dict(mmseqs2_output, taxdb, fraction=0.75):
    """
    Parse the MMSeqs2 output and assign a taxon to each contig based on the
    taxonomic assignment of its genes.

    Parameters
    ----------
    mmseqs2_output : Path
        Path object pointing to the MMSeqs2 output.
    taxdb : TaxDb
        A TaxDb object.
    fraction : float, default 0.75
        The taxon assigned to each contig will be shared by more than `fraction`
        of the contig genes. This value must be equal to or greater than 0.5 and
        less than 1.

    Returns
    -------
    dictionary
        Returns a nested dictionary. In the outer dictionary, the keys are the
        names of the genomes found in the MMSeqs2 output. In the inner
        dictionary the keys are the names of the contigs in the genome and the
        values are Taxon objects of the most specific taxon that is shared by
        more than the chosen fraction of the contig genes.
    """
    contig_taxids = defaultdict(lambda: defaultdict(list))
    with open(mmseqs2_output) as fin:
        for line in fin:
            taxid = line.split()[1]
            if taxid != "0":
                genome, contig, _ = line.split()[0].split("~")
                contig_taxids[genome][contig].append(taxid)
    taxonomy_dict = defaultdict(dict)
    for genome, contigs in contig_taxids.items():
        for contig in contigs:
            contig_taxon = [
                taxopy.Taxon(taxid, taxdb) for taxid in contig_taxids[genome][contig]
            ]
            if len(contig_taxon) > 1:
                majority_vote = taxopy.find_majority_vote(contig_taxon, taxdb, fraction)
            else:
                (majority_vote,) = contig_taxon
            taxonomy_dict[genome][contig] = majority_vote
    return taxonomy_dict


def create_embedding(
    data, min_dist, n_neighbors, set_op_mix_ratio, metric="euclidean", random_state=42,
):
    """
    Creates an UMAP embedding from the input data.

    Parameters
    ----------
    data : array-like
        Input data to be fit into an embedded space.
    min_dist : float
        The effective minimum distance between embedded points.
    n_neighbors : int
        The size of local neighborhood used for manifold approximation.
    set_op_mix_ratio : float
        Interpolate between (fuzzy) union and intersection as the set operation
        used to combine local fuzzy simplicial sets to obtain a global fuzzy
        simplicial sets.
    metric : str, default "euclidean"
        The metric to use to compute distances in high dimensional space.
    random_state : int, default 42
        The seed used by the random number generator.

    Returns
    -------
    ndarray
        Embedding of the training data in low-dimensional space.
    """
    os.environ["NUMBA_NUM_THREADS"] = "1"
    os.environ["THREADING_LAYER"] = "tbb"
    reducer = umap.UMAP(
        n_components=2,
        metric=metric,
        min_dist=min_dist,
        n_neighbors=n_neighbors,
        set_op_mix_ratio=set_op_mix_ratio,
        random_state=random_state,
    )
    return reducer.fit_transform(data)


def identify_contaminant_clusters(embedding, base_bandwidth, lengths, strictness):
    """
    Clusters the input data using the `MeanShift` algorithm and finds the main
    cluster using contigs lengths as weights. The reimaining clusters are
    labeled as contaminants.

    Parameters
    ----------
    embedding : array-like
        UMAP embedding of the contigs.
    base_bandwidth : float
        Base bandwidth used in the RBF kernel.
    lengths : list
        A list containing the lengths of the contigs in the input data.
    strictness : float
        Controls how strict the clustering is by adjusting the bandwidth with
        the following formula:
        `bandwidth = base_bandwidth + (strictness - 0.5) * 4`
        Accepts values ranging from 0 (more clusters) and 1 (less clusters).

    Returns
    -------
    ndarray
        A boolean array where `True` represents contaminants and `False`
        represents contigs contained in the main cluster.
    """
    bandwidth = base_bandwidth + (strictness - 0.5) * 4
    clusters = MeanShift(bandwidth=bandwidth, n_jobs=1).fit(embedding).labels_
    weights = defaultdict(int)
    for cluster, lenght in zip(clusters, lengths):
        weights[cluster] += lenght
    selected_cluster = max(weights, key=weights.get)
    return clusters != selected_cluster


def write_contamination_output(mag_contamination_list, contamination_output_file):
    """
    Write a file containing the information of which contigs of each genome were
    flagged as putative contaminants.

    Parameters
    ----------
    mag_contamination_list : list
        List of objects of the `Composition`, `Coverage` or `Taxonomy` classes.
    contamination_output_file : Path
        Path object pointing to the output file.
    """
    with open(contamination_output_file, "w") as fout:
        fout.write("genome\tcontig\tcontamination\n")
        for mag_contamination in mag_contamination_list:
            for contig, contaminant in mag_contamination:
                contaminant = "Yes" if contaminant else "No"
                fout.write(f"{mag_contamination.genome}\t{contig}\t{contaminant}\n")


def write_filtered_genome(mag, mags_contaminants, filtering, filtered_output_directory):
    """
    Writes a FASTA file containing contigs not flagged as contaminants.

    Parameters
    ----------
    mag : Mag
        A `Mag` object of a genome.
    mags_contaminants : dict
        A nested dictionary. In the outer dictionary, the keys are the names of
        genomes. In the inner dictionary the keys are the names of the contigs
        in the genome and the values are boolean lists. Each value in those
        lists indicates if the contig was flagged as a contaminant by one of the
        contaminant-detection approaches.
    filtering : str
        If `"any"`, the contig will be removed if it was flagged as a
        contaminant by any of the contaminant-detection approaches. If `"all"`,
        the contig will only be removed if it was flagged by all the approaches.
    filtered_output_directory : Path
        Path object pointing to the directory where the filtered genome will be
        written to.
    """
    output_fasta = filtered_output_directory.joinpath(mag.genome + ".filtered.fna")
    if mag.genome in mags_contaminants:
        with open(output_fasta, "w") as fout:
            for contig, description, sequence in mag:
                if filtering == "any":
                    if all(mags_contaminants[mag.genome][contig]):
                        fout.write(f">{description}\n")
                        fout.write(f"{textwrap.fill(sequence, 70)}\n")
                elif filtering == "all":
                    if any(mags_contaminants[mag.genome][contig]):
                        fout.write(f">{description}\n")
                        fout.write(f"{textwrap.fill(sequence, 70)}\n")
