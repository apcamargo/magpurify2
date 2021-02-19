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
import csv
import gzip
import hashlib
import logging
import lzma
import os
import sys
import textwrap
from collections import defaultdict
from enum import Enum, auto
from pathlib import Path

import hdbscan
import numpy as np
import umap
from Bio import SeqIO, bgzf

from magpurify2._codon import get_codon_usage_profile, get_cai, get_codon_index
from magpurify2._composition import get_gc, get_tnf
from magpurify2._coverage import get_bam_coverages
from magpurify2._mmseqs2 import get_mmseqs2

logger = logging.getLogger("timestamp")


class Compression(Enum):
    bzip2 = auto()
    gzip = auto()
    xz = auto()
    noncompressed = auto()


def validade_input(value, parameter_name, interval):
    """
    Checks the value of a numeric input. If it is between the provided interval,
    the same value is returned. If the value is outside of the interval the
    minimum or the maximum value will be returned (if the inpur is less than the
    minimum value or greater than the maximum value, respectively).

    Parameters
    ----------
    value : float
        The numerical value that will be validated.
    parameter_name : str
        Name of the input.
    interval : list
        A list containing the minimum and maximum allowerd values.

    Returns
    -------
    float
        A value thats falls within the provided interval.
    """
    interval = sorted(interval)
    if value > interval[1] or value < interval[0]:
        value = min(max(value, interval[0]), interval[1])
        logger.warning(
            f"The value of the `{parameter_name}` parameter must be between {interval[0]} and "
            f"{interval[1]}. Setting it to {value}."
        )
    return value


def check_output_directory(output_directory):
    """
    Checks if the directory exists. If it doesn't, the logger will give a
    warning and the directory will be created.

    Parameters
    ----------
    output_directory : Path
        Path object pointing to the output directory of the program.
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
        `bzip2`, `gzip`, `xz`, or `noncompressed`
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


def check_bam_files(bam_files):
    """
    Checks if the input BAM files exist and if they are sorted. Exits with an
    error otherwise.

    Parameters
    ----------
    bam_files : list
        List with the Path objects pointing to the input BAM files.
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


def get_tsv_coverages(filepath, contig_set=None):
    """
    Reads contig coverages from tabular text files. The first column of the
    coverage file must store the contigs names. The remaining columns should
    contain the coverage of each contig across multiple samples.

    Parameters
    ----------
    filepath : Path
        Path object pointing tabular text file.
    contig_set : set, optional
       If provided, only the coverages of the contigs within `contig_set` will
       returned. Default is None (return the coverages of all contigs).

    Returns
    -------
    tuple
        A tuple whose fist element is an numpy array of the contig names and the
        second one is a numpy matrix of contig coverages in the input BAM files.
    """
    with open(filepath) as fin:
        ncols = len(fin.readline().split("\t"))
    if ncols < 2:
        logger.error("The tabular coverage file must have at least two columns.")
        sys.exit(1)
    contig_names_vector = np.genfromtxt(
        filepath, dtype=str, comments=None, delimiter="\t", usecols=0
    )
    coverage_vector = np.genfromtxt(
        filepath, dtype=float, comments=None, delimiter="\t", usecols=range(1, ncols)
    )
    if coverage_vector.ndim == 1:
        coverage_vector = coverage_vector.reshape(-1, 1)
    if contig_set:
        mask = np.array([contig in contig_set for contig in contig_names_vector])
        contig_names_vector = contig_names_vector[mask]
        coverage_vector = coverage_vector[mask, :]
    return (contig_names_vector, coverage_vector)


def get_weighted_median(data, weights):
    """
    Computes the weighted median of the input data.

    Parameters
    ----------
    data : array-like
        Data vector whose weighted median will be computed.
    weights : array-like
        An array of weights associated with the values in the input data.

    Returns
    -------
    ndarray
        Weighted median of the input data.
    """
    data, weights = np.array(data).squeeze(), np.array(weights).squeeze()
    s_data, s_weights = map(np.array, zip(*sorted(zip(data, weights))))
    midpoint = 0.5 * sum(s_weights)
    if np.any(weights > midpoint):
        w_median = data[weights == np.max(weights)][0]
    else:
        cs_weights = np.cumsum(s_weights)
        idx = np.where(cs_weights <= midpoint)[0][-1]
        if cs_weights[idx] == midpoint:
            w_median = np.mean(s_data[idx : idx + 2])
        else:
            w_median = s_data[idx + 1]
    return w_median


def get_log_ratio_scores(data, weights, base, pseudo=1e-5, min=0.0, max=1.0):
    """
    Computes the deviation from the weighted median as the log of the ratio
    between each data point and the weighted median.

    Parameters
    ----------
    data : array-like
        Input data that will be used to compute the score.
    weights : array-like
        An array of weights associated with the values in the input data. Used
        to compute the weighted median.
    base : float
        Base of the log.
    pseudo : float
        Pseudocount added to the numerator and denominator in the ratio between
        the data and the weighted median.
    min : float
        Minimum allowed value within the resulting score array.
    max : float
        Maximum allowed value within the resulting score array.

    Returns
    -------
    ndarray
        Array of containing the score of each data point.
    """
    if data.ndim == 1:
        weighted_median = get_weighted_median(data, weights=weights)
    else:
        weighted_median = np.apply_along_axis(
            get_weighted_median, 0, data, weights=weights
        )
    scores = (data + pseudo) / (weighted_median + pseudo)
    scores = 1 - np.abs(np.log(scores) / np.log(base))
    if data.ndim != 1:
        scores = np.average(scores, axis=1)
    return np.clip(scores, min, max)


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
        if is_compressed(filepath) != Compression.noncompressed:
            genome_name = genome_name.rsplit(".", 1)[0]
        genome_faa = prodigal_output_directory.joinpath(f"{genome_name}_genes.faa")
        if not genome_faa.exists():
            genomes_without_prediction.append(filepath)
        else:
            genomes_with_prediction.append(filepath)
    return genomes_without_prediction, genomes_with_prediction


def write_mmseqs2_input(output_directory):
    """
    Prepare the FASTA input for MMseqs2 by the concatenation of Prodigal
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


def create_embedding(
    data,
    n_components,
    min_dist,
    n_neighbors,
    set_op_mix_ratio,
    metric="euclidean",
    random_state=42,
):
    """
    Creates an UMAP embedding from the input data.

    Parameters
    ----------
    data : array-like
        Input data to be fit into an embedded space.
    n_components : int
        The dimension of the space to embed into.
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
    init = (
        "random"
        if data.shape[1] <= n_components or data.shape[0] <= n_components + 1
        else "spectral"
    )
    reducer = umap.UMAP(
        n_components=n_components,
        metric=metric,
        min_dist=min_dist,
        n_neighbors=n_neighbors,
        set_op_mix_ratio=set_op_mix_ratio,
        init=init,
        random_state=random_state,
    )
    return reducer.fit_transform(data)


def get_cluster_score(data, allow_single_cluster, lengths):
    """
    Clusters the data with HDBSCAN, identify the main cluster as the one that
    contains the largest sum of contig lengths and returns the score of each
    contig as the measured membership to the main cluster.

    Parameters
    ----------
    data : array-like
        Input data to be clustered.
    allow_single_cluster : bool
        Allow HDBSCAN to identify a single cluster.
    lengths : list
        Lengths of the contigs.

    Returns
    -------
    ndarray
        The score of each contig, measured as the membership of each contig to
        the main cluster.
    """
    # Get min_cluster_size and min_samples values based on the size of the dataset.
    if len(lengths) >= 250:
        min_cluster_size = 5
        min_samples = len(lengths) // 50
    elif len(lengths) >= 150:
        min_cluster_size = 4
        min_samples = 4
    elif len(lengths) >= 50:
        min_cluster_size = 3
        min_samples = 3
    elif len(lengths) >= 6:
        min_cluster_size = 2
        min_samples = 2
    else:
        min_cluster_size = 2
        min_samples = 1
    weights = defaultdict(int)
    clusterer = hdbscan.HDBSCAN(
        min_cluster_size=min_cluster_size,
        min_samples=min_samples,
        prediction_data=True,
        allow_single_cluster=allow_single_cluster,
        core_dist_n_jobs=1,
    ).fit(data)
    clusters = clusterer.labels_
    try:
        soft_clusters = hdbscan.all_points_membership_vectors(clusterer)
        # If no clusters were identified, give a score of 1 to all contigs.
        if np.all(clusters == -1):
            scores = np.ones(len(lengths))
        else:
            for cluster, lenght in zip(clusters, lengths):
                weights[cluster] += lenght
            # Ignore unclustered contigs to identify the main cluster.
            weights.pop(-1, None)
            selected_cluster = max(weights, key=weights.get)
            scores = soft_clusters[:, selected_cluster]
        return scores
    except KeyError:
        return np.ones(len(lengths))


def get_cluster_score_from_embedding(
    data, lengths, n_iterations, n_components, min_dist, n_neighbors, set_op_mix_ratio
):
    """
    Iterativelly create UMAP embeddings and compute contigs score using
    different seeds.

    Parameters
    ----------
    data : array-like
        Input data to be fit into an embedded space.
    lengths : list
        Lengths of the contigs.
    n_iterations : int
        Number of iterations to execute using different seeds.
    n_components : int
        The dimension of the space to embed into.
    min_dist : float
        The effective minimum distance between embedded points.
    n_neighbors : int
        The size of local neighborhood used for manifold approximation.
    set_op_mix_ratio : float
        Interpolate between (fuzzy) union and intersection as the set operation
        used to combine local fuzzy simplicial sets to obtain a global fuzzy
        simplicial sets.

    Returns
    -------
    ndarray
        Embedding of the training data in low-dimensional space.
    """
    scores = np.zeros(len(lengths))
    n_neighbors = len(lengths) - 1 if len(lengths) <= n_neighbors else n_neighbors
    for i in range(n_iterations):
        embedding = create_embedding(
            data=data,
            n_components=n_components,
            min_dist=min_dist,
            n_neighbors=n_neighbors,
            set_op_mix_ratio=set_op_mix_ratio,
            random_state=i,
        )
        scores += get_cluster_score(
            data=embedding, allow_single_cluster=True, lengths=lengths,
        )
    scores = scores / max(scores)
    return scores


def write_contig_taxonomy_output(mag_taxonomy_list, taxonomy_output_file):
    """
    Write a file containing the genome and contig-level taxonomic lineages for
    each contig.

    Parameters
    ----------
    mag_taxonomy_list : list
        List containing the `Taxonomy` objects of the input genomes.
    taxonomy_output_file : Path
        Path object pointing to the file where the taxonomic lineages will be
        written to.
    """
    with open(taxonomy_output_file, "w") as fout:
        fout.write("genome\tcontig\tgenome_taxonomy\tcontig_taxonomy\n")
        for mag_taxonomy in mag_taxonomy_list:
            genome_taxonomy_str = ";".join(
                reversed(mag_taxonomy.genome_taxonomy.name_lineage)
            )
            for index, (genome, contig, _, _, _, _) in enumerate(mag_taxonomy):
                contig_taxonomy_str = ";".join(
                    reversed(mag_taxonomy.contig_taxonomy[index].name_lineage)
                )
                fout.write(
                    f"{genome}\t{contig}\t"
                    f"{genome_taxonomy_str}\t{contig_taxonomy_str}\n"
                )


def write_module_output(score_list, score_output_file):
    """
    Write a file containing the computed attributes for each contig of the input
    genomes.

    Parameters
    ----------
    score_list : list
        List of objects of the `CodonUsage`, `Composition`, `Coverage`,
        `Taxonomy` or `ContigClassifier` classes.
    score_output_file : Path
        Path object pointing to the output file.
    """
    with open(score_output_file, "w") as fout:
        header = "\t".join(score_list[0].attributes)
        fout.write(f"{header}\n")
        for mag_score in score_list:
            for attributes in mag_score:
                line = "\t".join(map(str, attributes))
                fout.write(f"{line}\n")


def get_checkm_contamination(filepath):
    """
    Reads a CheckM tabular file (i.e. the output of `checkm qa --tab_table …`)
    and returns a dictionary where the keys are genome names and the values are
    the estimated contamination values.

    Parameters
    ----------
    filepath : Path
        Path object pointing to a CheckM tabular file.

    Returns
    -------
    dict
        A dictionary where the keys are genome names and the values are the
        estimated contamination values.
    """
    if filepath.exists():
        with open(filepath) as fin:
            csvreader = csv.reader(fin, delimiter="\t")
            fields = next(csvreader)
            if set(["Bin Id", "Completeness", "Contamination"]).issubset(fields):
                bin_id_col = fields.index("Bin Id")
                completeness_col = fields.index("Completeness")
                contamination_col = fields.index("Contamination")
                checkm_cont_dict = {}
                for row in csvreader:
                    genome = row[bin_id_col]
                    checkm_cont_dict[genome] = float(row[contamination_col])
                return checkm_cont_dict
            else:
                logger.warning(
                    "The supplied CheckM tabular file is not correctly formatted. The "
                    "default probability threshold value will be used instead."
                )
                return None
    else:
        logger.warning(
            "The supplied CheckM tabular file was not found. The default probability "
            "threshold value will be used instead."
        )
        return None


def write_filtered_genome(mag, mags_contaminants_dict, filtered_output_directory, suffix):
    """
    Writes a FASTA file containing contigs not flagged as contaminants.

    Parameters
    ----------
    mag : Mag
        A `Mag` object of a genome.
    mags_contaminants_dict : dict
        A dictionary where keys are the names of the contigsand the values are
        booleans that indicate if the contig was classified as a contaminant.
    filtered_output_directory : Path
        Path object pointing to the directory where the filtered genome will be
        written to.
    suffix : str
        Suffix to be added to the filenames of the filtered genomes.
    """
    if mag.genome in mags_contaminants_dict:
        if suffix:
            output_fna = filtered_output_directory.joinpath(f"{mag.genome}.{suffix}.fna")
        else:
            output_fna = filtered_output_directory.joinpath(f"{mag.genome}.fna")
        with open(output_fna, "w") as fout:
            for contig, description, sequence in mag:
                if not mags_contaminants_dict[mag.genome][contig]:
                    fout.write(f">{description}\n")
                    fout.write(f"{textwrap.fill(sequence, 70)}\n")
