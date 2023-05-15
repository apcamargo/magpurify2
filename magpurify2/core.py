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
from collections import defaultdict

import numpy as np
import scipy.stats as ss
import taxopy
import xgboost as xgb
from scipy.signal import find_peaks

from magpurify2 import tools

logger = logging.getLogger("timestamp")


class Mag:
    def __init__(self, filepath, store_sequences=True):
        self.genome = filepath.stem
        if tools.is_compressed(filepath) != tools.Compression.noncompressed:
            self.genome = self.genome.rsplit(".", 1)[0]
        self.contigs = []
        self.descriptions = []
        self.sequences = []
        self.lengths = []
        for contig, description, sequence in tools.read_fasta(filepath):
            self.contigs.append(contig)
            self.descriptions.append(description)
            self.lengths.append(len(sequence))
            if store_sequences:
                self.sequences.append(sequence)
        self.lengths = np.array(self.lengths)

    def __len__(self):
        return len(self.contigs)

    def __iter__(self):
        if self.sequences:
            return zip(self.contigs, self.descriptions, self.sequences)
        else:
            return iter(self.contigs)

    def __repr__(self):
        return f"{self.genome}: {len(self)} contigs"


class CodonUsage:
    def __init__(
        self,
        mag,
        n_iterations,
        n_components,
        min_dist,
        n_neighbors,
        set_op_mix_ratio,
        prodigal_fna_filepath,
    ):
        self.attributes = [
            "genome",
            "contig",
            "codon_usage_cluster_score",
            "cai_score",
            "n_genes",
            "total_cds_length",
            "mean_strand_coding_density",
        ]
        self.genome = mag.genome
        self.contigs = mag.contigs
        self.lengths = mag.lengths
        self.cds = []
        self.cds_sequences = []
        self.n_genes = defaultdict(int)
        self.total_cds_length = defaultdict(int)
        for cds, _, cds_sequence in tools.read_fasta(prodigal_fna_filepath):
            self.cds.append(cds)
            self.cds_sequences.append(cds_sequence)
            contig, _ = cds.rsplit("_", 1)
            self.n_genes[contig] += 1
            self.total_cds_length[contig] += len(cds_sequence)
        self.n_genes = np.array([self.n_genes.get(contig, 0) for contig in self.contigs])
        self.total_cds_length = np.array(
            [self.total_cds_length.get(contig, 0) for contig in self.contigs]
        )
        self.mean_strand_coding_density = self.total_cds_length / (2 * self.lengths)
        if len(self) <= 2:
            self.codon_usage_cluster_scores = np.ones(len(self))
            self.cai_scores = np.ones(len(self))
        else:
            self.codon_usage_cluster_scores = self.compute_codon_usage_profile_scores(
                n_iterations, n_components, min_dist, n_neighbors, set_op_mix_ratio
            )
            self.delta_cai = self.get_delta_cai()
            self.cai_scores = self.compute_cai_scores()

    def compute_codon_usage_profile_scores(
        self, n_iterations, n_components, min_dist, n_neighbors, set_op_mix_ratio
    ):
        cds_codon_usage_profile = tools.get_codon_usage_profile(self.cds_sequences)
        cds_contig_name_list, cds_length_list = [], []
        for cds_name, cds_seq in zip(self.cds, self.cds_sequences):
            contig_name = cds_name.rsplit("_", 1)[0]
            cds_contig_name_list.append(contig_name)
            cds_length_list.append(len(cds_seq))
        cds_contig_name_list = np.array(cds_contig_name_list)
        cds_length_list = np.array(cds_length_list)
        contig_codon_usage_profile, contig_length_list, contig_name_list = [], [], []
        for contig_name, contig_length in zip(self.contigs, self.lengths):
            mask = cds_contig_name_list == contig_name
            if mask.sum():
                profile = np.average(
                    cds_codon_usage_profile[mask], weights=cds_length_list[mask], axis=0
                )
                contig_codon_usage_profile.append(profile)
                contig_length_list.append(contig_length)
                contig_name_list.append(contig_name)
        contig_codon_usage_profile = np.stack(contig_codon_usage_profile)
        codon_usage_profile_scores = tools.get_cluster_score_from_embedding(
            contig_codon_usage_profile,
            contig_length_list,
            n_iterations,
            n_components,
            min_dist,
            n_neighbors,
            set_op_mix_ratio,
        )
        codon_usage_profile_scores = dict(
            zip(contig_name_list, codon_usage_profile_scores)
        )
        codon_usage_profile_scores = np.array(
            [codon_usage_profile_scores.get(i, 1.0) for i in self.contigs]
        )
        return codon_usage_profile_scores

    def get_delta_cai(self, quantile=0.25):
        cds_sequences = np.array(self.cds_sequences)
        codon_index = tools.get_codon_index(cds_sequences)
        cai_list = tools.get_cai(cds_sequences, codon_index)
        cai_threshold = np.quantile(cai_list, quantile)
        upper_cds_sequences = cds_sequences[cai_list > cai_threshold]
        new_codon_index = tools.get_codon_index(upper_cds_sequences)
        new_cai_list = tools.get_cai(cds_sequences, new_codon_index)
        return new_cai_list - cai_list

    def compute_cai_scores(self):
        contig_delta_cai = defaultdict(list)
        for cds, delta_cai in zip(self.cds, self.delta_cai):
            contig, _ = cds.rsplit("_", 1)
            contig_delta_cai[contig].append(delta_cai)
        kept_contigs = np.array(self.contigs)[self.n_genes >= 1]
        kept_n_genes = self.n_genes[self.n_genes >= 1]
        mean_contig_delta_cai = [
            np.average(contig_delta_cai[contig]) for contig in kept_contigs
        ]
        # Scale the data and take the average between it and 0.5. This will bring the
        # points closer to the center and will give the KDE some space to breathe.
        mean_contig_delta_cai = np.array(mean_contig_delta_cai)
        mean_contig_delta_cai = (mean_contig_delta_cai - mean_contig_delta_cai.min()) / (
            mean_contig_delta_cai.max() - mean_contig_delta_cai.min()
        )
        mean_contig_delta_cai = (mean_contig_delta_cai + 0.5) / 2
        kernel = ss.gaussian_kde(mean_contig_delta_cai)
        contig_delta_cai_kde = kernel(np.linspace(0, 1, 1000))
        # Find the deepest valley in the KDE.
        valleys = find_peaks(-contig_delta_cai_kde, prominence=(0.035, None))
        if len(valleys[0]):
            max_valley = valleys[0][np.argmax(valleys[1]["prominences"])] / 1000
        else:
            return np.ones(len(self))
        # Compute the KDE again using the number of genes per contig as weight. Then, find
        # the highest peak, where the non-contaminant contigs are concentrated.
        kernel = ss.gaussian_kde(mean_contig_delta_cai, weights=kept_n_genes)
        contig_delta_cai_kde = kernel(np.linspace(0, 1, 1000))
        peaks = find_peaks(contig_delta_cai_kde, height=(None, None))
        if len(peaks[0]) >= 2:
            max_peak = peaks[0][np.argsort(peaks[1]["peak_heights"])[-1]] / 1000
            second_peak = peaks[0][np.argsort(peaks[1]["peak_heights"])[-2]] / 1000
            if np.abs(max_peak - second_peak) < 0.175:
                return np.ones(len(self))
        else:
            return np.ones(len(self))
        divisor = np.abs(max_peak - max_valley)
        if max_peak > max_valley:
            kept_scores = 1 - (max_valley - mean_contig_delta_cai) / divisor
        else:
            kept_scores = 1 - (mean_contig_delta_cai - max_valley) / divisor
        contig_score_dict = dict(zip(kept_contigs, kept_scores))
        scores = np.array([contig_score_dict.get(contig, 1.0) for contig in self.contigs])
        scores = np.clip(scores, 0, 1)
        return scores

    def __len__(self):
        return len(self.contigs)

    def __iter__(self):
        return zip(
            np.repeat(self.genome, len(self)),
            self.contigs,
            np.round(self.codon_usage_cluster_scores, 5),
            np.round(self.cai_scores, 5),
            self.n_genes,
            self.total_cds_length,
            np.round(self.mean_strand_coding_density, 5),
        )


class Composition:
    def __init__(
        self,
        mag,
        composition_dict,
        gc_content,
        n_iterations,
        n_components,
        min_dist,
        n_neighbors,
        set_op_mix_ratio,
    ):
        self.attributes = [
            "genome",
            "contig",
            "tnf_score",
            "gc_content_score",
            "contig_length",
        ]
        self.genome = mag.genome
        self.contigs = mag.contigs
        self.lengths = mag.lengths
        self.tnf = composition_dict[mag.genome]
        self.gc_content = gc_content
        if len(self) <= 2:
            self.tnf_scores = np.ones(len(self))
            self.gc_content_scores = np.ones(len(self))
        else:
            self.tnf_scores = tools.get_cluster_score_from_embedding(
                data=self.tnf,
                lengths=self.lengths,
                n_iterations=n_iterations,
                n_components=n_components,
                min_dist=min_dist,
                n_neighbors=n_neighbors,
                set_op_mix_ratio=set_op_mix_ratio,
            )
            self.gc_content_scores = tools.get_log_ratio_scores(
                self.gc_content, self.lengths, 2
            )

    def __len__(self):
        return len(self.contigs)

    def __iter__(self):
        return zip(
            np.repeat(self.genome, len(self)),
            self.contigs,
            np.round(self.tnf_scores, 5),
            np.round(self.gc_content_scores, 5),
            self.lengths,
        )


class Coverage:
    def __init__(
        self,
        mag,
        coverages,
        min_average_coverage,
        n_iterations,
        n_components,
        min_dist,
        n_neighbors,
        set_op_mix_ratio,
    ):
        self.attributes = [
            "genome",
            "contig",
            "coverage_score",
            "coverage_cluster_score",
            "n_samples",
        ]
        self.genome = mag.genome
        self.contigs = mag.contigs
        self.lengths = mag.lengths
        self.coverages = coverages
        self.selected_samples = (
            np.average(self.coverages, axis=0, weights=self.lengths)
            >= min_average_coverage
        )
        self.n_samples = self.selected_samples.sum()
        if len(self) <= 2 or self.n_samples == 0:
            self.scores = np.ones(len(self))
            self.cluster_scores = np.ones(len(self))
        else:
            self.scores = tools.get_log_ratio_scores(
                self.coverages[:, self.selected_samples], self.lengths, 25
            )
            if self.n_samples >= 3:
                self.cluster_scores = tools.get_cluster_score_from_embedding(
                    data=np.log1p(self.coverages),
                    lengths=self.lengths,
                    n_iterations=n_iterations,
                    n_components=n_components,
                    min_dist=min_dist,
                    n_neighbors=n_neighbors,
                    set_op_mix_ratio=set_op_mix_ratio,
                )
            else:
                self.cluster_scores = np.ones(len(self))

    def __len__(self):
        return len(self.contigs)

    def __iter__(self):
        return zip(
            np.repeat(self.genome, len(self)),
            self.contigs,
            np.round(self.scores, 5),
            np.round(self.cluster_scores, 5),
            np.repeat(self.n_samples, len(self)),
        )


class Taxonomy:
    def __init__(
        self,
        mag,
        mmseqs2_dict,
        contig_min_fraction,
        genome_min_fraction,
        min_genus_identity,
        taxdb,
    ):
        self.attributes = [
            "genome",
            "contig",
            "taxonomy_score",
            "genome_rank",
            "contig_rank",
            "average_protein_identity"
        ]
        self.genome = mag.genome
        self.contigs = mag.contigs
        self.lengths = mag.lengths
        self.taxdb = taxdb
        if mmseqs2_dict:
            self.contigs_in_mmseqs2 = np.array(
                [contig for contig in self.contigs if contig in mmseqs2_dict]
            )
            self.taxid_array = [
                mmseqs2_dict[contig][0] for contig in self.contigs_in_mmseqs2
            ]
            self.identity_array = [
                mmseqs2_dict[contig][1] for contig in self.contigs_in_mmseqs2
            ]
            self.bitscore_array = [
                mmseqs2_dict[contig][2] for contig in self.contigs_in_mmseqs2
            ]
            self.gene_taxonomy = self.get_gene_taxonomy(
                min_genus_identity=min_genus_identity
            )
            self.contig_taxonomy = self.get_contig_taxonomy(fraction=contig_min_fraction)
            self.genome_taxonomy = self.get_genome_taxonomy(fraction=genome_min_fraction)
            self.scores = self.compute_gene_agreement()
            self.average_protein_identity = self.get_average_protein_identity()
        else:
            self.contig_taxonomy = np.array([taxopy.Taxon(1, self.taxdb)] * len(self))
            self.genome_taxonomy = taxopy.Taxon(1, self.taxdb)
            self.scores = np.zeros(len(self))
            self.average_protein_identity = np.zeros(len(self))
        self.genome_rank = len(self.genome_taxonomy.taxid_lineage) - 1
        self.contig_rank = np.array(
            [len(i.taxid_lineage) - 1 for i in self.contig_taxonomy]
        )

    def get_gene_taxonomy(self, min_genus_identity):
        gene_taxonomy_array = []
        for i in range(len(self.taxid_array)):
            contig_taxa_list = []
            for j in range(len(self.taxid_array[i])):
                identity = self.identity_array[i][j]
                taxid = self.taxid_array[i][j]
                if taxid == 0:
                    taxid = 1
                taxon = taxopy.Taxon(taxid, self.taxdb)
                if identity >= min_genus_identity and taxon.rank == "species":
                    taxon = taxopy.Taxon(taxon.taxid_lineage[1], self.taxdb)
                elif identity >= min_genus_identity and taxon.rank == "genus":
                    pass
                elif taxon.rank == "species":
                    taxon = taxopy.Taxon(taxon.taxid_lineage[2], self.taxdb)
                elif taxon.rank == "genus":
                    taxon = taxopy.Taxon(taxon.taxid_lineage[1], self.taxdb)
                contig_taxa_list.append(taxon)
            gene_taxonomy_array.append(contig_taxa_list)
        return gene_taxonomy_array

    def get_contig_taxonomy(self, fraction):
        contig_taxonomy_array = []
        contigs_in_mmseqs2_set = set(self.contigs_in_mmseqs2)
        for contig in self.contigs:
            if contig in contigs_in_mmseqs2_set:
                gene_taxonomy = self.gene_taxonomy[
                    np.where(self.contigs_in_mmseqs2 == contig)[0][0]
                ]
                if len(gene_taxonomy) > 1:
                    genes_bitscore = self.bitscore_array[
                        np.where(self.contigs_in_mmseqs2 == contig)[0][0]
                    ]
                    contig_taxonomy = taxopy.find_majority_vote(
                        gene_taxonomy,
                        self.taxdb,
                        weights=genes_bitscore,
                        fraction=fraction,
                    )
                    # If the contig taxonomic rank is really high, get a new taxonomy with less
                    # conservative parameters
                    if len(contig_taxonomy.taxid_lineage) <= 2:
                        contig_taxonomy = taxopy.find_majority_vote(
                            gene_taxonomy,
                            self.taxdb,
                            weights=genes_bitscore,
                            fraction=fraction - 0.1,
                        )
                else:
                    contig_taxonomy = taxopy.Taxon(gene_taxonomy[0].taxid, self.taxdb)
            else:
                contig_taxonomy = taxopy.Taxon(1, self.taxdb)
            contig_taxonomy_array.append(contig_taxonomy)
        return np.array(contig_taxonomy_array)

    def get_genome_taxonomy(self, fraction):
        if len(self) > 1:
            genome_taxonomy = taxopy.find_majority_vote(
                self.contig_taxonomy,
                self.taxdb,
                weights=self.lengths.tolist(),
                fraction=fraction,
            )
            if len(genome_taxonomy.taxid_lineage) <= 2:
                genome_taxonomy = taxopy.find_majority_vote(
                    self.contig_taxonomy,
                    self.taxdb,
                    weights=self.lengths.tolist(),
                    fraction=fraction - 0.2,
                )
            if len(genome_taxonomy.taxid_lineage) == 3:
                genome_taxonomy = taxopy.find_majority_vote(
                    self.contig_taxonomy,
                    self.taxdb,
                    weights=self.lengths.tolist(),
                    fraction=fraction - 0.1,
                )
        elif len(self) == 1:
            genome_taxonomy = taxopy.Taxon(self.contig_taxonomy[0].taxid, self.taxdb)
        else:
            genome_taxonomy = taxopy.Taxon(1, self.taxdb)
        return genome_taxonomy

    def compute_gene_agreement(self):
        scores_array = []
        contigs_in_mmseqs2_set = set(self.contigs_in_mmseqs2)
        for contig in self.contigs:
            if contig in contigs_in_mmseqs2_set:
                score = 0
                gene_taxonomy = self.gene_taxonomy[
                    np.where(self.contigs_in_mmseqs2 == contig)[0][0]
                ]
                genes_bitscore = self.bitscore_array[
                    np.where(self.contigs_in_mmseqs2 == contig)[0][0]
                ]
                for i in range(len(gene_taxonomy)):
                    if len(self.genome_taxonomy.taxid_lineage) <= len(
                        gene_taxonomy[i].taxid_lineage
                    ):
                        if self.genome_taxonomy.taxid in gene_taxonomy[i].taxid_lineage:
                            score += genes_bitscore[i] / sum(genes_bitscore)
                    else:
                        if gene_taxonomy[i].taxid in self.genome_taxonomy.taxid_lineage:
                            score += genes_bitscore[i] / sum(genes_bitscore)
            else:
                score = 1.0
            scores_array.append(score)
        return np.array(scores_array)

    def get_average_protein_identity(self):
        average_protein_identity_dict = []
        for identity, bitscore in zip(self.identity_array, self.bitscore_array):
            average_protein_identity = np.average(identity, weights=bitscore)
            average_protein_identity_dict.append(average_protein_identity)
        average_protein_identity_dict = dict(
            zip(self.contigs_in_mmseqs2, average_protein_identity_dict)
        )
        return np.array(
            [average_protein_identity_dict.get(contig, 0.0) for contig in self.contigs]
        )

    def __len__(self):
        return len(self.contigs)

    def __iter__(self):
        return zip(
            np.repeat(self.genome, len(self)),
            self.contigs,
            np.round(self.scores, 5),
            np.repeat(self.genome_rank, len(self)),
            self.contig_rank,
            np.round(self.average_protein_identity, 5),
        )


class ContigClassifier:
    def __init__(
        self,
        genome_contig_matrix,
        feature_matrix,
        model_file,
        fast_mode,
        probability_threshold,
        checkm_file,
        threads,
    ):
        self.attributes = ["genome", "contig", "contaminant_probability"]
        self.fast_mode = fast_mode
        self.base_probability_threshold = probability_threshold
        feature_matrix = xgb.DMatrix(feature_matrix)
        model = xgb.Booster({"booster": "gbtree", "nthread": threads, "verbosity": 0})
        model.load_model(model_file)
        self.genomes = genome_contig_matrix[:, 0]
        self.contigs = genome_contig_matrix[:, 1]
        self.probabilities = model.predict(feature_matrix)
        if checkm_file:
            checkm_scores = tools.get_checkm_contamination(checkm_file)
            if checkm_scores:
                logger.info(
                    "Using completeness and contamination estimates to set dynamic thresholds."
                )
                if len(set(self.genomes)) > len(checkm_scores):
                    logging.warning(
                        "The input CheckM tabular file does not contain all the "
                        "genomes being filtered. The default probability threshold "
                        f"({self.base_probability_threshold}) will be used for the "
                        "genomes not found in it."
                    )
                self.probability_threshold = self.get_dynamic_thresholds(checkm_scores)
            else:
                self.probability_threshold = probability_threshold
        else:
            self.probability_threshold = probability_threshold
        self.flagged_contaminants = self.probabilities > self.probability_threshold
        self.mags_contaminants_dict = defaultdict(dict)
        for index, contaminant in enumerate(self.flagged_contaminants):
            genome = self.genomes[index]
            contig = self.contigs[index]
            self.mags_contaminants_dict[genome][contig] = contaminant

    def get_dynamic_thresholds(self, checkm_scores):
        f = lambda score: 0.85 + 0.72 * (np.exp(-0.141 * score) - 1)
        threshold_dict = {genome: f(score) for genome, score in checkm_scores.items()}
        probability_threshold = np.array(
            [
                threshold_dict.get(genome, self.base_probability_threshold)
                for genome in self.genomes
            ]
        )
        probability_threshold = np.clip(probability_threshold, 0.1, 0.85)
        return probability_threshold

    def __len__(self):
        return len(self.probabilities)

    def __iter__(self):
        return zip(self.genomes, self.contigs, np.round(self.probabilities, 5))
