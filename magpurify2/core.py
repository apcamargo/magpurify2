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

from collections import defaultdict

import numpy as np
import scipy.stats as ss
import taxopy
from scipy.signal import find_peaks

from magpurify2 import tools


class Mag:
    def __init__(self, filepath, store_sequences=True):
        self.genome = filepath.stem
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
    def __init__(self, mag, min_genes, prodigal_fna_filepath):
        self.genome = mag.genome
        self.contigs = mag.contigs
        self.lengths = mag.lengths
        self.cds = []
        self.cds_sequences = []
        if len(self) == 1:
            self.scores = np.array([1.0])
        else:
            n_genes_dict = defaultdict(int)
            for cds, _, cds_sequence in tools.read_fasta(prodigal_fna_filepath):
                self.cds.append(cds)
                self.cds_sequences.append(cds_sequence)
                contig, _ = cds.rsplit("_", 1)
                n_genes_dict[contig] += 1
            self.n_genes = np.array(
                [n_genes_dict.get(contig, 0) for contig in self.contigs]
            )
            self.delta_cai = self.get_delta_cai()
            self.scores = self.compute_codon_usage_scores(min_genes)

    def get_delta_cai(self, quantile=0.25):
        cds_sequences = np.array(self.cds_sequences)
        codon_index = tools.get_codon_index(cds_sequences)
        cai_list = tools.get_cai(cds_sequences, codon_index)
        cai_threshold = np.quantile(cai_list, quantile)
        upper_cds_sequences = cds_sequences[cai_list > cai_threshold]
        new_codon_index = tools.get_codon_index(upper_cds_sequences)
        new_cai_list = tools.get_cai(cds_sequences, new_codon_index)
        return new_cai_list - cai_list

    def compute_codon_usage_scores(self, min_genes):
        contig_delta_cai = defaultdict(list)
        for cds, delta_cai in zip(self.cds, self.delta_cai):
            contig, _ = cds.rsplit("_", 1)
            contig_delta_cai[contig].append(delta_cai)
        kept_contigs = np.array(self.contigs)[self.n_genes >= min_genes]
        kept_n_genes = self.n_genes[self.n_genes >= min_genes]
        mean_contig_delta_cai = [
            np.average(contig_delta_cai[contig]) for contig in kept_contigs
        ]

        # Scale the data and take the average between it and 0.5. This will bring the
        # points closer to the center and will give the KDE some space to breathe.
        # Scale gene-level delta cai:
        delta_cai = (self.delta_cai - self.delta_cai.min()) / (
            self.delta_cai.max() - self.delta_cai.min()
        )
        delta_cai = (delta_cai + 0.5) / 2
        # Scale contig-level delta cai:
        mean_contig_delta_cai = np.array(mean_contig_delta_cai)
        mean_contig_delta_cai = (mean_contig_delta_cai - mean_contig_delta_cai.min()) / (
            mean_contig_delta_cai.max() - mean_contig_delta_cai.min()
        )
        mean_contig_delta_cai = (mean_contig_delta_cai + 0.5) / 2

        kernel = ss.gaussian_kde(delta_cai)
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
        scores[scores > 1] = 1
        scores[scores < 0] = 0
        return scores

    def __len__(self):
        return len(self.contigs)

    def __iter__(self):
        return zip(self.contigs, self.scores)


class Composition:
    def __init__(
        self,
        mag,
        composition_dict,
        n_iterations,
        n_components,
        min_dist,
        n_neighbors,
        set_op_mix_ratio,
    ):
        self.genome = mag.genome
        self.contigs = mag.contigs
        self.lengths = mag.lengths
        self.tnf = composition_dict[mag.genome]
        if len(self) <= 2:
            self.scores = np.array([1.0] * len(self))
        else:
            self.scores = tools.get_cluster_score_from_embedding(
                data=self.tnf,
                lengths=self.lengths,
                n_iterations=n_iterations,
                n_components=n_components,
                min_dist=min_dist,
                n_neighbors=n_neighbors,
                set_op_mix_ratio=set_op_mix_ratio,
            )

    def __len__(self):
        return len(self.contigs)

    def __iter__(self):
        return zip(self.contigs, self.scores)


class Coverage:
    def __init__(
        self,
        mag,
        coverages,
        min_average_coverage,
        use_clustering,
        n_iterations,
        n_components,
        min_dist,
        n_neighbors,
        set_op_mix_ratio,
    ):
        # The `_iterate_cluster_scores` attribute controls whether the iterator will yield
        # scores generated by the relative error method or clustering methods.
        self._iterate_cluster_scores = False
        self.genome = mag.genome
        self.contigs = mag.contigs
        self.lengths = mag.lengths
        self.use_clustering = use_clustering
        self.coverages = coverages
        self.selected_samples = (
            np.average(self.coverages, axis=0, weights=self.lengths)
            >= min_average_coverage
        )
        if len(self) <= 2 or self.selected_samples.sum() == 0:
            self.scores = np.array([1.0] * len(self))
            self.cluster_scores = np.array([1.0] * len(self))
        else:
            self.scores = self.log_relative_error_scores()
            if self.use_clustering and self.selected_samples.sum() >= 3:
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

    def set_iterate_cluster_scores(self, iterate_cluster_scores=True):
        self._iterate_cluster_scores = iterate_cluster_scores

    def log_relative_error_scores(self):
        weighted_medians = np.apply_along_axis(
            tools.weighted_median, 0, self.coverages, weights=self.lengths
        )
        selected_data = self.coverages[:, self.selected_samples]
        weighted_medians = weighted_medians[self.selected_samples]
        deviation = selected_data / (weighted_medians + 1e-5)
        scores = 1 - np.abs(np.log(deviation + 1e-5) / np.log(25))
        scores = np.average(scores, axis=1)
        scores[scores < 0] = 0
        return scores

    def __len__(self):
        return len(self.contigs)

    def __iter__(self):
        if not self._iterate_cluster_scores:
            return zip(self.contigs, self.scores)
        else:
            return zip(self.contigs, self.cluster_scores)


class Taxonomy:
    def __init__(
        self,
        mag,
        taxonomy_dict,
        contig_min_fraction,
        genome_min_fraction,
        allow_genus,
        taxdb,
    ):
        self.genome = mag.genome
        self.contigs = mag.contigs
        self.lengths = mag.lengths
        self.taxdb = taxdb
        self.gene_taxonomy = [
            taxonomy_dict[mag.genome].get(contig, [taxopy.Taxon("1", self.taxdb)])
            for contig in self.contigs
        ]
        self.contig_taxonomy = self.get_contig_taxonomy(
            fraction=contig_min_fraction, allow_genus=allow_genus
        )
        self.genome_taxonomy = self.get_genome_taxonomy(fraction=genome_min_fraction)
        self.scores = self.compute_gene_agreement()

    def get_contig_taxonomy(self, fraction=0.75, allow_genus=False):
        contig_taxonomy = []
        for contig in self.gene_taxonomy:
            if len(contig) > 1:
                majority_vote = taxopy.find_majority_vote(contig, self.taxdb, fraction)
            else:
                (majority_vote,) = contig
            if allow_genus and majority_vote.rank == "species":
                majority_vote = taxopy.Taxon(majority_vote.taxid_lineage[1], self.taxdb)
            elif majority_vote.rank == "species":
                majority_vote = taxopy.Taxon(majority_vote.taxid_lineage[2], self.taxdb)
            elif majority_vote.rank == "genus":
                majority_vote = taxopy.Taxon(majority_vote.taxid_lineage[1], self.taxdb)
            contig_taxonomy.append(majority_vote)
        return contig_taxonomy

    def get_genome_taxonomy(self, fraction=0.75):
        selected_taxon = "1"
        threshold = fraction * sum(self.lengths)
        taxa_weights = defaultdict(lambda: defaultdict(int))
        for taxon, length in zip(self.contig_taxonomy, self.lengths):
            for rank, taxid in enumerate(reversed(taxon.taxid_lineage)):
                taxa_weights[rank][taxid] += length
        for rank in reversed(range(8)):
            if (
                taxa_weights[rank].values()
                and max(taxa_weights[rank].values()) > threshold
            ):
                selected_taxon = max(taxa_weights[rank], key=taxa_weights[rank].get)
                break
        return taxopy.Taxon(selected_taxon, self.taxdb)

    def compute_gene_agreement(self):
        gene_agreement = []
        for gene_taxonomy in self.gene_taxonomy:
            score = 0
            for gene in gene_taxonomy:
                lca = taxopy.find_lca([gene, self.genome_taxonomy], self.taxdb)
                if len(self.genome_taxonomy.taxid_lineage) <= len(gene.taxid_lineage):
                    if lca.taxid == self.genome_taxonomy.taxid:
                        score += 1
                else:
                    if lca.taxid == gene.taxid:
                        score += 1
            score /= len(gene_taxonomy)
            gene_agreement.append(score)
        return gene_agreement

    def __len__(self):
        return len(self.contigs)

    def __iter__(self):
        return zip(self.contigs, self.scores)
