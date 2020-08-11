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
    def __init__(self, mag, min_genes, prodigal_fna_filepath):
        self.attributes = [
            "genome",
            "contig",
            "codon_usage_score",
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
        self.mean_strand_coding_density = np.round(self.mean_strand_coding_density, 5)
        self.delta_cai = self.get_delta_cai()
        if len(self) == 1:
            self.scores = np.array([1.0])
        else:
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
            [self.genome] * len(self),
            self.contigs,
            np.round(self.scores, 5),
            self.n_genes,
            self.total_cds_length,
            self.mean_strand_coding_density,
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
        self.attributes = ["genome", "contig", "tnf_score", "gc_score", "length"]
        self.genome = mag.genome
        self.contigs = mag.contigs
        self.lengths = mag.lengths
        self.tnf = composition_dict[mag.genome]
        self.gc_content = gc_content
        if len(self) <= 2:
            self.tnf_scores = np.array([1.0] * len(self))
            self.gc_scores = np.array([1.0] * len(self))
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
            self.gc_scores = tools.get_log_ratio_scores(self.gc_content, self.lengths, 2)

    def __len__(self):
        return len(self.contigs)

    def __iter__(self):
        return zip(
            [self.genome] * len(self),
            self.contigs,
            np.round(self.tnf_scores, 5),
            np.round(self.gc_scores, 5),
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
            self.scores = np.array([1.0] * len(self))
            self.cluster_scores = np.array([1.0] * len(self))
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
            [self.genome] * len(self),
            self.contigs,
            np.round(self.scores, 5),
            np.round(self.cluster_scores, 5),
            [self.n_samples] * len(self),
        )


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
        self.attributes = ["genome", "contig", "taxonomy_score"]
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
        return zip([self.genome] * len(self), self.contigs, np.round(self.scores, 5))
