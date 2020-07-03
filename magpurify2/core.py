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
    def __init__(self, mag, min_genes, stringency, prodigal_fna_filepath):
        self.genome = mag.genome
        self.contigs = mag.contigs
        self.cds = []
        self.cds_sequences = []
        self.lengths = []
        for cds, _, cds_sequence in tools.read_fasta(prodigal_fna_filepath):
            self.cds.append(cds)
            self.cds_sequences.append(cds_sequence)
        self.delta_cai = self.get_delta_cai()
        self.scores = self.identify_codon_usage_outliers(min_genes, stringency)

    def get_delta_cai(self, quantile=0.25):
        cds_sequences = np.array(self.cds_sequences)
        codon_index = tools.get_codon_index(cds_sequences)
        cai_list = tools.get_cai(cds_sequences, codon_index)
        cai_threshold = np.quantile(cai_list, quantile)
        upper_cds_sequences = cds_sequences[cai_list > cai_threshold]
        new_codon_index = tools.get_codon_index(upper_cds_sequences)
        new_cai_list = tools.get_cai(cds_sequences, new_codon_index)
        return new_cai_list - cai_list

    def identify_codon_usage_outliers(self, min_genes, stringency):
        contig_delta_cai = defaultdict(list)
        for cds, delta_cai in zip(self.cds, self.delta_cai):
            contig, _ = cds.rsplit("_", 1)
            contig_delta_cai[contig].append(delta_cai)
        kept_contigs, n_genes, mean_contig_cai = [], [], []
        for contig, delta_cai_list in contig_delta_cai.items():
            if len(delta_cai_list) >= min_genes:
                kept_contigs.append(contig)
                mean_contig_cai.append(np.mean(delta_cai_list))
                n_genes.append(len(delta_cai_list))
        unit_mean_contig_cai = tools.zscore(mean_contig_cai, unit_interval=True)
        kernel = ss.gaussian_kde(unit_mean_contig_cai, weights=n_genes)
        contig_cai_kde = kernel(np.linspace(0, 1, 1000))
        # Find the deepest valley in the KDE.
        valleys = find_peaks(-contig_cai_kde, prominence=(0.05, None))
        if len(valleys[0]):
            max_valley = valleys[0][np.argmax(valleys[1]["prominences"])] / 1000
        else:
            return np.ones(len(self))
        # Find the highest peak, where the non-contaminant contigs are concentrated.
        peaks = find_peaks(contig_cai_kde, height=(None, None))
        if len(peaks[0]):
            max_peak = peaks[0][np.argmax(peaks[1]["peak_heights"])] / 1000
        else:
            return np.ones(len(self))
        if max_peak > max_valley:
            threshold = max_valley - np.abs(max_peak - max_valley) * stringency
            scores = (unit_mean_contig_cai >= threshold).astype(float)
        else:
            threshold = max_valley + np.abs(max_peak - max_valley) * stringency
            scores = (unit_mean_contig_cai <= threshold).astype(float)
        contig_score_dict = dict(zip(kept_contigs, scores))
        return np.array([contig_score_dict.get(contig, 0.0) for contig in self.contigs])

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
        coverage_dict,
        n_iterations,
        n_components,
        min_dist,
        n_neighbors,
        set_op_mix_ratio,
        max_deviation,
    ):
        self.genome = mag.genome
        self.contigs = mag.contigs
        self.lengths = mag.lengths
        self.coverages = np.array(
            [
                coverage_dict[contig]
                if contig in coverage_dict
                else [0.0] * len(list(coverage_dict.values())[0])
                for contig in self.contigs
            ]
        )
        selected_samples = self.coverages.mean(axis=0) >= 1
        self.coverages = self.coverages[:, selected_samples]
        if self.coverages.shape[1] > 1:
            self.scores = tools.get_cluster_score_from_embedding(
                data=np.log1p(self.coverages),
                lengths=self.lengths,
                n_iterations=n_iterations,
                n_components=n_components,
                min_dist=min_dist,
                n_neighbors=n_neighbors,
                set_op_mix_ratio=set_op_mix_ratio,
            )
        else:
            self.scores = self.identify_coverage_outliers(max_deviation)

    def identify_coverage_outliers(self, max_deviation=5.0):
        data = np.array(self.coverages).flatten()
        kernel = ss.gaussian_kde(data, weights=self.lengths)
        data_kde = kernel(np.linspace(0, np.max(data), 1000))
        peaks = find_peaks(data_kde, height=(None, None))
        if len(peaks[0]):
            threshold = (
                peaks[0][np.argmax(peaks[1]["peak_heights"])] / 1000 * np.max(data)
            )
            limits = sorted([max_deviation * threshold, threshold / max_deviation])
            scores = np.logical_and((data <= limits[1]), (data >= limits[0])).astype(
                float
            )
        else:
            scores = np.ones(len(data))
        return scores

    def __len__(self):
        return len(self.contigs)

    def __iter__(self):
        return zip(self.contigs, self.scores)


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
