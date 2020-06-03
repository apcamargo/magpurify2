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
import taxopy

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


class Composition:
    def __init__(self, mag, composition_dict):
        self.genome = mag.genome
        self.contigs = mag.contigs
        self.lengths = mag.lengths
        self.tnf = composition_dict[mag.genome]
        self.embedding = tools.create_embedding(
            data=self.tnf, n_components=3, min_dist=0.1, n_neighbors=15, set_op_mix_ratio=1,
        )
        self.scores = tools.compute_contig_cluster_score(
            data=self.embedding, allow_single_cluster=True, lengths=self.lengths,
        )

    def __len__(self):
        return len(self.contigs)

    def __iter__(self):
        return zip(self.contigs, self.scores)


class Coverage:
    def __init__(self, mag, coverage_dict):
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
        self.scores = tools.compute_contig_cluster_score(
            data=self.coverages, allow_single_cluster=True, lengths=self.lengths,
        )

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
        self.contig_taxonomy = self.get_contig_taxonomy(taxonomy_dict, contig_min_fraction, allow_genus)
        self.genome_taxonomy = self.get_genome_taxonomy(genome_min_fraction)
        self.scores = self.compute_gene_agreement()

    def get_contig_taxonomy(self, taxonomy_dict, fraction=0.75, allow_genus=False):
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
