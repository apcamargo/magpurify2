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

import taxopy

from magpurify2 import tools


class Mag:
    def __init__(self, filepath, store_sequences=True):
        self.genome = tools.get_genome_name(filepath)
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
    def __init__(self, mag, composition_dict, strictness, threads):
        self.genome = mag.genome
        self.contigs = mag.contigs
        self.lengths = mag.lengths
        self.strictness = strictness
        self.threads = threads
        self.tnf = composition_dict[mag.genome]
        self.embedding, self.contaminants = self.identify_contaminants()

    def identify_contaminants(self):
        embedding = tools.create_embedding(
            data=self.tnf, min_dist=0.1, n_neighbors=5, set_op_mix_ratio=1,
        )
        contaminants = tools.identify_contaminant_clusters(
            embedding=embedding,
            base_bandwidth=7,
            lengths=self.lengths,
            strictness=self.strictness,
            threads=self.threads,
        )
        return embedding, contaminants

    def __len__(self):
        return len(self.contigs)

    def __iter__(self):
        return zip(self.contigs, self.contaminants)


class Coverage:
    def __init__(self, mag, coverage_dict, strictness, threads):
        self.genome = mag.genome
        self.contigs = mag.contigs
        self.lengths = mag.lengths
        self.strictness = strictness
        self.threads = threads
        self.coverages = [
            coverage_dict[contig]
            if contig in coverage_dict
            else [0.0] * len(list(coverage_dict.values())[0])
            for contig in self.contigs
        ]
        self.embedding, self.contaminants = self.identify_contaminants()

    def identify_contaminants(self):
        embedding = tools.create_embedding(
            data=self.coverages, min_dist=0.05, n_neighbors=15, set_op_mix_ratio=0.25,
        )
        contaminants = tools.identify_contaminant_clusters(
            embedding=embedding,
            base_bandwidth=5.5,
            lengths=self.lengths,
            strictness=self.strictness,
            threads=self.threads,
        )
        return embedding, contaminants

    def __len__(self):
        return len(self.contigs)

    def __iter__(self):
        return zip(self.contigs, self.contaminants)


class Taxonomy:
    def __init__(self, mag, taxonomy_dict, strictness, taxdb):
        self.genome = mag.genome
        self.contigs = mag.contigs
        self.lengths = mag.lengths
        self.strictness = strictness
        self.contig_taxonomy = []
        self.taxdb = taxdb
        for contig in self.contigs:
            contig_taxon = taxonomy_dict[mag.genome].get(contig, taxopy.Taxon("1", taxdb))
            self.contig_taxonomy.append(contig_taxon)
        self.genome_taxonomy = self.get_genome_taxonomy()
        self.contaminants = self.identify_contaminants()

    def get_genome_taxonomy(self):
        selected_taxon = "1"
        threshold = (self.strictness + 0.5) / 2 * sum(self.lengths)
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

    def identify_contaminants(self):
        genome_taxonomy_rank = len(self.genome_taxonomy.taxid_lineage) - 1
        contaminants = []
        for taxon in self.contig_taxonomy:
            if len(taxon.taxid_lineage) >= genome_taxonomy_rank + 1:
                if (
                    list(reversed(taxon.taxid_lineage))[genome_taxonomy_rank]
                    != self.genome_taxonomy.taxid
                ):
                    contaminants.append(True)
                else:
                    contaminants.append(False)
            else:
                contaminants.append(False)
        return contaminants

    def __len__(self):
        return len(self.contigs)

    def __iter__(self):
        return zip(self.contigs, self.contaminants)
