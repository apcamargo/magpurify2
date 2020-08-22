# Introduction

**ADD ILLUSTRATION: genomes ➔ contigs ➔ binning ➔ magpurify2**

Metagenome-assembled genomes (MAGs) are genomic sequences extracted from broader metagenomic assemblies. This is acomplished through the binning of contigs that share similar sequence composition and exhibit correlated sequencing coverage across different environmental conditions.

Due to technical limitations of sequencing, metagenomic assembly and binning, MAGs are usually fragmented, incomplete and display varying degrees of contamination (i.e. sequences from a different organism). Under the typical assumption that a given MAG represents a single species, this last issue is especially concearning as contaminant sequences will lead to unreliable biological interpretations.

MAGpurify2 aims to improve the quality of MAGs through the identification and removal of contaminant sequences, improving the reliability of downstream analysis of the genomic sequences.