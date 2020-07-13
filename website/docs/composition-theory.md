# Composition module

Contigs assembled from reads derived from the same genome or from genomes of close organisms tend to display a similar sequence composition profile, represented as k-mer frequencies. Thus, metagenomic binners use 4-mer frequencies, or tetranucleotide frequencies (TNFs), to cluster contigs into putative genomic bins.

Whether or not two given contigs will be clustered into the same genomic bin does not depend exclusively on their TNF profiles. In most modern clustering algorithms local relationships are influenced by other data points, meaning that a given pair of contigs may end up in the same bin or not, depending on the full set contigs that is being clustered. Moreover, most binners also use sequencing coverage information in addition to TNF data to cluster contigs, which may lead to genomic bins that encompass contigs with distinct TNF profiles.

MAGpurify2 processes each genomic bin individually and finds potential contaminants with respect to TNF profile by identifying contigs that fall outside of the "core cluster" within the bin.

## How it works

To identify putative contaminants within a genomic bin, MAGpurify2: (1) computes the TNF profile of each contig, (2) embbeds data points into a low-dimentional space using a non-linear transformation, and (3) finds the "core cluster" and computes each contig score.

The four canonical DNA bases (A, T, C and G) can produce 256 distinct 4-mers, however in a strand-independent analysis reverse complement k-mers (eg.: `TTAC` and `GTAA`) are redundant and should be count as a single entity in order to reduce memory usage and data variance. Thus, MAGpurify2 counts the relative frequencies of 136 canonical 4-mers for each contig within the bin.

::: warning TNF profile of short sequences
Short contigs contain a reduced number of 4-mers and thus display greater data variance and statistical uncertainty than longer contigs. That is one of the reasons that most binners filter out contigs shorter than a set threshold (usually around 2,000 bp). MAGpurify2 currently doesn't take into account the length-dependent uncertainty of TNF estimation when identifying putative contaminants.
:::

