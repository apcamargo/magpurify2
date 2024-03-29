# Introduction

**ADD ILLUSTRATION: genomes ➔ contigs ➔ binning ➔ magpurify2**

Metagenome-assembled genomes (MAGs) are draft genome sequences that are recovered from metagenomes containing genomic fragments of multiple species. The retrieval of MAG sequences from metagenomes is accomplished through the clustering — in this context, also known as binning — of contigs that share similar sequence composition and exhibit correlated sequencing coverage across different environmental samples. Due to technical limitations inherent to sequencing, metagenomic assembly and binning, MAGs are usually fragmented, incomplete, and have varying degrees of contamination (I.e. sequences from a different organism). Under the typical assumption that a MAG represents a composite of a single species, contaminant sequences are especially concerning as they will cause unreliable biological interpretations.

MAGpurify2 aims to improve the quality of MAGs through the identification and removal of contaminant sequences, improving the reliability of downstream analysis of the genomic sequences. It can be be used to process thousands of genomes at once, making it especially useful for large-scale studies.