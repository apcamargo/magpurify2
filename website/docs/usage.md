# Using MAGpurify2

## Quick start

```bash
magpurify2 composition test_data/genomes/* output
magpurify2 coverage test_data/genomes/* output --bam_files test_data/bam_files/*
magpurify2 taxonomy test_data/genomes/* output magpurify2DB
magpurify2 filter test_data/genomes/* output filtered_genomes
```

## The `composition` module

`magpurify2 composition` is the command used to identify putative contaminants using tetranucleotide frequencies.

```
usage: magpurify2 composition [-h] [-s STRICTNESS] [-t THREADS] [--quiet] genomes [genomes ...] output_directory

positional arguments:
  genomes               Input genomes in the FASTA format.
  output_directory      Directory to write the output files to.

optional arguments:
  -h, --help            show this help message and exit
  -s STRICTNESS, --strictness STRICTNESS
                        Strictness of the contaminant detection algorithm. Must be a number between 0 (less strict) and 1 (more strict). (default: 0.5)
  -t THREADS, --threads THREADS
                        Number of threads to use. All by default. (default: 4)
  --quiet               Suppress the logger output (default: False)
```

## The `coverage` module

::: tip How to generate the BAM files
BAM files store reads that map to the target metagenome or MAG and are processed by MAGpurify2 to estimate the coverage of each contig. To generate the BAM required by MAGpurify2 you should first map your reads to the complete metagenome using a proper tool (such as [Bowtie 2](https://github.com/BenLangmead/bowtie2), [minimap2](https://github.com/lh3/minimap2) or [BWA-MEM2](https://github.com/bwa-mem2/bwa-mem2)) and then convert and sort the output using [samtools](https://github.com/samtools/samtools).

```bash
$ mkdir bt2
$ bowtie2-build --threads 4 metagenome.fna bt2/metagenome.fna
$ bowtie2 --threads 6 -x bt2/metagenome.fna \
  -1 sample1_R1.fastq.gz -2 sample1_R2.fastq.gz \
  | samtools sort -@ 6 -o sample1.bam -
```

We reccomend mapping the reads to the metagenome and not directly to the MAGs and this is because of two factors:
- When you map the reads into the MAG a read that was originated from the sequencing of a closely related genome might be erroneously aligned to the MAG (cross-mapping), introducing bias to the coverage estimation.
- Metagenome-wide mappings can be used to estimate the coverage of all the contigs in the metagenome, thus allowing MAGpurify2 to process multiple MAGs in a single execution.
:::