# Using MAGpurify2

## Generating bins

*Recommend binners and DAS Tool.*

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
BAM files store read alignment information to the target metagenome (or MAG) and used by MAGpurify2 to estimate the coverage of each contig. To generate the BAM inputs to MAGpurify2 you should first map your reads to the complete metagenome using a proper tool (such as [Bowtie 2](https://github.com/BenLangmead/bowtie2), [minimap2](https://github.com/lh3/minimap2) or [BWA-MEM2](https://github.com/bwa-mem2/bwa-mem2)) and then convert and sort the output using [samtools](https://github.com/samtools/samtools). For example:

```bash
# Create a Bowtie 2 index for your metagenome inside the 'bt2' directory:
$ mkdir bt2
$ bowtie2-build --threads 4 metagenome.fna bt2/metagenome
# Map the reads, sort the output and write it to 'sample1.bam':
$ bowtie2 --threads 4 -x bt2/metagenome \
  -1 sample1_R1.fastq.gz -2 sample1_R2.fastq.gz \
  | samtools sort -@ 4 -o sample1.bam -
```

We recommend mapping the reads to the metagenome (superset) and not directly to the MAGs retrieved from it (subsets). There are two main reasons for that:

- When you map the reads to the MAG there's a chance that reads that were originated from the sequencing of a closely related genome will be erroneously aligned to the MAG (cross-mapping), introducing bias to the coverage estimation.
- Metagenome-wide mappings can be used to estimate the coverage of all the contigs in the metagenome, thus allowing MAGpurify2 to process multiple MAGs in a single execution.

If the target MAGs are derived from multiple source metagenomes you need input BAM files containing read mappings to each one of them.
:::

If you don't have access to the raw sequencing data or to previously generated BAM files you can input coverage data stored in a tab-separated values (TSV) file. To do so, you should use the `--coverage_file` argument:

```bash
coverage genomes/* output --coverage_file contig_coverages.tsv --threads 4
```

The first column of the coverage file must store the contigs names. The remaining columns should contain the coverage of each contig across multiple samples, as shown in the example below:

```
contig_1     15.744    12.605    25.148    3.728    0.000
contig_2     34.466    48.019    18.222    3.707    4.195
contig_3      0.000    22.356    21.944    4.479    4.463
contig_4     14.201     9.993     0.000    0.925    4.608
contig_5     17.179    12.280    56.643    3.586    4.226
contig_6      5.239     8.430    5.2070    3.988    0.000
contig_7     17.737    16.005    29.692    4.039    4.190
contig_8      0.000    15.866    13.663    0.877    2.269
contig_9     19.129    15.145    21.249    0.000    2.342
contig_10    11.074     9.574    19.673    0.934    4.732
```