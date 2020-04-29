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