# MAGpurify2

Identify and remove contaminants from metagenome-assembled genomes.

## Quick start


Download test data:

```bash
fileId="1-Gf-FsVIcARqrUb-LHS_FZlb-sCGAmfo"
fileName="magpurify2_test_data.tar.gz"
curl -sc /tmp/cookie "https://drive.google.com/uc?export=download&id=${fileId}" > /dev/null
code="$(awk '/_warning_/ {print $NF}' /tmp/cookie)"
curl -Lb /tmp/cookie "https://drive.google.com/uc?export=download&confirm=${code}&id=${fileId}" -o ${fileName}
```

Download database:

```bash
fileId="14qAK-4NJNu0sOPldnZhO8vcAFpb1BcrB"
fileName="magpurify2DB.v1.0.tar.gz"
curl -sc /tmp/cookie "https://drive.google.com/uc?export=download&id=${fileId}" > /dev/null
code="$(awk '/_warning_/ {print $NF}' /tmp/cookie)"
curl -Lb /tmp/cookie "https://drive.google.com/uc?export=download&confirm=${code}&id=${fileId}" -o ${fileName}
```

Execute the pipeline:

```bash
magpurify2 composition test_data/genomes/* output
magpurify2 coverage test_data/genomes/* output --bam_files test_data/bam_files/*
magpurify2 taxonomy test_data/genomes/* output magpurify2DB
magpurify2 filter test_data/genomes/* output filtered_genomes
```
