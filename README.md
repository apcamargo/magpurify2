<img src="https://apcamargo.github.io/magpurify2/logo.svg" align="left" width="200" height="155px"/>
<img align="left" width="0" height="155px" hspace="20"/>

Improve the quality of metagenome-assembled genomes by identifying and removing contaminant sequences with an easy-to-use and modular command-line interface.

[![MIT License](https://img.shields.io/badge/license-MIT-007EC7.svg?style=flat-square)](/LICENSE) [![Fish Shell Version](https://img.shields.io/badge/fish-≥v2.2.0-007EC7.svg?style=flat-square)](http://fishshell.com) [![Travis Build Status](http://img.shields.io/travis/oh-my-fish/oh-my-fish.svg?style=flat-square)](https://travis-ci.org/oh-my-fish/oh-my-fish)

You can find MAGpurify2's full documentation at its [website](https://apcamargo.github.io/magpurify2/).

---

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
