# Installation

## Installing MAGpurify2

```bash
conda create -n magpurify2 -c conda-forge -c bioconda python=3.8 prodigal mmseqs2
conda activate magpurify2
pip install magpurify2-0.1.0-cp38-cp38-manylinux2010_x86_64.whl
```

## Download test data and database

- Download test data:

```bash
fileId="1-Gf-FsVIcARqrUb-LHS_FZlb-sCGAmfo"
fileName="magpurify2_test_data.tar.gz"
curl -sc /tmp/cookie "https://drive.google.com/uc?export=download&id=${fileId}" > /dev/null
code="$(awk '/_warning_/ {print $NF}' /tmp/cookie)"
curl -Lb /tmp/cookie "https://drive.google.com/uc?export=download&confirm=${code}&id=${fileId}" -o ${fileName}
tar zxfv magpurify2_test_data.tar.gz
```

- Download database:

```bash
fileId="1ooWiR3LplBy5GsY5wZ7o6dwswiCWVvmi"
fileName="magpurify2DB.v1.0.tar.gz"
curl -sc /tmp/cookie "https://drive.google.com/uc?export=download&id=${fileId}" > /dev/null
code="$(awk '/_warning_/ {print $NF}' /tmp/cookie)"
curl -Lb /tmp/cookie "https://drive.google.com/uc?export=download&confirm=${code}&id=${fileId}" -o ${fileName}
tar zxfv magpurify2DB.v1.0.tar.gz
```

## External dependencies

CheckV depends on…

If you installed MAGpurify2 via Conda the dependencies should have been installed automatically.
