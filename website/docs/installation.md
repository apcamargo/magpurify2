# Installation

## Installing MAGpurify2

MAGpurify2 can be installed using `pip` or `conda`. Alternatively, you can execute it via [Docker](#docker).

<code-group>
<code-block title="pip" active>
```bash
pip install magpurify2
```
</code-block>

<code-block title="conda">
```bash
conda install -c conda-forge -c bioconda magpurify2
```
</code-block>
</code-group>

If you choose install MAGpurify2 via `pip`, make sure that you have also installed the third-party dependencies: [Prodigal](https://github.com/hyattpd/Prodigal) and [MMseqs2](https://github.com/soedinglab/MMseqs2). The `conda` install method will automatically download and install these software.

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

## Docker

Foo.