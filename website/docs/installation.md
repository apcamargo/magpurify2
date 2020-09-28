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

### External dependencies

If you choose install MAGpurify2 via `pip`, make sure that you have also installed the third-party dependencies: [Prodigal](https://github.com/hyattpd/Prodigal) [^1] and [MMseqs2](https://github.com/soedinglab/MMseqs2) [^2]. The `conda` installation method will automatically download and install these software for you.

::: tip External dependencies are not needed in fast mode
If you are executing MAGpurify2 with the `--fast_mode` parameter you don't need to install Progidal or MMseqs2. These software are only required to run the `codon_usage` and `taxonomy` modules.
:::

## Docker

Foo.

## Test the installation

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

[^1]: Hyatt, Doug, et al. ["Prodigal: prokaryotic gene recognition and translation initiation site identification."](https://pubmed.ncbi.nlm.nih.gov/20211023/) *BMC Bioinformatics* 11.1 (2010): 119.

[^2]: Steinegger, Martin, and Johannes Söding. ["MMseqs2 enables sensitive protein sequence searching for the analysis of massive data sets."](https://pubmed.ncbi.nlm.nih.gov/29035372/) *Nature Biotechnology* 35.11 (2017): 1026-1028.