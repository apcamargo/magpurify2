# Codon usage module

Lorem ipsum dolor sit amet, consectetur adipiscing elit. Proin consectetur rutrum tellus sit amet egestas. Aliquam faucibus condimentum mauris vitae rhoncus. Praesent laoreet rutrum purus, sed sodales neque. Quisque ut tincidunt erat, at rhoncus orci. Aliquam condimentum eget libero sed suscipit. Duis nec eros nec augue porta malesuada. Nulla nec nisl justo. Mauris eget augue a leo laoreet aliquam.

Assumes the [standard genetic code](https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi#SG1), which is true for the majority of bacterial genomes/genes.

## How it works

To summarize the genome codon usage bias MAGpurify2 uses the the Relative Synonymous Codon Usage (RSCU)

and the codon usage adaptation of each gene relative to the rest of the genome is measured using the Codon Adaptation Index (CAI)

described by [Sharp PM and Li WH](https://pubmed.ncbi.nlm.nih.gov/3547335/).

$$
\mathit{RSCU}_{i,j} = \frac{x_{i,j}}{\frac{1}{n_i}\sum_{j=1}^{n_i}x_{i,j}\:}
$$

Aliquam nulla mi, commodo sit amet augue quis, mattis scelerisque nisl. Aliquam aliquam gravida odio ac facilisis. Curabitur in nisi egestas, molestie elit et, pellentesque lectus. Curabitur accumsan mollis pretium. In purus elit, bibendum a mauris vitae, finibus posuere justo. Maecenas non porta dui, vitae commodo lectus. Sed et elit non tortor tristique laoreet ac sed felis.

::: tip Codon usage bias in secondary chromosomes
[Genes on secondary chromosomes exhibited significantly weaker codon usage bias than those on primary chromosomes.](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1000732)
:::

$$
w_{i,j} = \frac{\mathit{RSCU}_{i,j}}{\max_{1\le j\le n}\mathit{RSCU}_{i,j}}
$$

Nullam rhoncus dui eu iaculis congue. Sed elit neque, ultricies eu venenatis nec, vehicula et nisl. Nullam tincidunt dapibus leo, faucibus venenatis justo rutrum et. In hac habitasse platea dictumst. Ut nunc risus, faucibus bibendum gravida vitae, bibendum vel ex.

$$
\mathit{CAI} = \left(\prod \:_{k=1}^Lw_k\right)^{\frac{1}{L}}
$$

Maecenas sed gravida lectus. Nulla aliquam purus ut justo bibendum ornare. Morbi vestibulum congue urna, quis dapibus urna. Nam ligula augue, rutrum a pellentesque non, ornare id risus. Maecenas malesuada, ipsum ac vehicula fermentum, turpis sapien sollicitudin arcu, a pulvinar dui magna vitae est. Aenean nec purus maximus, pulvinar lacus sed, fermentum magna. Nullam feugiat fringilla diam in auctor.

![gene-cai](./figures/gene-cai.svg)

Etiam placerat dui id mi fermentum, eget semper sem eleifend. In finibus tincidunt massa sit amet posuere. Integer ultricies nibh ac nunc iaculis, vel vehicula mi imperdiet. Nam consequat nunc convallis eros condimentum, non tincidunt sem varius.

![gene-delta-cai](./figures/gene-delta-cai.svg)

The ΔCAI values of genes within a same contig are averaged to obtain contig-level ΔCAI. Next, a kernel density estimation (KDE) is computed from the data in order to

![contig-delta-cai-kde](./figures/contig-delta-cai-kde.svg)

Pellentesque scelerisque nunc ligula, vitae congue tellus posuere accumsan. Nunc ac nibh sed ex semper tempus. Aliquam faucibus a magna sit amet hendrerit. Cras ac vestibulum mauris. In efficitur quam pretium massa aliquam maximus. Suspendisse at ornare augue, tincidunt gravida ante. Etiam sed dictum erat. Phasellus convallis ante ut turpis placerat convallis. Phasellus turpis lacus, molestie sed orci eu, sollicitudin tristique magna. Vivamus consequat eu nulla a consectetur. Duis in lacus bibendum, facilisis arcu eu, rhoncus lacus.

![contig-delta-cai-scores](./figures/contig-delta-cai-scores.svg)

Duis nisl ex, tempus ultrices mauris consectetur, convallis sagittis tellus. Nam nec lobortis dui. Curabitur a risus blandit erat luctus tristique a eget felis. In vestibulum, eros non auctor ultricies, ex leo eleifend lacus, at fermentum dui ligula id diam. Vestibulum arcu ex, finibus quis consequat fermentum, malesuada eu mi. In hendrerit ipsum vestibulum nibh elementum, suscipit faucibus augue viverra. Proin non condimentum orci. Ut sed metus ante. Nam pharetra interdum commodo.

::: tip Alternative codon usage bias encodings
The RSCU is not the only method to summarize the codon usage bias of an organism. [Hughes AL and Langley KJ](https://pubmed.ncbi.nlm.nih.gov/17000138/) proposed an alternative representation that encodes codon usage bias into five variables and this approach was succesfully applied to metagenomic binning by [Yu G et al](https://pubmed.ncbi.nlm.nih.gov/29947757/).
:::