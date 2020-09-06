# Coverage module

When read mappings are supplied by the user, MAGpurify2 processes the alignments to ensure that the coverage estimates are as accurate as possible. This is achieved by employing the following filters:

- Both edges of the contig are excluded from coverage computation (75 bp from each). These regions are unreliable for coverage estimation as read mappers have troubles partially aligning reads that would extend off edges.

- Reads that map to the contig with an alignment identity lower than a set threshold (97%, by default) are discarded. As species are typically defined as clusters of genomes that display an average nucleotide identity of at least 95%, this filter removes reads that are likely to be originated from a distinct species.

- Variations in sequence conservation across species, the presence of repeat regions, sequencing biases, and sampling biases can all cause different regions within the same contig to have very distinct read coverages. To prevent large within-contig coverage variations from shifting the estimated contig coverage, MAGpurify2 trims the lower and upper tails of the within-coverage distribution (5% of each, by default) before taking the mean.

![coverage-trimming](./figures/coverage-trimming.svg)

[Measurement of bacterial replication rates in microbial communities](https://www.nature.com/articles/nbt.3704)

Lorem ipsum dolor sit amet, consectetur adipiscing elit. Nam id nulla ac velit elementum tempus non eget urna. Etiam placerat leo ac risus semper feugiat ut id nulla. Fusce venenatis magna non feugiat convallis.

$a>n <=> a>=n+1$, if $a, n\in\Z$

$$
\Gamma(z) = \int_0^\infty t^{z-1}e^{-t}dt\,.
$$

Lorem ipsum dolor sit amet, consectetur adipiscing elit. Nam id nulla ac velit elementum tempus non eget urna. Etiam placerat leo ac risus semper feugiat ut id nulla. Fusce venenatis magna non feugiat convallis. Maecenas hendrerit orci quis elit pretium aliquet. Nulla lorem lectus, tempus ut sollicitudin nec, ultricies eu ex. Morbi at fringilla nunc, eget facilisis enim. Nullam a iaculis massa. Sed faucibus leo sed consectetur convallis. Aenean diam neque, imperdiet et mauris nec, mollis dignissim sapien. Vestibulum imperdiet magna a erat pellentesque, porta laoreet enim dapibus. Quisque vitae egestas enim. Donec in lacus volutpat dolor malesuada aliquet. Maecenas eget blandit ipsum. Vivamus eget posuere tortor.