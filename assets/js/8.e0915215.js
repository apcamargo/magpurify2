(window.webpackJsonp=window.webpackJsonp||[]).push([[8],{365:function(t,a,s){t.exports=s.p+"assets/img/coverage-trimming.7a02fc40.svg"},377:function(t,a,s){"use strict";s.r(a);var e=s(17),i=Object(e.a)({},(function(){var t=this,a=t.$createElement,e=t._self._c||a;return e("ContentSlotsDistributor",{attrs:{"slot-key":t.$parent.slotKey}},[e("h1",{attrs:{id:"coverage-module"}},[e("a",{staticClass:"header-anchor",attrs:{href:"#coverage-module"}},[t._v("#")]),t._v(" Coverage module")]),t._v(" "),e("p",[t._v("Contigs originated from the same genome are physically linked, so it is expected that they have a similar abundances (quantified as the average sequencing coverage) within a sample and correlated relative abundances across multiple samples. MAGpurify2 identifies putative contaminants by findings contigs with outlier coverage values in relation to the bin's profile.")]),t._v(" "),e("div",{staticClass:"custom-block tip"},[e("p",{staticClass:"custom-block-title"},[t._v("Uneven coverage across the genome")]),t._v(" "),e("p",[t._v("Even though metagenome binners assume a approximately uniform sequencing coverage across the genomes, there are factors of both technical and biological nature that lead to uneven coverage of contigs originated from the same genome. Three major causes of nonuniform genome coverage are: (1) GC content sequencing bias, that cause a coverage reduction in GC-rich regions; (2) repeat genomic regions, that exhibit copy number-dependent coverage; and differences in replication rates. Regarding the latter, regions of the genome that are closer to the replication origin tend to have higher coverage than origin-distal sequences and the magnitude of this gradient depends on the genomic replication rate "),e("sup",{staticClass:"footnote-ref"},[e("a",{attrs:{href:"#fn1",id:"fnref1"}},[t._v("[1]")])]),t._v(".")])]),t._v(" "),e("h2",{attrs:{id:"computing-the-relative-contig-average-coverage"}},[e("a",{staticClass:"header-anchor",attrs:{href:"#computing-the-relative-contig-average-coverage"}},[t._v("#")]),t._v(" Computing the relative contig average coverage")]),t._v(" "),e("p",[t._v("When read mappings are supplied by the user, MAGpurify2 processes the alignments to ensure that the coverage estimates are as accurate as possible. This is achieved by employing the following filters:")]),t._v(" "),e("ul",[e("li",[e("p",[t._v("Reads that map to the edges of the contig (75 bp) are excluded from coverage computation. These regions are unreliable for coverage estimation as read mappers have trouble to partially align reads that would extend off edges.")])]),t._v(" "),e("li",[e("p",[t._v("Reads that map to the contig with an alignment identity lower than a specified threshold (97%, by default) are discarded. As species are typically defined as genome clusters with an average nucleotide identity of at least 95%, this filter removes reads that are likely to be originated from a distinct species.")])]),t._v(" "),e("li",[e("p",[t._v("Variations in sequence conservation across species, the presence of repeat regions, sequencing biases, and sampling biases can all cause different regions within the same contig to have very distinct read coverages (see box above). To prevent large within-contig coverage variations from shifting the estimated contig coverage, MAGpurify2 trims the lower and upper tails of the within-coverage distribution (5% of each, by default) before taking the mean.")])])]),t._v(" "),e("p",[e("img",{attrs:{src:s(365),alt:"coverage-trimming"}})]),t._v(" "),e("p",[t._v("After performing read filtering, MAGpurify computes the absolute average coverage ("),e("eq",[e("span",{staticClass:"katex"},[e("span",{staticClass:"katex-mathml"},[e("math",{attrs:{xmlns:"http://www.w3.org/1998/Math/MathML"}},[e("semantics",[e("mrow",[e("mi",[t._v("a")])],1),e("annotation",{attrs:{encoding:"application/x-tex"}},[t._v("a")])],1)],1)],1),e("span",{staticClass:"katex-html",attrs:{"aria-hidden":"true"}},[e("span",{staticClass:"base"},[e("span",{staticClass:"strut",staticStyle:{height:"0.43056em","vertical-align":"0em"}}),e("span",{staticClass:"mord mathdefault"},[t._v("a")])])])])]),t._v(") of each contig ("),e("eq",[e("span",{staticClass:"katex"},[e("span",{staticClass:"katex-mathml"},[e("math",{attrs:{xmlns:"http://www.w3.org/1998/Math/MathML"}},[e("semantics",[e("mrow",[e("mi",[t._v("i")])],1),e("annotation",{attrs:{encoding:"application/x-tex"}},[t._v("i")])],1)],1)],1),e("span",{staticClass:"katex-html",attrs:{"aria-hidden":"true"}},[e("span",{staticClass:"base"},[e("span",{staticClass:"strut",staticStyle:{height:"0.65952em","vertical-align":"0em"}}),e("span",{staticClass:"mord mathdefault"},[t._v("i")])])])])]),t._v(") in each sample ("),e("eq",[e("span",{staticClass:"katex"},[e("span",{staticClass:"katex-mathml"},[e("math",{attrs:{xmlns:"http://www.w3.org/1998/Math/MathML"}},[e("semantics",[e("mrow",[e("mi",[t._v("j")])],1),e("annotation",{attrs:{encoding:"application/x-tex"}},[t._v("j")])],1)],1)],1),e("span",{staticClass:"katex-html",attrs:{"aria-hidden":"true"}},[e("span",{staticClass:"base"},[e("span",{staticClass:"strut",staticStyle:{height:"0.85396em","vertical-align":"-0.19444em"}}),e("span",{staticClass:"mord mathdefault",staticStyle:{"margin-right":"0.05724em"}},[t._v("j")])])])])]),t._v("):")],1),t._v(" "),e("section",[e("eqn",[e("span",{staticClass:"katex-display"},[e("span",{staticClass:"katex"},[e("span",{staticClass:"katex-mathml"},[e("math",{attrs:{xmlns:"http://www.w3.org/1998/Math/MathML"}},[e("semantics",[e("mrow",[e("msub",[e("mi",[t._v("a")]),e("mrow",[e("mi",[t._v("i")]),e("mo",{attrs:{separator:"true"}},[t._v(",")]),e("mi",[t._v("j")])],1)],1),e("mo",[t._v("=")]),e("mfrac",[e("msub",[e("mrow",[e("mi",[t._v("R")]),e("mi",[t._v("e")]),e("mi",[t._v("a")]),e("mi",[t._v("d")]),e("mi",[t._v("c")]),e("mi",[t._v("o")]),e("mi",[t._v("u")]),e("mi",[t._v("n")]),e("mi",[t._v("t")])],1),e("mi",[t._v("i")])],1),e("mrow",[e("msub",[e("mrow",[e("mi",[t._v("L")]),e("mi",[t._v("e")]),e("mi",[t._v("n")]),e("mi",[t._v("g")]),e("mi",[t._v("t")]),e("mi",[t._v("h")])],1),e("mi",[t._v("i")])],1),e("mo",[t._v("−")]),e("mn",[t._v("2")]),e("mo",[t._v("∗")]),e("mn",[t._v("75")])],1)],1)],1),e("annotation",{attrs:{encoding:"application/x-tex"}},[t._v("\na_{i,j} = \\frac{\\mathit{Read count}_i}{\\mathit{Length}_i - 2 * 75}\n")])],1)],1)],1),e("span",{staticClass:"katex-html",attrs:{"aria-hidden":"true"}},[e("span",{staticClass:"base"},[e("span",{staticClass:"strut",staticStyle:{height:"0.716668em","vertical-align":"-0.286108em"}}),e("span",{staticClass:"mord"},[e("span",{staticClass:"mord mathdefault"},[t._v("a")]),e("span",{staticClass:"msupsub"},[e("span",{staticClass:"vlist-t vlist-t2"},[e("span",{staticClass:"vlist-r"},[e("span",{staticClass:"vlist",staticStyle:{height:"0.311664em"}},[e("span",{staticStyle:{top:"-2.5500000000000003em","margin-left":"0em","margin-right":"0.05em"}},[e("span",{staticClass:"pstrut",staticStyle:{height:"2.7em"}}),e("span",{staticClass:"sizing reset-size6 size3 mtight"},[e("span",{staticClass:"mord mtight"},[e("span",{staticClass:"mord mathdefault mtight"},[t._v("i")]),e("span",{staticClass:"mpunct mtight"},[t._v(",")]),e("span",{staticClass:"mord mathdefault mtight",staticStyle:{"margin-right":"0.05724em"}},[t._v("j")])])])])]),e("span",{staticClass:"vlist-s"},[t._v("​")])]),e("span",{staticClass:"vlist-r"},[e("span",{staticClass:"vlist",staticStyle:{height:"0.286108em"}},[e("span")])])])])]),e("span",{staticClass:"mspace",staticStyle:{"margin-right":"0.2777777777777778em"}}),e("span",{staticClass:"mrel"},[t._v("=")]),e("span",{staticClass:"mspace",staticStyle:{"margin-right":"0.2777777777777778em"}})]),e("span",{staticClass:"base"},[e("span",{staticClass:"strut",staticStyle:{height:"2.30158em","vertical-align":"-0.9301400000000001em"}}),e("span",{staticClass:"mord"},[e("span",{staticClass:"mopen nulldelimiter"}),e("span",{staticClass:"mfrac"},[e("span",{staticClass:"vlist-t vlist-t2"},[e("span",{staticClass:"vlist-r"},[e("span",{staticClass:"vlist",staticStyle:{height:"1.37144em"}},[e("span",{staticStyle:{top:"-2.3139999999999996em"}},[e("span",{staticClass:"pstrut",staticStyle:{height:"3em"}}),e("span",{staticClass:"mord"},[e("span",{staticClass:"mord"},[e("span",{staticClass:"mord"},[e("span",{staticClass:"mord mathit"},[t._v("L")]),e("span",{staticClass:"mord mathit"},[t._v("e")]),e("span",{staticClass:"mord mathit"},[t._v("n")]),e("span",{staticClass:"mord mathit"},[t._v("g")]),e("span",{staticClass:"mord mathit"},[t._v("t")]),e("span",{staticClass:"mord mathit"},[t._v("h")])]),e("span",{staticClass:"msupsub"},[e("span",{staticClass:"vlist-t vlist-t2"},[e("span",{staticClass:"vlist-r"},[e("span",{staticClass:"vlist",staticStyle:{height:"0.21752399999999997em"}},[e("span",{staticStyle:{top:"-2.4558600000000004em","margin-right":"0.05em"}},[e("span",{staticClass:"pstrut",staticStyle:{height:"2.7em"}}),e("span",{staticClass:"sizing reset-size6 size3 mtight"},[e("span",{staticClass:"mord mathdefault mtight"},[t._v("i")])])])]),e("span",{staticClass:"vlist-s"},[t._v("​")])]),e("span",{staticClass:"vlist-r"},[e("span",{staticClass:"vlist",staticStyle:{height:"0.24414em"}},[e("span")])])])])]),e("span",{staticClass:"mspace",staticStyle:{"margin-right":"0.2222222222222222em"}}),e("span",{staticClass:"mbin"},[t._v("−")]),e("span",{staticClass:"mspace",staticStyle:{"margin-right":"0.2222222222222222em"}}),e("span",{staticClass:"mord"},[t._v("2")]),e("span",{staticClass:"mspace",staticStyle:{"margin-right":"0.2222222222222222em"}}),e("span",{staticClass:"mbin"},[t._v("∗")]),e("span",{staticClass:"mspace",staticStyle:{"margin-right":"0.2222222222222222em"}}),e("span",{staticClass:"mord"},[t._v("7")]),e("span",{staticClass:"mord"},[t._v("5")])])]),e("span",{staticStyle:{top:"-3.23em"}},[e("span",{staticClass:"pstrut",staticStyle:{height:"3em"}}),e("span",{staticClass:"frac-line",staticStyle:{"border-bottom-width":"0.04em"}})]),e("span",{staticStyle:{top:"-3.677em"}},[e("span",{staticClass:"pstrut",staticStyle:{height:"3em"}}),e("span",{staticClass:"mord"},[e("span",{staticClass:"mord"},[e("span",{staticClass:"mord"},[e("span",{staticClass:"mord mathit"},[t._v("R")]),e("span",{staticClass:"mord mathit"},[t._v("e")]),e("span",{staticClass:"mord mathit"},[t._v("a")]),e("span",{staticClass:"mord mathit"},[t._v("d")]),e("span",{staticClass:"mord mathit"},[t._v("c")]),e("span",{staticClass:"mord mathit"},[t._v("o")]),e("span",{staticClass:"mord mathit"},[t._v("u")]),e("span",{staticClass:"mord mathit"},[t._v("n")]),e("span",{staticClass:"mord mathit"},[t._v("t")])]),e("span",{staticClass:"msupsub"},[e("span",{staticClass:"vlist-t vlist-t2"},[e("span",{staticClass:"vlist-r"},[e("span",{staticClass:"vlist",staticStyle:{height:"0.31166399999999994em"}},[e("span",{staticStyle:{top:"-2.5500000000000003em","margin-right":"0.05em"}},[e("span",{staticClass:"pstrut",staticStyle:{height:"2.7em"}}),e("span",{staticClass:"sizing reset-size6 size3 mtight"},[e("span",{staticClass:"mord mathdefault mtight"},[t._v("i")])])])]),e("span",{staticClass:"vlist-s"},[t._v("​")])]),e("span",{staticClass:"vlist-r"},[e("span",{staticClass:"vlist",staticStyle:{height:"0.15em"}},[e("span")])])])])])])])]),e("span",{staticClass:"vlist-s"},[t._v("​")])]),e("span",{staticClass:"vlist-r"},[e("span",{staticClass:"vlist",staticStyle:{height:"0.9301400000000001em"}},[e("span")])])])]),e("span",{staticClass:"mclose nulldelimiter"})])])])])])])],1),e("p",[t._v("To account for differences in sequence depth between different samples, the relative coverage ("),e("eq",[e("span",{staticClass:"katex"},[e("span",{staticClass:"katex-mathml"},[e("math",{attrs:{xmlns:"http://www.w3.org/1998/Math/MathML"}},[e("semantics",[e("mrow",[e("mi",[t._v("c")])],1),e("annotation",{attrs:{encoding:"application/x-tex"}},[t._v("c")])],1)],1)],1),e("span",{staticClass:"katex-html",attrs:{"aria-hidden":"true"}},[e("span",{staticClass:"base"},[e("span",{staticClass:"strut",staticStyle:{height:"0.43056em","vertical-align":"0em"}}),e("span",{staticClass:"mord mathdefault"},[t._v("c")])])])])]),t._v(") of each contig in each sample computed by dividing its absolute coverage by the sample's total coverage:")],1),t._v(" "),e("section",[e("eqn",[e("span",{staticClass:"katex-error",staticStyle:{color:"#cc0000"},attrs:{title:"ParseError: KaTeX parse error: Expected group after &#x27;\\frac&#x27; at position 26: … \\frac{a_{i,j}}_̲i}{\\sum_{i=1}^{…"}},[t._v("\nc_{i,j} = \\frac{a_{i,j}}_i}{\\sum_{i=1}^{\\mathit{No contigs}_j}\\:a_{i,j}}\n")])])],1),e("p",[t._v("Lorem ipsum dolor sit amet, consectetur adipiscing elit. Nam id nulla ac velit elementum tempus non eget urna. Etiam placerat leo ac risus semper feugiat ut id nulla. Fusce venenatis magna non feugiat convallis. Maecenas hendrerit orci quis elit pretium aliquet. Nulla lorem lectus, tempus ut sollicitudin nec, ultricies eu ex. Morbi at fringilla nunc, eget facilisis enim. Nullam a iaculis massa. Sed faucibus leo sed consectetur convallis. Aenean diam neque, imperdiet et mauris nec, mollis dignissim sapien. Vestibulum imperdiet magna a erat pellentesque, porta laoreet enim dapibus. Quisque vitae egestas enim. Donec in lacus volutpat dolor malesuada aliquet. Maecenas eget blandit ipsum. Vivamus eget posuere tortor.")]),t._v(" "),e("hr",{staticClass:"footnotes-sep"}),t._v(" "),e("section",{staticClass:"footnotes"},[e("ol",{staticClass:"footnotes-list"},[e("li",{staticClass:"footnote-item",attrs:{id:"fn1"}},[e("p",[t._v("Skovgaard, Ole, et al. "),e("a",{attrs:{href:"https://genome.cshlp.org/content/21/8/1388",target:"_blank",rel:"noopener noreferrer"}},[t._v('"Genome-wide detection of chromosomal rearrangements, indels, and mutations in circular chromosomes by short read sequencing."'),e("OutboundLink")],1),t._v(" "),e("em",[t._v("Genome Research")]),t._v(" (2011). "),e("a",{staticClass:"footnote-backref",attrs:{href:"#fnref1"}},[t._v("↩︎")])])])])])])}),[],!1,null,null,null);a.default=i.exports}}]);