(window.webpackJsonp=window.webpackJsonp||[]).push([[4],{348:function(e,t,n){e.exports=n.p+"assets/img/tnf-embedding.094edc5e.svg"},359:function(e,t,n){"use strict";n.r(t);var i=n(42),o=Object(i.a)({},(function(){var e=this,t=e.$createElement,i=e._self._c||t;return i("ContentSlotsDistributor",{attrs:{"slot-key":e.$parent.slotKey}},[i("h1",{attrs:{id:"composition-module"}},[i("a",{staticClass:"header-anchor",attrs:{href:"#composition-module"}},[e._v("#")]),e._v(" Composition module")]),e._v(" "),i("p",[e._v("Contigs assembled from reads derived from the same genome or from genomes of closely related organisms tend to display a similar sequence composition profile, represented as k-mer frequencies. Metagenomic binners exploit this property to cluster contigs into putative genomic bins using 4-mer frequencies, or tetranucleotide frequencies (TNFs).")]),e._v(" "),i("div",{staticClass:"custom-block tip"},[i("p",{staticClass:"custom-block-title"},[e._v("Genomic islands and plasmids")]),e._v(" "),i("p",[e._v("Mobile genetic elements such as genomic islands and plasmids usually have a 4-mer composition that is distinct from the majority of the genome, making them problematic for binning algorithms.")])]),e._v(" "),i("p",[e._v("Whether or not two given contigs will be clustered into the same genomic bin does not depend exclusively on their TNF profiles. In most modern clustering algorithms local relationships are influenced by other data points, meaning that a given pair of contigs may end up in the same bin or not, depending on the full set contigs that is being clustered. Moreover, most binners also use sequencing coverage information in addition to TNF data to cluster contigs, which may lead to genomic bins that encompass contigs with distinct TNF profiles.")]),e._v(" "),i("p",[e._v('MAGpurify2 processes each genomic bin individually and finds potential contaminants with respect to TNF profile by identifying contigs that fall outside of the "core cluster" within the bin.')]),e._v(" "),i("h2",{attrs:{id:"how-it-works"}},[i("a",{staticClass:"header-anchor",attrs:{href:"#how-it-works"}},[e._v("#")]),e._v(" How it works")]),e._v(" "),i("p",[e._v('To identify putative contaminants within a genomic bin, MAGpurify2: (1) computes the TNF profile of each contig, (2) embbeds data points into a low-dimentional space using a non-linear transformation, and (3) finds the "core cluster" and computes each contig score.')]),e._v(" "),i("p",[i("img",{attrs:{src:n(348),alt:"tnf-embedding"}})]),e._v(" "),i("p",[e._v("The four DNA bases (A, T, C and G) can produce 256 distinct 4-mers, however, in a strand-independent analysis, reverse complement k-mers (eg.: "),i("code",[e._v("TTAC")]),e._v(" and "),i("code",[e._v("GTAA")]),e._v(") are redundant and should be counted as a single entity (a canonical k-mer) in order to reduce memory usage and data variance. Thus, MAGpurify2 counts the 136 canonical 4-mers for each contig within the bin and computes their relative frequencies as maximum-likelihood estimations of the underlying TNF profile of the sequence.")]),e._v(" "),i("p",[e._v("::: info TNF profile of short contigs\nShort contigs contain a reduced number of 4-mers and therefore provide less reliable estimations of the underlying genomic TNF profile than longer contigs. This is one of the reasons why most binners filter out contigs shorter than a set threshold (usually around 2,000 bp) before clustering. Currently, MAGpurify2 does not take into account the length-dependent statistical uncertainty of the TNF estimation when identifying putative contaminants.\n:::")]),e._v(" "),i("p",[e._v("The high dimensional 4-mer frequency data is then non-linearly projected into a three dimensional space using the "),i("a",{attrs:{href:"https://umap-learn.readthedocs.io/en/latest/",target:"_blank",rel:"noopener noreferrer"}},[e._v("UMAP"),i("OutboundLink")],1),e._v(" algorithm, which will bring similar data points together and distance contigs with distinct TNF profiles. Next, "),i("a",{attrs:{href:"https://hdbscan.readthedocs.io/en/latest/",target:"_blank",rel:"noopener noreferrer"}},[e._v("hdbscan"),i("OutboundLink")],1),e._v(' is used to identify clusters within the UMAP embedding and, if at least one cluster is found, compute the membership of each contig to the "core cluster". The "core cluster" is defined as the cluster that emcompassess the largest assembled fraction, that is, the sum of the lengths of all the contigs within the cluster.')]),e._v(" "),i("p",[e._v('As UMAP is a non-deterministic algorithm, MAGpurify2 executes multiple iterations of the dimension reduction and clustering steps. The final contig score correspond to the average of its membership to the "core cluster" across the iterations.')])])}),[],!1,null,null,null);t.default=o.exports}}]);