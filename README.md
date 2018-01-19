# Ballgown Extractor

Ballgown extractor is a R script able to manage transcriptome data obtained by a [pipeline](https://www.ncbi.nlm.nih.gov/pubmed/27560171) described by Petrea and collegues. This pipeline take advantage of [Hisat](https://ccb.jhu.edu/software/hisat2/index.shtml), [Stringtie](https://ccb.jhu.edu/software/stringtie/) and the R package [Ballgown](http://bioconductor.org/packages/release/bioc/html/ballgown.html).

In particular, these scripts will allow you to extract information about single genes, tissues and gene features.

# Dependencies

The following packages are needed:

- [ballgown](http://bioconductor.org/packages/release/bioc/html/ballgown.html)
- [dplyr](https://cran.r-project.org/web/packages/dplyr/README.html)
- [genefilter](http://bioconductor.org/packages/release/bioc/html/genefilter.html)
- [ggplot2](http://ggplot2.org/)