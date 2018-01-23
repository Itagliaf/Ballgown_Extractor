[![License: GPL v3](https://img.shields.io/badge/License-GPL%20v3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)

# Ballgown Extractor

Ballgown extractor is a R script able to manage transcriptome data
obtained by a [pipeline](https://www.ncbi.nlm.nih.gov/pubmed/27560171)
described by Petrea and collegues. This pipeline take advantage of
[Hisat](https://ccb.jhu.edu/software/hisat2/index.shtml),
[Stringtie](https://ccb.jhu.edu/software/stringtie/) and the R package
[Ballgown](http://bioconductor.org/packages/release/bioc/html/ballgown.html).

In particular, these scripts will allow you to extract information
about single genes, tissues and gene features.

# Dependencies

This script is written in [R](https://www.r-project.org/). It was
tested using the version 3.4.1 "Single Candle" but all distribution
supporting the packages required should be fine.


The following packages are needed:

- [ballgown](http://bioconductor.org/packages/release/bioc/html/ballgown.html)
- [dplyr](https://cran.r-project.org/web/packages/dplyr/README.html)
- [genefilter](http://bioconductor.org/packages/release/bioc/html/genefilter.html)
- [ggplot2](http://ggplot2.org/)

# Usage

This tool must be run from the command line using Rscript.

```bash
Rscript reader.R [arguments]
```
To show a synthetic help, run the above command without any argument.

To access the function, one should use the first positional argument as follows:

- 1 "Plotter": plots the FPKM values of one or more genes considered.
Arguments:
  - file containing a list of genes (separator=\n)

```bash
Rscript reader.R 1 input_genes.txt
```

- 2 "SearchByTissue": Searches a tissue and write  a file containing the FPKM values of all genes. If a third argument with the name of gene is given, the FPKM of this gene inside the chosen tissue. Arguments

  - tissue name;
  - gene name (optional);

```bash
#FPKM of all genes inside the tissue "liver"
Rscript reader.R 2 liver
#FPKM of the gene "gene1" in the tissue "liver"
Rscript reader.R 2 liver gene1
```

- 3 "SearchByGene": Searches the FPKM values of a single gene in all tissues. Arguments
  - Gene name to be analyzed

```bash
Rscript reader.R 3 gene1
```

- 4 "SearchByGeneFeature": allows to analyze the expression of gene features such as introns or exons. Arguments:
  - Gene name to be analyzed
  - Gene feature to be analyzed (Exon or Intron)

```bash
#results about the exons of gene1
Rscript reader.R 4 gene1 exon
#results about the introns of gene1
Rscript reader.R 4 gene1 intron
```

- 5 "SearchByDifferentialFold": Given 2 tissues, the differential fold (linear and log2) is calculated. To avoid "Inf" values, to all FPKM values was added 0.000001.
If a fourth argument is given reporting a gene name, the function will consider only the entries containing only that gene. Arguments
  - First tissue to be considered;
  - Second tissue to be considered;
  - Gene name (optional);

```bash
#Differential fold between liver and muscle
Rscript reader.R 5 liver muscle
#Differential fold between liver and muscle considering only the gene "gene1"
Rscript reader.R 5 liver muscle gene1
```
