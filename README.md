[![License: GPL v3](https://img.shields.io/badge/License-GPL%20v3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
|docs|



# Ballgown Extractor

Ballgown extractor is an R script able to manage transcriptome data
obtained by a [pipeline](https://www.ncbi.nlm.nih.gov/pubmed/27560171)
described by Petrea and collegues. This pipeline takes advantage of
[Hisat2](https://ccb.jhu.edu/software/hisat2/index.shtml),
[Stringtie](https://ccb.jhu.edu/software/stringtie/) and the R package
[Ballgown](http://bioconductor.org/packages/release/bioc/html/ballgown.html).

In particular, these scripts will allow you to extract information
about single genes, tissues and gene features.

The algorithmic part is contained in "definitons.R" script. This script 
can be run in 3 different ways:

- 1 using reader.R: running reader.R as script with the appropriate arguments
```
Rscript reader.R <function_code> [arguments]
```

- 2 using rpy_reader.py: this python scrypts uses [rpy2](https://pypi.python.org/pypi/rpy2)
to run the functions defined in definitions.R 

```
python3 rpy_reader.py <function_code> [aruments]
```


- 3 using App.R: this R scipt scrypts uses [shiny](https://www.rstudio.com/products/shiny/)
and opens a local shiny app.

```
Rscript App.R
```

# Dependencies

This script is written in [R](https://www.r-project.org/) and [pyhton 3](https://www.python.org/downloads/release/python-363/).
It wastested using the R version 3.4.1 "Single Candle" and the python 3.5.2.

The following R packages are required:

- [ballgown](http://bioconductor.org/packages/release/bioc/html/ballgown.html)
- [dplyr](https://cran.r-project.org/web/packages/dplyr/README.html)
- [genefilter](http://bioconductor.org/packages/release/bioc/html/genefilter.html)
- [ggplot2](http://ggplot2.org/)
- [igraph](http://igraph.org/r/)

The following python packages are required:

- [rpy2](https://pypi.python.org/pypi/rpy2)

All the dependecies are installed automatically by the script install_deps.sh


# Usage

## Command line
The "Launching scripts" must bu run from command line as follows

```
Rscript reader.R [arguments]

python3 rpy_reader.py [arguments]
```

To show a synthetic help, run the above command without any argument.

The arguments are a numerical code that identifies a single function and
their arguments.

The numeric code list to invoke the functions is reported below:

"Choose a function as follow:"

- 1 = Plotter Function

- 2 = Search by Tissue

- 3 = Search by Gene

- 4 = Search by gene feature

- 5 = Search by Differential Fold

- 6 = Search Coexpression Module

- 7 = Create Network Graph

- 8 = Differential Fold for a Gene in All Tissues

- 99 = Import Data From a Ballgonw Object

## Shiny app:

After running "App.R", a new browser window will open. All the funcions 
reported here are present as tabs inside the app.

# Functions

- 1 "Plotter": plots the FPKM values of all transcripts of a given gene
Arguments:
  - gene symbol wanted

```bash
Rscript reader.R 1 GENE1
```

- 2 "SearchByTissue": Searches a tissue and write a file containing the FPKM values of all genes.
If a second argument with the name of gene is given, the FPKM of this gene inside the chosen tissue.
Arguments:

  - tissue name;
  - gene name (optional);

```bash
#FPKM of all genes inside the tissue "liver"
Rscript reader.R 2 liver
#FPKM of the gene "gene1" in the tissue "liver"
Rscript reader.R 2 liver GENE1
```

- 3 "SearchByGene": Searches the FPKM values of a single gene in all tissues.
Arguments:
  - Gene name to be analyzed

```bash
Rscript reader.R 3 GENE1
```

- 4 "SearchByGeneFeature": allows to analyze the expression of gene features such as introns or exons. Arguments:
  - Gene symbol to be analyzed
  - Gene feature to be analyzed (Exon or Intron)

```bash
#results about the exons of GENE1
Rscript reader.R 4 GENE1 exon
#results about the introns of GENE1
Rscript reader.R 4 GENE1 intron
```

- 5 "SearchByDifferentialFold": Given 2 tissues, the differential fold (linear and log2) is calculated. To avoid "Inf" values, to all FPKM values was added 0.000001.
If a fourth argument is given reporting a gene name, the function will consider only the entries containing only that gene. Arguments
  - First tissue to be considered;
  - Second tissue to be considered;
  - Gene Symbol (optional);

```bash
#Differential fold between liver and muscle
Rscript reader.R 5 liver muscle
#Differential fold between liver and muscle considering only the gene "gene1"
Rscript reader.R 5 liver muscle GENE1
```

- 6 "SearchCoexpressionModule": Given a gene symbol, the coexpression module, found with [WGCNA](https://labs.genetics.ucla.edu/horvath/CoexpressionNetwork/Rpackages/WGCNA/)
is returned.
Arguments:

  - Gene Symbol;
  - File containing the coexpression module found with WGCNA;
  
```bash
Rscript reader.R 6 GENE1 Module_File
```

- 7 "Network": Given a gene symbol, from the coexpression module, found with [WGCNA](https://labs.genetics.ucla.edu/horvath/CoexpressionNetwork/Rpackages/WGCNA/)
a network graph containing only n nodes (genes) over a certain correlation value is produced
Arguments:

  - Gene Name or Symbol is used to identify the gene. This option will be deleted in a future version, just use "symbol";
  - Gene Symbol;
  - File containing the coexpression module found with WGCNA;
  - Correlation threashold value to show a connection between 2 genes
  - Number of nodes-genes to be shown
  
```bash
Rscript reader.R 7 symbol GENE1 Module_File 0.8 25
```


- 8 "Differential Fold for a Gene in All Tissues": Given a gene symbol, the differential fold between all genes is calculated considering
the oviduct as a refernce.
NB: up to now "oviduct" is hardcoded, this issue will be resolved soon

  - Gene Symbol;
  
```bash
Rscript reader.R 8 GENE1
```
