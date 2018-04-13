[![License: GPL v3](https://img.shields.io/badge/License-GPL%20v3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)

[![Documentation Status](https://readthedocs.org/projects/ballgown-extractor/badge/?version=master)](http://ballgown-extractor.readthedocs.io/en/latest/?badge=master)

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

The full documentation is available here: [documentation](http://ballgown-extractor.readthedocs.io/en/master/)

