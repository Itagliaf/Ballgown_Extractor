******************
Ballgown Extractor
******************

Ballgown extractor is a collection of scripts, written in R and
python, able to manage transcriptome data obtained by a `pipeline
<https://www.ncbi.nlm.nih.gov/pubmed/27560171>`_ described by Petrea
and collegues, which takes advantages of `Hisat2
<https://ccb.jhu.edu/software/hisat2/index.shtml>`_, `Stringtie
<https://ccb.jhu.edu/software/stringtie/>`_ and the R package
`Ballgown
<http://bioconductor.org/packages/release/bioc/html/ballgown.html>`_.
In particular these scripts will allow you to extract information
about single genes, tissues and gene features.

In this project there are different scripts that produces the same
result:


* A simple R script (reader.R);
  
* A python script (rpy_reader.py) that takes advantage of the `rpy2 <https://pypi.python.org/pypi/rpy2>`_
  library;
  
* A R script spawning a shiny app (App.R), that takes advantage of `shiny <https://www.rstudio.com/products/shiny/>`_

Installation
============

.. toctree::
   :maxdepth: 3

   Installation
   
Usage
=====
 
.. toctree::
   :maxdepth: 3

   Usage

Functions
=========

.. toctree::
   :maxdepth: 3

   Functions
   

Another Simple Header
=====================

Here is some text explaining some complicated stuff ::

  print 'hello'
  >>hello



..
   Indices and tables
   ==================

   * :ref:`genindex`
   * :ref:`modindex`
   * :ref:`search`
