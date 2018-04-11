=========
Functions
=========
--------------------------------------------
Import Data From Ballgown Object (code = 99)
--------------------------------------------

This function allows to create all the proper files to be used in the
other functions taking advantage of the `Ballgown
<http://bioconductor.org/packages/release/bioc/html/ballgown.html>`_
package (ballgown function). To gain a better understanding of this
function, check the Ballgown documentation
(`<https://bioconductor.org/packages/release/bioc/manuals/ballgown/man/ballgown.pdf>`_)

Syntax:
-------
::
   
  Rscript reader.R 99 <Mapping file> <Folder Containing The Data> <Common part of samples name>
  python rpy_reader.py 99 <Mapping file> <Folder Containing The Data> <Common part of samples name>

Arguments:

* Mapping File: data frame with rows corresponding to samples and
  columns corresponding to phenotipic variables. Up to now, this file
  was structured as follows:

  +-------+-----------+
  |  ids  | tissue    |
  +=======+===========+
  | id 1  | name 1    |
  +-------+-----------+
  | id 2  | name 2    |
  +-------+-----------+

  Reporting only the tissue name. If a file presents more columns,
  there should be no problem because they will be ignored. This file
  will be used as the argument pData in the "ballgown" function
  (`Ballgown documentation <https://bioconductor.org/packages/release/bioc/manuals/ballgown/man/ballgown.pdf>`_);

* Folder Containing The Data: folder containing the results from the
  `pipeline <https://www.ncbi.nlm.nih.gov/pubmed/27560171>`_ described
  by Petrea and collegues;

* Common part of samples name: a common root or part of samples name
  needed by ballgown to extract properl all the data needed.


----------------
Plotter (code=1)
----------------

Plots the FPKM values of all transcripts of a given gene. This
function takes advantage of the "plotTranscript" function in the
ballgown package. For further informations, check the `Ballgown
documentation
<https://bioconductor.org/packages/release/bioc/manuals/ballgown/man/ballgown.pdf>`_

Syntax:
-------

::
   
   Rscript reader.R 1 <gene>
   python3 rpy_reader.py 1 <gene>

* gene: the gene symbol (name) to be analyzed.

-------------------------
Search By Tissue (code=2)
-------------------------

Searches the a given tissue and writes a file containing the FPKM
vaules of all genes. If a second argument with the name of a gene is
given, the FPKM of this gene insider the chosen Tissue.

Syntax:
-------

::
   
   Rscript reader.R 2 <tissue> [gene]
   python3 rpy_reader.py 1 <tissue> [gene]

* tissue: the tissue name to be analyzed

* gene: [optional] the gene symbol (name) to be analyzed.
