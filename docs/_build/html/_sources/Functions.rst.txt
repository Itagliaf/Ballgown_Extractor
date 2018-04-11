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
  (`Ballgown documentation <https://bioconductor.org/packages/release/bioc/manuals/ballgown/man/ballgown.pdf>`_).

* Folder Containing The Data: folder containing the results from the
  `pipeline <https://www.ncbi.nlm.nih.gov/pubmed/27560171>`_ described
  by Petrea and collegues.

* Common part of samples name: a common root or part of samples name
  needed by ballgown to extract properl all the data needed.

**This function is not available for the shiny app: this step must be executed via command line**
  
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

* tissue: the tissue name to be analyzed.

* gene: [optional] the gene symbol (name) to be analyzed.


-------------------------
Search By Gene (code=3)
-------------------------

Searches the FPKM values of a single gene in all tissues.

Syntax:
-------

::
   
   Rscript reader.R 3 <gene>
   python3 rpy_reader.py 3 <gene> 

* gene: the gene symbol (name) to be analyzed.

-------------------------------
Search By Gene Feature (code=4)
-------------------------------

Allows to analyze the expression of gene feature (exons and introns).

Syntax:
-------

::
   
   Rscript reader.R 4 <gene> <feature>
   python3 rpy_reader.py 4 <gene> <feature>

* gene: the gene symbol (name) to be analyzed.

* feature: the gene feature to analyzed. Possible values: intron, exon.

------------------------------------
Search By Differential Fold (code=5)
------------------------------------

Given 2 tissues, the differential fold (linear and log2) is
calculated. To avoid "Inf" values, to all FPKM values was added
0.000001. If a fourth argument is given reporting a gene name, the
function will consider only the entries containing only that
gene.

Syntax:
-------

::
   
   Rscript reader.R 5 <tissue 1> <tissue 2> [gene]
   python3 rpy_reader.py 5 <tissue 1> <tissue 2> [gene]

* tissue 1: First tissue to be considered.

* tissue 2: Second tissue to be considered.

* gene: [optional] the gene symbol (name) to be analyzed.

------------------------------------
Search By Differential Fold (code=6)
------------------------------------

Given a gene symbol, the coexpression module, found with WGCNA is
returned.

Syntax:
-------

::
   
   Rscript reader.R 6 <gene> <ModuleFile>
   python3 rpy_reader.py 6 <gene> <ModuleFile>

* gene: the gene symbol (name) to be analyzed.
* File containing the coexpression module found with `WGCNA
  <https://labs.genetics.ucla.edu/horvath/CoexpressionNetwork/Rpackages/WGCNA/>`_.

  **AGGIUNGERE LO SCRIPT PER LA COESPRESSIONE NELLA REPO**


------------------------------------
Search By Differential Fold (code=6)
------------------------------------

Given a gene symbol, from the coexpression module, found with `WGCNA
<https://labs.genetics.ucla.edu/horvath/CoexpressionNetwork/Rpackages/WGCNA/>`_
a network graph containing only n nodes (genes) over a certain
correlation value is produced Arguments:

Syntax:
-------

::
   
   Rscript reader.R 7 <symbol> <gene> <ModuleFile> <corr> <nodes> 
   python3 rpy_reader.py 7 <symbol> <gene> <ModuleFile> <corr> <nodes> 

* symbol: a legacy argument. Use symbol as first argument (to removed
  in the next versions).
* gene: the gene symbol (name) to be analyzed.
* File containing the coexpression module found with `WGCNA
  <https://labs.genetics.ucla.edu/horvath/CoexpressionNetwork/Rpackages/WGCNA/>`_.
* Correlation threashold value to show a connection between 2 genes.
* Number of nodes-genes to be shown.

  **AGGIUNGERE LO SCRIPT PER LA COESPRESSIONE NELLA REPO**

  
