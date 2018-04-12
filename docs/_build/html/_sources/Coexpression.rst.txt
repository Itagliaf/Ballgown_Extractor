=============
Coexpression
=============

To run the functions "Search Coexpression Group" (code=6) and "Create
Network Graph" (code=7), a file containing the coexpression module
must be produced.

To create this file, is possible to use the "Coexpression.R" script,
which exetutes the `WGCNA part 1 tutorial
<https://labs.genetics.ucla.edu/horvath/CoexpressionNetwork/Rpackages/WGCNA/Tutorials/>`_
finding the coexpression group of each gene.

This step is not fully automated because a number of steps must be
supervised by the user to complete in the best way the analisys.

Moreover there are a few variables that must be modified accoringly to
your system and data (eg: number of cores). Up to now, no
"interactive" way to modify these variables are implemented: the user
must change the values of such variables in the code. This variables
are highlighted in the code.

To run the script, simply type::
  
  Rscript coexpression.R

*It's advisable to use this script as a roadmap and conduct the
analyses as reported in the official documentation of* `WGCNA
<https://labs.genetics.ucla.edu/horvath/CoexpressionNetwork/Rpackages/WGCNA/Tutorials/>`_
