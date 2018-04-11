As a R Script
=============

All operations can be conducted simply running the proper R script as
follows::

    Rscript reader.R [arguments]

To show a synthetic help, run the abouve command without any argument.

The arguments are a numerical code that identifies a single function
followed by it's arguments. To a more detailed description of the
functions, check the here `Functions <Functions>`_.

The arguments are the same for the python and R version

As a Python Script
==================

All operations can be conducted simply running the proper python
script as follows::
  
  python3 rpy_reader.py [arguments]

The arguments are a numerical code that identifies a single function
followed by it's arguments. To a more detailed description of the
functions, check the here `Functions <Functions>`_.

The arguments are the same for the python and R version

As a Shiny App
==============

A user-friendly `shiny <https://shiny.rstudio.com/>`_ interface will
be created using your computer as a local server. All the funtions
will be available in a GUI in your browser.

To use the interface, open a R interactive session::

  R

And source the environment.R file ::
  
  source("environment.R")

This file import all the necessaries libraries and data and creates a
shiny interface using your computer as server.

Warning: if you run the scripts in this way, consider that the server
is your own computer: in case of big data make shure you hane enough
memory.
