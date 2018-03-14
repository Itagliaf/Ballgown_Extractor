print(">>> Ballgown Extractor Installer <<<")

print(version)

print("Installing dplyr")
install.packages("dplyr",repos="https://cran.stat.unipd.it/")

print("Sourcing Bioconductor")
source("http://bioconductor.org/biocLite.R")

print("Installing genefilter")
biocLite("genefilter")

print("Installing ballgown")
biocLite("ballgown")

print("Installing ggplot2")
install.packages("ggplot2",,repos="https://cran.stat.unipd.it/")

print("Installing igraph")
install.packages("igraph",,repos="https://cran.stat.unipd.it/")



