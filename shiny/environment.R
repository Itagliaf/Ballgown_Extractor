library(shiny)
library(ggplot2)
library(ballgown)
library(Cairo)
options(shiny.usecairo=T)

##=== Loading functions ====
source("definitions.R")

##==== Loading and formatting data ====
load("data.fil.RData")
transcripts<-data.fil@expr$trans
head(transcripts)
phenodata<-read.csv("phenodata.csv","header"=TRUE)
colnames(transcripts)<-NameFormatter(transcripts,phenodata)
##==== END ====

runApp("app.R")
