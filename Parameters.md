---
title: "Parameters_Exploration_GQ"
output: html_document
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```

## R Markdown

This is an R Markdown document. Markdown is a simple formatting syntax for authoring HTML, PDF, and MS Word documents. For more details on using R Markdown see <http://rmarkdown.rstudio.com>.

When you click the **Knit** button a document will be generated that includes both content as well as the output of any embedded R code chunks within the document. You can embed an R code chunk like this:

The first part is exactly the same code as in the Readme and allows to perform the simulations. The second part allows exploring the parameters within the range indicated in the Schedule file.

```{r}

source("simu.R")

```
Now starts the parameters exploration.  
First, population size [5,1000].  
Graphs can be produced and are described below.   

```{r}

PopSize <- seq(5,1000, by=50)
for (i in 1:length(PopSize)) {
  current_PopSize <- PopSize[i]
  current_name <- paste0("PopSize_",current_PopSize)
  current_simu <- simulation(20,current_PopSize,10,1,1,1,3)
  assign(current_name,current_simu)
  rm(current_simu)
}

plot(NA,xlim=c(0,20),ylim=c(0,3.5),xlab="generations", ylab=("phen.mean"))
for (i in 1:length(PopSize)) 
  {
  ps<-get(paste("PopSize_",PopSize[i],sep=""))
  lines(ps$phen.mean,col=i)
  }
 
plot(NA,xlim=c(0,20),ylim=c(0,3.5),xlab="generations", ylab=("gen.mean"))
for (i in 1:length(PopSize)) 
{
  ps<-get(paste("PopSize_",PopSize[i],sep=""))
  lines(ps$gen.mean,col=i)
}

plot(NA,xlim=c(0,20),ylim=c(0,1),xlab="generations", ylab=("fit.mean"))
for (i in 1:length(PopSize)) 
{
  ps<-get(paste("PopSize_",PopSize[i],sep=""))
  lines(ps$fit.mean,col=i)
}

plot(NA,xlim=c(0,20),ylim=c(0,1),xlab="generations", ylab=("sel.diff"))
for (i in 1:length(PopSize)) 
{
  ps<-get(paste("PopSize_",PopSize[i],sep=""))
  lines(ps$sel.diff,col=i)
}

```

#Nb Locus
NbLocus <- seq(5,100, by=10)
for (i in 1:length(NbLocus)) {
  current_NbLocus <- NbLocus[i]
  current_name <- paste0("NbLocus_",current_NbLocus)
  current_simu <- simulation(20,20,current_NbLocus,1,1,1,3)
  assign(current_name,current_simu)
  rm(current_simu)
}

plot(NA,xlim=c(0,20),ylim=c(0,3.5),xlab="generations", ylab=("phen.mean"))
for (i in 1:length(NbLocus)) 
{
  ps<-get(paste("NbLocus_",NbLocus[i],sep=""))
  lines(ps$phen.mean,col=i)
}

plot(NA,xlim=c(0,20),ylim=c(0,3.5),xlab="generations", ylab=("gen.mean"))
for (i in 1:length(NbLocus)) 
{
  ps<-get(paste("NbLocus_",NbLocus[i],sep=""))
  lines(ps$gen.mean,col=i)
}

plot(NA,xlim=c(0,20),ylim=c(0,1),xlab="generations", ylab=("fit.mean"))
for (i in 1:length(NbLocus)) 
{
  ps<-get(paste("NbLocus_",NbLocus[i],sep=""))
  lines(ps$fit.mean,col=i)
}

