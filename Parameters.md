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

GPmap <- function(genotype) {
  sum(genotype)
}

get.phenotype <- function(genotype, var.env) {
  
  rnorm(1, mean=GPmap(genotype), sd=sqrt(var.env))
}

init.individual <- function(var.init, num.loci, var.env) {
  genotype <- matrix(
    rnorm(2*num.loci, mean=0, sd=sqrt(var.init/2/num.loci)), 
    ncol=2)
  list(
    genotype  = genotype, 
    genot.value= GPmap(genotype),
    phenotype = get.phenotype(genotype, var.env),
    fitness   = 1
  )
}

init.population <- function(pop.size, var.init, num.loci, var.env) {
  replicate(pop.size, init.individual(var.init, num.loci, var.env), simplify=FALSE)
}

make.gamete <- function(indiv) {
  indiv$genotype[cbind(1:nrow(indiv$genotype), sample(c(1,2), nrow(indiv$genotype), replace=TRUE))]
}

make.offspring <- function(mother, father, var.env) {
  genotype <- cbind(make.gamete(mother), make.gamete(father))
  list(
    genotype  = genotype, 
    genot.value= GPmap(genotype),
    phenotype = get.phenotype(genotype, var.env),
    fitness   = 1
  )
}

update.fitness <- function(population, sel.strength, sel.optimum) {
  lapply(population, function(indiv) {indiv$fitness <- exp(-(indiv$phenotype-sel.optimum)^2 /(2 * sel.strength)); indiv })
}

reproduction <- function(population, pop.size, var.env) {
  fitnesses <- sapply(population, "[[", "fitness")
  replicate(n=pop.size, 
            expr=make.offspring(
              mother=unlist(sample(population, 1, prob=fitnesses), recursive=FALSE), 
              father=unlist(sample(population, 1, prob=fitnesses), recursive=FALSE),
              var.env=var.env),
            simplify = FALSE)
}

summary.population <- function(population) {
  phenotypes <- sapply(population, "[[", "phenotype")
  genot.val  <- sapply(population, "[[", "genot.value")
  fitnesses  <- sapply(population, "[[", "fitness")
  data.frame(
    phen.mean = mean(phenotypes), 
    phen.var  = var (phenotypes),
    gen.mean  = mean(genot.val),
    gen.var   = var (genot.val),
    fit.mean  = mean(fitnesses),
    fit.var   = var (fitnesses),
    sel.diff  = mean(fitnesses*phenotypes)/mean(fitnesses) - mean(phenotypes)
  )
}

simulation <- function(generations, pop.size, num.loci, var.init, var.env, sel.strength, sel.optimum) {
  pop <- init.population(pop.size=pop.size, var.init=var.init, num.loci=num.loci, var.env=var.env)
  summ <- data.frame()
  for (gg in 1:generations) {
    pop <- update.fitness(pop, sel.strength, sel.optimum)
    summ <- rbind(summ, summary.population(pop))
    if (gg < generations)
      pop <- reproduction(pop, pop.size=pop.size, var.env=var.env)
  }
  summ
}

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

