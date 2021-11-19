---
title: "Simu_Arnaud_GQ"
output: html_document
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```

## R Markdown

This is an R Markdown document. Markdown is a simple formatting syntax for authoring HTML, PDF, and MS Word documents. For more details on using R Markdown see <http://rmarkdown.rstudio.com>.

When you click the **Knit** button a document will be generated that includes both content as well as the output of any embedded R code chunks within the document. You can embed an R code chunk like this:

The simulation program works with a certain number of functions that are described below before the main function that uses these fucntions.  

The GPmap function returns the genotypic value (mean phenotype) corresponding to a genotype.  

```{r}
GPmap <- function(genotype) {
  sum(genotype)
}
```
The get.phenotype function returns a phenotype value (1) from a Normal distribution (rnorm) of (mean) the genotypic value, and standard deviation (sd) the square root of the environmental variance as defined for the population.  

```{r}
get.phenotype <- function(genotype, var.env) {

	rnorm(1, mean=GPmap(genotype), sd=sqrt(var.env))
}
```
The init.individual function generates a random individual for the starting population. The genotype of the individual is defined as a matrix of 2 columns (2 alleles), the number of rows being equal to the number of loci. The value of each allele is drawn from a Normal distribution (rnorm) of (mean) 0 and (sd): the initial genetic variance divided by the number of alleles in the population (=2 times the number of loci). The individual is defined by its (genotype), its genotypic value (genot.value) which is the sum of the allelic values, its (phenotype), and its fitness. Note here that the fitness of the initial individuals is 1 whatever their phenotypes. Fitnesses are then updated in the simulations (see below). This is because we need all individuals before applying selection.  

```{r}
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
```
You can try to generate a single individual by typing:    
ind<-init.individual(1, 5, 1)  
Each individual is represented by a list of 4 elements containing the genotype, the genotypic value, the phenotype and the fitness.  
Access the genotype of that individual by typing:  
ind$genotype  
you can observe the values of the two alleles at the 5 loci.  
Now its genotypic value, by typing:  
ind$genot.value  
which is the sum of the genotypes (the 10 elements contained in ind$genot.value)  
and finally access its phenotype by:  
ind$phenotype  

The function init.population generates the initial population with as many individuals (init.individual) as in the population (pop.size) and returns those individuals as a list.   
```{r}
init.population <- function(pop.size, var.init, num.loci, var.env) {
	replicate(pop.size, init.individual(var.init, num.loci, var.env), simplify=FALSE)
}
```
Now generate your population by typing:  
pop<-init.population(100, 1, 5, 1)  
This population contains 100 individuals, each of which is a list of the 4 elements described above.  
The first individual can be accessed by  
pop[[1]] and its values pop[[1]]$genotype etc...  
The last is pop[[100]]  

The make.gamete function makes a haploid gamete out of an individual. Random drawing is done by the sample function which draws by chance either allele 1 or 2 (c(1,2)) as many times as the number of rows in the genotype of an individual (number of loci) with replacement.  
Recombination rate is 0.5 between loci (=free recombination), so that drawing allele 1 or 2 at a given locus does not depend on the preceding value drawn.  

```{r}
make.gamete <- function(indiv) {
	indiv$genotype[cbind(1:nrow(indiv$genotype), sample(c(1,2), nrow(indiv$genotype), replace=TRUE))]
}
```
To understand what the function does you can type:  
pop[[1]]$genotype to access the genotype of the first individual of the population  
and type:  
gam<-make.gamete(pop[[1]])  
gam is a random drawing at each locus of one of the two alleles of pop[[1]]$genotype.  
Visualize gam  

The make.offspring function makes an individual by binding the genotypes of two gametes coming from two individuals (mother and father).  

```{r}
make.offspring <- function(mother, father, var.env) {
	genotype <- cbind(make.gamete(mother), make.gamete(father))
	list(
		genotype  = genotype, 
		genot.value= GPmap(genotype),
		phenotype = get.phenotype(genotype, var.env),
		fitness   = 1
	)
}
```
For example type:  
offspr<-make.offspring(pop[[1]],pop[[2]], 1)  
Now verify which allele is coming from which parent at each locus.  

The update.fitness returns a new population object with updated fitnesses.  
The fitness function is a Gaussian distribution:  
Fitness = exp(- (phenotype - sel.optimum)^2 / (2 * sel.strength))  
The fitness of each individual depends on the distance between its phenotype and the optimum (sel.optimum) phenotype, as well as the strenghth of selection, a coefficient that determines the shape of the Gaussian distribution (a high value corresponds to a narrow Gausssian where a small proportion of individuals have a high fitness, while a small value corresponds to a wide Gaussian where a larger proportion of individuals have a higher fitness).  

```{r}
update.fitness <- function(population, sel.strength, sel.optimum) {
	lapply(population, function(indiv) {indiv$fitness <- exp(-(indiv$phenotype-sel.optimum)^2 /(2 * sel.strength)); indiv })
}
```
To look at the shape of the fitness function, you can generate phenotypes between -3 and 3 and visualize the function considering
sel.optimum<-0  
sel.strength<-1  
plot(seq(-3,3,0.01),exp(-(seq(-3,3,0.01) - sel.optimum)^2 /(2 * sel.strength)))  
now if you consider a greater value for sel.strength<-4, look what it changes:  
lines(seq(-3,3,0.01),exp(-(seq(-3,3,0.01) - sel.optimum)^2 /(2 * 4)))  
The function has a larger variance, and the selection is actually smoother than let's say for sel.strength<-0.1  
lines(seq(-3,3,0.01),exp(-(seq(-3,3,0.01) - sel.optimum)^2 /(2 * 0.1)))  

Now, for our simulated phenotypes, you can visualize the corresponding fitness function  
sel.optimum<-0  
sel.strength<-1  
phenotypes <- sapply(pop, "[[", "phenotype")  
and draw the fitness function:  
Fitness = exp(-(phenotypes - sel.optimum)^2 /(2 * sel.strength))  
plot(Fitness)  
You can modulate sel.strength, to make selection stronger (sel.strength=) or smoother (sel.strength=)  
Look at the shape of the fitness function.  

You can then type  
fitness<-update.fitness(pop, 1, 0)    
and check the update by looking at  
pop[[1]]$fitness  

The reproduction function returns the next generation, that is it generates individuals for the next generation population by sampling one mother and one father randomly from the previous population without replacement with a probability equals to their fitness. Mothers and fathers are generated during the make.offspring function (from gametes).  

```{r}
reproduction <- function(population, pop.size, var.env) {
	fitnesses <- sapply(population, "[[", "fitness")
	replicate(n=pop.size, 
		expr=make.offspring(
				mother=unlist(sample(population, 1, prob=fitnesses), recursive=FALSE), 
				father=unlist(sample(population, 1, prob=fitnesses), recursive=FALSE),
				var.env=var.env),
				simplify = FALSE)
}
```
if you type
newpop<-reproduction(pop,100,1)
you will see that a new population has been created.

The summary function computes summary statistics for the population. For phenotypes, genotypic values and fitness, it provides their mean and variance. 
```{r}
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
```
You can compare the summary of pop and new pop:
sumpop<-summary.population(pop)
sumnewpop<-summary.population(newpop)

We can plot the distribution of phenotypes, genotypic values and fitness of new pop by typing
hist(phenotypes)
hist(genot.val)
hist(fitnesses)


Below is the main function of the program: it initiates a population, and for each generation of simulations it creates a new pop and calculates the summary statitics.

```{r}
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
You can try running a first simulation with the following parameters values.
sim <- simulation(20, 100, 5, 1, 1, 1, 0)
You can plot several graphs to follow the evolution of phenotypes, genotypic values and fitness through time, as well as the selection differential.

