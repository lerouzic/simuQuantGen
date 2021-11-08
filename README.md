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

Here we define the parameters of the population (population size, number of (haploid) loci, initial genetic variance, environmental (residual) variance, strength of selection).
Truncation selection is defined as retaining the most extreme % of the individuals of the population to form the next generation. trunc.sel=1 indicates 100% of the individuals are retained (=no selection). Positive (resp. negative) values stand for selection towards higher (resp. lower) phenotypes. For example 0.1 means that you retain only 10% of the individuals (those with the highest values).
```{r}
global <- list(pop.size = 100, num.loci = 5, var.init = 1, var.env  = 1, trunc.sel= 1)

```
The GPmap function returns the genotypic value (mean phenotype) corresponding to a genotype
```{r}
GPmap <- function(genotype) {
  sum(genotype)
}
```
The get.phenotype function returns a phenotype value (1) from a Normal distribution (rnorm) of (mean) the genotypic value, and standard deviation (sd) the square root of the environmental variance as defined for the population.
```{r}
get.phenotype <- function(genotype, var.env) {

	rnorm(1, mean=GPmap(genotype), sd=sqrt(global$var.env))
}
```
The init.individual function generates a random individual for the starting population. The genotype of the individual is defined as a matrix of 2 columns (2 alleles), the number of rows being equal to the number of loci. The value of each allele is drawn from a Normal distribution (rnorm) of (mean) 0 and (sd): the initial genetic variance divided by the number of alleles in the population (=2 times the number of loci). The individual is defined by its (genotype), its genotypic value (genot.value) which is the sum of the allelic values, its (phenotype), and its fitness. Note here that the fitness of the initial individuals is 1 whatever their phenotypes. Fitnesses are then updated in the simulations (see below). This is because we need all individuals before applying truncation selection.
```{r}
 init.individual <- function(var.init, num.loci) {
	genotype <- matrix(
		rnorm(2*global$num.loci, mean=0, sd=sqrt(global$var.init/2/global$num.loci)), 
		ncol=2)
	list(
		genotype  = genotype, 
		genot.value= GPmap(genotype),
		phenotype = get.phenotype(genotype),
		fitness   = 1
	)
}
```
You can try to generate a single individual by typing 
ind<-init.individual()
Each individual is represented by a list of 4 elements containing the genotype, the genotypic value, the phenotype and the fitness.
Access the genotype of that individual by typing
ind$genotype
you can observe the values of the two alleles at the 5 loci.
Now its genotypic value, by typing
ind$genot.value
which is the sum of the genotypes (the 10 elements contained in ind$genot.value)
and finally access its phenotype by
ind$phenotype

The function init.population generates as many random individuals (init.individual) as in the population (pop.size) and returns those individuals as a list. 
```{r}
init.population <- function(pop.size, var.init, num.loci) {
	replicate(global$pop.size, init.individual(var.init, num.loci), simplify=FALSE)
}
```
Now generate your population by typing
pop<-init.population(100, 1, 5)
This population contains 100 individuals, each of which is a list of the 4 elements described above.
The first individual can be accessed by 
pop[[1]] and its values pop[[1]]$genotype etc...
The last is pop[[100]]

The make.gamete function makes a haploid gamete out of an individual. Random drawing is done by the sample function which draws by chance either allele 1 or 2 (c(1,2)) as many times as the number of rows in the genotype of an indiidual (number of loci) with replacement.
Recombination rate is 0.5 between loci (=free recombination), so that drawing allele 1 or 2 at a given locus does not depend on the preceding value drawn. 
```{r}
make.gamete <- function(indiv) {
	indiv$genotype[cbind(1:nrow(indiv$genotype), sample(c(1,2), nrow(indiv$genotype), replace=TRUE))]
}
```
To understand what the function does you can type:
pop[[1]]$genotype to access the genotype of the first individual of the population
and type 
gam<-make.gamete(pop[[1]])
gam is a random drawing at each locus of one of the two alleles of pop[[1]]$genotype.
Visualize gam

The make.offspring function makes an individual by binding the genotypes of two gametes coming from two individuals (mother and father).
```{r}
make.offspring <- function(mother, father) {
	genotype <- cbind(make.gamete(mother), make.gamete(father))
	list(
		genotype  = genotype, 
		genot.value= GPmap(genotype),
		phenotype = get.phenotype(genotype),
		fitness   = 1
	)
}
```
For example type:
offspr<-make.offspring(pop[[1]],pop[[2]])
Now verify which allele is coming from which parent at each locus.

The update.fitness returns a new population object with updated fitnesses. 
If trunc.sel is set to 1, there is no selection so that the fitness of eery individual is 1.

```{r}
update.fitness <- function(population, trunc.sel) {
	phenotypes <- sapply(population, "[[", "phenotype")
	keep.indiv <- if (global$trunc.sel > 0) phenotypes >= quantile(phenotypes, prob=1-global$trunc.sel)
				  else 						phenotypes <= quantile(phenotypes, prob= -global$trunc.sel)
	mapply(population, keep.indiv, FUN=function(indiv, keep) { indiv$fitness <- if (keep) 1 else 0; indiv }, SIMPLIFY=FALSE)
}
```
You can type
fitness<-update.fitness(pop)
and check the update by looking at
pop[[1]]$fitness

```{r}
reproduction <- function(population) {
	# Returns the next generation
	fitnesses <- sapply(population, "[[", "fitness")
	replicate(	n       = global$pop.size, 
				expr    = make.offspring(
							mother=unlist(sample(population, 1, prob=fitnesses)), 
							father=unlist(sample(population, 1, prob=fitnesses))),
				simplify = FALSE)
}
```

```{r}
reproduction <- function(population) {
	# Returns the next generation
	fitnesses <- sapply(population, "[[", "fitness")
	replicate(	n       = global$pop.size, 
				expr    = {
				      mother <- unlist(sample(population, 1, prob=fitnesses))
				      father <- if(runif(1) < self.rate) 
				                  mother 
				                else 
				                  unlist(sample(population, 1, prob=fitnesses))
				      make.offspring(mother, father)
				      },
				simplify = FALSE)
}
```

```{r}
summary.population <- function(population) {
	# Computes summary statistics for the population
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

```{r}
simulation <- function(generations=20, pop.size = 100, num.loci = 5, var.init = 1, var.env  = 1, trunc.sel= 1) {
	# Runs a simulation
	pop <- init.population(pop.size=pop.size)
	summ <- data.frame()
	for (gg in 1:generations) {
		pop <- update.fitness(pop)
		summ <- rbind(summ, summary.population(pop))
		if (gg < generations)
			pop <- reproduction(pop)
	}
	summ
}
```
