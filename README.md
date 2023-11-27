# Presentation of the Simulation program

## Parameters

The program defines a series of default parameters, that will be used in the simulations unless specified otherwise. 

```{r}
}
default <- list(
	pop.size     = 100, 
	num.loci     = 5, 
	var.init     = 1.0, 
	var.env      = 1.0, 
	sel.Vs       = 1.0, 
	sel.optimum  = 0.0, 
	sel.trunc    = 0.0,
	rate.mut     = 0.0000, 
	var.mut      = 1.0,
	rate.rec     = 0.5,
	rate.selfing = 0.0,
	rate.clonal  = 0.0, 
	num.pop      = 1,
	rate.migr    = 0.0,
	fitness      = "gaussian"
)
```

## Simulation objects

Simulations are individual-based, each individual is characterized by:
* Its genotype, a 2-column, num.loci-row matrix (one column per haplotype). Allelic values are real numbers.
* Its genotypic value (genot.value), the average phenotype corresponding to the genotype
* Its phenotype, accounting for the (random) effect of the environment
* Its fitness, proportional to the probability to reproduce. 

```{r}
$genotype
             [,1]       [,2]
[1,]  0.3327662   0.3327662
[2,] -0.3111255   0.1121897
[3,] -0.0073502   0.1664419
[4,] -0.1919904  -0.1919904
[5,] -0.0073234 - 0.4895492

$genot.value
[1] -0.2551651

$phenotype
[1] -0.7695404

$fitness
[1] 1
```

Populations are list of individuals. 

## Genotype - Phenotype

Individuals are diploid, the genotype-phenotype map is additive. 

The GPmap function returns the genotypic value (mean phenotype) corresponding to a genotype.  

```{r}
GPmap <- function(genotype) {
	# Returns the genotypic value (mean phenotype) corresponding to a genotype
	sum(genotype)
}
```
The get.phenotype function returns a phenotype value (1) from a Normal distribution (rnorm) of (mean) the genotypic value, and standard deviation (sd) the square root of the environmental variance as defined for the population.  

```{r}
get.phenotype <- function(
		genotype, 
		var.env = default$var.env) 
{
	# Returns a phenotype value corresponding to a specific genotype. Environmental effets are accounted for. 
	rnorm(1, mean=GPmap(genotype), sd=sqrt(var.env))
}
```
## Fitness

There are two families of fitness functions,
* fitness ==  "gaussian" defines fitness as a bell-shaped stabilizing function, defined by its optimum (the phenotype for which fitness == 1) and the selection strength Vs (the variance of the fitness curve, i.e. large values correspond to small selection). 

$$ w(z) = \exp ( -\frac{1}{2} \frac{(z - \text{sel.optimum})^2}{\text{sel.Vs}} ) $$

* fitness == "truncation" implements truncation selection, a procedure often associated with artificial selection. All individuals above or below a threshold are kept, all the others are discarded. The proportion of individuals discarded is given by sel.trunc (sel.trunc = 0: no selection, sel.trunc = 1: no one survives), the sign of sel.trunc defines the direction of selection (if negative, individuals with the lowest phenotype are kept). 

```{r}
update.fitness <- function(
		population, 
		sel.Vs       = default$sel.Vs, 
		sel.optimum  = default$sel.optimum,
		sel.trunc    = default$sel.trunc,
		fitness      = default$fitness) 
{
	# if fitness=="truncation", abs(sel.trunc) stands for the part of the population discarded (and the sign stands for the direction)
	
	# Returns a new population object with updated fitnesses. 

	if (fitness == "gaussian") {
		lapply(population, function(indiv) { indiv$fitness <- exp(-(indiv$phenotype-sel.optimum)^2 / 2 / sel.Vs); indiv })
	} else if (fitness == "truncation") {
		pp  <- sapply(population, "[[", "phenotype")
		if (sel.trunc > 0) {
			thr <- quantile(pp, probs=sel.trunc)
			lapply(population, function(indiv) { indiv$fitness <- if (indiv$phenotype >= thr) 1.0 else 0.0; indiv })
		} else {
			thr <- quantile(pp, probs=1-abs(sel.trunc))
			lapply(population, function(indiv) { indiv$fitness <- if (indiv$phenotype <= thr) 1.0 else 0.0; indiv })
		}
	}
}
``` 

## Simulation loop

The simulation loop consists in calling the reproduction() routine recursively up to the number of desired generations.

```{r}
	for (gg in 1:generations) {
		pop <- update.fitness(pop, sel.Vs, sel.optimum, sel.trunc, fitness)
		summ <- rbind(summ, summary.population(pop))
		if (gg < generations)
			pop <- reproduction(
						pop, 
						pop.size     = pop.size, 
						var.env      = var.env, 
						rate.mut     = rate.mut, 
						var.mut      = var.mut, 
						rate.rec     = rate.rec, 
						rate.selfing = rate.selfing, 
						rate.clonal  = rate.clonal, 
						optimization = optimization)
	}
```

Reproduction can be clonal or sexual. Sexual reproduction itself can rely on outcrossing (random mating) or selfing. This is set by two rates: rate.selfing and rate.clonal. When both are set to 0 (default), the population reproduces by random mating. The corresponding part of the code is below:

```{r}
num.outcros <- pop.size - num.clones - num.selfers

parent.outcros1 <- sample(population, num.outcros, prob=fitnesses, replace=TRUE)
parent.outcros2 <- sample(population, num.outcros, prob=fitnesses, replace=TRUE)
outcros <- mapply(parent.outcros1, parent.outcros2, FUN=function(p1, p2) 
	make.offspring(
		mother   = p1, 
		father   = p2,
		var.env  = var.env ,
		rate.mut = rate.mut, 
		var.mut  = var.mut, 
		rate.rec = rate.rec),
	SIMPLIFY=FALSE)
```

The routine make.offspring() gathers two gametes in a new individual object:

```{r}
make.offspring <- function(
		mother, 
		father, 
		var.env  = default$var.env, 
		rate.mut = default$rate.mut, 
		var.mut  = default$var.mut, 
		rate.rec = default$rate.rec) 
{
	# Makes an individual out of two parents. 
	genotype <- cbind(
		make.gamete(mother, rate.mut, var.mut, rate.rec),
		make.gamete(father, rate.mut, var.mut, rate.rec))

	list(
		genotype  = genotype, 
		genot.value= GPmap(genotype),
		phenotype = get.phenotype(genotype, var.env),
		fitness   = 1
	)
}
``` 

Finally, make.gamete() handles two important biological processes: recombination and mutation. 

```{r}
make.gamete <- function(
		indiv, 
		rate.mut = default$rate.mut, 
		var.mut  = default$var.mut,
		rate.rec = default$rate.rec) 
{
	# Recombination
	recs <- cumsum(runif(length(rate.rec)+1) < c(0.5, rate.rec))
	gam <- indiv$genotype[cbind(1:nrow(indiv$genotype), 1+(recs %% 2))]
#~ 	gam <- ifelse(recs %% 2 == 0, indiv$genotype[,1], indiv$genotype[,2])
	
	# Mutation
	if (rate.mut > 0 && runif(1) < rate.mut) {
		mut.loc <- sample(seq_along(gam), 1)
		gam[mut.loc] <- rnorm(1, mean=gam[mut.loc], sd=sqrt(var.mut))
	}
	gam
}
``` 

## Initialization

Two procedures to initialize the population.

* From the parameters. The routine init.population() is called, based on the following variables:
	* pop.size
	* var.init , the initial genetic variance
	* num.loci
	* var.env

* From a file (option input.file). This file had to be created from a previous simulation. 

## Simulation output

After having run for the requested number of generations, the simulation stops. The variable returned by the simulation routine is a data.frame with, for each generation, the average and the variance of a series of indicators:
	* phen.mean : the phenotypic mean
	* phen.var  : the phenotypic variance
	* gen.mean  : the genetic mean
	* gen.var   : the genetic variance (should be phen.var - var.env)
	* fit.mean  : the mean absolute fitness
	* fit.var   : the variance in absolute fitness
	* sel.diff  : the realized selection differential (measured from the individual fitnesses)

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
With the option output.pop=TRUE, the population from the last generation will be embedded as an attribute in the simulation output. It can then be used to start a new simulation, with the input.pop argument.

## Code examples

```{r}
# Run a simulation
s1 <- simulation(generations=20, pop.size=1000)

# s1 is a data.frame tracking the evolution of the population
head(s1)
    phen.mean phen.var   gen.mean   gen.var  fit.mean   fit.var      sel.diff
1  0.02388641 1.992135 0.02668777 1.0193897 0.5733783 0.1157158 -0.0078695807
2 -0.01749457 1.874539 0.03146603 0.8602849 0.5907936 0.1099869  0.0073106962
3  0.04127094 1.896728 0.03817701 0.8106886 0.6009933 0.1130402 -0.0480819267
4 -0.04555095 1.711153 0.02541822 0.7517501 0.6075296 0.1098529  0.0018137856
5  0.01933254 1.744701 0.04109299 0.8194891 0.6007816 0.1091351 -0.0007977639
6  0.02401994 1.848087 0.03743467 0.7537498 0.5960803 0.1087571 -0.0346293013

# Let's see how the population responds to selection
s2 <- simulation(generations=30, pop.size=1000, sel.optimum=2.0)
plot(s2$phen.mean, type="l", xlab="Generations", ylab="Phenotype")

# Reuse the population to start a new simulation
s3 <- simulation(generations=20, pop.size=1000, sel.optimum=1.0, output.pop=TRUE)
s4 <- simulation(generations=20, pop.size=100,  sel.optimum=0.0, input.pop=s3)
plot(s4$phen.mean, type="l", xlab="Generations", ylab="Phenotype")

# Replicate the evolutionary dynamics
s5 <- simulation(generations=20, pop.size=100, sel.optimum=c(3.0, 2.0, 1.0), num.pop=3, rate.migr=0.0)
plot(s5$phen.mean.1, type="l", xlab="Generations", ylab="Phenotype")
lines(s5$phen.mean.2, col="blue")
lines(s5$phen.mean.3, col="red")

# Ananlyse the final population
s6 <- simulation(generations=30, pop.size=1000, sel.optimum=2.0, output.pop=TRUE)
hist(sapply(attr(s6, "lastpop"), "[[", "fitness"), xlab="Fitness", main="Fitness distribution")
```
