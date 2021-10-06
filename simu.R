
global <- list(
	pop.size = 100,			# Population size every generation
	num.loci = 5,			# Number of (haploid) loci
	var.init = 1,			# Initial genetic variance
	var.env  = 1,			# Environmental (residual) variance
	trunc.sel= 1			# Strength of truncation selection
							#  (Part of the population kept for breeding, 1 = 100% => No selection)
							#  (negative values stand for selection towards smaller phenotypes)
)


GPmap <- function(genotype) {
	# Returns the genotypic value (mean phenotype) corresponding to a genotype
	sum(genotype)
}

get.phenotype <- function(genotype) {
	# Returns a phenotype value corresponding to a specific genotype. Environmental effets are accounted for. 
	rnorm(1, mean=GPmap(genotype), sd=sqrt(global$var.env))
}

init.individual <- function() {
	# Generates a random individual for the starting population
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

init.population <- function() {
	# Generates the initial population
	replicate(global$pop.size, init.individual(), simplify=FALSE)
}

make.gamete <- function(indiv) {
	# Makes a haploid gamete out of an individual. Recombination rate is 0.5 (free recombination)
	if (length(indiv) == 1) indiv <- indiv[[1]]
	indiv$genotype[cbind(1:nrow(indiv$genotype), sample(c(1,2), nrow(indiv$genotype), replace=TRUE))]
}

make.offspring <- function(mother, father) {
	# Makes an individual out of two parents. 
	genotype <- cbind(make.gamete(mother), make.gamete(father))
	list(
		genotype  = genotype, 
		genot.value= GPmap(genotype),
		phenotype = get.phenotype(genotype),
		fitness   = 1
	)
}

update.fitness <- function(population) {
	# Returns a new population object with updated fitnesses. 
	phenotypes <- sapply(population, "[[", "phenotype")
	keep.indiv <- if (global$trunc.sel > 0) phenotypes >= quantile(phenotypes, prob=1-global$trunc.sel)
				  else 						phenotypes <= quantile(phenotypes, prob= -global$trunc.sel)
	mapply(population, keep.indiv, FUN=function(indiv, keep) { indiv$fitness <- if (keep) 1 else 0; indiv }, SIMPLIFY=FALSE)
}

reproduction <- function(population) {
	# Returns the next generation
	fitnesses <- sapply(population, "[[", "fitness")
	replicate(	n       = global$pop.size, 
				expr    = make.offspring(
							mother=sample(population, 1, prob=fitnesses), 
							father=sample(population, 1, prob=fitnesses)),
				simplify = FALSE)
}

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

simulation <- function(generations=20) {
	# Runs a simulation
	pop <- init.population()
	summ <- data.frame()
	for (gg in 1:generations) {
		pop <- update.fitness(pop)
		summ <- rbind(summ, summary.population(pop))
		if (gg < generations)
			pop <- reproduction(pop)
	}
	summ
}
