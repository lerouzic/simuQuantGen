
global <- list(
	pop.size = 100,
	num.loci = 5,
	var.init = 1,
	var.env  = 1
)


GPmap <- function(genotype) {
	rnorm(1, mean=sum(genotype), sd=sqrt(global$var.env))
}

init.individual <- function() {
	genotype <- matrix(
		rnorm(2*global$num.loci, mean=0, sd=sqrt(global$var.init/2/global$num_loci)), 
		ncol=2)
	list(
		genotype  = genotype, 
		phenotype = GPmap(genotype),
		fitness   = 1
	)
}

init.population <- function() {
	replicate(global$pop.size, init.individual(), simplify="list")
}

make.gamete <- function(indiv) {
	indiv$genotype[1:nrow(indiv$genotype), sample(c(1,2), nrow(indiv$genotype), replace=TRUE)]
}

make.offspring <- function(mother, father) {
	genotype <- cbind(make.gamete(mother), make.gamete(father))
	list(
		genotype  = genotype, 
		phenotype = GPmap(genotype),
		fitness   = 1
	)
}

reproduction <- function(population) {
	fitnesses <- sapply(population, "[[", "fitness")
	replicate(	n       = global$pop.size, 
				expr    = make.offspring(
							mother=sample(population, 1, prob=fitnesses), 
							father=sample(population, 1, prob=fitnesses),
				simplify = "list")
}

