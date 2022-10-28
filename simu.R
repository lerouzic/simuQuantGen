library(abind)

# Function to get mean of repetitions
  list.sim.mean <- function(ll) {
    arr <- do.call(abind, c(ll, list(along=3)))
    ans <- as.data.frame(rowMeans(arr, dims=2))
    rownames(ans) <- rownames(ll[[1]])
    ans
  }


GPmap <- function(genotype) {
	# Returns the genotypic value (mean phenotype) corresponding to a genotype
	sum(genotype)
}

get.phenotype <- function(genotype, var.env) {
	# Returns a phenotype value corresponding to a specific genotype. Environmental effets are accounted for. 
	rnorm(1, mean=GPmap(genotype), sd=sqrt(var.env))
}

init.individual <- function(var.init, num.loci, var.env) {
	# Generates a random individual for the starting population
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
	# Generates the initial population	
	replicate(pop.size, init.individual(var.init, num.loci, var.env), simplify=FALSE)
}

make.gamete <- function(indiv, rate.mut, var.mut) {
	gam <- indiv$genotype[cbind(1:nrow(indiv$genotype), sample(c(1,2), nrow(indiv$genotype), replace=TRUE))]
	if (runif(1) < rate.mut) {
		mut.loc <- sample(seq_along(gam), 1)
		gam[mut.loc] <- rnorm(1, mean=gam[mut.loc], sd=sqrt(var.mut))
	}
	gam
}

make.offspring <- function(mother, father, var.env, rate.mut, var.mut) {
	# Makes an individual out of two parents. 
	genotype <- cbind(make.gamete(mother, rate.mut, var.mut), make.gamete(father, rate.mut, var.mut))
	list(
		genotype  = genotype, 
		genot.value= GPmap(genotype),
		phenotype = get.phenotype(genotype, var.env),
		fitness   = 1
	)
}

update.fitness <- function(population, sel.strength, sel.optimum) {
	# Returns a new population object with updated fitnesses. 
	# Fitness = exp(- (phenotype - sel.optimum)^2 / (2*sel.strength))
	lapply(population, function(indiv) { indiv$fitness <- exp(-(indiv$phenotype-sel.optimum)^2 / 2 / sel.strength); indiv })
}


reproduction <- function(population, pop.size, var.env, rate.mut, var.mut) {
	# Returns the next generation
	fitnesses <- sapply(population, "[[", "fitness")
	replicate(n=pop.size, 
		expr=make.offspring(
				mother=unlist(sample(population, 1, prob=fitnesses), recursive=FALSE), 
				father=unlist(sample(population, 1, prob=fitnesses), recursive=FALSE),
				var.env=var.env ,
				rate.mut=rate.mut, var.mut=var.mut),
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

simulation <- function(generations=20, pop.size = 100, num.loci = 5, var.init = 1, var.env = 1, sel.strength = 1, sel.optimum = 0, rate.mut=0.0001, var.mut=1, input.file=NULL, output.file=NULL) {
	# Runs a simulation
	pop <- if (!is.null(input.file)) {
			readRDS(input.file)
		} else {
			init.population(pop.size=pop.size, var.init=var.init, num.loci=num.loci, var.env=var.env)
		}
	summ <- data.frame()
	for (gg in 1:generations) {
		pop <- update.fitness(pop, sel.strength, sel.optimum)
		summ <- rbind(summ, summary.population(pop))
		if (gg < generations)
			pop <- reproduction(pop, pop.size=pop.size, var.env=var.env, rate.mut=rate.mut, var.mut=var.mut)
	}
	if (!is.null(output.file))
		saveRDS(pop, output.file)
	summ
}

# To run several repetitions, this command can be used
rep = 10 #number of repetitions
sims <- replicate(rep, simulation(), simplify=FALSE)
# Uses function list.sim.means() defined at the begining. The variable meansims is a dataframe with the mean over the repetitions of each variable over time
meansims<-list.sim.mean(sims)
