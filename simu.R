library(abind)

default <- list(
	pop.size     = 100, 
	num.loci     = 5, 
        mean.init    = c(0, 1),     # intercept, slope
	var.init     = 1.0*diag(2), # intercept, slope
	var.env      = 1.0, 
	sel.Vs       = 1.0, 
	sel.optimum  = c(0,1),      # intercept, slope 
	rate.mut     = 0.0, 
	var.mut      = 1.0*diag(2), # intercept, slope
	rate.rec     = 0.5,
        env.genFUN   = function(g) 0,
	rate.selfing = 0.0,
	rate.clonal  = 0.0, 
	fitness      = "gaussian"
)

# Function to get mean of repetitions
  list.sim.mean <- function(ll) {
    arr <- do.call(abind, c(ll, list(along=3)))
    ans <- as.data.frame(rowMeans(arr, dims=2))
    rownames(ans) <- rownames(ll[[1]])
    ans
  }

 rowVars <- function(x) {
    v <- var(t(x))
    setNames(v[upper.tri(v,diag=TRUE)], nm=outer(seq_len(nrow(x)),seq_len(nrow(x)),\(i,j) ifelse(i==j, as.character(i), paste(i,j,sep=".")))[upper.tri(v,diag=TRUE)])
 }

GPEmap <- function(genotype, environment) {
	# Returns the genotypic value (mean phenotype) corresponding to a genotype
	intercept <- sum(genotype[,1:2])
        slope     <- sum(genotype[,3:4])
        intercept + slope*environment
}

get.phenotype <- function(
		genotype, 
                environment,
		var.env = default$var.env) 
{
	# Returns a phenotype value corresponding to a specific genotype. Environmental variance + plasticity are accounted for. 
	rnorm(1, GPEmap(genotype,environment), sqrt(var.env))
}

init.individual <- function(
                mean.init= default$mean.init,
		var.init = default$var.init, 
		num.loci = default$num.loci, 
		var.env  = default$var.env,
                environment= default$env.genFUN(0)) 
{
	# Generates a random individual for the starting population
	genotype <- matrix(
			MASS::mvrnorm(2*num.loci, mu=mean.init/2/num.loci, Sigma=var.init/2/num.loci),
			ncol=2*2)
	list(
		genotype  = genotype, 
		genot.value= GPEmap(genotype, environment),
		phenotype = get.phenotype(genotype, environment, var.env),
		fitness   = 1
	)
}

init.population <- function(
		pop.size = default$pop.size, 
                mean.init= default$mean.init,
		var.init = default$var.init, 
		num.loci = default$num.loci, 
		var.env  = default$var.env,
                environment= default$env.genFUN(0)) 
{
	# Generates the initial population	
	replicate(pop.size, init.individual(mean.init, var.init, num.loci, var.env, environment), simplify=FALSE)
}



make.gamete.R <- function(
		indiv, 
		rate.mut = default$rate.mut, 
		var.mut  = default$var.mut,
		rate.rec = default$rate.rec) 
{
	num.traits <- ncol(indiv$genotype)/2
        num.loc    <- nrow(indiv$genotype)
	# Recombination
	recs <- cumsum(runif(length(rate.rec)+1) < c(0.5, rate.rec))
	gam <- indiv$genotype[cbind(rep(1:num.loc, num.traits), 2*rep(1:num.traits, each=num.loc)-1+(rep(recs %% 2, num.traits)))] |> 
                matrix(nrow=num.loc)
	# Mutation
	if (rate.mut > 0 && runif(1) < rate.mut) {
		mut.loc <- sample(seq_len(nrow(gam)), 1)
		gam[mut.loc,] <- MASS::mvrnorm(1, mu=gam[mut.loc,], Sigma=var.mut)
	}
	gam
}


make.gamete <- function(
		indiv, 
		rate.mut = default$rate.mut, 
		var.mut  = default$var.mut,
		rate.rec = default$rate.rec)
{
	FUN.gamete <- make.gamete.R
	make.gamete.R(indiv, rate.mut, var.mut, rate.rec)
} 	 


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
		make.gamete(father, rate.mut, var.mut, rate.rec))[,c(matrix(1:ncol(mother$genotype),ncol=2,byrow=TRUE))]

	list(
		genotype  = genotype, 
		genot.value= NA,
		phenotype = NA,
		fitness   = NA
	)
}

make.clonal.offspring <- function(
		parent, 
		var.env  = default$var.env, 
		rate.mut = default$rate.mut, 
		var.mut  = default$var.mut,
                environment)
{
	genotype <- parent$genotype
	num.traits <- ncol(genotype)/2
	# Not very clean, but in order to keep exactly the same shortcomings as sexual reproduction 
	# (when rate.mut is large for instance), the best is to repeat the mutation procedure twice
	sq1 <- seq(1,2*num.traits,by=2)
	sq2 <- seq(2,2*num.traits,by=2)
	if (rate.mut > 0 && runif(1) < rate.mut) {
		mut.loc <- sample(seq_len(nrow(genotype)), 1)
		genotype[mut.loc,sq1] <- MASS::mvrnorm(1, mu=genotype[mut.loc,sq1], Sigma=var.mut)
	}
	if (rate.mut > 0 && runif(1) < rate.mut) {
		mut.loc <- sample(seq_len(nrow(genotype)), 1)
		genotype[mut.loc,sq2] <- MASS::mvrnorm(1, mu=genotype[mut.loc,sq2], Sigma=var.mut)
	}
	
	list(
		genotype  = genotype, 
		genot.value= NA,
		phenotype = NA,
		fitness   = NA
	)
}

update.fitness <- function(
		population, 
		sel.Vs       = default$sel.Vs, 
		sel.optimum  = default$sel.optimum,
		sel.trunc    = default$sel.trunc,
		fitness      = default$fitness,
                environment) 
{
	# Returns a new population object with updated fitnesses. 

	if (fitness == "gaussian") {
		lapply(population, function(indiv) { 
                                indiv$genot.value <- GPEmap(indiv$genotype, environment)
                                indiv$phenotype   <- get.phenotype(indiv$genotype, environment)
				dd <- indiv$phenotype-(sel.optimum[1] - environment*sel.optimum[2])
				indiv$fitness     <- exp(-0.5*dd^2 / sel.Vs)
				indiv })
	} 
}


reproduction <- function(
		population, 
		pop.size     = default$pop.size, 
		var.env      = default$var.env, 
		rate.mut     = default$rate.mut, 
		var.mut      = default$var.mut,
		rate.rec     = default$rate.rec,
		rate.selfing = default$rate.selfing,
		rate.clonal  = default$rate.clonal) 
{
	# Returns the next generation
	fitnesses <- sapply(population, "[[", "fitness")
	num.clones  <- rbinom(1, pop.size, prob=rate.clonal)
	num.selfers <- rbinom(1, pop.size - num.clones, prob=if(rate.clonal == 1) 0 else rate.selfing/(1-rate.clonal))
	num.outcros <- pop.size - num.clones - num.selfers
	
	clones  <- sample(population, num.clones, prob=fitnesses, replace=TRUE)
	clones  <- lapply(clones, make.clonal.offspring, var.env=var.env, rate.mut=rate.mut, var.mut=var.mut)
	
	parent.selfers <- sample(population, num.selfers, prob=fitnesses, replace=TRUE)
	selfers <- lapply(parent.selfers, function(p) 
		make.offspring(
			mother   = p, 
			father   = p,
			var.env  = var.env ,
			rate.mut = rate.mut, 
			var.mut  = var.mut, 
			rate.rec = rate.rec)
		)
	
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
	
	return(c(clones, selfers, outcros)) # The order is not expected to matter
}

summary.population <- function(population) {
	# Computes summary statistics for the population
	phenotypes <- sapply(population, "[[", "phenotype")
	genot.val  <- sapply(population, "[[", "genot.value")
        genot.int  <- sapply(population, function(indiv) sum(indiv$genotype[,1:2]))
        genot.slp  <- sapply(population, function(indiv) sum(indiv$genotype[,3:4]))
	fitnesses  <- sapply(population, "[[", "fitness")
	htz        <- sapply(population, function(ind) mean(ind$genotype[,1] != ind$genotype[,2])) # not very precise in a multivariate context
	data.frame(
		phen.mean    = mean(phenotypes), 
		phen.var     = var (phenotypes),
		gen.mean     = mean(genot.val),
		gen.var      = var (genot.val),
                gen.int.mean = mean(genot.int),
                gen.int.var  = var (genot.int),
                gen.slp.mean = mean(genot.slp),
                gen.slp.var  = var (genot.slp),
                gen.int.slp.cov = cov(genot.int, genot.slp),
		fit.mean     = mean(fitnesses),
		fit.var      = var (fitnesses),
		htz.rate     = mean(htz),
		sel.diff     = mean(fitnesses*phenotypes)/mean(fitnesses) - mean(phenotypes)
	)
}


clean.inpop <- function(obj, n=1) {
	# Tries to find a proper population in object obj
	if (is.list(obj)) {
		if (is.data.frame(obj)) {
			pp <- attr(obj, "lastpop")
			if (is.null(pp)) {
				stop("Impossible to find a population in the provided object. Have you forgotten the output.pop=TRUE option?")
			}
                        if (is.list(pp[[1]][[1]])) { # pp is not a population, perhaps a list of populations?
                                return(pp[[n]])
                        }
			return(pp)
		}
		return(obj) # Hoping that the list can be interpreted as a population, otherwise the crash will happen later. 
	} else {
		stop("The object cannot be interpreted as a population")
	}
}

crosspopulations <- function(
		pop1, 
		pop2         = NULL, 
		numcross     = length(pop1), 
		var.env      = default$var.env,
		rate.rec     = default$rate.rec,
		inv.sel.Vs   = solve(default$sel.Vs),
		sel.trunc    = default$sel.trunc,
		sel.optimum  = default$sel.optimum,
		fitness      = default$fitness,
                environment  = default$env.genFUN(0),
		output.pop   = FALSE) 
{
        if (!is.null(pop2)) {
	        pop1 <- clean.inpop(pop1)
	        pop2 <- clean.inpop(pop2)
        } else {
                pop1 <- clean.inpop(pop1, 1)
                pop2 <- clean.inpop(pop1, 2)
        }
	stopifnot(length(pop1) > 0, length(pop2) > 0, 
	          nrow(pop1[[1]]$genotype) == nrow(pop2[[1]]$genotype))
	rate.rec <- rep_len(rate.rec, nrow(pop1[[1]]$genotype) - 1)
	cross <- replicate(numcross, 
		expr= {
				parent1 <- unlist(sample(pop1, 1), recursive=FALSE)
				parent2 <- unlist(sample(pop2, 1), recursive=FALSE)
				make.offspring(
					mother   = parent1, 
					father   = parent2,
					var.env  = var.env ,
					rate.mut = 0.0, 
					var.mut  = 0.0, 
					rate.rec = rate.rec)
				},
		simplify=FALSE)
	cross <- update.fitness(cross, sel.Vs, sel.optimum, sel.trunc, fitness, environment)
	
	summ <- summary.population(cross)
	if (output.pop)
		attr(summ, "lastpop") <- cross
	summ
}

simulation1pop <- function(
		generations  = 20, 
		pop.size     = default$pop.size, 
		num.loci     = default$num.loci, 
                mean.init    = default$mean.init,
		var.init     = default$var.init, 
		var.env      = default$var.env, 
		sel.Vs       = default$sel.Vs, 
		sel.optimum  = default$sel.optimum, 
		sel.trunc    = default$sel.trunc,
		rate.mut     = default$rate.mut, 
		var.mut      = default$var.mut, 
		rate.rec     = default$rate.rec,
		rate.selfing = default$rate.selfing,
		rate.clonal  = default$rate.clonal,
		fitness      = default$fitness,
                env.genFUN   = default$env.genFUN,
		input.pop    = NULL, 
		output.pop   = FALSE) 
{
	if (!is.null(input.pop)) {
		pop <- clean.inpop(input.pop)
		stopifnot(nrow(pop[[1]]$genotype) == num.loci, ncol(pop[[1]]$genotype == num.traits)) 
		# Number of loci and number of traits are the only parameters that cannot change
	} else {
		pop <- init.population(pop.size=pop.size, var.init=var.init, num.loci=num.loci, var.env=var.env, environment=env.genFUN(0))
	}
	summ <- if (is.data.frame(input.pop)) {attr(input.pop, "lastpop") <- NULL; input.pop} else data.frame()
	for (gg in 1:generations) {
                environment <- env.genFUN(gg)
		pop <- update.fitness(pop, sel.Vs, sel.optimum, sel.trunc, fitness, environment=environment)
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
						rate.clonal  = rate.clonal)
	}
	if (output.pop)
		attr(summ, "lastpop") <- pop
	summ
}



simulation <- function(
		generations  = 20, 
		pop.size     = default$pop.size, 
		num.loci     = default$num.loci, 
                mean.init    = default$mean.init,
		var.init     = default$var.init, 
		var.env      = default$var.env, 
		sel.Vs       = default$sel.Vs, 
		sel.optimum  = default$sel.optimum, 
		sel.trunc    = default$sel.trunc, 
		rate.mut     = default$rate.mut, 
		var.mut      = default$var.mut, 
		rate.rec     = default$rate.rec,
		rate.selfing = default$rate.selfing,
		rate.clonal  = default$rate.clonal,
		num.pop      = default$num.pop,
		rate.migr    = default$rate.migr,
		fitness      = default$fitness, # can be "gaussian" only
                env.genFUN   = default$env.genFUN,
		input.pop    = NULL,
		output.pop   = FALSE) 
{
	# Checks and adjust parameters
	stopifnot(
		generations >= 1,
		pop.size    >= 1,
		num.loci    >= 1,
                length(mean.init) == 2,
		is.matrix(var.init) && ncol(var.init) == 2 && nrow(var.init) == 2,
		var.env     >= 0,
		is.matrix(var.mut)  && ncol(var.mut)  == 2 && nrow(var.mut)  == 2,
		sel.Vs      > 0,
		rate.mut    >= 0.0, rate.mut <= 1.0,
		all(rate.rec >= 0.0), all(rate.rec <= 0.5),
		rate.selfing >= 0.0, rate.selfing <= 1.0,
		rate.clonal  >= 0.0, rate.clonal  <= 1.0,
		rate.selfing + rate.clonal <= 1.0,
		fitness %in% c("gaussian")
		) 
		
	rate.rec <- rep_len(rate.rec, num.loci - 1)
	
	simulation1pop(
			generations, pop.size, num.loci, 
			mean.init, var.init, var.env, sel.Vs, sel.optimum, sel.trunc, rate.mut, var.mut, rate.rec, rate.selfing, rate.clonal, 
			fitness, env.genFUN=env.genFUN, input.pop, output.pop)

}

