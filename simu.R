library(abind)

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
	fitness      = "gaussian",
	optimization = "none" # coulf be "cmpfun" or "c++"
)

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

get.phenotype <- function(
		genotype, 
		var.env = default$var.env) 
{
	# Returns a phenotype value corresponding to a specific genotype. Environmental effets are accounted for. 
	rnorm(1, mean=GPmap(genotype), sd=sqrt(var.env))
}

init.individual <- function(
		var.init = default$var.init, 
		num.loci = default$num.loci, 
		var.env  = default$var.env) 
{
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

init.population <- function(
		pop.size = default$pop.size, 
		var.init = default$var.init, 
		num.loci = default$num.loci, 
		var.env  = default$var.env) 
{
	# Generates the initial population	
	replicate(pop.size, init.individual(var.init, num.loci, var.env), simplify=FALSE)
}



make.gamete.R <- function(
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

if(require(compiler)) 
	make.gamete.cmpfun <- cmpfun(make.gamete.R)

if (require(Rcpp))
cppFunction('
NumericVector makeGameteCPP(
	const List&          indiv, 
	double               rateMut, 
	double               varMut, 
	const NumericVector& rateRec) {
	    unsigned int        nloc     = rateRec.size() + 1;
	    const NumericMatrix genotype = as<NumericMatrix>(indiv["genotype"]); 
		const NumericVector recrand  = Rcpp::runif(nloc, 0.0, 1.0);
		NumericVector       gam(nloc);
		
		unsigned int curpar = recrand[0] < 0.5 ? 1 : 0;
		for (unsigned int loc = 0; loc < nloc; loc++) {
			gam[loc] = genotype(loc, curpar);
			if (recrand[loc+1] < rateRec[loc]) 
				curpar = (curpar + 1 ) % 2;
		}
		if ( (rateMut > 0.0) && (R::runif(0.0, 1.0) < rateMut) ){
			unsigned int locmut = floor(R::runif(0.0, nloc));
			gam[locmut] += R::rnorm(0.0, sqrt(varMut));
		}
		
		return(gam);
	}
')

make.gamete <- function(
		indiv, 
		rate.mut = default$rate.mut, 
		var.mut  = default$var.mut,
		rate.rec = default$rate.rec, 
		optimization = default$optimization)
{
	FUN.gamete <- make.gamete.R
	if (optimization == "cmpfun") FUN.gamete <- make.gamete.cmpfun
	if (optimization == "c++")    FUN.gamete <- makeGameteCPP
	FUN.gamete(indiv, rate.mut, var.mut, rate.rec)
} 	 


make.offspring <- function(
		mother, 
		father, 
		var.env  = default$var.env, 
		rate.mut = default$rate.mut, 
		var.mut  = default$var.mut, 
		rate.rec = default$rate.rec,
		optimization = default$optimization) 
{
	# Makes an individual out of two parents. 
	genotype <- cbind(
		make.gamete(mother, rate.mut, var.mut, rate.rec, optimization),
		make.gamete(father, rate.mut, var.mut, rate.rec, optimization))

	list(
		genotype  = genotype, 
		genot.value= GPmap(genotype),
		phenotype = get.phenotype(genotype, var.env),
		fitness   = 1
	)
}

make.clonal.offspring <- function(
		parent, 
		var.env  = default$var.env, 
		rate.mut = default$rate.mut, 
		var.mut  = default$var.mut)
{
	genotype <- parent$genotype
	# Not very clean, but in order to keep exactly the same shortcomings as sexual reproduction 
	# (when rate.mut is large for instance), the best is to repeat the mutation procedure twice
	if (rate.mut > 0 && runif(1) < rate.mut) {
		mut.loc <- sample(seq_len(nrow(genotype)), 1)
		genotype[mut.loc,1] <- rnorm(1, mean=genotype[mut.loc,1], sd=sqrt(var.mut))
	}
	if (rate.mut > 0 && runif(1) < rate.mut) {
		mut.loc <- sample(seq_len(nrow(genotype)), 1)
		genotype[mut.loc,2] <- rnorm(1, mean=genotype[mut.loc,2], sd=sqrt(var.mut))
	}
	
	list(
		genotype  = genotype, 
		genot.value= GPmap(genotype),
		phenotype = get.phenotype(genotype, var.env),
		fitness   = 1
	)
}

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


reproduction <- function(
		population, 
		pop.size     = default$pop.size, 
		var.env      = default$var.env, 
		rate.mut     = default$rate.mut, 
		var.mut      = default$var.mut,
		rate.rec     = default$rate.rec,
		rate.selfing = default$rate.selfing,
		rate.clonal  = default$rate.clonal, 
		optimization = default$optimization) 
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
			rate.rec = rate.rec, 
			optimization = optimization),
		SIMPLIFY=FALSE)
	
	return(c(clones, selfers, outcros)) # The order is not expected to matter
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

migration <- function(pops, rate.migr = default$rate.migr) {
   npop <- length(pops)
   stopifnot(npop > 1)
   for (i in 1:(npop-1))
     for (j in (i+1):npop) {
       num.migr <- rbinom(1, min(length(pops[[i]]), length(pops[[j]])), rate.migr)
       if (num.migr == 0) next
       from.i <- sample(1:length(pops[[i]]), num.migr, replace=FALSE)
       from.j <- sample(1:length(pops[[j]]), num.migr, replace=FALSE)
       tmp.pop <- pops[[i]][from.i]
       pops[[i]][from.i] <- pops[[j]][from.j]
       pops[[j]][from.j] <- tmp.pop
     }
   pops
}

clean.inpop <- function(obj) {
	# Tries to find a proper population in object obj
	if (is.list(obj)) {
		if (is.data.frame(obj)) {
			pp <- attr(obj, "lastpop")
			if (is.null(pp)) {
				stop("Impossible to find a population in the provided object. Have you forgotten the output.pop=TRUE option?")
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
		pop2, 
		numcross     = length(pop1), 
		var.env      = default$var.env,
		rate.rec     = default$rate.rec,
		sel.Vs       = default$sel.Vs,
		sel.trunc    = default$sel.trunc,
		sel.optimum  = default$sel.optimum,
		fitness      = default$fitness,
		output.pop   = FALSE) 
{
	pop1 <- clean.inpop(pop1)
	pop2 <- clean.inpop(pop2)
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
	cross <- update.fitness(cross, sel.Vs, sel.optimum, sel.trunc, fitness)
	
	summ <- summary.population(cross)
	if (output.pop)
		attr(summ, "lastpop") <- cross
	summ
}

simulation1pop <- function(
		generations  = 20, 
		pop.size     = default$pop.size, 
		num.loci     = default$num.loci, 
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
		optimization = default$optimization,
		input.pop    = NULL, 
		output.pop   = FALSE) 
{
	if (!is.null(input.pop)) {
		pop <- clean.inpop(input.pop)
		stopifnot(nrow(pop[[1]]$genotype) == num.loci) # Number of loci is the only parameter that cannot change
	} else {
		pop <- init.population(pop.size=pop.size, var.init=var.init, num.loci=num.loci, var.env=var.env)
	}
	summ <- if (is.data.frame(input.pop)) {attr(input.pop, "lastpop") <- NULL; input.pop} else data.frame()
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
	if (output.pop)
		attr(summ, "lastpop") <- pop
	summ
}

simulationNpop <- function(
		generations  = 20, 
		pop.size     = default$pop.size, 
		num.loci     = default$num.loci, 
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
		fitness      = default$fitness,
		optimization = default$optimization,
		input.pop    = NULL,
		output.pop   = FALSE) 
{
	# Parameters that can be vectorized (one value for each population)
	pop.size    <- rep(pop.size,    length.out=num.pop)
	var.init    <- rep(var.init,    length.out=num.pop)
	var.env     <- rep(var.env,     length.out=num.pop)
	sel.Vs      <- rep(sel.Vs,      length.out=num.pop)
	sel.optimum <- rep(sel.optimum, length.out=num.pop)
	sel.trunc   <- rep(sel.trunc,   length.out=num.pop)
	rate.mut    <- rep(rate.mut,    length.out=num.pop)
	var.mut     <- rep(var.mut,     length.out=num.pop)
	rate.selfing<- rep(rate.selfing,length.out=num.pop)
	rate.clonal <- rep(rate.clonal, length.out=num.pop)
	fitness     <- rep(fitness,     length.out=num.pop)
	# rate.rec    <- rep(rate.rec,    length.out=num.pop) # Not possible yet, already a vector over loci
	 
	if (!is.null(input.pop)) {
		pops <- clean.inpop(input.pop)
	} else {
		pops <- lapply(seq_len(num.pop), function(i) init.population(pop.size=pop.size[i], var.init=var.init[i], num.loci=num.loci, var.env=var.env[i]))
	}
	
	summ <- if (is.data.frame(input.pop)) {attr(input.pop, "lastpop") <- NULL; input.pop} else data.frame()
	for (gg in 1:generations) {
		pops <- lapply(seq_len(num.pop), function(i) update.fitness(pops[[i]], sel.Vs[i], sel.optimum[i], sel.trunc[i], fitness[i]))
		summ <- rbind(summ, cbind(do.call(cbind, lapply(pops, summary.population)), summary.population(unlist(pops, recursive=FALSE))))
		if (gg < generations)
			pops <- lapply(seq_len(num.pop), function(i) reproduction(
						pops[[i]], 
						pop.size     = pop.size[i], 
						var.env      = var.env[i], 
						rate.mut     = rate.mut[i], 
						var.mut      = var.mut[i], 
						rate.rec     = rate.rec, 
						rate.selfing = rate.selfing[i], 
						rate.clonal  = rate.clonal[i], 
						optimization = optimization)
					)
			pops <- migration(pops, rate.migr)
	}
	colnames(summ) <- paste0(colnames(summ), ".", rep(c(as.character(1:num.pop), "All"), each=ncol(summ)/(num.pop+1)))
	if (output.pop)
		attr(summ, "lastpop") <- pops
	summ
}



simulation <- function(
		generations  = 20, 
		pop.size     = default$pop.size, 
		num.loci     = default$num.loci, 
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
		fitness      = default$fitness, # can be "gaussian" or "truncation"
		optimization = default$optimization,   # can be "none", "cmpfun", or "c++"
		input.pop    = NULL,
		output.pop   = FALSE) 
{
	# Checks and adjust parameters
	stopifnot(
		generations >= 1,
		pop.size    >= 1,
		num.loci    >= 1,
		var.init    >= 0.0,
		var.env     >= 0.0,
		rate.mut    >= 0.0, rate.mut <= 1.0,
		var.mut     >= 0.0,
		all(rate.rec >= 0.0), all(rate.rec <= 0.5),
		rate.selfing >= 0.0, rate.selfing <= 1.0,
		rate.clonal  >= 0.0, rate.clonal  <= 1.0,
		rate.selfing + rate.clonal <= 1.0,
		fitness %in% c("gaussian", "truncation"),
		fitness == "gaussian" | (sel.trunc >= -1.0 & sel.trunc <= 1.0),
		num.pop      >= 1, 
		rate.migr    >= 0.0, rate.migr <= 1.0) 
		
	rate.rec <- rep_len(rate.rec, num.loci - 1)
	
	if (num.pop == 1) {
		simulation1pop(
			generations, pop.size, num.loci, var.init, var.env, sel.Vs, sel.optimum, sel.trunc, rate.mut, var.mut, rate.rec, rate.selfing, rate.clonal, fitness, optimization, input.pop, output.pop)
	} else {
		simulationNpop(
			generations, pop.size, num.loci, var.init, var.env, sel.Vs, sel.optimum, sel.trunc, rate.mut, var.mut, rate.rec, rate.selfing, rate.clonal, num.pop, rate.migr, fitness, optimization, input.pop, output.pop)
	}
}

