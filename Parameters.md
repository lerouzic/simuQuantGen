# Exploration of the simulation parameters
  
First we upload the code found in simu.r 

```{r}
source("simu.R")
set.seed(2)
```

Let's explore Population Size. Here we choose 3 values, 10, 100, 1000. We consider 50 generations, 100 loci, and a selection optimum of 10. First, we consider no selection by setting the sel.Vs to Inf. The parameter sel.Vs actually represets the inverse of the strength of selection (i.e. the greater it is, the weaker the strength of selection). We can follow th response to selection by following the evolution of the phenotypic mean across generations. With sel.Vs = Inf, there should be no selection response.  


```{r}
# If you want to explore more: PopSize <- seq(5,1000, by=50)
PopSize<-c(10,100,1000)
for (i in 1:length(PopSize)) {
  current_PopSize <- PopSize[i]
  current_name <- paste0("PopSize_",current_PopSize)
  current_simu <- simulation(generations=50,pop.size=current_PopSize,num.loci=100,var.init=1,var.env=1,sel.Vs=Inf,sel.optimum=10)
  assign(current_name,current_simu)
  rm(current_simu)
}

plot(NA,xlim=c(0,50),ylim=c(-10,10),xlab="generations", ylab=("phen.mean"))
for (i in 1:length(PopSize)) 
  {
  ps<-get(paste("PopSize_",PopSize[i],sep=""))
  lines(ps$phen.mean,col=i)
}
legend("bottomright", lty=1, col=1:length(PopSize), legend=paste0("N=", PopSize))
```

In absence of selection, the phenotype remains around its initial value. We observe however the effect of drift at N=10. Let's look at the genetic variance.


```{r}
plot(NA,xlim=c(0,50),ylim=c(0,2),xlab="generations", ylab=("gen.var"))
for (i in 1:length(PopSize)) 
  {
  ps<-get(paste("PopSize_",PopSize[i],sep=""))
  lines(ps$gen.var,col=i)
}
legend("bottomright", lty=1, col=1:length(PopSize), legend=paste0("N=", PopSize))
```

The genetic variance remains around the initial value for N=100 and N=1000, but presents greater fluctuations for N=100 than for N=1000. As for N=10, drift exhausts the initial genetic variance. Indeed, we know from theory that the decrease of genetic variance over time can be written Va(t) = Va(0) (1-1/2N)^t (neutral case). Hence the greater N, the smaller the change in Va (the additive genetic variance) between time t and t+1.

Let's now explore the effect of Population Size with selection (sel.Vs set to 5). We expect two phases: first directional selection towards the optimum and subsequently stabilizing selection. We look at the selection response (phenotypic mean across generations) to depict these two phases.



```{r}
# If you want to explore more: PopSize <- seq(5,1000, by=50)
PopSize<-c(10,100,1000)
for (i in 1:length(PopSize)) {
  current_PopSize <- PopSize[i]
  current_name <- paste0("PopSize_",current_PopSize)
  current_simu <- simulation(generations=50,pop.size=current_PopSize,num.loci=100,var.init=1,var.env=1,sel.Vs=5,sel.optimum=10)
  assign(current_name,current_simu)
  rm(current_simu)
}

plot(NA,xlim=c(0,50),ylim=c(0,12),xlab="generations", ylab=("phen.mean"))
for (i in 1:length(PopSize)) 
  {
  ps<-get(paste("PopSize_",PopSize[i],sep=""))
  lines(ps$phen.mean,col=i)
}
legend("bottomright", lty=1, col=1:length(PopSize), legend=paste0("N=", PopSize))
```

In this graph you can see the two phases: one where the response to selection is strong, and one where the optimum is reached (around 30 generations at N=1000) and the response reaches a plateau. The smallest population (N=10) displays a response to selection that is dimished compared with the larger populations. In this simulation at N=10, the optimum is never reached. In addition, the variance around the response is important at N=10, less so at N=100, and there is little fluctuation at N=1000. This is partly due to sampling (since the variance is calculated using a smaller number of individuals), but even if sampling were to be corrected of (i.e. variance measured for 10 individuals for all values of N) this pattern would still be observable. You can run more simulations to verify this pattern.

Why is the response dimished at N=10? Let's look at the genetic variance.

```{r}
plot(NA,xlim=c(0,50),ylim=c(0,1.5),xlab="generations", ylab=("gen.var"))
for (i in 1:length(PopSize)) 
  {
  ps<-get(paste("PopSize_",PopSize[i],sep=""))
  lines(ps$gen.var,col=i)
}
legend("topright", lty=1, col=1:length(PopSize), legend=paste0("N=", PopSize))
```

You see in this graph that the genetic variance is rapidly exhausted at N=10: after 20 generations there is very little variance left. That is why after 20 generations, this population can no longer respond to selection. The speed of the exhaustion of the genetic variance depends on the population size. 

Note that increasing the selection of strength (sel.Vs = 0.1) depletes the genetic variance rapidly, limiting the response to selection and the optimum of 10 is never reached. The mechanisms behind this are discussed in more detail below.


```{r}
PopSize<-c(10,100,1000)
for (i in 1:length(PopSize)) {
  current_PopSize <- PopSize[i]
  current_name <- paste0("PopSize_",current_PopSize)
  current_simu <- simulation(generations=50,pop.size=current_PopSize,num.loci=100,var.init=1,var.env=1,sel.Vs=0.1,sel.optimum=10)
  assign(current_name,current_simu)
  rm(current_simu)
}

plot(NA,xlim=c(0,50),ylim=c(0,8),xlab="generations", ylab=("phen.mean"))
for (i in 1:length(PopSize)) 
  {
  ps<-get(paste("PopSize_",PopSize[i],sep=""))
  lines(ps$phen.mean,col=i)
}
legend("bottomright", lty=1, col=1:length(PopSize), legend=paste0("N=", PopSize))
```


We now explore the effect of the number of loci. Here we choose 3 values, 1, 50, 100 loci. We consider 50 generations, a population size of 100, and a selection optimum of 10, and sel.Vs = 5.

```{r}
NbLocus<-c(1,10,100)
for (i in 1:length(NbLocus)) {
  current_NbLocus <- NbLocus[i]
  current_name <- paste0("NbLocus_",current_NbLocus)
  current_simu <- simulation(generations=50,pop.size=100,current_NbLocus,var.init=1,var.env=1,sel.Vs=5,sel.optimum=10)
  assign(current_name,current_simu)
  rm(current_simu)
}

plot(NA,xlim=c(0,50),ylim=c(0,12),xlab="generations", ylab=("phen.mean"))
for (i in 1:length(NbLocus)) 
{
  ps<-get(paste("NbLocus_",NbLocus[i],sep=""))
  lines(ps$phen.mean,col=i)
}
legend("bottomright", lty=1, col=1:length(NbLocus), legend=paste0("NbLocus=", NbLocus))
```


The response to selection depends on the number of loci, this is because the genetic variance is rapidly exhausted with only a single locus and the same applies for a number of loci between 10 and 100 as a shown below.


```{r}
plot(NA,xlim=c(0,50),ylim=c(0,1.5),xlab="generations", ylab=("gen.var"))
for (i in 1:length(NbLocus)) 
{
  ps<-get(paste("NbLocus_",NbLocus[i],sep=""))
  lines(ps$gen.var,col=i)
}
legend("bottomright", lty=1, col=1:length(NbLocus), legend=paste0("NbLocus=", NbLocus))
```

With more loci, we approach the infinitesimal model: genetic variance is maintained through time.

Now let's have a look at the initial genetic variance. We again consider 50 generations, 100 individuals, 100 loci, sel.Vs = 5 and a selection optimum of 8. With no genetic variance, we expect no selection response. How does the amount of initial genetic variance impact the time needed to reach the optimum?


```{r}

VarInit<-c(0,1,2,3,4)
for (i in 1:length(VarInit)) {
  current_VarInit <- VarInit[i]
  current_name <- paste0("VarInit_",current_VarInit)
  current_simu <- simulation(generations=50,pop.size=100,num.loci=100,current_VarInit,var.env=1,sel.Vs=5,sel.optimum=8)
  assign(current_name,current_simu)
  rm(current_simu)
}

plot(NA,xlim=c(0,50),ylim=c(0,10),xlab="generations", ylab=("phen.mean"))
for (i in 1:length(VarInit)) 
{
  ps<-get(paste("VarInit_",VarInit[i],sep=""))
  lines(ps$phen.mean,col=i)
}
legend("bottomright", lty=1, col=1:length(VarInit), legend=paste0("VarInit=", VarInit))
```

As predicted, the optimum is never reached with 0 initial variation. It is reached around 40G, 22G, 20G and 10G for increasing values of var.init. The higher the initial variance, the more genetic variance maintained over time (see the graph below). The level of initial variance modifies primarily the slope of the initial response to selection (R = Va * beta; so that with constant selection (beta), the response between the generation 0 and 1 is proportional to Va).


```{r}
plot(NA,xlim=c(0,50),ylim=c(0,6),xlab="generations", ylab=("gen.var"))
for (i in 1:length(VarInit)) 
{
  ps<-get(paste("VarInit_",VarInit[i],sep=""))
  lines(ps$gen.var,col=i)
}
legend("topright", lty=1, col=1:length(VarInit), legend=paste0("VarInit=", VarInit))
```

Let's look at impact of the selection optimum.

```{r}

SelOpt<-c(0,2,4,6,8,10)
for (i in 1:length(SelOpt)) {
  current_SelOpt <- SelOpt[i]
  current_name <- paste0("SelOpt_",current_SelOpt)
  current_simu <- simulation(generations=50,pop.size=100,num.loci=100,var.init=1,var.env=1,sel.Vs=5,current_SelOpt)
  assign(current_name,current_simu)
  rm(current_simu)
}

plot(NA,xlim=c(0,50),ylim=c(0,10),xlab="generations", ylab=("phen.mean"))
for (i in 1:length(SelOpt)) 
{
  ps<-get(paste("SelOpt_",SelOpt[i],sep=""))
  lines(ps$phen.mean,col=i)
}
legend("bottomright", lty=1, col=1:length(SelOpt), legend=paste0("SelOpt=", SelOpt))
```

It takes longer to reach a far optimum: an optimum of 2 here is reached in 20 generations, while it takes 40 generations to reach an optimum of 10. We can predict that with an optimum which is further away the average initial fitness of the population would be lower than with an optimum that is close. Let's verify this. 




```{r}
plot(NA,xlim=c(0,50),ylim=c(0,1),xlab="generations", ylab=("fit.mean"))
for (i in 1:length(SelOpt)) 
{
  ps<-get(paste("SelOpt_",SelOpt[i],sep=""))
  lines(ps$fit.mean,col=i)
}
legend("bottomright", lty=1, col=1:length(SelOpt), legend=paste0("SelOpt=", SelOpt))
```
We indeed find that the initial fitness is lower for an optimum of 10 than for an optimum of 8, 6 etc.. 

Let's look now at the strength of selection through the parameter "sel.Vs".

```{r}

SelStr<-c(Inf,10,5,1,0.1)
for (i in 1:length(SelStr)) {
  current_SelStr <- SelStr[i]
  current_name <- paste0("SelStr_",current_SelStr)
  current_simu <- simulation(generations=50,pop.size=100,num.loci=100,var.init=1,var.env=1,current_SelStr,sel.optimum=10)
  assign(current_name,current_simu)
  rm(current_simu)
}

plot(NA,xlim=c(0,50),ylim=c(0,10),xlab="generations", ylab=("phen.mean"))
for (i in 1:length(SelStr)) 
{
  ps<-get(paste("SelStr_",SelStr[i],sep=""))
  lines(ps$phen.mean,col=i)
}
legend("bottomright", lty=1, col=1:length(SelStr), legend=paste0("SelStr=", SelStr))
```

We see here an interesting pattern where the optimum is reached only for intermediate values of sel.Vs (10 and 5). In absence of selection, the phenotype remains around its initial value and when the selection strength is too high, the best haplotypes present in the first generations are rapidly fixed, and no variability is left to create new haplotypes through recombination. In other words, we have here a trade off between the speed of the response and the capacity of the population to reach the optimum: with strong selection strength, the response is fast but limited; with mild selection strength, the response is slow but the population eventually reaches the optimum. The dynamics of the exhaustion of the genetic variance are illustrated below. 



```{r}
plot(NA,xlim=c(0,50),ylim=c(0,2),xlab="generations", ylab=("gen.var"))
for (i in 1:length(SelStr)) 
{
  ps<-get(paste("SelStr_",SelStr[i],sep=""))
  lines(ps$gen.var,col=i)
}
legend("bottomright", lty=1, col=1:length(SelStr), legend=paste0("SelStr=", SelStr))

```

As predicted we observe no genetic variance for strong selection (sel.Vs = 1 and 0.1), while variance is maintained for all other strengths of selection.






















