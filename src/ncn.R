###NCN = Node Contribution to Nestedenss
###MKLau 2 Jul 2014

###Setup simulations
nsim <- 100
require(bipartite)

###Load observed matrix
x <- dget('../data/soli.mod')
x[x!=0] <- 1
x <- x[apply(x,1,sum)!=0,apply(x,2,sum)!=0]

###Observed nestedness
N <- nestedtemp(x)$statistic

###Get multiple estimates of nestedness (Nstar) for each species after randomizing interactions

##For each species shuffle their interactions
##For each observation, shuffle their interactions
##Calculate nestedness
##sites
Nstar.i <- list()
for (i in 1:nrow(x)){
out <- 0
  for (k in 1:nsim){
    print(paste(i,k))
    null <- x
    y <- commsimulator(x,method='r1',thin=100)
    null[i,] <- y[i,]
    out[k] <- nestedtemp(null)$statistic
  }
  Nstar.i[[i]] <- out
}
##species
Nstar.j <- list()
for (j in 1:ncol(x)){
out <- 0
  for (k in 1:nsim){
    print(paste(j,k))
    null <- x
    y <- commsimulator(x,method='r1',thin=100)
    null[,j] <- y[,j]
    out[k] <- nestedtemp(null)$statistic
  }
  Nstar.j[[j]] <- out
}

###
Nstar.i.mu <- unlist(lapply(Nstar.i,mean))
Nstar.i.sd <- unlist(lapply(Nstar.i,sd))
Nstar.j.mu <- unlist(lapply(Nstar.j,mean))
Nstar.j.sd <- unlist(lapply(Nstar.j,sd))

###Calculate contributions
ci <- (N - Nstar.i.mu)/Nstar.i.sd #sites
cj <- (N - Nstar.j.mu)/Nstar.j.sd #species
names(ci) <- rownames(x)
names(cj) <- colnames(x)

###Write output
dput(list(obs=N,ci=ci,cj=cj),file='../data/ncn_out.rda')
