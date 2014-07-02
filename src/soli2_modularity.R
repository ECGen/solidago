###Testing the modularity of the solidago network
###MKLau 2 Jul 2014

#load packages
null.gen <- FALSE
nsim <- 1000
require(bipartite)
require(methods)

#load model 
x <- dget('../data/soli.mod')

#make binary
x[x!=0] <- 1

#compute observed modules
obs <- computeModules(x)

#generate null communities
if (null.gen){
  null.set <- list()
  for (i in 1:nsim){
    print(i)
    null.set[[i]] <- commsimulator(x,method='r1',thin=100)
    null.set[[i]] <- null.set[[i]][apply(null.set[[i]],1,sum)!=0,apply(null.set[[i]],2,sum)!=0]
  }
  dput(null.set,'../data/null_set.rda')
}else{
  null.set <- dget('../data/null_set.rda')
} 

###Get modularity for null communities
null <- lapply(null.set,function(x) slot(computeModules(x),name='likelihood'))
out <- c(obs=slot(obs,'likelihood'),unlist(null))

###Write output
dput(obs,'../data/mod_observed.rda')
dput(out,'../data/mod_like_obs_nul.rda')
