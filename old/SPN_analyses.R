##Solidago Pollinator Network Analyses
##Project Director: Dave Smith
##Data Recorder: Ryan ? (undergrad assistant)

##NOTE FROM DAVE ABOUT PREVIOUS ANALYSES: I did not include population FS in my analysis for elk / no-elk.  Although, they could be useful data for looking at heritability or perhaps something else.  Also, for my analysis, I only included families with more than one rep (I also did an analysis including only families with 3 or more reps).  And, I removed singletons from the data.

##Questions
##1. do more genetically similar individuals tend to have similar communities?
##2. does this pattern produce a characteristic pollinator network structure?
##3. do do more similar inidividuals tend to have different pollinator communities in order to avoid in-breeding?

require(ecodist)
require(sna)
require(vegan)
require(bipartite)
library(pbapply)
library(sna)
source('/Users/Aeolus/Documents/Active_Projects/CorNets/CorNets.R')
source('~/cor_nets/araujo_method/araujo_method.R')
bin.sum = function(x){x[x!=0]=1;sum(x)}
                                        #import data
data = read.csv('/Users/Aeolus/Documents/Active_Projects/SolidagoPollinationNetwork/SolidagoPollinators2010.csv')
summary(data)
colnames(data)
                                        #remove FS
data <- data[data[,4]!='FS',]
                                        #remove families with less than 1 rep
data <- data[data[,1]%in%names(table(data[,1]))[table(data[,1])>1],]
                                        #separate
com = data[14:ncol(data)]
                                        #remove singletons
com = com[,apply(com,2,sum)>1]
colnames(com)
ds = rep(1,nrow(com)) #dummy species
com. = cbind(com,ds) #Adjusted Bray-Curtis (Clarke et al. 2006)

##population factor
pop = factor(as.character(data[,4]))
pop

##family factor
fam <- factor(as.character(data[,1]))

###FIX FAMILY TYPE THROPE == THORPE
fam[fam == 'Thrope 5'] <- 'Thorpe 5'

##elk factor
elk = data[,5]
elk[elk=='']='No'#make all empty elk entries No's
elk <- factor(as.character(elk))

###Compositional Analyses
adonis(com.~pop%in%elk)
com.. <- apply(com.,2,function(x) x/max(x))
adonis(com..~pop%in%elk)

##Network Modeling
###JUST USE ELK AS A FACTOR
com.list <- split(com,elk)
                                        #co-occurrence patterns
coa.elk <- lapply(split(com,elk),CA.results,nmc=1000)
coa.pop <- lapply(split(com,pop),CA.results,nmc=1000)
do.call(rbind,lapply(coa.elk,function(x) x[[2]]))
do.call(rbind,lapply(coa.pop,function(x) x[[2]]))
                                        #network models
                                        #gams
nsp <- 5
library(gam)
gam.net <- lapply(com.list,function(com) gamNet(com[,apply(com,2,sum)>nsp]))
detach(package:gam)
no.v.cex <- com[elk=='No',]
no.v.cex <- apply(no.v.cex[,apply(no.v.cex,2,sum)>nsp],2,sum)
no.v.cex <- (no.v.cex/sum(no.v.cex)+1)^2
yes.v.cex <- com[elk=='Yes',]
yes.v.cex <- apply(yes.v.cex[,apply(yes.v.cex,2,sum)>nsp],2,sum)
yes.v.cex <- (yes.v.cex/sum(yes.v.cex)+1)^2
no.e.col <- gam.net$No
no.e.col[no.e.col<0] <- 'black'
no.e.col[no.e.col>0] <- 'grey'
yes.e.col <- gam.net$Yes
yes.e.col[yes.e.col<0] <- 'black'
yes.e.col[yes.e.col>0] <- 'grey'
                                        #plots
par(mfrow=c(1,2))
gplot(gam.net$No,displaylabels=TRUE,label.cex=0.5,gmode='graph',vertex.sides=50,edge.lwd=(gam.net$No)^2,vertex.cex=no.v.cex,main='Absent',edge.col=no.e.col,mode='circle')
gplot(gam.net$Yes,displaylabels=TRUE,label.cex=0.5,gmode='graph',vertex.sides=50,edge.lwd=(gam.net$Yes)^2,vertex.cex=yes.v.cex,main='Present',edge.col=yes.e.col,mode='circle')
                                        #structure
gam.d <- lapply(gam.net,function(x) ncol(x[apply(x,1,sum)!=0,apply(x,2,sum)!=0])) #degree
gam.tnc <- lapply(gam.net,function(x) sum(apply(x[apply(x,1,sum)!=0,apply(x,2,sum)!=0],1,bin.sum))) # total number of connections
gam.anc  <- lapply(gam.net,function(x) mean(apply(x[apply(x,1,sum)!=0,apply(x,2,sum)!=0],1,bin.sum))) # average number of connections
gam.c <- lapply(gam.net,function(x) centralization(x,degree)) #centralization
str.sum <- do.call(rbind,list(gam.d,gam.tnc,gam.anc,gam.c)) # structure summary
rownames(str.sum) <- c('Degree','Total Connections','Avg. Connections','Centrality')
str.sum
                                        #kendall
net.list <- lapply(com.list,function(x) kendall.pairs(x[,apply(x,2,sum)>nsp],adj.method='fdr',alpha=0.05,p.adj=FALSE))
                                        #network graph
par(mfrow=c(1,2))
gplot(net.list$No,displaylabels=TRUE,label.cex=0.75,gmode='graph',vertex.sides=50,mode='circle')
gplot(net.list$Yes,displaylabels=TRUE,label.cex=0.75,gmode='graph',vertex.sides=50,mode='circle')
                                        #araujo co-occurrence networks
com.elk <- split(com,elk)
elk.nets <- lapply(com.elk,araujoNet,min.abundance=nsp)
                                        #plots
par(mfrow=c(1,3))
gplot(araujoNet(com,min.abundance=nsp)$dp,displaylabels=TRUE,label.cex=0.65)
gplot(elk.nets[[1]]$dp,displaylabels=TRUE,label.cex=0.65)
gplot(elk.nets[[2]]$dp,displaylabels=TRUE,label.cex=0.65)
