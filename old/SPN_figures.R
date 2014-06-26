###################################################
### chunk number 1: 
###################################################
#line 54 "SPN_figures.Rnw"
require(ecodist)
require(sna)
require(vegan)
require(bipartite)
source('/Users/Aeolus/Documents/Active_Projects/CorNets/CorNets.R')


###################################################
### chunk number 2: 
###################################################
#line 62 "SPN_figures.Rnw"
bin.sum = function(x){x[x!=0]=1;sum(x)}


###################################################
### chunk number 3: 
###################################################
#line 68 "SPN_figures.Rnw"
data = read.csv('/Users/Aeolus/Documents/Active_Projects/SolidagoPollinationNetwork/SolidagoPollinators2010.csv')

summary(data)
colnames(data)

com = data[14:ncol(data)]
##remove singletons
com = com[,apply(com,2,sum)>1]
colnames(com)
ds = rep(1,nrow(com)) #dummy species
com. = cbind(com,ds) #Adjusted Bray-Curtis (Clarke et al. 2006)

##population factor
pop = data[,4]
pop

##family factor
fam = data[,1]
###FIX FAMILY TYPE THROPE == THORPE
fam[fam == 'Thrope 5'] <- 'Thorpe 5'

##elk factor
elk = data[,5]
elk[elk=='']='No'#make all empty elk entries No's




###################################################
### chunk number 4: 
###################################################
#line 99 "SPN_figures.Rnw"
##Network Modeling
###JUST USE ELK AS A FACTOR
com.list <- list()
for (i in (1:length(unique(elk)))){
com.list[[i]] <- com[elk == unique(elk)[i],]
}

names(com.list) <- unique(elk)

net.list <- lapply(com.list,kendall.pairs,adj.method='fdr',alpha=0.05,p.adj=FALSE)



###################################################
### chunk number 5: 
###################################################
#line 116 "SPN_figures.Rnw"
  #Graphs for dave
quartz('',15,9.5)
par(mfrow=c(1,2),mar=c(2.4,1.3,1.3,1.3),oma=c(0.1,0.1,0.1,0.1),cex=2,mar=c(2,1,1,1))
names(net.list)=c('No Elk','Elk') #rename the network graphs
net.list.reorder <- net.list[c(2,1)]
com.list.reorder <- com.list[c(2,1)]
for (i in (1:length(net.list))){
      v.col=apply(com.list.reorder[[i]],2,sum); v.col[v.col != 0] = 'black' #color the present species
      v.col[v.col == 0] <- 'lightgray' #color empty species
      gplot(abs(net.list.reorder[[i]]),gmode='graph',vertex.cex=3,vertex.sides=100,vertex.col=v.col,edge.lwd=(abs(net.list.reorder[[i]])+1)^5,edge.col=gray(0.1)[1],vertex.border='grey',mode='circle',displaylabels=FALSE,cex=2,main=names(net.list.reorder)[i]) #without titles
    }



