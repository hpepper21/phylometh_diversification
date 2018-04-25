getwd()
setwd("/Users/elipepper/Documents/GitHub")
#install.packages(c("ape", "TreeSim", "geiger", "diversitree", "devtools"))
library(ape)
library(TreeSim)
library(geiger)
library(diversitree)
devtools::install_github("thej022214/hisse")
library(hisse)
#simulating 30 taxon tree w/ extinction rate = 0
my.tree <- TreeSim::sim.bd.taxa(n=300, numbsim=1, lambda=0.1, mu=0)[[1]]
#plotting my.tree
plot(my.tree) #plot of phylogeny
ape::ltt.plot(my.tree) #species-time plot w/ exponential growth
#putting species-time plot on log scale
ape::ltt.plot(my.tree, log="y")
#looking at multiple trees
yule.trees <- TreeSim::sim.bd.taxa(n=300, numbsim=10, lambda=0.1, mu=0, complete=FALSE)
#plotting multiple trees
mltt.plot(yule.trees,legend = F)
#look at trees w/ birth & death rate
bd.trees <- TreeSim::sim.bd.taxa(n=300, numbsim=10, lambda=1, mu=.9, complete=FALSE)
ape::mltt.plot(bd.trees, log="y", legend=FALSE)
#comparing B-D trees w/ yule trees
depth.range <- range(unlist(lapply(yule.trees,ape::branching.times)), unlist(lapply(bd.trees,ape::branching.times)))
max.depth <- sum(abs(depth.range)) #ape rescales depths
plot(x=c(0, -1*max.depth), y=c(1, ape::Ntip(yule.trees[[1]])), log="y", type="n", bty="n", xlab="Time", ylab="N")
colors=c(rgb(1,0,0,0.5), rgb(0, 0, 0, 0.5))
list.of.both <- list(bd.trees, yule.trees)
for (i in sequence(2)) {
  tree.list <- list.of.both[[i]]
  for (j in sequence(length(tree.list))) {
    ape::ltt.lines(tree.list[[j]], col=colors[[i]])   
  }
}
legend("topleft", legend=c("Birth Death", "Yule"),cex = .5, fill=colors)
#zooming in on final part of the plot
depth.range <- range(unlist(lapply(yule.trees,ape::branching.times)), unlist(lapply(bd.trees,ape::branching.times)))
max.depth <- sum(abs(depth.range)) #ape rescales depths
plot(x=c(0, -5), y=c(200, ape::Ntip(yule.trees[[1]])), log="y", type="n", bty="n", xlab="Time", ylab="N")
colors=c(rgb(1,0,0,0.5), rgb(0, 0, 0, 0.5))
list.of.both <- list(bd.trees, yule.trees)
for (i in sequence(2)) {
  tree.list <- list.of.both[[i]]
  for (j in sequence(length(tree.list))) {
    ape::ltt.lines(tree.list[[j]], col=colors[[i]])   
  }
}
legend("topleft", legend=c("Birth Death", "Yule"),cex = .5, fill=colors)

###plotting trees w/ other diversification parameters
#what happens if speciation rate is much higher than extinction rate?
my.trees <- TreeSim::sim.bd.taxa(n=300, numbsim=10, lambda=1, mu=0, complete=FALSE)
ape::mltt.plot(my.trees, log="y", legend=FALSE)
my.trees1 <- TreeSim::sim.bd.taxa(n=300, numbsim=10, lambda=0.8, mu=0.2, complete=FALSE)
ape::mltt.plot(my.trees1, log="y", legend=FALSE)
my.trees2 <- TreeSim::sim.bd.taxa(n=300, numbsim=10, lambda=1, mu=0.4, complete=FALSE)
ape::mltt.plot(my.trees2, log="y", legend=FALSE)

###now going through HiSSE example
#simulation of tree and characters
speciation.rates <- c(0.1, 0.1, 0.1, 0.2) #0A, 1A, 0B, 1B
extinction.rates <- rep(0.03, 4)
transition.rates <- c(0.01,0.01,0, 0.01, 0, 0.01, 0.01,0,0.01, 0,0.01,0.01)
pars <- c(speciation.rates, extinction.rates, transition.rates)
phy <- tree.musse(pars, max.taxa=50, x0=1, include.extinct=FALSE)
sim.dat.true <- data.frame(names(phy$tip.state), phy$tip.state)
sim.dat <- sim.dat.true
# Now to hide the "hidden" state
sim.dat[sim.dat[,2]==3,2] = 1
sim.dat[sim.dat[,2]==4,2] = 2
# and convert states 1,2 to 0,1
sim.dat[,2] = sim.dat[,2] - 1
#now view plot
plot(phy)
knitr::kable(cbind(sim.dat, true.char=sim.dat.true$phy.tip.state))

#hisse examples
turnover.anc = c(1,1,0,0)
#turnover.anc has a single free parameter for both 0A and 1A state combos.;
#there is a single free parameter for extinction fraction...
#this is equivalent to a bisse model w/ a fixed turnover and extinction rates
#across the observed states 0 and 1
eps.anc = c(1,1,0,0)
#now, say we want to include separate turnover rates for both states 0A and 1A
turnover.anc = c(1,2,0,0)
#thus, a full hisse model would be 
turnover.anc = c(1,2,3,4)
#this corresponds to four separate net turnover rates for 1=0A, 2=1A, 3=0B, 4=1B.
#extinction fraction follows the same format, though including a zero for a state
#we want to include in the model corresponds to no extinction, which is the yule (pure birth) equivalent
eps.anc = c(0,0,0,0)

###setting up the transistion rate matrix
#generate the index matrix describing the free parameters in the transition model:
trans.rates = TransMatMaker(hidden.states=TRUE)
trans.rates
#remove the dual transitions b/w both the observed trait and the hidden trait from the model entirely:
trans.rates.nodual = ParDrop(trans.rates, c(3,5,8,10))
trans.rates.nodual

#run a model where we assume all transitions are equal to one another;
#using the example above, let's set parameter 1 and 6 to have the same rate:
trans.rates.nodual.equal16 = ParEqual(trans.rates.nodual, c(1,6))
trans.rates.nodual.equal16
#note that the rate for parameter 6 and the rate for parameter 1 will take the same value

#now setting all rates to be equal
trans.rates.nodual.allequal <- ParEqual(trans.rates.nodual, c(1,2,1,3,1,4,1,5,1,6,1,7,1,8))
trans.rates.nodual.allequal
#there are less confusing ways to get the same result - e.g.,
trans.rates.nodual.allequal <- trans.rates.nodual
trans.rates.nodual.allequal[!is.na(trans.rates.nodual.allequal) & !trans.rates.nodual.allequal == 0] = 1
trans.rates.nodual.allequal

#note that in order to run a BiSSE mod in HiSSE, the matrix would look like this:
trans.rates.bisse <- TransMatMaker(hidden.states = FALSE)
trans.rates.bisse

#input the transition matrix to the trans.rate = argument w/in the hisse() call
pp <- hisse(phy = phy, data = sim.dat, f=c(1,1), hidden.states = TRUE, turnover.anc = turnover.anc, eps.anc = eps.anc, trans.rate = trans.rates.nodual.allequal)

#example of a common mistake; it may be of interest to test a model
#where the hidden state is associated w only a single observed state,
#such that the model contains states 0A, 1A, & 1B
#diversification parameters might look something like this...
turnover.anc <- c(1,2,0,3)
eps.anc <- c(1,2,0,3)
#the 0 in the 3rd entry for state 0B designates that the parameter is 
#removed from the model;
#a common mistake is that the transitions to & from 0B are not removed
#from the transition matrix (this must be done manually):
trans.rates <- TransMatMaker(hidden.states = TRUE)
trans.rates.nodual.no0B <- ParDrop(trans.rates, c(2,3,5,7,8,9,10,12))
trans.rates.nodual.no0B

#the output.type = argument alters how the final parameters are printed after running hisse();
#there are 3 options for output: "turnover", "net.div", & "raw"; 
#whichever option you choose returns the results as estimates of 
#speciation (lambda) & extinction (mu);
#now, outputting net diversification:
pp <- hisse(phy = phy, sim.dat, f=c(1,1), hidden.states = TRUE, turnover.anc = turnover.anc, eps.anc = eps.anc, trans.rate = trans.rates.nodual.allequal, output.type = "net.div")

#the goal is to set up a mod where the diversification process is independent
#from the observed states (0 or 1) of the focal trait;
#diversification rates, if they exist, will only be associated w one of the
#hidden states (A or B) regardless of the state of the focal trait;
#the free parameters for this diversification model should look like this:
turnover.anc <- c(1,1,2,2)
eps.anc <- c(1,1,2,2)
#we are specifying that both 0A & 1A have one set of diversification rates,
#& 0B & 1B have another set of rates; this is the "null-two" mod.

#there are 3 ways to set up transition rates...
#the first is to assume the usual 8 transitions in the full hisse mod
trans.rates <- TransMatMaker(hidden.states = TRUE)
trans.rates.nodual <- ParDrop(trans.rates, c(3,5,8,10))
#we could also assume all rates are equal:
trans.rates.nodual.allequal <- ParEqual(trans.rates.nodual, c(1,2,1,3,1,4,1,5,1,6,1,7,1,8))
trans.rates.nodual.allequal
#a third option specifies 3 rates: one rate describing transitionss among the
#different hidden states (which could be interpreted as the rate by which shifts in diversification occur),
#& two rates for transitions b/w the observed character states:
#now we want 3 specific rates:
trans.rates.nodual.threerates <- trans.rates.nodual
#set all transitions from 0->1 to be governed by a single rate:
to.change <- cbind(c(1,3), c(2,4))
trans.rates.nodual.threerates[to.change] = 1
#now set all transitions from 1->0 to be governed by a single rate:
to.change <- cbind(c(2,4), c(1,3))
trans.rates.nodual.threerates[to.change] = 2
#finally, set all transitions b/w the hidden state to be a single rate,
#this essentially gives you an estimate of the rate by which shifts in
#diversification occur:
to.change <- cbind(c(1,3,2,4), c(3,1,4,2))
trans.rates.nodual.threerates[to.change] = 3
trans.rates.nodual.threerates

#now inputting these arugments into hisse():
pp <- hisse(phy = phy, sim.dat, f=c(1,1), hidden.states = TRUE, turnover.anc = turnover.anc, eps.anc = eps.anc, trans.rate = trans.rates.nodual.allequal)

#plotting hisse reconstructions:
knitr::opts_chunk$set(fig.width = 7, fig.height = 5)
#for this example, the mod assumes two diversificaiton rate parameters:
#turnover.anc=c(1,1,1,2) & eps.anc=c(1,1,1,1)...
load("testrecon1.Rsave")
class(pp.recon)
pp.recon
#plotting net diversification rates:
plot.hisse.states(pp.recon, rate.param = "net.div", show.tip.label = FALSE)
#you have to put contraints on the range, so that all models will have similar rates when comparing;
#a vector w/ the min & max rate across all mods can be passed to the visual:
plot.hisse.states(pp.recon, rate.param = "net.div", show.tip.label = FALSE, rate.range = c(0,0.072))

###modeling-averaging approach:
#first step is to make sure that the hisse.states objects contain the AIC from
#the mod fit embedded in it:
pp.recon$aic
#the AIC for the mod can be supplied as an argument in the MarginRecon() func.:
pp.recon <- MarginRecon(phy = phy, sim.dat, f=c(1,1), hidden.states = TRUE, pars = pp$solution, aic = pp$AIC, n.cores = 2)
#we are adding two additional hisse.states objects to the modeling-averaging approach;
#one reconstruction is based on the nul-two mod [turnover.anc=c(1,1,2,2)];
#the other assumes four free turnover rates [turnover.anc=c(1,2,3,4)];
#in all cases we assume equal transition rates & equal extinction fractions:
hisse.results.list = list()
load("testrecon1.Rsave")
hisse.results.list[[1]] = pp.recon
load("testrecon2.Rsave")
hisse.results.list[[2]] = pp.recon
load("testrecon3.Rsave")
#now supply the list the plotting func.
plot.hisse.states(hisse.results.list, rate.param = "net.div", show.tip.label = FALSE, rate.range = c(0,0.072))
#note that there is an even easier way to generate the list
#first, suck in all the files w .Rsave line ending in your working directory:
files <- system("ls -1 | grep .Rsave", intern = TRUE)
#create an empty list object:
hisse.results.list <- list()
#now loop through all files, adding the embedded pp.recon object in each
for(i in sequence(length(files))){
  load(files[i])
  hisse.results.list[[i]] = pp.recon
  rm(pp.recon)
}

###now run w my own data - using discrete dat on Gekkonids
##w/ 0 being diurnal & 1 being nocturnal
library(datelife)
gecko.trees <- datelife_search(input="Gekkonidae", get_spp_from_taxon=TRUE)
gek.tree <- gecko.trees[[which.max(lapply(gecko.trees, Ntip))]]
my.dat<-read.delim("BayesTraits_binary.data.txt", header=FALSE, stringsAsFactors=FALSE)
colnames(my.dat)<-c("taxon", "nocturnal")
CleanData<-function(phy,data){
  data.vector <- data$nocturnal
  names(data.vector) <- gsub("_", " ", data$taxon)
  
  cleaned <- treedata(phy, data.vector)
  return(cleaned)
}
clean.discrete<-CleanData(phy = gek.tree,data = my.dat)
pruned.dat <- clean.discrete$data
clean.phy<-clean.discrete$phy
VizualizeData<-function(phy,data){
  plot(phy)
  print(data)
}
VizualizeData(phy = clean.phy,data = pruned.dat)
treedata(clean.phy, pruned.dat)

gek.turnover.anc = c(1,2,3,4)
gek.eps.anc = c(1,2,3,4)

trans.rates <- TransMatMaker(hidden.states = TRUE)
trans.rates.nodual <- ParDrop(trans.rates, c(3,5,8,10))

trans.rates.nodual.threerates <- trans.rates.nodual

to.change <- cbind(c(1,3), c(2,4))
trans.rates.nodual.threerates[to.change] = 1

to.change <- cbind(c(2,4), c(1,3))
trans.rates.nodual.threerates[to.change] = 2

to.change <- cbind(c(1,3,2,4), c(3,1,4,2))
trans.rates.nodual.threerates[to.change] = 3
trans.rates.nodual.threerates

my.data <- data.frame(species=rownames(pruned.dat), trait=pruned.dat[,1], stringsAsFactors = FALSE)

gek.pp <- hisse(phy = clean.phy, my.data, f=c(1,1), hidden.states = TRUE, turnover.anc = gek.turnover.anc, eps.anc = gek.eps.anc, trans.rate = trans.rates.nodual.threerates)




