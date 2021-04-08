#R code for general nest category reconstructions
#data collected by Johnathan Hughes
#code by Jacob Berv, jsb439@cornell.edu
#heavily inspired by Liam Revell's phytools examples


#load requried packages
require(phytools)
require(plotrix)
require(ape)
require(stringr)
require(ratematrix)
require(parallel)
require(doSNOW)
require(doParallel)
require(devtools)
require(pbmcapply)
require(geiger)

## Revert to 'sequential' setup of PSOCK cluster in RStudio Console on macOS and R 4.0.0
if (Sys.getenv("RSTUDIO") == "1" && !nzchar(Sys.getenv("RSTUDIO_TERM")) && 
    Sys.info()["sysname"] == "Darwin" && getRversion() >= "4.0.0") {
  parallel:::setDefaultClusterOptions(setup_strategy = "sequential")
}


#set working directory
setwd("/Users/cotinga/jsb439@cornell.edu/Code/mammal_arboreality")

#load functions for SIMMAP in parallel
source("simmap_parallel.R")
source("match_phylo_nodes.R")
source("rate_through_time.R")
source("helpers.R")

#load the tree and clean up
TimeTree<-read.tree("upham_2019.tre")
TimeTree<-ladderize(TimeTree)

#relabel tips
genera<-as.character(as.data.frame(str_split_fixed(TimeTree$tip.label, "_",n=2))[,1])
TimeTree$tip.label<-genera

#read in posterior trees
# posterior_trees<-read.nexus(file="./tree-pruner-e1a1cb23-d757-4c59-b4cf-20ee0691a4ba/output.nex")
# posterior_trees.backup<-posterior_trees
# set.seed(seed=5)
# posterior_trees.sub<-sample(posterior_trees, size=1000)
# 
# #relabel tips
# for(i in 1:length(posterior_trees.sub)){
#   genera<-as.character(as.data.frame(str_split_fixed(posterior_trees.sub[[i]]$tip.label, "_",n=2))[,1])
#   posterior_trees.sub[[i]]$tip.label<-genera
# }
# 
# saveRDS(posterior_trees.sub, file="posterior_trees.sub.RDS")
posterior_trees.sub <- readRDS("posterior_trees.sub.RDS")

#import the data table
data<-read.table("mammal_data_v2_upham.txt", header=T)
rownames(data)<-as.character(as.data.frame(str_split_fixed(data$Taxon, "_",n=2))[,1])

#reorder the data frame to match the tip order
data<-data[TimeTree$tip.label,]

#remove the outgroups
#remove<-c('Danio', 'Xenopus', 'Anolis', 'Gallus', 'Taeniopygia')

#prune outgroups from datset
#data <-data[!data$Taxon %in% remove,]

#prune outgroups from tree
#TimeTree<-drop.tip(TimeTree, tip=remove)


#generating a subset vector for ASR #version 1 of the dataset
character1<-data$Nesting1
character1<-as.character(character1)
names(character1)<-rownames(data)

#generating a subset vector for ASR #version 2 of the dataset
character2<-data$Nesting2
character2<-as.character(character2)
names(character2)<-rownames(data)

#as.data.frame(character1==character2)


#begin ancestral state reconstructions
#begin parsimony reconstruction using phangorn
require(phangorn)

tmp<-c(as.character(character1))
names(tmp)<-rownames(data)

#set up the character data as phyDat object
character.phyDat<-phyDat(tmp, type="USER", levels=c('Arboreal','Non-arboreal','Semi-arboreal'))

#perform maximum parsimony reconstruction
pars<-ancestral.pars(TimeTree, character.phyDat, type='MPR')

#calculate consistency index using Nick Matzke's code
source("http://ib.berkeley.edu/courses/ib200b/scripts/_R_tree_functions_v1.R")
consistency <- CI_phyd(TimeTree, character.phyDat)

#calculate retention index using phangorn
retention <- RI(TimeTree, character.phyDat)

#calculate parsimony score (tree length)
score<-parsimony(TimeTree, character.phyDat)

#extract tip states
tips <- data.frame(matrix(unlist(pars[1:164]), nrow=length(pars[1:164]), byrow=T))
#fix column names
colnames(tips)<-c('Arboreal','Non-arboreal','Semi-arboreal')
#fix row names
rownames(tips) <- c(seq(1,164))
#convert to matrix
tips<-as.matrix(tips)

#extract reconstructed nodes
pars <- data.frame(matrix(unlist(pars[165:327]), nrow=length(pars[165:327]), byrow=T))
#fix column names
colnames(pars)<-c('Arboreal','Non-arboreal','Semi-arboreal')
#fix row names
rownames(pars) <- c(seq(165:327))
#convert to matrix
pars<-as.matrix(pars)


#setting up custom transition model
#one rate for transitions between arboreal and non-arboreal
#one rate for transitions that go through an intermediate state

transitions <- matrix(c(0, 1, 2, 
                        1, 0, 2, 
                        2, 2, 0), nrow=3)

#Rate index matrix:
#             Arboreal Non-arboreal Semi-arboreal
#Arboreal             .            1             2
#Non-arboreal         1            .             2
#Semi-arboreal        2            2             .


#setting up alternative custom transition model
#forward and reverse rates for non arboreal to semi-arboreal
#forward and revers rates for arboreal to semi-arboreal
#ie, forcing transitions to go through semi-arboreal

transitions.int <- matrix(c(0, 0, 1, 
                            0, 0, 2, 
                            3, 4, 0), nrow=3)


#running everything for character1 first
#fit with fitMk and with ace for comparison

#first, fit the first custom model, called "transitions"
fitCustom.character1 <- fitMk(tree = TimeTree, x = character1, model = transitions, pi="fitzjohn")
fitCustom.character1.ace <- ace(character1, TimeTree, model=transitions, type='discrete')
AIC(fitCustom.character1) #AIC = 244.2246
AIC(fitCustom.character1.ace) #AIC = 243.5212

#next, fit the second custom model, called "transitions.int"
fitCustom.int.character1 <- fitMk(tree = TimeTree, x = character1, model = transitions.int, pi="fitzjohn")
fitCustom.int.character1.ace <- ace(character1, TimeTree, model=transitions.int, type='discrete')
AIC(fitCustom.int.character1) #AIC = 234.3019
AIC(fitCustom.int.character1.ace) #AIC = 232.6041

#next, fit the ARD model
fitARD.character1 <- fitMk(tree = TimeTree, x = character1, model = "ARD", pi="fitzjohn")
fitARD.character1.ace <- ace(character1, TimeTree, model='ARD', type='discrete')
AIC(fitARD.character1) #AIC = 238.6664
AIC(fitARD.character1.ace) #AIC = 236.535

#Hidden Markov models with corHMM
#following the 2.1 vignette

require(corHMM)
#generate data frame appropriate for corHMM
character1.corhmm <- data.frame(sp=names(character1), d=character1)

# how many of each state do we have?
summary(as.factor(character1.corhmm[, 2])) 
#get legend and rat matrix from dataset
character1.LegendAndRate<-getStateMat4Dat(character1.corhmm)

#$legend
#              1               2               3 
#     "Arboreal"  "Non-arboreal" "Semi-arboreal"

#setting up 2 rate classes in this test

#rate class one is equivalent to the transition.int model from above
character1.corhmm_R1 <- dropStateMatPars(character1.LegendAndRate$rate.mat, c(1, 3))

#rate class two is equivalent to the ARD model from above
character1.corhmm_R2 <- character1.LegendAndRate$rate.mat

#combine the two rate classes into a list
character1.corhmm.ObsStateClasses <- list(character1.corhmm_R1, character1.corhmm_R2)

#assume transitions between rate classes occur at the same rate
character1.corhmm.RateClassMat <- getRateCatMat(2)
MFT_RateClassMat <- equateStateMatPars(character1.corhmm.RateClassMat, 1:4)

#generate corHMM compatible rate matrix
character1.corhmm_FullMat <- getFullMat(character1.corhmm.ObsStateClasses, character1.corhmm.RateClassMat)

#plotting to check it makes sense
#plotMKmodel(character1.corhmm_FullMat, rate.cat = 2, display = "square", text.scale = 0.9)

#running corHMM on the custom rate matrix
character1.corhmm.out <- corHMM(phy = TimeTree, data = character1.corhmm, rate.cat = 2, rate.mat = character1.corhmm_FullMat, node.states = "none", root.p="maddfitz")
character1.corhmm.out$AIC #AIC =  246.2115

#trying with ARD + 2 hidden categories
character1.corhmm.out.2 <- corHMM(phy = TimeTree, data = character1.corhmm, rate.cat = 2, node.states = "none", root.p="maddfitz")
character1.corhmm.out.2$AIC #AIC = 250.5569

#trying with ARD + 3 hidden categories
character1.corhmm.out.3 <- corHMM(phy = TimeTree, data = character1.corhmm, rate.cat = 3, node.states = "none", root.p="maddfitz")
character1.corhmm.out.3$AIC #AIC =  270.3725


#collect the AIC values to generate a table
AIC(fitCustom.character1)
AIC(fitCustom.int.character1)
AIC(fitARD.character1)
character1.corhmm.out$AIC
character1.corhmm.out.2$AIC
character1.corhmm.out.3$AIC


#stochastic character1 map with custom models

#this next line estimates 5000 maps on the TimeTree using the first custom model
TimeTree.simmap.character1 <- simmap.parallel(trees=TimeTree, model=transitions, sims=5000, characters=character1, ladderize=F, pi="fitzjohn", Q_cores=10, S_cores=10)
saveRDS(TimeTree.simmap.character1, file="TimeTree.simmap.character1.RDS")
#TimeTree.simmap.character1 <-readRDS("TimeTree.simmap.character1.RDS")

#the next line runs the same analysis on the posterior sample of trees
simmap.character1.posterior <- simmap.parallel(trees=posterior_trees.sub, model=transitions, sims=500, characters=character1, ladderize=F, pi="fitzjohn", Q_cores=63, S_cores=20)
saveRDS(simmap.character1.posterior, file="simmap.character1.posterior.RDS")
#simmap.character1.posterior <-readRDS("simmap.character1.posterior.RDS") 

Qs.pos<-Q_extract(simmap.character1.posterior)
saveRDS(Qs.pos, file="Qs.pos.RDS")
#Qs.pos<-readRDS("Qs.pos.RDS")


################

#this next line estimates 5000 maps on the TimeTree using the second custom model

#stochastic character1 map with custom model forcing intermediate states
#TimeTree.simmap.character1.int <- simmap.parallel(trees=TimeTree, model=transitions.int, sims=5000, characters=character1, ladderize=F, pi="fitzjohn", Q_cores=10, S_cores=10)
#saveRDS(TimeTree.simmap.character1.int, file="TimeTree.simmap.character1.int.RDS")
TimeTree.simmap.character1.int <- readRDS("TimeTree.simmap.character1.int.RDS")

#the next line runs the same analysis on the posterior sample of trees
#simmap.character1.posterior.int <- simmap.parallel(trees=posterior_trees.sub, model=transitions.int, sims=500, characters=character1, ladderize=F, pi="fitzjohn", Q_cores=63, S_cores=20)
#saveRDS(simmap.character1.posterior.int, file="simmap.character1.posterior.int.RDS")
simmap.character1.posterior.int<-readRDS("simmap.character1.posterior.int.RDS")

#Qs.pos.int<-Q_extract(simmap.character1.posterior.int)
Qs.pos.int.all<-Q_extract_all(simmap.character1.posterior.int)

test<-countQ_Zeros(simmap.character1.posterior.int)


#saveRDS(Qs.pos.int, file="Qs.pos.int.RDS")
Qs.pos.int<-readRDS("Qs.pos.int.RDS")

################

#stochastic character1 map with ARD model
#TimeTree.simmap.ARD.character1 <- simmap.parallel(trees=TimeTree, model='ARD', sims=5000, characters=character1, ladderize=F, pi="fitzjohn", Q_cores=10, S_cores=10)
#saveRDS(TimeTree.simmap.ARD.character1, file="TimeTree.simmap.ARD.character1.RDS")
TimeTree.simmap.ARD.character1<-readRDS("TimeTree.simmap.ARD.character1.RDS")

#the next line runs the same analysis on the posterior sample of trees
#simmap.ARD.character1.posterior <- simmap.parallel(trees=posterior_trees.sub, model='ARD', sims=500, characters=character1, ladderize=F, pi="fitzjohn", Q_cores=63, S_cores=20)
#saveRDS(simmap.ARD.character1.posterior, file="simmap.ARD.character1.posterior.RDS")
simmap.ARD.character1.posterior <- readRDS("simmap.ARD.character1.posterior.RDS")

#Qs.pos.ARD<-Q_extract(simmap.ARD.character1.posterior)
#saveRDS(Qs.pos.ARD, file="Qs.pos.ARD.RDS")
Qs.pos.ARD<-readRDS("Qs.pos.ARD.RDS")

################


#generate summaries
#posterior summary of stochastic mapping

#summarize stochastic maps for first custom model
pd.1 <- describe.simmap.alt(TimeTree.simmap.character1, plot=F, ref.tree=TimeTree)
saveRDS(pd.1, file="pd.1.RDS")
#pd.1 <- readRDS("pd.1.RDS")

pd.1.pos<-describe.simmap.alt(ref.tree=TimeTree, tree=simmap.character1.posterior)
saveRDS(pd.1.pos, file="pd.1.pos.RDS")
#pd.1.pos <- readRDS("pd.1.pos.RDS")

################

#summarize stochastic maps for the second custom model
pd.int.1 <- describe.simmap.alt(tree=TimeTree.simmap.character1.int, ref.tree=TimeTree)
saveRDS(pd.int.1, file="pd.int.1.RDS")
#pd.int.1 <- readRDS("pd.int.1.RDS")

pd.int.1.pos<-describe.simmap.alt(ref.tree=TimeTree, tree=simmap.character1.posterior.int)
saveRDS(pd.int.1.pos, file="pd.int.1.pos.RDS")
#pd.int.1.pos <- readRDS("pd.int.1.pos.RDS")

################

#summarize stochastic maps for the ARD model
pd.ARD.1 <- describe.simmap.alt(TimeTree.simmap.ARD.character1, plot=F, ref.tree=TimeTree)
saveRDS(pd.ARD.1, file="pd.ARD.1.RDS")
#pd.ARD.1 <- readRDS("pd.ARD.1.RDS")

pd.ARD.1.pos<-describe.simmap.alt(ref.tree=TimeTree, tree=simmap.ARD.character1.posterior)
saveRDS(pd.ARD.1.pos, file="pd.ARD.1.pos.RDS")
#pd.ARD.1.pos<-readRDS("pd.ARD.1.pos.RDS")


################



###setting up plotting

#setting up node labels
cols<-setNames(palette()[1:length(unique(character1))],sort(unique(character1)))

#modify colors to custom choices
cols[1]<- "#17B21F"
cols[2]<-"#BC9568"
cols[3]<-"#AEB4B3"

#temporary plotting SIMMMAP  reconstruction 
plot(TimeTree.simmap.character1[[1]], cols, type='fan', fsize=0.5, ftype='i', lwd=1.5)

#temporary plotting SIMMMAP  reconstruction 
plot(TimeTree.simmap.ARD.character1[[1]], cols, type='fan', fsize=0.5, ftype='i', lwd=1.5)

#setting up concentric circles based on geologic timescale
obj<-geo.legend() ## this is just to get the colors

#defining colors for the relevant periods
obj$colors[1]<-"#ffffff"
obj$colors[2]<-"#F0F0F0"
obj$colors[3]<-"#ffffff"
obj$colors[4]<-"#F0F0F0"
obj$colors[5]<-"#ffffff"
obj$colors[6]<-"#F0F0F0"
obj$colors[7]<-"#ffffff"
obj$colors[8]<-"#F0F0F0"
obj$colors[9]<-"#ffffff"
obj$colors[10]<-"#F0F0F0"

r<-max(obj$leg[,1])-obj$leg[,2] #calculates appropriate radius lengths




#generate PDF summary for the first custom model
pdfgenerator(filename="Upham_et_al_2019.SIMMAP-2_rate_V1.pdf", basetree=(TimeTree), simsum=pd.1, title="Two rate model", modelfit=fitCustom.character1)
#generate PDF summary for the second custom model
pdfgenerator(filename="Upham_et_al_2019.SIMMAP-4_rate_V1.pdf", basetree=(TimeTree), simsum=pd.int.1, title="Four rate model", modelfit=fitCustom.int.character1)
#generate PDF summary for the ARD.model
pdfgenerator(filename="Upham_et_al_2019.SIMMAP-6_rate_V1.pdf", basetree=(TimeTree), simsum=pd.ARD.1, title="Six rate model", modelfit=fitARD.character1)


#generate PDF summary for the first custom model + posterior
pdfgenerator(filename="Upham_et_al_2019.SIMMAP-2_rate_pos_V1.pdf", basetree=(TimeTree), simsum=pd.1.pos, title="Two rate model, posterior trees", modelfit=fitCustom.character1, posmodelfit=Qs.pos, posterior=T)
#generate PDF summary for the second custom model + posterior
pdfgenerator(filename="Upham_et_al_2019.SIMMAP-4_rate_pos_V1.pdf", basetree=(TimeTree), simsum=pd.int.1.pos, title="Four rate model, posterior trees", modelfit = fitCustom.int.character1, posmodelfit=Qs.pos.int, posterior=T)
#generate PDF summary for the ARD.model + posterior
pdfgenerator(filename="Upham_et_al_2019.SIMMAP-6_rate_pos_V1.pdf", basetree=(TimeTree), simsum=pd.ARD.1.pos, title="Six rate model, posterior trees", modelfit = fitARD.character1, posmodelfit=Qs.pos.ARD, posterior=T)


#writing out the posterior probabilities
write.table(pd.1$ace, file='Upham_et_al_2019.SIMMAP-2_rate_nodeprobability.SIMMAP.V1.txt', quote=F, sep="\t", row.names=T, col.names=T)
write.table(pd.int.1$ace, file='Upham_et_al_2019.SIMMAP-4_rate_nodeprobability.SIMMAP.V1.txt', quote=F, sep="\t", row.names=T, col.names=T)
write.table(pd.ARD.1$ace, file='Upham_et_al_2019.SIMMAP-6_rate_nodeprobability.SIMMAP.V1.txt', quote=F, sep="\t", row.names=T, col.names=T)
write.table(pd.1.pos$ace, file='Upham_et_al_2019.SIMMAP-2_rate_pos_nodeprobability.SIMMAP.V1.txt', quote=F, sep="\t", row.names=T, col.names=T)
write.table(pd.int.1.pos$ace, file='Upham_et_al_2019.SIMMAP-4_rate_pos_nodeprobability.SIMMAP.V1.txt', quote=F, sep="\t", row.names=T, col.names=T)
write.table(pd.ARD.1.pos$ace, file='Upham_et_al_2019.SIMMAP-6_rate_pos_nodeprobability.SIMMAP.V1.txt', quote=F, sep="\t", row.names=T, col.names=T)

#this chunk no longer used
# #comparing marginal ancestral states to posterior probabilities
# plot(fitCustom.character1$lik.anc, pd.1$ace, xlab='marginal ancestral states', ylab='posterior probabilites from stochastic mapping')
# lines(c(0,1), c(0,1), lty='dashed', col='red', lwd=2)
# 
# #regression summaries indicating marginal ancestral states and posterior probabilities are virtually identical
# summary(lm(fitCustom.character1$lik.anc~pd.1$ace))
# summary(lm(fitCustom.character1$lik.anc[,1]~pd.1$ace[,1]))
# summary(lm(fitCustom.character1$lik.anc[,2]~pd.1$ace[,2]))
# summary(lm(fitCustom.character1$lik.anc[,3]~pd.1$ace[,3]))
# #R^2 > 0.99
# 


#plotting posterior densities from SIMMAP
dd<-density(TimeTree.simmap.character1)
pdf(file="Upham_et_al_2019.2_rate_SIMMAP.pd-V1.pdf", width=10, height=10)
plot(dd)
dev.off()

#plotting posterior densities from SIMMAP
dd<-density(TimeTree.simmap.character1.int)
pdf(file="Upham_et_al_2019.4_rate_SIMMAP.pd-V1.pdf", width=10, height=10)
plot(dd)
dev.off()

#plotting posterior densities from SIMMAP
dd<-density(TimeTree.simmap.ARD.character1)
pdf(file="Upham_et_al_2019.6_rate_SIMMAP.pd-V1.pdf", width=10, height=10)
plot(dd)
dev.off()


#plotting posterior densities from SIMMAP
system.time(dd<-density(simmap.character1.posterior))
pdf(file="Upham_et_al_2019.2_rate_pos_SIMMAP.pd-V1.pdf", width=10, height=10)
plot(dd)
rm(dd)
dev.off()

#plotting posterior densities from SIMMAP
dd<-density(simmap.character1.posterior.int)
pdf(file="Upham_et_al_2019.4_rate_pos_SIMMAP.pd-V1.pdf", width=10, height=10)
plot(dd)
rm(dd)
dev.off()

#plotting posterior densities from SIMMAP
dd<-density(simmap.ARD.character1.posterior)
pdf(file="Upham_et_al_2019.6_rate_pos_SIMMAP.pd-V1.pdf", width=10, height=10)
plot(dd)
rm(dd)
dev.off()



##########################################
#start plotting parsimony reconstruction#
##########################################

pdf(file="Upham_et_al_2019.MPR-V1.pdf", width=10, height=10)
#change to your output directory

plotTree(TimeTree, cols, type='fan', fsize=0.5, ftype='i', lwd=2.5, offset=3)

for(i in 1:nrow(obj$leg)){
  color<-paste(strsplit(obj$colors[i],"")[[1]][1:7],collapse="")
  draw.circle(0,0,radius=r[i],col=color,border="transparent")
  draw.circle(0,0, radius=r[i:min(i)+1], col='transparent', border='black', lty=2, lwd=0.25)
}

draw.circle(0,0, radius=r[4], col='transparent', border='red', lty=1, lwd=3)


par(new=T)
plotTree(TimeTree, cols, type='fan', fsize=0.5, ftype='i', lwd=2.5, offset=3)
#addcladelabels()
par(lwd = 0.2)

#add node labels
nodelabels(pie=pars, piecol=cols, prompt=FALSE, x=0.9*par()$usr[1], y=-max(nodeHeights(TimeTree)), cex= 0.3)

#add tip labels
tiplabels(pie=tips, piecol=cols, prompt=FALSE, cex=0.2)

#add legend
add.simmap.legend(colors=cols,prompt=FALSE,x=0.9*par()$usr[1], y=-max(nodeHeights(TimeTree)),fsize=0.8)
addcladelabels()

#add parsimony statistics
legend('topleft', legend=paste(' Consistency Index = ', round(consistency,3), '\n', 'Retention Index = ', round(retention,3), '\n', 'steps = ', round(score,3)), bty='n', cex=0.75)

dev.off() #end plotting


write.table(pars, file='Upham_et_al_2019.parsimony.V1.txt', quote=F, sep="\t", row.names=T, col.names=T)




#######





#calculate rates through time for consensus tree
#execution commants commented out

#counts.all<-getrates(tree=TimeTree, bins=50, maps=TimeTree.simmap.character1.int, type="all")
#saveRDS(counts.all, file="counts.all.RDS")
counts.all<-readRDS("upham_out/counts.all.RDS")

#counts.NS<-getrates(tree=TimeTree, bins=50, maps=TimeTree.simmap.character1.int, type="Non-arboreal->Semi-arboreal")
#saveRDS(counts.NS, file="counts.NS.RDS")
counts.NS<-readRDS("upham_out/counts.NS.RDS")

#counts.SN<-getrates(tree=TimeTree, bins=50, maps=TimeTree.simmap.character1.int, type="Semi-arboreal->Non-arboreal")
#saveRDS(counts.SN, file="counts.SN.RDS")
counts.SN<-readRDS("upham_out/counts.SN.RDS")

#counts.SA<-getrates(tree=TimeTree, bins=50, maps=TimeTree.simmap.character1.int, type="Semi-arboreal->Arboreal")
#saveRDS(counts.SA, file="counts.SA.RDS")
counts.SA<-readRDS("upham_out/counts.SA.RDS")

#counts.AS<-getrates(tree=TimeTree, bins=50, maps=TimeTree.simmap.character1.int, type="Arboreal->Semi-arboreal")
#saveRDS(counts.AS, file="counts.AS.RDS")
counts.AS<-readRDS("upham_out/counts.AS.RDS")


#executes with no errors


#calculate rates through time for posterior trees
#first, reformat the posterior data object using the splitter function
split.simmap.character1.posterior.int<-splitter(maps=simmap.character1.posterior.int, length=500)

#now, calculate rates through time for the posterior trees
#all.rates.simmap.character1.posterior.int<-getratesloop(posteriortrees=posterior_trees.sub, bins=50, splitmaps=split.simmap.character1.posterior.int, type="all")
#saveRDS(all.rates.simmap.character1.posterior.int, file="all.rates.simmap.character1.posterior.int.RDS")
all.rates.simmap.character1.posterior.int<-readRDS("all.rates.simmap.character1.posterior.int.RDS")

#NS.rates.simmap.character1.posterior.int<-getratesloop(posteriortrees=posterior_trees.sub, bins=50, splitmaps=split.simmap.character1.posterior.int, type="Non-arboreal->Semi-arboreal")
#saveRDS(NS.rates.simmap.character1.posterior.int, file="NS.rates.simmap.character1.posterior.int.RDS")
NS.rates.simmap.character1.posterior.int<-readRDS("NS.rates.simmap.character1.posterior.int.RDS")

#SN.rates.simmap.character1.posterior.int<-getratesloop(posteriortrees=posterior_trees.sub, bins=50, splitmaps=split.simmap.character1.posterior.int, type="Semi-arboreal->Non-arboreal")
#saveRDS(SN.rates.simmap.character1.posterior.int, file="SN.rates.simmap.character1.posterior.int.RDS")
SN.rates.simmap.character1.posterior.int<-readRDS("SN.rates.simmap.character1.posterior.int.RDS")

#SA.rates.simmap.character1.posterior.int<-getratesloop(posteriortrees=posterior_trees.sub, bins=50, splitmaps=split.simmap.character1.posterior.int, type="Semi-arboreal->Arboreal")
#saveRDS(SA.rates.simmap.character1.posterior.int, file="SA.rates.simmap.character1.posterior.int.RDS")
SA.rates.simmap.character1.posterior.int<-readRDS("SA.rates.simmap.character1.posterior.int.RDS")

#AS.rates.simmap.character1.posterior.int<-getratesloop(posteriortrees=posterior_trees.sub, bins=50, splitmaps=split.simmap.character1.posterior.int, type="Arboreal->Semi-arboreal")
#saveRDS(AS.rates.simmap.character1.posterior.int, file="AS.rates.simmap.character1.posterior.int.RDS")
AS.rates.simmap.character1.posterior.int<-readRDS("AS.rates.simmap.character1.posterior.int.RDS")
 

#which nodes are identified to be kpg associated in upham and meredith
kpgNodes_upham<-c(166, 290, 291, 294, 300, 303, 302, 295, 170, 171, 180, 219, 238, 239, 242, 243, 255, 256, 260, 278)
kpgNodes_meridith<-c(166, 305, 311, 312, 290, 303, 295, 180, 219, 232, 242, 243, 256, 260, 278, 258)
kpgNodes_both<-intersect(kpgNodes_upham, kpgNodes_meridith)


#testing rate plots
rateplot(rates=counts.all, tree=TimeTree, ylim=c(0,0.035), spline=T, width=3, alpha=0.3)
rateplot(rates=counts.NS, tree=TimeTree, ylim=c(0,0.035), spline=T, width=3, alpha=0.3)
rateplot(rates=counts.SN, tree=TimeTree, ylim=c(0,0.035), spline=T, width=3, alpha=0.3)
rateplot(rates=counts.SA, tree=TimeTree, ylim=c(0,0.035), spline=T, width=3, alpha=0.3)
rateplot(rates=counts.AS, tree=TimeTree, ylim=c(0,0.035), spline=T, width=3, alpha=0.3)



#test plotting with the posterior lines

pdf(file="all.pdf", width=11, height=8.5)
#plot the line for all transitions
rateplot(rates=counts.all, tree=TimeTree, ylim=c(0,0.035), spline=T, width=3, alpha=0.3)

#plot the lines for posterior transitions
segplotloop(rates = all.rates.simmap.character1.posterior.int, ylim = c(0, 0.035), color = "#0000001A", spline = F, width = 0.001, xlim = c(max(nodeHeights(TimeTree)),  0), ax=F)
dev.off()


pdf(file="NS.pdf", width=11, height=8.5)
#plot the line for all transitions
rateplot(rates=counts.NS, tree=TimeTree, ylim=c(0,0.035), spline=T, width=3, alpha=0.3)

#plot the lines for posterior transitions
segplotloop(rates = NS.rates.simmap.character1.posterior.int, ylim = c(0,0.035), color = "#0000FF1A", spline = F, width = 0.001, xlim = c(max(nodeHeights(TimeTree)),  0), ax=F)
dev.off()


pdf(file="SN.pdf", width=11, height=8.5)
#plot the line for all transitions
rateplot(rates=counts.SN, tree=TimeTree, ylim=c(0,0.035), spline=T, width=3, alpha=0.3)

#plot the lines for posterior transitions
segplotloop(rates = SN.rates.simmap.character1.posterior.int, ylim = c(0, 0.035), color = "#8000801A", spline = F, width = 0.001, xlim = c(max(nodeHeights(TimeTree)),  0), ax=F)
dev.off()


pdf(file="SA.pdf", width=11, height=8.5)
#plot the line for all transitions
rateplot(rates=counts.SA, tree=TimeTree, ylim=c(0,0.035), spline=T, width=3, alpha=0.3)

#plot the lines for posterior transitions
segplotloop(rates = SA.rates.simmap.character1.posterior.int, ylim = c(0, 0.035), color = "#FFA5001A", spline = F, width = 0.001, xlim = c(max(nodeHeights(TimeTree)),  0), ax=F)
dev.off()


pdf(file="AS.pdf", width=11, height=8.5)
#plot the line for all transitions
rateplot(rates=counts.AS, tree=TimeTree, ylim=c(0,0.055), spline=T, width=3, alpha=0.3)

#plot the lines for posterior transitions
segplotloop(rates = AS.rates.simmap.character1.posterior.int, ylim = c(0, 0.035), color = "#00FF001A", spline = F, width = 0.001, xlim = c(max(nodeHeights(TimeTree)),  0), ax=F)
dev.off()




#get plotting data for the rectangle
lastPP <- get("last_plot.phylo", envir = .PlotPhyloEnv)
minx<-min(lastPP$xx[union(kpgNodes_meridith, kpgNodes_upham)])
maxx<-max(lastPP$xx[union(kpgNodes_meridith, kpgNodes_upham)])



pdf(file="consensusrates.pdf", width=11, height=8.5)

rateplot(rates=counts.all, tree=TimeTree, ylim=c(0,0.035), spline=T, width=3, alpha=0.3)
addcladelabels.flat(ln.offset=123.5, wing=0.5, size=0.65)

par(lwd=0.75)
nodelabels(text=rep("", length(kpgNodes_upham)), node=kpgNodes_upham, frame="circle", col="black", bg="black", cex=0.4)
nodelabels(text=rep("", length(kpgNodes_meridith)), node=kpgNodes_meridith, frame="circle", col="black", bg="white", cex=0.4)
nodelabels(text=rep("", length(kpgNodes_both)), node=kpgNodes_both, frame="circle", col="black", bg="grey", cex=0.4)
rect(xleft=minx, ybottom=-100, xright=maxx, ytop=200, col="#D2B48C1A", border=NA)

legend("topleft", legend=c("All transitions", "Non-Arboreal to Semi-Arboreal", "Semi-Arboreal to Non-Arboreal", "Semi-Arboreal to Arboreal","Arboreal to Semi-Arboreal"), col = c("black", "blue","purple","orange","green"), bty = 'n', lty=1, lwd=8)

par(new=T)
segplot(rates=counts.NS, ylim=c(0, 0.035), color="blue", spline=T, width=3, ax=F) #non to semi
par(new=T)
segplot(rates=counts.SN, ylim=c(0, 0.035), color="purple", spline=T, width=3, ax=F) #semi to non
par(new=T)
segplot(rates=counts.SA, ylim=c(0, 0.035), color="orange", spline=T, width=3, ax=F) #semi to ar
par(new=T)
segplot(rates=counts.AS, ylim=c(0, 0.035), color="green", spline=T, width=3, ax=F) #ar to semi

abline(v=66, lty=2, lwd=0.75)

dev.off()







pdf(file="consensus_posterior.pdf", width=11, height=8.5)
rateplot(rates=counts.all, tree=TimeTree, ylim=c(0,0.035), spline=T, width=3, alpha=0.3)
#addcladelabels.flat(ln.offset=161.5, wing=0.5, size=0.5)

#par(lwd=0.001)
#nodelabels(text=rep("", length(kpgNodes_upham)), node=kpgNodes_upham, frame="circle", col="black", bg="black", cex=0.5)
#nodelabels(text=rep("", length(kpgNodes_meridith)), node=kpgNodes_meridith, frame="circle", col="red", bg="red", cex=0.35)
rect(xleft=minx, ybottom=-100, xright=maxx, ytop=200, col="#D2B48C1A", border=NA)


#segplotloop(rates = all.rates.simmap.character1.posterior.int, ylim = c(0, 0.035), color = "#0000001A", spline = F, width = 0.001, xlim = c(max(nodeHeights(TimeTree)),  0), ax=F)
segplotloop(rates = AS.rates.simmap.character1.posterior.int, ylim = c(0, 0.035), color = "#00FF001F", spline = T, width = 0.001, xlim = c(max(nodeHeights(TimeTree)),  0), ax=F)
segplotloop(rates = SN.rates.simmap.character1.posterior.int, ylim = c(0, 0.035), color = "#8000801F", spline = T, width = 0.001, xlim = c(max(nodeHeights(TimeTree)),  0), ax=F)
segplotloop(rates = SA.rates.simmap.character1.posterior.int, ylim = c(0, 0.035), color = "#FFA50014", spline = T, width = 0.001, xlim = c(max(nodeHeights(TimeTree)),  0), ax=F)
segplotloop(rates = NS.rates.simmap.character1.posterior.int, ylim = c(0,0.035), color = "#0000FF0D", spline = T, width = 0.001, xlim = c(max(nodeHeights(TimeTree)),  0), ax=F)

par(new=T)
segplot(rates=counts.AS, ylim=c(0, 0.035), color="green", spline=T, width=3, ax=F, y.lab="", x.lab="")

par(new=T)
segplot(rates=counts.SN, ylim=c(0, 0.035), color="purple", spline=T, width=3, ax=F, y.lab="", x.lab="") #semi to non

par(new=T)
segplot(rates=counts.SA, ylim=c(0, 0.035), color="orange", spline=T, width=3, ax=F, y.lab="", x.lab="") #semi to ar

par(new=T)
segplot(rates=counts.NS, ylim=c(0, 0.035), color="blue", spline=T, width=3, ax=F, y.lab="", x.lab="") #non to semi

#legend("topleft", legend=c("All transitions", "Non-Arboreal to Semi-Arboreal", "Semi-Arboreal to Non-Arboreal", "Semi-Arboreal to Arboreal","Arboreal to Semi-Arboreal"), col = c("black", "blue","purple","orange","green"), bty = 'n', lty=1, lwd=8)

abline(v=66, lty=2, lwd=0.75)

dev.off()









######
#debugging
fitMk(tree=posterior_trees.sub[[26]], model=transitions.int, x=character1, pi="fitzjohn")


NS.rates.simmap.character1.posterior.int<-getratesloop(posteriortrees=posterior_trees.sub[26], bins=50, splitmaps=split.simmap.character1.posterior.int[26], type="Non-arboreal->Semi-arboreal")


describe.simmap(split.simmap.character1.posterior.int[[26]])
split.simmap.character1.posterior.int[[26]][[2]]$Q



plot(split.simmap.character1.posterior.int[[26]][[2]])



plotTree(TimeTree,add=F,ftype="i",lwd=1,color=make.transparent("blue", 0.5),
         mar=par()$mar,direction="rightwards", fsize=0.25)

nodelabels(cex=0.25, frame="circle")


