#R code for general nest category reconstructions
#data collected by Johnathan Hughes
#code by Jacob Berv, jsb439@cornell.edu
#heavily inspired by Liam Revell's phytools examples


#load requried packages
require(phytools)
require(plotrix)
require(ape)

#set working directory
setwd("~/jsb439@cornell.edu/Dan-Shared/the_trees_are_dead_part2/analyses")

#load the tree and clean up
TimeTree<-read.tree("AAtimetreeAUTOhardMEAN.tree")

#import the data table
data<-read.table("mammal_data_v2.txt", header=T)
rownames(data)<-data$Taxon

#reorder the data frame to match the tip order
data<-data[TimeTree$tip.label,]

#remove the outgroups
remove<-c('Danio', 'Xenopus', 'Anolis', 'Gallus', 'Taeniopygia')

#prune outgroups from datset
data <-data[!data$Taxon %in% remove,]

#prune outgroups from tree
TimeTree<-drop.tip(TimeTree, tip=remove)


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
transitions <- matrix(c(0, 1, 2, 1, 0, 2, 2, 2, 0), nrow=3)

#Rate index matrix:
#             Arboreal Non-arboreal Semi-arboreal
#Arboreal             .            1             2
#Non-arboreal         1            .             2
#Semi-arboreal        2            2             .

#maybe try setting non-arboreal to arboreal to it's own rate class


#running everything for character1 first

#fit custom transition model with ace, for maximum likelihood estimate, with different initial starting value
fitCustom.character1 <- ace(character1, TimeTree, model=transitions, type='discrete', ip=c(0.001, 0.001))
AIC(fitCustom.character1) #249.4975


fitARD.character1 <- ace(character1, TimeTree, model='ARD', type='discrete', ip=c(0.001, 0.001))
AIC(fitARD.character1) #242.3074, but runs with errors - revisit this

#trying fitMK instead of ace -- results are the same, so ignoring
test<-fitMk(TimeTree, character1, model='ARD', type='discrete')
AIC(test)

#setting up node labels
cols<-setNames(palette()[1:length(unique(character1))],sort(unique(character1)))

#modify colors to custom choices
cols[1]<- "#17B21F"
cols[2]<-"#BC9568"
cols[3]<-"#AEB4B3"

#stochastic character1 map with custom model
TimeTree.simmap.character1<-make.simmap(TimeTree, character1, model=transitions, nsim=1000, message=T, type='discrete')

#stochastic character1 map with ARD model
TimeTree.simmap.ARD.character1<-make.simmap(TimeTree, character1, model='ARD', nsim=1000, message=T, type='discrete')


#posterior summary of stochastic mapping
pd.1<-summary(TimeTree.simmap.character1, plot=F)
pd.ARD.1<- summary(TimeTree.simmap.ARD.character1, plot=F)


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

#generate PDF plot

#function for adding clade arc labels for mammalian groups
addcladelabels<-function(){
  arc.cladelabels(tree=TimeTree,'Marsupialia', node=168, ln.offset=1.20, lab.offset=1.25, cex=1, mark.node=F, lwd=3, lend=1, col='grey48')
  arc.cladelabels(tree=TimeTree,'Monotremata', node=166, ln.offset=1.20, lab.offset=1.25, cex=1, mark.node=F, lwd=3, lend=1, col='grey48')
  arc.cladelabels(tree=TimeTree,'Afrotheria', node=317, ln.offset=1.20, lab.offset=1.25, cex=1, mark.node=F, lwd=3, lend=1, col='grey48')
  arc.cladelabels(tree=TimeTree,'Xenarthra', node=313, ln.offset=1.20, lab.offset=1.25, cex=1, mark.node=F, lwd=3, lend=1, col='grey48')
  arc.cladelabels(tree=TimeTree,'Primatomorpha', node=298, ln.offset=1.20, lab.offset=1.25, cex=1, mark.node=F, lwd=3, lend=1, col='grey48')
  arc.cladelabels(tree=TimeTree,'Glires', node=264, ln.offset=1.20, lab.offset=1.25, cex=1, mark.node=F, lwd=3, lend=1, col='grey48')
  arc.cladelabels(tree=TimeTree,'Laurasiatheria', node=193, ln.offset=1.20, lab.offset=1.25, cex=1, mark.node=F, lwd=3, lend=1, col='grey48')
  arc.cladelabels(tree=TimeTree,'Scandentia', node=263, ln.offset=1.20, lab.offset=1.25, cex=1, mark.node=F, lwd=3, lend=1, col='grey48')
}

#start plotting
pdf(file="Meredith_etal2011.SIMMAP-V1.pdf", width=10, height=10)
#change to your output directory

plot(TimeTree.simmap.character1[[1]], cols, type='fan', fsize=0.5, ftype='i', lwd=2.5)

for(i in 1:nrow(obj$leg)){
  color<-paste(strsplit(obj$colors[i],"")[[1]][1:7],collapse="")
  draw.circle(0,0,radius=r[i],col=color,border="transparent")
  draw.circle(0,0, radius=r[i:min(i)+1], col='transparent', border='black', lty=2, lwd=0.25)
}

draw.circle(0,0, radius=r[4], col='transparent', border='red', lty=1, lwd=3)

par(new=T)
plot(TimeTree.simmap.character1[[1]], cols, type='fan', fsize=0.5, ftype='i', lwd=2.5)
addcladelabels()

par(lwd = 0.2)

#add node labels
nodelabels(pie=pd.1$ace, piecol=cols, prompt=FALSE, x=0.9*par()$usr[1], y=-max(nodeHeights(TimeTree)), cex= 0.3)

#add tip labels # can only be executed after pars recon for this sectino
tiplabels(pie=tips, piecol=cols, prompt=FALSE, cex=0.1)

#add legend
add.simmap.legend(colors=cols,prompt=FALSE,x=0.9*par()$usr[1], y=-max(nodeHeights(TimeTree)),fsize=0.8)


dev.off() #end plotting

write.table(pd.1$ace, file='Meredith_etal2011.posteriorprobability.SIMMAP.V1.txt', quote=F, sep="\t", row.names=T, col.names=T)

#comparing marginal ancestral states to posterior probabilities
plot(fitCustom.character1$lik.anc, pd.1$ace, xlab='marginal ancestral states', ylab='posterior probabilites from stochastic mapping')
lines(c(0,1), c(0,1), lty='dashed', col='red', lwd=2)

#regression summaries indicating marginal ancestral states and posterior probabilities are virtually identical
summary(lm(fitCustom.character1$lik.anc~pd.1$ace))
summary(lm(fitCustom.character1$lik.anc[,1]~pd.1$ace[,1]))
summary(lm(fitCustom.character1$lik.anc[,2]~pd.1$ace[,2]))
summary(lm(fitCustom.character1$lik.anc[,3]~pd.1$ace[,3]))
#R^2 > 0.99

#plotting posterior densities from SIMMAP
dd<-density(TimeTree.simmap.character1)
pdf(file="Meredith_etal2011.SIMMAP.pd-V1.pdf", width=10, height=10)
plot(dd)
dev.off()


#start plotting
pdf(file="Meredith_etal2011.SIMMAP.ARD-V1.pdf", width=10, height=10)
#change to your output directory

plot(TimeTree.simmap.ARD.character1[[1]], cols, type='fan', fsize=0.5, ftype='i', lwd=2.5)

for(i in 1:nrow(obj$leg)){
  color<-paste(strsplit(obj$colors[i],"")[[1]][1:7],collapse="")
  draw.circle(0,0,radius=r[i],col=color,border="transparent")
  draw.circle(0,0, radius=r[i:min(i)+1], col='transparent', border='black', lty=2, lwd=0.25)
}

draw.circle(0,0, radius=r[4], col='transparent', border='red', lty=1, lwd=3)

par(new=T)
plot(TimeTree.simmap.ARD.character1[[1]], cols, type='fan', fsize=0.5, ftype='i', lwd=2.5)
addcladelabels()

par(lwd = 0.2)

#add node labels
nodelabels(pie=pd.ARD.1$ace, piecol=cols, prompt=FALSE, x=0.9*par()$usr[1], y=-max(nodeHeights(TimeTree)), cex= 0.3)

#add tip labels # can only be executed after pars recon for this sectino
tiplabels(pie=tips, piecol=cols, prompt=FALSE, cex=0.1)


#add legend
add.simmap.legend(colors=cols,prompt=FALSE,x=0.9*par()$usr[1], y=-max(nodeHeights(TimeTree)),fsize=0.8)

dev.off() #end plotting

write.table(pd.ARD.1$ace, file='Meredith_etal2011.posteriorprobability.SIMMAP.ARD.V1.txt', quote=F, sep="\t", row.names=T, col.names=T)


#comparing marginal ancestral states to posterior probabilities
plot(fitARD.character1$lik.anc, pd.ARD.1$ace, xlab='marginal ancestral states', ylab='posterior probabilites from stochastic mapping')
lines(c(0,1), c(0,1), lty='dashed', col='red', lwd=2)

#regression summaries indicating marginal ancestral states and posterior probabilities are virtually identical across characters
summary(lm(fitARD.character1$lik.anc~pd.ARD.1$ace))
summary(lm(fitARD.character1$lik.anc[,1]~pd.ARD.1$ace[,1]))
summary(lm(fitARD.character1$lik.anc[,2]~pd.ARD.1$ace[,2]))
summary(lm(fitARD.character1$lik.anc[,3]~pd.ARD.1$ace[,3]))
#R^2 > 0.96 for all

#plotting posterior densities from SIMMAP
dd<-density(TimeTree.simmap.ARD.character1)
pdf(file="Meredith_etal2011.SIMMAP.pd.ARD-V1.pdf", width=10, height=10)
plot(dd)
dev.off()



#start plotting parsimony reconstruction
pdf(file="Meredith_etal2011.MPR-V1.pdf", width=10, height=10)
#change to your output directory

plotTree(TimeTree, cols, type='fan', fsize=0.5, ftype='i', lwd=2.5)

for(i in 1:nrow(obj$leg)){
  color<-paste(strsplit(obj$colors[i],"")[[1]][1:7],collapse="")
  draw.circle(0,0,radius=r[i],col=color,border="transparent")
  draw.circle(0,0, radius=r[i:min(i)+1], col='transparent', border='black', lty=2, lwd=0.25)
}

draw.circle(0,0, radius=r[4], col='transparent', border='red', lty=1, lwd=3)


par(new=T)
plotTree(TimeTree, cols, type='fan', fsize=0.5, ftype='i', lwd=2.5)
addcladelabels()
par(lwd = 0.2)

#add node labels
nodelabels(pie=pars, piecol=cols, prompt=FALSE, x=0.9*par()$usr[1], y=-max(nodeHeights(TimeTree)), cex= 0.3)

#add tip labels
tiplabels(pie=tips, piecol=cols, prompt=FALSE, cex=0.1)

#add legend
add.simmap.legend(colors=cols,prompt=FALSE,x=0.9*par()$usr[1], y=-max(nodeHeights(TimeTree)),fsize=0.8)

#add parsimony statistics
legend('topleft', legend=paste(' Consistency Index = ', round(consistency,3), '\n', 'Retention Index = ', round(retention,3), '\n', 'steps = ', round(score,3)), bty='n', cex=0.75)

dev.off() #end plotting


write.table(pars, file='Meredith_etal2011.parsimony.V1.txt', quote=F, sep="\t", row.names=T, col.names=T)



#redo everything for character2

#begin parsimony reconstruction using phangorn
require(phangorn)

tmp<-c(as.character(character2))
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



#fit custom transition model with ace, for maximum likelihood estimate, with different initial starting value
fitCustom.character2 <- ace(character2, TimeTree, model=transitions, type='discrete', ip=c(0.001, 0.001))
AIC(fitCustom.character2) #240.3987

fitARD.character2 <- ace(character2, TimeTree, model='ARD', type='discrete', ip=c(0.001, 0.001))
AIC(fitARD.character2) #233.7194, but runs with errors - revisit this

test<-fitMk(TimeTree, character2, model='ARD', type='discrete')
AIC(test)

#setting up node labels
cols<-setNames(palette()[1:length(unique(character2))],sort(unique(character2)))

#modify colors to custom choices
cols[1]<- "#17B21F"
cols[2]<-"#BC9568"
cols[3]<-"#AEB4B3"

#stochastic character2 map with custom model
TimeTree.simmap.character2<-make.simmap(TimeTree, character2, model=transitions, nsim=1000, message=T, type='discrete')

#stochastic character2 map with ARD model
TimeTree.simmap.ARD.character2<-make.simmap(TimeTree, character2, model='ARD', nsim=1000, message=T, type='discrete')

#posterior summary of stochastic mapping
pd.2<-summary(TimeTree.simmap.character2, plot=F)
pd.ARD.2<- summary(TimeTree.simmap.ARD.character2, plot=F)


#temporary plotting SIMMMAP  reconstruction 
plot(TimeTree.simmap.character2[[1]], cols, type='fan', fsize=0.5, ftype='i', lwd=1.5)

#temporary plotting SIMMMAP  reconstruction 
plot(TimeTree.simmap.ARD.character2[[1]], cols, type='fan', fsize=0.5, ftype='i', lwd=1.5)


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

#generate PDF plot

#start plotting
pdf(file="Meredith_etal2011.SIMMAP-V2.pdf", width=10, height=10)
#change to your output directory

plot(TimeTree.simmap.character2[[1]], cols, type='fan', fsize=0.5, ftype='i', lwd=2.5)

for(i in 1:nrow(obj$leg)){
  color<-paste(strsplit(obj$colors[i],"")[[1]][1:7],collapse="")
  draw.circle(0,0,radius=r[i],col=color,border="transparent")
  draw.circle(0,0, radius=r[i:min(i)+1], col='transparent', border='black', lty=2, lwd=0.25)
}

draw.circle(0,0, radius=r[4], col='transparent', border='red', lty=1, lwd=3)

par(new=T)
plot(TimeTree.simmap.character2[[1]], cols, type='fan', fsize=0.5, ftype='i', lwd=2.5)
addcladelabels()

par(lwd = 0.2)

#add node labels
nodelabels(pie=pd.2$ace, piecol=cols, prompt=FALSE, x=0.9*par()$usr[1], y=-max(nodeHeights(TimeTree)), cex= 0.3)

#add tip labels, only run after parsimony
tiplabels(pie=tips, piecol=cols, prompt=FALSE, cex=0.1)

#add legend
add.simmap.legend(colors=cols,prompt=FALSE,x=0.9*par()$usr[1], y=-max(nodeHeights(TimeTree)),fsize=0.8)


dev.off() #end plotting

write.table(pd.2$ace, file='Meredith_etal2011.posteriorprobability.SIMMAP.V2.txt', quote=F, sep="\t", row.names=T, col.names=T)


#comparing marginal ancestral states to posterior probabilities
plot(fitCustom.character2$lik.anc, pd.2$ace, xlab='marginal ancestral states', ylab='posterior probabilites from stochastic mapping')
lines(c(0,1), c(0,1), lty='dashed', col='red', lwd=2)

#regression summaries indicating marginal ancestral states and posterior probabilities are virtually identical
summary(lm(fitCustom.character2$lik.anc~pd.2$ace))
summary(lm(fitCustom.character2$lik.anc[,1]~pd.2$ace[,1]))
summary(lm(fitCustom.character2$lik.anc[,2]~pd.2$ace[,2]))
summary(lm(fitCustom.character2$lik.anc[,3]~pd.2$ace[,3]))
#R^2 > 0.99


#plotting posterior densities from SIMMAP
dd<-density(TimeTree.simmap.character2)
pdf(file="Meredith_etal2011.SIMMAP.pd-V2.pdf", width=10, height=10)
plot(dd)
dev.off()


#start plotting
pdf(file="Meredith_etal2011.SIMMAP.ARD-V2.pdf", width=10, height=10)
#change to your output directory

plot(TimeTree.simmap.ARD.character2[[1]], cols, type='fan', fsize=0.5, ftype='i', lwd=2.5)

for(i in 1:nrow(obj$leg)){
  color<-paste(strsplit(obj$colors[i],"")[[1]][1:7],collapse="")
  draw.circle(0,0,radius=r[i],col=color,border="transparent")
  draw.circle(0,0, radius=r[i:min(i)+1], col='transparent', border='black', lty=2, lwd=0.25)
}

draw.circle(0,0, radius=r[4], col='transparent', border='red', lty=1, lwd=3)


par(new=T)
plot(TimeTree.simmap.ARD.character2[[1]], cols, type='fan', fsize=0.5, ftype='i', lwd=2.5)
addcladelabels()

par(lwd = 0.2)

#add node labels
nodelabels(pie=pd.ARD.2$ace, piecol=cols, prompt=FALSE, x=0.9*par()$usr[1], y=-max(nodeHeights(TimeTree)), cex= 0.3)

#add legend
add.simmap.legend(colors=cols,prompt=FALSE,x=0.9*par()$usr[1], y=-max(nodeHeights(TimeTree)),fsize=0.8)

dev.off() #end plotting

write.table(pd.ARD.2$ace, file='Meredith_etal2011.posteriorprobability.SIMMAP.ARD.V2.txt', quote=F, sep="\t", row.names=T, col.names=T)


#comparing marginal ancestral states to posterior probabilities
plot(fitARD.character2$lik.anc, pd.ARD.2$ace, xlab='marginal ancestral states', ylab='posterior probabilites from stochastic mapping')
lines(c(0,1), c(0,1), lty='dashed', col='red', lwd=2)

#regression summaries indicating marginal ancestral states and posterior probabilities are virtually identical
summary(lm(fitARD.character2$lik.anc~pd.ARD.2$ace))
summary(lm(fitARD.character2$lik.anc[,1]~pd.ARD.2$ace[,1]))
summary(lm(fitARD.character2$lik.anc[,2]~pd.ARD.2$ace[,2]))
summary(lm(fitARD.character2$lik.anc[,3]~pd.ARD.2$ace[,3]))
#R^2 > 0.98

#plotting posterior densities from SIMMAP
dd<-density(TimeTree.simmap.ARD.character2)
pdf(file="Meredith_etal2011.SIMMAP.pd.ARD-V2.pdf", width=10, height=10)
plot(dd)
dev.off()


#start plotting parsimony reconstruction
pdf(file="Meredith_etal2011.MPR-V2.pdf", width=10, height=10)
#change to your output directory

plotTree(TimeTree, cols, type='fan', fsize=0.5, ftype='i', lwd=2.5)

for(i in 1:nrow(obj$leg)){
  color<-paste(strsplit(obj$colors[i],"")[[1]][1:7],collapse="")
  draw.circle(0,0,radius=r[i],col=color,border="transparent")
  draw.circle(0,0, radius=r[i:min(i)+1], col='transparent', border='black', lty=2, lwd=0.25)
}

draw.circle(0,0, radius=r[4], col='transparent', border='red', lty=1, lwd=3)

par(new=T)
plotTree(TimeTree, cols, type='fan', fsize=0.5, ftype='i', lwd=2.5)
addcladelabels()
par(lwd = 0.2)

#add node labels
nodelabels(pie=pars, piecol=cols, prompt=FALSE, x=0.9*par()$usr[1], y=-max(nodeHeights(TimeTree)), cex= 0.3)

#add tip labels
tiplabels(pie=tips, piecol=cols, prompt=FALSE, cex=0.1)

#add legend
add.simmap.legend(colors=cols,prompt=FALSE,x=0.9*par()$usr[1], y=-max(nodeHeights(TimeTree)),fsize=0.8)

#add parsimony statistics
legend('topleft', legend=paste(' Consistency Index = ', round(consistency,3), '\n', 'Retention Index = ', round(retention,3), '\n', 'steps = ', round(score,3)), bty='n', cex=0.75)

dev.off() #end plotting

write.table(pars, file='Meredith_etal2011.parsimony.V2.txt', quote=F, sep="\t", row.names=T, col.names=T)

