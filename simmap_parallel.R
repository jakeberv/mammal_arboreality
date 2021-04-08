#some dependencies include
# require(phytools)
# require(plotrix)
# require(ape)
# require(stringr)
# require(ratematrix)
# require(parallel)
# library(doSNOW)
# library(doParallel)
# library(devtools)
# library(pbmcapply)


#fitQ fits a Q matrix under ML for a given tree, character, model, and root state assumption
#tree must be class phylo
#model can be ARD, ER, SYM, or matrix
#pi can be equal or fitzjohn (fitzjohn 2009)
fitQ<-function(tree, characters, model, pi){
  fit<-fitMk(tree=tree, x=characters, model=model, pi=pi)
  ## extracted the fitted transition matrix:
  fittedQ<-matrix(NA,length(fit$states),length(fit$states))
  fittedQ[]<-c(0,fit$rates)[fit$index.matrix+1]
  diag(fittedQ)<-0
  diag(fittedQ)<--rowSums(fittedQ)
  colnames(fittedQ)<-rownames(fittedQ)<-fit$states
  return(fittedQ)
}

# Qs<-function(tree, characters, model, pi){
# 
#   if (class(tree)=="multiPhylo") {
#     
#     
#   } else {
#     fit<-fitMk(tree=tree, x=characters, model=model, pi=pi)
#     ## extracted the fitted transition matrix:
#     fittedQ<-matrix(NA,length(fit$states),length(fit$states))
#     fittedQ[]<-c(0,fit$rates)[fit$index.matrix+1]
#     diag(fittedQ)<-0
#     diag(fittedQ)<--rowSums(fittedQ)
#     colnames(fittedQ)<-rownames(fittedQ)<-fit$states
#     
#   }
#   
#   return(fittedQ)
# }

#function to estimate the Q matrix in parallel across posterior trees
# Qs_loop<-function(i, trees, characters, mod){
#   if (class(trees)=="multiPhylo") {
#     Qs<-make.simmap(tree=trees[[i]], x=characters, model=mod, nsim=1)
#   } else {
#     Qs<-make.simmap(tree=trees, x=characters, model=mod, nsim=1)
#   }
#   #print(paste("parsed ", i, "trees out of ", length(targets), "or ", round(i/length(targets)*100, digits=3), " %"))
#   return(Qs)
# } #backup
# Qs_loop.mc<-function(trees, characters, mod, fun, cores){
# #set up cluster
# #define number of cores to use (max cores - 1)
# no_cores <- cores
# #how long is your data
# n <- length(trees)
# #spool up the cluster
# cl <- makeCluster(no_cores)
# registerDoSNOW(cl)
# 
# #for the length of data, run fun_parallel on each in parallel
# result<-foreach(i=1:n, .export=ls(globalenv()), .packages=c('ape', 'phytools', 'geiger'), .verbose=T) %dopar% fun(i, trees, characters, mod)
# 
# if(!("multiSimmap"%in%class(result))) class(result)<-c("multiSimmap","multiPhylo", class(result))
# 
# #kill the cluster
# stopCluster(cl)
# 
# #return data
# return(result)
# } #backup

#this funciton is a helper to estimate the Q matrix in parallel
Qs_loop<-function(i, trees, characters, mod, pi){
  if (class(trees)=="multiPhylo") {
    Qs<-fitQ(tree=trees[[i]], characters, model = mod, pi)
  } else {
    Qs<-fitQ(tree = trees, characters, model = mod, pi)
  }
  #print(paste("parsed ", i, "trees out of ", length(targets), "or ", round(i/length(targets)*100, digits=3), " %"))
  return(Qs)
}

#this functions generates a local cluster to loop the Q matrix estimation
Qs_loop.mc<-function(trees, characters, mod, pi, fun, cores){
  #set up cluster
  #define number of cores to use (max cores - 1)
  no_cores <- cores
  #how long is your data
  n <- length(trees)
  #spool up the cluster
  cl <- makeCluster(no_cores)
  print(paste("cluster started", "with", no_cores, "cores"))
  registerDoSNOW(cl)
  
  #for the length of data, run fun_parallel on each in parallel
  result<-foreach(i=1:n, .export=ls(globalenv()), .packages=c('ape', 'phytools', 'geiger'), .verbose=T) %dopar% fun(i, trees, characters, mod, pi)
  
  #if(!("multiSimmap"%in%class(result))) class(result)<-c("multiSimmap","multiPhylo", class(result))
  
  #kill the cluster
  stopCluster(cl)
  
  #return data
  return(result)
}


#function to estimate Q matrix for each posterior tree
#then simulate nsim stochastic maps on each, and combine into
#one multiSimmap object
# simmap.parallel<-function(trees, model, sims, characters, pi, cores){
#   #generate a distribution of Q matricies from the posterior distribution
#   #Qs <- make.simmap(trees, characters, model, nsim=1)
#   if (class(trees)=="multiPhylo"){
#     Qs <- Qs_loop.mc(trees=trees, characters = characters, mod=model, pi=pi, fun=Qs_loop, cores = cores)
#   } else {
#     Qs <- Qs_loop(i=1, trees=trees, characters = characters, mod=model, pi=pi)
#   }
#   
#   #print(Qs)
#   
#   if (class(trees) == "multiPhylo"){
#     #generating stochastic maps in parallel
#     maps<-list()
#     tmp<-list()
#     for (i in 1:length(Qs)){
#       tmp <- mclapply(1:sims, function(x) fastSimmap(tree = trees[[i]], x = characters, Q = Qs[[i]]$Q), mc.cores = cores)
#       if(!("multiSimmap"%in%class(tmp))) class(tmp)<-c("multiSimmap","multiPhylo", class(tmp))
#       maps[[i]]<-tmp
#       print(paste("iterating through", round((i/length(Qs))*100, 2), "% of the tree distribution"))
#     }
#     #code from http://blog.phytools.org/2017/11/running-makesimmap-in-parallel.html
#     maps <- do.call(c, maps)
#     if(!("multiSimmap"%in%class(maps))) class(maps)<-c("multiSimmap","multiPhylo", class(maps))
#   } else {
#     #generating stochastic maps in parallel
#     maps<-list()
#     #tmp<-list()
#     #for (i in 1:length(Qs)){
#     maps <- mclapply(1:sims, function(x) fastSimmap(tree = trees, x = characters, Q = Qs$Q), mc.cores = cores)
#     if(!("multiSimmap"%in%class(maps))) class(maps)<-c("multiSimmap","multiPhylo", class(maps))
#     #maps<-tmp
#   }
#   
#   return(maps)
# }


#alternative describe.simmmap function that 
#uses the matchNodes code from Eliot Miller
#much faster than matchNodes
#also parallelizes the summary function when using
#a reference tree
describe.simmap.parallel<-function (tree, ncpus, ...) {
  if (hasArg(plot)) 
    plot <- list(...)$plot
  else plot <- FALSE
  if (hasArg(check.equal)) 
    check.equal <- list(...)$check.equal
  else check.equal <- FALSE
  if (hasArg(message)) 
    message <- list(...)$message
  else message <- FALSE
  if (hasArg(ref.tree)) 
    ref.tree <- list(...)$ref.tree
  else ref.tree <- NULL
  if (inherits(tree, "multiPhylo")) {
    if (check.equal) {
      TT <- sapply(tree, function(x, y) sapply(y, all.equal.phylo, 
                                               x), y = tree)
      check <- all(TT)
      if (!check) 
        cat("Note: Some trees are not equal.\nA \"reference\" tree will be computed if none was provided.\n\n")
    }
    else check <- TRUE
    YY <- getStates(tree)
    states <- sort(unique(as.vector(YY)))
    if (is.null(ref.tree) && check) 
      ZZ <- t(apply(YY, 1, function(x, levels, Nsim) summary(factor(x, 
                                                                    levels))/Nsim, levels = states, Nsim = length(tree)))
    else {
      if (is.null(ref.tree)) {
        cat("No reference tree provided & some trees are unequal.\nComputing majority-rule consensus tree.\n")
        ref.tree <- consensus(tree, p = 0.5)
      }
      YYp <- matrix(NA, ref.tree$Nnode, length(tree), dimnames = list(1:ref.tree$Nnode + 
                                                                        Ntip(ref.tree), NULL))
      registerDoParallel(cores=ncpus)
      foreach(i=seq_along(1:length(tree)), .combine=c) %dopar% {
        M <- match_phylo_nodes(ref.tree, tree[[i]])
        jj <- sapply(M[, 2], function(x, y) if (x %in% 
                                                y) 
          which(as.numeric(y) == x)
          else NA, y = as.numeric(rownames(YY)))
        YYp[, i] <- YY[jj, i]
      }
      
      # for (i in 1:length(tree)) {
      #   M <- match_phylo_nodes(ref.tree, tree[[i]])
      #   jj <- sapply(M[, 2], function(x, y) if (x %in% 
      #                                           y) 
      #     which(as.numeric(y) == x)
      #     else NA, y = as.numeric(rownames(YY)))
      #   YYp[, i] <- YY[jj, i]
      # }
      ZZ <- t(apply(YYp, 1, function(x, levels) summary(factor(x[!is.na(x)], 
                                                               levels))/sum(!is.na(x)), levels = states))
    }
    XX <- countSimmap(tree, states, FALSE)
    XX <- XX[, -(which(as.vector(diag(-1, length(states))) == 
                         -1) + 1)]
    AA <- t(sapply(unclass(tree), function(x) c(colSums(x$mapped.edge), 
                                                sum(x$edge.length))))
    colnames(AA)[ncol(AA)] <- "total"
    BB <- getStates(tree, type = "tips")
    CC <- t(apply(BB, 1, function(x, levels, Nsim) summary(factor(x, 
                                                                  levels))/Nsim, levels = states, Nsim = length(tree)))
    x <- list(count = XX, times = AA, ace = ZZ, tips = CC, 
              tree = tree, ref.tree = if (!is.null(ref.tree)) ref.tree else NULL)
    class(x) <- "describe.simmap"
  }
  else if (inherits(tree, "phylo")) {
    XX <- countSimmap(tree, message = FALSE)
    YY <- getStates(tree)
    states <- sort(unique(YY))
    AA <- setNames(c(colSums(tree$mapped.edge), sum(tree$edge.length)), 
                   c(colnames(tree$mapped.edge), "total"))
    AA <- rbind(AA, AA/AA[length(AA)])
    rownames(AA) <- c("raw", "prop")
    x <- list(N = XX$N, Tr = XX$Tr, times = AA, states = YY, 
              tree = tree)
    class(x) <- "describe.simmap"
  }
  if (message) 
    print(x)
  if (plot) 
    plot(x)
  x
}

describe.simmap.alt<-function (tree, ...) {
  if (hasArg(plot)) 
    plot <- list(...)$plot
  else plot <- FALSE
  if (hasArg(check.equal)) 
    check.equal <- list(...)$check.equal
  else check.equal <- FALSE
  if (hasArg(message)) 
    message <- list(...)$message
  else message <- FALSE
  if (hasArg(ref.tree)) 
    ref.tree <- list(...)$ref.tree
  else ref.tree <- NULL
  if (inherits(tree, "multiPhylo")) {
    if (check.equal) {
      TT <- sapply(tree, function(x, y) sapply(y, all.equal.phylo, 
                                               x), y = tree)
      check <- all(TT)
      if (!check) 
        cat("Note: Some trees are not equal.\nA \"reference\" tree will be computed if none was provided.\n\n")
    }
    else check <- TRUE
    YY <- getStates(tree)
    states <- sort(unique(as.vector(YY)))
    if (is.null(ref.tree) && check) 
      ZZ <- t(apply(YY, 1, function(x, levels, Nsim) summary(factor(x, 
                                                                    levels))/Nsim, levels = states, Nsim = length(tree)))
    else {
      if (is.null(ref.tree)) {
        cat("No reference tree provided & some trees are unequal.\nComputing majority-rule consensus tree.\n")
        ref.tree <- consensus(tree, p = 0.5)
      }
      YYp <- matrix(NA, ref.tree$Nnode, length(tree), dimnames = list(1:ref.tree$Nnode + 
                                                                        Ntip(ref.tree), NULL))
      for (i in 1:length(tree)) {
        M <- match_phylo_nodes(ref.tree, tree[[i]])
        jj <- sapply(M[, 2], function(x, y) if (x %in% 
                                                y) 
          which(as.numeric(y) == x)
          else NA, y = as.numeric(rownames(YY)))
        YYp[, i] <- YY[jj, i]
      }
      ZZ <- t(apply(YYp, 1, function(x, levels) summary(factor(x[!is.na(x)], 
                                                               levels))/sum(!is.na(x)), levels = states))
    }
    XX <- countSimmap(tree, states, FALSE)
    XX <- XX[, -(which(as.vector(diag(-1, length(states))) == 
                         -1) + 1)]
    AA <- t(sapply(unclass(tree), function(x) c(colSums(x$mapped.edge), 
                                                sum(x$edge.length))))
    colnames(AA)[ncol(AA)] <- "total"
    BB <- getStates(tree, type = "tips")
    CC <- t(apply(BB, 1, function(x, levels, Nsim) summary(factor(x, 
                                                                  levels))/Nsim, levels = states, Nsim = length(tree)))
    x <- list(count = XX, times = AA, ace = ZZ, tips = CC, 
              tree = tree, ref.tree = if (!is.null(ref.tree)) ref.tree else NULL)
    class(x) <- "describe.simmap"
  }
  else if (inherits(tree, "phylo")) {
    XX <- countSimmap(tree, message = FALSE)
    YY <- getStates(tree)
    states <- sort(unique(YY))
    AA <- setNames(c(colSums(tree$mapped.edge), sum(tree$edge.length)), 
                   c(colnames(tree$mapped.edge), "total"))
    AA <- rbind(AA, AA/AA[length(AA)])
    rownames(AA) <- c("raw", "prop")
    x <- list(N = XX$N, Tr = XX$Tr, times = AA, states = YY, 
              tree = tree)
    class(x) <- "describe.simmap"
  }
  if (message) 
    print(x)
  if (plot) 
    plot(x)
  x
}



#run simmap in parallel
simmap.parallel<-function(trees, model, sims, characters, pi, ladderize, Q_cores, S_cores, preschedule=T){
  #generate a distribution of Q matricies from the posterior distribution
  if (class(trees)=="multiPhylo"){
    print(paste("detected", length(trees), "input trees, estimating Q matricies"))
    Qs <- Qs_loop.mc(trees=trees, characters = characters, mod=model, pi=pi, fun=Qs_loop, cores = Q_cores)
  } else {
    print("detected single input tree, estimating Q matrix")
    Qs <- Qs_loop(trees=trees, characters = characters, mod=model, pi=pi)
  }
  
  
  if (class(trees) == "multiPhylo"){
    print(paste("detected", length(trees), "input trees, estimating", sims, "stochastic maps on each"))
    if(pi=="equal"){
      print("using equal root state prior")
      #generating stochastic maps in parallel
      maps<-list()
      tmp<-list()
      for (i in 1:length(Qs)){
        tmp <- pbmclapply(1:sims, function(x) fastSimmap(tree = trees[[i]], x = characters, Q = Qs[[i]]), mc.cores = S_cores, mc.preschedule = preschedule)
        if(!("multiSimmap"%in%class(tmp))) class(tmp)<-c("multiSimmap","multiPhylo", class(tmp))
        maps[[i]]<-tmp
        tmp<-NULL
        print(paste("iterating through", round((i/length(Qs))*100, 2), "% of the tree distribution"))
        invisible(capture.output(gc(full=T)))
      }
      #code from http://blog.phytools.org/2017/11/running-makesimmap-in-parallel.html
      maps <- do.call(c, maps)
      #print("ladderizing output")
      #maps <- pbmclapply(maps, ladderize.simmap, mc.cores=cores)
      #if(!("multiSimmap"%in%class(maps))) class(maps)<-c("multiSimmap","multiPhylo", class(maps))
    }
    
    if(pi=="fitzjohn") {
      print("using Fitzjohn 2009 root state prior")
      #generating stochastic maps in parallel
      maps<-list()
      tmp<-list()
      for (i in 1:length(Qs)){
        tmp <- pbmclapply(1:sims, function(x) fastSimmap(tree = trees[[i]], x = characters, Q = Qs[[i]], pi="madfitz"), mc.cores = S_cores, mc.preschedule = preschedule)
        if(!("multiSimmap"%in%class(tmp))) class(tmp)<-c("multiSimmap","multiPhylo", class(tmp))
        maps[[i]]<-tmp
        tmp<-NULL
        print(paste("iterating through", round((i/length(Qs))*100, 2), "% of the tree distribution"))
        invisible(capture.output(gc(full=T)))
      }
      
      #code from http://blog.phytools.org/2017/11/running-makesimmap-in-parallel.html
      maps <- do.call(c, maps)
      #print("ladderizing output")
      #maps <- pbmclapply(maps, ladderize.simmap, mc.cores=cores)
      #if(!("multiSimmap"%in%class(maps))) class(maps)<-c("multiSimmap","multiPhylo", class(maps))
    }
    
  } else {
    print(paste("estimating", sims, "stochastic maps on single input tree"))
    if(pi=="equal"){
      print("using equal root state prior")
      #generating stochastic maps in parallel
      maps<-list()
      maps <- pbmclapply(1:sims, function(x) fastSimmap(tree = trees, x = characters, Q = Qs), mc.cores = S_cores, mc.preschedule = preschedule)
      #invisible(capture.output(gc()))
      #print("ladderizing output")
      #maps <- pbmclapply(maps, ladderize.simmap, mc.cores=cores)
      #if(!("multiSimmap"%in%class(maps))) class(maps)<-c("multiSimmap","multiPhylo", class(maps))
    }
    if(pi=="fitzjohn"){
      print("using Fitzjohn 2009 root state prior")
      #generating stochastic maps in parallel
      maps<-list()
      maps <- pbmclapply(1:sims, function(x) fastSimmap(tree = trees, x = characters, Q = Qs, pi="madfitz"), mc.cores = S_cores, mc.preschedule = preschedule)
      #invisible(capture.output(gc()))
      #print("ladderizing output")
      #maps <- pbmclapply(maps, ladderize.simmap, mc.cores=cores)
      #if(!("multiSimmap"%in%class(maps))) class(maps)<-c("multiSimmap","multiPhylo", class(maps))
    }
  }
  
  if(ladderize==T){
    print("ladderizing output")
    maps <- pbmclapply(maps, ladderize.simmap, mc.cores=Q_cores, mc.preschedule = preschedule)
    if(!("multiSimmap"%in%class(maps))) class(maps)<-c("multiSimmap","multiPhylo", class(maps))
    
  } else{
    if(!("multiSimmap"%in%class(maps))) class(maps)<-c("multiSimmap","multiPhylo", class(maps))
  }
  
  #print("ladderizing output")
  #maps <- mclapply(maps, ladderize.simmap, mc.cores=cores)
  #if(!("multiSimmap"%in%class(maps))) class(maps)<-c("multiSimmap","multiPhylo", class(maps))
  print("analysis complete, have a nice day")
  #browseURL(url="http://mrwgifs.com/wp-content/uploads/2013/05/Denzel-Washington-Boom-Gif.gif")
  return(maps)
  
}


#function for extracting an average Q matrix from a big run (posterior)
Q_extract<-function(simmap, every=500){
  result<-list()
  for(i in seq(1, length(simmap), every)){
    print(paste(round(i/length(simmap), 3)*100, "% complete"))
    result[[1+length(result)]]<-(simmap[[i]]$Q)
  }
  
  result<-apply(simplify2array(result), 1:2, mean)
  return(result)
}

#function for extracting all Qmatricies separately for inspection
Q_extract_all<-function(simmap, every=500){
  result<-list()
  for(i in seq(1, length(simmap), every)){
    print(paste(round(i/length(simmap), 3)*100, "% complete"))
    result[[1+length(result)]]<-(simmap[[i]]$Q)
  }
  
  #result<-apply(simplify2array(result), 1:2, mean)
  return(result)
}



#count number of zeros in vector
countZeros<-function(dat){
sum(dat==0)
}


#function for counting the number of appearances of zero
countQ_Zeros<-function(simmap, every=500){
  result<-list()
  for(i in seq(1, length(simmap), every)){
    #print(paste(round(i/length(simmap), 3)*100, "% complete"))
    result[[1+length(result)]]<-(simmap[[i]]$Q)
  }
  
  result<-apply(simplify2array(result), 1:2, countZeros)
  return(result)
  
}



#example
#simmap.parallel(trees=multiPhylo_object, model="ARD, ER, SYM OR MATRIX", sims=NUMBER OF SIMS, characters = single character vector, pi=roor prior assumption, cores=number of cores)
