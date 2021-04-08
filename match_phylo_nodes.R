
#function from Eliot Miller to generate the matches
#depends on geiger
tableSetup <- function(tree1, tree2){
  output <- data.frame(matrix(nrow=min(c(length(tree1$tip.label),length(tree2$tip.label)))-1, ncol=2))
  names(output) <- c("tree1", "tree2")
  output
}
equivCut <- function(tree1, tree2, results, least.inclusive){
  #identify nodes that appear more than once
  dups <- unique(results[,2][duplicated(results[,2])])
  
  #if there aren't any, return the original results
  if(length(dups)==0)
  {
    return(results)
  }
  
  #loop through, add details on size of clade to data frame
  for(i in 1:length(dups))
  {
    #set aside all the results that pertain to that node
    setAside <- results[results[,2]==dups[i],]
    
    #go through each node and figure out how many taxa descend from that node
    temp <- data.frame(matrix(nrow=dim(setAside)[1], ncol=2))
    names(temp) <- c("node","no.taxa")
    for(j in 1:dim(setAside)[1])
    {
      temp[j,"node"] <- setAside[j,"tree1"]
      temp[j,"no.taxa"] <- length(extract.clade(tree1, setAside[j,"tree1"])$tip.label)
    }
    
    #now you can decide which of these to keep
    if(least.inclusive)
    {
      keep <- temp$node[temp$no.taxa==min(temp$no.taxa)]
      toDrop <- setdiff(temp$node, keep)
      results <- results[!(results$tree1 %in% toDrop),]
    }
    
    else
    {
      keep <- temp$node[temp$no.taxa==max(temp$no.taxa)]
      toDrop <- setdiff(temp$node, keep)
      results <- results[!(results$tree1 %in% toDrop),]
    }
  }
  
  results
}
filter<-function(ref,results){
  output<-matrix(nrow=ref$Nnode, ncol=2)
  output[,1]<-(length(ref$tip.label)+1):(length(ref$tip.label)+length(ref$tip.label)-1)
  rownames(output)<-output[,1]
  rownames(results)<-results[,1]
  
  merged<-merge(output, results, by=0, all=T)
  merged$Row.names<-NULL
  merged$V2<-NULL
  merged$tree1<-NULL
  colnames(merged)<-c("tree1", "tree2")
  merged
}
match_phylo_nodes <- function(tree1, tree2, least.inclusive=TRUE){
  #set the results table up
  results <- tableSetup(tree1, tree2)
  
  #define tree1 nodes here
  nodes1 <- (length(tree1$tip.label)+1):max(tree1$edge)
  
  for(i in 1:length(nodes1))
  {
    #set the correct row in the tree1 column to the node in question
    results[i,1] <- nodes1[i]
    
    #find the descendants of this node
    tips1 <- tips(tree1, nodes1[i])
    
    #drop these to just tips that are in tree 2
    tips1 <- tips1[tips1 %in% tree2$tip.label]
    
    #if there's nothing left, set to no match and move on
    if(length(tips1)==0)
    {
      results[i,2] <- "no.match"
      next()
    }
    
    #find the MRCA of those tips in tree 2
    #if there's only a single taxon left, pull the node it descends from (getMRCA will fail)
    if(length(tips1==1))
    {
      edge2 <- which.edge(tree2, tips1)
      mrca2 <- tree2$edge[edge2,][1]
    }
    else
    {
      mrca2 <- getMRCA(tree2, tips1)
    }
    
    #find the descendants of that node
    tips2 <- tips(tree2, mrca2)
    
    #drop to just taxa that are in tree1
    tips2 <- tips2[tips2 %in% tree1$tip.label]
    
    #if these tips are equivalent, the nodes match
    if(setequal(tips1, tips2))
    {
      results[i,2] <- mrca2
    }
    
    else
    {
      results[i,2] <- "no.match"
    }
  }
  #print(results)
  #consider what you want to do next. do you want to, e.g., add a row for
  #every not-yet-mentioned tree2 node and add "no.match" to the first column,
  #or do you want to return as is, or do you want to only return matches? for
  #now, drop to only matches and re-class both as integer so can demonstrate
  #equality no matter which tree is first or second
  results <- results[results[,2] != "no.match",]
  results[,2] <- as.integer(results[,2])
  
  #run the equivCut function
  results <- equivCut(tree1, tree2, results, least.inclusive)
  
  results <- filter(ref=tree1, results)
  results
}
