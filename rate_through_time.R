

##this function is a modified version of markChanges that only returns the values 
markChanges.alt<-function (tree, colors = NULL, cex = 1, lwd = 2, plot = F) {
  states <- sort(unique(getStates(tree)))
  #if (is.null(colors)) 
    #colors <- setNames(palette()[1:length(states)], states)
  obj <- get("last_plot.phylo", envir = .PlotPhyloEnv)
  nc <- sapply(tree$maps, length) - 1
  ii <- which(nc > 0)
  nc <- nc[ii]
  xx <- yy <- vector()
  for (i in 1:length(ii)) {
    for (j in 1:nc[i]) {
      ss <- names(tree$maps[[ii[i]]])[j + 1]
      mm <- tree$edge[ii[i], 1]
      dd <- tree$edge[ii[i], 2]
      x <- rep(obj$xx[mm] + cumsum(tree$maps[[ii[i]]])[j], 2)
      y <- c(obj$yy[dd] - 0.5 * mean(strheight(LETTERS) * cex), obj$yy[dd] + 0.5 * mean(strheight(LETTERS) * cex))
      #if (plot) 
      #  lines(x, y, lwd = lwd, col = colors[ss], lend = 2)
      xx <- c(xx, setNames(x[1], paste(names(tree$maps[[ii[i]]])[j:(j + 1)], collapse = "->")))
      yy <- c(yy, mean(y))
    }
  }
  XY <- cbind(xx, yy)
  colnames(XY) <- c("x", "y")
  return(XY)
}

#this function computes the total number of changes through the list of stochastic maps
countchanges<-function(tree, maps){
  #use markChanges to identify which changes occur at what timepoints
  #print("compute bins")
  tmp<-plotTree(tree, plot=F)
  #par(fg="transparent")
  #tiplabels(pie=to.matrix(character1,c("Arboreal", "Non-arboreal", "Semi-arboreal")),piecol=colors,cex=0.2)
  #par(fg="black")
  plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n')
  text(x = 0.5, y = 0.5, paste("BRB"), cex=5)
  #print("counting changes")
  changes<-sapply(maps, markChanges.alt)
  dev.off()
  return(changes)
  #set up the temporal segments
}

#parse out changes of particular types
parsetypes<-function(summary, type) {
  tmp<-list()
  for(i in 1:length(summary)){
    tmp[[i]]<-summary[[i]][which(rownames(summary[[i]])==type),,drop = FALSE]
  }
  return(tmp)
}

#this function generates a new object containing the average across bins
getrates<-function(tree, bins, maps, type, log=F) { 
	#mostly from liam's blog
	#http://blog.phytools.org/2017/11/visualizing-rate-of-change-in-discrete.html
	#use markChanges to identify which changes occur at what timepoints
  #print("compute bins")
  #plotTree(tree,ftype="off",lwd=1)
  #par(fg="transparent")
  #tiplabels(pie=to.matrix(character1,c("Arboreal", "Non-arboreal", "Semi-arboreal")),piecol=colors,cex=0.2)
  #par(fg="black")
  #changes<-sapply(maps, markChanges, plot=F)
  #dev.off()
  if (type == "all"){
    changes<-countchanges(tree, maps)
  } else {
    changes<-countchanges(tree, maps)
    changes<-parsetypes(changes, type)
  }
  #return(changes)
	#set up the temporal segments
	h<-max(nodeHeights(tree))
	b<-bins
	segs<-cbind(seq(0,h-h/b,h/b), seq(1/b*h,h,h/b))
	#print(segs)
	
	#print("counting changes in bins")
	#compute the mean number of sampled changes for each segment
	nchanges<-rep(0,b)
	#print(nchanges)

	
	for(i in 1:length(changes)){
	    #print(paste("i is", i))
	    #print(length(changes[[i]]))
	    #return(changes[[136]])
	      
	    if(nrow(changes[[i]]) > 0) {
	    
	    #if(nrow(changes[[i]]) == 0) {
	    #  stop()
	    #}
	  
	      for(j in 1:nrow(changes[[i]])){
	        #print(paste("i is", i))
	        #print(paste("j is", j))
	        #print(1:nrow(changes[[145]]))
	        #print(changes[[i]])
	        ind<-which((changes[[i]][j,1]>segs[,1])+(changes[[i]][j,1]<=segs[,2])==2)
	        nchanges[ind]<-nchanges[ind]+1/length(changes)
	      }
	      
	    } else {
        next #do nothing, nchanges will remain vector of zeros
	    }
	    
	}
	
	#print(nchanges)
	
	#print("accounting for LTT")
	#control for the total edge length present in a segment using LTT computation
	LTT<-ltt(tree,plot=FALSE)
	
	LTT<-cbind(LTT$ltt[2:(length(LTT$ltt)-1)],
    LTT$times[2:(length(LTT$ltt)-1)],
    LTT$times[3:length(LTT$ltt)])
	
	ii<-1
	edge.length<-rep(0,b)
	
	#print(segs)
	for(i in 1:nrow(segs)){
		done.seg<-FALSE
		while(LTT[ii,2]<=segs[i,2]&&done.seg==FALSE){
			edge.length[i]<-edge.length[i]+
			LTT[ii,1]*(min(segs[i,2],LTT[ii,3])-
			max(segs[i,1],LTT[ii,2]))
		if(LTT[ii,3]>=segs[i,2]) done.seg<-TRUE
		if(LTT[ii,3]<=segs[i,2]) ii<-if(ii<nrow(LTT)) ii+1 else ii
    	}
	}
	
	#print bins thru time
	timesegs<-h-as.vector(t(segs))
	if(log==T) {
	  magnitude<-rbind(nchanges/log(edge.length),nchanges/log(edge.length))
	} else {
	  magnitude<-rbind(nchanges/edge.length,nchanges/edge.length)
	}
	
	
	
	out<-list()
	out[[1]]<-segs
  out[[2]]<-timesegs
  out[[3]]<-magnitude
	out[[4]]<-type
  names(out)<-c("segs", "timesegs", "magnitude", "type")
  
  return(out)

}

#this function generates the summary plot
rateplot<-function(rates, tree, ylim, spline, width, alpha){
  if (spline == T) {
    smoothed <- smooth.spline(x=jitter(rates$timesegs, factor=0.1), y=rates$magnitude, cv=T)
    newx<- seq(from = max(rates$timesegs), to=0, length.out = 100)
    pred<-predict(smoothed, x = newx)
    #print(c(max(rates$segs), min(rates$segs)))
    plot(x=newx, y=pred$y, 
         lwd=width, type='l', 
         xlim=c(max(rates$segs), min(rates$segs)), 
         ylab="mean number of changes / total edge length", 
         xlab="time since the present", 
         #main=rates$type, 
         ylim)
    
  } else {
    
    plot(rates$timesegs, rates$magnitude, 
         lwd=width, type="l", 
         xlim=c(max(rates$segs), min(rates$segs)),
         lend=0,xlab="time since the present",
         ylab="mean number of changes / total edge length", 
         #main=rates$type, 
         ylim)

  }
  
  plotTree(tree,add=TRUE,ftype="off",lwd=1,color=make.transparent("blue", alpha),
           mar=par()$mar,direction="leftwards",xlim=c(max(rates$segs),min(rates$segs)))
  
  #abline(v=66, col="red", lwd=2)
}

#only plot the segments
segplot<-function(rates, ylim, color, spline, width, ax=T, y.lab="mean number of changes / total edge length", x.lab="time since the present") {
  if (spline == T) {
    smoothed <- smooth.spline(x=jitter(rates$timesegs, factor=0.1), y=rates$magnitude, cv=T)
    newx<- seq(from = max(rates$timesegs), to=0, length.out = 100)
    pred<-predict(smoothed, x = newx)
    plot(x=newx, y=pred$y, 
         lwd=width, type='l', 
         xlim=c(max(rates$segs), min(rates$segs)), 
         ylab=y.lab, 
         xlab=x.lab, 
         #main=rates$type, 
         ylim,
         col=color,
         axes=ax)
    
  } else {
    
    plot(rates$timesegs, rates$magnitude, 
         lwd=width, type="l", 
         xlim=c(max(rates$segs), min(rates$segs)),
         lend=0,xlab="time since the present",
         ylab="mean number of changes / total edge length", 
         #main=rates$type, 
         ylim,
         col=color,
         axes=ax)
    
  }
}

#function for adding posterior splines to pre-existing plot
segplotloop<-function(rates, ylim, color, spline, width, xlim, ax=F) {
  if (spline == T) {
    for (i in 1:length(rates)){
      print(paste("plotting spline", i))
      par(new=T)
      #plot.window(xlim, ylim)
      smoothed <- smooth.spline(x=jitter(rates[[i]]$timesegs, factor=0.1), y=rates[[i]]$magnitude, cv=T)
      #print("smoothed")
      newx<- seq(from = max(rates[[i]]$timesegs), to=0, length.out = 100)
      #print("newx")
      pred<-predict(smoothed, x = newx)

      plot(x=newx, y=pred$y,
           lwd=width, type='l',
           xlim,
           ylab="",
           xlab="",
           ylim,
           col=color,
           bty="n",
           axes = ax,
           #xaxs='i'
           )
      
      # points(x=newx, y=pred$y, 
      #      lwd=width, type='l', 
      #      #xlim, 
      #      #ylab="", 
      #      #xlab="",
      #      #ylim,
      #      col=color#,
      #      #bty="n",
      #      #axes = T#,
      #      #xaxs='i'
      #     )
    }

  } else {
    for (i in 1:length(rates)){
      print(paste("plotting rates", i))
      par(new=T)
      plot(rates[[i]]$timesegs, rates[[i]]$magnitude, 
           lwd=width, type="l", 
           xlim, 
           ylab="", 
           xlab="",
           ylim,
           col=color,
           bty="n",
           axes = ax)
    }
      
  }
  
}

#function for splitting posterior map object into sub-lists
#this is used to prepare an object to pass to getratesloop
splitter<-function(maps, length){
  #https://stackoverflow.com/questions/18857275/r-given-a-list-return-a-list-of-equal-length-sublists
  n <- length(maps)
  k <- length ## your LEN
  tmp <- split(maps, rep(1:ceiling(n/k), each=k)[1:n])
}

#get the binned rates across the posterior trees
getratesloop<-function(posteriortrees, bins, splitmaps, type, log=F){
  tmp<-list()
  
  for (i in 1:length(posteriortrees)){
    print(paste("operating on tree", i))
    tmp[[i]] <- getrates(tree=force.ultrametric(as.phylo(posteriortrees[[i]])), bins, maps=splitmaps[[i]], type, log=F)

  }
  
  return(tmp)
}


# 
# https://archive.ph/c29Mw
# 
# {{cite web
#   | title       = Phylogenetic Tools for Comparative Biology: Visualizing the rate of c…
#   | url         = http://blog.phytools.org/2017/11/visualizing-rate-of-change-in-discrete.html
#   | date        = 2021-04-08
#   | archiveurl  = http://archive.today/c29Mw
#   | archivedate = 2021-04-08 }}