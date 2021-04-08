
#function for adding clade arc labels to upham tree
addcladelabels<-function(ln.offset=1.175,lab.offset=1.225){
  arc.cladelabels(tree=TimeTree,'Marsupialia', node=findMRCA(TimeTree, tips=c('Acrobates','Rhyncholestes')), ln.offset=ln.offset, lab.offset=lab.offset, cex=1, mark.node=F, lwd=3, lend=1, col='grey48')
  arc.cladelabels(tree=TimeTree,'Monotremata', node=findMRCA(TimeTree, tips=c('Ornithorhynchus','Tachyglossus')), ln.offset=ln.offset, lab.offset=lab.offset, cex=1, mark.node=F, lwd=3, lend=1, col='grey48')
  arc.cladelabels(tree=TimeTree,'Afrotheria', node=findMRCA(TimeTree, tips=c('Amblysomus','Trichechus')), ln.offset=ln.offset, lab.offset=lab.offset, cex=1, mark.node=F, lwd=3, lend=1, col='grey48')
  arc.cladelabels(tree=TimeTree,'Xenarthra', node=findMRCA(TimeTree, tips=c('Bradypus','Dasypus')), ln.offset=ln.offset, lab.offset=lab.offset, cex=1, mark.node=F, lwd=3, lend=1, col='grey48')
  arc.cladelabels(tree=TimeTree,'Euarchonta', node=findMRCA(TimeTree, tips=c('Macaca','Galeopterus')), ln.offset=ln.offset, lab.offset=lab.offset, cex=1, mark.node=F, lwd=3, lend=1, col='grey48')
  arc.cladelabels(tree=TimeTree,'Glires', node=findMRCA(TimeTree, tips=c('Oryctolagus','Glis')), ln.offset=ln.offset, lab.offset=lab.offset, cex=1, mark.node=F, lwd=3, lend=1, col='grey48')
  arc.cladelabels(tree=TimeTree,'Laurasiatheria', node=findMRCA(TimeTree, tips=c('Pteropus','Solenodon')), ln.offset=ln.offset, lab.offset=lab.offset, cex=1, mark.node=F, lwd=3, lend=1, col='grey48')
}

#function for adding clade labels to rate plot tree upham
addcladelabels.flat<-function(ln.offset=1.175, wing=0.5, size=0.5){
  cladelabels(tree=TimeTree,'Mo.', node=findMRCA(TimeTree, tips=c('Ornithorhynchus','Tachyglossus')), offset=ln.offset, wing.length = wing, cex=size, orientation="horizontal")
  cladelabels(tree=TimeTree,'Ma.', node=findMRCA(TimeTree, tips=c('Acrobates','Rhyncholestes')), offset=ln.offset, wing.length = wing, cex=size, orientation="horizontal")
  cladelabels(tree=TimeTree,'Xe.', node=findMRCA(TimeTree, tips=c('Bradypus','Dasypus')), offset=ln.offset, wing.length = wing, cex=size, orientation="horizontal")
  cladelabels(tree=TimeTree,'Af.', node=findMRCA(TimeTree, tips=c('Amblysomus','Trichechus')), offset=ln.offset, wing.length = wing, cex=size, orientation="horizontal")
  cladelabels(tree=TimeTree,'Eu.', node=findMRCA(TimeTree, tips=c('Macaca','Galeopterus')), offset=ln.offset, wing.length = wing, cex=size, orientation="horizontal")
  cladelabels(tree=TimeTree,'Gl.', node=findMRCA(TimeTree, tips=c('Oryctolagus','Glis')), offset=ln.offset, wing.length = wing, cex=size, orientation="horizontal")
  cladelabels(tree=TimeTree,'La.', node=findMRCA(TimeTree, tips=c('Pteropus','Solenodon')), offset=ln.offset, wing.length = wing, cex=size, orientation="horizontal")
}

#function for adding clade arc labels for mammalian groups on the meredith tree
addcladelabels.m<-function(){
  arc.cladelabels(tree=TimeTree,'Marsupialia', node=168, ln.offset=1.20, lab.offset=1.25, cex=1, mark.node=F, lwd=3, lend=1, col='grey48')
  arc.cladelabels(tree=TimeTree,'Monotremata', node=166, ln.offset=1.20, lab.offset=1.25, cex=1, mark.node=F, lwd=3, lend=1, col='grey48')
  arc.cladelabels(tree=TimeTree,'Afrotheria', node=317, ln.offset=1.20, lab.offset=1.25, cex=1, mark.node=F, lwd=3, lend=1, col='grey48')
  arc.cladelabels(tree=TimeTree,'Xenarthra', node=313, ln.offset=1.20, lab.offset=1.25, cex=1, mark.node=F, lwd=3, lend=1, col='grey48')
  arc.cladelabels(tree=TimeTree,'Primatomorpha', node=298, ln.offset=1.20, lab.offset=1.25, cex=1, mark.node=F, lwd=3, lend=1, col='grey48')
  arc.cladelabels(tree=TimeTree,'Glires', node=264, ln.offset=1.20, lab.offset=1.25, cex=1, mark.node=F, lwd=3, lend=1, col='grey48')
  arc.cladelabels(tree=TimeTree,'Laurasiatheria', node=193, ln.offset=1.20, lab.offset=1.25, cex=1, mark.node=F, lwd=3, lend=1, col='grey48')
  arc.cladelabels(tree=TimeTree,'Scandentia', node=263, ln.offset=1.20, lab.offset=1.25, cex=1, mark.node=F, lwd=3, lend=1, col='grey48')
}

#function for adding clade labels to rate meredith tree
addcladelabels.flat.m<-function(ln.offset=1.175, wing=0.25, size=0.5){
  cladelabels(tree=TimeTree,'Ma.', node=168, offset=ln.offset, wing.length = wing, cex=size, orientation="horizontal")
  cladelabels(tree=TimeTree,'Mo.', node=166, offset=ln.offset, wing.length = wing, cex=size, orientation="horizontal")
  cladelabels(tree=TimeTree,'Af.', node=317, offset=ln.offset, wing.length = wing, cex=size, orientation="horizontal")
  cladelabels(tree=TimeTree,'Xe.', node=313, offset=ln.offset, wing.length = wing, cex=size, orientation="horizontal")
  cladelabels(tree=TimeTree,'Pr.', node=298, offset=ln.offset, wing.length = wing, cex=size, orientation="horizontal")
  cladelabels(tree=TimeTree,'Gl.', node=264, offset=ln.offset, wing.length = wing, cex=size, orientation="horizontal")
  cladelabels(tree=TimeTree,'La.', node=193, offset=ln.offset, wing.length = wing, cex=size, orientation="horizontal")
  cladelabels(tree=TimeTree,'Sc.', node=263, offset=ln.offset, wing.length = wing, cex=size, orientation="horizontal")
}

#modified plotting function for the Q matrix inset
plot.Qmatrix<-function (x, ...) {
  Q <- unclass(x)
  if (hasArg(signif)) 
    signif <- list(...)$signif
  else signif <- 3
  if (hasArg(main)) 
    main <- list(...)$main
  else main <- NULL
  if (hasArg(cex.main)) 
    cex.main <- list(...)$cex.main
  else cex.main <- 1.2
  if (hasArg(cex.traits)) 
    cex.traits <- list(...)$cex.traits
  else cex.traits <- 1
  if (hasArg(cex.rates)) 
    cex.rates <- list(...)$cex.rates
  else cex.rates <- 0.6
  if (hasArg(show.zeros)) 
    show.zeros <- list(...)$show.zeros
  else show.zeros <- TRUE
  if (hasArg(tol)) 
    tol <- list(...)$tol
  else tol <- 1e-06
  if (hasArg(mar)) 
    mar <- list(...)$mar
  #else mar <- c(1.1, 1.1, 3.1, 1.1)
  else mar <- c(1, 0.5, 0, 1)
  
  if (hasArg(lwd)) 
    lwd <- list(...)$lwd
  else lwd <- 1
  if (hasArg(umbral)) 
    umbral <- list(...)$umbral
  else umbral <- FALSE
  if (hasArg(ncat)) 
    ncat <- list(...)$ncat
  else ncat <- NULL
  if (hasArg(spacer)) 
    spacer <- list(...)$spacer
  else spacer <- 0.2
  plot.new()
  par(mar = mar, oma=c(0,0,0,0))
  xylim <- c(-1.2, 1.2)
  plot.window(xlim = xylim, ylim = xylim, asp = 1)
  #plot.window(xlim = c(0, 1.2), ylim = c(0, 1.2), asp = 1)
  if (!is.null(main)) 
    title(main = main, cex.main = cex.main)
  nstates <- nrow(Q)
  if (!umbral || is.null(ncat)) {
    step <- 360/nstates
    angles <- seq(0, 360 - step, by = step)/180 * pi
    if (nstates == 2) 
      angles <- angles + pi/2
    v.x <- cos(angles)
    v.y <- sin(angles)
  }
  else {
    v.x <- v.y <- vector()
    for (i in 1:length(ncat)) {
      Q <- Q[sort(rownames(Q)), sort(colnames(Q))]
      xp <- -1 + 2 * (i - 1)/(length(ncat) - 1)
      v.x <- c(v.x, rep(xp, ncat[i]))
      yp <- seq(1, -1, length.out = max(ncat))[1:ncat[i]]
      v.y <- c(v.y, yp)
    }
  }
  for (i in 1:nstates) for (j in 1:nstates) if (if (!isSymmetric(Q)) 
    i != j
    else i > j) {
    dx <- v.x[j] - v.x[i]
    dy <- v.y[j] - v.y[i]
    slope <- abs(dy/dx)
    shift.x <- 0.05 * sin(atan(dy/dx)) * sign(j - i) * if (dy/dx > 
                                                           0) 
      1
    else -1
    shift.y <- 0.05 * cos(atan(dy/dx)) * sign(j - i) * if (dy/dx > 
                                                           0) 
      -1
    else 1
    s <- c(v.x[i] + spacer * cos(atan(slope)) * sign(dx) + 
             if (isSymmetric(Q)) 0 else shift.x, v.y[i] + spacer * 
             sin(atan(slope)) * sign(dy) + if (isSymmetric(Q)) 0 else shift.y)
    e <- c(v.x[j] + spacer * cos(atan(slope)) * sign(-dx) + 
             if (isSymmetric(Q)) 0 else shift.x, v.y[j] + spacer * 
             sin(atan(slope)) * sign(-dy) + if (isSymmetric(Q)) 0 else shift.y)
    if (show.zeros || Q[i, j] > tol) {
      if (abs(diff(c(i, j))) == 1 || abs(diff(c(i, j))) == 
          (nstates - 1)) 
        text(mean(c(s[1], e[1])) + 1.5 * shift.x, mean(c(s[2], 
                                                         e[2])) + 1.5 * shift.y, round(Q[i, j], signif), 
             cex = cex.rates, srt = atan(dy/dx) * 180/pi)
      else text(mean(c(s[1], e[1])) + 0.3 * diff(c(s[1], 
                                                   e[1])) + 1.5 * shift.x, mean(c(s[2], e[2])) + 
                  0.3 * diff(c(s[2], e[2])) + 1.5 * shift.y, round(Q[i, 
                                                                     j], signif), cex = cex.rates, srt = atan(dy/dx) * 
                  180/pi)
      arrows(s[1], s[2], e[1], e[2], length = 0.05, code = if (isSymmetric(Q)) 
        3
        else 2, lwd = lwd)
    }
  }
  text(v.x, v.y, rownames(Q), cex = cex.traits, col = make.transparent("black", 
                                                                       0.7))
  object <- data.frame(states = rownames(Q), x = v.x, y = v.y)
  invisible(object)
}

#helper function to plot PDF trees for upham
pdfgenerator<-function(filename,
                       width=10,
                       height=10,
                       basetree=TimeTree,
                       geo=obj,
                       rad=r,
                       simsum,
                       tps=tips,
                       title,
                       modelfit,
                       posmodelfit,
                       cls=cols,
                       posterior=F) {
  
  #start plotting
  pdf(file=filename, width=width, height=height)
  
  #generate initial base tree plot
  plotTree(basetree, type='fan', fsize=0.5, ftype='i', lwd=2.5, offset=3)
  
  #apply the time periods 
  for(i in 1:nrow(geo$leg)){
    color<-paste(strsplit(geo$colors[i],"")[[1]][1:7],collapse="")
    draw.circle(0,0,radius=rad[i],col=color,border="transparent")
    draw.circle(0,0, radius=rad[i:min(i)+1], col='transparent', border='black', lty=2, lwd=0.25)
  }
  
  #add a circle for kpg
  draw.circle(0,0, radius=rad[4], col='transparent', border='red', lty=1, lwd=3)
  
  #start a new plot
  par(new=T)
  plotTree(basetree, type='fan', fsize=0.5, ftype='i', lwd=2.5, offset=3)
  addcladelabels()
  
  #adjust pie circle outline width
  par(lwd = 0.2)
  #add node labels
  nodelabels(pie=simsum$ace, piecol=cls, prompt=FALSE, x=0.9*par()$usr[1], y=-max(nodeHeights(basetree)), cex= 0.3)
  
  #add tip labels # can only be executed after pars recon for this sectino
  tiplabels(pie=tps, piecol=cls, prompt=FALSE, cex=0.2)
  
  #add legend
  add.simmap.legend(colors=cls,prompt=FALSE,x=0.9*par()$usr[1], y=-max(nodeHeights(basetree)), fsize=0.9)
  
  #add labels
  legend("topleft", legend = paste(title,"\nAIC = ", round(AIC(modelfit), 2)), bty='n', cex=1)
  
  #add Q matrix plot
  par(new=T)
  par(fig = c(0.8, 1.0, 0.8, 1.0)) 
  if (posterior==F) {
    plot.Qmatrix(as.Qmatrix(modelfit), cex.rates=0.5, cex.traits=0.6, show.zeros=T, spacer=0.25)
    title("Q matrix", line = -1, cex.main=0.5)
  }
  if (posterior==T) {
    plot.Qmatrix((posmodelfit), cex.rates=0.5, cex.traits=0.6, show.zeros=T, spacer=0.25)
    title("Average Q matrix", line = -1, cex.main=0.5)
  }
  
  dev.off() 
  
}




#helper function to plot PDF trees for meredith
pdfgenerator.m<-function(filename,
                       width=10,
                       height=10,
                       basetree=TimeTree,
                       geo=obj,
                       rad=r,
                       simsum,
                       tps=tips,
                       title,
                       modelfit,
                       posmodelfit,
                       cls=cols,
                       posterior=F) {
  
  #start plotting
  pdf(file=filename, width=width, height=height)
  
  #generate initial base tree plot
  plotTree(basetree, type='fan', fsize=0.5, ftype='i', lwd=2.5, offset=3)
  
  #apply the time periods 
  for(i in 1:nrow(geo$leg)){
    color<-paste(strsplit(geo$colors[i],"")[[1]][1:7],collapse="")
    draw.circle(0,0,radius=rad[i],col=color,border="transparent")
    draw.circle(0,0, radius=rad[i:min(i)+1], col='transparent', border='black', lty=2, lwd=0.25)
  }
  
  #add a circle for kpg
  draw.circle(0,0, radius=rad[4], col='transparent', border='red', lty=1, lwd=3)
  
  #start a new plot
  par(new=T)
  plotTree(basetree, type='fan', fsize=0.5, ftype='i', lwd=2.5, offset=3)
  addcladelabels.m()
  
  #adjust pie circle outline width
  par(lwd = 0.2)
  #add node labels
  nodelabels(pie=simsum$ace, piecol=cls, prompt=FALSE, x=0.9*par()$usr[1], y=-max(nodeHeights(basetree)), cex= 0.3)
  
  #add tip labels # can only be executed after pars recon for this sectino
  tiplabels(pie=tps, piecol=cls, prompt=FALSE, cex=0.2)
  
  #add legend
  add.simmap.legend(colors=cls,prompt=FALSE,x=0.9*par()$usr[1], y=-max(nodeHeights(basetree)), fsize=0.9)
  
  #add labels
  legend("topleft", legend = paste(title,"\nAIC = ", round(AIC(modelfit), 2)), bty='n', cex=1)
  
  #add Q matrix plot
  par(new=T)
  par(fig = c(0.8, 1.0, 0.8, 1.0)) 
  if (posterior==F) {
    plot.Qmatrix(as.Qmatrix(modelfit), cex.rates=0.5, cex.traits=0.6, show.zeros=T, spacer=0.25)
    title("Q matrix", line = -1, cex.main=0.5)
  }
  if (posterior==T) {
    plot.Qmatrix((posmodelfit), cex.rates=0.5, cex.traits=0.6, show.zeros=T, spacer=0.25)
    title("Average Q matrix", line = -1, cex.main=0.5)
  }
  
  dev.off() 
  
}


