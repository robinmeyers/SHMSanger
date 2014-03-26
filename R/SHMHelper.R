getBases <- function () {
  bases <- c("A","C","G","T","N")
}
getAscii <- function () {
  ascii <- c(65,67,71,84,78)
}
getBasecolors <- function () {
  basecolors <- brewer.pal(7,"Set1")
}

tictactoePlot <- function (subs, dels, blocks, ref, tstart, tend, plotrows, cloneIDs) {
  
  bases <- getBases()
  ascii <- getAscii()
  basecolors <- getBasecolors()
  refseq <- paste(ref$Base,collapse="")
  par(mai=c(0.2,0.75,0.2,0.75),omi=c(0.5,0,0,0))
  
  layout(as.matrix(1:plotrows,ncol=1,nrow=plotrows))
  
  
  ymax <- max(5.5,length(cloneIDs))
  
  
  dels$y <- match(dels$Clone,cloneIDs)
  
  
  subs$y <- match(subs$Clone,cloneIDs)
  subs$color <- match(subs$To,bases)
  subs$pch <- match(subs$To,bases)
  
  
  
  

  
  agct <- unlist(gregexpr("AGCT",as.character(refseq)))
  rgyw <- unique(c(unlist(gregexpr("[AG]G[CT][AT]",as.character(refseq))),unlist(gregexpr("[AT][AG]C[CT]",as.character(refseq)))))
  
  rowwidth <- ceiling((tend-tstart+1)/plotrows)
  
  tstarts <- tstart+rowwidth*0:(plotrows-1)
  tends <- tstarts+rowwidth-1
  
  
  
  for (i in 1:plotrows) {
    
    plot(c(),c(),ylab="",xlab="",xaxt="n",xlim=c(max(1,tstarts[i]-2),min(nchar(refseq),tends[i]+2)),ylim=c(0,ymax),xaxs="i",bty="o")
    axis(1,lwd=0,lwd.ticks=1)
    
    rect(xleft=rgyw-0.5,ybottom=-1,xright=rgyw+3.5,ytop=ymax,col=rgb(254,217,142,max=255),border=F)
    rect(xleft=agct-0.5,ybottom=-1,xright=agct+3.5,ytop=ymax,col=rgb(254,153,41,max=255),border=F)
    
    
    grid(ny=0,col=grey(0.5),lty=3)
    points(1:nrow(ref),rep(0,nrow(ref)),col=basecolors[ref$color],pch=ascii[ref$pch],cex=0.6)
    if (length(cloneIDs) > 0) {
      rect(xleft=blocks$Start-0.5,ybottom=blocks$Clone-0.5,xright=blocks$End+0.5,ytop=blocks$Clone+0.5,col=grey(0.1,0.25),border=F)
      points(subs$Pos,subs$y,col=basecolors[subs$color],pch=ascii[subs$pch],cex=0.5)
      segments(x0=dels$Pos-0.5,y0=dels$y,x1=dels$End+0.5,col=basecolors[7],lwd=1)
    } else {
      text(x=(tstarts[i]+tends[i])/2,y=ymax/2-0.5,"No Mutations To Display")
    }
  }
  
}

calculateBlocks <- function(clones,tstart,tend) {
  blocks <- data.frame(Clone=character(),Start=integer(),End=integer(),stringsAsFactors=F)

  if (nrow(clones) > 0) {
    for (i in 1:nrow(clones)) {
      coordlist <- unlist(strsplit(clones[i,"Coords"],","))
      
      for (j in 1:length(coordlist)) {
        coords <- as.integer(unlist(strsplit(coordlist[j],"-")))
        if (j == 1 && coords[1] > tstart) {
          blocks[nrow(blocks)+1,] <- c(clones[i,"Clone"],tstart,coords[1]-1)
        }
        if (j > 1) {
          coordslast <- as.integer(unlist(strsplit(coordlist[j-1],"-")))
          blocks[nrow(blocks)+1,] <- c(clones[i,"Clone"],coordslast[2]+1,coords[1]-1)
        }
        if (j == length(coordlist) && coords[2] < tend) {
          blocks[nrow(blocks)+1,] <- c(clones[i,"Clone"],coords[2]+1,tend)
        }
      }
      
    }
  }
  blocks$Start <- as.integer(blocks$Start)
  blocks$End <- as.integer(blocks$End)
  return(blocks)
}

calculateProfile <- function (subs, clones, refseq, tstart, tend) {
  
  blocks <- calculateBlocks(clones,tstart,tend)
  profile <- data.frame(Pos=tstart:tend,Base=strsplit(substr(as.character(refseq),tstart,tend),""))
  colnames(profile) <- c("Pos","Base")
  profile$Clones <- 0
  profile$Subs <- 0
  profile$Y <- 0

  if (nrow(subs) > 0) {
    profile$Subs <- hist(subs$Pos,breaks=seq(tstart-0.5,tend+0.5,by=1),plot=F)$counts
  }

  profile$Clones <- nrow(clones) - unlist(lapply(tstart:tend,function(x) {
      return(nrow(blocks[blocks$Start <= x & blocks$End >= x,]))
    }))

  profile$Y <- ifelse(profile$Clones == 0, 0, profile$Subs/profile$Clones)

  return(profile)
}

connectfourSubPlot <- function (subs, blocks, ref, tstart, tend, plotrows, cloneIDs) {
  
  bases <- getBases()
  ascii <- getAscii()
  basecolors <- getBasecolors()
  refseq <- paste(ref$Base,collapse="")
  par(mai=c(0.2,0.75,0.2,0.75),omi=c(0.5,0,0,0))
  
  layout(as.matrix(1:plotrows,ncol=1,nrow=plotrows))
  
  agct <- unlist(gregexpr("AGCT",as.character(refseq)))
  rgyw <- unique(c(unlist(gregexpr("[AG]G[CT][AT]",as.character(refseq))),unlist(gregexpr("[AT][AG]C[CT]",as.character(refseq)))))
  rowwidth <- ceiling((tend-tstart+1)/plotrows)
  
  tstarts <- tstart+rowwidth*0:(plotrows-1)
  tends <- tstarts+rowwidth-1
  
  
  ref$dens_denom <- 0
  ref$dens_numer <- 0
  ref$dens <- 0
  
  if (nrow(subs) > 0) {
    ref$dens_numer <- hist(subs$Pos,0:nrow(ref),plot=F)$counts
  }
  ref$dens_denom <- length(cloneIDs) - unlist(lapply(1:nrow(ref),function(x) {
    return(nrow(blocks[blocks$Start <= ref$Pos[x] & blocks$End >= ref$Pos[x],]))
  }))
  
  ref$dens <- ifelse(ref$dens_denom==0,0,ref$dens_numer/ref$dens_denom)
  
  dens <- data.frame(x=c(ref$Pos-0.49,ref$Pos+0.49),y=c(ref$dens,ref$dens))
  dens <- dens[order(dens$x),]
  ymax <- max(0.05,1.25*dens$y)
  
  refy <- -ymax/20
  
  #   rgyw_tops <- ymax
  #   agct_tops <- ymax
  #   rgyw_bottoms <- 0.82*ymax
  #   agct_bottoms <- 0.82*ymax
  
  #   Bottom of the plot - below the sequence  
  rgyw_bottoms <- 2*refy
  agct_bottoms <- 2*refy
  
  #   Adjust to height of the peak
  rgyw_tops <- unlist(lapply(rgyw,function(x,dens) {max(dens[x:(x+3)])},ref$dens)) + 0.025*ymax
  agct_tops <- unlist(lapply(agct,function(x,dens) {max(dens[x:(x+3)])},ref$dens)) + 0.025*ymax
  
  
  
  for (i in 1:plotrows) {
    plot(c(),c(),ylab="",xaxt="n",xlab="",xlim=c(max(1,tstarts[i]-2),min(nchar(refseq),tends[i]+2)),ylim=c(refy,ymax),xaxs="i",bty="o")
    axis(1,lwd=0,lwd.ticks=1)
    rect(xleft=rgyw-0.5,ybottom=rgyw_bottoms,xright=rgyw+3.5,ytop=rgyw_tops,col=rgb(254,217,142,max=255),border=rgb(254,217,142,max=255),lwd=2)
    rect(xleft=agct-0.5,ybottom=agct_bottoms,xright=agct+3.5,ytop=agct_tops,col=rgb(254,153,41,max=255),border=rgb(254,153,41,max=255),lwd=2)
    
    grid(col=grey(0.5))
    points(1:nrow(ref),rep(refy,nrow(ref)),col=basecolors[ref$color],pch=ascii[ref$pch],cex=0.6)
    lines(dens$x,dens$y)
  }

  return(ref)
}

connectfourDelPlot <- function (dels, blocks, ref, tstart, tend, plotrows, cloneIDs) {
  
  bases <- getBases()
  ascii <- getAscii()
  basecolors <- getBasecolors()
  refseq <- paste(ref$Base,collapse="")
  par(mai=c(0.2,0.75,0.2,0.75),omi=c(0.5,0,0,0))
  
  layout(as.matrix(1:plotrows,ncol=1,nrow=plotrows))
  
  agct <- unlist(gregexpr("AGCT",as.character(refseq)))
  rgyw <- unique(c(unlist(gregexpr("[AG]G[CT][AT]",as.character(refseq))),unlist(gregexpr("[AT][AG]C[CT]",as.character(refseq)))))
  rowwidth <- ceiling((tend-tstart+1)/plotrows)
  
  tstarts <- tstart+rowwidth*0:(plotrows-1)
  tends <- tstarts+rowwidth-1
  
  
  ref$dens_denom <- 0
  ref$dens_numer <- 0
  ref$dens <- 0
  
  if (nrow(dels) > 0) {
    dels_expand <- c()
    for (i in 1:nrow(dels)) {
      dels_expand <- c(dels_expand,seq(dels$Pos[i],dels$Pos[i]+dels$Size[i]-1))
    }
    ref$dens_numer <- hist(dels_expand,c(0,1:nrow(ref)),plot=F)$counts
  }
  
  
  ref$dens_denom <- length(cloneIDs) - unlist(lapply(1:nrow(ref),function(x) {
    return(nrow(blocks[blocks$Start <= ref$Pos[x] & blocks$End >= ref$Pos[x],]))
  }))
  
  ref$dens <- ifelse(ref$dens_denom==0,0,ref$dens_numer/ref$dens_denom)
  
  dens <- data.frame(x=c(ref$Pos-0.49,ref$Pos+0.49),y=c(ref$dens,ref$dens))
  dens <- dens[order(dens$x),]
  ymax <- max(0.01,1.25*dens$y)
  
  refy <- -ymax/20
  
  #   rgyw_tops <- ymax
  #   agct_tops <- ymax
  #   rgyw_bottoms <- 0.82*ymax
  #   agct_bottoms <- 0.82*ymax
  
  #   Bottom of the plot - below the sequence  
  rgyw_bottoms <- 2*refy
  agct_bottoms <- 2*refy
  
  #   Adjust to height of the peak
  rgyw_tops <- unlist(lapply(rgyw,function(x,dens) {max(dens[x:(x+3)])},ref$dens)) + 0.025*ymax
  agct_tops <- unlist(lapply(agct,function(x,dens) {max(dens[x:(x+3)])},ref$dens)) + 0.025*ymax
  
  
  
  for (i in 1:plotrows) {
    plot(c(),c(),ylab="",xaxt="n",xlab="",xlim=c(max(1,tstarts[i]-2),min(nchar(refseq),tends[i]+2)),ylim=c(refy,ymax),xaxs="i",bty="o")
    axis(1,lwd=0,lwd.ticks=1)
    rect(xleft=rgyw-0.5,ybottom=rgyw_bottoms,xright=rgyw+3.5,ytop=rgyw_tops,col=rgb(254,217,142,max=255),border=rgb(254,217,142,max=255),lwd=2)
    rect(xleft=agct-0.5,ybottom=agct_bottoms,xright=agct+3.5,ytop=agct_tops,col=rgb(254,153,41,max=255),border=rgb(254,153,41,max=255),lwd=2)
    
    grid(col=grey(0.5))
    points(1:nrow(ref),rep(refy,nrow(ref)),col=basecolors[ref$color],pch=ascii[ref$pch],cex=0.6)
    lines(dens$x,dens$y)
  }
  
  
}