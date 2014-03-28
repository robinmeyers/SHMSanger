ARGS <- c(
  "datafile","character","file path of data file - columns: Pos, Base, Y, Err",
  "output","character","file path of plot pdf"
  
)

OPTS <- c(
  "tstart","numeric",0,"Start of reference to view",
  "tend","numeric",0,"End of reference to view",
  "plotrows","numeric",4,"Rows on plot",
  "ymax","numeric",0,"Maximum y-axis height"
)

source_local <- function(fname){
  argv <- commandArgs(trailingOnly = FALSE)
  base_dir <- dirname(substring(argv[grep("--file=", argv)], 8))
  source(paste(base_dir, fname, sep="/"))
}

source_local("Rsub.R")
source_local("SHMHelper.R")

parseArgs("SHMPlot.R", ARGS, OPTS)

library(RColorBrewer)
bases <- getBases()
basecolors <- getBasecolors()
ascii <- getAscii()

data <- read.delim(datafile, header=T, as.is=T)
refseq <- paste(data$Base,collapse="")

if (tstart < 1) {
  tstart <- 1
}
if (tend == 0 || tend > nrow(data)) {
  tend <- nrow(data)
}

data$Style <- match(data$Base,bases)

agct <- unlist(gregexpr("AGCT",as.character(refseq)))
rgyw <- unique(c(unlist(gregexpr("[AG]G[CT][AT]",as.character(refseq))),unlist(gregexpr("[AT][AG]C[CT]",as.character(refseq)))))
rowwidth <- ceiling((tend-tstart+1)/plotrows)

tstarts <- tstart+rowwidth*0:(plotrows-1)
tends <- tstarts+rowwidth-1


plotline <- data.frame(x=c(data$Pos-0.49,data$Pos+0.49),y=c(data$Y,data$Y))
plotline <- plotline[order(plotline$x),]
if (ymax==0) { ymax <- 1.1*max(plotline$y[plotline$x >= tstart & plotline$x <= tend]) }
refy <- -ymax/20

if (rgyw[1] > 0) {
  rgyw_bottoms <- 2*refy
  rgyw_tops <- unlist(lapply(rgyw,function(x,y) { max(y[x:(x+3)]) }, data$Y)) + 0.025*ymax
}
if (agct[1] > 0) {
  agct_bottoms <- 2*refy
  agct_tops <- unlist(lapply(agct,function(x,y) { max(y[x:(x+3)]) }, data$Y)) + 0.025*ymax
}




pdf(output,height=8,width=11)

  par(mai=c(0.2,0.75,0.2,0.75),omi=c(0.5,0,0,0))
  layout(as.matrix(1:plotrows,ncol=1,nrow=plotrows))

  for (i in 1:plotrows) {
    plot(c(),c(),ylab="",xaxt="n",xlab="",xlim=c(max(1,tstarts[i]-2),min(nchar(refseq),tends[i]+2)),ylim=c(refy,ymax),xaxs="i",bty="o")
    axis(1,lwd=0,lwd.ticks=1)
    if (rgyw[1] > 0) rect(xleft=rgyw-0.5,ybottom=rgyw_bottoms,xright=rgyw+3.5,ytop=rgyw_tops,col=rgb(254,217,142,max=255),border=rgb(254,217,142,max=255),lwd=2)
    if (agct[1] > 0) rect(xleft=agct-0.5,ybottom=agct_bottoms,xright=agct+3.5,ytop=agct_tops,col=rgb(254,153,41,max=255),border=rgb(254,153,41,max=255),lwd=2)
    
    if ("Err" %in% colnames(data)) {

      toperrline <- data.frame(x=c(data$Pos-0.49,data$Pos+0.49),y=c(data$Y+data$Err,data$Y+data$Err))
      toperrline <- toperrline[order(toperrline$x),]
      boterrline <- data.frame(x=c(data$Pos-0.49,data$Pos+0.49),y=c(data$Y-data$Err,data$Y-data$Err))
      boterrline <- boterrline[order(boterrline$x),]
      boterrline$y <- unlist(lapply(boterrline$y,function(y) {max(0,y)}))
      polygon(c(toperrline$x,rev(boterrline$x)),c(toperrline$y,rev(boterrline$y)),col=grey(0.5,0.5),border=F)

    }


    grid(col=grey(0.5))
    points(1:nrow(data),rep(refy,nrow(data)),col=basecolors[data$Style],pch=ascii[data$Style],cex=0.6)
    lines(plotline$x,plotline$y)
  }



dev.off()