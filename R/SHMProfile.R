ARGS <- c(
  "mutfile","character","file path of mutation file",
  "clonefile","character","file path of summary file",
  "refseqfile","character","reference sequence file",
  "outstub","character","file to write profile to"
  
)

OPTS <- c(
  "tstart","numeric",0,"Start of reference to include",
  "tend","numeric",0,"End of reference to include",
  "minsubs","numeric",0,"minimum number of substitutions for a clone to be included",
  "maxsubs","numeric",0,"maximum number of substitutions for a clone to be included",
  "rmdups","logical",TRUE,"remove clones marked as duplicates"
)

source_local <- function(fname){
  argv <- commandArgs(trailingOnly = FALSE)
  base_dir <- dirname(substring(argv[grep("--file=", argv)], 8))
  source(paste(base_dir, fname, sep="/"))
}

source_local("Rsub.R")
source_local("SHMHelper.R")

parseArgs("SHMProfile.R", ARGS, OPTS)


suppressPackageStartupMessages(library(Biostrings, quietly=TRUE))

refseq <- readDNAStringSet(refseqfile)

clones <- read.delim(clonefile,header=T,as.is=T)
clones <- clones[clones$Coords != "",]

if (!all(grepl("([[:digit:]]+-[[:digit:]]+)(,[[:digit:]]+-[[:digit:]]+)*",clones$Coords))) {
  stop ("Clone coordinates not in correct format")
}

muts <- read.delim(mutfile,header=F,as.is=T,col.names=c("Expt","Clone","Pos","Type","From","To","Size","End","Ins"))

subs <- muts[ muts$Type == "sub",]

#### Exit if substitutions exist that don't have a clone in clone data.frame
if (nrow(subs[! subs$Clone %in% clones$Clone, ]) > 0 ) {
  stop ("Mutations from clones not included in clonefile")
}

if (anyDuplicated(clones$Clone) > 0) {
  stop ("Duplicate Clone IDs in clonefile")
}

if (tstart == 0) {
  tstart <- 1
}
if (tend == 0) {
  tend <- nchar(refseq)
}

subs <- subs[subs$Pos <= tend,]
subs <- subs[subs$Pos >= tstart,]


subtable <- table(subs$Clone)
emptyclones <- rep(0,nrow(clones[! clones$Clone %in% names(subtable), ]))
names(emptyclones) <- clones[! clones$Clone %in% names(subtable), "Clone"]
subtable <- c(subtable,emptyclones)


if (rmdups) {
  clones <- clones[ clones$Bp > 0, ]
  subs <- subs[ subs$Clone %in% clones$Clone, ]
}


if (minsubs > 0) {
  subs <- subs[ subtable[subs$Clone] >= minsubs, ]
  clones <- clones[ subtable[clones$Clone] >= minsubs, ]
}

if (maxsubs > 0) {
  subs <- subs[ subtable[subs$Clone] <= maxsubs, ]
  clones <- clones[ subtable[clones$Clone] <= maxsubs, ]
}



profile <- calculateProfile(subs,clones,refseq,tstart,tend)



write.table(profile,paste(outstub,"_profile.txt",sep=""),quote=F,sep="\t",na="",row.names=F,col.names=T)

