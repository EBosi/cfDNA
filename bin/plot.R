library("ggplot2")
source("https://bioconductor.org/biocLite.R")
biocLite(c("ggbio"))
library("ggbio")
library(GenomicRanges)
library(grid)
library(gridExtra)

slidingWindow = function(i,v,wlength=3){
  # Function to compute a sliding window of a vector
  # nb wlength should be an odd number
  windowTailLength = (wlength - 1)/2
  left = max(1,i-windowTailLength)
  right = min(length(v),i+windowTailLength)
  return(mean(v[left:right]))
}

loadKaryoData = function(chr="all"){
  # Function to instanciate karyogram data as a data.frame
  # needed for a grand-scheme plot
  library("ggbio")
  library(GenomicRanges)
  if(chr=="all"){chr=paste("chr",c(1:21,"X","Y"),sep="")}
  data(hg19IdeogramCyto, package = "biovizBase")
  hg19 <- keepSeqlevels(hg19IdeogramCyto,chr,pruning.mode = "coarse")
  return(hg19)
}

setwd("~/hdd/DNA_circolante")
df = read.table("bon-ali.prova.tsv",header=T)
v = df$bon.alice
df$smoothed_v = sapply((1:length(v)), {function (x) slidingWindow(x,v)})

p1 <- ggplot(df,aes(x=1:length(v),bon.alice)) +
  geom_line()


gr <- GRanges(c("chr1", "chr1", "chr2", "chr3"), IRanges(1:4, width=3))
seqlevels(gr)
## Keep only 'chr1'
chr1 <- keepSeqlevels(gr, "chr1",pruning.mode = "coarse")
## Drop 'chr1'. Both 'chr2' and 'chr3' are kept.
chr2 <- dropSeqlevels(gr, "chr1")


data(hg19IdeogramCyto, package = "biovizBase")
hg19 <- keepSeqlevels(hg19IdeogramCyto,paste("chr",c(1:21,"X","Y"),sep=""),pruning.mode = "coarse")
p2 <- ggplot(hg19) + layout_karyogram(cytoband = TRUE)
  + theme_clear()


p1 = qplot(1:10, rnorm(10))
p2 = qplot(1:10, rnorm(10))
grid.arrange(p1, p2, ncol = 1)



p1 + p2
hg19 <- keepSeqlevels(hg19IdeogramCyto, paste0("chr", c(1:22, "X", "Y")))
p <- ggplot(hg19) + layout_karyogram(cytoband = TRUE)
p

