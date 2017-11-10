###############
# Modulo Plot #
###############

########
# Import

# kayoploteR
library(karyoploteR)
# kpDataBackground(kp, data.panel = 1)
# kpText(kp, chr="chrX", x=60e6, y=0.5, labels="data.panel=1", data.panel = 1)
# kpAxis(kp,ymin=0.5,ymax=1.5)
# kp$chromosome.lengths["chr12"]


########
# Fnxs

slidingWindow = function(i,v,wlength=3){
  # Function to compute a sliding window of a vector
  # nb wlength should be an odd number
  windowTailLength = (wlength - 1)/2
  left = max(1,i-windowTailLength)
  right = min(length(v),i+windowTailLength)
  return(mean(v[left:right]))
}

createSlidingWindowArray = function(v,wlength=3){
  # wrapper function of slidingWindow
  out = sapply((1:length(v)), {function (x) slidingWindow(x,v,wlength)})
  return(out)
}

getXsBins = function(bins,binSize,chrLength){
  # get coordinate values for simple plotting
  # chrLength can be obtained from kp object
  # e.g. kp$chromosome.lengths
  out = bins * binSize - (binSize/2)
  if(out[length(out)] > chrLength){out[length(out)] = chrLength}
  return(out)
}

getBinsRanges = function(df){
  # get coordinate values for simple plotting
  # chrLength can be obtained from kp object
  # e.g. kp$chromosome.lengths
  starts = df
  ends
  out = bins * binSize - (binSize/2)
  if(out[length(out)] > chrLength){out[length(out)] = chrLength}
  return(out)
}

########
# Inputs

myChr="12"
binSize=50000
table="195-17.binned.50000.counts.reads.csv"

#######
# Script

setwd("~/hdd/DNA_circolante")

karyoPlotter = function(table,myChr,normCounts,binSize=50000){
  
  # read DF
  df = read.table(table,header=T,sep=",")
  df_ = df[df$CHR == myChr,]
  
  # get chromosome info
  kp = plotKaryotype(genome = "hg19",chromosomes = paste("chr",myChr,sep=""), plot.type=1, main="Titolo plot?")
  myChrLength = kp$chromosome.lengths
  
  # get bins (Xs)
  bins = getXsBins(1:length(df_$CHR),50000,myChrLength)
  
  # get norm counts (Ys)
  # bonAli = read.table("bon-ali.prova.tsv",header=T)[,1] # NB This should be computed from df!!!
  # smoothed = createSlidingWindowArray(bonAli,wlength=19)
  smoothed = createSlidingWindowArray(normCounts,wlength=19)  
  
  # Plot
  kpAxis(kp,ymin=0.5,ymax=1.5)
  kpLines(kp, chr=paste("chr",myChr,sep=""), x=bins, y=normCounts -0.5)
  kpLines(kp, chr=paste("chr",myChr,sep=""), x=bins, y=smoothed - 0.5,col="#FF0000")
  kpAddBaseNumbers(kp)
}
