setwd("~/hdd/DNA_circolante/references/")
source("../loess_norm.r")
source("../normalizeR.r")

table = read.csv("campioni_alt.csv")
sample = read.csv("../195-17.binned.50000.counts.reads.csv")

makeControlGenderList = function(){
  table = read.csv("/hdd/DNA_circolante/references/campioni_alt.csv")
  cNames = colnames(table)[grep("GC.10M",colnames(table),fixed=T)]
  males = sub(".COUNT.GC.10M","",cNames[22:42])
  females = sub(".COUNT.GC.10M","",cNames[1:21])
  return(list(males=males,females=females))
  }

makeOldControlTable = function(){
  wd = "/hdd/DNA_circolante/references/controlli"
  core_df = read.csv(paste(wd,"110-17.binned.50000.counts.reads.csv",sep="/"))[,1:6]
  for(f in system(paste("ls",wd,sep=" "),intern=T)){
    my_df = read.csv(paste(wd,f,sep="/"))[,7:12]
    core_df = cbind(core_df,my_df)
  }
  return(core_df)
}


makeMedianVector = function(table,alt_col=F){
  allGC = colnames(table)[grep("GC.10M",colnames(table),fixed=T)]
  autosomesGC = colnames(table)[grep("GC.autosomes.10M",colnames(table),fixed=T)]
  if(alt_col){
    allGC = colnames(table)[grep("GC.1GR0M",colnames(table),fixed=T)]
    autosomesGC = colnames(table)[grep("GC.AUTO.GR10M",colnames(table),fixed=T)]
  }
  allMedian = apply(table[table$CHR %in% c("X","Y"),allGC],1,median)
  autosomeMedian = apply(table[!table$CHR %in% c("X","Y"),autosomesGC],1,median)
  countMedian = c(autosomeMedian,allMedian)
  return(countMedian)
}

makeSampleValues = function(df){
  extractor = sample$CHR %in% c("X","Y")
  autosomes = normalize10M(loessGCNormalize_alt(sample,autosomal = T))[!extractor]
  all = normalize10M(loessGCNormalize_alt(sample,autosomal = F))[extractor]
  return(c(autosomes,all))
}

extractor = sample$CHR %in% c("X","Y")
sampleOldValues = c(sample$X195.17.GC.AUTO.GR10M[!extractor],sample$X195.17.GC.GR10M[extractor])
sampleNewValues = c(makeSampleValues(sample))

controlMedian = makeMedianVector(table)

finalValuesOld = sampleOldValues/controlMedian
finalValuesNew = sampleNewValues/controlMedian

bonAli = read.table("/hdd/DNA_circolante/bon-ali.prova.tsv",header=T)[,1]

karyoPlotter("campioni_alt.csv","1",finalValuesOld[sample$CHR=="1"])



df = read.csv("~/hdd/DNA_circolante/references/counts-controllli.csv")

controlMedian_old = makeMedianVector(old_tbl,alt_col =T)

finalValuesOld = sampleOldValues/controlMedian_old

head(df)
colnames(df)[7:158]

plot(finalValuesOld,finalValuesNew)
finalValuesNew>1000

plot(finalValuesNew[which(!finalValuesNew>10)],finalValuesOld[which(!finalValuesNew>10)])
