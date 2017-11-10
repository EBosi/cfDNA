# source("loess_norm.r")

normalize10M = function(v){
  return(v*(10000000/sum(v)))
}

# makeNormColumnNames(prefix){
#   suffixes=c(".GC",".GC.10M",".GC.autosomes",".GC.autosomes.10M")
#   sapply(suffixes,function(x){pas})
# }

computeReferenceCounts = function(csv_control="references/counts-controllli_.csv",sep=",",norm_fxn = "loessGCNormalize"){
  ref_df = read.csv(csv_control,sep=sep)
  referenceIDs = colnames(ref_df)[6:dim(ref_df)[2]]
  core_df = ref_df[,1:5]
  out_df = ref_df[,1:5]
  norm_FUN = match.fun(norm_fxn)
  for(id in referenceIDs){
    # print(id)
    gcNormCounts = norm_FUN(cbind(1,core_df,ref_df[,id]),autosomal=F)
    gcNormCounts10M = normalize10M(gcNormCounts)
    gcNormCounts.Autosome = norm_FUN(cbind(1,core_df,ref_df[,id]),autosomal=T)
    gcNormCounts10M.Autosome = normalize10M(gcNormCounts.Autosome)
    tmp_out = cbind(gcNormCounts,gcNormCounts10M,gcNormCounts.Autosome,gcNormCounts10M.Autosome)
    colnames(tmp_out) = paste(id,c(".GC",".GC.10M",".GC.autosomes",".GC.autosomes.10M"),sep="")
    out_df = cbind(out_df,tmp_out)
  }
  return(out_df)
}


# rapporto campione.gc.10m.tutti / 
# mediana popolazione per bin gc.10m.tutti

# normalizzato con la prima funzione
# asd = computeReferenceCounts("references/prova.gc_.csv",sep=";")
# write.csv(asd,"references/campioni.csv")

# normalizzato con la seconda funzione
# prova
 asd_ = computeReferenceCounts("references/prova.gc_.csv",sep=";",norm_fxn = "loessGCNormalize_alt")
 write.csv(asd_,"references/campioni_alt_prova.csv")
# # actual reference
 asd_ = computeReferenceCounts("../references/counts-controllli_.csv",norm_fxn = "loessGCNormalize_alt")
 write.csv(asd_,"../references/campioni_alt.csv")

