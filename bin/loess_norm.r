getLoessInp = function(df,autosomal=T,return.mask=F){
  # Fxn to extract relevant data for loess
  if(autosomal){mask = df[,6]==0 & df[,7] > 0 & !(df[,2] %in% c("X","Y"))}
  else{mask = df[,6]==0 & df[,7] > 0}
  df_out = df[mask,c(5,7)]
  if(return.mask){return(list(new_df=df_out,mask=mask))}
  return(df_out)
}

loessRegr = function(vars,span=0.2,return_average.count=F){
  # Fxn to compute Loess-based regression
  # var1 and var2 should be GC count and read count of masked bins
  # that is, bins with Ns or with no reads are excluded
  var1 = vars[,1]
  var2 = vars[,2]
  data = data.frame(x=var1,y=var2)
  data_sorted = data[order(data$x),]
  data_sorted$GC=sapply(data_sorted$x/50000,function(x){round(x,2)})
  if(return_average.count){return(mean(data_sorted$y))}

  # compute average values
  GC_bin_counts = as.numeric(names(table(data_sorted$GC)))
  GC_bins_means = sapply(GC_bin_counts,function(x){mean(data_sorted[data_sorted$GC == x,]$y)})
  average_GC = data.frame(GC=GC_bin_counts,averageReadCount=GC_bins_means)

  # compute regression with loess
  regress.fit = loess(averageReadCount ~ GC, data=average_GC, span=span)
  regress.ys = predict(regress.fit)
  average_GC$e = regress.ys
  return(average_GC)  
  }


loessGCNormalize = function(mydf,autosomal=T,span=0.2){
  # get df and do loess
  df = mydf
  regress_df = loessRegr(getLoessInp(df,autosomal),span)
  average_count = loessRegr(getLoessInp(df,autosomal),span,return_average.count = T)
  
  # define normalization fnxs
  g = function(x,df_=df){
    # return GC value for x
    return(round(df_$BIN.GC.CONTENT[x]/50000,2))
  }
  
  y = function(x,df_=df){
    # return read count value for x
    return(df_[x,7])
  }
  
  f = function(g,df=regress_df){
    # return average read count for GC value g 
    return(df[df$GC==g,2])
  }
  
  e = function(x,df=regress_df){
    # return the expected value of x
    return(df[df$GC == x,3])
  }
  
  e_alt = function(x,df_=regress_df){
    # return the average value of count over all bins
    return(average_count)
  }
  
  correctionFxn = function(x,df_=df,alt=F){
    # correct x-th element as described in 
    # y'(x) = y(x) - [f(x) - e(x)]
    if(df_[x,6] >0 || df_[x,7]==0){return(0)}
    if(y(x)==0 || g(x)==0){return(0.0)}
    y_x = y(x)
    f_x = f(g(x))
    if(alt){e_x=e_alt()}
    else{e_x = e(g(x))}
    corrected = y_x - (f_x - e_x)
    return(corrected)}

  corrected_ys = sapply(1:dim(df)[1],correctionFxn)
  return(corrected_ys)
}

####################################################
## normalize exactly as said in Alkan et al - 2009 #
####################################################

loessGCNormalize_alt = function(mydf,autosomal=T){
  # get df and do loess
  filtered_inp = getLoessInp(mydf,autosomal,return.mask = T)
  norm_inp = filtered_inp$new_df
  my_mask = filtered_inp$mask
  colnames(norm_inp) = c("x","y")

  getResponseLM = function(lm,xs){
    coeff = lm$coefficients
    return((xs*coeff[2]) + coeff[1]) }
  
  getC = function(df){
    expected = mean(df$y)
    linearFit = lm(y ~ x,df)
    fitted = getResponseLM(linearFit,df$x)
    return(fitted - expected) }

  expected = mean(norm_inp$y)
  linearFit = lm(y ~ round(x/50000,1),norm_inp)
  fitted = expected * !my_mask
  fitted[my_mask] = getResponseLM(linearFit,round(mydf[my_mask,5]/50000,1))
  C.x = fitted - expected
  
  out_vec = mydf[,7] - C.x
  out_vec[!my_mask] = 0
  return(out_vec)
}

# prova
df="/home/eb/hdd/DNA_circolante/195-17.binned.50000.counts.reads.csv"
table = read.csv(df)
loessGCNormalize_alt(table,autosomal = T)
dim(table)

head(table)


prova = function(df="/home/eb/hdd/DNA_circolante/195-17.binned.50000.counts.reads.csv"){
  table = read.csv(df)
  return(loessGCNormalize(table))
}

