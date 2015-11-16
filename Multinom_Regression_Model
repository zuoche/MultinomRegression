allele_subg <- function (dat, sSG, eSG, sSnp, eSnp, grnam, ...)
{
  snp<-c()
  SNP<-c()
  rowname<-c()
  pval<-c()
    
  da<-dat[,c(sSG:eSG,sSnp:eSnp)]
  lend<-eSG-sSG+1
  lend2<-eSnp-sSnp+1
  
  for (i in (lend+1):(lend2+lend)) {
    ssnp<-names(da)[i]
    formul<-reformulate(termlabels = grnam, response = ssnp)
    analysis<-summary(multinom(formul, data=da))
    snp<-append(snp,ssnp,after=length(snp))
    
    pval<-rbind(pval, pt(abs(analysis$coefficients/analysis$standard.errors), df=nrow(da)-16, lower=F))
    
    
  }
  
  for (i in 1:length(snp)) SNP[i] <- substr(snp[i], 1, nchar(as.character(snp[i]))-2)
  for (i in 1:length(SNP)) rowname <- append(rowname, rep(SNP[i],2))
  total<-cbind(rowname,pval)
  
  return(total)
  
  
  
  
}

test<-allele_subg(ccp2,3,14,21,1187,"KL_biogrp1+KL_biogrp2+KL_biogrp3+KL_biogrp4+KL_biogrp5+KL_biogrp6+KL_biogrp7+KL_biogrp8+KL_biogrp9+KL_biogrp10+KL_biogrp11+KL_biogrp12")
