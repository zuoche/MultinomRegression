allele_subg <- function (dat, sSG, eSG, sSnp, eSnp, grnam, ...)
{
  snp<-c()
  SNP<-c()
  rsq<-c()
  adjr2<-c()
  cp<-c()
  bic<-c()
  
  da<-dat[,c(sSG:eSG,sSnp:eSnp)]
  lend<-eSG-sSG+1
  lend2<-eSnp-sSnp+1
  
  for (i in (lend+1):(lend2+lend)) {
    ssnp<-names(da)[i]
    formul<-reformulate(termlabels = grnam, response = ssnp)
    analysis<-summary(regsubsets(formul, data=da, nbest=1, nvmax=12))
    snp<-append(snp,ssnp,after=length(snp))
    
    rsq<-rbind(rsq,analysis$rsq)   #R square
  
    adjr2<-rbind(adjr2,analysis$adjr2)   #adjusted R square
    
    cp<-rbind(cp,analysis$cp)   #Mallow's Cp
    
    bic<-rbind(bic,analysis$bic)   #BIC
    
  }
  
  for (i in 1:length(snp)) SNP[i] <- substr(snp[i], 1, nchar(as.character(snp[i]))-2)
  total<-cbind(SNP, adjr2)  #choose method=c("rsq","adjr2","cp","bic")
  
  return(total)
  
}
    
test<-allele_subg(ccp2, 3, 14, 21, 1187, "KL_biogrp1+KL_biogrp2+KL_biogrp3+KL_biogrp4+KL_biogrp5+KL_biogrp6+KL_biogrp7+KL_biogrp8+KL_biogrp9+KL_biogrp10+KL_biogrp11+KL_biogrp12")



