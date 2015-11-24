##################Regression analysis#################
library(nnet)
pass
analysis<-summary(multinom(formul, data=da))
# multi-level independent variables (# of Ab in each group)
str(analysis)
pval<-rbind(pval, pt(abs(analysis$coefficients/analysis$standard.errors), df=nrow(da)-16, lower=F))

################Backward stepwise########################
#Whether a group remains in the final model
#Decrease in AIC
library(nnet)
mod <- multinom(formul, data=da, model=T)
library(MASS)
aic <- stepAIC(mod, direction=”backward”)
AICvalue <- matrix(nrow=1167, ncol=12)
names(AICvalue) <- c("KL_biogrp1","KL_biogrp2","KL_biogrp3","KL_biogrp4","KL_biogrp5","KL_biogrp6","KL_biogrp7","KL_biogrp8","KL_biogrp9","KL_biogrp10","KL_biogrp11","KL_biogrp12")
for (j in 1:12) ifelse(names(AICvalue)[j] %in% summary(aic)$coefnames, AICvalue[i,j] <- summary(aic)$AIC, AICvalue[i,j] <- AICvalue[i,j])
# if a group included in the model, return AIC; else NA
# compare AIC in fitted model with full model

#####################All-subsets regression###################
#Adjusted 𝑅^2 for each SNP
#Mallow’s Cp
library(leaps)
modall <- regsubsets(formul, data=, nbest=1, nvmax=12)    
#nbest - # of models for each subset size
#nvmax – max subset size to be included
plot(modall, scale="adjr2")
library(car)
subsets(modall, statistic="cp") 

################Multicollinearity############################
#Generalized VIF among 12 biological groups as predictors
#No big differences across the chromosome
#Collinearity in Group1&5 worth concern (VIF>2.5 => 𝑅^2>0.6), but not yet problematic (VIF<10.0)
library(car)
vif(lm())      #applicable for lm() or glm()

#################GRS###############################
#1167 SNPs in relation to RA risk
#Sample size drop: ~70% missing
#Marker selection: >50% missing
#ANOVA within each bio-group
#Paire-wise t-test across 12 groups
aovX <- aov(GRS ~ sapply(KL_biogrpX, as.character), data=grs) 
library(multcomp)
K <- diag(length(coef(aovX)))[-1,]
rownames(K) <- names(coef(aovX))[-1]
ci <- glht(aov, linfct=K)  
#ncol(K)=length of coef(aovX)
plot(ci)








