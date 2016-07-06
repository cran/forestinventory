

# This function "merges" the reduced- and full-model coefficients of a
# three-phase global estimation for a more compact output in the
# summary()
#

# ---------------------------------------------------------------------------- #

convert_coefs_table.global3p<- function(global_obj){

# place reduced-model coefficients beside full-model counterparts:

coefs<- global_obj$coefficients

coef.merge<- data.frame(matrix(data = NA, nrow = 2, ncol = length(coefs$beta)))
colnames(coef.merge)<- names(coefs$beta)
row.names(coef.merge)[1]<- ""
row.names(coef.merge)[2]<- "*"

# paste in full-model coefficients:
coef.merge[1,]<- coefs$beta

# paste in reduced-model coefficients:
ind.red.coefs<- which(colnames(coef.merge) %in% names(coefs$alpha))
coef.merge[2, ind.red.coefs]<- coefs$alpha

coef.merge

}

# ---------------------------------------------------------------------------- #

convert_coefs_table.smallarea3p<- function(sae_obj){

  coefs<- sae_obj$coefficients

  coef.merge<- data.frame(matrix(data = NA, nrow = 2*nrow(coefs$alpha), ncol = ncol(coefs$beta)-1))
  colnames(coef.merge)<- colnames(coefs$beta)[-1]
  row.names(coef.merge)[seq(1,nrow(coef.merge),by=2)]<- coefs$alpha$area
  row.names(coef.merge)[seq(2,nrow(coef.merge),by=2)]<- paste(coefs$alpha$area,"*", sep="")

  # paste in full-model coefficients:
  ind.beta<- which(row.names(coef.merge) %in% coefs$beta$area)
  coefs$beta<- coefs$beta[-(colnames(coefs$beta) %in% "area")]
  for(i in 1:length(ind.beta)){
    coef.merge[ind.beta[i], ]<- coefs$beta[i,]
  }

  # paste in reduced-model coefficients:
  ind.alpha<- which(row.names(coef.merge) %in% paste(coefs$alpha$area,"*", sep=""))
  coefs$alpha<- coefs$alpha[-(colnames(coefs$alpha) %in% "area")]
  ind.red.coefs<- which(colnames(coef.merge) %in% colnames(coefs$alpha))

  for(i in 1:length(ind.alpha)){
    coef.merge[ind.alpha[i], ind.red.coefs]<- coefs$alpha[i,]
  }

  coef.merge

}











