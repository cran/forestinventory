#Helper function that overrides response variable with non-NA (e.g. 1) so that the returned design.matrix for all of s1 is returned

#returns design matrix for s1 with complete.obs for auxiliary variables
design_matrix.s1_return <- function(formula, data){
  data[,all.vars(formula)[1]] <- 1 #change response values non-NA to trick lm function
  model.matrix(formula, data=data)}