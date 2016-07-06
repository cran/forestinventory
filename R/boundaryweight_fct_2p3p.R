

# help-function to calculate weighted means (non-cluster sampling):
boundaryweight_fct_2p3p<- function(formula, model.object, data_select_from, boundary_weights){

  # extend formula.s0 by term "boundary-weights" to also extract the weights:
  formula_plus_boundary.weight <- as.formula(paste(paste(all.vars(formula)[1],"~"), paste(c(attr(terms(model.object),"term.labels"), boundary_weights), collapse= "+")))

  # extract design_matrix:
  design_matrix<- design_matrix.s1_return(formula=formula_plus_boundary.weight, data=data_select_from)

  # calculate the weighted mean of auxiliary variables (Z_bar):
  Z1_bar_weighted<- apply(design_matrix, 2, weighted.mean, w=design_matrix[,boundary_weights])

  Z1_bar_weighted[- which(names(Z1_bar_weighted) %in% boundary_weights)]

}
