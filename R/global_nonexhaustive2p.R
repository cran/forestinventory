# This function helps the "two phase" estimator function when the first phase sample size is finite.
# This corresponds double sampling in global estimation
# The external and g-weight variances are calulated

# INPUTS:
#   formula = model formula
#   data = data.frame
#   boundary_weights = vector of weights representing the proportion of the forest area in the s1 interpretation area (vector of 1s if NA)
#
# OUTPUTS:
#   estimate = point estimate
#   ext_variance = external variance
#   g_variance = g-weight variance
#   n1 = ...
#   n2 = sample size
#   r.squared = R squared of model

global_nonexhaustive2p <- function(formula, data, phase_id, boundary_weights,...){

  # retrieve phase.columnname and indicator of s2-terrestrial-grid:
  phase.col<- phase_id[["phase.col"]]
  s2.ind<- phase_id[["terrgrid.id"]]

  model.object <- lm(formula, data=data[data[,phase.col] == s2.ind,], x=TRUE, y=TRUE)
  design_matrix.s2 <- model.object$x
  design_matrix.s1 <- design_matrix.s1_return(formula=formula, data=data)
  Y <- model.object$y
  n1 <- nrow(design_matrix.s1)
  n2 <- nrow(design_matrix.s2)
  As2inv <- solve( t(design_matrix.s2) %*% design_matrix.s2/n2 )
  Z_bar <- apply(design_matrix.s1, 2, mean)

  # correct Z_bar boundary weights if defined
  if(!is.na(boundary_weights)){
    formula_plus_boundary.weight <- as.formula(paste(paste(all.vars(formula)[1],"~"),
                                                     paste(c(attr(terms(model.object),  "term.labels"), boundary_weights), collapse= "+")))
    #return design matrix with boundary_weights to treat NAs but only save weight vector
    w <- data.frame(design_matrix.s1_return(formula=formula_plus_boundary.weight, data=data))[,boundary_weights]
    Z_bar <- apply(design_matrix.s1, 2, weighted.mean, w=w)
  } #partial interpretation intersect forest weight adjustment

  g <- Z_bar %*% As2inv %*% t(design_matrix.s2) # g-weights
  resid <- model.object$residuals
  estimate <- mean(g*Y)
  ext_variance <- (1/n1)*var(design_matrix.s1%*%coef(model.object)) + (1/n2)*var(resid) #This is the better external formula available only when Y_hat and R_hat are orthogonal
  g_variance <- (1/n2^2)*sum((g^2)*(resid^2)) + (1/n1)*var(design_matrix.s1%*%coef(model.object))


  ## ------- create outputs ------------------------------------------------- ##

  # summarize sample size info:
  samplesizes<- data.frame(cbind (n1, n2))
  colnames(samplesizes)<- c("n1", "n2")
  rownames(samplesizes)<- "plots"

  estimation<- data.frame(estimate=estimate, ext_variance=ext_variance, g_variance=g_variance,
                 n1=samplesizes$n1, n2=samplesizes$n2, r.squared=summary(model.object)$r.squared)

  # ... to store inputs used:
  inputs<- list()
  inputs[["data"]]<- data
  inputs[["formula"]]<- formula
  inputs[["boundary_weights"]]<- boundary_weights
  inputs[["method"]]<- "non-exhaustive"
  inputs[["cluster"]]<- FALSE
  inputs[["exhaustive"]]<- FALSE


  # save warning-messages:
  warn.messages<- NA

  result<- list(input=inputs,
                estimation=estimation,
                samplesizes=samplesizes,
                coefficients=coef(model.object),
                cov_coef=NA,   #WARNING:  This is calculable.  We just don't because we used the g-weight representation
                Z_bar_1=Z_bar,
                cov_Z_bar_1G=NA,
                Rc_x_hat=resid,
                mean_Rc_x_hat=mean(resid),
                warn.messages=warn.messages)

  class(result)<- c("global", "twophase")

  return(result)

}


















