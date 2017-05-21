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

  # -------------------------------------------------------------------------- #

  # retrieve phase.columnname and indicator of s2-terrestrial-grid:
  phase.col<- phase_id[["phase.col"]]
  s2.ind<- phase_id[["terrgrid.id"]]


  # -------------------------------------------------------------------------- #

  # fit models with lm()-function, based on terrestrial sample s2:
  model.object <- lm(formula, data=data[data[,phase.col] == s2.ind,], x=TRUE, y=TRUE)

  # -------------------------------------------------------------------------- #

  # derive design-matrices restricted to terrestrial sample s2:
  design_matrix.s2 <- model.object$x

  # derive design-matrices restricted to reduced aux.sample s1:
  design_matrix.s1 <- design_matrix.s1_return(formula=formula, data=data)
  Y <- model.object$y


  # -------------------------------------------------------------------------- #
  # determine sample sizes:
  n1 <- nrow(design_matrix.s1)
  n2 <- nrow(design_matrix.s2)

  # -------------------------------------------------------------------------- #

  # calculate the Z_bars...:
  Z_bar <- apply(design_matrix.s1, 2, mean)

  # correct Z_bar boundary weights if defined (replace by new function...)
  if(!is.na(boundary_weights)){
    formula_plus_boundary.weight <- as.formula(paste(paste(all.vars(formula)[1],"~"),
                                                     paste(c(attr(terms(model.object),  "term.labels"), boundary_weights), collapse= "+")))
    w <- data.frame(design_matrix.s1_return(formula=formula_plus_boundary.weight, data=data))[,boundary_weights]
    Z_bar <- apply(design_matrix.s1, 2, weighted.mean, w=w)
  }

  # calculate covariance matrix of Z_bar:
  design_matrix.s1_centered <- t(apply(design_matrix.s1, 1, function(x){x-Z_bar}))
  cov_Z_bar <- t(design_matrix.s1_centered) %*% design_matrix.s1_centered / (n1*(n1-1))

  # -------------------------------------------------------------------------- #

  # calculate the g-weights:
  As2inv <- solve( t(design_matrix.s2) %*% design_matrix.s2/n2 )
  g <- Z_bar %*% As2inv %*% t(design_matrix.s2) # g-weights

  # get regression coefficients:
  beta_s2<-  model.object$coefficients
  # beta_s2<- t(As2inv %*% t(Y %*% design_matrix.s2 / n2)) same as extracted by lm()

  # get residuals:
  resid<- model.object$residuals


  # asymp. consistent design-based covariance-matrices  of regression coefficients:
  cov_beta_s2<- As2inv %*% ((t(resid*design_matrix.s2) %*% (resid*design_matrix.s2)) / n2^2) %*% As2inv



  # -------------------------------------------------------------------------- #

  # calculate estimates:

  # estimate <- mean(g*Y) # 1) using g-weights
  estimate<- Z_bar %*% beta_s2 # 2) equivalent to 1)


  # External Variance: This is the better external formula available only when Y_hat and R_hat are orthogonal:
  ext_variance <- (1/n1)*var(design_matrix.s1%*%coef(model.object)) + (1/n2)*var(resid)

  # g.weight Variance:

  # g_variance <- (1/n2^2)*sum((g^2)*(resid^2)) + (1/n1)*var(design_matrix.s1%*%coef(model.object)) # see .... ?
  g_variance <- (t(Z_bar) %*% cov_beta_s2 %*% Z_bar) + (t(beta_s2) %*% cov_Z_bar %*% beta_s2) # see ... ? (according Y_G_psynth in Report 1, page 13 [28])


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
                cov_coef=cov_beta_s2,
                Z_bar_1=Z_bar,
                cov_Z_bar_1G=cov_Z_bar,
                Rc_x_hat=resid,
                mean_Rc_x_hat=mean(resid),
                warn.messages=warn.messages)

  class(result)<- c("global", "twophase")

  return(result)

}


