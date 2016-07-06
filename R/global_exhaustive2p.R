# This function helps the "two phase" estimator function when the first phase sample size goes to infinity.
# This corresponds to the situation where wall2wall info is known for post-stratification
# The external and g-weight variances are calulated

# INPUTS:
#   formula = model formula
#   data = data.frame
#   exhaustive = the vector of auxiliary variable means known exhaustively (should ideally be weighted mean adjusted for proportion of pixels in interpretation are in the forest)
#
# OUTPUTS:
#   estimate = point estimate
#   ext_variance = external variance
#   g_variance = g-weight variance
#   n1 = Inf
#   n2 = sample size
#   r.squared = R squared of model

global_exhaustive2p <- function(formula, data, exhaustive){


  model.object <- lm(formula, data=data, x=TRUE, y=TRUE)
  design_matrix.s2 <- model.object$x
  n1 <- Inf
  n2 <- nrow(design_matrix.s2)
  As2inv <- solve( t(design_matrix.s2)%*%design_matrix.s2/n2 )
  g <- exhaustive %*% As2inv %*% t(design_matrix.s2)
  resid <- model.object$residuals
  estimate <- mean(g*model.object$y)
  ext_variance <- var(resid)/n2
  g_variance <- (1/n2^2)*sum((g^2)*(resid^2))





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
  inputs[["boundary_weights"]]<- NA
  inputs[["method"]]<- "exhaustive"
  inputs[["cluster"]]<- FALSE
  inputs[["exhaustive"]]<- TRUE


  # save warning-messages:
  warn.messages<- NA

  result<- list(input=inputs,
                estimation=estimation,
                samplesizes=samplesizes,
                coefficients=coef(model.object),
                cov_coef=NA,
                Z_bar_1=NA,
                cov_Z_bar_1G=NA,
                Rc_x_hat=resid,
                mean_Rc_x_hat=mean(resid),
                warn.messages=warn.messages)

  class(result)<- c("global", "twophase")

  return(result)

}
