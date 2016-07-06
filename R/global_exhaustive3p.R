# This function helps the "three phase" estimator function when the first phase sample size goes to infinity.
# --> also refered to as "generalized two-phase estimator with partially exhaustive information"
# The external and g-weight variances are calulated

# INPUTS:
#   formula = model formula
#   data = data.frame
#   exhaustive = the vector of auxiliary variable means known exhaustively for the reduced model
#                 (should ideally be weighted-mean-adjusted for proportion of pixels in interpretation are in the forest)
#
# OUTPUTS:
#   estimate = point estimate
#   ext_variance = external variance
#   g_variance = g-weight variance
#   n0 = Inf
#   n1 = sample size s1
#   n2 = sample size s2
#   r.squared = R squared of model

global_exhaustive3p <- function(formula.s0, formula.s1, data, phase_id, boundary_weights, exhaustive,...){

  # retrieve phase.columnname and indicator of s1 grid id and terrestrial-grid id:
  phase.col<- phase_id[["phase.col"]]
  s1.ind<- phase_id[["s1.id"]]        # s1.id identifies the small sample of the auxiliary vars (s1-sample, refers to formula.s1, ie. the full-model)
  s2.ind<- phase_id[["terrgrid.id"]]  # identifies the terrestrial sample (s2-sample)

  # fit models with lm()-function, based on terrestrial sample s2:
  model.object.reduced<- lm(formula.s0, data=data[data[,phase.col] == s2.ind,], x=TRUE, y=TRUE)     # "reduced" model
  model.object.full<-    lm(formula.s1, data=data[data[,phase.col] == s2.ind,], x=TRUE, y=TRUE)     # "full" model


  # derive design-matrices restricted to terrestrial sample s2:
  design_matrix_1.s2<- model.object.reduced$x  #"reduced" design matrix defined at terrestrial points s2
  design_matrix.s2<-   model.object.full$x      #"full" design matrix defined at terrestrial points s2


  # derive design-matrices restricted to reduced aux.sample s1:
  design_matrix_1.s1<- design_matrix.s1_return(formula=formula.s0, data=data[data[,phase.col] %in% c(s1.ind, s2.ind),]) #"reduced" design matrix defined at s1
  design_matrix.s1<-   design_matrix.s1_return(formula=formula.s1, data=data[data[,phase.col] %in% c(s1.ind, s2.ind),]) #"full" design matrix defined at s1
  # design_matrix1.s1 <- design_matrix.s1_return(formula=formula.s1, data=data) #"reduced" design matrix defined at s1

  # get sample-sizes:
  n0 <- Inf
  n1 <- nrow(design_matrix.s1)
  n2 <- nrow(design_matrix.s2)

  # calculate the Z_bars...:
  Z1_bar_s1<- apply(design_matrix_1.s1, 2, mean) # that's Z_bar(1)_1 in Daniels Report...
  Z_bar_s1<-  apply(design_matrix.s1, 2, mean)   # that's Z_bar_1 in Daniels Report...


  # calculate the Z_bars as weighted auxiliary means if boundary_weights are used:
  if(!is.na(boundary_weights)){

    Z1_bar_s1<- boundaryweight_fct_2p3p(formula=formula.s0,
                                        model.object=model.object.reduced,
                                        data_select_from=data[data[,phase.col] %in% c(s1.ind, s2.ind),],
                                        boundary_weights)

    Z_bar_s1<- boundaryweight_fct_2p3p(formula=formula.s1,
                                       model.object=model.object.full,
                                       data_select_from=data[data[,phase.col] %in% c(s1.ind, s2.ind),],
                                       boundary_weights)
  }


  # calculate the g-weights:
  A1_s1_inv<- solve( (t(design_matrix_1.s1)%*%design_matrix_1.s1) / n1 )
  A_s2_inv<-  solve( (t(design_matrix.s2)%*%design_matrix.s2) / n2 )
  g1_1<- exhaustive %*% A1_s1_inv %*% t(design_matrix_1.s2)   # that's g(1)_1(x) in the report, only needed for s2-sample
  g2<-  Z_bar_s1 %*% A_s2_inv %*% t(design_matrix.s2)         # that's g2(x) in the report


  # get regression coefficients:
  alpha<- coef(model.object.reduced)
  beta<-  coef(model.object.full)

  # get residuals:
  resid_reduced<- model.object.reduced$residuals # that's R_1 in the report
  resid_full<-    model.object.full$residuals    # that's R in the report


  # calculate estimates:
  estimate <- ((exhaustive - Z1_bar_s1) %*% alpha) + (Z_bar_s1 %*% beta)
  # ext_variance <- (var(resid_reduced)/n1) + ((1-(n2/n1)) * var(resid_full)/n2)
  ext_variance<- ((1/(n1*n2))*sum(resid_reduced^2)) + ((1-(n2/n1)) * (1/n2^2) * sum(resid_full^2))
  g_variance<- (1/(n1*n2))*sum((g1_1^2)*(resid_reduced^2)) + (1-(n2/n1)) * (1/n2^2) * sum((g2^2)*(resid_full^2))


  ## ------- create outputs ------------------------------------------------- ##

  # summarize sample size info:
  samplesizes<- data.frame(cbind (n0, n1, n2))
  colnames(samplesizes)<- c("n0", "n1", "n2")
  rownames(samplesizes)<- "plots"

  estimation<- data.frame(estimate=estimate, ext_variance=ext_variance, g_variance=g_variance,
                          n0=samplesizes$n0, n1=samplesizes$n1, n2=samplesizes$n2, r.squared_reduced=summary(model.object.reduced)$r.squared,
                          r.squared_full=summary(model.object.full)$r.squared)

  # ... to store inputs used:
  inputs<- list()
  inputs[["data"]]<- data
  inputs[["formula.s0"]]<- formula.s0
  inputs[["formula.s1"]]<- formula.s1
  inputs[["boundary_weights"]]<- boundary_weights
  inputs[["exhaustive"]]<- exhaustive
  inputs[["method"]]<- "exhaustive"
  inputs[["cluster"]]<- FALSE
  inputs[["exhaustive"]]<- TRUE

  # save warning-messages:
  warn.messages<- NA

  result<- list(input=inputs,
                estimation=estimation,
                samplesizes=samplesizes,
                coefficients=list(alpha=alpha, beta=beta),
                g_reduced=g1_1,
                g_full=g2,
                Z1_bar_s1=Z1_bar_s1,
                Z_bar_s1=Z_bar_s1,
                resid_reduced=resid_reduced,
                resid_full=resid_full,
                warn.messages=warn.messages)

  class(result)<- c("global", "threephase")

  result

}

