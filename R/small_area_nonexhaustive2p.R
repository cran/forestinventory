# This function helps the "two phase" estimator function when the first phase sample size is finite.
# This corresponds double sampling in small area estimation
# The external and g-weight variances with extended indicator model are calculated

# INPUTS:
#   formula = model formula
#   data = data.frame
#   boundary_weights = char string indicating weight variable representing the proportion of the forest area in the s1 interpretation area=
#   small_area = list indicating variable which plots occur in target small area, the value of the variable, whether to use extended model (i.e. unbiased) (small area must equal 1)
#
# OUTPUTS:
#   estimate = point estimate
#   ext_variance = external variance
#   g_variance = g-weight variance
#   n1 = sample size of all first phase plots
#   n2 = sample size of all second phase plots
#   n1G = sample size of all first phase plots in the small area G
#   n2G = sample size of all second phase plots in the small area G
#   r.squared = R squared of model


small_area_nonexhaustive2p <- function(formula, data, phase_id, small_area, boundary_weights, thresh.n2g){
  
  # intialize warning-message vector (i.e. collector of warnings produced during function execution):
  w1<- NA; w2<- NA; w3<- NA; w4<- NA
  
  # model fit applied to entire forest area
  model.object <- lm(formula, data=data[data[,phase_id[["phase.col"]]]==phase_id[["terrgrid.id"]],], x=TRUE, y=TRUE) 
  # local density:
  Yx<-model.object$y
  # plot level R square:
  r.squared <- summary(model.object)$r.squared
  
  # design matrix for s1 is only in the small area
  design_matrix.s1G <- design_matrix.s1_return(formula=formula, data=data[data[,small_area[["sa.col"]]]==1,])
  # design matrix for s1 is in entire area (all areas in dataset):
  design_matrix.s1 <- design_matrix.s1_return(formula=formula, data=data)

  # determine sample sizes:
  n1 <- nrow(design_matrix.s1)
  n1G <- nrow(design_matrix.s1G)
  n2 <- nrow(model.object$x)
  n2G <- sum(!is.na(data[,all.vars(formula)[1]]) & data[,phase_id[["phase.col"]]]==phase_id[["terrgrid.id"]] & data[,small_area[["sa.col"]]]==1) 

  ## ------- calculate estimations ----------------------------------------- ##
  
  design_matrix.s2 <- model.object$x
  p <- ncol(design_matrix.s2)
  Rx_s2 <- model.object$residuals
  
  # in case we don`t have any s2-points in small area G:
  if(n2G==0){ 
    if(!small_area[["unbiased"]]){ 
      # if "unbiased" == FALSE, we have no indicator introduced and hence we can invert the design-matrix although n2G == 0:
      As2inv <- solve( t(design_matrix.s2) %*% design_matrix.s2/n2 )
      Rx_s2G <- NA
      Yx_s2G <- NA
    } else { 
      Rx_s2G <- NA
      Yx_s2G <- NA
      As2inv <- matrix(NA, nrow=p, ncol=p) # ... we can't invert the design-matrix (singularity)
      r.squared<- NA # so we also set r.squared to NA
      w1<- ifelse(small_area[["unbiased"]], 
                  warning(paste("Unbiased estimation for Small Area", small_area[["areas"]], "not possible, set to 'NA': Small Area",small_area[["areas"]] ,"does not contain any terrestrial data", sep = " "), call. = F), NA)
    }
  }
  
  # in case we have only one s2-points in small area G:
  if(n2G==1){ Yx_s2G <- data[data[,small_area[["sa.col"]]]==1, all.vars(formula)[1]] #residuals only in small area are needed in external variance
              Rx_s2G <- Yx_s2G - predict.lm(model.object, data[data[,small_area[["sa.col"]]]==1,])
              Rx_s2G <- Rx_s2G[!is.na(Rx_s2G)]
              As2inv <- solve( t(design_matrix.s2) %*% design_matrix.s2/n2 ) # p. 17 Mandallaz Small Area technical report
              w2<- warning(paste("Terrestrial sample size in small area",small_area[["areas"]],"is 1. Computation of external variance not possible", sep = " "), call. = F)
#                if(small_area[["unbiased"]]){
#                   w4<- warning(paste("Terrestrial sample size in small area",small_area[["areas"]],"is considerably small. Consider true synthetic estimator instead", sep = " "), call. = F)
#                }
  }
  
  # in case we have more than one s2-point in small area G:
  if(n2G >= 2){Yx_s2G <- data[data[,small_area[["sa.col"]]]==1, all.vars(formula)[1]] #residuals only in small area are needed in external variance
               Rx_s2G <- Yx_s2G - predict.lm(model.object, data[data[,small_area[["sa.col"]]]==1,])
               Rx_s2G <- Rx_s2G[!is.na(Rx_s2G)]
               As2inv <- solve( t(design_matrix.s2) %*% design_matrix.s2/n2 ) # p. 17 Mandallaz Small Area technical report
#                if(!small_area[["unbiased"]] & n2G > thresh.n2g){
#                  w3<- warning(paste("Terrestrial sample size in small area",small_area[["areas"]],"sufficient for unbiased estimation. Consider extended pseudosynthetic estimator instead.",sep = " "), call. = F)
#                  }
#                if(small_area[["unbiased"]] & n2G <= thresh.n2g){
#                  w4<- warning(paste("Terrestrial sample size in small area",small_area[["areas"]],"is considerably small. Consider true synthetic estimator instead", sep = " "), call. = F)
#                  }
  }

  #asymp. consistent design-based cov matrix of beta_s2
  cov_beta_s2 <-   As2inv %*% ((t(Rx_s2*design_matrix.s2) %*% (Rx_s2*design_matrix.s2))/n2^2) %*% As2inv # eq [43] p. 18 Mandallaz Small Area technical report

  # compute unweighted auxiliary variable mean vector
  Z_bar_1G <- apply(design_matrix.s1G, 2, mean)
  
  
  # correct Z_bar_1G boundary weights if defined
  if(!is.na(boundary_weights)){
      formula_plus_boundary.weight <- as.formula(paste(paste(all.vars(formula)[1],"~"), paste(c(attr(terms(model.object),  "term.labels"), boundary_weights), collapse= "+")))
      #return design matrix with boundary_weights to treat NAs but only save weight vector
      w <- data.frame(design_matrix.s1_return(formula=formula_plus_boundary.weight, data=data[data[,small_area[["sa.col"]]]==1,]))[,boundary_weights] 
      Z_bar_1G <- apply(design_matrix.s1G, 2, weighted.mean, w=w)
  } 
  
  
  #store model coeff with indicator model
  beta_s2<- (Yx %*% design_matrix.s2 / n2) %*% As2inv 
  
  # calculate covariance matrix of Z_bar_1G: 
  Z1G_centered <- t(apply(design_matrix.s1G, 1, function(Zx,...){Zx - Z_bar_1G}))
  cov_Z_bar_1G <- (t(Z1G_centered) %*% Z1G_centered)/(n1G*(n1G-1))# eq [52] p. 20 Mandallaz Small Area technical report
  
  # calculate estimations:
  estimate <- Z_bar_1G %*% t(beta_s2) # eq [48] p. 19 Mandallaz Small Area technical report (if "unbiases" == TRUE); eq [18] p. 9  (if "unbiases" == FALSE)
  
  g_variance <- (Z_bar_1G %*% cov_beta_s2 %*% Z_bar_1G) + (beta_s2 %*% cov_Z_bar_1G %*% t(beta_s2)) # if "unbiased" == TRUE, we apply eq. [51] p. 19 Mandallaz Small Area technical report
                                                                                                    # if "unbiased" == FALSE, we apply eq. [28] p. 13 Mandallaz Small Area technical report
  ext_variance <- (1/n1G)*var(Yx_s2G, na.rm=TRUE) + (1-(n2G/n1G))*(1/n2G)*var(Rx_s2G, na.rm=TRUE)  

 
   ## ------- create outputs ------------------------------------------------- ##
  
  # summarize sample size info:
  samplesizes<- data.frame(cbind (n1G, n2G, n1, n2))
  colnames(samplesizes)<- c("n1G", "n2G", "n1", "n2")
  rownames(samplesizes)<- "plots"
  
  estimation<- c(estimate=estimate, ext_variance=ext_variance, g_variance=g_variance, 
                 n1=n1, n2=n2, n1G=n1G, n2G=n2G, r.squared=r.squared)
  
  # est.type can be "unbiased" in case of using the indicator, and "biased" in case of the true synthetic estimator:
  est.type<- ifelse(small_area[["unbiased"]], "psynth extended", "psynth")
  
  # save warning-messages:
  warn.messages<- c(w1, w2, w3, w4)[which(!is.na(c(w1, w2, w3, w4)))]

  return(
    list(regmodel=formula,
         estimation=estimation, 
         samplesizes=samplesizes, 
         coef=beta_s2, 
         cov_coef=cov_beta_s2, 
         Z_bar_1G=Z_bar_1G, 
         cov_Z_bar_1G=cov_Z_bar_1G, 
         Rc_x_hat_G = Rx_s2G,
         Mx_s2G = 1, # just a trick for computing the psmall ...
         mean_Rc_x_hat_G=mean(Rx_s2G),
         Yx_s2G=Yx_s2G, 
         boundary_weights=boundary_weights,
         estimator=est.type,
         warn.messages=warn.messages)
  )
  
  
}
