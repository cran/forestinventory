

# small_area_looper executes sae-estimation for ervery small area existing in the dataset
# and returns result-matrix
# we control the critical minimal number of n2g-points here ... parameter "thresh.n2g"

small_area_looper<- function(formula, data, phase_id, cluster, small_area, boundary_weights, exhaustive, progressbar, psmall){

  # *** for the package authors: *** #
  # --> set critical threshold of n2g for unbiased estimations:
  thresh.n2g = 6


  # -----------------------------------------------------------------------------#
  # create enclosure:
  # -----------------------------------------------------------------------------#


  # This enclosure function treats each small area successively and converts the variable to 1 0 indicator.
  # If model formula needs to be extended it does that, otherwise ( if "unbiased" == FALSE), the indicator variable is just not used (as if not introduced at all).
  extended_model <- function(){

    data[[small_area[["sa.col"]]]]<- as.character(data[[small_area[["sa.col"]]]])
    data[data[[small_area[["sa.col"]]]] != small_area[["areas"]][i], small_area[["sa.col"]]] <- 0
    data[data[[small_area[["sa.col"]]]] == small_area[["areas"]][i], small_area[["sa.col"]]] <- 1
    data[[small_area[["sa.col"]]]] <- as.numeric(data[[small_area[["sa.col"]]]])
    small_area[["areas"]]<- small_area[["areas"]][i]

    # extending model-formula by indicator variable for small area G if unbiased =TRUE
    if(small_area[["unbiased"]]){formula <- as.formula(paste(paste(all.vars(formula)[1],"~"), paste(c(all.vars(formula)[-1], small_area[["sa.col"]]), collapse= "+")))}

    # if small_area_nonexhaustive2p is required, then this function is # # sourced and applied:
     if(is.na(cluster) & all(is.na(exhaustive))){
      # source("small_area_nonexhaustive2p.R")
      return(small_area_nonexhaustive2p(formula=formula, data=data, phase_id=phase_id, small_area=small_area,
                                        boundary_weights=boundary_weights, thresh.n2g))
     }

     # if small_area_nonexhaustive2p_cluster is required, then this function is # sourced and applied:
     if(!is.na(cluster) & all(is.na(exhaustive))){
       # source("small_area_nonexhaustive2p_cluster.R")
       return(small_area_nonexhaustive2p_cluster(formula=formula, data=data, phase_id=phase_id, cluster=cluster, small_area=small_area,
                                                 boundary_weights=boundary_weights, thresh.n2g, psmall))
     }

    # if small_area_exhaustive2p is required, then this function is # sourced and applied:
    if(is.na(cluster) & all(!is.na(exhaustive))){
      # source("small_area_exhaustive2p.R")
      return(small_area_exhaustive2p(formula=formula, data=data, phase_id=phase_id, small_area=small_area,
                                     exhaustive=exhaustive[which(rownames(exhaustive) == small_area[["areas"]]), ], thresh.n2g))
    }

    # if small_area_exhaustive2p_cluster is required, then this function is # sourced and applied:
    if(!is.na(cluster) & all(!is.na(exhaustive))){
        # source("small_area_exhaustive2p_cluster.R")
        return(small_area_exhaustive2p_cluster(formula=formula, data=data, phase_id=phase_id, cluster=cluster, small_area=small_area,
                                               boundary_weights=boundary_weights,
                                               exhaustive=exhaustive[which(rownames(exhaustive) == small_area[["areas"]]), ], thresh.n2g, psmall))
    }
     } # end of extended_model function



  # -----------------------------------------------------------------------------#
  # initialize result matrices and lists:
  # -----------------------------------------------------------------------------#

  # ... to store inputs used:
  inputs<- list()

  # ... to store the estimates:
  # sa.col <- small_area[["sa.col"]]
  sa.areas <- small_area[["areas"]]
  no.sa <- length(sa.areas)
  if (progressbar) {pb <- txtProgressBar(min = 0, max = no.sa, style = 3)}

  result.sa <- data.frame(area=rep(NA_character_, no.sa), estimate=NA_real_, ext_variance=NA_real_,
                          g_variance=NA_real_, n1=NA_integer_, n2=NA_integer_, n1G=NA_integer_,
                          n2G=NA_integer_, r.squared=NA_real_)
  result.sa[,"area"]<-as.character(result.sa[,"area"])

  # ... to store the sample sizes:
  samplesizes<- list()

  # ... to store the regression coefficients:
  no.col.designmatrix<- ncol(design_matrix.s1_return(formula=formula, data=data)) + ifelse(small_area[["unbiased"]], 1, 0)
  coef<- data.frame(matrix(data = NA, nrow = no.sa, ncol = 1 + no.col.designmatrix))
  colnames(coef)[1]<- "area"
  coef[,"area"]<-as.character(coef[,"area"])

  # ... to store variance_covariance matrix of regression coefficients:
  cov_coef<- list()

  # ... to store variance_covariance matrix of auxiliary variables on cluster level:
  cov_Z_bar_1G<- list()

  # ... to store mean of auxilary variables:
  Z_bar_1G<- data.frame(matrix(data = NA, nrow = no.sa, ncol = 1 + no.col.designmatrix))
  colnames(Z_bar_1G)[1]<- "area"
  Z_bar_1G[,"area"]<-as.character(Z_bar_1G[,"area"])

  # ... to store Residuals of the small area G:
  Rc_x_hat_G<- list()

  # ... to store local densities of the small area G:
  Yx_s2G<- list()

  # ... to store M(x) of the small area G:
  Mx_s2G<- list()

  # ... to store mean of Residuals of the small area G:
  mean_Rc_x_hat_G<- data.frame(area=rep(NA_character_, no.sa), mean_Rc_x_hat_G=NA_real_)
  mean_Rc_x_hat_G[,"area"]<-as.character(mean_Rc_x_hat_G[,"area"])

  # ... to store estimation method:
  estimator<- data.frame(area=rep(NA_character_, no.sa), method=NA_integer_)
  estimator[,"area"]<-as.character(estimator[,"area"])

  # ... to store if boundary adjustment was used:
  # boundadj<- NA

  # ... to store warnings:
  warn.messages<- list()



  # -----------------------------------------------------------------------------#
  # loop over small areas:                                                       #
  # -----------------------------------------------------------------------------#

  for (i in 1:no.sa){
    #run enclosure to temporarily overwrite smallarea column with indicator variable for loop index i_th small variable
    # (the model is extended if necessary in this function)
    result.temp<- extended_model()

    # store inputs used (run only once for first loop):
    if (i==1){
    inputs[["data"]]<- data
    inputs[["smallareas"]]<- small_area[["areas"]]
    inputs[["formula"]]<- result.temp$regmodel
    inputs[["boundary_weights"]]<- result.temp$boundary_weights
    inputs[["method"]]<- result.temp$estimator
    inputs[["cluster"]]<- !is.na(cluster)
    inputs[["exhaustive"]]<- ifelse(all(is.na(exhaustive)), FALSE, TRUE)
    }

    # we should not give the external variance for the synthetic estimations ... (even if its possible, its not consistent)
    if(!small_area[["unbiased"]] & !psmall) {result.temp$estimation[["ext_variance"]]<- NA}

    # store estimation results
    result.sa[i,names(result.sa)[-1]] <- result.temp$estimation
    result.sa[i,names(result.sa)[1]] <- sa.areas[i] #add name of small area

    # store samplesize of small area:
    samplesizes[[sa.areas[i]]]<- result.temp$samplesizes

    # store variance_covariance matrix of auxiliary variables on cluster level:
    cov_Z_bar_1G[[sa.areas[i]]]<- result.temp$cov_Z_bar_1G

    # store mean of auxilary variables in G:
    Z_bar_1G[i,names(Z_bar_1G)[-1]] <- result.temp$Z_bar_1G
    Z_bar_1G[i,names(Z_bar_1G)[1]] <- sa.areas[i] #add name of small area
    if (i==1) colnames(Z_bar_1G)[-1]<- names(result.temp$Z_bar_1G)

    # store regression coefficients:
    coef[i,names(coef)[-1]] <- result.temp$coef
    coef[i,names(coef)[1]] <- sa.areas[i] #add name of small area
    if (i==1) colnames(coef)[-1]<- names(result.temp$Z_bar_1G)

    # store variance_covariance matrix of regression coefficients:
    cov_coef[[sa.areas[i]]] <- result.temp$cov_coef

    # Residuals in G:
    Rc_x_hat_G[[sa.areas[i]]]<- result.temp$Rc_x_hat_G

    # Local densities in G:
    Yx_s2G[[sa.areas[i]]]<- result.temp$Yx_s2G

    # Mx_s2G in G:
    Mx_s2G[[sa.areas[i]]]<- result.temp$Mx_s2G

    # mean of Residuals in G:
    mean_Rc_x_hat_G[i,names(mean_Rc_x_hat_G)[-1]] <- result.temp$mean_Rc_x_hat_G
    mean_Rc_x_hat_G[i,names(mean_Rc_x_hat_G)[1]] <- sa.areas[i] #add name of small area

    # warnings messages:
    warn.messages[[sa.areas[i]]] <- result.temp$warn.messages

    # visualize progress:
    if(progressbar){setTxtProgressBar(pb, i)}

  } # end of loop over small areas

  # close progress:
  if(progressbar){close(pb)}

  # --------------------------------------------------------------------- #


  # -------------------------------#
  # return results:                #
  # -------------------------------#

  result<- list(input=inputs,
                estimation=result.sa,
                samplesizes=samplesizes,
                coefficients =coef,
                cov_coef=cov_coef,
                Z_bar_1G=Z_bar_1G,
                cov_Z_bar_1G=cov_Z_bar_1G,
                Rc_x_hat_G=Rc_x_hat_G,
                Yx_s2G=Yx_s2G, #I just added this
                Mx_s2G = Mx_s2G,
                mean_Rc_x_hat_G=mean_Rc_x_hat_G,
                warn.messages=warn.messages)

  class(result)<- c("smallarea", "twophase")

  return(result)

}

