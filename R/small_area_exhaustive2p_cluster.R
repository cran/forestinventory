
# INPUTS:
#   formula = model formula
#   data = data.frame
#   boundary_weights = char string indicating weight variable representing the proportion of the forest area in the s1 interpretation area=
#   small_area = list indicating variable which plots occur in target small area, the value of the variable, whether to use extended model (i.e. unbiased) (small area must equal 1)
#   cluster
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
# 
# plot(Y_TERR ~ X_TERR, data=data, pch=19, cex=0.7, col="grey")
# points(y= data[data$ismallg2==1,"Y_TERR"], x= data[data$ismallg2==1,"X_TERR"], pch=1, cex=0.7)

small_area_exhaustive2p_cluster<- function(formula, data, phase_id, cluster, small_area, boundary_weights, exhaustive, thresh.n2g, psmall){ # !!! ugly but true: we might also need boundary_weights in this case ...
  
  # intialize warning-message vector (i.e. collector of warnings produced during function execution):
  w1<- NA; w2<- NA; w3<- NA; w4<- NA
  
  s2 <- data[data[[phase_id[["phase.col"]]]]%in%phase_id[["terrgrid.id"]],]
  model.object <- lm(formula=formula, data=s2, x=TRUE, y=TRUE)
  # data_ext = design-matrix + infos on PLOT LEVEL
  data_ext<- data.frame(cbind(model.object[["x"]], model.object[["y"]]))
  colnames(data_ext)[ncol(data_ext)] <- all.vars(formula)[1]
  colnames.designmatrix<- colnames(data.frame(model.object[["x"]])) # trick: here, data_ext contains the colnames of the design-matrix, including the indicator-column only if "unbiased" == T
  data_ext[,cluster] <- data[rownames(data_ext), cluster] 
  data_ext[,small_area[["sa.col"]]]<- data[rownames(data_ext), small_area[["sa.col"]]] 
  if(!is.na(boundary_weights)){
    data_ext[,boundary_weights]<- data[data[[phase_id[["phase.col"]]]]%in%phase_id[["terrgrid.id"]],boundary_weights]
  }
  
  
  # plot level R square
  r.squared <- summary(lm(formula, data, y=TRUE))$r.squared 
  
  # info which cluster is s2 and s1 -> "mean" function here returns the phase_id on cluster level
#   gridind_clust<- aggregate(as.numeric(data_ext[ ,phase_id[["phase.col"]] ]), list(cluster = data_ext[ ,cluster]), mean)  #PROBLEMATIC!!! for factor phase IDs
#   colnames(gridind_clust)[2]<- phase_id[["phase.col"]]
#   
  # M(x) for entire s2-sample:
  cluster_weights <- aggregate(data_ext[,all.vars(formula)[1]], list(cluster = data_ext[,cluster]), length) 
  colnames(cluster_weights)[2]<- "Mx"

  IcG <- aggregate(data_ext[,small_area[["sa.col"]]], list(data_ext[,cluster]), mean)
  #colnames(IcG)[2]<- "IcG"
  colnames(IcG)<- c(cluster,small_area[["sa.col"]])
  

  # merge cluster-info together:
  clust_info<- merge(cluster_weights, IcG, by=cluster)

  ## ------- calculate everything for samples s2:---------------------------- ##
  
  ### --> retrieve Zc(x) for ENTIRE s2-sample:
  cluster_means <- aggregate(data_ext[,-c(which(names(data_ext)==cluster) ,which(names(data_ext)==small_area[["sa.col"]]))], # has to be without indicator variable!
                             list(cluster = data_ext[,cluster]), mean) 
  

  # for weighted cluster-means:
  if(!is.na(boundary_weights)){
    # library(plyr) #* load as dependency ...
    # weight the auxiliary information, but not the response:
    clustmeans.temp<- ddply(data_ext[,-c(which(names(data_ext) == small_area[["sa.col"]]), which(names(data_ext) == all.vars(formula)[1]))], 
                            .(cluster), function(x) colwise(weighted.mean, w = x[[boundary_weights]]) (x))
    clustmeans.temp<- clustmeans.temp[, -which(names(clustmeans.temp) == boundary_weights)] # remove column with boundary-weights, we don't need them anymore
    # merge unweighted cluster-means of response:
    cluster_means<- merge(clustmeans.temp, cluster_means[,c(cluster, all.vars(formula)[1])], by=cluster)
  }
  
  
  data_clust_s2<- merge(cluster_means, clust_info, by="cluster") # attach clusterinfo incl. indicator variable
  rm(cluster_means)
  
  # get cluster local density Yc(x):
  Yc_x <- data_clust_s2[, all.vars(formula)[1] ]

  # get design_matrices for s1 and s2:
  design_matrix.s2<- as.matrix(data_clust_s2[ ,colnames.designmatrix ]) # Trick: indicator column in data_clust_sx only gets into the design-matrix.sx, if "unbiased" == T
  
  # get cluster-weights:
  M_x.s2<- data_clust_s2[["Mx"]]
  
  # determine sample sizes:
  n1_clusters<- Inf
  n2_clusters<- nrow(design_matrix.s2)
  n1_plots<- Inf
  n2_plots<- nrow(data_ext)
  
  
  ## ------- calculate everything for samples s2_G:-------------------------- ##
  
  # get design_matrices for s1 and s2 (has to be calculated individually since s1/s2-cluster may not lie entirely in G):
  data_ext_s2G<- data_ext[ data_ext[small_area[["sa.col"]]] == 1, ]
  
  cluster_means_s2G<- aggregate(data_ext_s2G[ ,-c(which(names(data_ext_s2G)==cluster),which(names(data_ext_s2G)==small_area[["sa.col"]]))], 
                                list(cluster = data_ext_s2G[,cluster]), mean) 
  
  # for weighted cluster-means:
  if(!is.na(boundary_weights)){
    # library(plyr) #* load as dependency ...
    # weight the auxiliary information, but not the response:
    clustmeans_s2G.temp<- ddply(data_ext_s2G[,-c(which(names(data_ext_s2G) == small_area[["sa.col"]]), which(names(data_ext_s2G) == all.vars(formula)[1]))], 
                                .(cluster), function(x) colwise(weighted.mean, w = x[[boundary_weights]]) (x))
    clustmeans_s2G.temp<- clustmeans_s2G.temp[, -which(names(clustmeans_s2G.temp) == boundary_weights)] # remove column with boundary-weights, we don't need them anymore
    # merge unweighted cluster-means of response:
    cluster_means_s2G<- merge(clustmeans_s2G.temp, cluster_means_s2G[,c(cluster, all.vars(formula)[1])], by=cluster)
  }
  
  # local density on cluster-level of G: 
  Yc_x_G<- cluster_means_s2G[ ,all.vars(formula)[1]]
  
  # M(x)_s2G for  s2G-sample:
  cluster_weights.s2G <- aggregate(data_ext_s2G[,all.vars(formula)[1]], list(cluster = data_ext_s2G[,cluster]), length) 
  colnames(cluster_weights.s2G)[2]<- "Mx_s2G"
  
  # re-calculate indicator variable on cluster level IcG --> should now be "1" for all observations
  IcG.G<- aggregate(data_ext_s2G[[ small_area[["sa.col"]] ]],list(data_ext_s2G[,cluster]), mean)
  colnames(IcG.G)<- c(cluster,small_area[["sa.col"]])
  
  data_clust_s2G<- merge(cluster_means_s2G, IcG.G, by= cluster) # attach clusterinfo incl. indicator variable
  
  # the cluster-weights have to be merged AFTER calculating data_clust_s2G, since the merge()-function in line 117
  # can change the row-order of cluster_means_s2G --> if this happens, cluster_weights.s2G do not match the right data_clust_s1G-rows:
  Mx_s2G<- merge(data_clust_s2G, cluster_weights.s2G, by = cluster)[["Mx_s2G"]]
  
  # get design_matrices for s2_G:
  design_matrix.s2G<- as.matrix(data_clust_s2G[ ,colnames.designmatrix ])

  # determine sample sizes:
  n1G_clusters<- Inf
  n2G_clusters<- nrow(design_matrix.s2G)
  n1G_plots<- Inf
  n2G_plots<- nrow(data_ext_s2G)

  
  ## ------- calculate estimations ----------------------------------------- ##

  if(n2G_clusters==0){ 
    if(!small_area[["unbiased"]]) { # is "unbiased" == FALSE" ? -> then we have no indicator introduced
          As2inv <- solve( (t(design_matrix.s2) %*% (M_x.s2*design_matrix.s2)) / n2_clusters ) # ... and we can invert the design-matrix
          design_matrix.s2G<- matrix(NA, nrow=1, ncol= ncol(design_matrix.s2))
          Yc_x_G<- NA
          Mx_s2G<- NA
    } else { # if we have the indicator introduced ...
          As2inv <- matrix(NA, nrow=ncol(design_matrix.s2), ncol= ncol(design_matrix.s2)) # since we cannot invert the design-matrix due to singularity
          # plot level R square
          r.squared <- NA
          design_matrix.s2G<- matrix(NA, nrow=1, ncol= ncol(design_matrix.s2))
          Yc_x_G<- NA
          Rc_x_hat_G<- NA
          Mx_s2G<- NA
          w1<- warning(paste("Unbiased estimation for Small Area", small_area[["areas"]], "not possible, set to 'NA': Small Area",small_area[["areas"]],"does not contain any terrestrial data"), call. = F)
    }
  }
    
  if(n2G_clusters!=0){ 
    As2inv <- solve( (t(design_matrix.s2) %*% (M_x.s2*design_matrix.s2)) / n2_clusters )
#     if(!small_area[["unbiased"]] & n2G_plots > thresh.n2g){
#       w2<- warning(paste("Terrestrial sample size in small area",small_area[["areas"]],"sufficient for unbiased estimation. Consider extended pseudosynthetic estimator instead."),call. = F)
#     }
    if(small_area[["unbiased"]]){
#       if(n2G_plots <= thresh.n2g){
#         w3<- warning(paste("Terrestrial sample size in small area",small_area[["areas"]],"is considerably small. Consider true synthetic estimation instead"), call. = F)
#       }
      if (any(IcG[[small_area[["sa.col"]]]] > 0 & IcG[[small_area[["sa.col"]]]] < 1)){ 
        if(!psmall){
        w4<- warning(
          paste("At least one terrestrial cluster not entirely included within the small area ",small_area[["areas"]],".\n",
                "Zero mean residual assumption for small area maybe violated.","\n", 
                "Check mean_Rc_x_hat_G and consider alternative estimator 'psmall'", sep = ""), call. = F)
        }
      }
    }
 }
  

  # regression coefficient theta: eq [63] p. 24 Mandallaz Small Area technical report
  theta_s2<- ((M_x.s2 * Yc_x) %*% design_matrix.s2 / n2_clusters) %*% As2inv
  
  # (mean) residuals in F:
  Yc_x_hat_s2<- design_matrix.s2 %*% t(theta_s2)
  Rc_x_hat<- Yc_x - Yc_x_hat_s2
  mean_Rc_x_hat<- weighted.mean(Rc_x_hat, w = M_x.s2)
  
  # (mean) residuals in G (only possible if we have s2-data in G!)
  Yc_x_hat_s2_G<- design_matrix.s2G %*% t(theta_s2)
  Rc_x_hat_G<- Yc_x_G - Yc_x_hat_s2_G
  mean_Rc_x_hat_G<- weighted.mean(Rc_x_hat_G, w = Mx_s2G)

  # cov.matrix of theta_s2: eq [64] p. 24 Mandallaz Small Area technical report
  MR_square <- as.vector((M_x.s2^2)*(Rc_x_hat^2))
  middle_term<- ((t(design_matrix.s2) %*%  (MR_square * design_matrix.s2)) / n2_clusters^2)
  cov_theta_s2<- As2inv %*% middle_term %*% As2inv   
  
  #define unweighted auxiliary variable mean vector (according to design matrix!!! Could be tricky for categorical variables)
  if(small_area[["unbiased"]]) {Z_bar_c1G <- c(unlist(exhaustive),1)} else {Z_bar_c1G<- unlist(exhaustive)}
  names(Z_bar_c1G)<- colnames(design_matrix.s2)
  
  # estimates:
  estimate<- Z_bar_c1G %*% t(theta_s2)
  g_variance<- (t(Z_bar_c1G) %*% cov_theta_s2 %*% Z_bar_c1G)
  ext_variance<- (1/n2G_clusters) * (1/(n2G_clusters - 1)) * sum( ((Mx_s2G / mean(Mx_s2G))^2) * ((Rc_x_hat_G - mean_Rc_x_hat_G)^2) ) 

  
  
  ## ------- create outputs ------------------------------------------------- ##
 
  # summarize sample size info:
  samplesizes<- data.frame( rbind(cbind (n1G_clusters, n2G_clusters, n1_clusters, n2_clusters),
                                  cbind (n1G_plots, n2G_plots, n1_plots, n2_plots)))
  colnames(samplesizes)<- c("n1G", "n2G", "n1", "n2")
  rownames(samplesizes)<- c("clusters", "plots")
  
  estimation<- c(estimate=estimate, ext_variance=ext_variance, g_variance=g_variance, 
                 n1=n1_clusters, n2=n2_clusters, n1G=n1G_clusters, n2G=n2G_clusters, 
                 r.squared=r.squared)
  
  # est.type can be "unbiased" in case of using the indicator, and "biased" in case of the true synthetic estimator:
  est.type<- ifelse(small_area[["unbiased"]], "synth extended", "synth")
  
  # save warning-messages:
  warn.messages<- c(w1, w2, w3, w4)[which(!is.na(c(w1, w2, w3, w4)))]
  
  return(
    list(regmodel=formula,
         estimation=estimation, 
         samplesizes=samplesizes, 
         coef=theta_s2, 
         cov_coef=cov_theta_s2, 
         Z_bar_1G=Z_bar_c1G, 
         cov_Z_bar_1G=NA, # we dont have this term in the exhaustive case
         Rc_x_hat_G=Rc_x_hat_G, 
         Mx_s2G = Mx_s2G,
         mean_Rc_x_hat_G=mean_Rc_x_hat_G,
         boundary_weights=boundary_weights,
         estimator=est.type,
         warn.messages=warn.messages)
  )
  

} # end of "small_area_exhaustive2p_cluster"

