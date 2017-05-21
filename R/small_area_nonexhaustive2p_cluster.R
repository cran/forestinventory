
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

small_area_nonexhaustive2p_cluster<- function(formula, data, phase_id, cluster, small_area, boundary_weights, thresh.n2g, psmall){


  # intialize warning-message vector (i.e. collector of warnings produced during function execution):
  w1<- NA; w2<- NA; w3<- NA; w4<- NA

  ## ------- setting up the upper-level dataset for the small area:  ------- ##

  # data_ext = design-matrix + infos on PLOT LEVEL
  data_ext<- data.frame(design_matrix.s1_return(formula=formula, data=data))
  colnames.designmatrix<- colnames(data_ext) # trick: here, data_ext contains the colnames of the design-matrix, including the indicator-column only if "unbiased" == T
  data_ext[,all.vars(formula)[1]] <- data[,all.vars(formula)[1]]
  data_ext[,cluster] <- data[,cluster]
  data_ext[,small_area[["sa.col"]]]<- data[,small_area[["sa.col"]]]
  data_ext[,phase_id[["phase.col"]]]<- data[,phase_id[["phase.col"]]]
  if(!is.na(boundary_weights)){
  data_ext[,boundary_weights]<- data[,boundary_weights]
  }

  # plot level R square
  r.squared <- summary(lm(formula, data, y=TRUE))$r.squared

  # info which cluster is s2 and s1 -> "mean" function here returns the phase_id on cluster level
  gridind_clust<- aggregate(as.numeric(data[ ,phase_id[["phase.col"]] ]), list(cluster = data[ ,cluster]), mean)
  colnames(gridind_clust)[2]<- phase_id[["phase.col"]]

  # M(x) for entire s1-sample:
  cluster_weights <- aggregate(data_ext[,all.vars(formula)[1]], list(cluster = data_ext[,cluster]), length)
  colnames(cluster_weights)[2]<- "Mx"

  # re-calculate indicator variable on cluster level IcG: p. 24 in text, Mandallaz Small Area technical report.
  IcG<- aggregate(data_ext[[ small_area[["sa.col"]] ]],list(data_ext[,cluster]), mean)
  colnames(IcG)<- c(cluster,small_area[["sa.col"]])

    # merge cluster-info together:
  clust_info<- merge(merge(cluster_weights, gridind_clust, by = cluster), IcG, by=cluster)
  IcG.s2<- clust_info [clust_info[phase_id[["phase.col"]]]  == phase_id[["terrgrid.id"]], small_area[["sa.col"]] ]


  ## ------- calculate everything for samples s1 and s2: ------------------- ##

  ### --> retrieve Zc(x) for ENTIRE s1-sample (has to be WITHOUT indicator variable!):
  cluster_means <- aggregate(data_ext[,-c(which(names(data_ext)==cluster) ,which(names(data_ext)==small_area[["sa.col"]]),
                                          which(names(data_ext)==phase_id[["phase.col"]]))],
                             list(cluster = data_ext[,cluster]), mean)


  # for weighted cluster-means:
  if(!is.na(boundary_weights)){
    # library(plyr) #* load as dependency ...
    # weight the auxiliary information, but not the response:
    clustmeans.temp<- ddply(data_ext[,-c(which(names(data_ext) == small_area[["sa.col"]]), which(names(data_ext) == all.vars(formula)[1]),
                                         which(names(data_ext)==phase_id[["phase.col"]]))],
                            .(cluster), function(x) colwise(weighted.mean, w = x[[boundary_weights]]) (x))
    clustmeans.temp<- clustmeans.temp[, -which(names(clustmeans.temp) == boundary_weights)] # remove column with boundary-weights, we don't need them anymore
    # merge unweighted cluster-means of response:
    cluster_means<- merge(clustmeans.temp, cluster_means[,c(cluster, all.vars(formula)[1])], by=cluster)
  }


  data_clust_s1<- merge(cluster_means, clust_info, by="cluster") # attach clusterinfo incl. indicator variable
  data_clust_s2<- data_clust_s1 [data_clust_s1[phase_id[["phase.col"]]]  == phase_id[["terrgrid.id"]], ]
  rm(cluster_means)

  # get cluster local density Yc(x):
  Yc_x <- data_clust_s2[, all.vars(formula)[1] ]

  # get design_matrices for s1 and s2:
  design_matrix.s1<- as.matrix(data_clust_s1[ ,colnames.designmatrix ]) # Trick: indicator column in data_clust_sx only gets into the design-matrix.sx, if "unbiased" == T
  design_matrix.s2<- as.matrix(data_clust_s2[ ,colnames.designmatrix ])

  # get cluster-weights:
  M_x.s1<- data_clust_s1[["Mx"]]
  M_x.s2<- data_clust_s2[["Mx"]]

  # determine sample sizes:
  n1_clusters<- nrow(design_matrix.s1)
  n2_clusters<- nrow(design_matrix.s2)
  n1_plots<- nrow(data)
  n2_plots<- sum( data[[ phase_id[["phase.col"]] ]] == phase_id[["terrgrid.id"]] )


  ## ------- calculate everything for samples s1_G and s2_G: ---------------- ##

  # get design_matrices for s1 and s2 (has to be calculated individually since s1/s2-cluster may not lie entirely in G):
  data_ext_s1G<- data_ext[ data_ext[small_area[["sa.col"]]] == 1, ]

  # --> retrieve ZcG(x) for s1-sample IN G (has to be without indicator variable!):
  cluster_means_s1G<- aggregate(data_ext_s1G[ ,-c(which(names(data_ext_s1G)==cluster),which(names(data_ext_s1G)==small_area[["sa.col"]]),
                                                  which(names(data_ext_s1G)==phase_id[["phase.col"]]))],
                                list(cluster = data_ext_s1G[,cluster]), mean)


  # for weighted cluster-means:
  if(!is.na(boundary_weights)){
    # # library(plyr) #* load as dependency ...
    # weight the auxiliary information, but not the response:
    clustmeans_s1G.temp<- ddply(data_ext_s1G[,-c(which(names(data_ext_s1G) == small_area[["sa.col"]]), which(names(data_ext_s1G) == all.vars(formula)[1]),
                                                 which(names(data_ext_s1G)==phase_id[["phase.col"]]))],
                            .(cluster), function(x) colwise(weighted.mean, w = x[[boundary_weights]]) (x))
    clustmeans_s1G.temp<- clustmeans_s1G.temp[, -which(names(clustmeans_s1G.temp) == boundary_weights)] # remove column with boundary-weights, we don't need them anymore
    # merge unweighted cluster-means of response:
    cluster_means_s1G<- merge(clustmeans_s1G.temp, cluster_means_s1G[,c(cluster, all.vars(formula)[1])], by=cluster)
  }



  cluster_means_s1G<- merge(cluster_means_s1G, gridind_clust, by = cluster)
  cluster_means_s2G<- cluster_means_s1G[ cluster_means_s1G[phase_id[["phase.col"]]] == phase_id[["terrgrid.id"]], ]

  # local density on cluster-level of G:
  Yc_x_G<- cluster_means_s2G[ ,all.vars(formula)[1]]

  # M(x)_s1G for  s1G-sample:
  cluster_weights.s1G<- aggregate(data_ext_s1G[,all.vars(formula)[1]], list(cluster = data_ext_s1G[,cluster]), length)
  colnames(cluster_weights.s1G)[2]<- "Mx_s1G"

  # re-calculate indicator variable on cluster level IcG --> should now be "1" for all observations
  IcG.G<- aggregate(data_ext_s1G[[ small_area[["sa.col"]] ]],list(data_ext_s1G[,cluster]), mean)
  colnames(IcG.G)<- c(cluster,small_area[["sa.col"]])

  data_clust_s1G<- merge(cluster_means_s1G, IcG.G, by= cluster) # attach clusterinfo incl. indicator variable
  data_clust_s2G<- merge(cluster_means_s2G, IcG.G, by= cluster)

  # the cluster-weights have to be merged AFTER calculating data_clust_s1(s2)G, since the merge()-function in line 133/134
  # can change the row-order of cluster_means_s1G --> if this happens, cluster_weights.s1G do not match the right data_clust_s1G-rows:
  Mx_s1G<- merge(data_clust_s1G, cluster_weights.s1G, by = cluster)[["Mx_s1G"]]
  Mx_s2G<- merge(data_clust_s2G, cluster_weights.s1G, by = cluster)[["Mx_s1G"]]

  # get design_matrices for s1_G and s2_G:
  design_matrix.s1G<- as.matrix(data_clust_s1G[ ,colnames.designmatrix ])
  design_matrix.s2G<- as.matrix(data_clust_s2G[ ,colnames.designmatrix ])

  # determine sample sizes:
  n1G_clusters<- nrow(design_matrix.s1G)
  n2G_clusters<- nrow(design_matrix.s2G)
  n1G_plots<- nrow(data_ext_s1G)
  n2G_plots<- sum( data_ext_s1G[[ phase_id[["phase.col"]] ]] == phase_id[["terrgrid.id"]] )


  # ---------------------------- #
  # Catch special cases and "adjust" accordingly:

  # 1) if n1G_clusters and n2G_clusters have the same size, they are identical, and thus
  #     two-phase est. method cannot be applied (and doesnt make sense at all by the way)
  if (n1G_clusters == n2G_clusters){

    warning(paste("Two-phase estimation methods for Small Area", small_area[["areas"]],"not possible, set to 'NA': only terrestrial data available in Small Area"), call. = F)
    As2inv <- matrix(NA, nrow=ncol(design_matrix.s2), ncol= ncol(design_matrix.s2)) # since we cannot invert the design-matrix due to singularity
    r.squared<- NA # so we also set r.squared to NA
    design_matrix.s2G<- matrix(NA, nrow=1, ncol= ncol(design_matrix.s2))
    Yc_x_G<- NA
    Rc_x_hat_G<- NA
    Mx_s2G<- NA

  } else {


    if(n2G_clusters==0){
      # if "unbiased" == FALSE, we have no indicator introduced and hence we can invert the design-matrix although n2G == 0:
      if(!small_area[["unbiased"]]) {
        As2inv <- solve( (t(design_matrix.s2) %*% (M_x.s2*design_matrix.s2)) / n2_clusters ) # ... and we can invert the design-matrix
        design_matrix.s2G<- matrix(NA, nrow=1, ncol= ncol(design_matrix.s2))
        Yc_x_G<- NA
        Mx_s2G<- NA
      } else { # if we have the indicator introduced ...
        As2inv <- matrix(NA, nrow=ncol(design_matrix.s2), ncol= ncol(design_matrix.s2)) # since we cannot invert the design-matrix due to singularity
        r.squared<- NA # so we also set r.squared to NA
        design_matrix.s2G<- matrix(NA, nrow=1, ncol= ncol(design_matrix.s2))
        Yc_x_G<- NA
        Rc_x_hat_G<- NA
        Mx_s2G<- NA
        w1<- warning(paste("Unbiased estimation for Small Area", small_area[["areas"]],"not possible, set to 'NA': Small Area",small_area[["areas"]],"does not contain any terrestrial data"), call. = F)
      }
    }

    if(n2G_clusters!=0){
      As2inv <- solve( (t(design_matrix.s2) %*% (M_x.s2*design_matrix.s2)) / n2_clusters )
      #     if(!small_area[["unbiased"]] & n2G_plots > thresh.n2g){
      #         w2<- warning(paste("Terrestrial sample size in small area",small_area[["areas"]],"sufficient for unbiased estimation. Consider extended pseudosynthetic estimator instead."), call. = F)
      #       }
      if(small_area[["unbiased"]]){
        #       if(n2G_plots <= thresh.n2g){
        #        w3<- warning(paste("Terrestrial sample size in small area",small_area[["areas"]],"is considerably small. Consider true synthetic estimation instead"), call. = F)
        #       }
        if (any(IcG.s2 > 0 & IcG.s2 < 1)){
          if(!psmall){
            w4<- warning(
              paste("At least one terrestrial cluster not entirely included within the small area ",small_area[["areas"]],".\n",
                    "Zero mean residual assumption for small area maybe violated.","\n",
                    "Check mean_Rc_x_hat_G and consider alternative estimator 'psmall'", sep = ""), call. = F)
          }
        }
      }
    }

  }

  ## ------- calculate estimations ----------------------------------------- ##

  # regression coefficient theta: eq [63] p. 24 Mandallaz Small Area technical report
  theta_s2<- ((M_x.s2 * Yc_x) %*% design_matrix.s2 / n2_clusters) %*% As2inv

  # (mean) residuals in F:
  Yc_x_hat_s2<- design_matrix.s2 %*% t(theta_s2)
  Rc_x_hat<- Yc_x - Yc_x_hat_s2
  mean_Rc_x_hat<- weighted.mean(Rc_x_hat, w = M_x.s2)

  # (mean) residuals in G (only possible if we have s2-data in G!)
  Yc_x_hat_s1_G<- design_matrix.s1G %*% t(theta_s2)
  Yc_x_hat_s2_G<- design_matrix.s2G %*% t(theta_s2)
  Rc_x_hat_G<- Yc_x_G - Yc_x_hat_s2_G
  mean_Rc_x_hat_G<- weighted.mean(Rc_x_hat_G, w = Mx_s2G)

  # cov.matrix of theta_s2: eq [64] p. 24 Mandallaz Small Area technical report
  MR_square <- as.vector((M_x.s2^2)*(Rc_x_hat^2))
  middle_term<- ((t(design_matrix.s2) %*%  (MR_square * design_matrix.s2)) / n2_clusters^2)
  cov_theta_s2<- As2inv %*% middle_term %*% As2inv

  Z_bar_c1G<- apply(design_matrix.s1G, 2, weighted.mean, w=Mx_s1G)

  # prepare matrix calculation for cov_Z_bar_c1G by centering design matrix and scaling it with M_x ...eq [67] p. 25 Mandallaz Small Area technical report
  design_matrix.s1G_centered <- t(apply(design_matrix.s1G, 1, function(x){x - Z_bar_c1G}))
  design_matrix.s1G_centered_weighted <- apply(design_matrix.s1G_centered, 2, function(col){as.vector(Mx_s1G / mean(Mx_s1G))*col})
  cov_Z_bar_c1G <- t(design_matrix.s1G_centered_weighted) %*% design_matrix.s1G_centered_weighted / (n1G_clusters*(n1G_clusters-1))

  # estimates:
  estimate<- Z_bar_c1G %*% t(theta_s2)

  # Variances (only calculable if n1G_clusters > 1, set to NA otherwise!)
  if (n1G_clusters == 1) {

    warning(paste("Variance estimation for Small Area", small_area[["areas"]],"not possible, set to 'NA': Sample Size (n1G) in Small Area is one"), call. = F)
    g_variance<- NA
    ext_variance<- NA

  } else {

    g_variance<- (t(Z_bar_c1G) %*% cov_theta_s2 %*% Z_bar_c1G) + (theta_s2 %*% cov_Z_bar_c1G %*% t(theta_s2))

    ext_variance<- ( (1 - (n2G_clusters / n1G_clusters)) * (1/n2G_clusters) * (1/(n2G_clusters - 1)) * sum( ((Mx_s2G / mean(Mx_s2G))^2) * ((Rc_x_hat_G - mean_Rc_x_hat_G)^2) ) ) +
      ( (1 / n1G_clusters) * (1 / (n2G_clusters - 1)) * sum( ((Mx_s2G / mean(Mx_s2G))^2) * ((Yc_x_G - Yc_x_hat_s2_G)^2) ) )

  }


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
  est.type<- ifelse(small_area[["unbiased"]], "psynth extended", "psynth")

  # save warning-messages:
  warn.messages<- c(w1, w2, w3, w4)[which(!is.na(c(w1, w2, w3, w4)))]

return(
 list(regmodel=formula,
      estimation=estimation,
      samplesizes=samplesizes,
      coef=theta_s2,
      cov_coef=cov_theta_s2,
      Z_bar_1G=Z_bar_c1G,
      cov_Z_bar_1G=cov_Z_bar_c1G,
      Rc_x_hat_G = Rc_x_hat_G,
      Mx_s2G = Mx_s2G,
      mean_Rc_x_hat_G=mean_Rc_x_hat_G,
      boundary_weights=boundary_weights,
      estimator=est.type,
      warn.messages=warn.messages)
)


} # end of "small_area_nonexhaustive2p_cluster"


