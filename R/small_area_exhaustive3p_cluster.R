



small_area_exhaustive3p_cluster<- function(formula.s0, formula.s1, data, phase_id, cluster, small_area, boundary_weights, exhaustive, thresh.n2g, psmall){


  # intialize warning-message vector (i.e. collector of warnings produced during function execution):
  w1<- NA; w2<- NA; w3<- NA; w4<- NA

  # -------------------------------------------------------------------------- #
  # retrieve phase.columnname and indicator of s1 grid id and terrestrial-grid id:
  phase.col<- phase_id[["phase.col"]]
  sa.col<- small_area[["sa.col"]]
  s1.ind<- phase_id[["s1.id"]]        # s1.id identifies the small sample of the auxiliary vars (s1-sample, refers to formula.s1, ie. the full-model)
  s2.ind<- phase_id[["terrgrid.id"]]  # identifies the terrestrial sample (s2-sample)


  # calculate the r-squares for both models (on plot-level):
  r.squared_reduced<- summary(lm(formula.s0, data=data[data[[phase.col]] %in% s2.ind,]))$r.squared
  r.squared_full<- summary(lm(formula.s1, data=data[data[[phase.col]] %in% s2.ind,]))$r.squared

  # -------------------------------------------------------------------------- #


  # -------------------------------------------------------------------------- #
  ## ------- calculate everything for samples s1 and s2: -------------------- ##
  # -------------------------------------------------------------------------- #


  # Prepare data on plot level #

  # compute data_ext_1.s1 i.e dataset of reduced model for s1-sample:
  design_matrix_1.s1_plot_level <- design_matrix.s1_return(formula=formula.s0, data=data[data[[phase.col]] %in% c(s1.ind, s2.ind),]) # plot-level des.Mat of reduced model for s1-sample
  data_ext_1.s1 <- data.frame(design_matrix_1.s1_plot_level)
  data_ext_1.s1[,all.vars(formula.s0)[1]] <- data[data[[phase.col]] %in% c(s1.ind, s2.ind),all.vars(formula.s0)[1]]
  data_ext_1.s1[,cluster] <- data[data[[phase.col]] %in% c(s1.ind, s2.ind),cluster]
  data_ext_1.s1[,phase.col] <- data[data[[phase.col]] %in% c(s1.ind, s2.ind),phase.col]
  if(!is.na(boundary_weights)){
    data_ext_1.s1[,boundary_weights]<- data[data[[phase.col]] %in% c(s1.ind, s2.ind),boundary_weights]
  }

  # compute data_ext.s1 i.e dataset of full model for s1-sample:
  design_matrix.s1_plot_level <-  design_matrix.s1_return(formula=formula.s1, data=data[data[[phase.col]] %in% c(s1.ind, s2.ind),]) # plot-level des.Mat of full model for s1
  data_ext.s1 <- data.frame(design_matrix.s1_plot_level)
  data_ext.s1[,all.vars(formula.s0)[1]] <- data[data[[phase.col]] %in% c(s1.ind, s2.ind),all.vars(formula.s0)[1]]
  data_ext.s1[,cluster] <- data[data[[phase.col]] %in% c(s1.ind, s2.ind),cluster]
  data_ext.s1[,phase.col] <- data[data[[phase.col]] %in% c(s1.ind, s2.ind),phase.col]
  if(!is.na(boundary_weights)){
    data_ext.s1[,boundary_weights]<- data[data[[phase.col]] %in% c(s1.ind, s2.ind),boundary_weights]
  }

  # -------------------------------------------------------------------------- #

  # compute cluster weights M(x) for entire sample in F:
  cluster_weights <- aggregate(data_ext.s1[,all.vars(formula.s0)[1]], list(cluster = data_ext.s1[,cluster]), length) # the M(x)
  colnames(cluster_weights)[2]<- "Mx"

  # -------------------------------------------------------------------------- #


  # compute s1- and s2- datasets on cluster level:

  ## ... of full model for s1:
  cluster_means.s1 <- aggregate(data_ext.s1[,-which(names(data_ext.s1) == cluster)], list(cluster = data_ext.s1[,cluster]), mean)
  # for weighted cluster-means:
  if(!is.na(boundary_weights)){cluster_means.s1<- boundaryweight_fct_3p_clust(formula.s1, data_ext.s1, boundary_weights, cluster)}
  data_clust.s1 <- merge(cluster_means.s1, cluster_weights, by=cluster)  # cluster_weights is M(x)


  ## ... of reduced model for s1:
  cluster_means_1.s1 <- aggregate(data_ext_1.s1[,-which(names(data_ext_1.s1)==cluster)], list(cluster = data_ext_1.s1[,cluster]), mean)
  # for weighted cluster-means:
  if(!is.na(boundary_weights)){cluster_means_1.s1<- boundaryweight_fct_3p_clust(formula.s0, data_ext_1.s1, boundary_weights, cluster)}
  data_clust_1.s1 <- merge(cluster_means_1.s1, cluster_weights, by=cluster)  # right hand column is M(x)


  ## ... of full model for s2:
  data_clust.s2<- data_clust.s1[ data_clust.s1[[phase.col]] == s2.ind, ]
  Yc_x<- data_clust.s2[, all.vars(formula.s0)[1]]


  ## of reduced model for s2:
  data_clust_1.s2<- data_clust_1.s1[ data_clust_1.s1[[phase.col]] == s2.ind, ]


  # -------------------------------------------------------------------------- #

  # extract design-matrices (Zc's):

  ## ... of full model for s1:
  design_matrix.s1<- as.matrix(data_clust.s1[ , -c(which(names(data_clust.s1) %in% c(cluster,all.vars(formula.s0)[1], "Mx", phase.col)))]) # Zc(x)

  ## ... of reduced model for s1:
  design_matrix_1.s1<- as.matrix(data_clust_1.s1[ , -c(which(names(data_clust_1.s1) %in% c(cluster,all.vars(formula.s0)[1], "Mx", phase.col)))])# Zc_1(x)

  ## ... of full model for s2:
  design_matrix.s2<- as.matrix(data_clust.s2[ , -c(which(names(data_clust.s2) %in% c(cluster,all.vars(formula.s0)[1], "Mx", phase.col)))]) # Zc(x) for s2-sample

  ## ... of reduced model for s2:
  design_matrix_1.s2<- as.matrix(data_clust_1.s2[ , -c(which(names(data_clust_1.s2) %in% c(cluster,all.vars(formula.s0)[1], "Mx", phase.col)))]) # Zc_1(x) for s2-sample


  # extract and store the M(x)
  M_x.s1<- data_clust.s1[["Mx"]]
  M_x.s2<- data_clust.s2[["Mx"]]


  # get samplesizes for phases in F:
  n2_plots<- sum(data[[phase.col]] == s2.ind)
  n1_plots<- sum(data[,phase.col] %in% c(s1.ind, s2.ind))
  n0_plots<- Inf
  n2_clusters <- nrow(design_matrix.s2)
  n1_clusters <- nrow(design_matrix.s1)
  n0_clusters <- Inf

  # -------------------------------------------------------------------------- #


  # ---------------------------------------------------------------------------#
  ## ------- calculate everything for samples s1G and s2G: ------------------ ##
  # ---------------------------------------------------------------------------#


  # Prepare data on plot level #

  # compute data_ext_1.s1G i.e dataset of reduced model for s1G-sample:
  design_matrix_1.s1G_plot_level <- design_matrix.s1_return(formula=formula.s0, data=data[data[,phase.col] %in% c(s1.ind, s2.ind) & data[,sa.col] == 1 ,])
  data_ext_1.s1G <- data.frame(design_matrix_1.s1G_plot_level)
  data_ext_1.s1G[,all.vars(formula.s0)[1]] <- data[data[,phase.col] %in% c(s1.ind, s2.ind) & data[,sa.col] == 1, all.vars(formula.s0)[1]]
  data_ext_1.s1G[,cluster] <- data[data[,phase.col] %in% c(s1.ind, s2.ind) & data[,sa.col] == 1, cluster]
  data_ext_1.s1G[,phase.col] <- data[data[,phase.col] %in% c(s1.ind, s2.ind) & data[,sa.col] == 1,phase.col]
  if(!is.na(boundary_weights)){
    data_ext_1.s1G[,boundary_weights]<- data[data[,phase.col] %in% c(s1.ind, s2.ind) & data[,sa.col] == 1, boundary_weights]
  }


  # compute data_ext.s1 i.e dataset of full model for s1G-sample:
  design_matrix.s1G_plot_level <- design_matrix.s1_return(formula=formula.s1, data=data[data[,phase.col] %in% c(s1.ind, s2.ind) & data[,sa.col] == 1 ,])
  data_ext.s1G <- data.frame(design_matrix.s1G_plot_level)
  data_ext.s1G[,all.vars(formula.s0)[1]] <- data[data[,phase.col] %in% c(s1.ind, s2.ind) & data[,sa.col] == 1, all.vars(formula.s0)[1]]
  data_ext.s1G[,cluster] <- data[data[,phase.col] %in% c(s1.ind, s2.ind) & data[,sa.col] == 1, cluster]
  data_ext.s1G[,phase.col] <- data[data[,phase.col] %in% c(s1.ind, s2.ind) & data[,sa.col] == 1,phase.col]
  if(!is.na(boundary_weights)){
    data_ext.s1G[,boundary_weights]<- data[data[,phase.col] %in% c(s1.ind, s2.ind) & data[,sa.col] == 1, boundary_weights]
  }

  ## compute data_ext_1.s2G i.e dataset of reduced model for s2G-sample:
  data_ext_1.s2G<- data_ext_1.s1G[data_ext_1.s1G[[phase.col]]==s2.ind, ]

  ## compute data_ext.s2 i.e dataset of full model for s2G-sample:
  data_ext.s2G<- data_ext.s1G[data_ext.s1G[[phase.col]]==s2.ind, ]


  # -------------------------------------------------------------------------- #

  # compute s1G- and s2G- datasets on cluster level:

  # Mx_s1G (same for full and reduced model):
  cluster_weights.s1G<- aggregate(data_ext.s1G[,all.vars(formula.s0)[1]], list(cluster = data_ext.s1G[,cluster]), length) # Mx_s1G
  colnames(cluster_weights.s1G)[2]<- "Mx_s1G"

  # Mx_s2G (same for full and reduced model):
  cluster_weights.s2G<- aggregate(data_ext.s2G[,all.vars(formula.s0)[1]], list(cluster = data_ext.s2G[,cluster]), length) # Mx_s1G
  colnames(cluster_weights.s2G)[2]<- "Mx_s2G"


  ## ... of full model for s1G:
  data_clust.s1G <- aggregate(data_ext.s1G[,-which(names(data_ext.s1G) == cluster)], list(cluster = data_ext.s1G[,cluster]), mean)
  # for weighted cluster-means:
  if(!is.na(boundary_weights)){data_clust.s1G<- boundaryweight_fct_3p_clust(formula.s1, data_ext.s1G, boundary_weights, cluster)}
  data_clust.s1G<- merge(data_clust.s1G, cluster_weights.s1G, by=cluster)


  ## ... of reduced model for s1G:
  data_clust_1.s1G <- aggregate(data_ext_1.s1G[,-which(names(data_ext_1.s1G) == cluster)], list(cluster = data_ext_1.s1G[,cluster]), mean)
  # for weighted cluster-means:
  if(!is.na(boundary_weights)){data_clust_1.s1G<- boundaryweight_fct_3p_clust(formula.s0, data_ext_1.s1G, boundary_weights, cluster)}
  data_clust_1.s1G<- merge(data_clust_1.s1G, cluster_weights.s1G, by=cluster)


  ## ... of full model for s2G:
  data_clust.s2G <- aggregate(data_ext.s2G[,-which(names(data_ext.s2G) == cluster)], list(cluster = data_ext.s2G[,cluster]), mean)
  # for weighted cluster-means:
  if(!is.na(boundary_weights)){data_clust.s2G<- boundaryweight_fct_3p_clust(formula.s1, data_ext.s2G, boundary_weights, cluster)}
  data_clust.s2G<- merge(data_clust.s2G, cluster_weights.s2G, by=cluster)


  ## ... of reduced model for s2G:
  data_clust_1.s2G <- aggregate(data_ext_1.s2G[,-which(names(data_ext_1.s2G) == cluster)], list(cluster = data_ext_1.s2G[,cluster]), mean)
  # for weighted cluster-means:
  if(!is.na(boundary_weights)){data_clust_1.s2G<- boundaryweight_fct_3p_clust(formula.s0, data_ext_1.s2G, boundary_weights, cluster)}
  data_clust_1.s2G<- merge(data_clust_1.s2G, cluster_weights.s2G, by=cluster)

  Yc_x_G<- data_clust.s2G[, all.vars(formula.s0)[1]]

  # -------------------------------------------------------------------------- #

  # extract design-matrices (Zc's):

  ## ... of full model for s1G:
  design_matrix.s1G<- as.matrix(data_clust.s1G[ , -c(which(names(data_clust.s1G) %in% c(cluster,all.vars(formula.s0)[1], "Mx_s1G", phase.col)))]) # Zc(x)

  ## ... of reduced model for s1G:
  design_matrix_1.s1G<- as.matrix(data_clust_1.s1G[ , -c(which(names(data_clust_1.s1G) %in% c(cluster,all.vars(formula.s0)[1], "Mx_s1G", phase.col)))])# Zc_1(x)

  ## ... of full model for s2G:
  design_matrix.s2G<- as.matrix(data_clust.s2G[ , -c(which(names(data_clust.s2G) %in% c(cluster,all.vars(formula.s0)[1], "Mx_s2G", phase.col)))])# Zc_1(x)

  ## ... of reduced model for s2G:
  design_matrix_1.s2G<- as.matrix(data_clust_1.s2G[ , -c(which(names(data_clust_1.s2G) %in% c(cluster,all.vars(formula.s0)[1], "Mx_s2G", phase.col)))])# Zc_1(x)


  # extract and store the M(x)_G
  M_x.s1G<- data_clust.s1G[["Mx_s1G"]]
  M_x.s2G<- data_clust.s2G[["Mx_s2G"]]


  # get sample size:
  n2G_plots<-  sum(data_ext.s1G[[phase.col]] == s2.ind)
  n1G_plots<-  nrow(data_ext.s1G)
  n0G_plots<-  Inf
  n0G_clusters <- Inf
  n1G_clusters <- nrow(design_matrix.s1G)
  n2G_clusters <- sum(data_clust_1.s1G[[phase.col]] == s2.ind)


  # -------------------------------------------------------------------------- #


  # calculate the A's: --> depends on cases "no. n2G in G"

  # ------------------------------------ #
  # 1) case "no n2g in G":

  if(n2G_clusters == 0){

    if(!small_area[["unbiased"]]) {
      # if "unbiased" == FALSE (i.e. if synthetic est. applied):

      # calculate the A's: --> we have no indicator introduced and hence we can invert the design-matrix although n2G == 0:
      A_1_s2_inv<- solve( (t(design_matrix_1.s2) %*%  (data_clust_1.s2[["Mx"]] * design_matrix_1.s2)) / n2_clusters )
      A_1_s1_inv<- solve( (t(design_matrix_1.s1) %*%  (data_clust_1.s1[["Mx"]] * design_matrix_1.s1)) / n1_clusters )
      A_s2_inv<-   solve( (t(design_matrix.s2)   %*%  (data_clust.s2[["Mx"]] * design_matrix.s2))   / n2_clusters )

      # set design-matrix to NA:
      design_matrix.s2G<-   matrix(NA, nrow=1, ncol= ncol(design_matrix.s2))
      design_matrix_1.s2G<- matrix(NA, nrow=1, ncol= ncol(design_matrix_1.s2))

      # set local densities and Mx_s2G for s2G to NA:
      Yc_x_G<- NA
      M_x.s2G<- NA

    } else { # if "unbiased" == TRUE:

      # calculate the A's: not possible, since we cannot invert the design-matrix:
      A_1_s2_inv<- matrix(NA, nrow=ncol(design_matrix_1.s2), ncol=ncol(design_matrix_1.s2))
      A_1_s1_inv<- matrix(NA, nrow=ncol(design_matrix_1.s1), ncol=ncol(design_matrix_1.s1))
      A_s2_inv<-   matrix(NA, nrow=ncol(design_matrix.s2), ncol=ncol(design_matrix.s2))

      # set design-matrix to NA:
      design_matrix.s2G<-   matrix(NA, nrow=1, ncol= ncol(design_matrix.s2))
      design_matrix_1.s2G<- matrix(NA, nrow=1, ncol= ncol(design_matrix_1.s2))

      # set local densities and Mx_s2G for s2G to NA:
      Yc_x_G<- NA
      M_x.s2G<- NA

      # set r-squares to NA to force entire output to appear as 'NA':
      r.squared_reduced<- NA
      r.squared_full<- NA

      w1<- warning(paste("Unbiased estimation for Small Area", small_area[["areas"]],"not possible, set to 'NA': Small Area",small_area[["areas"]],"does not contain any terrestrial data"), call. = F)

    }
  }


  # ------------------------------------ #
  # 2) case "only 1 n2g in G":

  if(n2G_clusters ==1){

    # calculate the A's:
    A_1_s2_inv<- solve( (t(design_matrix_1.s2) %*%  (data_clust_1.s2[["Mx"]] * design_matrix_1.s2)) / n2_clusters )
    A_1_s1_inv<- solve( (t(design_matrix_1.s1) %*%  (data_clust_1.s1[["Mx"]] * design_matrix_1.s1)) / n1_clusters )
    A_s2_inv<-   solve( (t(design_matrix.s2)   %*%  (data_clust.s2[["Mx"]] * design_matrix.s2))   / n2_clusters )

#     A_1_s2_inv<- matrix(NA, nrow=ncol(design_matrix_1.s2), ncol=ncol(design_matrix_1.s2))
#     A_1_s1_inv<- matrix(NA, nrow=ncol(design_matrix_1.s1), ncol=ncol(design_matrix_1.s1))
#     A_s2_inv<-   matrix(NA, nrow=ncol(design_matrix.s2), ncol=ncol(design_matrix.s2))

    w2<- warning(paste("Terrestrial sample size in small area",small_area[["areas"]],"is 1. Computation of external variance not possible"), call. = F)

  }


  # ------------------------------------ #
  # 3) "normal" case (few n2g available):

  if(n2G_clusters != 0){

    # calculate the A's:
    A_1_s2_inv<- solve( (t(design_matrix_1.s2) %*%  (data_clust_1.s2[["Mx"]] * design_matrix_1.s2)) / n2_clusters )
    A_1_s1_inv<- solve( (t(design_matrix_1.s1) %*%  (data_clust_1.s1[["Mx"]] * design_matrix_1.s1)) / n1_clusters )
    A_s2_inv<-   solve( (t(design_matrix.s2)   %*%  (data_clust.s2[["Mx"]] * design_matrix.s2))   / n2_clusters )

    if(small_area[["unbiased"]]){

      if (any(data_clust.s2G[[sa.col]] > 0 & data_clust.s2G[[sa.col]] < 1) | any(data_clust_1.s2G[[sa.col]] > 0 & data_clust_1.s2G[[sa.col]] < 1)){ # check with example

        if(!psmall){
          w4<- warning(
            paste("At least one terrestrial cluster not entirely included within the small area ",small_area[["areas"]],".\n",
                  "Zero mean residual assumption for small area maybe violated.","\n",
                  "Check 'mean_Rc_1_x_hat_G', 'mean_Rc_x_hat_G' and consider alternative estimator 'psmall'", sep = ""), call. = F)
        }

      }
    }
  }

  # ------------------------------------ #


  # -------------------------------------------------------------------------- #

  # regression coefficients:
  beta<-  A_s2_inv   %*% t((M_x.s2 * Yc_x) %*% design_matrix.s2   / n2_clusters) # reg.coef of full model
  alpha<- A_1_s2_inv %*% t((M_x.s2 * Yc_x) %*% design_matrix_1.s2 / n2_clusters) # reg.coef of reduced model

  # (mean) residuals in F:
  Yc_x_hat<-   design_matrix.s2 %*% beta    # predictions of full model
  Yc_1_x_hat<- design_matrix_1.s2 %*% alpha # predictions of reduced model
  Rc_x_hat<-   Yc_x - Yc_x_hat   # residuals of full model
  Rc_1_x_hat<- Yc_x - Yc_1_x_hat # residuals of reduced model
  mean_Rc_x_hat<-   weighted.mean(Rc_x_hat, w = M_x.s2)
  mean_Rc_1_x_hat<- weighted.mean(Rc_1_x_hat, w = M_x.s2)


  # (mean) residuals in G (only possible if we have s2-data in G!)
  Yc_x_hat_s2G<-    design_matrix.s2G %*% beta     # predictions of full model
  Yc_1_x_hat_s2G<-  design_matrix_1.s2G %*% alpha  # predictions of reduced model
  Rc_x_hat_G<-   Yc_x_G - Yc_x_hat_s2G   # residuals of full model
  Rc_1_x_hat_G<- Yc_x_G - Yc_1_x_hat_s2G # residuals of reduced model
  mean_Rc_x_hat_G<-   weighted.mean(Rc_x_hat_G, w = M_x.s2G)
  mean_Rc_1_x_hat_G<- weighted.mean(Rc_1_x_hat_G, w = M_x.s2G)

  # variance-covariance-matrix of beta:
  MR_square <- as.vector((M_x.s2^2)*(Rc_x_hat^2))
  middle_term <- ((t(design_matrix.s2)) %*% (MR_square*(design_matrix.s2))) / n2_clusters^2
  cov_beta_s2 <- A_s2_inv %*% middle_term %*% A_s2_inv

  # variance-covariance-matrix of alpha:
  MR_square_1<- as.vector((M_x.s2^2)*(Rc_1_x_hat^2))
  middle_term_1<- ((t(design_matrix_1.s2)) %*% (MR_square_1*(design_matrix_1.s2))) / n2_clusters^2
  cov_alpha_s2 <- A_1_s1_inv %*% middle_term_1 %*% A_1_s1_inv

  # -------------------------------------------------------------------------- #

  # calculate the Z_bars...:
  Z_bar_s1G<-   apply(design_matrix.s1G,   2, weighted.mean, w=M_x.s1G)   # that's Z_bar_1 in Daniels Report...
  Z_1_bar_s1G<- apply(design_matrix_1.s1G, 2, weighted.mean, w=M_x.s1G) # that's Z_bar(1)_1 in Daniels Report...
  if(small_area[["unbiased"]]){
  Z1_G <- c(unlist(exhaustive),1)
  } else {
  Z1_G<- unlist(exhaustive)}

  # -------------------------------------------------------------------------- #


  # ---------------------------------------------------------------------------#
  ## ------- estimations: --------------------------------------------------- ##
  # ---------------------------------------------------------------------------#

  # db-estimates:
  estimate<- ((Z1_G - Z_1_bar_s1G) %*% alpha) + (Z_bar_s1G %*% beta)

  g_variance<- (n2_clusters / n1_clusters) * (Z1_G %*% cov_alpha_s2 %*% Z1_G) +
    (1 - (n2_clusters / n1_clusters)) * (Z_bar_s1G %*% cov_beta_s2 %*% Z_bar_s1G)


  # external variance:
  ext_variance<- ( 1 / (n1G_clusters)) * ( 1 / (n2G_clusters - 1)) * sum( ((M_x.s2G / mean(M_x.s2G))^2) * ((Rc_1_x_hat_G - mean_Rc_1_x_hat_G)^2) ) +
                 ( 1 - (n2G_clusters/n1G_clusters)) * ( 1 / (n2G_clusters * (n2G_clusters - 1))) * sum( ((M_x.s2G / mean(M_x.s2G))^2) * ((Rc_x_hat_G - mean_Rc_x_hat_G)^2) )


  # ---------------------------------------------------------------------------#
  ## ------- create outputs:------------------------------------------------- ##
  # ---------------------------------------------------------------------------#



  # summarize sample size info:
  samplesizes<- data.frame( rbind(cbind (n0G_clusters, n1G_clusters, n2G_clusters, n0_clusters, n1_clusters, n2_clusters),
                                  cbind (n0G_plots, n1G_plots, n2G_plots, n0_plots, n1_plots, n2_plots)))
  colnames(samplesizes)<- c("n0G", "n1G", "n2G", "n0", "n1", "n2")
  rownames(samplesizes)<- c("clusters", "plots")


  estimation<- data.frame(estimate=estimate, ext_variance=ext_variance, g_variance=g_variance,
                          n0=n0_clusters, n1=n1_clusters, n2=n2_clusters,
                          n0G=n0G_clusters, n1G=n1G_clusters, n2G=n2G_clusters,
                          r.squared_reduced=r.squared_reduced,
                          r.squared_full=r.squared_full)



  # est.type can be "unbiased" in case of using the indicator, and "biased" in case of the true synthetic estimator:
  est.type<- ifelse(small_area[["unbiased"]], "synth extended", "synth")

  # save warning-messages:
  warn.messages<- c(w1, w2, w3)[which(!is.na(c(w1, w2, w3, w4)))]

  result<- list(regmodel.s0=formula.s0,
                regmodel.s1=formula.s1,
                estimation=estimation,
                samplesizes=samplesizes,
                coefficients=list(alpha=alpha, beta=beta),
                cov_coef=list(cov_beta_s2=cov_beta_s2, cov_alpha_s2=cov_alpha_s2),
                cov_Z1_bar_0G=NA,  # there is no uncertainty in the aux. variables in the exhaustive case
                Z1_bar_0G=Z1_G,    # in the exhaustive case, Z1_bar_0G is the vector of the known exhaustive means
                Z1_bar_1G= Z_1_bar_s1G,
                Z_bar_1G= Z_bar_s1G,
                Rc_x_hat_G=list(resid_reduced_G=Rc_1_x_hat_G, resid_full_G=Rc_x_hat_G),
                mean_Rc_x_hat_G=list(mean_resid_reduced_G=mean_Rc_1_x_hat_G, mean_resid_full_G=mean_Rc_x_hat_G),
                Mx_s2G = M_x.s2G,
                boundary_weights=boundary_weights,
                estimator=est.type,
                warn.messages=warn.messages)


  result


}
