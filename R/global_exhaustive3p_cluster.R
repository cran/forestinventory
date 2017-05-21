# This function helps the "two phase" cluster sampling estimator function when the first phase sample size goes to infinity.
# This corresponds to the situation where wall2wall info is known for post-stratification
# The external and g-weight variances are calulated


# INPUTS:
#   formula.s0 = formula for "reduced" model calculable for all of s0
#   formula.s1 = formula for "full" model calculable for all of s1
#   data = data.frame
#   boundary_weights = vector of weights representing the proportion of the forest area in the s1 interpretation area
#                      (vector of 1s if NA)
#
# OUTPUTS:
#   estimate = point estimate
#   ext_variance = external variance
#   g_variance = g-weight variance
#   n0 = Inf
#   n1 = sample size s1
#   n2 = sample size s2
#   r.squared_reduced = R squared of reduced model
#   r.squared_full = R squared of fullmodel


global_exhaustive3p_cluster <- function(formula.s0, formula.s1, data, phase_id, cluster, boundary_weights, exhaustive, ...){


  # retrieve phase.columnname and indicator of s1 grid id and terrestrial-grid id:
  phase.col<- phase_id[["phase.col"]]
  s1.ind<- phase_id[["s1.id"]]        # s1.id identifies the small sample of the auxiliary vars (s1-sample, refers to formula.s1, ie. the full-model)
  s2.ind<- phase_id[["terrgrid.id"]]  # identifies the terrestrial sample (s2-sample)

  # get number of plots for phases:
  n2<- sum(data[[phase.col]] == s2.ind)
  n1<- sum(data[,phase.col] %in% c(s1.ind, s2.ind))
  n0<- Inf


  # -------------------------------------------------------------------------- #
  # Prepare data on plot level #

  # compute data_ext_1.s1 i.e dataset of reduced model for s1:
  design_matrix_1.s1_plot_level <- design_matrix.s1_return(formula=formula.s0,
                                    data=data[data[[phase.col]] %in% c(s1.ind, s2.ind),]) # plot-level des.Mat of reduced model for s1
  data_ext_1.s1 <- data.frame(design_matrix_1.s1_plot_level)
  data_ext_1.s1[,all.vars(formula.s0)[1]] <- data[data[[phase.col]] %in% c(s1.ind, s2.ind),all.vars(formula.s0)[1]]
  data_ext_1.s1[,cluster] <- data[data[[phase.col]] %in% c(s1.ind, s2.ind),cluster]
  data_ext_1.s1[,phase.col] <- data[data[[phase.col]] %in% c(s1.ind, s2.ind),phase.col] # append phase-info
  if(!is.na(boundary_weights)){
    data_ext_1.s1[,boundary_weights]<- data[data[[phase.col]] %in% c(s1.ind, s2.ind),boundary_weights]
  }


  # compute data_ext.s1 i.e dataset of full model for s1:
  design_matrix.s1_plot_level <-  design_matrix.s1_return(
    formula=formula.s1, data=data[data[[phase.col]] %in% c(s1.ind, s2.ind),]) # plot-level des.Mat of full model for s1
  data_ext.s1 <- data.frame(design_matrix.s1_plot_level)
  data_ext.s1[,all.vars(formula.s0)[1]] <- data[data[[phase.col]] %in% c(s1.ind, s2.ind),all.vars(formula.s0)[1]]
  data_ext.s1[,cluster] <- data[data[[phase.col]] %in% c(s1.ind, s2.ind),cluster]
  data_ext.s1[,phase.col] <- data[data[[phase.col]] %in% c(s1.ind, s2.ind),phase.col] # append phase-info
  if(!is.na(boundary_weights)){
    data_ext.s1[,boundary_weights]<- data[data[[phase.col]] %in% c(s1.ind, s2.ind),boundary_weights]
  }


  # -------------------------------------------------------------------------- #

  # compute cluster weights M(x):
  cluster_weights <- aggregate(data_ext.s1[,all.vars(formula.s0)[1]],
                               list(cluster = data_ext.s1[,cluster]), length) # the M(x)

  # * Remark * #
  # --> cluster weights M(x) do not depend on reduced or full model. They're just computed for every cluster in the s1-sample,
  #     which are included in both the reduced and the full model dataset...

  # -------------------------------------------------------------------------- #

  # compute datasets on cluster level:

  ## ... of full model for s1:
  cluster_means.s1 <- aggregate(data_ext.s1[,-which(names(data_ext.s1)==cluster)],
                                list(cluster = data_ext.s1[,cluster]), mean)
  # for weighted cluster-means:
  if(!is.na(boundary_weights)){cluster_means.s1<- boundaryweight_fct_3p_clust(formula.s1, data_ext.s1, boundary_weights, cluster)}
  data_clust.s1 <- merge(cluster_means.s1, cluster_weights, by=cluster)  #right hand column is M(x)


  ## ... of reduced model for s1:
  cluster_means_1.s1 <- aggregate(data_ext_1.s1[,-which(names(data_ext_1.s1)==cluster)],
                                  list(cluster = data_ext_1.s1[,cluster]), mean)
  # for weighted cluster-means:
  if(!is.na(boundary_weights)){cluster_means_1.s1<- boundaryweight_fct_3p_clust(formula.s0, data_ext_1.s1, boundary_weights, cluster)}
  data_clust_1.s1 <- merge(cluster_means_1.s1, cluster_weights, by=cluster)  #right hand column is M(x)


  ## ... of full model for s2:
  data_clust.s2<- data_clust.s1[ data_clust.s1[[phase.col]] == s2.ind, ]
  Yc_x<- data_clust.s2[, all.vars(formula.s0)[1]]


  ## of reduced model for s2:
  data_clust_1.s2<- data_clust_1.s1[ data_clust_1.s1[[phase.col]] == s2.ind, ]


  # -------------------------------------------------------------------------- #

  # extract design-matrices (Zc's):

  ## ... of full model for s1:
  design_matrix.s1<- as.matrix(data_clust.s1[ , -c(which(names(data_clust.s1) %in% c(cluster,all.vars(formula.s0)[1], "x", phase.col)))]) # Zc(x)

  ## ... of reduced model for s1:
  design_matrix_1.s1<- as.matrix(data_clust_1.s1[ , -c(which(names(data_clust_1.s1) %in% c(cluster,all.vars(formula.s0)[1], "x", phase.col)))])

  ## ... of full model for s2:
  design_matrix.s2<- as.matrix(data_clust.s2[ , -c(which(names(data_clust.s2) %in% c(cluster,all.vars(formula.s0)[1], "x", phase.col)))])

  ## ... of reduced model for s2:
  design_matrix_1.s2<- as.matrix(data_clust_1.s2[ , -c(which(names(data_clust_1.s2) %in% c(cluster,all.vars(formula.s0)[1], "x", phase.col)))])


  # extract and store the M(x)
  M_x.s1 <- data_clust.s1[,"x"]
  M_x.s2 <- data_clust.s2[,"x"]

  # get sample size:
  n0_clusters <- Inf
  n1_clusters <- nrow(data_clust.s1)
  n2_clusters <- nrow(data_clust.s2)

  # -------------------------------------------------------------------------- #

  # calculate the Z_bars...:
  Z_bar_s1<-   apply(design_matrix.s1,   2, weighted.mean, w=M_x.s1)   # that's Z_bar_1 in Daniels Report...
  Z_1_bar_s1<- apply(design_matrix_1.s1, 2, weighted.mean, w=M_x.s1) # that's Z_bar(1)_1 in Daniels Report...


  # the A's:
  A_1_s2_inv<- solve( (t(design_matrix_1.s2) %*%  (M_x.s2*design_matrix_1.s2)) / n2_clusters )
  A_s2_inv<-   solve( (t(design_matrix.s2)   %*%  (M_x.s2*design_matrix.s2))   / n2_clusters )

  # -------------------------------------------------------------------------- #

  # regression coefficients:
  beta<-  A_s2_inv   %*% t((M_x.s2 * Yc_x) %*% design_matrix.s2   / n2_clusters) # reg.coef of full model
  alpha<- A_1_s2_inv %*% t((M_x.s2 * Yc_x) %*% design_matrix_1.s2 / n2_clusters) # reg.coef of reduced model

  Yc_x_hat<-   design_matrix.s2 %*% beta    # predictions of full model
  Yc_1_x_hat<- design_matrix_1.s2 %*% alpha # predictions of reduced model

  Rc_x_hat<-   Yc_x - Yc_x_hat   # residuals of full model
  Rc_1_x_hat<- Yc_x - Yc_1_x_hat # residuals of reduced model
  mean_Rc_x_hat<-   weighted.mean(Rc_x_hat, w = M_x.s2)
  mean_Rc_1_x_hat<- weighted.mean(Rc_1_x_hat, w = M_x.s2)

  # variance-covariance-matrix of beta:
  MR_square <- as.vector((M_x.s2^2)*(Rc_x_hat^2))
  middle_term <- ((t(design_matrix.s2)) %*% (MR_square*(design_matrix.s2))) / n2_clusters^2
  cov_beta_s2 <- A_s2_inv %*% middle_term %*% A_s2_inv

  # variance-covariance-matrix of alpha:
  MR_square_1<- as.vector((M_x.s2^2)*(Rc_1_x_hat^2))
  middle_term_1<- ((t(design_matrix_1.s2)) %*% (MR_square_1*(design_matrix_1.s2))) / n2_clusters^2
  cov_alpha_s2 <- A_1_s2_inv %*% middle_term_1 %*% A_1_s2_inv


  # -------------------------------------------------------------------------- #
  # estimates:

  estimate<- ((exhaustive - Z_1_bar_s1) %*% alpha) + (Z_bar_s1 %*% beta)

  # external variance is an extension of equation[9], page 6 in two-phase with exh.information (Mandallaz) to cluster sampling:
  ext_variance<- (1 / (n1_clusters * (n2_clusters - 1))) * sum( (M_x.s2 / mean(M_x.s2))^2 * (Rc_1_x_hat - mean_Rc_1_x_hat)^2) +
                 (1 - (n2_clusters / n1_clusters)) * (1 / (n2_clusters * (n2_clusters-1) )) * sum( (M_x.s2 / mean(M_x.s2))^2 * (Rc_x_hat - mean_Rc_x_hat)^2)

  g_variance<- (n2_clusters / n1_clusters) * (exhaustive %*% cov_alpha_s2 %*% exhaustive) +
               (1 - (n2_clusters / n1_clusters)) * (Z_bar_s1 %*% cov_beta_s2 %*% Z_bar_s1)


  ## ------- create outputs ------------------------------------------------- ##


  # summarize sample size info:
  samplesizes<- data.frame(cbind (n0_clusters,n1_clusters, n2_clusters, n0, n1, n2))
  colnames(samplesizes)<- c("n0_clust", "n1_clust", "n2_clust", "n0", "n1", "n2")
  rownames(samplesizes)<- "plots"

  estimation<- data.frame(estimate=estimate, ext_variance=ext_variance, g_variance=g_variance,
                          n0=samplesizes$n0, n1=samplesizes$n1, n2=samplesizes$n2,
                          n0_clust=samplesizes$n0_clust, n1_clust=samplesizes$n1_clust, n2_clust=samplesizes$n2_clust,
                          r.squared_reduced=summary(lm(formula.s0, data=data[data[[phase.col]] %in% s2.ind,]))$r.squared,
                          r.squared_full=summary(lm(formula.s1, data=data[data[[phase.col]] %in% s2.ind,]))$r.squared)


  # ... to store inputs used:
  inputs<- list()
  inputs[["data"]]<- data
  inputs[["formula.s0"]]<- formula.s0
  inputs[["formula.s1"]]<- formula.s1
  inputs[["boundary_weights"]]<- boundary_weights
  inputs[["exhaustive"]]<- exhaustive
  inputs[["method"]]<- "exhaustive"
  inputs[["cluster"]]<- TRUE
  inputs[["exhaustive"]]<- TRUE

  # save warning-messages:
  warn.messages<- NA

  result<- list(input=inputs,
                estimation=estimation,
                samplesizes=samplesizes,
                coefficients=list(alpha=t(alpha)[1,], beta=t(beta)[1,]),
                cov_coef=list(cov_beta_s2=cov_beta_s2, cov_alpha_s2=cov_alpha_s2),
                Rc_x_hat=Rc_x_hat,
                Rc_1_x_hat = Rc_1_x_hat,
                mean_Rc_x_hat=mean_Rc_x_hat,
                mean_Rc_1_x_hat = mean_Rc_1_x_hat,
                warn.messages=warn.messages)

  class(result)<- c("global" ,"threephase")

  return(result)

}


