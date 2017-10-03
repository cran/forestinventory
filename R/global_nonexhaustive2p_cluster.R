# This function helps the "two phase" cluster sampling estimator function when the first phase sample size is finite.
# This corresponds to the situation where wall2wall info is known for post-stratification
# The external (large sample version for linear models) and g-weight variances are calulated

# INPUTS:
#   formula = model formula
#   data = data.frame (plot level)
#   exhaustive = the vector of auxiliary variable means known exhaustively (should ideally be weighted mean adjusted for proportion of pixels in interpretation are in the forest)
#
# OUTPUTS:
#   estimate = point estimate
#   ext_variance = external variance
#   g_variance = g-weight variance
#   n1 = Inf
#   n2 = sample size
#   r.squared = R squared of model

global_nonexhaustive2p_cluster <- function(formula, data, phase_id, cluster, boundary_weights){

  n2<- sum(data[[phase_id[["phase.col"]]]] == phase_id[["terrgrid.id"]])
  n1<- nrow(data)

  design_matrix.s1_plot_level <- design_matrix.s1_return(formula=formula, data=data)
  data_ext <- data.frame(design_matrix.s1_plot_level)
  data_ext[,all.vars(formula)[1]] <- data[,all.vars(formula)[1]]
  data_ext[,cluster] <- data[,cluster]
  if(!is.na(boundary_weights)){
    data_ext[,boundary_weights]<- data[,boundary_weights]
  }

  cluster_weights <- aggregate(data_ext[,all.vars(formula)[1]], list(cluster = data_ext[,cluster]), length) # the M(x)
  cluster_means <- aggregate(data_ext[,-which(names(data_ext)==cluster)], list(cluster = data_ext[,cluster]), mean)


  # for weighted cluster-means:
  if(!is.na(boundary_weights)){
    # weight the auxiliary information, but not the response:
    clustmeans.temp<- ddply(data_ext[,-c(which(names(data_ext) == all.vars(formula)[1]))],
                            .(cluster), function(x) colwise(weighted.mean, w = x[[boundary_weights]]) (x))
    clustmeans.temp<- clustmeans.temp[, -which(names(clustmeans.temp) == boundary_weights)] # remove column with boundary-weights, we don't need them anymore
    # merge unweighted cluster-means of response:
    cluster_means<- merge(clustmeans.temp, cluster_means[,c(cluster, all.vars(formula)[1])], by=cluster)
  }

  data_clust_s1 <- merge(cluster_means, cluster_weights, by=cluster)  #right hand column is M(x)


  data_clust_s1_groundindicator <- merge(data_clust_s1, aggregate(data[,c(cluster, phase_id[["phase.col"]])], list(data[,cluster]), unique)[-1], by=cluster, all.x=TRUE)



  # we need the s2 design matrix from the cluster-level dataset "data_clust"
  data_clust_s2 <- data_clust_s1_groundindicator[data_clust_s1_groundindicator[[phase_id[["phase.col"]]]]==phase_id[["terrgrid.id"]],] #THIS SHOULD HAVE PHASE_ID AT SOME POINT
  Yc_x <- data_clust_s2[, all.vars(formula)[1]]
  design_matrix.s1 <- data_clust_s1[, -1*c(which(names(data_clust_s1) %in% c(cluster,all.vars(formula)[1])), ncol(data_clust_s1))]
  design_matrix.s2 <- as.matrix(data_clust_s2[, -1*c(which(names(data_clust_s2) %in% c(cluster,all.vars(formula)[1],phase_id[["phase.col"]])), ncol(data_clust_s2)-1)])

  M_x.s1 <- data_clust_s1[, ncol(data_clust_s1)]
  M_x.s2 <- data_clust_s2[, ncol(data_clust_s2)-1] # the number of columns minus one location of M_x due to the merge order

  n1_clusters <- length(M_x.s1)
  n2_clusters <- length(Yc_x)

  As2inv <- solve( (t(design_matrix.s2) %*% (M_x.s2*design_matrix.s2)) / n2_clusters )  # (1/n2)SUMx[Mc(x)Zc(x)t(Zc(x))]
  beta_s2 <- t(As2inv %*% t(t(M_x.s2*Yc_x) %*% design_matrix.s2 / n2_clusters))

  Yc_x_hat_s2 <- design_matrix.s2 %*% t(beta_s2)
  Yc_x_hat_s1 <- as.matrix(design_matrix.s1) %*% t(beta_s2)
  Rc_x_hat <- Yc_x - Yc_x_hat_s2
  MR_square <- as.vector((M_x.s2^2)*(Rc_x_hat^2))  # eq [64] p. 24 Mandallaz Small Area technical report
  middle_term <- ((t(as.matrix(design_matrix.s2)) %*% ( MR_square*(as.matrix(design_matrix.s2)))) / n2_clusters^2) # eq [64] p. 24 Mandallaz Small Area technical report
  cov_beta_s2 <- As2inv %*% middle_term %*% As2inv  # eq [64] p. 24 Mandallaz Small Area technical report

  Z_bar_1 <- apply(design_matrix.s1, 2, weighted.mean, w=M_x.s1) # eq [65] p. 25 Mandallaz Small Area technical report

  #prepare matrix calculation for cov_Z_bar_1 by centering design matrix and scaling it with M_x ...eq [67] p. 25 Mandallaz Small Area technical report
  design_matrix.s1_centered <- t(apply(design_matrix.s1, 1, function(x){x-Z_bar_1}))
  design_matrix.s1_centered_weighted <- apply(design_matrix.s1_centered, 2, function(col){as.vector(M_x.s1 / mean(M_x.s1))*col})
  cov_Z_bar_1 <- t(design_matrix.s1_centered_weighted) %*% design_matrix.s1_centered_weighted / (n1_clusters*(n1_clusters-1))

  w_s2 <- (M_x.s2 / mean(M_x.s2))
  w_s1 <- (M_x.s1 / mean(M_x.s1))

  weighted_mean_Yc_x_hat <- sum(M_x.s1*Yc_x_hat_s1)/sum(M_x.s1)
  weighted_mean_Rc_x_hat <- sum(M_x.s2*Rc_x_hat)/sum(M_x.s2)

  estimate <- Z_bar_1 %*% t(beta_s2)
  ext_variance <- ((1/n1_clusters)*(1/(n1_clusters-1))*sum((w_s1^2)*(Yc_x_hat_s1 - weighted_mean_Yc_x_hat)^2)) + ((1/n2_clusters)*(1/(n2_clusters-1))*sum((w_s2^2)*(Rc_x_hat - weighted_mean_Rc_x_hat)^2))  #approximation in cluster case
  g_variance <- (t(Z_bar_1) %*% cov_beta_s2 %*% Z_bar_1) + (beta_s2 %*% cov_Z_bar_1 %*% t(beta_s2))  #eq [67] p. 25 Mandallaz Small Area technical report

  r.squared <- summary(lm(formula, data, y=TRUE))$r.squared # plot level R square, normally worse than R-square on cluster level


  ## ------- create outputs ------------------------------------------------- ##


  # summarize sample size info:
  samplesizes<- data.frame(cbind (n1_clusters, n2_clusters, n1, n2))
  colnames(samplesizes)<- c("n1_clust", "n2_clust", "n1", "n2")
  rownames(samplesizes)<- "plots"

  estimation<- data.frame(estimate=estimate, ext_variance=ext_variance, g_variance=g_variance,
                          n1=samplesizes$n1_clust, n2=samplesizes$n2_clust,
                          r.squared=r.squared)

  # ... to store inputs used:
  inputs<- list()
  inputs[["data"]]<- data
  inputs[["formula"]]<- formula
  inputs[["boundary_weights"]]<- boundary_weights
  inputs[["method"]]<- "non-exhaustive"
  inputs[["cluster"]]<- TRUE
  inputs[["exhaustive"]]<- FALSE


  # save warning-messages:
  warn.messages<- NA

  result<- list(input=inputs,
                estimation=estimation,
                samplesizes=samplesizes,
                coefficients=beta_s2,
                cov_coef=cov_beta_s2,
                Z_bar_1=Z_bar_1,
                cov_Z_bar_1G=cov_Z_bar_1,
                Rc_x_hat=Rc_x_hat,
                mean_Rc_x_hat=weighted_mean_Rc_x_hat,
                warn.messages=warn.messages)

  class(result)<- c("global", "twophase")

  return(result)

}

