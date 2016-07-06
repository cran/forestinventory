# This function helps the "two phase" cluster sampling estimator function when the first phase sample size goes to infinity.
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

global_exhaustive2p_cluster <- function(formula, data, phase_id, cluster, boundary_weights, exhaustive){

  n2<- sum(data[[phase_id[["phase.col"]]]] == phase_id[["terrgrid.id"]])
  n1<- Inf

  # extract data only for s2 sample:
  model.object_plot_level <- lm(formula, data=data[data[[phase_id[["phase.col"]]]]%in%phase_id[["terrgrid.id"]],], x=TRUE, y=TRUE)
  data_ext <- as.data.frame(cbind(model.object_plot_level$y, model.object_plot_level$x))
  names(data_ext)[1] <- all.vars(formula)[1]  #add response name
  data_ext[, cluster] <- data[,cluster][match(rownames(data_ext),rownames(data))] #add cluster variable to plot level s2
  if(!is.na(boundary_weights)){
    data_ext[,boundary_weights]<- data  [data[[phase_id[["phase.col"]]]]%in%phase_id[["terrgrid.id"]],boundary_weights]
  }

#   # alternative way to compute data_ext:
#   design_matrix.s1_plot_level <- design_matrix.s1_return(formula=formula, data=data[data[[phase_id[["phase.col"]]]]%in%phase_id[["terrgrid.id"]],])
#   data_ext <- data.frame(design_matrix.s1_plot_level)
#   data_ext[,all.vars(formula)[1]] <- data[data[[phase_id[["phase.col"]]]]%in%phase_id[["terrgrid.id"]],all.vars(formula)[1]]
#   data_ext[,cluster] <- data[data[[phase_id[["phase.col"]]]]%in%phase_id[["terrgrid.id"]],cluster]
#   if(!is.na(boundary_weights)){
#     data_ext[,boundary_weights]<- data[data[[phase_id[["phase.col"]]]]%in%phase_id[["terrgrid.id"]],boundary_weights]
#   }
#

  # calculate cluster means:
  cluster_weights <- aggregate(data_ext[,all.vars(formula)[1]], list(data_ext[,cluster]), length)
  cluster_means <- aggregate(data_ext[,-which(names(data_ext)==cluster)], list(cluster = data_ext[,cluster]), mean) #means of relevant variables at cluster level


  # for weighted cluster-means: --> boundary-weights have to be used when averaging the plots over clusters...
  if(!is.na(boundary_weights)){
    # library(plyr) #* load as dependency ...
    # weight the auxiliary information, but not the response:
    clustmeans.temp<- ddply(data_ext[,-which(names(data_ext) == all.vars(formula)[1])], .(cluster), function(x) colwise(weighted.mean, w = x[[boundary_weights]]) (x))
    clustmeans.temp<- clustmeans.temp[, -which(names(clustmeans.temp) == boundary_weights)] # remove column with boundary-weights, we don't need them anymore
    # merge unweighted cluster-means of response:
    cluster_means<- merge(clustmeans.temp, cluster_means[,c(cluster, all.vars(formula)[1])], by=cluster)
  }


  names(cluster_weights)[1] <- cluster
  names(cluster_means)[1] <- cluster
  data_clust_s2 <- merge(cluster_means, cluster_weights, by=cluster)

  Yc_x <- data_clust_s2[, all.vars(formula)[1]]
  design_matrix.s2 <- as.matrix(data_clust_s2[, -1*c(which(names(data_clust_s2) %in% c(cluster,all.vars(formula)[1])), ncol(data_clust_s2))])
  M_x.s2 <- data_clust_s2[, ncol(data_clust_s2)]  #this is clunky but the number of columns is the right hand one due to the merge order

  n1_clusters <- Inf
  n2_clusters <- length(Yc_x)

  As2inv <- solve( (t(design_matrix.s2) %*% (M_x.s2*design_matrix.s2)) / n2_clusters )
  beta_s2 <- As2inv %*% t(t(M_x.s2*Yc_x) %*% design_matrix.s2 / n2_clusters)

  Yc_x_hat <- design_matrix.s2 %*% beta_s2
  Rc_x_hat <- Yc_x - Yc_x_hat
  MR_square <- as.vector((M_x.s2^2)*(Rc_x_hat^2))
  middle_term <- ((t(as.matrix(design_matrix.s2)) %*% ( MR_square*(as.matrix(design_matrix.s2)))) / n2_clusters^2)
  cov_beta_s2 <- As2inv %*% middle_term %*% As2inv

  w <- (M_x.s2 / mean(M_x.s2))
  weighted_mean_Yc_x <- sum(M_x.s2*Yc_x)/sum(M_x.s2)
  weighted_mean_Rc_x_hat <- sum(M_x.s2*Rc_x_hat)/sum(M_x.s2)

  estimate <- exhaustive %*% beta_s2
  ext_variance <- (1/n2_clusters)*(1/(n2_clusters-1))*sum((w^2)*(Rc_x_hat - weighted_mean_Rc_x_hat)^2)  #approximation in cluster case
  g_variance <- (t(exhaustive) %*% cov_beta_s2 %*% exhaustive)


  ## ------- create outputs ------------------------------------------------- ##


  # summarize sample size info:
  samplesizes<- data.frame(cbind (n1_clusters, n2_clusters, n1, n2))
  colnames(samplesizes)<- c("n1_clust", "n2_clust", "n1", "n2")
  rownames(samplesizes)<- "plots"

  estimation<- data.frame(estimate=estimate, ext_variance=ext_variance, g_variance=g_variance,
                 n1=samplesizes$n1, n2=samplesizes$n2, n1_clust=samplesizes$n1_clust, n2_clust=samplesizes$n2_clust,
                 r.squared=summary(model.object_plot_level)$r.squared)

  # ... to store inputs used:
  inputs<- list()
  inputs[["data"]]<- data
  inputs[["formula"]]<- formula
  inputs[["boundary_weights"]]<- boundary_weights
  inputs[["method"]]<- "exhaustive"
  inputs[["cluster"]]<- TRUE
  inputs[["exhaustive"]]<- TRUE


  # save warning-messages:
  warn.messages<- NA

  result<- list(input=inputs,
                estimation=estimation,
                samplesizes=samplesizes,
                coefficients=t(beta_s2),
                cov_coef=cov_beta_s2,
                Rc_x_hat=Rc_x_hat,
                mean_Rc_x_hat=weighted_mean_Rc_x_hat,
                warn.messages=warn.messages)

  class(result)<- c("global", "twophase")

  return(result)


}






