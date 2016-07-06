

# help-function to calculate weighted cluster means:
boundaryweight_fct_3p_clust<- function(formula, data_ext, boundary_weights, cluster){

  # weight the auxiliary information, but NOT the response:
  clustmeans.temp<- ddply(data_ext[,-c(which(names(data_ext) == all.vars(formula)[1]))],
                          .(cluster), function(x) colwise(weighted.mean, w = x[[boundary_weights]]) (x))

  # calc. cluster-means of response:
  clustmeans.response<- aggregate(data_ext[,-which(names(data_ext)==cluster)],
                                  list(cluster = data_ext[,cluster]), mean)

  # merge unweighted cluster-means of response:
  cluster_means<- merge(clustmeans.temp, clustmeans.response[,c(cluster, all.vars(formula)[1])], by=cluster)
  cluster_means[,- which(colnames(cluster_means) %in% boundary_weights)] # remove boundary_weights- column
}
