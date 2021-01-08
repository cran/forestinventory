#' onephase
#'
#' \code{onephase} is used to calculate estimations exclusively based on
#' terrestrial observations of a forest inventory (i.e. the \emph{local densities}).
#' The estimation method is available for simple and cluster-sampling
#' and provides point estimates of the sample mean and their variances.
#'
#' @param formula an object of class "\code{\link[stats]{formula}}" that
#'                must be of the form \code{Y ~ 1}, where Y is the terrestrial response
#'                value of interest provided in \code{data}.
#'
#' @param data  a data frame or vector containing the response value Y.
#'              Specifications are given under 'Details'.
#'
#' @param phase_id an object of class "\code{\link[base]{list}}" containing two elements:
#'             \itemize{
#'                  \item \code{phase.col}: the column name in \code{data} that specifies the
#'                                          phase membership of each observation
#'                  \item \code{terrgrid.id}: the indicator identifying the the terrestrial
#'                                            (a.k.a. "ground truth") phase for that column
#'                                            (must be of type "\code{\link[base]{numeric}}")
#'                     }
#'        \strong{Note:} Only has to be specified if \code{data} is of class \code{data.frame}.
#'
#'
#' @param cluster Specifies the column name in \code{data}
#'                containing the cluster identification. Only used in case of
#'                cluster sampling.
#'
#' @param area (\emph{Optional}) an object of class "\code{\link[base]{list}}" containing two elements:
#'          \itemize{
#'                \item \code{sa.col}: the column name in \code{data} containing
#'                                      domain identification
#'                \item \code{areas}: vector of desired domains for which the estimation
#'                                    should be computed. If estimations for multiple domains should be computed,
#'                                    the domains have to be defined within a \code{character} vector using \code{c()}
#'                  }
#'
#'            Further details of the parameter-specifications are given under \emph{'Details'}.
#'
#' @details
#'
#'   \code{data} can either be a vector only containing the observations of the
#'               response variable Y,
#'               \emph{or} a data frame containing a column for the response variable and
#'               a column for the sample-grid indication that has to be further specified
#'               by argument \code{phase_id}.
#'               Additional \emph{optional} columns include a cluster identification in case of
#'               cluster sampling, as well as a column that specifies a domain (e.g. a forest district)
#'               the respective terrestrial observation falls into.
#'               The latter allows to compute onephase-estimations
#'               for multiple domains at a time (see \emph{'Examples'}).
#'
#'
#' @return \code{onephase} returns an object of class \code{"onephase"}.
#'
#' The functions \code{summary} and \code{confint} can be used to obtain a summary of the
#' estimation results (point estimations, variances and sample sizes) and the confidence intervals
#' for the respective point estimates.
#'
#' An object of class \code{"onephase"} returns a \code{list} of the following components:
#'
#'  \item{input}{a \code{list} containing the function inputs}
#'  \item{estimation}{a data frame containing the following components:
#'                   \itemize{
#'                    \item \code{area:} the domain (only present if argument \code{area} has been used)
#'                    \item \code{estimate:} the point estimate
#'                    \item \code{variance:} the variance of the point estimate
#'                    \item \code{n2:} the terrestrial sample size
#'                    }}
#'  \item{samplesizes}{a named numeric vector giving the terrestrial samplesize}
#'
#' @references Hill, A., Massey, A. F. (2021). \emph{The R Package forestinventory: Design-Based Global and Small Area Estimations for Multiphase Forest Inventories.} Journal of Statistical Software, 97(4), 1-40.
#' @references Mandallaz, D. (2007). \emph{Sampling techniques for forest inventories.} Chapter 4. CRC Press.

#'
#' @example examples/example_onephase_estimations.R
#'
#'
#' @export
onephase <- function(formula, data,
                     phase_id = list(phase.col=NA,terrgrid.id=NA),
                     cluster= NA, area= list(sa.col = NA, areas = NA)){

  # -------------------------------------------------------------------------- #
  # -------------------------------------------------------------------------- #
  # initial error-checking of mandatory input parameters:

    # ------------------------------------- #
    # 1) check if data is of type data.frame or vector:
    if (!is.data.frame(data) & !is.vector(data)) { stop("data must be a vector or a data.frame")}

    # --------------------------------------------------- #
    # 2) check if variables used in formula(s) exist in data:

    if (is.data.frame(data)) {

      if ( prod( all.vars(formula) %in% colnames(data) ) != 1){
        missing.var<- all.vars(formula) [which(all.vars(formula) %in% colnames(data) == FALSE)]
        stop(paste("Variable ", missing.var," used in regression-formula does not exist in data", sep = ""))
      }

    }

    # --------------------------------- #
    # 3) check for input-type "phase_id":

  if (all(!is.na(unlist(phase_id)))){

    # check if "phase.col" is character:
    if (!is.character(phase_id[["phase.col"]])) {
      stop("phase.col argument in input-parameter phase_id must be of type character")
    }
    # check if "phase.col" - name exists in data:
    if(phase_id[["phase.col"]] %in% colnames(data) == FALSE){ stop("specified phase column does not exist in data")}
    # check if terrgrid in "phase_id" does not exist in data:
    if(!all(phase_id[["terrgrid.id"]] %in% sapply(unique(data[phase_id[["phase.col"]]]), as.character))) {
      stop("the specified terrestrial grid indicator does not exist in data")
    }

  }



  # -------------------------------------------------------------------------- #
  # -------------------------------------------------------------------------- #
  # store fucntion-call:

  call<- match.call()

  # -------------------------------------------------------------------------- #
  # -------------------------------------------------------------------------- #
  # create one-phase estimation function as a closure:

  onephase_estimation<- function(){

    # -----------------------#

    if (all(!is.na(unlist(area)))){ # if a list of areas is given ...

      sa_col<- area[["sa.col"]]
      areas<-  area[["areas"]]

      # ... shrink dataset to current area:
      data<- data[data[[sa_col]] == areas[i], ]

    }
    # if area is not given, dataset remains as it is ...


    # -----------------------#
    # extract terrestrial dataset:

    # get name of response value:
    response<- all.vars(as.formula(formula))

    # retrieve phase.columnname and indicator of s1 grid id and terrestrial-grid id:

    if (all(!is.na(unlist(phase_id)))){ # if phase_id is specified ...

      phase.col<- phase_id[["phase.col"]]
      s2.ind<- phase_id[["terrgrid.id"]] # identifies the terrestrial sample (s2-sample)

    } else { # else assume that dataset only contains of terrestrial data

      # --> introduce an artificial terrestrial ground ID (==1):
      data<- data.frame(data, phase.col = rep(1, times= ifelse(is.data.frame(data), nrow(data), length(data))))
      if(!(response %in% colnames(data))){colnames(data)[1]<- response}
      phase.col<- "phase.col"
      s2.ind<- 1

    }

    # extract terrestrial data:
    data.terr<- data[data[[phase.col]] == s2.ind, ]


    # -----------------------#
    # NA-treatment:

    # rows to be deleted due to missing response value:
    deleted.s2<- !complete.cases(data.terr [,response]) # logical vector returning rows with missing response
    sum.NA_omitted<- sum(deleted.s2)

    # delete missing rows in entire dataset and produce message:
    if(sum.NA_omitted != 0) {
      data.terr<- data.terr[- which(deleted.s2),]
      message(paste(sum.NA_omitted," rows deleted due to missingness in the response value",sep = ""))
    }

    # -----------------------#


    # ------ non-cluster -----#

    if (is.na(cluster)){

      # sample size:
      n2<- nrow(data.terr)

      # check if any terrestrial data available:
      if (n2 == 0){

        warning("Estimation not possible, set to 'NA': Domain does not contain any terrestrial data", call. = F)

        estimation<- data.frame(estimate=NA, variance=NA, n2=n2)

        # summarize sample size info:
        samplesizes<- n2
        names(samplesizes)<- "n2"

      } else {

        # estimations:
        estimate<- mean(data.terr[[response]])
        variance<- (1/n2)*var(data.terr[[response]])

        if(n2 < 2){
          warning("Variance estimation not possible, set to 'NA': Domain contains only one cluster", call. = F)
          variance<- NA
        }

        # summarize sample size info:
        samplesizes<- n2
        names(samplesizes)<- "n2"

        estimation<- data.frame(estimate=estimate, variance=variance, n2=n2)

      }

    }

    # ------- cluster --------#

    if (!is.na(cluster)){

      # sample size:
      n2<- nrow(data.terr)

      if (n2 == 0){

        warning("Estimation not possible, set to 'NA': Domain does not contain any terrestrial data", call. = F)

        estimation<- data.frame(estimate=NA, variance=NA, n2=n2)

        # summarize sample size info:
        samplesizes<- data.frame(cbind ( 0,  n2))
        colnames(samplesizes)<- c("n2_clust", "n2")
        rownames(samplesizes)<- "plots"

      } else {

        # get on cluster level:
        cluster_weights<- aggregate(data.terr[[response]], list(cluster = data.terr[,cluster]), length) # the M(x) for sample s2
        names(cluster_weights)[2]<- "Mx"

        # plot sample size:
        n2_clusters<- nrow(cluster_weights)

        # local densities on cluster level:
        cluster_locdens<- aggregate(data.terr[[response]], list(cluster = data.terr[,cluster]), mean)
        names(cluster_locdens)[2]<- response

        # merge locdens and Mx:
        est.dat<- merge(cluster_locdens, cluster_weights, by=cluster)

        # estimations:
        estimate<- weighted.mean(est.dat[[response]], w = est.dat[["Mx"]])
        variance<- (1/(n2_clusters*(n2_clusters-1))) * sum( ((est.dat[["Mx"]] / mean(est.dat[["Mx"]]))^2) * ((est.dat[[response]] - estimate)^2) )

        if(n2_clusters < 2){
          warning("Variance estimation not possible, set to 'NA': Domain contains only one cluster", call. = F)
          variance<- NA
        }

        # summarize sample size info:
        samplesizes<- data.frame(cbind ( n2_clusters,  n2))
        colnames(samplesizes)<- c("n2_clust", "n2")
        rownames(samplesizes)<- "plots"

        estimation<- data.frame(estimate=estimate, variance=variance, n2=samplesizes$n2_clust)
      }

    }


    # ---- create outputs ------#

    # ... to store inputs used:
    inputs<- list()
    inputs[["data"]]<- data
    inputs[["method"]]<- "onephase"
    inputs[["cluster"]]<- !is.na(cluster)
    inputs[["call"]]<- call

    result<- list(input=inputs,
                  estimation=estimation,
                  samplesizes=samplesizes)

    return(result)

    } # end of onephase_estimation-function



  # -------------------------------------------------------------------------- #
  # -------------------------------------------------------------------------- #
  # run estimations:


  if (all(!is.na(unlist(area)))){ # if a list of areas is given ...


    # ------------------------------------ #

    # ... to store the estimates:
    areas<-  area[["areas"]]
    no.sa <- length(areas)
    results <- data.frame(area=rep(NA_character_, no.sa),
                         estimate=NA_real_, variance=NA_real_,
                         n2=NA_integer_)
    results[,"area"]<-as.character(results[,"area"])

    # ... to store the sample sizes:
    samplesizes<- list()

    # ------------------------------------ #

    # loop over areas:
    for (i in 1:no.sa) {

      result.temp<- onephase_estimation()

      # store estimation results
      results[i,names(results)[-1]] <- result.temp$estimation
      results[i,names(results)[1]]  <- areas[i] # add name of small area

      # store samplesize of small area:
      samplesizes[[areas[i]]]<- result.temp$samplesizes
    }


    result.temp[["input"]][["data"]]<- data
    result<- list(input=result.temp[["input"]],
                  estimation=results,
                  samplesizes=samplesizes)

    # ------------------------------------ #


  } else { # if no areas are given... --> just make onephase est. for entire dataset

  result.temp<- onephase_estimation()

  result<- list(input=result.temp[["input"]],
                estimation=result.temp$estimation,
                samplesizes=result.temp$samplesizes)

  }

  # -------------------------------------------------------------------------- #
  # -------------------------------------------------------------------------- #

    class(result)<- "onephase"

    return(result)


} # end of one-phase function


