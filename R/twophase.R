#' twophase
#'
#' \code{twophase} is used to calculate estimations based on double sampling under the
#' \emph{model-assisted Monte Carlo approach}. A \emph{first phase} of auxiliary information
#' (e.g. taken from remote sensing data) is used to generate model predictions based on multiple linear
#' regression  using the method of ordinary least squares. A subsample of the first phase comprises
#' the \emph{second phase} which contains terrestrial observations (i.e. the \emph{local densities}
#' of the ground truth) that is used to correct for bias in the design-based sense.
#' The estimation method is available for \emph{simple} and \emph{cluster sampling} and includes
#' the special case where the first phase is based on an \emph{exhaustive} sample (i.e. a census).
#' \emph{Small-area applications} are supported for synthetic estimation as well as two varieties
#' of bias-corrected estimators: the traditional small-area estimator and an asymptotically
#' equivalent version derived under Mandallaz' extended model approach.
#'
#' @param formula an object of class "\code{\link[stats]{formula}}" as would be used in the function \code{\link[stats]{lm}}
#'
#' @param data  a data frame containing all variables contained in \code{formula} and a column indexing
#'                  phase membership.  Additional columns designating small-area membership, cluster ID and
#'                  boundary weights should also be contained in the data frame if they are
#'                  requested in the function.
#'
#' @param phase_id an object of class "\code{\link[base]{list}}" containing two elements:
#'             \itemize{
#'                  \item \code{phase.col}: the column name in \code{data} that specifies the
#'                                          phase membership of each observation
#'                  \item \code{terrgrid.id}: the indicator identifying the terrestrial
#'                                            (a.k.a. "ground truth") phase for that column
#'                                             (must be of type "\code{\link[base]{numeric}}")
#'                     }
#'
#' @param cluster (\emph{Optional}) Specifies the column name in \code{data}
#'                containing the cluster ID. Only used in case of
#'                cluster sampling.
#'
#' @param small_area (\emph{Optional}) a list that if containing three elements:
#'             \itemize{
#'                  \item \code{sa.col}: the column name in \code{data} containing
#'                                       domain identification
#'                  \item \code{areas}: vector of desired small-area domain identifiers
#'                  \item \code{unbiased}: an object of type "\code{\link[base]{logical}}"
#'                                         that when FALSE designates that the estimator is allowed to be
#'                                         biased (i.e. the synthetic estimator) and when TRUE forces
#'                                         it to be design-unbiased. See \emph{'Details'}.
#'                     }
#'
#'             \strong{Note}: If \code{small_area} is left unchanged then \code{twophase} defaults to global estimation.
#'
#' @param boundary_weights (\emph{Optional}) Specifies the column name in \code{data}
#'                containing the weights for boundary adjustment.  See \emph{'Details'}
#'
#' @param exhaustive (\emph{Optional}) For global estimation, a vector of true auxiliary means corresponding to
#'                   an exhaustive first phase.
#'                   The vector must be input in the same order that \code{lm} processes a \code{formula} object
#'                   and include the intercept term.
#'                   For small area estimation, \code{exhaustive} is a \code{data.frame} containing column names
#'                   (\code{\link[base]{colnames}}) for every variable appearing in the parameter \code{formula} including
#'                   the variable "Intercept".Rownames (\code{\link[base]{row.names}}) have to be used and must correspond
#'                   to the names of the small areas. See \emph{'Details'}.
#'
#' @param progressbar (\emph{Optional}) an object a type "\code{\link[base]{logical}}" that when TRUE prints
#'            the progress of the calculation in the console (recommended for large amount of small areas).  Defaults to FALSE.
#'
#' @param psmall (\emph{Optional}) an object a type "\code{\link[base]{logical}}" used for small area estimations
#'            that only works when \code{unbiased} in the parameter \code{small_area} is set to TRUE. See \emph{'Details'}.
#'
#'
#' @details If estimations for multiple small-area domains should be computed, the domains have to be
#'          defined within a \code{character} vector using \code{c()}. Using \code{small_area(..., unbiased=FALSE)}
#'          calculates design-based estimates with the synthetic estimator and may be design-biased if
#'          the model is biased in that small area.  The default, \code{small_area(..., unbiased=TRUE)}, allows for a residual
#'          correction by one of two asymptotically equivalent methods to create design-unbiased estimates:
#'          \itemize{
#'              \item Mandallaz' extended model approach calculates the residual correction by extending the
#'                    model formula with an indicator variable in the small area.  It is the default method
#'                    \code{psmall}=FALSE.
#'              \item the traditional small area estimator calculates the residual correction by taking the
#'                    synthetic estimator and adding the mean residual observed in the small area.  It is activated
#'                    when \code{psmall}=TRUE.
#'                  }
#'
#'          Missing values (\code{NA}) in the auxiliary variables (i.e. at least one auxiliary variable cannot be observed at
#'          an inventory location) are automatically removed from the dataset \emph{before} the estimations are computed.
#'          Note that missingness in the auxiliary variables is only allowed if we assume that they are \emph{missing at random},
#'          since the unbiasedness of the estimates is based on the sampling design.
#'
#'          The boundary weight adjustment is pertinent for auxiliary information derived from remote sensing and
#'          is equal to the percentage of forested area (e.g. as defined by a forest mask) in the interpretation area.
#'
#'          Exhaustive estimation refers to when the true means of certain auxiliary variables are known
#'          and an exhaustive first phase (i.e. a census).  For global estimation, the vector must be input
#'          in the same order that \code{lm} processes a \code{formula} object including the intercept term whose
#'          true mean will always be one.  For small area estimation, \code{exhaustive} is a \code{data.frame} containing column names for every variable appearing in
#'          the parameter \code{formula} including the variable "Intercept".  The observations of the data.frame
#'          must represent the true auxiliary means in the same order as was presented in \code{areas} from the
#'          parameter \code{small_area}.  See \emph{'Examples'}.
#'
#'
#' @return \code{twophase} returns an object of class \code{"twophase"}.
#'
#' An object of class \code{"twophase"} returns a \code{list} of the following components:
#'
#'  \item{input}{a \code{list} containing the function's inputs}
#'  \item{estimation}{a data frame containing the following components:
#'                   \itemize{
#'                    \item \code{area:} the domain (only present if argument \code{areas} has been used)
#'                    \item \code{estimate:} the point estimate
#'                    \item \code{ext_variance:} the external variance of the point estimate that doesn't account for
#'                                               fitting the model from the current inventory
#'                    \item \code{g_variance:} the internal (g-weight) variance that accounts for
#'                                               fitting the model from the current inventory
#'                    \item \code{n1} the first phase sample size of plots
#'                    \item \code{n2} the second phase (i.e. terrestrial) sample size of plots
#'                    \item \code{n1G} the first phase sample size in the small area
#'                    \item \code{n2G} the second phase (i.e. terrestrial) sample size in the small area
#'                    \item \code{r.squared} the R squared of the linear model
#'                    }}
#'  \item{samplesizes}{a \code{\link[base]{data.frame}} summarizing all samplesizes: in case of cluster sampling both,
#'                     the number of individual plots and the number of clusters is reported.}
#'  \item{coefficients}{the linear model coefficients}
#'  \item{cov_coef}{the design-based covariance matrix of the model coefficients}
#'  \item{Z_bar_1G}{the estimated auxiliary means of \code{formula} based on the first phase.
#'                  If the first phase is exhaustive, these are the true auxiliary means specified in the input-argument \code{exhaustive}.}
#'  \item{cov_Z_bar_1G}{the covariance matrix of \code{Z_bar_1G}}
#'  \item{Rc_x_hat_G}{the small-area residuals at either the plot level or cluster level depending on the call}
#'  \item{Rc_x_hat}{the residuals at either the plot level or cluster level depending on the call}
#'  \item{Yx_s2G}{the local densities in the small area}
#'  \item{Mx_s2G}{the cluster weights in the small area}
#'  \item{mean_Rc_x_hat_G}{the mean residual (weighted mean in the case of cluster sampling) in the small area}
#'  \item{mean_Rc_x_hat}{the mean residual (weighted mean in the case of cluster sampling)}
#'  \item{warn.messages}{logical indicating if warning messages were issued}
#'
#' @note
#' In the special case of cluster sampling, the reported sample sizes in \code{estimation} are the number of clusters.
#' The \code{samplesize}-object also provides the respective number of single plot units for cluster sampling.
#' The reported \code{r.squared} describe the model fit of the applied linear regression
#' model (i.e. on \emph{plot-level}, not on \emph{cluster level}).
#'
#' @references Hill, A., Massey, A. F. (2021). \emph{The R Package forestinventory: Design-Based Global and Small Area Estimations for Multiphase Forest Inventories.} Journal of Statistical Software, 97(4), 1-40.
#' @references Mandallaz, D. (2007). \emph{Sampling techniques for forest inventories.} Chapter 4. CRC Press.
#' @references Mandallaz, D. (2013). \emph{Design-based properties of some small-area estimators in forest inventory with two-phase sampling.} Can. J. For. Res. 43: 441-449
#' @references Mandallaz, D. and Hill, A. and Massey, A. (2016). \emph{Design-based properties of some small-area estimators in forest inventory with two-phase sampling.} ETH Zurich, Department of Environmental Systems Science,Tech. rep. Available from \url{http://e-collection.library.ethz.ch}.
#'
#' @example examples/example_twophase_estimations_long.R
#'
#' @import plyr
#' @import stats
#' @import utils
#' @export

# This is the master two-phase function that calls helper functions for the following estimates
#     -global (exhaustive/nonexhaustive)
#     -small area (exhaustive/nonexhaustive)


# -----------------------------------------------------------------------------#
# SUPER-FUNCTION for two-phase estimations: (Documentation draft)

# Function does...
# functions takes all arguments, performs error-checking of input parameters
# and calls respective estimator functions according to function-input

# MANDATORY INPUT:
# formula:            a formula object specifying the regression-formula as for lm-function
# data:               dataset containing the inventory-data (large sample containing auxiliary information
#                     for each plot, reponse available for small (terrestrial) sample, otherwise set to NA
# phase_id:           list of characters, specifying
#                               1) the columnname where information of sample-grid membership is stored
#                               2) the indicator of the second phase grid (terrestrial sample)
#
# OPTIONAL INPUT:
# cluster:            character indicating the columnname where cluster-membership for each plot is stored
# boundary_weights:   character indicating the columnname where boundary-weights for each plot is stored
#                     to calculate accurate weighted means of auxiliary information
#
# small_area:         list of characters, specifying
#                               1) the columnname where information of small area membership is stored
#                               2) vector with small area names for which sae-estimations are executed
#                               3) a extended logical value (TRUE or FALSE or BOTH):
#
#                                  - if TRUE, we expand the design-matrix and the model-formula by the
#                                             indicator-variable to enforce zero-mean residuals over F and G.
#                                             Introducing the indicator variable requires the presents of terrestrial data in the small area G (rule of thumb: at least 5).
#                                             If the small area G does not contain any s2-points, the results will be signed "NA", since the calculation
#                                             of the point estimate and the design-based and external variance is not possible.
#                                             We refer to this estimator as not-synthetic but design-"unbiased",
#                                             although Mandallaz calls it synthetic due to the vanishing residual-terms in the formulas).
#
#                                  - if FALSE, we do NOT introduce the indicator variable to enforce zero-mean residuals over F and G.
#                                              In this case, the residuals in G are not insured to be zero in general.
#                                              This is actually the only chance to get an estimation if small area G does not contain any terrestrial data, but
#                                              NO external variance can be given in this case. The method is also be advisable if small area G does contain only very few
#                                              terrestrial data, since fitting an additional intercept for G is only based on few data points.
#                                              Since the design-based point estimate and design-based variance is hence calculated using the globally-calculated
#                                              regression coefficients,the estimations are potentially biased and the design-based variance will tend to be dramatically
#                                              smaller than their "unbiased" counterpart (we have to believe that the global model will also be appropriate for the small area).
#
#                                  - if BOTH, the decision whether to apply the unbiased small area estimator or the potentially biased sysnthetic estimator is made for each
#                                             small area individually based on a user defined threshold (thresh.n2g) of minimal number of n2g within a small area
#
#
# exhaustive:         vector containing the true mean for each auxiliary information
#
#
# progressbar:      default=FALSE; if TRUE, progress of estimation is displayed in console (recommended for large amount of small areas)


# -----------------------------------------------------------------------------#
# FUNCTION STARTS HERE:

twophase <- function(formula, data, phase_id, cluster=NA,
                     small_area= list(sa.col = NA, areas = NA, unbiased=TRUE),
                     boundary_weights=NA, exhaustive=NA, progressbar=FALSE, psmall=FALSE){


  # -------------------------------------------------------------------------- #
  # -------------------------------------------------------------------------- #
  # # source everything that we need "globally":

  # source("design_matrix.s1_return.R")
  # source("errorchecking.R")
  # source("print_methods.R")
  # source("summary_methods.R")
  # source("confint_methods.R")
  # ...

  call<- match.call()

  # -------------------------------------------------------------------------- #
  # -------------------------------------------------------------------------- #
  # initial error-checking of mandatory input parameters:

  check.mandatoryInputs(formula, data, phase_id)

  # -------------------------------------------------------------------------- #
  # -------------------------------------------------------------------------- #
  # Checking the nesting of the sample-design:
  # --> each s2-point muss have the complete set of auxvars (s1-info) available

  nest.violation<- sum(!complete.cases(data [ data[[phase_id[["phase.col"]]]] == phase_id[["terrgrid.id"]] , which(colnames(data) %in% all.vars(formula)[-1]) ]))

  if(nest.violation > 0){
    warning(paste("Sample design not nested: for",nest.violation,"terrestrial plots at least one auxiliary parameter is missing"))
  }


  # -------------------------------------------------------------------------- #
  # -------------------------------------------------------------------------- #
  # NA-treatment:

  # rows to be deleted due to missingness of expl. variables in s1 (i.e. set of auxiliary informations):
  deleted.s1<- !complete.cases(data [, which(colnames(data) %in% all.vars(formula)[-1])]) # logical vector returning rows with missing entries
  sum.NA_omitted<- sum(deleted.s1)

  # delete missing rows in entire dataset and produce message:
  if(sum.NA_omitted != 0) {
    data<- data[- which(deleted.s1),]
    warning(paste(sum.NA_omitted," rows deleted due to missingness in the set of auxiliary parameters (",
                  nest.violation," terrestrial plots affected by deletion)",sep = ""))
  }

  # check if there is an s2 point without response value, but with all auxiliary info:
  change.s2.to.s1<- data[[phase_id[["phase.col"]]]] == phase_id[["terrgrid.id"]] & !complete.cases(data[,all.vars(formula)[1]]) & complete.cases(data[,all.vars(formula)[-1]])
  sum.NA_change.s2.to.s1<- sum(change.s2.to.s1)
  # ... and in that case change missing reponse information for s2-grid to s1-grid:
  if(sum.NA_change.s2.to.s1 > 0) {
   s1.id<- unique(data [[ phase_id[["phase.col"]] ]])[!(unique(data [[ phase_id[["phase.col"]] ]]) %in% phase_id[["terrgrid.id"]])]
   data[[phase_id[["phase.col"]]]] [change.s2.to.s1]<- s1.id
   warning(paste("Changed the phase_id for ",sum.NA_change.s2.to.s1," rows to the first phase (s1) due to missing value for the response variable in the second phase (s2)",sep = ""))
  }

  # check if every remaining terrestrial plot (s2) has a response-value assigned:
  deleted.s2 <- data[[phase_id[["phase.col"]]]] == phase_id[["terrgrid.id"]] & !complete.cases(data[,all.vars(formula)[1]]) & !complete.cases(data[,all.vars(formula)[-1]])
  sum.deleted.s2<- sum(deleted.s2)
  # ... and if not, delete these rows from the dataset:
  if(sum.deleted.s2 != 0) {
    data<- data[- which(deleted.s2),]
    warning(paste(sum.deleted.s2," rows deleted due to missing value for the response variable AND at least one of the auxiliary variables", sep = ""))
  }


  # -------------------------------------------------------------------------- #
  # -------------------------------------------------------------------------- #
  # NA-treatment:

  # # rows to be deleted due to missingness in s1 (i.e. set of auxiliary informations):
  # deleted.s1<- !complete.cases(data [, which(colnames(data) %in% all.vars(formula)[-1])]) # logical vector returning rows with missing entries
  # sum.NA_omitted<- sum(deleted.s1)
  #
  # # delete missing rows in entire dataset and produce message:
  # if(sum.NA_omitted != 0) {
  #   data<- data[- which(deleted.s1),]
  #   message(paste(sum.NA_omitted," rows deleted due to missingness in the auxiliary parameters",sep = ""))
  #   }
  #
  # # check if every entry marked as a terrestrial plot (s2) has a response-value assigned:
  # deleted.s2 <- data[[phase_id[["phase.col"]]]] == phase_id[["terrgrid.id"]] & !complete.cases(data[,all.vars(formula)[1]])
  # sum.deleted.s2<- sum(deleted.s2)
  #
  # # change missing reponse information for s2-grid to s1-grid and produce message:
  # if(sum.deleted.s2 != 0) {
  #   data<- data[- which(deleted.s2),]
  #   message(paste(sum.deleted.s2," rows deleted due to missing value for the response variable", sep = ""))
  #   }


  # -------------------------------------------------------------------------- #
  # -------------------------------------------------------------------------- #

  #---------------------------------#
  # call for non-cluster functions: #
  # --------------------------------#

  if(is.na(cluster)) {

  # ---------------------------------------------------------------------#
  # check if --> "global 2-phase non-exhaustive" is required and apply estimator:

    if(is.na(small_area[["sa.col"]]) & all(is.na(exhaustive))) {

      # --- error checking -- :
      if(!is.na(boundary_weights)){ check.boundary_weightsInput(data, boundary_weights)}

      # -- # source global_nonexhaustive2p.R - function --:
      # source("global_nonexhaustive2p.R")

      # -- call function --:
      result<- global_nonexhaustive2p(formula, data, phase_id, boundary_weights)

    }

  # ---------------------------------------------------------------------#
  # check if --> "global 2-phase exhaustive" is required and apply estimator:

    if(is.na(small_area[["sa.col"]]) & all(!is.na(exhaustive))) {

      # --- error checking -- :
      check.exhaustive(formula, boundary_weights, exhaustive, data)

      # -- # source global_exhaustive2p.R - function --:
      # source("global_exhaustive2p.R")

      # -- call function -- :
      result<- global_exhaustive2p(formula, data, exhaustive)

    }

    # ---------------------------------------------------------------------#
    # check if --> "SMALL AREA 2-phase non-exhaustive" is required and apply estimator:

    if(!is.na(small_area[["sa.col"]]) & all(is.na(exhaustive))) {

      # --- error checking -- :
      check.smallareaInput(data, small_area)
      if(!is.na(boundary_weights)){ check.boundary_weightsInput(data, boundary_weights)}

      # -- # source global_exhaustive2p.R - function --:
      # source("small_area_looper.R")

      # -- call function -- :
      if(!psmall){
        result <- small_area_looper(formula, data, phase_id, cluster, small_area, boundary_weights, exhaustive, progressbar, psmall)
      } else {
        result<- psmall_fct(formula, data, phase_id, cluster, small_area, boundary_weights, exhaustive, progressbar, psmall)
      }


    }

    # ---------------------------------------------------------------------#
    # check if --> "SMALL AREA 2-phase exhaustive" is required and apply estimator:

    if(!is.na(small_area[["sa.col"]]) & all(!is.na(exhaustive))) {

      # --- error checking -- :
      check.smallareaInput(data, small_area)
      check.exhaustive(formula, boundary_weights, exhaustive, data)

      # -- # source global_exhaustive2p.R - function --:
      # source("small_area_looper.R")

      # -- call function -- :
      if(!psmall){
        result <- small_area_looper(formula, data, phase_id, cluster, small_area, boundary_weights, exhaustive, progressbar, psmall)
      } else {
        result<- psmall_fct(formula, data, phase_id, cluster, small_area, boundary_weights, exhaustive, progressbar, psmall)
      }

    }

  } # end of non-cluster function calls


  # -------------------------------------------------------------------------- #
  # -------------------------------------------------------------------------- #

  #---------------------------------#
  # call for cluster functions:     #
  # --------------------------------#


  if(!is.na(cluster)) {

    # --- error checking -- :
    check.clusterInput(data, cluster)

    # --- rename cluster colname to "cluster" --- #
    colnames(data)[which(colnames(data) %in% cluster)]<- "cluster"
    cluster.orig<- cluster
    cluster<- "cluster"

    if(!all(!is.na(data[["cluster"]]))){print(paste("WARNING: NAs removed from ",cluster))} #warning for NAs in cluster id
    data <- data[!is.na(data[["cluster"]]),] #strip NA clusterIDS


    # ---------------------------------------------------------------------#
    # check if --> "global 2-phase non-exhaustive_cluster" is required and apply estimator:

    if(is.na(small_area[["sa.col"]]) & all(is.na(exhaustive))) {

      # --- error checking -- :
      if(!is.na(boundary_weights)){ check.boundary_weightsInput(data, boundary_weights)}

      # -- # source ....R - function --:
      # source("global_nonexhaustive2p_cluster.R")

      # -- call function -- :
      result<- global_nonexhaustive2p_cluster(formula, data, phase_id, cluster, boundary_weights)

    }

    # ---------------------------------------------------------------------#
    # check if --> "global 2-phase exhaustive_cluster" is required and apply estimator:

    if(is.na(small_area[["sa.col"]]) & all(!is.na(exhaustive))) {

      # --- error checking -- :
      check.exhaustive.cluster(formula, boundary_weights, exhaustive, data)

      # -- # source ....R - function --:
      # source("global_exhaustive2p_cluster.R")

      # -- call function -- :
      result<- global_exhaustive2p_cluster(formula, data, phase_id, cluster, boundary_weights, exhaustive)

    }

    # ---------------------------------------------------------------------#
    # check if --> "SMALL AREA 2-phase non-exhaustive cluster" is required and apply estimator:

    if(!is.na(small_area[["sa.col"]]) & all(is.na(exhaustive))) {

      # --- error checking -- :
      check.smallareaInput(data, small_area)
      if(!is.na(boundary_weights)){ check.boundary_weightsInput(data, boundary_weights)}

      # -- # source ....p.R - function --:
      # source("small_area_looper.R")

      # -- call function -- :
      if(!psmall){
        result <- small_area_looper(formula, data, phase_id, cluster, small_area, boundary_weights, exhaustive, progressbar, psmall)
      } else {
        result<- psmall_fct(formula, data, phase_id, cluster, small_area, boundary_weights, exhaustive, progressbar, psmall)
      }

    }


    # ---------------------------------------------------------------------#
    # check if --> "SMALL AREA 2-phase exhaustive cluster" is required and apply estimator:

    if(!is.na(small_area[["sa.col"]]) & all(!is.na(exhaustive))) {

      # --- error checking -- :
      check.smallareaInput(data, small_area)
      check.exhaustive.cluster(formula, boundary_weights, exhaustive, data)

      # -- # source ....p.R - function --:
      # source("small_area_looper.R")

      # -- call function -- :
      if(!psmall){
        result <- small_area_looper(formula, data, phase_id, cluster, small_area, boundary_weights, exhaustive, progressbar, psmall)
      } else {
        result<- psmall_fct(formula, data, phase_id, cluster, small_area, boundary_weights, exhaustive, progressbar, psmall)
      }

    }


    # ---------------------------------------------------------------------#

    # rename cluster colname to original #
    colnames(result$input$data)[which(colnames(result$input$data) %in% "cluster")]<- cluster.orig


  } # end of cluster function calls


  # -------------------------------------------------------------------------- #
  # -------------------------------------------------------------------------- #

  # add function call to returned-list:
  result[["input"]]<- c(result[["input"]], call=call)


  result

} # end of Super-Function


  # -------------------------------------------------------------------------- #
  # -------------------------------------------------------------------------- #
