#' threephase
#'
#' \code{threephase} is used to calculate estimations based on triple sampling under the
#' \emph{model-assisted Monte Carlo approach}. A \emph{first phase} of auxiliary information
#' (e.g. taken from remote sensing data) is used to generate model predictions based on multiple linear
#' regression using the method of ordinary least squares. A subsample of the first phase comprises
#' a \emph{second phase} which contains further auxiliary information that produces another set of model predictions.
#' A further subsample produces a \emph{third final phase} based on terrestrial observations
#' (i.e. the \emph{local densities} of the ground truth) and is used to correct for bias in the design-based sense.
#' The estimation method is available for \emph{simple} and \emph{cluster sampling} and includes
#' the special case where the first phase is based on an \emph{exhaustive} sample (i.e. a census).
#' \emph{Small-area applications} are supported for synthetic estimation as well as two varieties
#' of bias-corrected estimators: the traditional small-area estimator and an asymptotically
#' equivalent version derived under Mandallaz's extended model approach.
#'
#' @param formula.s0 an object of class "\code{\link[stats]{formula}}" as would be used in the function \code{\link[stats]{lm}}
#'                   that contains a reduced set of auxiliary variables available for all first phase plots
#'
#' @param formula.s1 an object of class "\code{\link[stats]{formula}}" as would be used in the function \code{\link[stats]{lm}}
#'                   that contains the predictors from \code{formula.s0} as well as further ancilliary predictors available
#'                   for all second phase plots (i.e. \code{formula.s0} is \strong{nested} in \code{formula.s1})
#'
#' @param data  a data frame containing all variables contained in \code{formula} and a column indexing
#'                  phase membership.  Additional columns designating small-area membership, cluster ID and
#'                  boundary weights should also be contained in the data frame if they are
#'                  requested in the function.
#'
#' @param phase_id an object of class "\code{\link[base]{list}}" containing three elements:
#'             \itemize{
#'                  \item \code{phase.col}: the column name in \code{data} that specifies the
#'                                          phase membership of each observation
#'                  \item \code{s1.id}: the indicator identifying the "second phase only" plots
#'                                            for that column  (must be of type "\code{\link[base]{numeric}}")
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
#' @details \code{s1.id} identifies "second phase only" plots because the terrestrial phase is
#'          known to be part of the second phase by the construction of the subsampling.
#'
#'          If estimations for multiple small-area domains should be computed, the domains have to be
#'          defined within a \code{character} vector using \code{c()}. Using \code{small_area(..., unbiased=FALSE)}
#'          calculates design-based estimates with the synthetic estimator and may be design-biased if
#'          the model is biased in that small area.  The default, \code{small_area(..., unbiased=TRUE)}, allows for a residual
#'          correction by one of two asympototically equivalent methods to create design-unbiased estimates:
#'          \itemize{
#'              \item Mandallaz's extended model approach calculates the residual correction by extending the
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
#' @return \code{threephase} returns an object of class \code{"threephase"}.
#'
#' An object of class \code{"threephase"} returns a \code{list} of the following components:
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
#'                    \item \code{n0} the first phase sample size of plots
#'                    \item \code{n1} the second phase sample size of plots
#'                    \item \code{n2} the third phase (i.e. terrestrial) sample size of plots
#'                    \item \code{n0G} the first phase sample size in the small area
#'                    \item \code{n1G} the second phase sample size in the small area
#'                    \item \code{n2G} the third phase (i.e. terrestrial) sample size in the small area
#'                    \item \code{r.squared_reduced} the R-squared of the linear model based on \code{formula.s0} (i.e. the reduced model)
#'                    \item \code{r.squared_full} the R-squared of the linear model based on \code{formula.s1} (i.e. the full model)
#'                    }}
#'  \item{samplesizes}{a \code{\link[base]{data.frame}} summarizing all samplesizes: in case of cluster sampling both,
#'                     the number of individual plots and the number of clusters is reported.}
#'  \item{coefficients}{the coefficients of the two linear models:
#'                   \itemize{
#'                     \item \code{alpha:} the reduced model coefficients
#'                    \item \code{beta:} the full model coefficients
#'                    }}
#'  \item{cov_alpha_s2}{the design-based covariance matrix of the reduced model coefficients}
#'  \item{cov_beta_s2}{the design-based covariance matrix of the full model coefficients}
#'  \item{Z_bar_1_s0}{the estimated auxiliary means of \code{formula.s0} based on the first phase.
#'                    If the first phase is exhaustive, these are the true auxiliary means specified in the input-argument \code{exhaustive}.}
#'  \item{Z1_bar_s1}{the estimated auxiliary means of \code{formula.s0} based on the second phase}
#'  \item{Z_bar_s1}{the estimated auxiliary means of \code{formula.s1} based on the second phase}
#'  \item{cov_Z_bar_1_s0}{the covariance matrix for \code{Z_bar_1_s0}}
#'  \item{resid_reduced}{the reduced model residuals at either the plot level or cluster level depending on the call}
#'  \item{resid_full}{the full model residuals at either the plot level or cluster level depending on the call}
#'  \item{warn.messages}{logical indicating if warning messages were issued}
#'
#' @note
#' In the special case of cluster sampling, the reported sample sizes in \code{estimation} are the number of clusters.
#' The \code{samplesize}-object also provides the respective number of single plot units for cluster sampling.
#' The reported \code{r.squared_reduced} and \code{r.squared_full} describe the model fit of the applied linear regression
#' models (i.e. on \emph{plot-level}, not on \emph{cluster level}).
#'
#' @references Mandallaz, D., Breschan, J., & Hill, A. (2013). \emph{New regression estimators in forest inventories
#' with two-phase sampling and partially exhaustive information: a design-based monte carlo approach
#' with applications to small-area estimation.} Canadian Journal of Forest Research, 43(11), 1023-1031.
#' @references Mandallaz, D. (2014). \emph{A three-phase sampling extension of the generalized regression estimator with partially exhaustive information.} Can. J. For. Res. 44: 383-388
#' @references Massey, A. and Mandallaz, D. and Lanz, A. (2014). \emph{Integrating remote sensing and past inventory data under the new annual design of the Swiss National Forest Inventory using three-phase design-based regression estimation.} Can. J. For. Res. 44(10): 1177-1186
#' @references Mandallaz, D. (2013). \emph{Regression estimators in forest inventories with three-phase sampling and two multivariate components of auxiliary information.} ETH Zurich, Department of Environmental Systems Science,Tech. rep. Available from \url{http://e-collection.library.ethz.ch}.
#'
#' @example examples/example_threephase_estimations_long.R
#'
#' @import plyr
#' @import stats
#' @import utils
#' @export

# This is the master three-phase function that calls helper functions for the following estimates
#     -global (exhaustive/nonexhaustive)
#     -small area (exhaustive/nonexhaustive)


# -----------------------------------------------------------------------------#
# SUPER-FUNCTION for three-phase estimations: (Documentation draft)

# Function does...
# functions takes all arguments, performs error-checking of input parameters
# and calls respective estimator functions according to function-input

# MANDATORY INPUT:
# formula:            a formula object specifying the regression-formula as for lm-function
# data:               dataset containing the inventory-data (large sample containing auxiliary information
#                     for each plot, reponse available for small (terrestrial) sample, otherwise set to NA
# phase_id:           list of characters, specifying
#                               1) the columnname where information of sample-grid membership is stored (phase.col)
#                               2) the indicator value of the grid containing the all auxiliary variables (s2grid.id)
#                               3) the indicator value of the grid containing the ground truth (terrgrid.id)
#
# OPTIONAL INPUT:
# cluster:            character indicating the columnname where cluster-membership for each plot is stored
# boundary_weights:   character indicating the columnname where boundary-weights for each plot is stored
#                     to calculate accurate weighted means of auxiliary information
#
# small_area:         list of characters, specifying
#                               1) the columnname where information of small area membership is stored
#                               2) vector with small area names for which sae-estimations are executed
#                               3) a logical value (TRUE or FALSE):
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
# exhaustive:         vector containing the true mean for each first phase auxiliary variable


# -----------------------------------------------------------------------------#
# FUNCTION STARTS HERE:

threephase <- function(formula.s0, formula.s1, data, phase_id, cluster=NA,
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
  # source("convert_coefs_table_global3p.R")
  # ...

  # -------------------------------------------------------------------------- #
  # -------------------------------------------------------------------------- #

  # save function call for output:
  call<- match.call()

  # -------------------------------------------------------------------------- #
  # -------------------------------------------------------------------------- #
  # initial error-checking of mandatory input parameters:

  check.mandatoryInputs3p(formula.s0, formula.s1, data=data, phase_id=phase_id)


  # -------------------------------------------------------------------------- #
  # -------------------------------------------------------------------------- #
  # Checking the nesting of the sample-design:
  # --> each s2-point muss have the complete set of s1-auxvars (s1-info) available
  # --> each s2-point muss have the complete set of s0-auxvars (s0-info) available
  # --> each s1-point muss have the complete set of s0-auxvars (s0-info) available

  # test 1: s2 c s1 ?
  s2_in_s1.nest.violation<- sum(is.na(data [ data[[phase_id[["phase.col"]]]] == phase_id[["terrgrid.id"]] , which(colnames(data) %in% all.vars(formula.s1)[-1])]))

  # test 1: s2 c s0 ?
  s2_in_s0.nest.violation<- sum(is.na(data [ data[[phase_id[["phase.col"]]]] == phase_id[["terrgrid.id"]] , which(colnames(data) %in% all.vars(formula.s0)[-1])]))

  # test 1: s1 c s0 ?
  s1_in_s0.nest.violation<- sum(is.na(data [ data[[phase_id[["phase.col"]]]] == phase_id[["s1.id"]] , which(colnames(data) %in% all.vars(formula.s0)[-1])]))


  if(s2_in_s1.nest.violation > 0){ # read: "s2 with no s1-info"
    warning(paste("Sample design not nested: for",s2_in_s1.nest.violation,"terrestrial plots at least one auxiliary parameter of the second phase (s1) is missing"))
  }

  if(s2_in_s0.nest.violation > 0){
    warning(paste("Sample design not nested: for",s2_in_s0.nest.violation,"terrestrial plots at least one auxiliary parameter of the first phase (s0) is missing"))
  }

  if(s1_in_s0.nest.violation > 0){
    warning(paste("Sample design not nested: for",s1_in_s0.nest.violation,"second phase (s1) plots at least one auxiliary parameter of the first phase (s0) is missing"))
  }


  # -------------------------------------------------------------------------- #
  # -------------------------------------------------------------------------- #
  # NA-treatment:

  # rows to be deleted due to missing auxiliary information in s0: (s0 in indicated by s0-id AND s1-id in our datasets!)
  deleted.s0<- !complete.cases(data [ , which(colnames(data) %in% all.vars(formula.s0)[-1])]) # logical vector returning rows with missing entries
  sum.NA_omitted.s0<- sum(deleted.s0)

  # delete missing rows in s0 of entire dataset and produce message:
  if(sum.NA_omitted.s0 != 0) {
    data<- data[- which(deleted.s0),]
    m0<- message(paste(sum.NA_omitted.s0," rows deleted due to missingness in the set of auxiliary parameters for the first phase (s0) (",
                       s2_in_s0.nest.violation," terrestrial plots affected by deletion)",sep = ""))
  }


  # get NA-auxvars-entries for s1-phase (indicated by 's1.id' AND 'terrgrid.id' and turn them into s0-plots (i.e. change phase id).
  # If they were s2, we should also delete the terrestrial information for safety and clarity:
  change.s1.to.s0<- !complete.cases(data [, which(colnames(data) %in% all.vars(formula.s1)[-1])]) & data[[phase_id[["phase.col"]]]] %in% c(phase_id[["s1.id"]] , phase_id[["terrgrid.id"]])
  sum.NA_change.s1.to.s0<- sum(change.s1.to.s0)

  if(sum.NA_change.s1.to.s0 > 0) {
    data[[phase_id[["phase.col"]]]] [change.s1.to.s0]<- 0
    data[change.s1.to.s0, which(colnames(data) %in% all.vars(formula.s1)[1])]<- NA
    m1<- message(paste("Changed the phase_id for ",sum.NA_change.s1.to.s0," rows to the first phase (s0) due to missingness in the set of auxiliary parameters
                       for the second phase (s1) (",s2_in_s1.nest.violation," terrestrial information no longer usable by this change)",sep = ""))
  }


  # check if every terrestrial plot has a response-value:
  deleted.terr <- data[[phase_id[["phase.col"]]]] == phase_id[["terrgrid.id"]] & !complete.cases(data[,all.vars(formula.s1)[1]])
  sum.deleted.terr<- sum(deleted.terr)

  # delete items with missing reponse information for s3-grid and produce message:
  if(sum.deleted.terr != 0) {
    data<- data[- which(deleted.terr),]
    m2<- message(paste("Additional ",sum.deleted.terr," rows deleted due to missing value for the response variable", sep = ""))
  }


  # -------------------------------------------------------------------------- #
  # -------------------------------------------------------------------------- #
  # NA-treatment:

  # # rows to be deleted due to missing auxiliary information or any input parameters:
  # deleted.s0<- !complete.cases(data [, which(colnames(data) %in% all.vars(formula.s0)[-1])]) # logical vector returning rows with missing entries
  # sum.NA_omitted.s0<- sum(deleted.s0)
  #
  # # delete missing rows in s0 of entire dataset and produce message:
  # if(sum.NA_omitted.s0 != 0) {
  #   data<- data[- which(deleted.s0),]
  #   m0<- message(paste(sum.NA_omitted.s0," rows deleted due to missingness in the auxiliary parameters or any of the input parameters in s1",sep = ""))
  # }
  #
  # # delete missing rows in s1-phase and produce message:
  # deleted.s1<- !complete.cases(data [, which(colnames(data) %in% all.vars(formula.s1)[-1])]) & data[[phase_id[["phase.col"]]]]==phase_id[["s1.id"]] # logical vector returning rows with missing entries
  # sum.NA_omitted.s1<- sum(deleted.s1)
  # if(sum.NA_omitted.s1 != 0) {
  #   data<- data[- which(deleted.s1),]
  #   m1<- message(paste(sum.NA_omitted.s1," rows deleted due to missingness in the auxiliary parameters or any of the input parameters in s2",sep = ""))
  # }
  #
  # # check if  every terrestrial plot has a response-value:
  # deleted.terr <- data[[phase_id[["phase.col"]]]] == phase_id[["terrgrid.id"]] & !complete.cases(data[,all.vars(formula.s1)[1]])
  # sum.deleted.terr<- sum(deleted.terr)
  #
  # # delete items with missing reponse information for s3-grid and produce message:
  # if(sum.deleted.terr != 0) {
  #   data<- data[- which(deleted.terr),]
  #   m2<- message(paste(sum.deleted.terr," rows deleted due to missing value for the response variable", sep = ""))
  # }


  # -------------------------------------------------------------------------- #
  # -------------------------------------------------------------------------- #

  #---------------------------------#
  # call for non-cluster functions: #
  # --------------------------------#

  if(is.na(cluster)) {

    # ---------------------------------------------------------------------#
    # check if --> "global 3-phase non-exhaustive" is required and apply estimator:

    if(is.na(small_area[["sa.col"]]) & all(is.na(exhaustive))) {

      # --- error checking -- :
      if(!is.na(boundary_weights)){ check.boundary_weightsInput(data, boundary_weights)}

      # -- # source global_nonexhaustive2p.R - function --:
      # source("global_nonexhaustive3p.R")

      # -- call function --:
      result<- global_nonexhaustive3p(formula.s0, formula.s1, data, phase_id, boundary_weights)

    }

    # ---------------------------------------------------------------------#
    # check if --> "global 3-phase exhaustive" is required and apply estimator:

    if(is.na(small_area[["sa.col"]]) & all(!is.na(exhaustive))) {

      # --- error checking -- :
      check.exhaustive3p(formula.s0, exhaustive, data)

      # -- # source global_exhaustive2p.R - function --:
      # source("global_exhaustive3p.R")

      # -- call function -- :
      result<- global_exhaustive3p(formula.s0, formula.s1, data, phase_id, boundary_weights, exhaustive)

    }

    # ---------------------------------------------------------------------#
    # check if --> "SMALL AREA 3-phase non-exhaustive" is required and apply estimator:

    if(!is.na(small_area[["sa.col"]]) & all(is.na(exhaustive))) {

      # --- error checking -- :
      check.smallareaInput(data, small_area)
      if(!is.na(boundary_weights)){ check.boundary_weightsInput(data, boundary_weights)}

      # -- # source small_area_looper_3p.R - function --:
      # source("small_area_looper_3p.R")

      # -- call function -- :
      if(!psmall){
        result <-  small_area_looper_3p(formula.s0, formula.s1, data, phase_id, cluster, small_area, boundary_weights, exhaustive, progressbar, psmall)
      } else {
        result<- psmall_fct3p(formula.s0, formula.s1, data, phase_id, cluster, small_area, boundary_weights, exhaustive, progressbar, psmall)
      }

    }

    # ---------------------------------------------------------------------#
    # check if --> "SMALL AREA 3-phase exhaustive" is required and apply estimator:

    if(!is.na(small_area[["sa.col"]]) & all(!is.na(exhaustive))) {

      # --- error checking -- :
      check.smallareaInput(data, small_area)
      check.exhaustive3p(formula.s0, exhaustive, data)

      # -- # source small_area_looper_3p.R - function --:
      # source("small_area_looper_3p.R")

      # -- call function -- :
      if(!psmall){
        result <- small_area_looper_3p(formula.s0, formula.s1, data, phase_id, cluster, small_area, boundary_weights, exhaustive, progressbar, psmall)
      } else {
        result<- psmall_fct3p(formula.s0, formula.s1, data, phase_id, cluster, small_area, boundary_weights, exhaustive, progressbar, psmall)
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
    if(!all(!is.na(data[["cluster"]]))){print(paste("WARNING: NAs removed from ",cluster))} #warning for NAs in cluster id
    data <- data[!is.na(data[["cluster"]]),] #strip NA clusterIDS


    # ---------------------------------------------------------------------#
    # check if --> "global 3-phase non-exhaustive_cluster" is required and apply estimator:

    if(is.na(small_area[["sa.col"]]) & all(is.na(exhaustive))) {

      # --- error checking -- :
      if(!is.na(boundary_weights)){ check.boundary_weightsInput(data, boundary_weights)}

      # -- # source ....R - function --:
      # source("global_nonexhaustive3p_cluster.R")

      # -- call function -- :
      result<- global_nonexhaustive3p_cluster(formula.s0, formula.s1, data, phase_id, cluster, boundary_weights)

    }

    # ---------------------------------------------------------------------#
    # check if --> "global 3-phase exhaustive_cluster" is required and apply estimator:

    if(is.na(small_area[["sa.col"]]) & all(!is.na(exhaustive))) {

      # --- error checking -- :
      check.exhaustive.cluster3p(formula.s0, exhaustive, data)

      # -- # source ....R - function --:
      # source("global_exhaustive3p_cluster.R")

      # -- call function -- :
      result<- global_exhaustive3p_cluster(formula.s0, formula.s1, data, phase_id, cluster, boundary_weights, exhaustive)

    }

    # ---------------------------------------------------------------------#
    # check if --> "SMALL AREA 3-phase non-exhaustive cluster" is required and apply estimator:

    if(!is.na(small_area[["sa.col"]]) & all(is.na(exhaustive))) {

      # --- error checking -- :
      check.smallareaInput(data, small_area)
      if(!is.na(boundary_weights)){ check.boundary_weightsInput(data, boundary_weights)}

      # -- # source ....p.R - function --:
      # source("small_area_looper_3p.R")

      # -- call function -- :
      if(!psmall){
        result <- small_area_looper_3p(formula.s0, formula.s1, data, phase_id, cluster, small_area, boundary_weights, exhaustive, progressbar, psmall)
      } else {
        result<- psmall_fct3p(formula.s0, formula.s1, data, phase_id, cluster, small_area, boundary_weights, exhaustive, progressbar, psmall)
      }


    }


    # ---------------------------------------------------------------------#
    # check if --> "SMALL AREA 3-phase exhaustive cluster" is required and apply estimator:

    if(!is.na(small_area[["sa.col"]]) & all(!is.na(exhaustive))) {

      # --- error checking -- :
      check.smallareaInput(data, small_area)
      check.exhaustive.cluster3p(formula.s0, exhaustive, data) ###

      # -- # source ....p.R - function --:
      # source("small_area_looper_3p.R")

      # -- call function -- :
      if(!psmall){
        result <- small_area_looper_3p(formula.s0, formula.s1, data, phase_id, cluster, small_area, boundary_weights, exhaustive, progressbar, psmall)
      } else {
        result<- psmall_fct3p(formula.s0, formula.s1, data, phase_id, cluster, small_area, boundary_weights, exhaustive, progressbar, psmall)
      }

    }


  } # end of non-cluster function calls


  # -------------------------------------------------------------------------- #
  # -------------------------------------------------------------------------- #

  # add function call to returned-list:

  result[["input"]]<- c(result[["input"]], call=call)

  result

} # end of Super-Function


# -------------------------------------------------------------------------- #
# -------------------------------------------------------------------------- #


