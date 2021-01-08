# --------------------------------------------------- #
# Takes as input an esttable - object and returns a validation of
# which multiphase-methods performed BEST in comparison to the onephase
# (baseline) std-error.
#
# We define the "gain" as the reduction ("+") or possibly also the increase ("-")
# of the best multiphase-standard-error compared to the onephase-std-error
#
# We also give the percentage of the best multiphase-standard-error
# compared to the onephase std-error AND the relative efficiency
#
# --------------------------------------------------- #

#' mphase.gain
#'
#' \code{mphase.gain} takes as input an object created by the \code{\link{estTable}} function
#' and returns a validation of which multiphase method and estimator performed best in comparison
#' to the onephase estimation (baseline) in terms of estimation precision.
#'
#'
#' @param esttable.obj an object of class \code{esttable} created by the \code{\link{estTable}} function
#'
#' @param pref.vartype preferred type of multiphase variance that should be compared to the \code{onephase} variance,
#'                     if more then one variance type has been calculated in the multiphase estimation object(s) stored in
#'                     \code{esttable}. Valid input values are \code{"g_variance"} (default) and \code{"ext_variance"}.
#'
#' @param exclude.synth \code{logical}. If set to \code{TRUE} (default), synthetic estimations are not considered in the validation.
#'
#'
#' @return \code{mphase.gain} returns a \code{data.frame} containing the following components:
#'
#'  \itemize{
#'     \item \code{area:} in case of small area estimation: the name of the small area
#'     \item \code{var_onephase:} standard error of the \code{\link{onephase}} estimation
#'     \item \code{var_multiphase:} smallest variance among the (set of) multiphase estimations stored in \code{esttable.obj}
#'     \item \code{estimator:} multiphase estimator with the smallest variance
#'     \item \code{method:} estimation Method of the multiphase estimator with the smallest variance
#'     \item \code{gain:} the \emph{gain} is the reduction (if value is positive) or possibly also the increase (if value is negative)
#'                        in variance when applying the multiphase as alternative to the onephase estimation
#'     \item \code{rel.eff:} the \emph{relative efficiency} defined as the ratio between the onephase variance and the multiphase variance
#'     %\item \code{perc.of.onephase:} ratio between the smallest multiphase standard error and the onephase standard error
#'  }
#'
#' @note
#'
#' The \emph{gain} can be interpreted as: "The multiphase estimation procedure leads to a \code{gain} \% reduction in variance compared to the
#' onephase procedure".
#'
#' The \emph{relative efficiency} can be interpreted as: "Using the onephase estimation procedure, the terrestrial sample size would have to be \code{rel.eff} times larger in order to achieve the same precision (in terms of variance) as the mutiphase estimation procedure".
#'
#'
#' % @example examples/example_mphasegain.R
#'
#' @import methods
#' @export


# -----------------------------------------------------------------------------#
# FUNCTION STARTS HERE:


mphase.gain<- function(esttable.obj, pref.vartype = "g_variance", exclude.synth = TRUE){

  # check input:
  if(!is(esttable.obj, "esttable")){stop("'mphase.gain()' expects an 'esttable' object created by 'estTable()'")}

  if(is(esttable.obj, "global")){
    etype<- "global"
  }

  if(is(esttable.obj, "smallarea")){
    etype<- "smallarea"
  }


  # convert esttable-object into data.frame:
  esttable.obj<- as.data.frame(esttable.obj)

  # -------- #
  # closure:
  prec.gain<- function(est.tab){

    ind.1ph<- est.tab$vartype %in% "variance"

    if(exclude.synth){
      ind.not.1ph<- est.tab$vartype %in% c(pref.vartype) & !(est.tab$estimator %in% c("psynth", "synth"))
    } else {
      ind.not.1ph<- est.tab$vartype %in% c(pref.vartype)
    }

    var.1ph<- est.tab[ind.1ph, "variance"]
    est.tab.not1ph<- est.tab[ind.not.1ph,]

    ind.best.multiph<- which.min(est.tab.not1ph[["variance"]])[1]

    best.mphase<- est.tab.not1ph[ind.best.multiph, c("estimator", "method", "variance", "vartype")]

    if(all(!ind.1ph)){ # if no onephase is there

      d<- data.frame(var_onephase = NA, var_multiphase = best.mphase$variance,
                     estimator = best.mphase$estimator, method = best.mphase$method,
                     gain = NA, rel.eff = NA)
    } else {

      # gain:
      red.to.1ph<- round(100* (1 - (best.mphase[["variance"]] / var.1ph)), digits = 1)

      perc.of.1phvar<- round(100*(best.mphase[["variance"]] / var.1ph), digits = 1)

      rel.eff<- var.1ph / best.mphase[["variance"]]

      d<-
        data.frame(var_onephase = var.1ph, var_multiphase = best.mphase$variance,
                   method = best.mphase$method, estimator = best.mphase$estimator,
                   gain = red.to.1ph, rel.eff = rel.eff)
    }

    return(d)

  }
  # -------- #

  if(etype == "global"){
    # apply closure:
      return(prec.gain(esttable.obj))
  }

  if(etype == "smallarea"){
    # apply closure:
      return(ddply(esttable.obj,"area", prec.gain))
  }


} # end of function
