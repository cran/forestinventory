#' Calculates Confidence Intervals for Global and Small-Area Estimations
#'
#' @param object object of class \code{onephase}, \code{twophase} or \code{threephase},
#'        containing estimation results of the respective estimation method.
#' @param parm  ignored.
#' @param level the confidence level required.
#' @param ... additional arguments, so far ignored.
#' @param adjust.method correction method to obtain simultaneous confidence intervals
#'                      for a set of estimates (thus restricted to objects of class \code{"onephase"},
#'                      \code{c("smallarea", "twophase")} and  \code{c("smallarea", "threephase")}).
#'                      Available correction methods are \code{c("none","bonferroni")}. Defaults to \code{"none"}.
#'
#' @details
#'
#' Depending on the estimation method specified, \code{confint()} computes confidence intervals as follows:
#'
#' \code{\link{onephase}}:
#'
#' Two-sided confidence intervals are computed based on the t-distribution with \code{n2 - p} \emph{degrees of freedom},
#' where \code{n2} is the number of terrestrial data in the respective inventory domain.
#'
#' \code{\link{twophase}}:
#'
#' The calculation of the two-sided confidence intervals for \emph{global} twophase estimates
#' (objects of class \code{global}) are calculated based on the quantiles of the \emph{t}-distribution
#' with \code{n2 - p} \emph{degrees of freedom}, where \code{p} is the number of parameters used in
#' the regression model, and \code{n2} is the number of terrestrial observations (i.e. \emph{local densities})
#' in the inventory domain.
#'
#' The calculation of the two-sided confidence intervals for \emph{smallarea} twophase estimates
#' (objects of class \code{smallarea}) are calculated based on the quantiles of the \emph{t}-distribution
#' with \code{n2G - 1} \emph{degrees of freedom}, where \code{n2G} is the number of
#' terrestrial observations (i.e. \emph{local densities}) in the smallarea.
#'
#' \code{\link{threephase}}:
#'
#' The calculation of the two-sided confidence intervals for \emph{global} threephase estimates
#' (objects of class \code{global}) are calculated based on the quantiles of the \emph{t}-distribution
#' with \code{n2 - p} \emph{degrees of freedom}, where \code{p} is the number of parameters used in
#' the \strong{full} regression model, and \code{n2} is the number of terrestrial observations
#' (i.e. \emph{local densities}) in the inventory domain (note: in notation used here n0, n1 and n2
#' correspond to the first, second and third phase sample sizes respectively).
#'
#' The calculation of the two-sided confidence intervals for \emph{smallarea} theephase estimates
#' (objects of class \code{smallarea}) are calculated based on the quantiles of the \emph{t}-distribution
#' with \code{n2G - 1} \emph{degrees of freedom}, where \code{n2G} is the number of
#' terrestrial observations (i.e. \emph{local densities}) in the smallarea.
#'
#'@return \code{confint} returns a list of  the following 3 components:
#'
#'  \item{ci}{a \code{data.frame} containing the columns:
#'                    \itemize{
#'                    \item \code{area} the domain, i.e. small area
#'                    \item \code{ci_lower_ext} the lower confidence limit based on the external variance
#'                    \item \code{ci_upper_ext} the upper confidence limit based on the external variance
#'                    \item \code{ci_lower_g}   the lower confidence limit based on the g-weight variance
#'                    \item \code{ci_upper_g}   the upper confidence limit based on the g-weight variance
#'                    }}
#'  \item{level}{the applied confidence level}
#'  \item{adjust.method}{the adjustment method applied to retrieve simultaneous confidence intervals}
#'
#'
#'
#' @note
#' In the special case of \emph{synthetic} smallarea estimations, the two-sided confidence intervals
#' are calculated based on the quantiles of the \emph{t}-distribution
#' with \code{n2 - p} \emph{degrees of freedom}, i.e. based on the global sample size.
#'
#' The confidence intervals for \emph{synthetic} smallarea estimations do not account for the potential
#' bias of a linear model that was fit in a large forest area and applied to a small area.  Thus,
#' the coverage rates for confidence intervals produced by synthetic estimators may be less than
#' the nominal level of confidence.
#'
#' In case of cluster-sampling, \code{n2G} is the number of terrestrial clusters
#' (a cluster constitutes the sample unit). This is automatically considered by \code{confint}.
#'
#' The adjustment methods passed to \code{adjust.method} are designed to achieve
#' \emph{simultaneous} confidence intervals by correcting the confidence level given by \code{level}.
#' The use of this option is recommended if a set of estimates contained in a \code{onephase}- or
#' \code{smallarea}-object should be compared by their confidence intervals. It ensures that the
#' percentage of confidence intervals containing the true value will correspond to
#' the nominal confidence level.
#'
#'
#' @references
#'
#' Mandallaz, D. (2013). \emph{Design-based properties of some small-area estimators in forest inventory
#' with two-phase sampling.} Canadian Journal of Forest Research, 43(5), 441-449.
#'
#' Mandallaz, D., Breschan, J., & Hill, A. (2013). \emph{New regression estimators in forest inventories
#' with two-phase sampling and partially exhaustive information: a design-based monte carlo approach
#' with applications to small-area estimation.} Canadian Journal of Forest Research,
#' 43(11), 1023-1031.
#'
#' Mandallaz, D. (2013). \emph{A three-phase sampling extension of the generalized regression estimator
#' with partially exhaustive information.} Canadian Journal of Forest Research,
#' 44(4), 383-388.
#'
#' Benjamini, Y., and Hochberg, Y. (1995). \emph{Controlling the false discovery rate:
#' a practical and powerful approach to multiple testing.}
#' Journal of the Royal Statistical Society Series B 57, 289-300.
#'
#' @example examples/example_confint.R
#'
#' @name confint
NULL
#>

# -----------------------------------------------------------------------------#
# -----------------------------------------------------------------------------#

#' @rdname confint
#' @export
confint.onephase<- function(object, parm, level = 0.95, adjust.method="none",...){
  # confint-method for one-phase outputs:

  orig.level<- level

  # ----------#
  # adjust confidence level if required:
  if(adjust.method %in% c("none", "bonferroni")){

    # bonferroni-correction:
    if (adjust.method == "bonferroni"){
      if(nrow(object$estimation)==1){
        stop("'bonferroni-correction only reasonable for objects that contain a set of estimates'")
      } else {
         level<- 1 - ((1 - level) / nrow(object$estimation))
        }
    }

  } else {
    stop("invalid argument for 'adjust.method'")
  }
  # ----------#

  # calculating the quantile of the t-distribution with n2 - 1 df:
  t.val<- suppressWarnings(qt(p = 1-( (1-level)/2 ), df = object$estimation$n2 - 1))

  # calculating the lower and upper confidence intervals based on the design-based variance:
  ci_lower_op<- object$estimation$estimate - (t.val * sqrt(object$estimation$variance))
  ci_upper_op<- object$estimation$estimate + (t.val * sqrt(object$estimation$variance))

  if ("area"  %in% names(object$estimation)){
    ci<- list(ci=data.frame(area=object$estimation$area, estimate=object$estimation$estimate, ci_lower_op, ci_upper_op), level=orig.level, adjust.method=adjust.method)
  } else {
    ci<- list(ci=data.frame(estimate=object$estimation$estimate, ci_lower_op, ci_upper_op), level=orig.level, adjust.method=adjust.method)
  }

  class(ci)<- c("confint.global", "onephase")

  return(ci)
} # end of confint.onephase



# -----------------------------------------------------------------------------#
# -----------------------------------------------------------------------------#

#' @rdname confint
#' @export
confint.twophase<- function(object, parm, level = 0.95, adjust.method="none",...){
  # confint-method for for twophase estimations:

  orig.level<- level

  # -------------------------------------------#
  # confidence-intervals for twophase-smallarea:

  if(class(object)[1]=="smallarea"){ # if class(object) is c("smallarea", "twophase")

    # ----------#
    # adjust confidence level if required:
    if(adjust.method %in% c("none", "bonferroni")){

      # bonferroni-correction:
      if (adjust.method == "bonferroni"){
        if(nrow(object$estimation)==1){
          stop("'bonferroni-correction only reasonable for objects that contain a set of estimates'")
        } else {
          level<- 1 - ((1 - level) / nrow(object$estimation))
        }
      }

    } else {
      stop("invalid argument for 'adjust.method'")
    }

    # ----------#
    # calulating the quantile of the t-distribution:

    # if synthetic-estimation was applied, the df's are calculated as n2 - # parameters of regression model:
    if(object$input$method %in% c("synth", "psynth")){

      t.val<- suppressWarnings(qt(p = 1-( (1-level)/2 ), df = unique(object$estimation[["n2"]]) - length(all.vars(object$input$formula)[-1])))

      # if non-synthetic-estimation was applied, the df's are calculated as n2G - 1:
    } else {
      t.val<- suppressWarnings(qt(p = 1-( (1-level)/2 ), df = object$estimation$n2G - 1))
    }

    #-----------#
    # calculating the lower and upper confidence intervals based on the external variance:
    ci_lower_g<- object$estimation$estimate - (t.val * sqrt(object$estimation$g_variance))
    ci_upper_g<- object$estimation$estimate + (t.val * sqrt(object$estimation$g_variance))

    # calculating the lower and upper confidence intervals based on the design-based variance:
    ci_lower_ext<- object$estimation$estimate - (t.val * sqrt(object$estimation$ext_variance))
    ci_upper_ext<- object$estimation$estimate + (t.val * sqrt(object$estimation$ext_variance))

    ci<- list(ci=data.frame(area=object$estimation$area, estimate=object$estimation$estimate, ci_lower_ext, ci_upper_ext, ci_lower_g, ci_upper_g),
              level=orig.level, adjust.method=adjust.method)

    class(ci)<- c("confint.smallarea", "twophase")

    return(ci)

  } # end of smallarea-CI

  # -------------------------------------------#
  # confidence-intervals for twophase-global:

  if(class(object)[1]=="global"){ # if class(object) is c("global", "twophase")

    if(adjust.method != "none"){
      stop("'bonferroni-correction only reasonable for objects that contain a set of estimates'")
    }

    n2<- ifelse(object$input$cluster, object$samplesizes$n2_clust, object$samplesizes$n2)

    t.val<- qt(p = 1-( (1-level)/2 ), df = n2 - length(all.vars(object$input$formula)[-1]))

    ci_lower_g<- object$estimation[["estimate"]] - (t.val * sqrt(object$estimation[["g_variance"]]))
    ci_upper_g<- object$estimation[["estimate"]] + (t.val * sqrt(object$estimation[["g_variance"]]))

    ci_lower_ext<- object$estimation[["estimate"]] - (t.val * sqrt(object$estimation[["ext_variance"]]))
    ci_upper_ext<- object$estimation[["estimate"]] + (t.val * sqrt(object$estimation[["ext_variance"]]))

    ci<- list(ci=data.frame(estimate=object$estimation[["estimate"]], ci_lower_ext, ci_upper_ext, ci_lower_g, ci_upper_g),
              level=orig.level, adjust.method=adjust.method)

    class(ci)<- c("confint.global", "twophase")

    return(ci)
  }# end of global-CI

} ## end of confint.twophase



# -----------------------------------------------------------------------------#
# -----------------------------------------------------------------------------#

#' @rdname confint
#' @export
confint.threephase<- function(object, parm, level = 0.95, adjust.method="none", ...){
  # confint-method for for threephase estimations:

  orig.level<- level

  # --------------------------------#
  # confint for threephase-smallarea:

  if(class(object)[1]=="smallarea"){

    #-----------#
    # adjust confidence level if required:
    if(adjust.method %in% c("none", "bonferroni")){

      # bonferroni-correction:
      if (adjust.method == "bonferroni"){
        if(nrow(object$estimation)==1){
          stop("'bonferroni-correction only reasonable for objects that contain a set of estimates'")
        } else {
          level<- 1 - ((1 - level) / nrow(object$estimation))
        }
      }

    } else {
      stop("invalid argument for 'adjust.method'")
    }

    #-----------#
    # calulating the quantile of the t-distribution:

    # if synthetic-estimation was applied, the df's are calculated as n2 - # parameters of regression model:
    if(object$input$method %in% c("synth", "psynth")){

      t.val<- suppressWarnings(qt(p = 1-( (1-level)/2 ), df = unique(object$estimation[["n2"]]) - length(all.vars(object$input$formula.s1)[-1])))

      # if non-synthetic-estimation was applied, the df's are calculated as n2G - 1:
    } else {

      t.val<- suppressWarnings(qt(p = 1-( (1-level)/2 ), df = object$estimation$n2G - 1))
    }
    #-----------#

    # calculating the lower and upper confidence intervals based on the external variance:
    ci_lower_g<- object$estimation$estimate - (t.val * sqrt(object$estimation$g_variance))
    ci_upper_g<- object$estimation$estimate + (t.val * sqrt(object$estimation$g_variance))

    # calculating the lower and upper confidence intervals based on the design-based variance:
    ci_lower_ext<- object$estimation$estimate - (t.val * sqrt(object$estimation$ext_variance))
    ci_upper_ext<- object$estimation$estimate + (t.val * sqrt(object$estimation$ext_variance))

    ci<- list(ci=data.frame(area=object$estimation$area, estimate=object$estimation$estimate, ci_lower_ext, ci_upper_ext, ci_lower_g, ci_upper_g),
              level=orig.level, adjust.method=adjust.method)

    class(ci)<- c("confint.smallarea", "threephase")

    return(ci)
  } # end of smallarea-CI


  # --------------------------------#
  # summary for threephase-global:

  if(class(object)[1]=="global"){

    if(adjust.method != "none"){
      stop("'bonferroni-correction only reasonable for objects that contain a set of estimates'")
    }

    n2<- ifelse(object$input$cluster, object$samplesizes$n2_clust, object$samplesizes$n2)

    t.val<- qt(p = 1-( (1-level)/2 ), df = n2 - length(all.vars(object$input$formula.s1)[-1])) # that's however not recommended by Daniel,
                                                                                                 # see Article Crfr, Threephase, page 387
    ci_lower_g<- object$estimation[["estimate"]] - (t.val * sqrt(object$estimation[["g_variance"]]))
    ci_upper_g<- object$estimation[["estimate"]] + (t.val * sqrt(object$estimation[["g_variance"]]))

    ci_lower_ext<- object$estimation[["estimate"]] - (t.val * sqrt(object$estimation[["ext_variance"]]))
    ci_upper_ext<- object$estimation[["estimate"]] + (t.val * sqrt(object$estimation[["ext_variance"]]))

    ci<- list(ci=data.frame(estimate=object$estimation[["estimate"]], ci_lower_ext, ci_upper_ext, ci_lower_g, ci_upper_g),
              level=orig.level,adjust.method=adjust.method)

    class(ci)<- c("confint.global", "threepase")

    return(ci)
  }# end of global-CI

}## end of confint.threephase


