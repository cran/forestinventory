#' Plotting Estimation Results
#'
#' Function plots the estimation results of an object created by the \code{\link{estTable}} function.
#' Provides the possibility to visualize and compare the point estimates and their estimation errors
#' differentiated by the applied estimation method and estimator.
#'
#' @param x object of class \code{"list" "esttable"} created by the \code{\link{estTable}} function.
#'
#' @param yvar if set to \code{"error"} (default), the estimation error is plotted on the y-axis. If set to \code{"estimate"},
#'        point estimates with their confidence intervals are plotted.
#' @param ncol number of columns to plot small area estimations.
#' @param yscale.free \code{logical}: should y-axis scales be free (default) or fixed.
#' @param ... ignored.
#'
#'
#' @example examples/example_plot_estTable.R
#'
#' @import ggplot2
#' @import methods
#' @export

plot.esttable<- function(x, yvar="error", ncol=5, yscale.free=TRUE,...){

  # check input:
  if(!is(x, "esttable")){stop("'mphase.gain()' expects an 'esttable' object created by 'estTable()'")}

  if(is(x, "global")){
    etype<- "global"
  }

  if(is(x, "smallarea")){
    etype<- "smallarea"
  }


  dat<- as.data.frame(x)

  # add missing factor levels and reorder levels for plotting:
  levs.estimators<- c("onephase", "non-exhaustive", "exhaustive", "psynth extended", "synth extended", "psmall", "small", "psynth", "synth")
  levels(dat$estimator)<- append(levels(dat$estimator), levs.estimators[!(levs.estimators %in% levels(dat$estimator))])
  dat$estimator<- factor(dat$estimator, levels = levs.estimators)

  # add missing factor levels and reorder levels for plotting:
  levs.methest<- c("onephase(variance)", "non-exhaustive(ext_variance)", "non-exhaustive(g_variance)",
                   "psmall(ext_variance)", "psmall(g_variance)", "small(ext_variance)", "small(g_variance)",
                   "psynth extended(ext_variance)", "psynth extended(g_variance)", "synth extended(ext_variance)", "synth extended(g_variance)",
                   "psynth(g_variance)", "synth(g_variance)")
  dat$methest<-  factor(paste(dat$estimator,"(",dat$vartype,")",sep = ""))
  levels(dat$methest)<- append(levels(dat$methest), levs.methest[!(levs.methest %in% levels(dat$methest))])
  dat$methest<- factor(dat$methest, levels = levs.methest)


  # create color lookup-table and extract colorvector for current plotting:
  colors.tab<- data.frame( methest= c("onephase(variance)", "non-exhaustive(ext_variance)", "non-exhaustive(g_variance)",
                                      "psmall(ext_variance)", "psmall(g_variance)", "small(ext_variance)", "small(g_variance)",
                                      "psynth extended(ext_variance)", "psynth extended(g_variance)", "synth extended(ext_variance)", "synth extended(g_variance)",
                                      "psynth(g_variance)", "synth(g_variance)"),
                           colors=  c("chartreuse4", "deepskyblue1", "deepskyblue3",
                                      "deepskyblue1", "deepskyblue3", "lightskyblue1", "lightskyblue3",
                                      "royalblue1", "royalblue3", "lightsteelblue1", "lightsteelblue3",
                                      "red3","firebrick1"))
  colors.temp<- as.character(colors.tab$colors[colors.tab$methest %in% unique(dat$methest)])




  # ----------------------------------------------------------------------------- #
  if(yvar=="error"){


    # ************************* #
    if(etype =="global"){

      p<-  ggplot(data = dat, aes_q(x=quote(method), y=quote(error), fill=quote(methest))) +
        geom_bar(colour="black", stat="identity", position=position_dodge()) +
        labs(x="Estimation Method", y="Estimation Error [%]") +
        scale_fill_manual("Estimator", values=colors.temp)
    }


    # ************************* #
    if(etype =="smallarea"){

      if(yscale.free){scaleset<- "free_y"}
      if(!yscale.free){scaleset<- "fixed"}

      p<-  ggplot(data=dat, aes_q(x=quote(method), y=quote(error), fill=quote(methest))) +
        geom_bar(colour="black", stat="identity", position=position_dodge()) +
        labs(x="Estimation Method", y="Estimation Error [%]") +
        facet_wrap( ~ area, ncol=ncol, scales = scaleset) +
        scale_fill_manual("Estimator", values=colors.temp)
    }

  }



  # ----------------------------------------------------------------------------- #
  if(yvar=="estimate"){

    # ************************* #
    if(etype == "global"){

      p<- ggplot(data=dat, aes_q(x=quote(method), y=quote(estimate), fill=quote(methest))) +
        geom_bar(colour="black", stat="identity", position=position_dodge()) +
        geom_errorbar(aes_q(ymax=quote(ci_upper), ymin=quote(ci_lower)), position=position_dodge(width=0.9), width=0.25, lwd=0.75) +
        xlab("Estimation Method") +
        scale_fill_manual("Estimator", values=colors.temp)
    }


    # ************************* #
    if(etype == "smallarea"){

      if(yscale.free){scaleset<- "free_y"}
      if(!yscale.free){scaleset<- "fixed"}

      p<- ggplot(data=dat, aes_q(x=quote(method), y=quote(estimate), fill=quote(methest))) +
        geom_bar(colour="black", stat="identity", position=position_dodge()) +
        geom_errorbar(aes_q(ymax=quote(ci_upper), ymin=quote(ci_lower)), position=position_dodge(width=0.9), width=0.25, lwd=0.75) +
        xlab("Estimation Method") +
        facet_wrap( ~ area, ncol=ncol, scales = scaleset) +
        scale_fill_manual("Estimator", values=colors.temp)
    }

  }


  return(p)

} # end of plotting function
