#' Summarizing Global and Small-Area Estimation Results
#'
#' @param object object of class \code{onephase}, \code{twophase} or \code{threephase},
#'        containing estimation results of the respective estimation method.
#' @param coefs of type "\code{\link[base]{logical}}". If set to \code{TRUE}, also
#'               gives the regression coefficients of \code{\link{twophase}} and
#'               \code{\link{threephase}} estimations. Defaults to \code{FALSE}.
#' @param ... additional arguments, so far ignored.
#' @name summary
NULL
#>


#' @rdname summary
#' @export
summary.onephase<- function(object, coefs=FALSE, ...){
  # print-method for one-phase small area outputs

  stopifnot(inherits(object, "onephase"))

  cat("\n")
  cat("One-phase estimation")
  cat("\n \n")
  cat("Call: ")
  cat("\n")
  print(object$input$call)
  cat("\n")
  cat("Method used:")
  cat("\n")

  if (object$input$cluster){
    cat("One-phase estimator for cluster sampling")
  } else {
    cat("One-phase estimator")
  }

  cat("\n", "\n")
  cat("Estimation results:")
  cat("\n")
  print(object$estimation, row.names = FALSE)
  cat("\n\n")

}


#' @rdname summary
#' @export
summary.twophase<- function(object, coefs=FALSE, ...){
  # summary for small area estimations:

  stopifnot(inherits(object, "twophase"))

  # --------------------------------#
  # summary for twophase-smallarea:

  if(class(object)[1]=="smallarea"){ # if class(sae_obj) is c("smallarea", "twophase")

    cat("\n")
    cat("Two-phase small area estimation")
    cat("\n \n")
    cat("Call: ")
    cat("\n")
    print(object$input$call)
    cat("\n")
    cat("Method used:")
    cat("\n")

    s<- FALSE # indicator for displaying only once the global calculation for coefficients

    if(object$input$exhaustive & !object$input$cluster){ # exhaustive & non-cluster
      if(object$input$method == "synth") { cat("Synthetic small area estimator"); s<- TRUE}
      if(object$input$method == "synth extended") { cat("Extended synthetic small area estimator")}
      if(object$input$method == "psmall"){ cat("Small area estimator"); s<- TRUE}
    }


    if(object$input$exhaustive & object$input$cluster){ # exhaustive & cluster
      if(object$input$method == "synth") { cat("Synthetic small area estimator for cluster sampling"); s<- TRUE}
      if(object$input$method == "synth extended") { cat("Extended synthetic small area estimator for cluster sampling")}
      if(object$input$method == "psmall"){ cat("Small area estimator for cluster sampling"); s<- TRUE}
    }


    if(!object$input$exhaustive & !object$input$cluster){ # non-exhaustive & non-cluster
      if(object$input$method == "psynth") { cat("Pseudosynthetic small area estimator"); s<- TRUE}
      if(object$input$method == "psynth extended") { cat("Extended pseudosynthetic small area estimator")}
      if(object$input$method == "psmall"){ cat("Pseudo small area estimator"); s<- TRUE}
    }


    if (!object$input$exhaustive & object$input$cluster){ # non-exhaustive & cluster
      if(object$input$method == "psynth") { cat("Pseudosynthetic small area estimator for cluster sampling"); s<- TRUE}
      if(object$input$method == "psynth extended") { cat("Extended pseudosynthetic small area estimator for cluster sampling")}
      if(object$input$method == "psmall"){ cat("Pseudo small area estimator for cluster sampling"); s<- TRUE}
    }

    cat("\n", "\n")
    cat("Regression Model:")
    cat("\n")
    print(object$input$formula, showEnv=FALSE)
    cat("\n")
    if(coefs){
      cat("Regression Coefficients:")
      cat("\n")
      ifelse(s, print(object$coefficients[1,-1],row.names = FALSE), print(object$coefficients, row.names = FALSE))
      cat("\n")
    }
    cat("Estimation results:")
    cat("\n")

    print(object$estimation, row.names = FALSE)
    cat("\n")

    if(!is.na(object$input$boundary_weights)){
      cat("'boundary_weight'- option was used to calculate weighted means of auxiliary variables")
      cat("\n\n")
    }
  } # end of smallarea-summary


  # ------------------------------#
  # summary for twophase-global:

  if(class(object)[1]=="global"){ # if class(sae_obj) is c("global", "twophase")

    cat("\n")
    cat("Two-Phase global estimation")
    cat("\n \n")
    cat("Call: ")
    cat("\n")
    print(object$input$call)
    cat("\n")
    cat("Method used:")
    cat("\n")

    if (object$input$exhaustive & !object$input$cluster){ # exhaustive & non-cluster
      cat("Exhaustive global estimator")
    }


    if (object$input$exhaustive & object$input$cluster){ # exhaustive & cluster
      cat("Exhaustive global estimator for cluster sampling")
    }


    if (!object$input$exhaustive & !object$input$cluster){ # non-exhaustive & non-cluster
      cat("Non-exhaustive global estimator")
    }


    if (!object$input$exhaustive & object$input$cluster){ # non-exhaustive & cluster
      cat("Non-exhaustive global estimator for cluster sampling")
    }


    cat("\n", "\n")
    cat("Regression Model:")
    cat("\n")
    print(object$input$formula, showEnv=FALSE)
    cat("\n")
    if (coefs){
      cat("Regression Coefficients:")
      cat("\n")
      print(object$coefficients, row.names = FALSE)
      cat("\n")
    }
    cat("Estimation results:")
    cat("\n")

    print(object$estimation, row.names = FALSE)
    cat("\n")

    if(!is.na(object$input$boundary_weights)){
      cat("'boundary_weight'- option was used to calculate weighted means of auxiliary variables")
      cat("\n\n")
    }
  }# end of global-summary

} # end of summary.twophase


#' @rdname summary
#' @export
summary.threephase<- function(object, coefs=FALSE, ...){
  # summary for threephase estimations:

  stopifnot(inherits(object, "threephase"))

  # --------------------------------#
  # summary for threephase-smallarea:

  if(class(object)[1]=="smallarea"){

    cat("\n")
    cat("Three-phase small area estimation")
    cat("\n \n")
    cat("Call: ")
    cat("\n")
    print(object$input$call)
    cat("\n")
    cat("Method used:")
    cat("\n")

    s<- FALSE # indicator for displaying only once the global calculation for coefficients

    if (object$input$exhaustive & !object$input$cluster){ # exhaustive & non-cluster
      if(object$input$method == "synth") { cat("Synthetic small area estimator"); s<- TRUE}
      if(object$input$method == "synth extended") { cat("Extended synthetic small area estimator")}
      if(object$input$method == "psmall"){ cat("Small area estimator"); s<- TRUE}
    }


    if (object$input$exhaustive & object$input$cluster){ # exhaustive & cluster
      if(object$input$method == "synth") { cat("Synthetic small area estimator for cluster sampling"); s<- TRUE}
      if(object$input$method == "synth extended") { cat("Extended synthetic small area estimator for cluster sampling")}
      if(object$input$method == "psmall"){ cat("Small area estimator for cluster sampling"); s<- TRUE}
    }


    if (!object$input$exhaustive & !object$input$cluster){ # non-exhaustive & non-cluster
      if(object$input$method == "psynth") { cat("Pseudosynthetic small area estimator"); s<- TRUE}
      if(object$input$method == "psynth extended") { cat("Extended pseudosynthetic small area estimator")}
      if(object$input$method == "psmall"){ cat("Pseudo small area estimator"); s<- TRUE}
    }


    if (!object$input$exhaustive & object$input$cluster){ # non-exhaustive & cluster
      if(object$input$method == "psynth") { cat("Pseudosynthetic small area estimator for cluster sampling"); s<- TRUE}
      if(object$input$method == "psynth extended") { cat("Extended pseudosynthetic small area estimator for cluster sampling")}
      if(object$input$method == "psmall"){ cat("Pseudo small area estimator for cluster sampling"); s<- TRUE}
    }

    cat("\n", "\n")
    cat("Full Regression Model:")
    cat("\n")
    print(object$input$formula.s1, showEnv=FALSE)
    cat("\n")
    cat("Reduced Regression Model:")
    cat("\n")
    print(object$input$formula.s0, showEnv=FALSE)
    cat("\n")
    if(coefs){
      cat("Summary of Coefficients:")
      cat("\n")
        if(s){
          coefs.cut<- convert_coefs_table.smallarea3p(object)[c(1,2), ]; row.names(coefs.cut)<- c("", "*")
          print(as.matrix(coefs.cut), na.print = "", quote = F)
        } else {
          print(as.matrix(convert_coefs_table.smallarea3p(object)), na.print = "", quote = F)
        }
      cat("\n")
      cat(" Coefficients of reduced model indicated by '*'")
      cat("\n", "\n")
    }

    cat("Estimation results:")
    cat("\n")

    print(object$estimation, row.names = FALSE)
    cat("\n")

    if(!is.na(object$input$boundary_weights)){
      cat("'boundary_weight'- option was used to calculate weighted means of auxiliary variables")
      cat("\n\n")
    }
  }# end of smallarea-summary


  # --------------------------------#
  # summary for threephase-global:

  if(class(object)[1]=="global"){

    cat("\n")
    cat("Three-phase global estimation")
    cat("\n \n")
    cat("Call: ")
    cat("\n")
    print(object$input$call)
    cat("\n")
    cat("Method used:")
    cat("\n")

    if (object$input$exhaustive & !object$input$cluster){ # exhaustive & non-cluster
      cat("Exhaustive global estimator")
    }


    if (object$input$exhaustive & object$input$cluster){ # exhaustive & cluster
      cat("Exhaustive global estimator for cluster sampling")
    }


    if (!object$input$exhaustive & !object$input$cluster){ # non-exhaustive & non-cluster
      cat("Non-exhaustive global estimator")
    }


    if (!object$input$exhaustive & object$input$cluster){ # non-exhaustive & cluster
      cat("Non-exhaustive global estimator for cluster sampling")
    }

    cat("\n", "\n")
    cat("Full Regression Model:")
    cat("\n")
    print(object$input$formula.s1, showEnv=FALSE)
    cat("\n")
    cat("Reduced Regression Model:")
    cat("\n")
    print(object$input$formula.s0, showEnv=FALSE)
    cat("\n")
    if(coefs){
      cat("Summary of Coefficients:")
      cat("\n")
      print(as.matrix(convert_coefs_table.global3p(object)), na.print = "", quote = F)
      cat("\n")
      cat(" Coefficients of reduced model indicated by '*'")
      cat("\n", "\n")
    }

    cat("Estimation results:")
    cat("\n")

    print(object$estimation, row.names = FALSE)
    cat("\n")

    if(!is.na(object$input$boundary_weights)){
      cat("'boundary_weight'- option was used to calculate weighted means of auxiliary variables")
      cat("\n\n")
    }
  } # end of global-summary

}# end of summary.threephase
