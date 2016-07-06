


#' @export
print.onephase<- function(x, ...){
  # print-method for one-phase outputs:

  cat("\n")
  cat("One-phase estimation")
  cat("\n \n")
  cat("Call: ")
  cat("\n")
  print(x$input$call)
  cat("\n")
  cat("Method used:")
  cat("\n")

  if (x$input$cluster){
    cat("One-phase estimator for cluster sampling")
  } else {
    cat("One-phase estimator")
  }

  cat("\n", "\n")

  # cat("\n")
  # cat("Number of areas calculated: ",nrow(sae_obj$estimation))
  # cat("\n \n")

}


#' @export
print.twophase<- function(x, ...){
  # print-method for all two-phase outputs:

  # --------------------------------#
  # summary for twophase-smallarea:

  if(class(x)[1]=="smallarea"){ # if class(x) is c("smallarea", "twophase")

    cat("\n")
    cat("Two-phase small area estimation")
    cat("\n \n")
    cat("Call: ")
    cat("\n")
    print(x$input$call)
    cat("\n")
    cat("Method used:")
    cat("\n")


    if (x$input$exhaustive & !x$input$cluster){ # exhaustive & non-cluster
      if(x$input$method == "synth") { cat("Synthetic small area estimator")}
      if(x$input$method == "synth extended") { cat("Extended synthetic small area estimator")}
      if(x$input$method == "psmall"){ cat("Small area estimator")}
    }


    if (x$input$exhaustive & x$input$cluster){ # exhaustive & cluster
      if(x$input$method == "synth") { cat("Synthetic small area estimator for cluster sampling")}
      if(x$input$method == "synth extended") { cat("Extended synthetic small area estimator for cluster sampling")}
      if(x$input$method == "psmall"){ cat("Small area estimator for cluster sampling")}
    }


    if (!x$input$exhaustive & !x$input$cluster){ # non-exhaustive & non-cluster
      if(x$input$method == "psynth") { cat("Pseudosynthetic small area estimator")}
      if(x$input$method == "psynth extended") { cat("Extended pseudosynthetic small area estimator")}
      if(x$input$method == "psmall"){ cat("Pseudo small area estimator")}
    }


    if (!x$input$exhaustive & x$input$cluster){ # non-exhaustive & cluster
      if(x$input$method == "psynth") { cat("Pseudosynthetic small area estimator for cluster sampling")}
      if(x$input$method == "psynth extended") { cat("Extended pseudosynthetic small area estimator for cluster sampling")}
      if(x$input$method == "psmall"){ cat("Pseudo small area estimator for cluster sampling")}
    }


    cat("\n", "\n")
    #   cat("Regression Model:")
    #   cat("\n")
    #   print(x$input$formula, showEnv=FALSE)
    cat("\n")
    cat("Number of small areas calculated: ",nrow(x$estimation))
    cat("\n \n")
  } # end of smallarea-print


  # ------------------------------#
  # print for twophase-global:

  if(class(x)[1]=="global"){ # if class(x) is c("global", "twophase")

    cat("\n")
    cat("Two-phase global estimation")
    cat("\n \n")
    cat("Call: ")
    cat("\n")
    print(x$input$call)
    cat("\n")
    cat("Method used:")
    cat("\n")


    if (x$input$exhaustive & !x$input$cluster){ # exhaustive & non-cluster
      cat("Exhaustive global estimator")
    }


    if (x$input$exhaustive & x$input$cluster){ # exhaustive & cluster
      cat("Exhaustive global estimator for cluster sampling")
    }


    if (!x$input$exhaustive & !x$input$cluster){ # non-exhaustive & non-cluster
      cat("Non-exhaustive global estimator")
    }


    if (!x$input$exhaustive & x$input$cluster){ # non-exhaustive & cluster
      cat("Non-exhaustive global estimator for cluster sampling")
    }

    #   cat("\n", "\n")
    #   cat("Regression Model:")
    #   cat("\n")
    #   print(x$input$formula, showEnv=FALSE)
    cat("\n \n")
  }

}# end of print.twophase



#' @export
print.threephase<- function(x, ...){
  # print-method for all three-phase outputs:

  # --------------------------------#
  # summary for threephase-smallarea:

  if(class(x)[1]=="smallarea"){

    cat("\n")
    cat("Three-phase small area estimation")
    cat("\n \n")
    cat("Call: ")
    cat("\n")
    print(x$input$call)
    cat("\n")
    cat("Method used:")
    cat("\n")


    if (x$input$exhaustive & !x$input$cluster){ # exhaustive & non-cluster
      if(x$input$method == "synth") { cat("Synthetic small area estimator")}
      if(x$input$method == "synth extended") { cat("Extended synthetic small area estimator")}
      if(x$input$method == "psmall"){ cat("Small area estimator")}
    }


    if (x$input$exhaustive & x$input$cluster){ # exhaustive & cluster
      if(x$input$method == "synth") { cat("Synthetic small area estimator for cluster sampling")}
      if(x$input$method == "synth extended") { cat("Extended synthetic small area estimator for cluster sampling")}
      if(x$input$method == "psmall"){ cat("Small area estimator for cluster sampling")}
    }


    if (!x$input$exhaustive & !x$input$cluster){ # non-exhaustive & non-cluster
      if(x$input$method == "psynth") { cat("Pseudosynthetic small area estimator")}
      if(x$input$method == "psynth extended") { cat("Extended pseudosynthetic small area estimator")}
      if(x$input$method == "psmall"){ cat("Pseudo small area estimator")}
    }


    if (!x$input$exhaustive & x$input$cluster){ # non-exhaustive & cluster
      if(x$input$method == "psynth") { cat("Pseudosynthetic small area estimator for cluster sampling")}
      if(x$input$method == "psynth extended") { cat("Extended pseudosynthetic small area estimator for cluster sampling")}
      if(x$input$method == "psmall"){ cat("Pseudo small area estimator for cluster sampling")}
    }


    cat("\n")
    #   cat("Regression Model:")
    #   cat("\n")
    #   print(x$input$formula, showEnv=FALSE)
    cat("\n")
    cat("Number of small areas calculated: ",nrow(x$estimation))
    cat("\n \n")

  }# end of smallarea-summary

  # --------------------------------#
  # summary for threephase-global:

  if(class(x)[1]=="global"){

    cat("\n")
    cat("Three-phase global estimation")
    cat("\n \n")
    cat("Call: ")
    cat("\n")
    print(x$input$call)
    cat("\n")
    cat("Method used:")
    cat("\n")


    if (x$input$exhaustive & !x$input$cluster){ # exhaustive & non-cluster
      cat("Exhaustive global estimator")
    }


    if (x$input$exhaustive & x$input$cluster){ # exhaustive & cluster
      cat("Exhaustive global estimator for cluster sampling")
    }


    if (!x$input$exhaustive & !x$input$cluster){ # non-exhaustive & non-cluster
      cat("Non-exhaustive global estimator")
    }


    if (!x$input$exhaustive & x$input$cluster){ # non-exhaustive & cluster
      cat("Non-exhaustive global estimator for cluster sampling")
    }

    #   cat("\n", "\n")
    #   cat("Regression Model:")
    #   cat("\n")
    #   print(x$input$formula, showEnv=FALSE)
    cat("\n \n")
  }# end of global-print

} # end of print.threephase



#' @method print confint.smallarea
#' @export
print.confint.smallarea<- function(x, ...){
  # print-method for small area confint-objects:

  cat("\n")
  if(x$adjust.method != "none"){
    cat(paste(x[["level"]]*100, "% Simultaneous Confidence Intervals for ", class(x)[2] ," small area estimation" ,sep = ""))
  } else {
    cat(paste(x[["level"]]*100, "% Confidence Intervals for ", class(x)[2] ," small area estimation" ,sep = ""))
  }
  cat("\n \n")
  print(x[[1]])
  cat("\n")

  if(x$adjust.method != "none"){
    cat(paste("Confidence Interval adjustment by method: ", x$adjust.method))
    cat("\n \n")
  }

}



#' @method print confint.global
#' @export
print.confint.global<- function(x, ...){
  # print-method for global confint-objects:

  cat("\n")
 if(x$adjust.method != "none"){
  cat(paste(x[["level"]]*100, "% Simultaneous Confidence Intervals for ", class(x)[2] ," global estimation" ,sep = ""))
 } else{
  cat(paste(x[["level"]]*100, "% Confidence Intervals for ", class(x)[2] ," global estimation" ,sep = ""))
 }
  cat("\n \n")
  print(x[[1]])
  cat("\n")

  if(x$adjust.method != "none"){
    cat(paste("Confidence Interval adjustment by method: ", x$adjust.method))
    cat("\n \n")
  }

}




