# Error checking functions called in two-phase and three-phase superfunctions:


# ------------------------------------- #
# check for mandatory input parameters: #
# ------------------------------------- #

# for twophase
check.mandatoryInputs<- function(formula, data, phase_id){

  # ------------------------------------- #
  # 1) check if data is of type data.frame:
  if (!is.data.frame(data)) { stop("data must be of class data.frame")}

  # --------------------------------------------------- #
  # 2) check if variables used in formula(s) exist in data:
  if ( prod( all.vars(formula) %in% colnames(data) ) != 1){
    missing.var<- all.vars(formula) [which(all.vars(formula) %in% colnames(data) == FALSE)]
    stop(paste("Variable ", missing.var," used in regression-formula does not exist in data", sep = ""))
  }

  # --------------------------------- #
  # 3) check for input-type "phase_id":

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

 if (!is.numeric( data [[ phase_id[["phase.col"]] ]])){
   stop("'phase_id'-argument must be of type numeric")
 }


 idset<- unique(data [[ phase_id[["phase.col"]] ]])

 if (any(is.na(idset))){
   stop("Invalid NA value(s) found as 'phase_id'")
 }

 if (length(idset) > 2){
   stop(paste(length(idset), "phase ids found. Only 2 different phase ids possible for twophase sampling."))
 }


  # --------------------------------- #
  # 4) check for availability of categorical vars in s2-sample:

    # get column-indices (ind.catvars) of categorical variables in dataset used in model-formula:
    names.catvars<- names(which(sapply(data[,all.vars(formula)[-1]], is.factor)))
    ind.catvars<- which(colnames(data) %in% names.catvars)
    error.list<- list()

    if (length(ind.catvars) > 0) {

      for(i in 1:length(ind.catvars)){

        tab.g<- table(data[,ind.catvars[i]], data[[ phase_id[["phase.col"]] ]])
        ind.s2col<- which(colnames(tab.g) == as.character(phase_id[["terrgrid.id"]]))
        ind.s1col<- which(colnames(tab.g) != as.character(phase_id[["terrgrid.id"]]))

        log.ind<- tab.g[ ,ind.s2col] == 0 & tab.g[ ,ind.s1col] > 0

        if(any(log.ind)){
          ind<- which(tab.g[ ,ind.s2col] == 0 & tab.g[ ,ind.s1col] > 0)

          error.list[i]<- paste("\n Level '", rownames(tab.g)[ind],"' of factor variable '",
                                colnames(data)[ind.catvars[i]],
                                "' existing in s1(s0)- but not in s2 sample. Calculation of coefficient not feasible.",
                                sep="")
        }
      }

      if(length(error.list) > 0){stop(unlist(error.list))}

    }


} #end check.mandatoryInputs


# for threephase
check.mandatoryInputs3p<- function(formula.s0, formula.s1, data, phase_id){

  # ------------------------------------- #
  # 1) check if data is of type data.frame:
  if (!is.data.frame(data)) { stop("data must be of class data.frame")}

  # --------------------------------------------------- #
  # 2) check if variables used in formula(s) exist in data:
  if ( prod( all.vars(formula.s1) %in% colnames(data) ) != 1){
    missing.var<- all.vars(formula.s1) [which(all.vars(formula.s1) %in% colnames(data) == FALSE)]
    stop(paste("Variable '", missing.var,"' used in 'formula.s1' does not exist in data", sep = ""))
  }

  if ( prod( all.vars(formula.s0) %in% colnames(data) ) != 1){
    missing.var<- all.vars(formula.s0) [which(all.vars(formula.s0) %in% colnames(data) == FALSE)]
    stop(paste("Variable '", missing.var,"' used in 'formula.s0' does not exist in data", sep = ""))
  }

  # ------------------------------------- #
  # 3) check if reduced model is a nested in the full model:

  if(!all(paste(c(attr(terms(formula.s0),"term.labels")))  %in%  paste(c(attr(terms(formula.s1),"term.labels"))))) {
    warning("reduced- and full model not nested: variables in 'formula.s0' not a subset of the variables used in 'formula.s1'")
  }

  # --------------------------------- #
  # 4) check for input-type "phase_id":

  # check if "phase.col" is character:
  if (!is.character(phase_id[["phase.col"]])) {
    stop("phase.col argument in input-parameter phase_id must be of type character")
  }
  # check if "phase.col" - name exists in data:
  if(phase_id[["phase.col"]] %in% colnames(data) == FALSE){ stop("specified phase column does not exist in data")}
  # check if terrgrid in "phase_id" does exist in data:
  if(!all(phase_id[["terrgrid.id"]] %in% sapply(unique(data[phase_id[["phase.col"]]]), as.character))) {
    stop("the specified terrestrial grid indicator does not exist in data")
  }
  # check if s1grid.id in "phase_id" does not exist in data:
  if(!all(phase_id[["s2grid.id"]] %in% sapply(unique(data[phase_id[["phase.col"]]]), as.character))) {
    stop("the specified first-phase grid indicator does not exist in data")
  }

  idset<- unique(data [[ phase_id[["phase.col"]] ]])

  if (any(is.na(idset))){
    stop("Invalid NA value(s) found as 'phase_id'")
  }

  if (length(idset) > 3){
    stop(paste(length(idset), "phase ids found. Only 3 different phase ids possible for twophase sampling."))
  }


  # --------------------------------- #
  # 5) check for availability of categorical vars in s2-sample:

  # get column-indices (ind.catvars) of categorical variables in dataset used in full-model formula:
  names.catvars<- names(which(sapply(data[,all.vars(formula.s1)[-1]], is.factor)))
  ind.catvars<- which(colnames(data) %in% names.catvars)
  error.list<- list()

  if (length(ind.catvars) > 0) {

    for(i in 1:length(ind.catvars)){

      tab.g<- table(data[,ind.catvars[i]], data[[ phase_id[["phase.col"]] ]])
      ind.s2col<- which(colnames(tab.g) == as.character(phase_id[["terrgrid.id"]]))
      ind.s1col<- which(colnames(tab.g) != as.character(phase_id[["terrgrid.id"]]))

      log.ind<- tab.g[ ,ind.s2col] == 0 & tab.g[ ,ind.s1col] > 0

      if(any(log.ind)){
        ind<- which(tab.g[ ,ind.s2col] == 0 & tab.g[ ,ind.s1col] > 0)

        error.list[i]<- paste("\n Level '", rownames(tab.g)[ind],"' of factor variable '",
                              colnames(data)[ind.catvars[i]],
                              "' existing in s1(s0)- but not in s2 sample. Calculation of coefficient not feasible.",
                              sep="")
      }
    }

    if(length(error.list) > 0){stop(unlist(error.list))}

  }

} #end check.mandatoryInputs3p


# ---------------------------------------------------------------------------- #
# ---------------------------------------------------------------------------- #


# -------------------------------- #
# check optional input parameters: #
# -------------------------------- #


# ---------------------------------------------------------------------------- #
# check cluster Input:

check.clusterInput<- function(data, cluster){

  # check if input-type "cluster" is character / if it exists in data:
  if(!is.character(cluster)){ stop("input-parameter cluster must be of type character")}
  if(!(cluster %in% colnames(data))){ stop(paste("specified cluster-column '",cluster,"' does not exist in data", sep = ""))}
}


# ---------------------------------------------------------------------------- #
# check boundary weight Input:

check.boundary_weightsInput<- function(data, boundary_weights){

  # check if boundary_weights input is of type character:
  if(!is.character(boundary_weights)){
    stop("input-parameter 'boundary_weights' must be of type 'character' describing the variable name")
  }

  # check if boundary_weights input is of type character:
  if(any(data[[boundary_weights]] < 0) | any(data[[boundary_weights]] > 1)){
    stop("invalid values for 'boundary_weights' used. Only values between 0 and 1 applicable")
  }

  # check if boundary_weights-column does exist in dataset:
  if(boundary_weights %in% colnames(data) == FALSE){
    stop(paste("specified boundary_weights-column '",boundary_weights,"' does not exist in data", sep = ""))
  }
}


# ---------------------------------------------------------------------------- #
# check small area Input:

check.smallareaInput<- function(data, small_area){

  # check if "sa.col" is character:

  if (!is.logical(small_area[["unbiased"]])) {
    stop("unbiased argument in input-parameter small_area must be of type logical")
  }


  if (!is.character(small_area[["sa.col"]])) {
    stop("sa.col argument in input-parameter small_area must be of type character")
  }

  # check if "sa.col" - name exists in data:
  if(small_area[["sa.col"]] %in% colnames(data) == FALSE){ stop("specified small_area column does not exist in data")}

  # check if any specified areas in "areas" does not exist in data:
  if(!all(small_area[["areas"]] %in% sapply(unique(data[small_area[["sa.col"]]]), as.character))) {
    nomatch.sa<- small_area[["areas"]] [- which(small_area[["areas"]] %in% sapply(unique(data[small_area[["sa.col"]]]), as.character))]
    stop(paste("the specified small areas ", paste0(nomatch.sa, collapse = ","), " do not exist in data",sep = ""))
  }
}


# ---------------------------------------------------------------------------- #
# check exhaustive Input:

# for twophase
check.exhaustive<- function(formula, boundary_weights, exhaustive, data){

  # check if exhaustive AND boundary-weights are not given at the same time:
  if(!is.na(boundary_weights)){
    stop("boundary_weights and exhaustive true means cannot be given at the same time,
          exhaustive true means must already include boundary weight adjustment")
  }

  # check if exhaustive means are provided for every auxiliary variable:
  if (!all(!is.na(exhaustive)) & ncol(data.frame(exhaustive)) != ncol(model.matrix(formula, data))){
    stop("number of true auxiliary means and number of auxiliary variables used in formula differ")
  }

} # end check.exhaustive for two-phase


# for threephase
check.exhaustive3p<- function(formula.s0, exhaustive, data){

  # check if exhaustive means are provided for every auxiliary variable in formula.s0:
  if (!all(!is.na(exhaustive)) & ncol(data.frame(exhaustive)) != ncol(model.matrix(formula.s0, data))){
    stop("number of true auxiliary means and number of auxiliary variables used in 'formula.s0' differ")
  }

} # end check.exhaustive for three-phase


# ---------------------------------------------------------------------------- #

# check exhaustive Input for cluster:

# for twophase
check.exhaustive.cluster<- function(formula, boundary_weights, exhaustive, data){

  # boundary-weights can be given to calculate the appropriate WEIGHTED cluster-means!

  # check if exhaustive means are provided for every auxiliary variable:
  if (!all(!is.na(exhaustive)) & ncol(data.frame(exhaustive)) != ncol(model.matrix(formula, data))){
    stop("number of true auxiliary means and number of auxiliary variables used in formula differ")
  }

}

# for threephase
check.exhaustive.cluster3p<- function(formula.s0, exhaustive, data){

  # boundary-weights can be given to calculate the appropriate WEIGHTED cluster-means!

  # check if exhaustive means are provided for every auxiliary variable:
  if (!all(!is.na(exhaustive)) & ncol(data.frame(exhaustive)) != ncol(model.matrix(formula.s0, data))){
    stop("number of true auxiliary means and number of auxiliary variables used in 'formula.s0' differ")
  }

}

# ---------------------------------------------------------------------------- #


