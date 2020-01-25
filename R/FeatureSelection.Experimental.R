# FS.greedy.heuristic.reduct.RST.positive <- function(decision.table,
#                                                     attrDescriptions = attr(decision.table, "desc.attrs"),
#                                                     decisionIdx = attr(decision.table, "decision.attr"),
#                                                     qualityF = X.nOfConflicts, nAttrs = NULL,
#                                                     epsilon = 0.0, inconsistentDecisionTable = NULL)  {
#   toRmVec = decisionIdx
#   attrIdxVec = (1:ncol(decision.table))[-toRmVec]
#   
#   if (!is.null(nAttrs)) {
#     if (nAttrs == 0 || nAttrs > ncol(decision.table) - 1) {
#       stop("There is something wrong with data (too little attributes?) or the parameter nAttrs has a wrong value.")
#     }
#     tmpAttrSub = sample(attrIdxVec, min(nAttrs, length(attrIdxVec)))
#   }  else tmpAttrSub = attrIdxVec
#   
#   if (epsilon >= 1 || epsilon < 0) {
#     stop("Wrong value of the parameter epsilon. It must be within [0,1) interval.")
#   }
#   
#   INDrelation = list(1:nrow(decision.table))
#   INDsizes = nrow(decision.table)
#   decisionChaos = compute_chaos(INDrelation, decision.table[[decisionIdx]],
#                                 attrDescriptions[[decisionIdx]])
#   decisionChaos = qualityF(decisionChaos[[1]])
#   attrScoresVec = mapply(qualityGain, decision.table[tmpAttrSub], attrDescriptions[tmpAttrSub],
#                          MoreArgs = list(decisionVec = decision.table[[decisionIdx]],
#                                          uniqueDecisions = attrDescriptions[[decisionIdx]],
#                                          INDclassesList = INDrelation,
#                                          INDclassesSizes = INDsizes,
#                                          baseChaos = decisionChaos,
#                                          chaosFunction = qualityF),
#                          SIMPLIFY = TRUE, USE.NAMES = FALSE)
#   tmpBestIdx = which.max(attrScoresVec)
#   
#   selectedAttrIdxVec  = tmpAttrSub[tmpBestIdx]
#   INDrelation = compute_indiscernibility(INDrelation,
#                                          as.character(decision.table[[selectedAttrIdxVec]]),
#                                          attrDescriptions[[selectedAttrIdxVec]])
#   attrIdxVec = (1:ncol(decision.table))[-c(selectedAttrIdxVec, toRmVec)]
#   
#   if(is.null(inconsistentDecisionTable)) {
#     if(sum(duplicated(decision.table)) == sum(duplicated(decision.table[-decisionIdx]))) {
#       inconsistentDecisionTable = FALSE
#     } else {
#       inconsistentDecisionTable = TRUE
#     }
#   }
#   
#   endFlag = FALSE
#   iteration = 1
#   if(inconsistentDecisionTable) {
#     tmpIND = split(1:nrow(decision.table),
#                    do.call(paste, decision.table[-decisionIdx]), drop = TRUE)
#     totalChaos = compute_chaos(tmpIND,
#                                as.character(decision.table[[decisionIdx]]),
#                                attrDescriptions[[decisionIdx]])
#     totalChaos = sum(sapply(totalChaos, qualityF) * sapply(tmpIND, length) / nrow(decision.table))
#     rm(tmpIND)
#   } else totalChaos = 0
#   totalDependencyInData = decisionChaos - totalChaos - 10^(-16)
#   approxThereshold = (1 - epsilon)*totalDependencyInData
#   
#   while (!endFlag) {
#     contingencyTabs = lapply(INDrelation,
#                              function(x,y) table(y[x]),
#                              decision.table[[decisionIdx]])
#     chaosVec = sapply(contingencyTabs, qualityF)
#     if(any(chaosVec == 0)) {
#       tmpIdx = chaosVec != 0
#       INDrelation = INDrelation[tmpIdx]
#       contingencyTabs = contingencyTabs[tmpIdx]
#       chaosVec = chaosVec[tmpIdx]
#       rm(tmpIdx)
#     }
#     sumsVec = sapply(contingencyTabs, sum)
#     if(length(INDrelation) > 0) tmpChaos = sum(chaosVec*(sumsVec/nrow(decision.table)))
#     else tmpChaos = 0
#     tmpDependencyInData = decisionChaos - tmpChaos
#     rm(contingencyTabs, sumsVec, chaosVec)
#     
#     if (approxThereshold <= tmpDependencyInData) {
#       endFlag = TRUE
#     }	else {
#       if (!is.null(nAttrs)) {
#         tmpAttrSub = sample(attrIdxVec, min(nAttrs, length(attrIdxVec)))
#       }	else tmpAttrSub = attrIdxVec
#       INDsizes = sapply(INDrelation, length)
#       attrScoresVec = mapply(qualityGain, decision.table[tmpAttrSub], attrDescriptions[tmpAttrSub],
#                              MoreArgs = list(decisionVec = decision.table[[decisionIdx]],
#                                              uniqueDecisions = attrDescriptions[[decisionIdx]],
#                                              INDclasses = INDrelation,
#                                              INDclassesSizes = INDsizes,
#                                              baseChaos = tmpChaos,
#                                              chaosFunction = qualityF),
#                              SIMPLIFY = TRUE, USE.NAMES = FALSE)
#       tmpBestIdx = which.max(attrScoresVec)
#       selectedAttrIdxVec[iteration + 1] = tmpAttrSub[tmpBestIdx]
#       INDrelation = compute_indiscernibility(INDrelation,
#                                              as.character(decision.table[[tmpAttrSub[tmpBestIdx]]]),
#                                              attrDescriptions[[tmpAttrSub[tmpBestIdx]]])
#       if(length(INDrelation) == 0) endFlag = TRUE
#       
#       attrIdxVec = (1:ncol(decision.table))[-c(selectedAttrIdxVec, toRmVec)]
#       iteration = iteration + 1
#     }
#   }
#   
#   if (iteration > 1) {
#     endFlag = FALSE
#     iteration = iteration - 1
#     while (!endFlag) {
#       clsContingencyTab = as.matrix(table(do.call(paste, decision.table[selectedAttrIdxVec[-iteration]]), decision.table[[decisionIdx]]))
#       tmpChaos = sum(apply(clsContingencyTab, 1, qualityF)*(rowSums(clsContingencyTab)/nrow(decision.table)))
#       tmpDependencyInData = decisionChaos - tmpChaos
#       if (approxThereshold <= tmpDependencyInData) {
#         selectedAttrIdxVec = selectedAttrIdxVec[-iteration]
#       }
#       iteration = iteration - 1
#       if (iteration == 0) endFlag = TRUE
#     }
#   }
#   
#   reduct <- selectedAttrIdxVec[order(selectedAttrIdxVec)]
#   names(reduct) <- colnames(decision.table)[reduct]
#   mod <- list(reduct = reduct, type.method = "greedy.heuristic",
#               epsilon = epsilon,
#               type.task = "feature selection", model = "RST")
#   
#   class.mod <- ObjectFactory(mod, classname = "FeatureSubset")
#   return(class.mod)
# }

# qualityGain.positive <- function(vec, uniqueValues, decisionVec, uniqueDecisions,
#                                  INDclassesList, INDclassesSizes, baseChaos, chaosFunction = X.gini)  {
#   
#   classCounts = .C("computeIndiscernibilityAndChaos",
#                    INDclasses = as.integer(unlist(INDclassesList)),
#                    INDsizes = as.integer(INDclassesSizes),
#                    NOfINDClasses = as.integer(length(INDclassesSizes)),
#                    attrValues = as.integer(vec),
#                    NOfAttrValues = as.integer(length(uniqueValues)),
#                    decValues = as.integer(decisionVec),
#                    NOfDecs = as.integer(length(uniqueDecisions)),
#                    output = as.integer(rep(0,length(INDclassesList)*length(uniqueValues)*length(uniqueDecisions))), PACKAGE="RoughSets")
#   
#   classCounts = matrix(classCounts$output,
#                        nrow = length(INDclassesList)*length(uniqueValues),
#                        ncol = length(uniqueDecisions), byrow=TRUE)
#   newINDsizes = rowSums(classCounts)
#   tmpInd = newINDsizes > 1
#   if(sum(tmpInd) > 0) {
#     remainingChaos = apply(classCounts[tmpInd, ,drop=FALSE], 1, chaosFunction)
#     remainingChaos = sum(remainingChaos * newINDsizes[tmpInd]) / length(decisionVec)
#   } else remainingChaos = 0
#   
#   return(as.numeric(baseChaos - remainingChaos))
# }



#' This function implements a positive decision-relative discernibility matrix. This notion
#' was proposed by (Sikora et al.) as a middle-step in many RST algorithms for computaion of reducts, 
#' discretization and rule induction in a case when the discernibility of objects from the positive class 
#' by positive attribute values is more desirable than by the negative ones. The implementation currently 
#' works only for binary decision system (all attributes, including the decision must be binary
#' and the positive value is marked by "1").
#'
#' @title Computation of a positive decision-relative discernibility matrix based on the rough set theory
#' @author Andrzej Janusz and Dominik Slezak
#'
#' @param decision.table an object inheriting from the \code{DecisionTable} class, which represents a decision system. 
#'        See \code{\link{SF.asDecisionTable}}.
#' @param return.matrix a logical value determining whether the discernibility matrix should be retunred in the output. 
#'        If it is set to FALSE (the default) only a list containing unique clauses from the CNF representation 
#'        of the discernibility function is returned.
#' @param attach.data a logical indicating whether the original decision table should be attached as 
#'        an additional element of the resulting list named as \code{dec.table}.
#' 
#' @return An object of a class \code{DiscernibilityMatrix} which has the following components: 
#' \itemize{
#' \item \code{disc.mat}: the decision-relative discernibility matrix which for pairs of objects from different 
#'       decision classes stores names of attributes which can be used to discern between them. Only pairs of 
#'       objects from different decision classes are considered. For other pairs the \code{disc.mat} contains
#'       \code{NA} values. Moreover, since the classical discernibility matrix is symmetric only the pairs 
#'       from the lower triangular part are considered.
#' \item \code{disc.list}: a list containing unique clauses from the CNF representation of the discernibility 
#'       function,
#' \item \code{dec.table}: an object of a class \code{DecisionTable}, which was used to compute the
#'       discernibility matrix,
#' \item \code{discernibility.type}: a type of discernibility relation used in the computations,
#' \item \code{type.model}: a character vector identifying the type of model which was used. 
#'                In this case, it is \code{"RST"} which means the rough set theory.
#' }
#' 
#' @seealso \code{\link{BC.IND.relation.RST}}, \code{\link{BC.LU.approximation.RST}}, \code{\link{BC.LU.approximation.FRST}}
#'          and \code{\link{BC.discernibility.mat.FRST}}
#' 
#' @examples
#' ###############################################################################
#' ## Example 1: Constructing the positive decision-relative discernibility matrix
#' ###############################################################################
#' data(RoughSetData)
#' binary.dt <- RoughSetData$binary.dt
#'
#' ## building the decision-relation discernibility matrix
#' disc.matrix <- BC.discernibility.positive.mat.RST(binary.dt, return.matrix = TRUE)
#' disc.matrix$disc.mat
#' 
#' ## compute all classical reducts
#' classic.reducts <- FS.all.reducts.computation(BC.discernibility.mat.RST(binary.dt))
#' head(classic.reducts$decision.reduct)
#' cat("A total number of reducts found: ", 
#'     length(classic.reducts$decision.reduct), "\n", sep = "")
#' classic.reducts$core
#' 
#' ## compute all positive reducts
#' positive.reducts <- FS.all.reducts.computation(disc.matrix)
#' head(positive.reducts$decision.reduct)
#' cat("A total number of positive reducts found: ", 
#'     length(positive.reducts$decision.reduct), "\n", sep = "")
#' print("The core:")
#' positive.reducts$core
#'
#' @references
#' TO BE ADDED
#' 
#' @export
BC.discernibility.positive.mat.RST <- function(decision.table, 
                                               return.matrix = FALSE, attach.data = FALSE){
  
  if(!inherits(decision.table, "DecisionTable")) {
    stop("Provided data should inherit from the \'DecisionTable\' class.")
  }
  
  ## get the data
  objects <- decision.table
  desc.attrs <- attr(decision.table, "desc.attrs")
  nominal.att <- attr(decision.table, "nominal.attrs")
  if (is.null(attr(decision.table, "decision.attr"))){
    stop("A decision attribute is not indicated.")
  }	else {
    decision.attr <- attr(decision.table, "decision.attr")
    ## check if index of decision attribute is not in last index, then change their position
    if (decision.attr != ncol(objects)){
      objects <- cbind(objects[, -decision.attr, drop = FALSE], objects[, decision.attr, drop = FALSE])
      nominal.att <- c(nominal.att[-decision.attr], nominal.att[decision.attr])
    }
  }
  if(!all(sapply(desc.attrs, length) == 2) | !all(sapply(desc.attrs, function(x) "1" %in% x))) {
    stop("This is not a binary decision table - make sure that '1' is among possible values of each attribute")
  }
  
  num.object <- nrow(objects)
  names.attr <- t(colnames(objects)[-decision.attr])
  
  ## initialize the discernibility matrix
  disc.mat <- array(list(NA), dim = c(num.object, num.object, 1))
  
  decVector = as.character(objects[, ncol(objects)])
  dataMatrix = objects[, -ncol(objects)] == '1'
  
  ## construct the positive discernibility matrix 
  ones_idxs = which(decVector == "1")
  zero_idxs = (1:length(decVector))[-ones_idxs]
  if(length(zero_idxs) > 0)  {
    for (i in ones_idxs) {
      for (j in zero_idxs) {
        disc.attr <- names.attr[dataMatrix[i,] & !dataMatrix[j,]]
        if(length(disc.attr) == 0) {
          disc.attr <- names.attr[!dataMatrix[i,] & dataMatrix[j,]]
        }
        disc.mat[j, i, 1] <- list(disc.attr)
      }
    }
  }
  disc.mat = as.data.frame(disc.mat)
  disc.list = unique(do.call(c, disc.mat))[-1]
  
  ## build class
  if (return.matrix){
    discernibilityMatrix = list(disc.mat = disc.mat, disc.list = disc.list, 
                                names.attr = colnames(decision.table), type.discernibility = "RST", type.model = "RST")
  }
  else {
    discernibilityMatrix = list(disc.list = disc.list, 
                                names.attr = colnames(decision.table), type.discernibility = "RST", type.model = "RST")
  }
  if(attach.data) {
    discernibilityMatrix = c(list(dec.table = decision.table), discernibilityMatrix)
  }
  
  discernibilityMatrix = ObjectFactory(discernibilityMatrix, classname = "DiscernibilityMatrix")
  return(discernibilityMatrix)
}
