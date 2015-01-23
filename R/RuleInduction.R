#############################################################################
#
#  This file is a part of the R package "RoughSets".
#
#  Author: Lala Septem Riza and Andrzej Janusz
#  Supervisors: Chris Cornelis, Francisco Herrera, Dominik Slezak and Jose Manuel Benitez
#  Copyright (c):
#       DiCITS Lab, Sci2s group, DECSAI, University of Granada and
#       Institute of Mathematics, University of Warsaw
#
#  This package is free software: you can redistribute it and/or modify it under
#  the terms of the GNU General Public License as published by the Free Software
#  Foundation, either version 2 of the License, or (at your option) any later version.
#
#  This package is distributed in the hope that it will be useful, but WITHOUT
#  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
#  A PARTICULAR PURPOSE. See the GNU General Public License for more details.
#
#############################################################################
#' It is a function for generating rules based on hybrid fuzzy-rough rule induction and feature selection.
#' It allows for classification and regression tasks.
#' 
#' It was proposed by (Jensen et al, 2009) attempting to combine rule induction and feature selection
#' at the same time. Basically this algorithm inserts some steps to generate rules
#' into the fuzzy QuickReduct algorithm (see \code{\link{FS.quickreduct.FRST}}. 
#' Furthermore, by introducing the degree of coverage, this algorithm selects proper rules. 
#'
#' This function allows not only for classification but also for regression problems. After obtaining the rules,
#' predicting can be done by calling \code{predict} or \code{\link{predict.RuleSetFRST}}.
#' Additionally, to get better representation we can execute \code{\link{summary}}.
#'
#' @title Hybrid fuzzy-rough rule and induction and feature selection
#' 
#' @param decision.table a \code{"DecisionTable"} class representing the decision table. See \code{\link{SF.asDecisionTable}}. 
#' @param control a list of other parameters which consist of 
#'         \itemize{
#'         \item \code{type.aggregation} a list representing the type of aggregation. The default value is \code{type.aggregation = c("t.tnorm", "lukasiewicz")}.
#' 
#'                             See \code{\link{BC.IND.relation.FRST}}.
#'         \item \code{type.relation} the type of indiscernibility relation. The default value is \code{type.relation = c("tolerance", "eq.3")}.
#'                             See \code{\link{BC.IND.relation.FRST}}.
#'         \item \code{t.implicator} the type of implication function. The default value is \code{"lukasiewicz"}. 
#'                             See BC.LU.approximation.FRST.
#'         }
#' @seealso \code{\link{RI.indiscernibilityBasedRules.RST}}, \code{\link{predict.RuleSetFRST}}, and \code{\link{RI.GFRS.FRST}}.
#' @return A class \code{"RuleSetFRST"} which has similar components as \code{\link{RI.GFRS.FRST}} 
#'         but in this case the \code{threshold} component is not included.
#' @references
#' R. Jensen, C. Cornelis, and Q. Shen, "Hybrid Fuzzy-rough Rule Induction and Feature Selection", 
#' in: IEEE International Conference on Fuzzy Systems (FUZZ-IEEE), p. 1151 - 1156 (2009).
#' 
#' @examples
#' ###########################################################
#' ## Example 1: Regression problem
#' ###########################################################
#' data(RoughSetData)
#' decision.table <- RoughSetData$housing7.dt
#'
#' control <- list(type.aggregation = c("t.tnorm", "lukasiewicz"), type.relation = 
#'                 c("tolerance", "eq.3"), t.implicator = "lukasiewicz")
#' res.1 <- RI.hybridFS.FRST(decision.table, control)
#'
#' ###########################################################
#' ## Example 2: Classification problem
#' ##############################################################
#' data(RoughSetData)
#' decision.table <- RoughSetData$pima7.dt
#' 
#' control <- list(type.aggregation = c("t.tnorm", "lukasiewicz"), type.relation =  
#'                 c("tolerance", "eq.3"), t.implicator = "lukasiewicz")
#' res.2 <- RI.hybridFS.FRST(decision.table, control)
#'
#' @export
RI.hybridFS.FRST <- function(decision.table, control = list()){	

	return(quickreduct.alg(decision.table, type.method = "hybrid.rule.induction", control = control))
}

#' It is a function generating rules in classification tasks using the fuzzy variable precision rough sets (FVPRS) approach (see \code{\link{BC.LU.approximation.FRST}}). 
#' 
#' The method proposed by (Zhao, 2010) consists of three steps as follows.
#' First, it builds a general lower approximation that is able to deal with misclassification and perturbation.
#' In this case, the fuzzy variable precision rough sets (FVPRS)
#' is used to calculate the lower approximation (see \code{\link{BC.LU.approximation.FRST}}). 
#' Secondly, a discernibility matrix considering a consistence degree is constructed for obtaining rules.
#' The details about the matrix can be seen in \code{\link{BC.discernibility.mat.FRST}}.
#' Then, we calculate attribute value reduction of every object and perform near-minimal rule set. 
#' The final step is to construct rules considering the consistence degree of associated objects. 
#' 
#' It should be noted that this function only allows classification problems. After obtaining the rules,
#' predicting can be done by calling \code{predict} or \code{\link{predict.RuleSetFRST}}.
#' Additionally, to get better representation we can execute \code{\link{summary}}.
#' 
#' @title Generalized fuzzy rough set rule induction based on FRST
#' 
#' @param decision.table a \code{"DecisionTable"} class representing the decision table. See \code{\link{SF.asDecisionTable}}. 
#' @param control a list of other parameters which consist of 
#'         \itemize{
#'         \item \code{alpha.precision}: a numeric value representing variable precision of FVPRS. 
#'                            The default value is 0.05. See \code{\link{BC.LU.approximation.FRST}}.
#'         \item \code{type.aggregation}: a list of a type of aggregations. The default value is \code{type.aggregation = c("t.tnorm", "lukasiewicz")}. 
#'
#'                             See \code{\link{BC.IND.relation.FRST}}.
#'         \item \code{type.relation}: a type of indiscernibility relation. The default value is \code{type.relation = c("tolerance", "eq.1")}.
#'                             See \code{\link{BC.IND.relation.FRST}}.
#'         \item \code{t.implicator}: a type of implication function. The default value is \code{"lukasiewicz"}. 
#'                             See \code{\link{BC.LU.approximation.FRST}}.
#'         }
#' @seealso \code{\link{RI.indiscernibilityBasedRules.RST}}, \code{\link{predict.RuleSetFRST}}, and \code{\link{RI.hybridFS.FRST}}.
#' @return A class \code{"RuleSetFRST"} which has following components:
#' \itemize{
#' \item \code{rules}: It is a list containing two elements which are \code{rules} and \code{threshold}. 
#'              The \code{rules} represent knowledge in data set that can be expressed as an IF-THEN form.
#'              For example, we got the rule as follows: \code{90 8 2} and its colnames: \code{pres}, \code{preg}, and \code{class}. 
#'              It refers to the following rule: \code{IF pres is about 90 and preg is about 8 THEN class is 2}.
#'              In other words, while the last column represents the consequent part, 
#'              the rest expresses the antecedent part.
#'              The second part of this object is \code{threshold} representing a value used to predict new data. 
#'              In order to change IF-THEN form, we can use \code{\link{summary}}.
#' \item \code{type.model}: It is the type of the theory whether \code{"FRST"} or \code{"RST"}. In this case, it is \code{FRST}.
#' \item \code{type.method}: It is the considered method. In this case, it is \code{RI.GFRS.FRST}.
#' \item \code{type.task}: It is the type of task. In this case, it is \code{"classification"}.
#' \item \code{t.similariy}: It is the type of similarity equation. See \code{\link{BC.IND.relation.FRST}}.
#' \item \code{t.tnorm}: It is the type of triangular operator. See \code{\link{BC.IND.relation.FRST}}.
#' \item \code{variance.data}: It represents the variance of the data set. It has \code{NA} values when the associated attributes are nominal values.
#' \item \code{range.data}: It represents the range of the data set. It has \code{NA} values when the associated attributes are nominal values.
#' \item \code{antecedent.attr}: It is a list of attribute names involved in the antecedent part.
#' \item \code{consequent.attr}: It is the attribute in the consequent part.
#' \item \code{nominal.att}: It is a list of boolean that represent whether a attribute is nominal or not. 
#' } 
#' @references
#' S. Y. Zhao, E. C. C. Tsang, D. G. Chen, and X. Z. Wang, "Building a Rule-based
#' Classifier -- A Fuzzy-rough Set Approach",
#' IEEE Trans. on Knowledge and Data Engineering, vol. 22, no. 5, p. 624 - 638 (2010).
#' 
#' @examples
#' ###########################################################
#' ## Example 
#' ##############################################################
#' data(RoughSetData)
#' decision.table <- RoughSetData$pima7.dt							 
#'	
#' control <- list(alpha.precision = 0.01, type.aggregation = c("t.tnorm", "lukasiewicz"), 
#'                 type.relation = c("tolerance", "eq.3"), t.implicator = "lukasiewicz")						 
#' rules <- RI.GFRS.FRST(decision.table, control)
#'
#' @export
RI.GFRS.FRST <- function(decision.table, control = list()){	

	if (is.null(attr(decision.table, "decision.attr"))){
		stop("A decision attribute is not indicated.")
	}
	
	if (attr(decision.table, "nominal.attrs")[attr(decision.table, "decision.attr")] == FALSE){
		stop("The decision attribute must be nominal values")
	}
	
	## set default values of all parameters
	control <- setDefaultParametersIfMissing(control, list(alpha.precision = 0,05, type.aggregation = c("t.tnorm", "lukasiewicz"), type.relation = c("tolerance", "eq.1"), 
														t.implicator = "lukasiewicz"))   
										
	## get the data
	objects <- decision.table
	desc.attrs <- attr(decision.table, "desc.attrs")
	nominal.att <- attr(decision.table, "nominal.attrs")
	num.att <- ncol(objects)
	indx.univ <- c(seq(1, nrow(objects)))
	num.object <- nrow(objects)
	names.attr <- t(colnames(objects))
	alpha.precision <- control$alpha.precision
	type.aggregation <- ch.typeAggregation(control$type.aggregation)
	type.relation <- control$type.relation
	t.similarity <- type.relation[2]
	t.implicator <- control$t.implicator
	t.tnorm <- type.aggregation[[2]]

	###############################
	#### build decision-relative discernibility matrix
	##################################
	res.temp <- build.discMatrix.FRST(decision.table = decision.table, type.discernibility = "alpha.red", num.red = "all", alpha.precision = alpha.precision, 
	                                  type.relation = type.relation, t.implicator = t.implicator, type.aggregation = type.aggregation, show.discernibilityMatrix = TRUE)
	disc.mat <- res.temp$disc.mat
	
	###################################
	#### get one near optimal attribute value reduction for EACH OBJECTS
	####################################	
	att.val.red <- c()
	for (i in 1 : ncol(disc.mat)){
		att.val.core <- c()
		disc.list <- disc.mat[, i]
		
		## delete redundant of discernibility matrix
		disc.list <- unique(disc.list)		
		disc.list <- disc.list[which(!is.na(disc.list))]		
		att.val.core <- disc.list[which(lapply(disc.list, length) == 1)]
		
		if (length(att.val.core) != 0){
			## delete element containing core (for each core) 
			disc.list <- disc.list[which(lapply(disc.list, length) > 1)]
			
			for (ii in 1 : length(att.val.core)){
				disc.list <- disc.list[which(lapply(disc.list, function(x){att.val.core[ii] %in% x}) == FALSE)] 	
			}
		}
		else {
			att.val.core <- disc.list[which.min(lapply(disc.list, length))]
			disc.list <- disc.list[which(lapply(disc.list, function(x){att.val.core %in% x}) == FALSE)] 	
		}
	
		## initialization
		temp.red <- c()
		while (length(disc.list) != 0){
		
			## get the max of freq.
			new.red <- names(which.max(table(unlist(disc.list))))
	
			## add new.red into reducts
			temp.red <- c(temp.red, new.red)

			## delete element containing new.red
			disc.list <- disc.list[which(lapply(disc.list, function(x){new.red %in% x}) == FALSE)] 
		}
		
		temp.red <- c(temp.red, unlist(att.val.core))
		
		## collect attribute value reduction
		att.val.red <- append(att.val.red, list(unique(temp.red)))
	}
	
	####################################
	#### Procedure for calculating degree of consistence (which is positive region)
	####################################	
    #### Perform fuzzy indiscernibility relation ####
	attributes <- c(seq(1, (num.att - 1)))
    control.ind <- list(type.aggregation = type.aggregation, type.relation = type.relation)
    IND <- BC.IND.relation.FRST(decision.table, attributes = attributes, control = control.ind)	
	IND.decAttr <- BC.IND.relation.FRST(decision.table, attributes = c(num.att), control = control.ind)	
	control <- list(t.implicator = t.implicator, t.tnorm = t.tnorm, alpha = alpha.precision)
	
	#### perform lower approximation
    FRST.fvprs <- BC.LU.approximation.FRST(decision.table, IND, IND.decAttr, 
               type.LU = "fvprs", control = control)
		
	#### Determine fuzzy regions ####
	cons.deg <- BC.positive.reg.FRST(decision.table, FRST.fvprs)$positive.freg
	
	###################################
	#### near-minimal rule set
	######################################
	## calculate range of data
	res <- cal.var.range(objects, attributes = seq(1, (num.att - 1)), nominal.att)
	range.dt <- res$range.data
	variance.dt <- res$variance.data
	
	all.rules <- att.val.red
	min.rules <- c()
	object.rules <- list()
	exit <- FALSE
	
	while (exit == FALSE && length(all.rules) >= 1){
		cover.rule <- c()
		num.rl <- matrix()
		indx.seq <- seq(1, length(all.rules))		
		indx.logical.0 <- which(lapply(all.rules, length) == 0)
		if (length(indx.logical.0) != 0){
			indx.seq <- indx.seq[-indx.logical.0] 
		}
		
		for (j in indx.seq){
			cov.rl <- c()
			counter <- 0
			attrs <- all.rules[[j]]
			
			## get object associated with attributes in attrs
			rl.obj <- objects[j, sapply(attrs, function(x) which(colnames(objects) == x))]
			
			## get type of attributes according to selected features
			state.att <- nominal.att[sapply(attrs, function(x) which(colnames(objects) == x))]
			
			## calculate range of objects which is associated with attributes in attrs
			range.data.temp <- matrix(range.dt[1, sapply(attrs, function(x) which(colnames(objects) == x))], nrow = 1) 
			variance.data.temp <- matrix(variance.dt[1, sapply(attrs, function(x) which(colnames(objects) == x))], nrow = 1) 
			
			for (k in indx.seq){
				if (j != k){
					## get object associated with attrs to be compared
					temp.rl.obj <- objects[k, sapply(attrs, function(x) which(colnames(objects) == x))]
					
					## calculate the similarity between two objects based on type of similarity
					IND.rel <- matrix(nrow = 1, ncol = length(state.att))
					for (ii in 1 : length(state.att)){
						if (state.att[ii] == TRUE){
							## overwrite type of similarity when crisp
							t.similarity <- "boolean"
						}						
						IND.rel[1, ii] <- unlist(similarity.equation(rl.obj[ii], temp.rl.obj[ii], t.similarity, delta = 0.2, 
						                         range.data = range.data.temp[1, ii], variance.data = variance.data.temp[1, ii]))
						
					}

					##	check whether object k is covered by j				
					if ((1 - min(IND.rel)) < cons.deg[j]){
						cov.rl <- c(cov.rl, k)
						counter <- counter + 1
					}
				}
			}	

			## list that shows a collection of all covered rules
			cover.rule <- append(cover.rule, list(cov.rl))
			num.rl[j] <- counter
		}	
			
		##	Initialization matrix list.obj showing how many objects are covered. 
		if (length(num.rl) == num.object){
			list.obj <- rbind(seq(1, num.object), num.rl)
		}
			
		## get max coverage
		indx.max.cov <- which.max(num.rl)
		
		## get one rule that covers max objects
		selected.obj <- list.obj[1, indx.max.cov]	
		min.rules <- append(min.rules, list(all.rules[[indx.max.cov]]))
		
		## update object corresponding with attributes in rules
		object.rules <- append(object.rules, list(selected.obj))
	
		## update all.rules
		all.rules <- all.rules[-c(cover.rule[[indx.max.cov]], indx.max.cov)]
		list.obj <- list.obj[, -c(cover.rule[[indx.max.cov]], indx.max.cov), drop = FALSE]
		
		## when remaining rule = 1
		if (length(all.rules) == 1){
			min.rules <- append(min.rules, all.rules) 
			object.rules <- append(object.rules, list(list.obj[1, 1]))
			exit <- TRUE
		}
		
	}
	indx.logical.0 <- which(lapply(min.rules, length) == 0)
	if (length(indx.logical.0) != 0){
		min.rules <- min.rules[-indx.logical.0]
		object.rules <- object.rules[-indx.logical.0]
	}
	
	rules <- list(attributes = min.rules, objects = object.rules, threshold = cons.deg)
	return(std.rules(rules, type.method = "RI.GFRS.FRST", decision.table, t.similarity, t.tnorm))
}

#' It is a function generating rule induction based on indiscernibility classes.  
#' 
#' Basically, this function uses the output of feature selection step as input data. 
#' So, before calling this function, we execute a particular function implementing feature selection, 
#'
#' e.g. \code{\link{FS.permutation.heuristic.reduct.RST}}. 
#'
#' It should be noted that this function only allows classification problems. After obtaining the rules,
#' predicting can be done by calling \code{predict} or \code{\link{predict.RuleSetRST}}.
#'
#' @title Rule induction based on RST
#' 
#' @param decision.table a \code{"DecisionTable"} class representing the decision table. See \code{\link{SF.asDecisionTable}}. 
#' @param feature.set a \code{"FeatureSubset"} class which is produced by feature selection functions based on RST 
#'                  e.g. \code{\link{FS.permutation.heuristic.reduct.RST}}, etc. 
#'
#'                  See \code{\link{FS.reduct.computation}}, \code{\link{FS.feature.subset.computation}}, and 
#'
#'                  \code{\link{FS.all.reducts.computation}} based on RST.
#' @seealso \code{\link{predict.RuleSetFRST}}, \code{\link{RI.GFRS.FRST}}, and \code{\link{RI.hybridFS.FRST}}.
#' @return A class \code{"RuleSetRST"} containing the following components/attributes:
#' \itemize{
#'	\item \code{rules}: every rule which is a list consisting of five parts as follows:
#'             \itemize{
#'             \item \code{idx}: it is a vector of indexes of attributes involved in each rule.
#'             \item \code{values}: it is a vector of values of each attributes used as the antecedent part.
#'             \item \code{consequent}: it is a vector of a value of decision attribute as the consequent part.
#'             \item \code{support}: it is a numeric value representing the support value of each rule.
#'             \item \code{laplace}: it is a numeric value representing the laplace value of the consequent part.
#'             }
#'  \item \code{uniqueCls}: it is a vector containing the decision classes.
#'  \item \code{clsProbs}: it is vector showing the probability of the decision classes.
#'  \item \code{majorityCls}: it is a value representing the majority value in the rules.
#'  \item \code{method}: it is the type of used method. In this case, it is \code{"indiscernibilityBasedRules"}.
#' }
# @references
# J. Stefanowski and S. Wilk, "Rough Sets for Handling Imbalanced Data: Combining Filtering and Rule-based Classifiers",
# Fundamenta Informaticae, vol. 72, no. 1 - 3, p. 379 - 391 (2006).
# 
#' @examples
#' ###########################################################
#' ## Example 
#' ##############################################################
#' data(RoughSetData)
#' decision.table <- RoughSetData$hiring.dt 							 
#'
#' ## determine feature subset/reduct 	
#' reduct <- FS.permutation.heuristic.reduct.RST(decision.table,  permutation = NULL)
#'						 
#' rules <- RI.indiscernibilityBasedRules.RST(decision.table, reduct)
#' @export
RI.indiscernibilityBasedRules.RST <- function(decision.table, feature.set) {
  
	if (!inherits(decision.table, "DecisionTable")) {
		stop("Provided data should inherit from the \'DecisionTable\' class.")
	}
  
	if (!inherits(feature.set, "FeatureSubset")) {
		stop("Provided feature subset should inherit from the \'FeatureSubset\' class.")
	}
  
	if(is.null(attr(decision.table, "decision.attr"))) {
    stop("A decision attribute is not indicated.")
	} 
	else decisionIdx <- attr(decision.table, "decision.attr")
  
	if(!all(attr(decision.table, "nominal.attrs"))) {
	  stop("Some of the attributes are numerical. 
         Discretize attributes before calling RST-based rule induction methods.")
	}
  
  reduct <- feature.set$reduct
  
	if (length(reduct) >= 1) {
		INDrelation <- BC.IND.relation.RST(decision.table, reduct)$IND.relation
	} 
	else stop("empty feature subset")
  
	clsFreqs <- table(decision.table[[decisionIdx]])
	uniqueCls <- names(clsFreqs)
	ruleSet <- lapply(INDrelation,
                   function(idxs, rule, dataS, clsVec, uniqueCls) laplaceEstimate(list(idx = reduct, values = as.vector(dataS[idxs[1],reduct])), 
                                                                                  dataS, clsVec, uniqueCls, suppIdx = idxs), 
                   reduct, decision.table, decision.table[[decisionIdx]], uniqueCls)
  
	attr(ruleSet, "uniqueCls") <- as.character(uniqueCls)
	attr(ruleSet, "clsProbs") <- clsFreqs/sum(clsFreqs)
	attr(ruleSet, "majorityCls") <- as.character(uniqueCls[which.max(table(decision.table[[decisionIdx]]))])
	attr(ruleSet, "method") <- "indiscernibilityBasedRules"
	attr(ruleSet, "dec.attr") <- colnames(decision.table[decisionIdx])
	attr(ruleSet, "colnames") <- colnames(decision.table)[-decisionIdx]
	
	ruleSet = ObjectFactory(ruleSet, classname = "RuleSetRST")	
	return(ruleSet)  
}

#' It is
#' 
#' @title Rule induction using the CN2 algorithm
#' 
#' @param decision.table a \code{"DecisionTable"} class representing the decision table. See \code{\link{SF.asDecisionTable}}. 
#' @param K 
#' @seealso \code{\link{predict.RuleSetFRST}}.
#' @return A class \code{"RuleSetRST"}. See \code{\link{RI.indiscernibilityBasedRules.RST}}.
#' @references
#' P.E. Clark and T. Niblett, "The CN2 Induction algorithm",
#' Machine Learning, 3, p. 261 - 284 (1986).
#' 
#' @examples
#' ###########################################################
#' ## Example 
#' ##############################################################
#' data(RoughSetData)
#' decision.table <- RoughSetData$hiring.dt 
#' 
#' rules <- RI.CN2Rules.RST(decision.table, K = 3)
#' 
#' @export
RI.CN2Rules.RST <- function(decision.table, K = 3)  {
  
  if (!inherits(decision.table, "DecisionTable")) {
    stop("Provided data should inherit from the \'DecisionTable\' class.")
  }
  
  if(is.null(attr(decision.table, "decision.attr"))) {
    stop("A decision attribute is not indicated.")
	} 
	else decIdx = attr(decision.table, "decision.attr")
  
  if(!all(attr(decision.table, "nominal.attrs"))) {
    stop("Some of the attributes are numerical.
         Discretize attributes before calling RST-based rule induction methods.")
  }
  
  clsVec <- decision.table[,decIdx]
  uniqueCls <- unique(clsVec)
  decisionName = colnames(decision.table)[decIdx]
  decision.table <- decision.table[, -decIdx]
  clsFreqs <- table(clsVec)
  uncoveredObjIdx <- 1:nrow(decision.table)
  
  rules = list()
  
  while (length(uncoveredObjIdx) > 0) {
    rules[[length(rules)+1]] = findBestCN2Rule(decision.table[uncoveredObjIdx,], clsVec[uncoveredObjIdx], uniqueCls, K)
    uncoveredObjIdx = uncoveredObjIdx[setdiff(1:length(uncoveredObjIdx), rules[[length(rules)]]$support)]
  }
  
  attr(rules, "uniqueCls") <- as.character(sort(uniqueCls))
  attr(rules, "clsProbs") <- clsFreqs/sum(clsFreqs)
  attr(rules, "majorityCls") <- as.character(sort(uniqueCls)[which.max(clsFreqs)])
  attr(rules, "method") <- "CN2Rules"
  attr(rules, "dec.attr") <- decisionName
  attr(rules, "colnames") <- colnames(decision.table)[-decIdx]
  
  rules = ObjectFactory(rules, classname = "RuleSetRST")
  
  return(rules);
}

#' It is
#' 
#' @title Rule induction using the LEM2 algorithm
#' @param decision.table a \code{"DecisionTable"} class representing the decision table. See \code{\link{SF.asDecisionTable}}. 
#' @seealso \code{\link{predict.RuleSetFRST}}.
#' @return A class \code{"RuleSetRST"}. See \code{\link{RI.indiscernibilityBasedRules.RST}}.
#' @references
#' J.W. Grzymala-Busse, "A New Version of the Rule Induction System LERS",
#' Fundamenta Informaticae, 31, p. 27 - 39 (1997).
#' 
#' @examples
#' ###########################################################
#' ## Example 
#' ##############################################################
#' data(RoughSetData)
#' decision.table <- RoughSetData$hiring.dt 
#' 
#' rules <- RI.LEM2Rules.RST(decision.table)
#' 
#' @export
RI.LEM2Rules.RST <- function(decision.table)  {
  
  if (!inherits(decision.table, "DecisionTable")) {
    stop("Provided data should inherit from the \'DecisionTable\' class.")
  }
  
  if(is.null(attr(decision.table, "decision.attr"))) {
    stop("A decision attribute is not indicated.")
  } else decIdx = attr(decision.table, "decision.attr")
  
  if(!all(attr(decision.table, "nominal.attrs"))) {
    stop("Some of the attributes are numerical.
         Discretize attributes before calling RST-based rule induction methods.")
  }
  
  clsVec <- decision.table[,decIdx]
  uniqueCls <- unique(clsVec)
  decisionName = colnames(decision.table)[decIdx]
  clsFreqs <- table(clsVec)
  
  INDrelation = BC.IND.relation.RST(decision.table, (1:ncol(decision.table))[-decIdx])
  approximations = BC.LU.approximation.RST(decision.table, INDrelation)
  lowerApproximations = approximations$lower.approximation
  rm(INDrelation, approximations)
  
  descriptorsList = attr(decision.table, "desc.attrs")[-decIdx]
  descriptorCandidates = list()
	for (i in 1:length(descriptorsList)) {
    descriptorCandidates = c(descriptorCandidates,
                             lapply(descriptorsList[[i]],
                                    function(v, x) return(list(idx = x, values = v)), i))
  }
  
  attributeValuePairs = lapply(descriptorCandidates, laplaceEstimate, 
                               decision.table[,-decIdx], clsVec, uniqueCls)
  rm(descriptorsList, descriptorCandidates)
  
  rules = list()
  for(i in 1:length(lowerApproximations)) {
    rules[[i]] = computeLEM2covering(as.integer(lowerApproximations[[i]]), attributeValuePairs, 
                                     clsVec, uniqueCls)
  }
  
  rules = unlist(rules, recursive = FALSE)
  rules = lapply(rules, function(x) laplaceEstimate(list(idx = x$idx, values = x$values), 
                                                    decision.table, clsVec, uniqueCls, suppIdx = x$support))
  
  attr(rules, "uniqueCls") <- as.character(sort(uniqueCls))
  attr(rules, "clsProbs") <- clsFreqs/sum(clsFreqs)
  attr(rules, "majorityCls") <- as.character(sort(uniqueCls)[which.max(clsFreqs)])
  attr(rules, "method") <- "LEM2Rules"
  attr(rules, "dec.attr") <- decisionName
  attr(rules, "colnames") <- colnames(decision.table)[-decIdx]
  
  rules = ObjectFactory(rules, classname = "RuleSetRST")
  
  return(rules);
}


#' It is
#' 
#' @title Rule induction using the AQ-style algorithm
#' @param decision.table a \code{"DecisionTable"} class representing the decision table. See \code{\link{SF.asDecisionTable}}. 
#' @param epsilon a numeric value representing small margin. 
#' @param timesCovered an integer value representing the number of coverage. 
#' @seealso \code{\link{predict.RuleSetFRST}}.
#' @return A class \code{"RuleSetRST"}. See \code{\link{RI.indiscernibilityBasedRules.RST}}.
# @references
# J.W. Grzymala-Busse, "A New Version of the Rule Induction System LERS",
# Fundamenta Informaticae, 31, p. 27 - 39 (1997).
# 
#' @examples
#' ###########################################################
#' ## Example 
#' ##############################################################
#' data(RoughSetData)
#' decision.table <- RoughSetData$hiring.dt 
#'
#' rules <- RI.AQRules.RST(decision.table)
#' 
#' @export
RI.AQRules.RST <- function(decision.table, epsilon = 0.0, timesCovered = 1)  {
	
  if (!inherits(decision.table, "DecisionTable")) {
		stop("Provided data should inherit from the \'DecisionTable\' class.")
	}
  
	if (is.null(attr(decision.table, "decision.attr"))) {
		stop("A decision attribute is not indicated.")
	} else decIdx = attr(decision.table, "decision.attr")
  
	if (!all(attr(decision.table, "nominal.attrs"))) {
		stop("Some of the attributes are numerical.
         Discretize attributes before calling RST-based rule induction methods.")
	}
  
	clsVec <- decision.table[,decIdx]
	uniqueCls <- unique(clsVec)
	decisionName = colnames(decision.table)[decIdx]
	clsFreqs <- table(clsVec)
  
	INDrelation = BC.IND.relation.RST(decision.table, (1:ncol(decision.table))[-decIdx])
	approximations = BC.LU.approximation.RST(decision.table, INDrelation)
	lowerApproximations = approximations$lower.approximation
	rm(INDrelation, approximations)
  
	descriptorsList = attr(decision.table, "desc.attrs")[-decIdx]
	descriptorCandidates = list()
	for (i in 1:length(descriptorsList)) {
		tmpDescriptors = lapply(descriptorsList[[i]],
                            function(v, x) return(list(idx = x, values = v)), i)
		tmpDescriptors = lapply(tmpDescriptors, laplaceEstimate, 
		                        decision.table[,-decIdx], clsVec, uniqueCls)
		names(tmpDescriptors) = descriptorsList[[i]]
		descriptorCandidates[[length(descriptorCandidates) + 1]] = tmpDescriptors
	}
	rm(descriptorsList, tmpDescriptors)
  
	rules = list()
	set.seed(123)
	for(i in 1:length(lowerApproximations)) {
		rules[[i]] = computeAQcovering(as.integer(lowerApproximations[[i]]),
                                   descriptorCandidates, 
		                               decision.table[,-decIdx], 
                                   epsilon = epsilon, K = timesCovered)
	}
  
	rules = unlist(rules, recursive = FALSE)
	rules = lapply(rules, function(x) laplaceEstimate(list(idx = x$idx, values = x$values), 
	                                                  decision.table, clsVec, uniqueCls, suppIdx = x$support))
  
	attr(rules, "uniqueCls") <- as.character(sort(uniqueCls))
	attr(rules, "clsProbs") <- clsFreqs/sum(clsFreqs)
	attr(rules, "majorityCls") <- as.character(sort(uniqueCls)[which.max(clsFreqs)])
	attr(rules, "method") <- "AQRules"
	attr(rules, "dec.attr") <- decisionName
	attr(rules, "colnames") <- colnames(decision.table)[-decIdx]
  
	rules = ObjectFactory(rules, classname = "RuleSetRST")
  
	return(rules) 
}


#' It is a function used to obtain predicted values after obtaining rules by using rule induction methods. 
#' We have provided the functions \code{\link{RI.GFRS.FRST}} and \code{\link{RI.hybridFS.FRST}} to generate rules based on FRST.
#' 
#' @title The predicting function for rule induction methods based on FRST
#'
#' @param object a \code{"RuleSetFRST"} class resulted by \code{\link{RI.GFRS.FRST}} and \code{\link{RI.hybridFS.FRST}}.
#' @param newdata a \code{"DecisionTable"} class containing a data frame or matrix (m x n) of data for the prediction process, where m is the number of instances and 
#' n is the number of input attributes. It should be noted that this data must have \code{colnames} on each attributes.
#' @param ... the other parameters.
#' @seealso \code{\link{RI.indiscernibilityBasedRules.RST}}, \code{\link{RI.GFRS.FRST}} and \code{\link{RI.hybridFS.FRST}}
#' @return The predicted values. 
#' @examples
#' ##############################################
#' ## Example: Classification Task
#' ##############################################
#' data(RoughSetData)
#' decision.table <- RoughSetData$pima7.dt 							 
#'		
#' ## using RI.hybrid.FRST for generating rules
#' control <- list(type.aggregation = c("t.tnorm", "lukasiewicz"), 
#'            type.relation = c("tolerance", "eq.1"), t.implicator = "lukasiewicz")
#' rules.hybrid <- RI.hybridFS.FRST(decision.table, control)
#'
#' ## in this case, we are using the same data set as the training data
#' res.1 <- predict(rules.hybrid, decision.table[, -ncol(decision.table)])
#'
#' ## using RI.GFRS.FRST for generating rules
#' control <- list(alpha.precision = 0.01, type.aggregation = c("t.tnorm", "lukasiewicz"), 
#'                 type.relation = c("tolerance", "eq.3"), t.implicator = "lukasiewicz")						 
#' rules.gfrs <- RI.GFRS.FRST(decision.table, control)
#'
#' ## in this case, we are using the same data set as the training data
#' res.2 <- predict(rules.gfrs, decision.table[, -ncol(decision.table)])
#'
#' ##############################################
#' ## Example: Regression Task
#' ##############################################
#' data(RoughSetData)
#' decision.table <- RoughSetData$housing7.dt 							 
#'							 
#' ## using RI.hybrid.FRST for generating rules
#' control <- list(type.aggregation = c("t.tnorm", "lukasiewicz"), 
#'            type.relation = c("tolerance", "eq.1"), t.implicator = "lukasiewicz")
#' rules <- RI.hybridFS.FRST(decision.table, control)
#'
#' ## in this case, we are using the same data set as the training data
#' res.1 <- predict(rules, decision.table[, -ncol(decision.table)])
#'
#' @aliases predict.FRST
#' @export  
#' @method predict RuleSetFRST
predict.RuleSetFRST <- function(object, newdata, ...) {
	if(!inherits(object, "RuleSetFRST")) stop("not a legitimate rules based on FRST model")
	
	if (any(object$type.method == c("RI.hybridFS.FRST", "RI.GFRS.FRST")) && object$type.model == "FRST"){
		## get data
		rules <- object$rules
		t.similarity <- object$t.similarity
		t.tnorm <- object$t.tnorm
		variance.data <- object$variance.data
		range.data <- object$range.data
		antecedent.attr <- object$antecedent.attr
		
		## filter newdata according to the antecedent part of rules
		testData <- newdata[, c(colnames(newdata) %in% antecedent.attr), drop = FALSE]
		
		## make results not to convert into numeric
		res <- data.frame(stringsAsFactors = FALSE)
		
		for (i in 1 : nrow(newdata)){	
			miu.Ra <- data.frame(stringsAsFactors = FALSE)
			
			for (j in 1 : length(rules$rules)){
				
				## newdata considering attributes on rules
				namesAtt <- colnames(rules$rules[[j]])
				test.dt <- testData[i, match(c(namesAtt[-length(namesAtt)]), names(testData)), drop = FALSE]
				valRules <- rules$rules[[j]]				
				tra.dt <- valRules[-ncol(valRules)]
				
				## build 2 rows matrix for comparison among them
				mat.calc <- rbind(test.dt, tra.dt)	
				
				## get attributes based on attributes in the rules
				P <- which(object$antecedent.attr %in% c(namesAtt[-length(namesAtt)]))
				temp.relation <- NA
				for (ii in 1 : length(P)){
					 if (object$nominal.att[P[ii]] == TRUE){
						 ## overwrite type of similarity when crisp
						 t.similarity <- "boolean"
					 }	

					 temp <- unlist(similarity.equation(test.dt[1, ii], tra.dt[1, ii], t.similarity, delta = 0.2, range.data = range.data[1, P[ii]], 
					 								variance.data = variance.data[1, P[ii]]))
					 if (is.na(temp.relation)){
						 temp.relation <- temp
					 }
					 else {
						 temp.relation <- func.tnorm(temp.relation, temp, t.tnorm = t.tnorm)
					 }
				}
				
				miu.Ra[j, 1] <- valRules[1, ncol(valRules)]
				
				## get indiscernibility relation or membership function
				miu.Ra[j, 2] <- temp.relation	
			}
			
			if (object$type.task == "classification"){
				options(stringsAsFactors = FALSE)
				
				## get fired rules
				if (object$type.method == c("RI.hybridFS.FRST")){
					fired.rules <- which.max(miu.Ra[, 2])
				}
				else if (object$type.method == "RI.GFRS.FRST"){
					miu.Ra[, 2] <- 1 - miu.Ra[, 2]				
					fired.rules <- which(miu.Ra[, 2] <= rules$threshold)					
					if (length(fired.rules) > 1){
						fired.rules <- which.min(miu.Ra[fired.rules, 2])
					}
				}
				
				## get results
				res <- rbind(res, miu.Ra[fired.rules, 1, drop = FALSE])
			}
			else {
				if (object$type.method == c("RI.hybridFS.FRST")){
					miu.Ra <- as.matrix(miu.Ra)
					result <- (miu.Ra[, 2] %*% miu.Ra[, 1, drop = FALSE]) / sum(miu.Ra[, 2])
					res <- rbind(res, result)
				}
			}
		}
	}
	
	rownames(res) <- NULL
	colnames(res) <- object$consequent.attr
	return(res)
}


#' It is a function used to obtain predicted values after obtaining rules by using rule induction methods. 
#' We have provided the functions \code{\link{RI.indiscernibilityBasedRules.RST}} to generate rules based on RST.
#' 
#' @title The predicting function for rule induction methods based on RST
#'
#' @param object a \code{"RuleSetRST"} class resulted by \code{\link{RI.indiscernibilityBasedRules.RST}}.
#' @param newdata a \code{"DecisionTable"} class containing a data frame or matrix (m x n) of data for the prediction process, where m is the number of instances and 
#' n is the number of input attributes. It should be noted that this data must have \code{colnames} on each attribute.
#' @param ... other parameters of a voting strategy (can be applied only to methods which classify objects by voting).
#' @seealso \code{\link{RI.GFRS.FRST}} and \code{\link{RI.hybridFS.FRST}}.
#' @return The predicted values. 
#' @examples
#' ##############################################
#' ## Example: Classification Task
#' ##############################################
#' data(RoughSetData)
#' decision.table <- RoughSetData$hiring.dt 							 
#'		
#' #### Rule induction based on RST ####
#' ## determine feature subset/reduct 	
#' reduct <- FS.permutation.heuristic.reduct.RST(decision.table,  permutation = NULL)
#'
#' ## generate rules						 
#' rules.rst <- RI.indiscernibilityBasedRules.RST(decision.table, reduct)
#'					
#' ## predicting newdata
#' ## in this case, we are using the same data set as the training data
#' res.1 <- predict(rules.rst, decision.table[, -ncol(decision.table)])
#'
#' @aliases predict.RST
#' @export  
#' @method predict RuleSetRST
predict.RuleSetRST <- function(object, newdata, ...) {
  
  if(!inherits(object, "RuleSetRST")) stop("Not a legitimate rules based on RST model.")
  
  if(!inherits(newdata, "DecisionTable")) stop("Provided data should inherit from the \'DecisionTable\' class.")
  
  method = attr(object, "method")
  if (!(method %in% c("indiscernibilityBasedRules", "CN2Rules", "LEM2Rules", "AQRules"))) {
    stop("Unrecognized classification method")
  }
  
  majorityCls = attr(object, "majorityCls")
  uniqueCls = attr(object, "uniqueCls")
  clsProbs = attr(object, "clsProbs")
  
  predVec = switch(method, 
                   indiscernibilityBasedRules = {
                     ruleSet = object[order(sapply(object, function(X) X$laplace), decreasing = TRUE)]
                     INDclasses = sapply(ruleSet, function(x) paste(unlist(x$values), collapse = " "))
                     consequents = sapply(ruleSet, function(x) x$consequent)
                     
                     newdata = do.call(paste, newdata[,ruleSet[[1]]$idx, drop = FALSE])
						sapply(newdata, bestFirst, INDclasses, consequents, majorityCls, uniqueCls, clsProbs) },
                   
                   CN2Rules = apply(as.matrix(newdata), 1, firstWin, ruleSet = object, ...), 
                   
                   LEM2Rules = apply(as.matrix(newdata), 1, rulesVoting, ruleSet = object, ...),
                   
                   AQRules = apply(as.matrix(newdata), 1, rulesVoting, ruleSet = object, ...) )
  
  predVec <- as.data.frame(predVec)
  rownames(predVec) <- NULL
  colnames(predVec) <- "class"
  
  return(predVec)
}

