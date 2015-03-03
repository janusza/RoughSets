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
#' This is a function implementing a fundamental part of RST: the indiscernibility relation.
#' The indiscernibility relation is a binary relation showing whether two objects can be discerned. The detailed description based on theoritical point of view
#' can be seen in \code{\link{A.Introduction-RoughSets}}.
#'  
#' This function is used as a basic function and is needed by other functions such as \code{\link{BC.LU.approximation.RST}}, \code{\link{BC.positive.reg.RST}} for calculating
#' lower and upper approximations and determining the positive region. The formula of the indiscernibility relation has been explained in \code{\link{A.Introduction-RoughSets}}.
#'
#' @title Indiscernibility relation based on rough set theory
#'
#' @param decision.table a \code{"DecisionTable"} class representing a decision table. See \code{\link{SF.asDecisionTable}}. 
#' @param attribute a numerical vector expressing indexes of subsets of attributes to be considered. The default value is \code{NULL} which means that 
#'                 all condition attributes will be considered. It should be noted that in this case, all attributes considered should be nominal attributes,  
#'                 otherwise discretization must be performed first. 
#' @seealso \code{\link{BC.LU.approximation.RST}}, \code{\link{BC.LU.approximation.RST}}
#' @return A class \code{"IndiscernibilityRelation"} which contains
#'          \itemize{
#'          \item \code{IND.relation}: a list representing indiscernibility relation over all objects. 
#'          \item \code{type.relation}: it is \code{"equivalence"}. 
#'          \item \code{type.model}: a string showing the type of model which is used. In this case, it is \code{"RST"} which means rough set theory.
#'          }
#' @references
#' Z. Pawlak, "Rough Sets", International Journal of Computer and Information Sciences, 
#' vol. 11, no. 5, p. 341 - 356 (1982).
#'
#' @examples
#' #############################################
#' ## Example 1: Using simple data set
#' ## Objects must be nominal/symbolic values
#' ## Otherwise, we must use discretization first 
#' #############################################
#' ## Construct decision table as data frame
#' dt.ex1 <- data.frame(c(1,0,2,1,1,2,2,0), c(0, 1,0, 1,0,2,1,1), 
#'                         c(2,1,0,0,2,0,1,1), c(2,1,1,2,0,1,1,0), c(0,2,1,2,1,1,2,1))
#' colnames(dt.ex1) <- c("aa", "bb", "cc", "dd", "ee")
#' decision.table <- SF.asDecisionTable(dataset = dt.ex1, decision.attr = 5, 
#'                                      indx.nominal = c(1:5))
#'
#' ## In this case, we only consider the second and third attributes.
#' P <- c(2,3)
#' 
#' ####### Perform indiscernibility relation #######
#' IND <- BC.IND.relation.RST(decision.table, attribute = P)
#' @export
BC.IND.relation.RST <- function(decision.table, attribute = NULL){
	
	## get data
	objects <- decision.table
	desc.attrs <- attr(decision.table, "desc.attrs")
	nominal.att <- attr(decision.table, "nominal.attrs")
	decision.attr <- attr(decision.table, "decision.attr")
	
	## initialize
	if (is.null(attribute)){
		if(length(decision.attr) > 0) attribute <- (1:ncol(objects))[-decision.attr]
    else attribute <- 1:ncol(objects)
	}

	## check for non nominal attribute
	if (!all(nominal.att[c(attribute)])){
		stop("please discretize attributes before computing an equivalence-based indiscernibility relation")
	}
	
	#compute the indiscernibility classes
	if (length(attribute) == 1){
		IND = split(1:nrow(objects), do.call(paste, list(objects[,attribute])))
	} else {
		IND = split(1:nrow(objects), do.call(paste, objects[,attribute]))
	}
	
	## construct class
	mod <- list(IND.relation = IND, type.relation = "equivalence", type.model = "RST")	
	class.mod <- ObjectFactory(mod, classname = "IndiscernibilityRelation")
	
	return(class.mod)
}

#' This is a function implementing a fundamental part of rough set theory: 
#' lower and upper approximations. The lower and upper approximations determine whether the objects can be certainty or possibly classified in a particular class based on the basis of available knowledge.
#' The detailed theoretical description 
#' can be seen in \code{\link{A.Introduction-RoughSets}}.
#' 
#' This function depends on \code{\link{BC.IND.relation.RST}} which calculates the equivalence classes of the indiscernibility relation. So, it is obvious that 
#' before performing this function, users must execute \code{\link{BC.IND.relation.RST}} first. Furthermore, we provide parameter \code{decision.attr} representing
#' a column index of the decision attribute, so actually, users may choose any index to be considered as the decision attribute. 
#'
#' @title The lower and upper approximations based on rough set
#'
#' @param decision.table a \code{"DecisionTable"} class representing the decision table. See \code{\link{SF.asDecisionTable}}. 
#'                 It should be noted that in this case, attributes considered must be nominal,  
#'                 otherwise the discretization task must be performed first. 
#' @param IND an \code{"IndiscernibilityRelation"} class representing the partitions of the indiscernibility relation. 
#' @param ... other parameters
#' @seealso \code{\link{BC.IND.relation.RST}}, \code{\link{BC.LU.approximation.FRST}}
#' @return A class \code{"LowerUpperApproximation"} representing rough set (the lower and upper approximations). It contains the following components:
#'         \itemize{
#'          \item \code{lower.approximation}: a list containing object indexes included in the lower approximation based on decision concepts. 
#'          \item \code{upper.approximation}: a list containing object indexes included in the upper approximation based on decision concepts.
#'          \item \code{type.model}: a string showing type of the used model. In this case, it is \code{"RST"} means rough set theory.
#'          } 
#' @references
#' Z. Pawlak, "Rough Sets", International Journal of Computer and Information Sciences, 
#' vol. 11, no. 5, p. 341 - 356 (1982).
#'
#' @examples
#' #######################################
#' ## Example: Using simple data set
#' #######################################
#' dt.ex1 <- data.frame(c(1,0,2,1,1,2,2,0), c(0, 1,0, 1,0,2,1,1), 
#'                         c(2,1,0,0,2,0,1,1), c(2,1,1,2,0,1,1,0), c(0,2,1,2,1,1,2,1))
#' colnames(dt.ex1) <- c("aa", "bb", "cc", "dd", "ee")
#' decision.table <- SF.asDecisionTable(dataset = dt.ex1, decision.attr = 5, 
#'                                      indx.nominal = c(1:5))
#'
#' ## Define considered attributes
#' P <- c(2,3)
#' 
#' ####### Compute indiscernibility relation #######
#' IND <- BC.IND.relation.RST(decision.table, attribute = P)
#'
#' ####### Compute lower and upper approximation #####
#' decision.attr <- c(5)
#' roughset <- BC.LU.approximation.RST(decision.table, IND)
#' @export
BC.LU.approximation.RST <- function(decision.table, IND, ...){
	
	## get the data
	objects <- decision.table
	desc.attrs <- attr(decision.table, "desc.attrs")
	nominal.att <- attr(decision.table, "nominal.attrs")
	IND <- IND$IND.relation
	
	if (is.null(attr(decision.table, "decision.attr"))){
		decision.attr <- ncol(objects)
	}
	else {
		decision.attr <- attr(decision.table, "decision.attr")
		## check if index of decision attribute is not in last index, then change their position
		if (decision.attr != ncol(objects)){
			objects <- cbind(objects[, -decision.attr, drop = FALSE], objects[, decision.attr, drop = FALSE])
			nominal.att <- c(nominal.att[-decision.attr], nominal.att[decision.attr])
			decision.attr = ncol(objects)
		}
	}	
	num.att <- ncol(objects)

	## get indiscernibility of decision attribute
	if (!all(nominal.att[-decision.attr])){
		stop("please, discretize attributes before computing an equivalence-based indiscernibility relation")
	} else {
		## get unique decisions
		uniqueDecisions <- as.character(unique(objects[[decision.attr]]))
	}
	
	IND.decision.attr <- lapply(IND, function(x, splitVec) split(x, splitVec[x]), as.character(objects[[decision.attr]]))
	
	## initialization
	lower.appr <- list()
	upper.appr <- list()
	for(i in 1:length(uniqueDecisions)) {
		tmpIdx1 = which(sapply(IND.decision.attr, function(x) (uniqueDecisions[i] %in% names(x))))
		upper.appr[[i]] = unlist(IND[tmpIdx1])
		tmpIdx2 = which(sapply(IND.decision.attr[tmpIdx1], function(x) (length(x) == 1)))
		if(length(tmpIdx2) > 0) {
			lower.appr[[i]] = unlist(IND[tmpIdx1][tmpIdx2])
		}
		else lower.appr[[i]] = integer()
		colnames(lower.appr[[i]]) <- NULL
		colnames(upper.appr[[i]]) <- NULL
	}
	rm(tmpIdx1, tmpIdx2)
  	
	## give the names of list
	names(lower.appr) <- uniqueDecisions
	names(upper.appr) <- uniqueDecisions
	
	## build class
	res <- list(lower.approximation = lower.appr, upper.approximation = upper.appr, type.model = "RST")
	class.mod <- ObjectFactory(res, classname = "LowerUpperApproximation")
	
	return(class.mod)
}


#' It is a function used to implement a fundamental concept of rough set theory: positive region and
#' degree of dependency. The explanation about this concept can be seen in \code{\link{A.Introduction-RoughSets}}.
#' 
#' In order to compute the function, we need to calculate the indiscernibility relation by executing \code{\link{BC.IND.relation.RST}} 
#' and the lower and upper approximations by calling \code{\link{BC.LU.approximation.RST}}.
#'
#' @title Regions based on rough set theory
#'
#' @param decision.table a \code{"DecisionTable"} class representing the decision table. See \code{\link{SF.asDecisionTable}}.
#' @param roughset a \code{"LowerUpperApproximation"} class representing rough sets that are produced by \code{\link{BC.LU.approximation.RST}}.
#' @seealso \code{\link{BC.IND.relation.RST}}, \code{\link{BC.LU.approximation.RST}}, \code{\link{BC.LU.approximation.FRST}}, 
#'
#'          and \code{\link{BC.positive.reg.FRST}}
#' @return A class \code{"PositiveRegion"} containing the following components:
#'         \itemize{
#'         \item \code{positive.reg}: a vector containing indexes of objects belonging to the positive region.
#'         \item \code{degree.dependency}: a value of degree of dependency. 
#'          \item \code{type.model}: a string showing type of the used model. In this case, it is \code{"RST"} means rough set theory.       
#'         } 
#' @references
#' Z. Pawlak, "Rough Sets", International Journal of Computer and Information Sciences, 
#' vol. 11, no. 5, p. 341 - 356 (1982).
#'
#' @examples
#' dt.ex1 <- data.frame(c(1,0,2,1,1,2,2,0), c(0, 1,0, 1,0,2,1,1), 
#'                         c(2,1,0,0,2,0,1,1), c(2,1,1,2,0,1,1,0), c(0,2,1,2,1,1,2,1))
#' colnames(dt.ex1) <- c("aa", "bb", "cc", "dd", "ee")
#' decision.table <- SF.asDecisionTable(dataset = dt.ex1, decision.attr = 5, 
#'                                     indx.nominal = c(1:5))
#'
#' ## in this case, we consider second and third attributes only
#' P <- c(2,3)
#' 
#' ####### Perform indiscernibility relation #######
#' IND <- BC.IND.relation.RST(decision.table, attribute = P)
#'
#' ####### Perform lower and upper approximations #####
#' roughset <- BC.LU.approximation.RST(decision.table, IND)
#' 
#' ####### Determine the positive region ######
#' region <- BC.positive.reg.RST(decision.table, roughset) 
#' 
#' @export
BC.positive.reg.RST <- function(decision.table, roughset){
	
	## get all objects from the lower approximations of decision classes
	positive.reg <- unlist(roughset$lower.approximation)
	names(positive.reg) <- NULL
	
	## get degree of dependecy 
	degree.depend <- length(positive.reg)/nrow(decision.table)
	
	res <- list(positive.reg = positive.reg[order(positive.reg)],
	            degree.dependency = degree.depend, type.model = "RST")
	
	## build class
	class.mod <- ObjectFactory(res, classname = "PositiveRegion")

	return(class.mod)
}

#' This is a function that builds the decision-relative discernibility matrix based on rough set theory.
#' 
#' It was proposed by (Skowron and Rauszer, 1992), and is used to find all reducts. The detailed explanation can be seen in \code{\link{A.Introduction-RoughSets}}.
#'
#' @title The decision-relative discernibility matrix based on rough set theory
#'
#' @param decision.table a \code{"DecisionTable"} class representing a decision table. See \code{\link{SF.asDecisionTable}}. 
#' @param range.object a vector representing considered objects to construct the \eqn{k}-relative discernibility matrix. 
#'                The default value is \code{NULL} which means that we are using all objects in the decision table.
#' @param show.discernibilityMatrix a boolean value determining whether the discernibility matrix will be shown or not. The default value is \code{FALSE}. 
#' @return A class \code{"DiscernibilityMatrix"} containing the following components: 
#' \itemize{
#' \item \code{disc.mat}: a matrix showing the decision-relative discernibility matrix \eqn{M(\mathcal{A})} 
#'        which contains \eqn{n \times n} where \eqn{n} is the number of objects.
#' \item \code{disc.list}: it refers to the decision-relation discernibility matrix in a list format.
#' \item \code{discernibility.type}: it is \code{"RST"}.
#' \item \code{type.model}: in this case, it is \code{"RST"}.
#' }
#' @examples
#' #######################################################################
#' ## Example 1: Constructing the decision-relative discernibility matrix
#' #######################################################################
#' dt.ex1 <- data.frame(c(1,0,2,1,1,2,2,0), c(0, 1,0, 1,0,2,1,1), 
#'                         c(2,1,0,0,2,0,1,1), c(2,1,1,2,0,1,1,0), c(0,2,1,2,1,1,2,1))
#' colnames(dt.ex1) <- c("aa", "bb", "cc", "dd", "ee")
#' decision.table <- SF.asDecisionTable(dataset = dt.ex1, decision.attr = 5, 
#'                                      indx.nominal = c(1:5))
#'
#' ## build the decision-relation discernibility matrix
#' res.1 <- BC.discernibility.mat.RST(decision.table, range.object = NULL)
#'
#' @seealso \code{\link{BC.IND.relation.RST}}, \code{\link{BC.LU.approximation.RST}}, and \code{\link{BC.LU.approximation.FRST}}
#' @references
#' A. Skowron and C. Rauszer,  
#' "The Discernibility Matrices and Functions in Information Systems", 
#' in: R. Słowinski (Ed.), Intelligent Decision Support: Handbook of Applications and
#' Advances of Rough Sets Theory, Kluwer Academic Publishers, Dordrecht, Netherland,  
#' p. 331 - 362 (1992).
#' @export
BC.discernibility.mat.RST <- function(decision.table, range.object = NULL, show.discernibilityMatrix = FALSE){

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
	num.object <- nrow(objects)
	names.attr <- t(colnames(objects)[-decision.attr])
	
	## replace if range.object = NULL
	if (is.null(range.object)){
		range.object <- matrix(c(1, nrow(objects)), nrow = 1)
		num.object <- range.object[1, 2] - range.object[1, 1] + 1
	}
	
	## initialize the discernibility matrix
  disc.mat <- array(list(NA), dim = c(num.object, num.object, 1))
	
	decVector = as.character(objects[, ncol(objects)])
	dataMatrix = as.matrix(objects[, -ncol(objects)])
  
   ## construct the discernibility matrix
	for (i in 1 : (num.object-1)){
		tmpIdx = which(decVector[(i+1) : num.object] != decVector[i])
		if(length(tmpIdx) > 0)  {
			tmpIdx = tmpIdx + i
			for (j in tmpIdx){
    			## select different values on decision attribute only
  				## construct one element has multi value/list which object is discerned.
  				disc.attr <- names.attr[which(dataMatrix[i,] != dataMatrix[j,])]
  				disc.mat[j, i, 1] <- list(disc.attr)
			}
		}
		rm(tmpIdx)
	}
	disc.mat = as.data.frame(disc.mat)
	disc.list = unique(do.call(c, disc.mat))[-1]
	
	## build class
	if (show.discernibilityMatrix){
		discernibilityMatrix = list(disc.mat = disc.mat, disc.list = disc.list, 
                              names.attr = colnames(decision.table), type.discernibility = "RST", type.model = "RST")
	}
	else {
		discernibilityMatrix = list(disc.list = disc.list, 
                              names.attr = colnames(decision.table), type.discernibility = "RST", type.model = "RST")
	}
	discernibilityMatrix = ObjectFactory(discernibilityMatrix, classname = "DiscernibilityMatrix")
	return(discernibilityMatrix)
}

