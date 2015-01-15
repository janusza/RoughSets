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
# It is used to get discernibility matrix for discretization and the function W
# @param data a matrix
# @param cut.val a vector contains 2 values: index and a cut value of an attribute
# @param type a type of computation whether "table" or "func.w"
# @param nrow.data number of row of data
# @param ncol.data number of coloumn of data
def.discern.mat <- function(data, cut.val, type = "table"){
######################################
## Reference:
## Jan G. Bazan, Hung Son Nguyen, Sinh Hoa Nguyen, Piotr Synak, Jakub Wroblewski, 
## "Rough Set Algorithms in Classification Problem", Chapter 2
##  In: L. Polkowski, S. Tsumoto and T.Y. Lin (eds.): Rough Set Methods and Applications
## Physica-Verlag, Heidelberg, New York 2000, pp. 49 - 88.
######################################
	if (type == "table"){	
		temp.test <- matrix(1, nrow = nrow(data), ncol = 1)
		temp.test[which(data[, cut.val[1]] <= cut.val[2])] <- 0
		dec.data <- data[, ncol(data), drop = FALSE]
				
		temp.all.test <- matrix(ncol = 1)	
		
		# it is a function for vectorizing and used to construct discernibility matrix which is pairs of objects toward cut values
		# @param i index 
		# @param j index
		# @param dec.data values of decision attributes
		# @param temp.test
		func.def.discern <- function(i, j, dec.data, discrete.val){
			if ((j > i) && (dec.data[i, 1] != dec.data[j, 1])){
				if (temp.test[i] != temp.test[j]){
					temp.all.test <<- rbind(temp.all.test, 1)
				}
				else {
					temp.all.test <<- rbind(temp.all.test, 0)
				}
			}
			else {
				temp.all.test <<-  rbind(temp.all.test, NA)
			}
		}

		## vectorize func.def.discern
		VecFun <- Vectorize(func.def.discern, vectorize.args=list("i","j"))
		outer(1 : nrow(dec.data), 1 : nrow(dec.data), VecFun, dec.data, temp.test)
		temp.all.test <- na.omit(temp.all.test)
			
		return(temp.all.test)				
	} 
	else if (type == "func.w"){
		## get indx which element values are less then cut value
		indx.left <- which(data[, cut.val[1]] < cut.val[2])
		
		## get indx which element values are greater then cut value
		indx.right <- which(data[, cut.val[1]] > cut.val[2])
		
		## calculate number of object which is less/greater than cut value
		lx.ac <- length(indx.left)
		rx.ac <- length(indx.right)
		
		## get decision class
		dec.class <- unique(c(data[, ncol(data)]))
		
		## initialization
		left.temp <- matrix()
		right.temp <- matrix()
		
		## compute left hand side and right hand side considering their decision class
		for (i in 1 : length(dec.class)){
			left.temp[i] <- length(which(data[indx.left, ncol(data)] == dec.class[i]))
			right.temp[i] <- length(which(data[indx.right, ncol(data)] == dec.class[i]))
		}
		
		## calculate total number
		sum.temp <- sum(left.temp * right.temp)
		
		## function W
		return(lx.ac * rx.ac - sum.temp)
	}
}
		
# a function for discretizing a single attribute
# @param vec a vector of data
# @param cuts a vector of cut values
applyDiscretization <- function(vec, cuts, isNominal) {
   if(!isNominal) vec = cut(vec, c(-Inf,cuts,Inf), right=TRUE, include.lowest = TRUE)
   return(as.character(vec))
}

# a function for computing cuts of the "quantile-based" discretization of a single attribute into a given number of intervals
# @param vec a vector of data
# @param n a number 
discretize.quantiles <- function(vec, n) {
  uniqueValues = unique(vec)
  if (length(uniqueValues) < n) {
		if(length(uniqueValues) == 1) cutVec = NULL
		else cutVec = uniqueValues
	} 
	else {
		cutVec = quantile(vec, (1:(n-1))/n)
		cutVec = unique(cutVec)
	}
  
  return(cutVec)
}

# a function for computing cuts of the "equal interval size" discretization of a single attribute into a given number of intervals
# @param vec a vector of data
# @param n a number
discretize.equal.intervals <- function(vec, n) {
  attrRange = range(vec)
  if(diff(attrRange) > 0 && n > 1) {
		cutVec = seq(attrRange[1], attrRange[2], length.out = n + 1)
  }
  else {
		cutVec = NULL
  }
  
  return(cutVec[-c(1,length(cutVec))])
}

# a function for computing cuts using the maximum discernibility heuristic (the global approach)
global.discernibility <- function(vecList, cutCandidatesVecList, decVec, nOfCuts, 
                                 nAttrs = length(vecList), minIntSupport = 0, ...) {
  attrCount = length(vecList)
  cutVecList = list()
  rmVecList = list()
  notRemovedVecIndicator = sapply(cutCandidatesVecList, function(x) return(length(x) > 0))
  cutVecList[1:attrCount] = list(numeric())
  
  INDclasses = list(1:length(decVec))
  minIntervalSize = ceiling(length(decVec)*minIntSupport)
  scrHistVec = rep(0, attrCount)
  maxDiscernPairs = conflictsCouner(decVec)
  i = 0
  numOfChosenCuts = 0
  nDecisions = length(unique(decVec))
  endFlag = F
  
  while(numOfChosenCuts < nOfCuts && !endFlag) {
    i = i + 1
    rmVecList[1:attrCount] = list(integer())
    candidateVecIdx = which(notRemovedVecIndicator)
    attrSampleIdx = sample(candidateVecIdx, min(nAttrs, length(candidateVecIdx)), replace = F)
    tmpINDclassesVec = unlist(INDclasses)
    tmpObjIdxLengths = sapply(INDclasses, length)
    
    bestCutsList = mapply(evaluateCuts, vecList[attrSampleIdx], cutCandidatesVecList[attrSampleIdx], 
                          MoreArgs = list(decVec = decVec[tmpINDclassesVec],
                                          nOfDec = nDecisions,
                                          INDclassesVec = tmpINDclassesVec, 
                                          INDclassesSizes = tmpObjIdxLengths, 
                                          minIntervalSize = minIntervalSize),
                          SIMPLIFY = F)
    
    maxCutScoreVec = sapply(bestCutsList, function(x) x$maxTPtoFP)
    maxCutIdxVec = sapply(bestCutsList, function(x) x$maxIdx)
    rmVecList[attrSampleIdx] = lapply(bestCutsList, function(x) x$rmVec) 
    rm(tmpINDclassesVec, tmpObjIdxLengths, bestCutsList)
    
    maxScr = max(maxCutScoreVec)
    tmpIdxVec = which(maxCutScoreVec == maxScr)
    bestAttrIdx = tmpIdxVec[which.max(scrHistVec[attrSampleIdx[tmpIdxVec]])]
    numOfChosenCuts = numOfChosenCuts + length(bestAttrIdx)
    scrHistVec[attrSampleIdx] = maxCutScoreVec
    chosenCutIdx = maxCutIdxVec[bestAttrIdx]
    rm(maxCutScoreVec, maxCutIdxVec, tmpIdxVec)
    
    if(maxScr == 0 || numOfChosenCuts >= nOfCuts) {
      endFlag = TRUE
      if(numOfChosenCuts >= nOfCuts) 
        cutVecList[[attrSampleIdx[bestAttrIdx]]] = c(cutVecList[[attrSampleIdx[bestAttrIdx]]],
                                                     cutCandidatesVecList[[attrSampleIdx[bestAttrIdx]]][chosenCutIdx])
    }
    else  { 
      cutVecList[[attrSampleIdx[bestAttrIdx]]] = c(cutVecList[[attrSampleIdx[bestAttrIdx]]],
                                                   cutCandidatesVecList[[attrSampleIdx[bestAttrIdx]]][chosenCutIdx])
      rmVecList[[attrSampleIdx[bestAttrIdx]]] = c(rmVecList[[attrSampleIdx[bestAttrIdx]]], 
                                                  chosenCutIdx)
      
      tmpDiscretizedVec = as.integer(vecList[[attrSampleIdx[bestAttrIdx]]] >= 
                                       cutCandidatesVecList[[attrSampleIdx[bestAttrIdx]]][chosenCutIdx])
      newINDclasses = unlist(lapply(INDclasses, splitINDclass, tmpDiscretizedVec), recursive = F)
      tmpClassesToRmIdx = which(sapply(newINDclasses, function(x, decV) length(unique(decV[x])) == 1,
                                       as.character(decVec)))
      if(length(tmpClassesToRmIdx) > 0) newINDclasses = newINDclasses[-tmpClassesToRmIdx]
      
      INDclasses = newINDclasses
      if(length(INDclasses) > 0)  maxDiscernPairs = sum(sapply(INDclasses, function(x, decV) conflictsCouner(decV[x]), decVec))
      else maxDiscernPairs = 0
      rm(newINDclasses, tmpClassesToRmIdx, tmpDiscretizedVec)
      
      for(j in attrSampleIdx) {
        if(length(rmVecList[[j]]) > 0)  {
          cutCandidatesVecList[[j]] = cutCandidatesVecList[[j]][-rmVecList[[j]]]
          if(length(cutCandidatesVecList[[j]]) == 0) notRemovedVecIndicator[j] = F
        }
      }
      
      if(maxDiscernPairs == 0 || max(sapply(INDclasses,length)) < 2*minIntervalSize)  endFlag = TRUE
    }
  }
  
  return(cutVecList)
}

#  an auxiliary function for computation of a number of conflicts with a decision vector
conflictsCouner <- function(decisionVector)  {
   decisionDistrib = as.numeric(table(decisionVector))
   return(as.numeric(sum(as.numeric(sum(decisionDistrib) - decisionDistrib) * decisionDistrib)))
}

# an auxiliary function for evaluation of a set of candidate cuts using the number of conflicts
# it uses a C++ code to speed up the computations
evaluateCuts <- function(numVec, cutsCandidates, decVec, nOfDec,
                        INDclassesVec, INDclassesSizes, minIntervalSize) {
  
  bestCut = .C("chooseBestCutC", k = as.integer(length(cutsCandidates)), 
               cutCandidates = as.double(cutsCandidates), 
               N = as.integer(length(INDclassesVec)), 
               vec = as.double(numVec),
               objectsIdx = as.integer(INDclassesVec), 
               objectsIdxLengths = as.integer(INDclassesSizes),
               numOfInt = as.integer(length(INDclassesSizes)), 
               decVec = as.integer(decVec),
               nOfDec = as.integer(nOfDec),
               attrType = as.integer(T), 
               minIntervalSize = as.integer(minIntervalSize),
               rmVec = as.integer(rep(0, length(cutsCandidates))), 
               idxVec = as.integer(0), 
               maxTPtoFP = as.double(0.0))
  
  return(list(maxTPtoFP = as.numeric(bestCut$maxTPtoFP),
              maxIdx = as.integer(bestCut$idxVec),
              rmVec = which(bestCut$rmVec > 0)))
}

# an auxiliary function for construction of sets of candidate cuts using for a given attribute
# it uses a C++ code to speed up the computations
chooseCutCandidates <- function(attrVec, decVec)  {
  
  tmpIdx = order(attrVec)
  tmpCutCandidates = .C("chooseCutCandidatesC", vec = as.double(attrVec[tmpIdx]), 
                        decVec = as.integer(decVec[tmpIdx]), 
                        N = as.integer(length(decVec)),
                        candidatesIdx = as.integer(rep(FALSE,length(decVec)-1)), 
                        candidates = as.double(rep(0.0,length(decVec)-1)))
  return(unique(tmpCutCandidates$candidates[which(as.logical(tmpCutCandidates$candidatesIdx))]))
}

# It is used to generate cut values
# @param dataset a data which contains non nominal values
gen.cut.values <- function(dataset){
	## initialization
	cut.values <- matrix(NA, nrow = 1, ncol = 2)
	
	## loop for all continuous attributes
	for (i in 1:ncol(dataset)){
	
		## initialization
		cut.att.i <- matrix()
		
		## sort values of continuous attributes
		non.nominal.obj.sort <- sort(unique(dataset[, i]))
		
		## calculate cut values
		cut.att.i <- non.nominal.obj.sort[-length(non.nominal.obj.sort)] + 0.5 * diff(non.nominal.obj.sort)
		cut.values <- rbind(cut.values, cbind(cut.att.i, i))
	}

	cut.values <- na.omit(cut.values)
	colnames(cut.values) <- NULL
	
	return(cut.values)
}