#include "R.h"
#include "Rmath.h"

#include <vector>

extern "C" {
void  chooseBestCutC(int *k, double *cutCandidates, int *N, double *vec, int *objectsIdx, int *objectsIdxLengths,
                     int *numOfInt, int *decVec, int *nOfDec, int *attrType, int *minIntervalSize,
                     int *rmVec, int *idxVec, double *maxTPtoFP)	{

  int maxScore = 0;
  int maxIndIdx = 0;
  int curScore = 0;
  int flag = false;
  int sFlag = false;


  if(k[0] >= 1)  {
  int posAttr = 0;
  int negAttr = 0;

  std::vector<int> posAttrDecisionCounts(nOfDec[0], 0);
  std::vector<int> negAttrDecisionCounts(nOfDec[0], 0);

  int tmpObjIdx = 0;
  double tmpIntLength = 0.0;
  for(int j = 0; j <= k[0] - 1;  ++j) {
    if(attrType[0] > 0) {
      tmpObjIdx = 0;
      flag = false;
      curScore = 0;
      for(int m = 0; m <= numOfInt[0] - 1; ++m)  {
        posAttr = 0;
        negAttr = 0;
        for(int i = 0; i <= nOfDec[0] - 1; ++i) {
          posAttrDecisionCounts[i] = 0;
          negAttrDecisionCounts[i] = 0;
        }
        for(int i = tmpObjIdx; i <= tmpObjIdx + objectsIdxLengths[m] - 1; ++i) {
          if(vec[objectsIdx[i] - 1] >= cutCandidates[j]) {
            posAttr++;
            posAttrDecisionCounts[decVec[i] - 1]++;
          }
          else {
            negAttr++;
            negAttrDecisionCounts[decVec[i] - 1]++;
          }
        }
        if( posAttr >= 0.5*objectsIdxLengths[m] ) tmpIntLength = posAttr - 0.5*objectsIdxLengths[m];
        else tmpIntLength = 0.5*objectsIdxLengths[m] - posAttr;
        if( 0.5*objectsIdxLengths[m] - tmpIntLength >= minIntervalSize[0])  {
          if(flag == false) {
            curScore = 0;
            flag = true;
            sFlag = true;
          }
          for(int i = 0; i <= nOfDec[0] - 1; ++i) {
             curScore = curScore + (negAttr - negAttrDecisionCounts[i])*posAttrDecisionCounts[i] + (posAttr - posAttrDecisionCounts[i])*negAttrDecisionCounts[i];
          }
        }
        tmpObjIdx = tmpObjIdx + objectsIdxLengths[m];
      }
      if(curScore > maxScore) {
        maxScore = curScore;
        maxIndIdx = j + 1;
      }
    }
    if(flag == false) {
      rmVec[j] = 1;
    }
  }
  }

  if(sFlag == false) {
    idxVec[0] = 0;
    maxTPtoFP[0] = 0;
  }
  else  {
    idxVec[0] = maxIndIdx;
    maxTPtoFP[0] = maxScore;
  }

  return;
}
}

extern "C" {
void  chooseCutCandidatesC(double *vec, int *decVec, int *N, int *candidatesIdx, double *candidates) {

  for(int i = 0; i < N[0] - 1; ++i)	{
    if(decVec[i] != decVec[i+1])	{
      candidates[i] = (vec[i] + vec[i+1])*0.5;
      candidatesIdx[i] = true;
    }
  }
  return;
}
}

extern "C" {

void  computeIndiscernibilityAndChaos(int *INDclasses, int *INDsizes, int *NOfINDClasses,
                                      int *attrValues, int *NOfAttrValues,
                                      int *decValues, int *NOfDecs,
                                      int *output)  {

  int objectIterator = 0;
  for (int i=0; i<NOfINDClasses[0]; ++i) {
    for (int j=0; j<INDsizes[i]; ++j) {
      ++output[i*NOfAttrValues[0]*NOfDecs[0] +
               (attrValues[INDclasses[objectIterator] - 1] - 1)*NOfDecs[0] +
               decValues[INDclasses[objectIterator] - 1] - 1];
      ++objectIterator;
    }
  }

  return;
}
}

