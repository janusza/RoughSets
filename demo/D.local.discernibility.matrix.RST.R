 #################################################################
 ## Example: Determine cut values and generate new decision table
 #################################################################
 library(RoughSets)
 
 dt.ex1 <- data.frame(c(1, 2, 3, 3, 4, 5, 6, 7, 7, 8), c(2,5, 7, 6, 6, 6, 1, 8, 1, 1),
                              c(3, 5, 1, 1, 3, 5, 8, 8, 1, 1), c(0, 1, 2, 1, 0, 1, 2, 2, 0, 0))
 colnames(dt.ex1) <- c("a1", "a2", "a3", "d")
 decision.table <- SF.asDecisionTable(dataset = dt.ex1, decision.attr = 4, indx.nominal = c(4)) 

 cut.values <- D.local.discernibility.matrix.RST(decision.table)
  
 ## generate new decision table
 new.decTable <- SF.applyDecTable(decision.table, cut.values)