library(data.table)

##compare data.tables of two predicitions. Inputs must only have two columns: ID and BGC_pred (names don't matter)
##returns spatial similarity, aspatial similarity, and optionally the confusion matrix
simple_comparison <- function(preds1, preds2, return_confmat = FALSE){ ##columns: ID, BGC prediction. IDs must match
  setnames(preds1, c("ID","BGC1"))
  setnames(preds2, c("ID","BGC2"))
  preds1[preds2, BGC2 := i.BGC2, on = "ID"]
  
  preds1[,Same := BGC1 == BGC2]
  preds_sum <- preds1[,.(Num = .N), by = .(Same)]
  spatial_acc <- preds_sum[Same == T, Num]/sum(preds_sum$Num)
  
  sm1 <- as.data.table(table(preds1$BGC1))
  sm2 <- as.data.table(table(preds1$BGC2))
  sm1[sm2, N_Pred2 := i.N, on = "V1"]
  sm1[,sum_diff := abs(N - N_Pred2)]
  aspatial_acc <- (sum(sm1$N) - sum(sm1$sum_diff))/sum(sm1$N)
  
  if(return_confmat){
    confMat <- dcast(preds1, BGC1 ~ BGC2, fun.aggregate = length)
    out <- list(spatial_acc = spatial_acc, aspatial_acc = aspatial_acc, confMat = confMat)
  } else {
    out <- list(spatial_acc = spatial_acc, aspatial_acc = aspatial_acc)
  }
  return(out)
}
