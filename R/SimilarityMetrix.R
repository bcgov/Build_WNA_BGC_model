library(data.table)
library(sf)
library(terra)
library(catsim)

mod1 <- vect("Y:\\WNA_BGC\\USA_BGC_comparisons/WA_BGC_4Jan2024_eliminated_33var_ClimNA.gpkg")
mod2 <- vect("Y:\\WNA_BGC\\USA_BGC_comparisons/WA_BGC_9Jan2024_eliminated_33var_Colin.gpkg")

compare_maps <- function(mod1, mod2){
  cwtab <- data.table(BGC = unique(c(mod1$BGC, mod2$BGC)))
  cwtab[, bgc_id := 1:nrow(cwtab)]
  cwtab[, Subzone := gsub("_.*","",BGC)]
  cwtab[, Zone := gsub("[a-z]*_[A-Z]*","",BGC)]
  cwtab[,zone_id := as.numeric(as.factor(Zone))]
  cwtab[,subzone_id := as.numeric(as.factor(Subzone))]
  
  mod1 <- merge(mod1, cwtab, by = "BGC")
  mod2 <- merge(mod2, cwtab, by = "BGC")
  
  ref_rast <- rast(mod1, resolution = 200)
  ##subzone variants
  rast1 <- rasterize(mod1, ref_rast, field = "bgc_id")
  rast2 <- rasterize(mod2, ref_rast, field = "bgc_id")
  rast1_sz <- subst(rast1, from = cwtab$bgc_id, to = cwtab$subzone_id)
  rast2_sz <- subst(rast2, from = cwtab$bgc_id, to = cwtab$subzone_id)
  rast1_z <- subst(rast1, from = cwtab$bgc_id, to = cwtab$zone_id)
  rast2_z <- subst(rast2, from = cwtab$bgc_id, to = cwtab$zone_id)
  
  rast_out <- copy(ref_rast)
  
  values(rast_out) <- NA
  rast_out[!is.na(rast2)] <- 0
  rast_out[rast1 != rast2] <- 1  
  rast_out[rast1_sz != rast2_sz] <- 2
  rast_out[rast1_z != rast2_z] <- 3
  #plot(rast_out)
  
  dat_freq <- as.data.table(freq(rast_out))
  
  dat_freq[,Total := sum(count)]
  dat_freq[,Prop := count/Total]
  return(list(summary = dat_freq, map = rast_out))
}

test1 <- compare_maps(mod1,mod2)

mod1 <- vect("Y:\\WNA_BGC\\USA_BGC_comparisons/WA_BGC_9Jan2024_eliminated_46var_Colin.gpkg")
mod2 <- vect("Y:\\WNA_BGC\\USA_BGC_comparisons/WA_BGC_9Jan2024_eliminated_33var_Colin.gpkg")

test2 <- compare_maps(mod1, mod2)

mod1 <- vect("Y:\\WNA_BGC\\USA_BGC_comparisons/WA_BGC_9Jan2024_eliminated_allvar.gpkg")
mod2 <- vect("Y:\\WNA_BGC\\USA_BGC_comparisons/WA_BGC_9Jan2024_eliminated_33var_Colin.gpkg")

test3 <- compare_maps(mod1, mod2)
plot(test3$map)
test3$summary

mod1 <- vect("Y:\\WNA_BGC\\USA_BGC_comparisons/WA_BGC_9Jan2024_eliminated_46var_Colin.gpkg")
mod2 <- vect("Y:\\WNA_BGC\\USA_BGC_comparisons/WA_BGC_9Jan2024_eliminated_33var_Colin.gpkg")

test4 <- compare_maps(mod1, mod2)
plot(test4$map)
print(test4$summary)

##balanced vs not
mod1 <- vect("Y:\\WNA_BGC\\USA_BGC_comparisons/WA_BGC_9Jan2024_eliminated_allvar.gpkg")
mod2 <- vect("Y:\\WNA_BGC\\USA_BGC_comparisons/WA_BGC_9Jan2024_eliminated_allvar_balanced.gpkg")

test5 <- compare_maps(mod1, mod2)
plot(test5$map)
print(test4$summary)

# temp1 <- as.matrix(rast1, wide = T)
# temp2 <- as.matrix(rast2, wide = T)
# catmssim_2d(temp1,temp2,levels = 3)

writeRaster(rast_out, "ModelDiff.tif")

rast_z_all <- rast1_z - rast2_z
rast_all <- rast1 - rast2


plot(rast2_z)

r1 <- rast(mod2)
colnames(mod1)[1] <- "BGC1"
colnames(mod2)[1] <- "BGC2"
mod_comb <- st_intersection(mod1, mod2)
mod_comb$diff <- fifelse(mod_comb$BGC1 == mod_comb$BGC2, 0, 
                         fifelse(gsub("_*","",mod_comb$BGC1) == gsub("_*","",mod_comb$BGC2), 1, 
                                 fifelse(gsub("[a-z]*_[A-Z]*","",mod_comb$BGC1) == gsub("[a-z]*_[A-Z]*","",mod_comb$BGC2),2,5)))

st_write(mod_comb, "WA_BGC_ClimNA_vs_Colin.gpkg", overwrite = T)
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
