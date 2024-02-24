library(data.table)
library(ranger)
library(foreach)
library(terra)
library(sf)
library(dplyr)
library(climr)
library(vip)
library(tictoc)

addVars <- function(dat){ ##this function modifies everything inplace, so no need for return value
  dat[,`:=`(PPT_MJ = PPT05+PPT06,
            PPT_JAS = PPT07+PPT08+PPT09,
            PPT.dormant = PPT_at+PPT_wt)]
  dat[,`:=`(CMD.def = 500-PPT.dormant)]
  dat[CMD.def < 0, CMD.def := 0]
  dat[,`:=`(CMDMax = CMD07,
            CMD.total = CMD.def + CMD)]
  dat[,`:=`(CMD.grow = CMD05+CMD06+CMD07+CMD08+CMD09,
            DD5.grow = DD5_05+DD5_06+DD5_07+DD5_08+DD5_09,
            #DDgood = DD5 - DD18,
            #DDnew = (DD5_05+DD5_06+DD5_07+DD5_08)-(DD18_05+DD18_06+DD18_07+DD18_08),
            TmaxJuly = Tmax07)]
}


bgc_info <- fread("../Common_Files/BGC_Info/WNA_BGCs_Info_v12_15.csv")
BC_BGCs <- bgc_info[grep("BGC.*",Source),BGC]



# bgc_map$ID <- seq_along(bgc_map$BGC)
# 
# neighbours_ls <- list()
# for(bgc in BC_BGCs[-(1:180)]){
#   cat(".")
#   focal <- bgc_map[bgc_map$BGC == bgc,]
#   if(nrow(focal) > 0){
#     focal <- st_union(focal$geom)
#     neighbours <- st_intersects(focal, bgc_map)
#     neighbours_ls[[bgc]] <- neighbours[[1]]
#   }
# }
# 
# saveRDS(neighbours_ls,"BC_BGC_NeighbourList.rds",)
# st_write(bgc_map, "BC_BGCs_with_ID.gpkg")
bgc_map <- st_read("BC_BGCs_with_ID.gpkg")
dem <- terra::rast("D:/DEMs/WNA_DEM_4326_clipped.tif")# %>% terra::project("epsg:4326")
#writeRaster(dem, "D:/DEMs/WNA_DEM_4326_clipped.tif")
#dem <- rast("//objectstore2.nrs.bcgov/ffec/DEM/DEM_NorAm/NA_Elevation/data/northamerica/northamerica_elevation_cec_2023.tif")
#dem <- dem.noram

neighbours_ls <- readRDS("BC_BGC_NeighbourList.rds")
remove = c("PPT_at", "PPT_wt", "CMD.def")
vars.selected <- fread("no_month_vars.csv") %>% dplyr::filter(!vars %in% remove)

bgc_list <- names(neighbours_ls)
bgc_list <- bgc_list[!grepl("BAFA.* | *.un.", bgc_list)]
bgc_list <- bgc_list[-1]

tic()

vars <- c(vars.selected$vars,"BGC")
bgc_list_short <- c("BGxh3", "BWBSdk", "CDFmm", "ICHmc1", "CWHmm1", "ESSFwk1","ICHmw1","IDFdk3",
                    "IDFxx2","MSdw","SBSdk","MSxv")
bgc_list <- bgc_list[-"CWHvm2"]
tic()
res_list <- list()
for(bgc in bgc_list_short){
  cat("Processing",bgc,"\n")
  out <- bgc_map[bgc_map$ID %in% neighbours_ls[[bgc]],]
  out_union <- group_by(out, BGC) %>% 
    summarize(geom = st_union(geom),
              BGC = BGC[1])
  pnts <- st_sample(out_union, size = rep(150, nrow(out_union)), type = "random", by_polygon = T)
  pnts_all <- st_as_sf(data.frame(BGC = rep(out_union$BGC, each = 150), geom = pnts))
  pnts_all <- st_transform(pnts_all, 4326)
  coords <- st_coordinates(pnts_all)
  temp_elev <- terra::extract(dem, coords)
  coords <- data.frame(coords, elev = temp_elev$WNA_DEM_3005_clipped,
                       ID = 1:nrow(pnts_all), BGC = pnts_all$BGC)
  
  # coords <- data.frame(coords, elev = temp_elev$WNA_DEM_SRT_30m, 
  #                      ID = 1:nrow(pnts_all), BGC = pnts_all$BGC)
  
  clim_vars <- suppressMessages(climr_downscale(coords, which_normal = "composite_normal", 
                                                vars = c(list_variables(),"CMI"), return_normal = T))  
  clim_vars <- data.table:::na.omit.data.table(clim_vars)
  addVars(clim_vars)
  clim_vars <- clim_vars[,..vars]
  clim_vars[,BGC := as.factor(BGC)]
  rf_mod <- ranger(BGC ~ ., data = clim_vars, num.trees = 101, importance = "impurity", splitrule = "gini")
  varimp <- sort(importance(rf_mod),decreasing = T)[1:10]
  res_list[[bgc]] <- data.table(Focal = bgc, Var = names(varimp), 
                                Importance = unname(varimp), 
                                OOB = rf_mod$prediction.error,
                                NumberBGCs = nrow(out_union))
}


dat_all2 <- rbindlist(res_list)
fwrite(dat_all, "Focal_Variable_Importance_v2.csv")
test <- dat_all[,.(Num = .N), by = .(Var)]
setorder(test, -Num)
toc()
# 
# 
# st_write(out, "Test_IDFdk3_Neighbours.gpkg")
# st_write(pnts_all, "Test_pnts.gpkg")
# 
# bgc_vect <- vect(bgc_map)
# focal_vect <- vect(focal)
# test <- terra::nearby(focal_vect,bgc_vect,distance = 0.1)
# plot(bgc_vect[test[,2]])