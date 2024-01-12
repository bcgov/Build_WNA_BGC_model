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
dem <- terra::rast("E:/WNA_DEM_SRT_30m_cropped.tif")
dem.noram <- rast("//objectstore2.nrs.bcgov/ffec/DEM/DEM_NorAm/NA_Elevation/data/northamerica/northamerica_elevation_cec_2023.tif")
dem <- dem.noram

neighbours_ls <- readRDS("BC_BGC_NeighbourList.rds")
vars.selected <- fread("no_month_vars.csv")

bgc_list <- names(neighbours_ls)
bgc_list <- bgc_list[!grepl("BAFA.* | *.un.", bgc_list)]
bgc_list <- bgc_list[-1]

tic()
res_list <- list()

bgc_list <- bgc_list[-"CWHvm2"]
for(bgc in bgc_list){
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
  # coords <- data.frame(coords, elev = temp_elev$WNA_DEM_SRT_30m, 
  #                      ID = 1:nrow(pnts_all), BGC = pnts_all$BGC)
  
  coords <- data.frame(coords, elev = temp_elev$northamerica_elevation_cec_2023, 
                       ID = 1:nrow(pnts_all), BGC = pnts_all$BGC)
  
  clim_vars <- suppressMessages(climr_downscale(coords, which_normal = "composite_normal", 
                                                vars = c(list_variables(),"CMI"), return_normal = T))  
  clim_vars <- data.table:::na.omit.data.table(clim_vars)
  clim_vars2 <- addVars(clim_vars)
  vars <- vars.selected$vars
  clim_vars3 <- clim_vars2 %>% dplyr::select(BGC, all_of(vars))
  clim_vars2[,BGC := as.factor(BGC)]
  rf_mod <- ranger(BGC ~ ., data = clim_vars2, num.trees = 201, importance = "impurity", splitrule = "gini")
  varimp <- sort(importance(rf_mod),decreasing = T)[1:10]
  res_list[[bgc]] <- data.table(Focal = bgc, Var = names(varimp), Importance = unname(varimp))
}


dat_all <- rbindlist(res_list)
fwrite(dat_all, "Focal_Variable_Importance.csv")

# toc()
# 
# 
# st_write(out, "Test_IDFdk3_Neighbours.gpkg")
# st_write(pnts_all, "Test_pnts.gpkg")
# 
# bgc_vect <- vect(bgc_map)
# focal_vect <- vect(focal)
# test <- terra::nearby(focal_vect,bgc_vect,distance = 0.1)
# plot(bgc_vect[test[,2]])