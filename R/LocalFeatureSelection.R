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

bgc_map <- st_read("D:/OneDrive - Government of BC/WNA_BGC - Future Forest Ecosystems Centre/ccissv13_files/WNA_BGC_v13_15Nov2024.gpkg") %>% dplyr::filter(!BGC == "(None)")
bgc_map$ID <- seq_along(bgc_map$BGC)
 bgcs <- unique(bgc_map$BGC) %>% data.frame
 fwrite(bgcs, "wna_bgcs.csv")
 
bgc_info <- fread("D:/OneDrive - Government of BC/WNA_BGC - Future Forest Ecosystems Centre/ccissv13_files/WNA_BGCs_Info_v13_1.csv") %>% filter(BGC %in% bgcs$.) #%>% mutate(bgc = as.factor(BGC))
BC_BGCs <- bgc_info[grep("BGC.*",Source),BGC]

#bgc_map <- st_read("BC_BGCs_with_ID.gpkg")

neighbours_ls <- list()
#for(bgc in BC_BGCs[-(1:180)]){
for(bgc in bgc_info$BGC){
  cat(".")
  focal <- bgc_map[bgc_map$BGC == bgc,]
  if(nrow(focal) > 0){
    focal <- st_union(focal$geom)
    neighbours <- st_intersects(focal, bgc_map)
    neighbours_ls[[bgc]] <- neighbours[[1]]
  }
}

saveRDS(neighbours_ls,"WNA_BGC_NeighbourList.rds",)

st_write(bgc_map, "WNA_BGCs_with_ID.gpkg", append = FALSE)

#dem <- terra::rast("D:/CommonTables/DEMs/WNA_DEM_4326_clipped.tif")
saveRDS(dem, "//objectstore2.nrs.bcgov/ffec/CCISS_Working/WNA_DEM_4326_clippedDEM.rds")
# %>% terra::project("epsg:4326")
# writeRaster(dem2, "D:/CommonTables/DEMs/WNA_DEM_4326_clipped.tif")
#dem <- rast("//objectstore2.nrs.bcgov/ffec/DEM/DEM_NorAm/NA_Elevation/data/northamerica/northamerica_elevation_cec_2023.tif")

# dem <- terra::rast("D:/CommonTables/DEMs/WNA_DEM_SRT_30m_cropped.tif")
# dem2 <- crop(dem, ext(st_transform(bgc_map, 4326)))
#dem.noram <- rast("//objectstore2.nrs.bcgov/ffec/DEM/DEM_NorAm/NA_Elevation/data/northamerica/northamerica_elevation_cec_2023.tif")
#dem <- dem.noram

neighbours_ls <- readRDS("WNA_BGC_NeighbourList.rds")
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
bgc = "CDFmm"
for(bgc in bgc_list){
  cat("Processing",bgc,"\n")
  out <- bgc_map[bgc_map$ID %in% neighbours_ls[[bgc]],]
  out_union <- group_by(out, BGC) %>% 
    summarize(geom = st_union(geom),
              BGC = BGC[1])
  pnts <- st_sample(out_union, size = rep(150, nrow(out_union)), type = "random", by_polygon = T)
  pnts_all <- st_as_sf(data.frame(BGC = rep(out_union$BGC, each = 150), geometry = pnts))
  pnts_all <- st_transform(pnts_all, 4326)
  coords <- st_coordinates(pnts_all)
  temp_elev <- terra::extract(dem, coords)
  coords.bgc <- data.frame(coords, elev = temp_elev$WNA_DEM_SRT_30m,
                       ID = 1:nrow(pnts_all), BGC = pnts_all$BGC) %>% rename(lat = Y, lon = X, id = ID) %>% 
    dplyr::select(lon, lat, elev,id, BGC)
  coords <- coords.bgc %>%   dplyr::select(lon, lat, elev,id)
  coords.bgc <- coords.bgc %>% dplyr::select(id, BGC)
  # coords <- data.frame(coords, elev = temp_elev$WNA_DEM_SRT_30m, 
  #                      ID = 1:nrow(pnts_all), BGC = pnts_all$BGC)
  
  clim_vars <- suppressMessages(climr_downscale(coords, which_normal = "auto", 
                                                vars = list_variables(), return_normal = T))  
  clim_vars <- data.table:::na.omit.data.table(clim_vars)
  addVars(clim_vars)
  clim_vars <- left_join(clim_vars, coords.bgc)
  
  clim_vars <- setDT(clim_vars)[,..vars]
  
  clim_vars[,BGC := as.factor(BGC)]
  rf_mod <- ranger(BGC ~ ., data = clim_vars, num.trees = 101, importance = "impurity", splitrule = "gini")
  varimp <- sort(importance(rf_mod),decreasing = T)[1:6]
  res_list[[bgc]] <- data.table(Focal = bgc, Var = names(varimp), 
                                Importance = unname(varimp), 
                                OOB = rf_mod$prediction.error,
                                NumberBGCs = nrow(out_union))
}


dat_all2 <- rbindlist(res_list)
fwrite(dat_all2, "Focal_Variable_Importance_WNAv1.csv")
var.count <- dat_all2[,.(Num = .N), by = .(Var)]
setorder(test, -Num)
toc()
fwrite(var.count, "Count_of_Variables_v2.csv")
#### Add type of variable by 3 categories and look by BGC
### GS/NGS, Temp?Precip?interact, Extreme/Mean

mean.error <- dat_all2 %>% select(Focal, OOB ) %>% distinct %>% mutate(mean.error = mean(OOB))
mean.error2 <- mean.error %>% filter(!str_detect(Focal, "p$"))%>% filter(!str_detect(Focal, "w$"))
require(ggplot2);require(stringr)


ggplot(mean.error2, aes(x="Y" , y=OOB))+
  geom_violin()+
  stat_summary(fun.data = "mean_sdl",  geom="crossbar", width=0.05)

var.type <- fread("Variable_Types.csv")
dat_all2 <- fread("Focal_Variable_Importance_WNAv1.csv") %>% left_join(var.type)
BGC_using_group <- dat_all2 %>% group_by(Focal, var.group, var.season) %>% count() %>% select(-n) %>%  ungroup() %>% group_by(var.group, var.season) %>% count()
var.count <- dat_all2[,.(Num = .N), by = .(var.type)]
setorder(test, -Num)

# 
# st_write(out, "Test_IDFdk3_Neighbours.gpkg")
# st_write(pnts_all, "Test_pnts.gpkg")
# 
# bgc_vect <- vect(bgc_map)
# focal_vect <- vect(focal)
# test <- terra::nearby(focal_vect,bgc_vect,distance = 0.1)
# plot(bgc_vect[test[,2]])