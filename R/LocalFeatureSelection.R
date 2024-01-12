library(data.table)
library(ranger)
library(foreach)
library(terra)
library(sf)
library(dplyr)
library(climr)
library(vip)

###Make 
bgc_map <- st_read("../Common_Files/WNA_BGC_v12_5Apr2022.gpkg")
bgc_info <- fread("../Common_Files/BGC_Info/WNA_BGCs_Info_v12_15.csv")
BC_BGCs <- bgc_info[grep("BGC.*",Source),BGC]
bgc_map <- bgc_map[bgc_map$BGC %in% BC_BGCs,]
#plot(bgc_map)
bgc_map <- st_cast(bgc_map, "POLYGON")
dem <- rast("../Common_Files/WNA_DEM_SRT_30m_cropped.tif")
library(tictoc)

bgc_map$ID <- seq_along(bgc_map$BGC)

neighbours_ls <- list()
for(bgc in BC_BGCs){
  cat(".")
  focal <- bgc_map[bgc_map$BGC == bgc,]
  focal <- st_union(focal$geom)
  neighbours <- st_intersects(focal, bgc_map)
  neighbours_ls[[bgc]] <- neighbours[[1]]
}

saveRDS(neighbours_ls,"BC_BGC_NeighbourList.rds")
st_write(bgc_map, "BC_BGCs_with_ID.gpkg")

require(doParallel)
coreNum <- as.numeric(detectCores()-1)
coreNo <- makeCluster(coreNum)
registerDoParallel(coreNo, cores = coreNum)

tic()
res <- foreach(bgc = names(neighbours_ls)[30:40], .combine = rbind) %do% {
  out <- bgc_map[bgc_map$ID %in% neighbours_ls[[bgc]],]
  out_union <- group_by(out, BGC) %>% 
    summarize(geom = st_union(geom),
              BGC = BGC[1])
  pnts <- st_sample(out_union, size = rep(150, nrow(out_union)), type = "random", by_polygon = T)
  pnts_all <- st_as_sf(data.frame(BGC = rep(out_union$BGC, each = 150), geom = pnts))
  pnts_all <- st_transform(pnts_all, 4326)
  coords <- st_coordinates(pnts_all)
  temp_elev <- terra::extract(dem, coords)
  coords <- data.frame(coords, elev = temp_elev$WNA_DEM_SRT_30m, 
                       ID = 1:nrow(pnts_all), BGC = pnts_all$BGC)
  clim_vars <- climr_downscale(coords, which_normal = "composite_normal", 
                               vars = c(list_variables(),"CMI"), return_normal = T)
  clim_vars2 <- clim_vars[,c("BGC",c(list_variables(),"CMI")), with = T]
  clim_vars2[,BGC := as.factor(BGC)]
  rf_mod <- ranger(BGC ~ ., data = clim_vars2, num.trees = 201, importance = "impurity", splitrule = "gini")
  varimp <- sort(importance(rf_mod),decreasing = T)[1:10]
  data.table(Focal = bgc, Var = names(varimp), Importance = unname(varimp))
}
toc()


st_write(out, "Test_IDFdk3_Neighbours.gpkg")
st_write(pnts_all, "Test_pnts.gpkg")

bgc_vect <- vect(bgc_map)
focal_vect <- vect(focal)
test <- terra::nearby(focal_vect,bgc_vect,distance = 0.1)
plot(bgc_vect[test[,2]])