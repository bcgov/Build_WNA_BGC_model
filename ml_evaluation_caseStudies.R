
# Case study for evaluating BGC model sensitivities

library(terra)
library(climr)
library(ranger)
library(reproducible)

options(reproducible.cachePath = "./reproducible.cache")

# BGC color scheme
# TODO modify for direct download from CCISS_ShinyApp repo
BGCcolors.BC <- prepInputs(targetFile = "BGCzone_Colorscheme.csv", 
                        url = "https://github.com/bcgov/CCISS_ShinyApp/raw/main/data-raw/data_tables/BGCzone_Colorscheme.csv",
                        fun = "utils::read.csv")
BGCcolors <- prepInputs(targetFile = "WNAv11_Zone_Colours.csv", 
                        url = "https://github.com/bcgov/CCISS_ShinyApp/raw/main/data-raw/data_tables/WNAv11_Zone_Colours.csv",
                        fun = "utils::read.csv")
BGCcolors.subzone <- prepInputs(targetFile = "WNAv12_3_SubzoneCols.csv", 
                        url = "https://github.com/bcgov/CCISS_ShinyApp/raw/main/data-raw/data_tables/WNAv12_3_SubzoneCols.csv",
                        fun = "utils::read.csv")
BGCcolors$colour <- as.character(BGCcolors$colour)
BGCcolors$colour[match(BGCcolors.BC$zone, BGCcolors$classification)] <- as.character(BGCcolors.BC$HEX) # reset BC colors to traditional BGC zone colors
ColScheme <- BGCcolors$colour
levels.bgc <- BGCcolors.subzone[,1]
levels.zone <- BGCcolors[,1]
zone.lookup <- levels.bgc
for(i in levels.zone){ zone.lookup[grep(i,levels.bgc)] <- i }

# DEMs
dem.noram <- rast("//objectstore2.nrs.bcgov/ffec/DEM/DEM_NorAm/NA_Elevation/data/northamerica/northamerica_elevation_cec_2023.tif")
dem.bc <- rast("//objectstore2.nrs.bcgov/ffec/Climatologies/PRISM_BC/PRISM_dem/PRISM_dem.asc")
dem.wna <- rast("//objectstore2.nrs.bcgov/ffec/DEM/DEM_Composite_WNA_800m/composite_WNA_dem.tif")
dem.noram <- project(dem.noram, dem.bc)

# define the study area and dem
studyarea <- ext(c(-123, -117, 49, 52.5))

plot(dem.bc)
plot(studyarea, add=T)

# choose a larger area for truncation of training points
plot(dem.wna)
trainingarea <- ext(c(-125, -112, 43, 55))
plot(trainingarea, add=T)

# create the study area DEM
dem <- crop(dem.noram, studyarea)
plot(dem)

# select subsets of the study area for hold-outs (gaps)
centre <- c(mean(studyarea[1:2]), mean(studyarea[3:4]))
range <- c(diff(studyarea[1:2]), diff(studyarea[3:4]))
gapcentre <- matrix(c(-1,1,1,1,1,-1,-1,-1), 4, byrow=T)
gapextents <- ext(c(centre[1]+c(-1,1)/8*range[1], centre[2]+c(-1,1)/8*range[2]))
plot(gapextents, add=T)
for(i in 1:4){
  temp <- ext(c(centre[1]+sort(gapcentre[i,1]*c(1,3)/8)*range[1], centre[2]+sort(gapcentre[i,2]*c(1,3)/8)*range[2]))
  gapextents <- list(gapextents, temp)
  plot(temp,add=T)
}


