### var comparison
library(data.table)
library(sf)
library(RPostgreSQL)
library(dplyr)
library(dbplyr)
library(foreach)
library(ggplot2)

##connect to local postgres instance
drv <- dbDriver("PostgreSQL")
##kiri
con <- dbConnect(drv, user = "postgres", host = "localhost",
                 password = "Kiriliny41", port = 5432, dbname = "cciss_data") 
##will
con <- dbConnect(drv, user = "postgres", host = "FLNRServer",
                 password = "Kiriliny41", port = 5432, dbname = "cciss_data") 

hexArea <- 138564.1 ##area of 1 hex poly
convFact <- hexArea/1e9 ##convert to units of 1000km^2

src_dbi(con)
datCon <- tbl(con, "varset_current") ##table

##summarise and pull predicted data
dat <- datCon %>%
  group_by(varset,bgc_pred) %>%
  summarise(Num = n()) %>%
  collect()

##summarise and pull actual (historic) data
datHist <- datCon %>%
  filter(varset == "All") %>%
  group_by(varset,bgc) %>%
  summarise(Num = n()) %>%
  collect()

datHist <- as.data.table(datHist)
datHist[,Num := Num * convFact]
datHist[,varset := "Actual"]

dat <- as.data.table(dat)
dat[,Num := Num * convFact]
setnames(dat, c("varset","bgc","Num"))
dat <- rbind(dat, datHist) ##combine
##create table by subzone
dt1 <- data.table::dcast(dat, bgc ~ varset, value.var = "Num")
fwrite(dt1,"VarComp_BC_Sz.csv")

##create table by zone
dat[,bgc := gsub("[[:lower:]]|[[:digit:]]", "",bgc)]
dat[,bgc := gsub("_.*", "",bgc)]
dt2 <- data.table::dcast(dat, bgc ~ varset, value.var = "Num", fun.aggregate = sum)
fwrite(dt2,"VarComp_BC_Zone.csv")

##graph change aspatially
dat2 <- dat[,.(Num = sum(Num)), by = .(varset,bgc)]
cutoff <- 8 ##create two graphs - cutoff for which graph
dat2[,Small := if(any(Num > cutoff)) T else F, by = .(bgc)]
datSmall <- dat2[Small == F,!"Small"]
datBig <- dat2[Small == T,!"Small"]

###graph of big units
gDat <- datBig
datAct <- gDat[varset == "Actual",]
datAct[,varset := NULL]
setnames(datAct, c("bgc","Actual"))
gDat <- gDat[varset != "Actual",]
gDat <- datAct[gDat, on = "bgc"]
gDat[,Change := Num - Actual]

ggplot(gDat, aes(x = bgc, y = Change, group = varset, fill = varset))+
  geom_bar(stat = "identity", position = position_dodge())

###graph of small units - mostly coming in from Ab or US
gDat <- datSmall
datAct <- gDat[varset == "Actual",]
datAct[,varset := NULL]
setnames(datAct, c("bgc","Actual"))
gDat <- gDat[varset != "Actual",]
gDat <- datAct[gDat, on = "bgc"]
gDat[is.na(Actual), Actual := 0]
gDat[,Change := Num - Actual]

ggplot(gDat, aes(x = bgc, y = Change, group = varset, fill = varset))+
  geom_bar(stat = "identity", position = position_dodge())

### Spatial change
dat <- datCon %>% ##pull data
  group_by(varset,bgc,bgc_pred) %>%
  summarise(Num = n()) %>%
  collect()

dat <- as.data.table(dat)
## only use zones
dat[,bgc := gsub("[[:lower:]]|[[:digit:]]", "",bgc)]
dat[,bgc_pred := gsub("[[:lower:]]|[[:digit:]]", "",bgc_pred)]
dat[,bgc_pred := gsub("_.*", "",bgc_pred)]
dt2 <- dat[,.(Num = sum(Num)), by = .(varset,bgc,bgc_pred)]
dt2[,Num := Num * convFact]
dt <- data.table::dcast(dt2, bgc + bgc_pred ~ varset, value.var = "Num")
fwrite(dt,"VarComp_All_ZoneChange.csv")

### spatial change, just for one district
District <- 'DRM'

datCon <- tbl(con, "varset_current_att")
dat2 <- datCon %>%
  filter(dist_code == District) %>%
  group_by(varset,bgc,bgc_pred) %>%
  summarise(Num = n()) %>%
  collect()

dat2 <- as.data.table(dat2)
dat2 <- dat2[Num > 100,] ## remove small
dt <- data.table::dcast(dat2, bgc + bgc_pred ~ varset, value.var = "Num")
fwrite(dt, paste0("./outputs/VarComp_",District, "_Change.csv"))


###Now spatial
datsf <- st_read(con, query = paste0("SELECT siteno,varset,bgc_pred,geom FROM vartest_sf WHERE dist_code = '", District, "'"))

fname = paste0("./outputs/",District, "_VarTests.gpkg")
for(var in unique(datsf$varset)){
  cat("Processing",var,"\n")
  temp <- datsf[datsf$varset == var,]
  mapComb <- aggregate(temp[,"geom"],  by = list(temp$bgc_pred), FUN = mean)
  colnames(mapComb)[1] <- "BGC"
  st_crs(mapComb) <- 3005
  st_write(mapComb,dsn = fname, layer = var, append = T)
}


###old####
dat <- st_drop_geometry(test) %>% as.data.table()
Num <- dat[,.(Size = .N), by = .(varset,bgc_pred)]
N2 <- Num[,.(bgc_pred,SizeProp = Size/sum(Size)), by = .(varset)]
mat <- data.table::dcast(N2, bgc_pred ~ varset)
mapComb <- aggregate(test[,"geom"],  by = list(test$bgc_pred), FUN = mean)
colnames(mapComb)[1] <- "BGC"
st_crs(mapComb) <- 3005
st_write(mapComb, dsn = "TestMap.gpkg", layer = "BioClim", append = T)

library(fasterize)
library(raster)
map <- as.data.table(mapComb) %>% st_as_sf()
##have to convert BGC labels to ints
legend <- data.table(BGC = unique(map$BGC))
legend[,Value := seq_along(BGC)]
map <- legend[map, on = "BGC"]
map <- st_as_sf(map) %>% st_cast("MULTIPOLYGON")
rast <- raster(map, resolution = 250) ##covert to raster with twice the resolution
rast <- fasterize(map,rast, field = "Value")

tMat <- as.matrix(values(rast))
Rcpp::sourceCpp("_RasterDespeckle.cpp")
tMat <- pem_focal(tMat, wrow = 7, wcol = 7) ##7*7 focal window to finding mode
tMat[tMat < 0] <- NA
values(rast) <- tMat
library(stars)
crs(rast) <- CRS("+init=EPSG:3005")
rS <- st_as_stars(rast)
poly <- st_as_sf(rS, use_integer = T, merge = T)
colnames(poly)[1] <- "Value"
poly <- merge(poly, legend, by = "Value")
st_write(poly, dsn = "TestMap.gpkg", layer = "BioClim_Despeckle2", append = T)
