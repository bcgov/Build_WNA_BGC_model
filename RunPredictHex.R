###This script is used to predict all of BC in 400m resolution and send it to the database
## Kiri Daust, 2020

require(data.table)
require(doBy)
require(randomForest)
require(foreach)
require(dplyr)
require(reshape2)
library(doParallel)
library(tidyr)
require(ranger)
require(DBI)
require(sf)
library(here)
require(RPostgreSQL)

addVars <- function(dat){
  dat[,`:=`(PPT_MJ = PPT05+PPT06,
            PPT_JAS = PPT07+PPT08+PPT09,
            PPT.dormant = PPT_at+PPT_wt)]
  dat[,`:=`(CMD.def = 500-PPT.dormant)]
  dat[CMD.def < 0, CMD.def := 0]
  dat[,`:=`(CMDMax = CMD07,
            CMD.total = CMD.def +CMD)]
  dat[,`:=`(CMD.grow = CMD05+CMD06+CMD07+CMD08+CMD09,
            DD5.grow = DD5_05+DD5_06+DD5_07+DD5_08+DD5_09,
            DDgood = DD5 - DD18,
            DDnew = (DD5_05+DD5_06+DD5_07+DD5_08)-(DD18_05+DD18_06+DD18_07+DD18_08),
            TmaxJuly = Tmax07)]
  return(dat)
}

origVars <- c("PPT05","PPT06","PPT07","PPT08","PPT09","PPT_at","PPT_wt",
              "CMD05","CMD06","CMD07","CMD08","CMD09","CMD",
              "DD5_05","DD5_06","DD5_07","DD5_08","DD5_09","DD5","DD18",
              "DD18_05","DD18_06","DD18_07","DD18_08","Tmax07")
newVars <- c("PPT_MJ","PPT_JAS","PPT.dormant","CMD.def","CMDMax","CMD.total",
             "CMD.grow","DD5.grow","DDgood","DDnew","TmaxJuly")
##connect to database
drv <- dbDriver("PostgreSQL")
con <- dbConnect(drv, user = "postgres", host = "localhost",password = "Kiriliny41", port = 5432, dbname = "cciss_data")
##climate bc data dir
datDir <- "~/ClimBC_Tiles/CurrentData/"
load("WNAv11_35_VAR_SubZone_ranger.Rdata")

vars <- BGCmodel[["forest"]][["independent.variable.names"]]
vars2 <- vars[!vars %in% newVars]
varImport <- c("ID1","ID2",vars2,origVars) %>% unique()

##if the tiles are too big, have to split them up
tile_predict <- function(Y1, maxSize = 6000000){
  n = nrow(Y1)
  brks <- seq(1,n,by = maxSize)
  brks <- c(brks,n)
  Y1[,BGC.pred := NA_character_]
  for(j in 1:(length(brks)-1)){
    Y1[brks[j]:brks[j+1],BGC.pred := predict(BGCmodel, Y1[brks[j]:brks[j+1],-c(1:3)])[['predictions']]]
  }
  TRUE
}

## for predicting current period
pers <- c("1991_2000MSY","2001_2010MSY","2011_2019MSY")

##loop through each tile
for(i in 1:13){
  cat("Processing tile",i,"... \n")
  dat <- foreach(per = pers, .combine = rbind) %do% {
    dt <- fread(paste0(datDir,"Tile",i,"_In_Decade_",per,".csv"),select = varImport)
    dt[,Period := per]
    dt
  }
  dat <- dat[,lapply(.SD, mean), by = .(ID1,ID2), .SD = -"Period"]
  #dat <- fread(paste0(datDir,"Tile",i,"_In_Normal_1961_1990MSY.csv"),select = varImport)
  temp <- rep("Current91", nrow(dat))
  dat <- cbind(temp,dat)
  Y1 <- addVars(dat)
  ##Y1 <- dat
  
  varList = c("Model", "SiteNo", "BGC", vars)
  colnames (Y1) [1:3] = c("Model", "SiteNo", "BGC")
  Y1=Y1[,..varList]
  Y1 <- Y1[complete.cases(Y1),]
  
  ##Predict future subzones######
  tile_predict(Y1)
  gc()
  #Y1 <- separate(Y1, Model, into = c("GCM","Scenario","FuturePeriod"), sep = "_", remove = T)
  #gc()
  #Y1$FuturePeriod <- gsub(".gcm","",Y1$FuturePeriod)
  setnames(Y1, old = "Model", new = "Period")
  #Y1 <- Y1[,c("GCM","Scenario","FuturePeriod","SiteNo","BGC","BGC.pred")]
  Y1 <- Y1[,c("Period","SiteNo","BGC","BGC.pred")]
  #setnames(Y1, c("gcm","scenario","futureperiod","siteno","bgc","bgc_pred"))
  setnames(Y1, c("period","siteno","bgc","bgc_pred"))
  Y1[,varset := "Gini16"]
  dbWriteTable(con, "varset_current", Y1,row.names = F, append = T)
  rm(Y1,dat)
  gc()
    
}



####OLD####
grd <- st_read(dsn = "../BCGrid/HexGrd400.gpkg")
pts <-st_read(dsn = "../BCGrid/HexPts400.gpkg")
colnames(grd)[1] <- "id"
st_write(grd, con, drop = T)
test <- st_read(con,query = "SELECT * FROM grd WHERE id IN (5,6);")

districts <- st_read(dsn = "../CommonTables/ForestRegions.gpkg", layer = "ForestRegions_clipped")
districts <- districts[,c("REGION_NAM","ORG_UNIT")]

distJn <- st_join(pts,districts)
distJn <- distJn[!is.na(distJn$ORG_UNIT),]
dj <- st_drop_geometry(distJn) %>% as.data.table()
dj <- dj[!is.na(ORG_UNIT),]

districts <- st_read(dsn = "../CommonTables/ForestRegions.gpkg", layer = "ForestDistricts")
districts <- districts[,c("DISTRICT_N","ORG_UNIT")]
distJn2 <- st_join(pts,districts)
dj2 <- st_drop_geometry(distJn2) %>% as.data.table()
dj2 <- dj2[!is.na(ORG_UNIT),]

tID <- fread(file.choose())

dj3 <- dj2[dj, on = "ID"]
setnames(tID, c("ID","TileNum"))
dj3 <- tID[dj3, on = "ID"]
setnames(dj3,c("siteno","tileno","district","dist_code","region","reg_code"))

dbWriteTable(con,"id_atts", dj3,row.names = F)
dbGetQuery(con, "select distinct region from id_atts")

  sk <- districts[4,]

sites <- dbGetQuery(con, "select distinct siteno from cciss_400m")
sites <- dbGetQuery(con, "select * from cciss_400m 
                    inner join id_atts on cciss_400m.siteno = id_atts.siteno 
                    where reg_code = 'RSK'")
sno <- dbGetQuery(con, "select siteno from id_atts where reg_code = 'RSK'")

q <- paste0("select * from cciss_400m where siteno in (",paste(sno$siteno, collapse = ","),")")      
test <- dbGetQuery(con,q)
