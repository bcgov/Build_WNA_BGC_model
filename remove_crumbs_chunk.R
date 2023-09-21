
```{r clean crumbs}
###now cleanup and remove crumbs
library(units)
t3 <- st_cast(t2, "MULTIPOLYGON") %>% st_cast("POLYGON")
t3 <- t3 %>%
  mutate(Area = st_area(.)) %>%
  mutate(ID = seq_along(BGC))
#unique(t3$Area) ### find the size of polygon to eliminate

size <- 300000
size <- set_units(size, "m^2")
t3$Area <- set_units(size, "m^2")
tSmall <- t3[t3$Area <= size,]
t3$BGC <- as.character(t3$BGC)

require(doParallel)
coreNum <- as.numeric(detectCores()-1)
coreNo <- makeCluster(coreNum)
registerDoParallel(coreNo, cores = coreNum)

###loop through each polygon < size, determine intersects, and assign to zone with most edge touching
###all the built in functions I found only dealt with holes in the middle of polygons
i = 1
new <- foreach(i = 1:length(tSmall$ID), .combine = rbind, .packages = c("foreach","sf")) %dopar% {
  ID <- tSmall$ID[i]
  nbrs <- st_intersects(tSmall[i,],t3)[[1]]
  nbrs <- nbrs[!nbrs %in% ID]
  if(length(nbrs) == 0){return(NULL)}
  lines <- st_intersection(t3[ID,],t3[nbrs,])
  lines <- st_cast(lines)
  l.len <- st_length(lines)
  names(l.len) <- lines$BGC.1
  zn <- names(l.len)[l.len == max(l.len)][1]
  newDat <- t3[ID,]
  newDat$BGC <- zn
  newDat
}

stopCluster(coreNo)
gc()
temp <- t3[!t3$ID %in% new$ID,]
t3 <- rbind(temp, new) %>%
  mutate(Zone = as.factor(BGC))

###now have to combine crumbs with existing large polygons
temp2 <- t3
st_precision(temp2) <- 0.5
t3 <- temp2 %>%
  group_by(BGC) %>%
  summarise(geometry = sf::st_union(geometry)) %>%
  ungroup()

#mapview(t2, zcol = "BGC")
t3 <- st_zm(t3, drop=T, what='ZM')
t3 <- t3 %>% st_buffer (0) ### need st_buffer(0) to get things right
st_write(t3, dsn = paste0("./outputs/", region, "_SubZoneMap_1991_2019_eliminated.gpkg"), driver = "GPKG")

```