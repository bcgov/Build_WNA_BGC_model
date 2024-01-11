getClimate <- function(coords, bgcs, ...) {
  dots <- list(...)
  
  coords_bgc <- st_join(coords_sf, bgcs)
  coords_bgc <- data.table(coords_bgc[,c("id","BGC")])
  coords_bgc[,geometry := NULL]
  coords_bgc <- coords_bgc[!is.na(BGC),]
  
  #coords <- fread("WNA_2km_grid_WHM.csv")
  #setcolorder(coords, c("long","lat","elev","id"))
  coords <- as.data.frame(coords)# %>% dplyr::rename(long = 1, lat = 2)
  
  ##turned it into one function
  args <- append(list(coords = coords, coords_bgc = coords_bgc), dots)
  out <- do.call(.getClimVars, args) |>
    Cache(.)

  return(out)
}

.getClimVars <- function(coords, coords_bgc, ...) {
  clim_vars <- climr_downscale(coords, ...)
  setDT(clim_vars)
  clim_vars <- clim_vars[!is.nan(PPT05),] ##lots of points in the ocean
  clim_vars[coords_bgc, BGC := i.BGC, on = "ID==id"]
  clim_vars <- clim_vars[!is.na(BGC), ]
  clim_vars[, PERIOD := NULL]
  clim_vars[, ID := NULL]
  
  return(clim_vars)
}

addVars <- function(dat) {
  dat[, PPT_MJ := PPT05 + PPT06]
  dat[, PPT_JAS := PPT07 + PPT08 + PPT09]
  dat[, PPT.dormant := PPT_at + PPT_wt]
  dat[, CMD.def := 500 - PPT.dormant]
  dat[CMD.def < 0, CMD.def := 0]
  dat[, CMDMax := CMD07]
  dat[, CMD.total := CMD.def + CMD]
  dat[, DD_delayed := ((DD_0_at + DD_0_wt)*0.0238) - 1.8386]
  dat[DD_delayed < 0, DD_delayed := 0]
}

# Function to log-transform fields in a data.table (e.g., climr output)
logVars <- function(dat,
                    elements = c("AHM", "DD", "Eref", "FFP", "NFFD", "PAS", "PPT", "SHM", "CMI"),
                    base = exp(1),
                    add.fields = FALSE,
                    zero_adjust = FALSE) {
  
  dat <- copy(dat)
  
  # Fields to operate on (generally these should be ratio (zero-limited) variable)
  logFields <- grep(paste(elements, collapse = "|"), names(dat), value = TRUE)
  dat.log <- dat[, .SD, .SDcols = logFields]
  
  # If specified by the user, give zero values a positive value that is one order of magnitude less than the minimum positive value
  if (zero_adjust) {
    dat.log[, lapply(.SD, function(x) {
      x[x <= 0] <- base^(log(min(x[x > 0], na.rm = TRUE), base = base) - 1)
      return(x)
    })]
  }
  
  # Perform log transformation
  dat.log <- dat.log[, lapply(.SD, function(x) log(x, base = base))]
  
  # Add 
  if(add.fields){
    setnames(dat.log, logFields, paste0(logFields, "_log"))
    dat <- cbind(dat, dat.log)
  } else {
    dat[, (logFields) := Map(x =.SD, xname = logFields, f = function(x, xname) {
      x <- dat.log[[xname]]
      return(x)
    }), .SDcols = logFields]
  }
  return(dat)
}


removeOutlier <- function(dat, alpha,numIDvars){
  out <- foreach(curr = unique(as.character(dat$BGC)), .combine = rbind) %do% {
    temp <- dat[dat$BGC == curr,]
    md <- tryCatch(mahalanobis(temp[,-c(1:numIDvars)],
                               center = colMeans(temp[,-c(1:numIDvars)]),
                               cov = cov(temp[,-c(1:numIDvars)])), error = function(e) e)
    if(!inherits(md,"error")){
      ctf <- qchisq(1-alpha, df = ncol(temp)-1)
      outl <- which(md > ctf)
      message("Removing", length(outl), "outliers from",curr, "; ")
      if(length(outl) > 0){
        temp <- temp[-outl,]
      }
    }
    temp
    
  }
  return(out)
}

##remove BGCs with low numbers of points
rmLowSampleBGCs <- function(X2, cutoff = 30) {
  XAll <- as.data.table(X2)
  BGC_Nums <- XAll[,.(Num = .N), by = BGC]
  BGC_good <- XAll[!BGC %in% BGC_Nums[Num < cutoff, BGC],] 
  return(BGC_good)
}