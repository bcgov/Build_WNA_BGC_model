#' Prepares coordinates and obtains climate normals
#'  using `climr_downscale`
#'
#' @param coords a `data.table` with point coordinates ("x" = longitude,
#'   "y" = latitude), elevation ("elev") and point IDs ("id").
#' @param bgcs a polygon `sf` object of biogeoclimatic zones to train the model. 
#' @param ... further arguments passed to [climr_downscale()].
#' 
#' @seealso [climr_downscale()]
#'
#' @return climate normals as a `data.table`
#' @export
#'  
#' @importFrom sf st_as_sf  st_join
#' @importFrom data.table data.table
getClimate <- function(coords, bgcs, ...) {
  dots <- list(...)
  
  coords_bgc <- st_join(coords_sf, bgcs)
  coords_bgc <- data.table(coords_bgc[,c("id","BGC")])
  coords_bgc[,geometry := NULL]
  coords_bgc <- coords_bgc[!is.na(BGC),]
  
  coords <- as.data.frame(coords)
  
  args <- append(list(coords = coords, coords_bgc = coords_bgc), dots)
  out <- do.call(.getClimVars, args) |>
    Cache(.)

  return(out)
}


#' Wrapper for `climr_downscale`
#'
#' @inheritParams getClimate 
#'
#' @return climate normals as a `data.table`.
#' 
#' @importFrom climr climr_downscale
#' @importFrom data.table setDT
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

#' Add extra climate variables
#'
#' @param dat a `data.table` with columns:
#'    * PPT05, PPT06, PPT07, PPT08, PPT09 (May, ..., September precip), 
#'    PPT_at (autumn precip), PPT_wt (winter precip)
#'    * CMD07 (July climate moisture deficit), CMD (annual CMD)
#'    * DD_0_at (autumn degree-days below 0 deg), DD_0_wt (winter degree-days below 0 deg)
#'    
#' @details This function calculates more climate variables derived from those
#'   output by `climr_downscale`. Presently it adds the following:
#'   *   May to June precip.: \eqn{PPT_MJ = PPT05 + PPT06}
#'   *   July to September precip \eqn{PPT_JAS = PPT07 + PPT08 + PPT09}
#'   *   Precipitation during vegetation dormant period: \eqn{PPT.dormant = PPT_at + PPT_wt}
#'   *   CMD deficit \eqn{CMD.def = 500 - PPT.dormant} (bounded to 0)
#'   *   \eqn{CMDMax = CMD07}
#'   *   \eqn{CMD.total = CMD.def + CMD}
#'   *   \eqn{DD_delayed = ((DD_0_at + DD_0_wt)*0.0238) - 1.8386} bounded to 0)
#'  
#' @return `dat` with all of the above added climate variables.
#' @export
#'
#' @examples
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

#' Log-transform climate variables
#'
#' @param dat a `data.table` with columns of climate variables corresponding to 
#'   the selected climate `elements`.
#' @param elements character. Climate elements to search for in `dat`.
#' @param base numeric. Log base.
#' @param add.vars logical. If `TRUE`, the new logged variables are added to `dat` 
#'   (TRUE). Otherwise, original column values are replaced with the logs (FALSE).
#' @param zero_adjust logical. If `TRUE` adjusts zeroes in raw data as:
#'   \eqn{base^{\log_base{x_min} - 1}.
#'   where \eqn{x_min} is the minimum non-zero, non-NA value.
#'
#' @details
#'   All column names that partially match strings in `elements` will be 
#'   log-transformed.
#' 
#' @return
#' @export
#'
#' @examples
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
#'
#' @param dat a `data.table` target data to "clean"
#' @param alpha numeric. The alpha value used to determine the cutoff for outliers.
#' @param IDvars character. Names of columns from which outliers should be excluded.
#' 
#' @details TODO. Parallelizes computations internally with `foreach`.
#'
#' @return
#' @seealso [foreach::foreach()]
#' 
#' @importFrom foreach foreach %do%
#' @importFrom stats mahalanobis qchisq cov
#' 
#' @export
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

rmLowSampleBGCs <- function(dat, cutoff = 30) {
  dat <- as.data.table(dat)
  BGC_Nums <- dat[,.(Num = .N), by = BGC]
  BGC_good <- dat[!BGC %in% BGC_Nums[Num < cutoff, BGC],] 
  return(BGC_good)
}