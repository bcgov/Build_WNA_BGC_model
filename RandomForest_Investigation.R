##figure out random forest vars
library(ranger)
library(randomForestExplainer)
library(data.table)
library(foreach)
load("./outputs/WNA_Subzone_NoMonth.Rdata")
mod1 <- BGCmodel
load("./outputs/WNA_Subzone_BioClim.Rdata")
mod2 <- BGCmodel
load("./outputs/WNA_Subzone_17Var_extratree.Rdata")
mod3 <- BGCmodel

minDepth <- function(rfMod){
  allTrees <- foreach(t = 1:101, .combine = rbind) %do% {
    temp <- treeInfo(rfMod,tree = t)
    dat <- as.data.table(temp)
    dat2 <- dat[!is.na(splitvarName),]
    dat2[,Depth := calcDepth(nodeID,leftChild,rightChild)]
    dat2 <- dat2[,.(minDepth = min(Depth),meanDepth = mean(unique(Depth))), by = .(splitvarName)]
    dat2[,Tree := t]
    dat2
  }
  
  treeSumm <- allTrees[,.(meanMin = mean(minDepth),q5 = quantile(minDepth,0.05), 
                          q95 = quantile(minDepth,0.95)),by = .(splitvarName)]
  return(treeSumm)
}

meanSplit <- function(rfMod){
  allTreesSplit <- foreach(t = 1:101, .combine = rbind) %do% {
    temp <- treeInfo(rfMod,tree = t)
    dat <- as.data.table(temp)
    dat2 <- dat[!is.na(splitvarName),]
    dat2[,Depth := calcDepth(nodeID,leftChild,rightChild)]
    dat2 <- dat2[,.(splitval = splitval[Depth == min(Depth)]), by = .(splitvarName)]
    dat2[,Tree := t]
    dat2
  }
  
  valSumm <- allTreesSplit[,.(meanSplit = mean(splitval),q5 = quantile(splitval,0.05), 
                              q95 = quantile(splitval,0.95)),by = .(splitvarName)]
  return(valSumm)
}

m1Depth <- minDepth(mod1)
m1Split <- meanSplit(mod1)
m2Depth <- minDepth(mod2)
m2Split <- meanSplit(mod2)
m3Depth <- minDepth(mod3)
m3Split <- meanSplit(mod3)

bioClimStats <- m1Depth[m1Split, on = "splitvarName"]
bioClimStats[,Model := "BioClim"]
noMonthStats <- m2Depth[m2Split, on = "splitvarName"]
noMonthStats[,Model := "NoMonth"]
new17Stats <- m3Depth[m3Split, on = "splitvarName"]
new17Stats[,Model := "17Var"]
allStats <- rbind(new17Stats,noMonthStats,bioClimStats)
allStats[,`:=`(i.q5 = NULL,i.q95 = NULL)]
fwrite(allStats,"RangerMinVarDepth.csv")

splitFactor <- function(rfMod,focalUnit){
  idfOut <- foreach(t = 1:101, .combine = rbind) %do% {
    cat(".")
    temp <- treeInfo(rfMod,tree = t)
    temp$prediction <- as.character(temp$prediction)
    starts <- temp$nodeID[grep(focalUnit,temp$prediction)]
    rnodes <- findSplit(temp$nodeID,temp$leftChild,temp$rightChild,initNode = starts,cutoff = length(starts))
    if(length(rnodes) > 0){
      rnodes <- unique(rnodes)
      byRoot <- foreach(rnode = rnodes, .combine = rbind) %do% {
        roots = c(temp$leftChild[rnode+1],temp$rightChild[rnode+1])
        tempres <- 1
        for(r in roots){
          testout <- capture.output(subtreePreds(temp$nodeID,temp$leftChild,temp$rightChild,temp$prediction,root = r))
          testout <- strsplit(testout," +")
          testout <- testout[[1]]
          tab <- table(testout)
          tab <- tab[names(tab) %in% c("SBSmc2","IDFdk1")]
          tempres <- tempres*(min(tab)/max(tab))
        }
        t1 <- temp[temp$nodeID == rnode,]
        t1$Clean <- tempres
        t1
      }
      
      out <- temp[temp$nodeID %in% rnodes,]
      out$Tree <- t
      out
    }else{
      NULL
    }
  }
  freq <- table(idfOut$splitvarName)
  freq <- freq[freq > quantile(freq, 0.5)]
  return(freq)
}
sf17var <- splitFactor(mod3,focalUnit = "IDFdk1")
sf17var <- as.data.table(sf17var)
setnames(sf17var, old = "N", new = "Var17")
sfBioClim <- splitFactor(mod2, focalUnit = "IDFdk1")
sfBioClim <- as.data.table(sfBioClim)
setnames(sfBioClim, old = "N", new = "BioClim")
sfNoMonth <- splitFactor(mod1, focalUnit = "IDFdk1")
sfNoMonth <- as.data.table(sfNoMonth)
setnames(sfNoMonth, old = "N", new = "NoMonth")

sfall <- merge(sf17var, sfBioClim, by = "V1", all = T)
sfall <- merge(sfall, sfNoMonth, by = "V1", all = T)
fwrite(sfallCWH, "SplittingNodesCWHws2.csv")
fwrite(sfall, "SplittingNodesIDFdk1.csv")

noMonthTest <- table(idfOut$splitvarName)
noMonthTest <- noMonthTest[noMonthTest > 2]
bioClimTest <- bioClimTest[bioClimTest > 7]

temp <- treeInfo(mod3,tree = 1)
temp$prediction <- as.character(temp$prediction)
starts <- temp$nodeID[grep("CWHws2",temp$prediction)]
rnodes <- findSplit(temp$nodeID,temp$leftChild,temp$rightChild,initNode = starts,cutoff = length(starts))
rnodes <- unique(rnodes)

roots = c(temp$leftChild[rnodes[2]+1],temp$rightChild[rnodes[2]+1])
#roots <- c(309,310)
for(r in roots){
  testout <- capture.output(subtreePreds(temp$nodeID,temp$leftChild,temp$rightChild,temp$prediction,root = r))
  testout <- strsplit(testout," +")
  testout <- testout[[1]]
  tab <- table(testout)
  tab <- tab[names(tab) %in% c("IDFdk1","CWHws2")]
  print(tab)
}

testout <- capture.output(subtreePreds(temp$nodeID,temp$leftChild,temp$rightChild,temp$prediction,root = 0))
testout <- strsplit(testout," +")
testout <- testout[[1]]
tab <- table(testout)
tab <- tab[names(tab) %in% c("IDFdk1","CWHws2")]


rf <- ranger(Species ~ ., data = iris)
temp <- treeInfo(rf, 1)
temp$prediction <- as.character(temp$prediction)
starts <- temp$nodeID[grep("virginica",temp$prediction)]
starts <- sort(starts,decreasing = T)
test <- findSplit(temp$nodeID,temp$leftChild,temp$rightChild,initNode = starts,cutoff = 5)
