library(raster)
library(rgdal)
library(dismo)
library(maptools)
library(ggplot2)

setwd("USER_WORKING_DIRECTORY/data/")

#Writen by Hannah Lois Owens <hannah.owens@sund.ku.dk>

# Get sites
sites <- readOGR("~/Dropbox/PaleoMegafauna/data/arctic_sites/arctic.shp")
sites <- spTransform(sites, 
                     CRS("+proj=laea +lat_0=90 +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m"))

templateForDownsamplingOccs <- raster("Inputs/PaleoClimateRasters/X005/bio_1.asc")

# Extract data for each site and calculate mean ----
humanBinary <- stack(list.files("Outputs/SDM_Redo/PostProcessedLayers/Binary/", 
                                pattern = ".asc", full.names = T))
crs(humanBinary) <- CRS("+proj=longlat +datum=WGS84")
humanBinary <- projectRaster(humanBinary, 
                             crs = CRS("+proj=laea +lat_0=90 +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m"),
                             method = "ngb")
humanBinaryExt <- extract(x = humanBinary, y = sites, fun = mean)

humanContinuous <- stack(list.files("Outputs/SDM_Redo/PostProcessedLayers/Continuous/", pattern = ".asc", full.names = T))
crs(humanContinuous) <- CRS("+proj=longlat +datum=WGS84")
humanContinuous <- projectRaster(humanContinuous, 
                             crs = CRS("+proj=laea +lat_0=90 +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m"))
humanContinuousExt <- extract(x = humanContinuous, y = sites, fun = mean)

# Collating and writing data ----
cols <- paste0(gsub(colnames(humanBinaryExt), pattern = "ensProjNoClampBinary_", replacement = ""),"bp")

humanBinaryExt <- round(humanBinaryExt)
colnames(humanBinaryExt) <- cols
humanBinaryExt <- cbind(sites@data, humanBinaryExt)
write.csv(humanBinaryExt, file = "Outputs/SDM_Redo/humanPresenceAbsence.csv", row.names = F)

colnames(humanContinuousExt) <- cols
humanContinuousExt <- cbind(sites@data, humanContinuousExt)
write.csv(humanContinuousExt, file = "Outputs/SDM_Redo/humanContinuousSuitability.csv", row.names = F)

# Mapping results ----
occurrences <- read.csv(paste0("USER_WORKING_DIRECTORY/data/Inputs/",
                               "HomoSapiensOccurrences_FinalForModeling_FINAL.csv"), 
                        header = T, stringsAsFactors = F)
occurrences <- occurrences[-4340:-4342,] # These points were breaking spTransform.
occurrences <- SpatialPointsDataFrame(occurrences[,c("Longitude", "Latitude")], 
                                      data = as.data.frame(occurrences$TimeStepIndex), 
                                      proj4string = CRS("+proj=longlat +datum=WGS84"))
occurrences <- spTransform(occurrences, 
                           CRS("+proj=laea +lat_0=90 +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m"))
names(occurrences) <- "Interval"
sliceIndex <- as.numeric(gsub(x = gsub(x = names(humanBinary), 
                                       pattern = "ensProjNoClampBinary_", 
                                       replacement = ""), 
                              pattern = "k", 
                              replacement = ""))
gridlineReferencePolygon <- rasterToPolygons(humanBinary[[1]], dissolve = T)

pdf("Outputs/SDM_Redo/PresenceAbsencePlot.pdf", width = 6, height = 6)
for(index in 1:nlayers(humanBinary)){
  sliceOccs <- occurrences[occurrences$Interval==sliceIndex[index],]
  sliceOccs <- gridSample(sliceOccs, humanBinary[[index]])
  sliceOccs <- sliceOccs[!is.na(extract(x = humanBinary[[index]], sliceOccs)),]
  plot(humanBinary[[index]], col = c("gray", "gold"), 
       main = paste0("Human Presence/Absence, ", gsub(names(humanBinary[[index]]), 
                                                      pattern = "ensProjNoClampBinary_", 
                                                      replacement = ""), "ybp"),
       legend = F, axes = F, useRaster=FALSE)
  llgridlines(gridlineReferencePolygon, lty = 1, 
              easts = c(-120, -60, 0, 60, 120), 
              norths = c(45, 55, 65, 75))
  points(sliceOccs, pch = 20, col = "black", cex = 0.75)
  plot(sites, border = "firebrick", cex = 2, add = T)
  legend(x = 'topright', legend = c("Present", "Absent", "Occurrence", "Paleo Site"), 
         fill = c("gold", "gray", "black", "firebrick"), bty = "n", cex = .75, horiz = T)
}
dev.off()

pdf("Outputs/SDM_Redo/ContinuousPlot.pdf", width = 6, height = 6)
for(index in 1:nlayers(humanContinuous)){
  sliceOccs <- occurrences[occurrences$Interval==sliceIndex[index],]
  sliceOccs <- gridSample(sliceOccs, humanBinary[[index]], n = 1)
  sliceOccs <- sliceOccs[!is.na(extract(x = humanBinary[[index]], sliceOccs)),]
  plot(humanContinuous[[index]], col = viridis::cividis(10), 
       main = paste0("Relative human suitability, ", gsub(names(humanBinary[[index]]), 
                                                      pattern = "ensProjNoClampBinary_", 
                                                      replacement = ""), "ybp"),
       legend = F, axes = F, useRaster=T)
  llgridlines(gridlineReferencePolygon, lty = 1,
              easts = c(-120, -60, 0, 60, 120), 
              norths = c(45, 55, 65, 75))
  points(sliceOccs, pch = 20, col = "black", cex = 0.75)
  plot(sites, border = "red", add = T)
  legend(x = 'topright', legend = c("Occurrence", "Paleo Site"), 
         fill = c("black", "red"), bty= "n", cex = .75, horiz = T)
}
dev.off()

# Occurrence count plot
#Load continents
cont <- readOGR("~/Dropbox/PaleoMegafauna/data/continent/continent.shp")
occurrences <- read.csv(paste0("USER_WORKING_DIRECTORY/data/Inputs/",
                               "HomoSapiensOccurrences_FinalForModeling_FINAL.csv"), 
                        header = T, stringsAsFactors = F)
occs <- cbind(occurrences, extract(cont, occurrences[,c("Longitude", "Latitude")]))
humanBinary <- stack(list.files("Outputs/SDM_Redo/PostProcessedLayers/Binary/", pattern = ".asc", full.names = T))
crs(humanBinary) <- CRS("+proj=longlat +datum=WGS84")

allOccs <- as.vector(NULL)
uniqueOccs <- as.vector(NULL)
naOccs <- as.vector(NULL)
eOccs <- as.vector(NULL)
aOccs <- as.vector(NULL)
bins <- sort(unique(occurrences$TimeStepIndex)) 
for (x in bins){
  print(x)
  occ_per_bin <- occurrences[occurrences$TimeStepIndex==x,]	 # subset data.frame for each time bin
  allOccs <- c(allOccs, nrow(occ_per_bin)) # number of Occurrences for each time interval
  occ_per_bin <- gridSample(occ_per_bin[,c("Longitude", "Latitude")], 
                            templateForDownsamplingOccs, n = 1)
  uniqueOccs <- c(uniqueOccs, nrow(occ_per_bin))
  if (nrow(occ_per_bin) > 0) {
    contOccs <- extract(cont, occ_per_bin)
    naOccs <- c(naOccs, nrow(contOccs[contOccs$CONTINENT=="North America",]))
    eOccs <- c(eOccs, nrow(contOccs[contOccs$CONTINENT=="Europe",]))
    aOccs <- c(aOccs, nrow(contOccs[contOccs$CONTINENT=="Asia",]))
    print(nrow(occ_per_bin))
  }
}
pdf("Outputs/SDM_Redo/HumanOccurrenceCounts.pdf", width = 6, height = 6)
barplot(allOccs, names.arg = bins, main = "Paleo Human Occurrences")
barplot(uniqueOccs, add = T, col = "black")
abline(a = 8, b = 0, col = "red")
legend(x = 'topright', legend = c("Raw Occurrences", "Unique Occurrences"), 
       fill = c("gray", "black"), bty= "n", cex = 1, horiz = F)
dev.off()
pdf("Outputs/SDM_Redo/occurrenceLocation.pdf")
# Stacked
data <- rbind(naOccs, eOccs, aOccs)
colnames(data) <- bins[1:ncol(data)]
barplot(data, 
        col=c("gray50","gray25","black"), 
        border="white", 
        xlab="KBP",
        main = "Unique Points By Continent")
legend(x = "topright", legend = c("North America", "Europe", "Asia"), 
       fill=c("gray50","gray25","black"), 
       bty= "n")
dev.off()