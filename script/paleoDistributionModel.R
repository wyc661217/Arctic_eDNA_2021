library(biomod2)
library(dismo)
library(raster)
library(rangeBuilder)
library(dplyr)

#Writen by Hannah Lois Owens <hannah.owens@sund.ku.dk>


workingDirectory <- "USER_WORKING_DIRECTORY/data/Inputs/"
wd <- workingDirectory
templateForDownsamplingOccs <- raster(paste0(workingDirectory, "/PaleoClimateRasters/X005/bio_1.asc"))
# Load and process occurrence data ----
# Load, calculate, and filter uncertainty and calculate time step index
occurrences <- read.csv(paste0(workingDirectory,"HomoSapiensOccurrences_FinalForModeling_Redo95.csv"), 
                        header = T, stringsAsFactors = F, sep = ';')
occurrences <- occurrences[occurrences$Latitude != 0,]
occurrences$uncertaintyRange <- occurrences$Max_of_the_probability_interval - 
  occurrences$Min_of_the_probability_interval
occurrences$TimeStepIndex <- round(((occurrences$Max_of_the_probability_interval - 
                                       occurrences$Min_of_the_probability_interval)/2 + 
                                      occurrences$Min_of_the_probability_interval - 500)/1000,
                                   0)
occurrences <- filter(occurrences, (uncertaintyRange < 1000 & TimeStepIndex < 32) | (uncertaintyRange < 2000 & TimeStepIndex >= 32))
occurrences$TimeStepIndex <- gsub(occurrences$TimeStepIndex, pattern = "^32", replacement = 32)
occurrences$TimeStepIndex <- gsub(occurrences$TimeStepIndex, pattern = "^33", replacement = 32)
occurrences$TimeStepIndex <- gsub(occurrences$TimeStepIndex, pattern = "^34", replacement = 34)
occurrences$TimeStepIndex <- gsub(occurrences$TimeStepIndex, pattern = "^35", replacement = 34)
occurrences$TimeStepIndex <- gsub(occurrences$TimeStepIndex, pattern = "^36", replacement = 36)
occurrences$TimeStepIndex <- gsub(occurrences$TimeStepIndex, pattern = "^37", replacement = 36)
occurrences$TimeStepIndex <- gsub(occurrences$TimeStepIndex, pattern = "^38", replacement = 38)
occurrences$TimeStepIndex <- gsub(occurrences$TimeStepIndex, pattern = "^39", replacement = 38)
occurrences$TimeStepIndex <- gsub(occurrences$TimeStepIndex, pattern = "^40", replacement = 40)
occurrences$TimeStepIndex <- gsub(occurrences$TimeStepIndex, pattern = "^41", replacement = 40)
occurrences$TimeStepIndex <- gsub(occurrences$TimeStepIndex, pattern = "^42", replacement = 42)
occurrences$TimeStepIndex <- gsub(occurrences$TimeStepIndex, pattern = "^43", replacement = 42)
occurrences$TimeStepIndex <- gsub(occurrences$TimeStepIndex, pattern = "^44", replacement = 44)
occurrences$TimeStepIndex <- gsub(occurrences$TimeStepIndex, pattern = "^45", replacement = 44)
occurrences$TimeStepIndex <- gsub(occurrences$TimeStepIndex, pattern = "^46", replacement = 46)
occurrences$TimeStepIndex <- gsub(occurrences$TimeStepIndex, pattern = "^47", replacement = 46)

occurrences$TimeStepIndex <- as.numeric(occurrences$TimeStepIndex)

write.csv(occurrences, paste0(workingDirectory,"HomoSapiensOccurrences_FinalForModeling_FINAL.csv"))

allOccs <- as.vector(NULL)
uniqueOccs <- as.vector(NULL)
bins <- sort(unique(occurrences$TimeStepIndex)) 
for (x in bins){ 
	 occ_per_bin <- occurrences[occurrences$TimeStepIndex==x,]	 # subset data.frame for each time bin
	 allOccs <- c(allOccs, nrow(occ_per_bin)) # number of Occurrences for each time TimeStepIndices
	 occ_per_bin <- gridSample(occ_per_bin[,c("Longitude", "Latitude")], 
	                           templateForDownsamplingOccs, n = 1)
	 uniqueOccs <- c(uniqueOccs, nrow(occ_per_bin))
	 rm(occ_per_bin)
}
occStats <- data.frame("Time"= paste(bins,"k", sep=""), "All"=allOccs, "Unique"=uniqueOccs)
n_var <- 4 # Set number of explanatory variables
list_len <- length(uniqueOccs[uniqueOccs >= 2*n_var + 1]) # count time TimeStepIndices w/ at least 2*n_var + 1 occurrences
sdm_bins <- bins[uniqueOccs >= 2*n_var] # Select time TimeStepIndices w/ enough occurrences to run analysis
sdm_bins <- sdm_bins[sdm_bins > 11]
used_occ <- uniqueOccs[uniqueOccs >= 2*n_var] # Get counts for

# Plot points by time slice as a first check ----
pdf(paste0(workingDirectory, "pointPlots_Redo.pdf"))
for (x in 1:length(bins)){
  occs <- occurrences[occurrences$TimeStepIndex==bins[[x]],c("Longitude", "Latitude")]
  if (nrow(occs) > 0){
    plot(gshhs, col = "yellow", xlim = c(-180,180), ylim = c(40,90), main = paste0("Time: ", bins[x], "k, n = ", 
                                                                                 nrow(occs)))
    points(occs, pch = 20, col = "red")
  }
}
dev.off()

# Input settings for biomod ----
myBiomodOption <- Print_Default_ModelingOptions()
myBiomodOption@MAXENT.Phillips$path_to_maxent.jar = paste(system.file(package="dismo"), 
                                                          "/java", sep='');
myBiomodOption@MAXENT.Phillips$memory_allocated <- 1024
myBiomodOption@MAXENT.Phillips$maximumiterations = 10000;
myBiomodOption@MAXENT.Phillips$threshold = F;
myBiomodOption@MAXENT.Phillips$hinge = F;
myBiomodOption@MAXENT.Phillips$visible = F;
myBiomodOption@MAXENT.Phillips$beta_lqp = .95;

myBiomodOption@GAM$algo <- "BAM_mgcv"
myBiomodOption@GAM$control$maxit <- 100

# Model each time slice that has enough occurrences ----
setwd(paste0(workingDirectory, "/PaleoClimateRasters/"))
dirs <- list.dirs()[-1]
dirs <- dirs[grepl(dirs, pattern = "X")]
sdm_dirs <- unlist(lapply(as.character(sdm_bins), FUN = function(x) dirs[grepl(pattern = x, x = dirs)]))
for (index in 1:length(sdm_dirs)) {
  print(paste0("Now modeling ", sdm_dirs[[index]], "k."))
  print(Sys.time())
  clList <- list.files(path = sdm_dirs[index], pattern = "bio", full.names = T)
  climate <- stack(clList)
  raster::crs(climate) <- raster::crs("+proj=longlat +datum=WGS84")

  #Get the occurrences per time TimeStepIndex
  sliceOccs <- occurrences[occurrences$TimeStepIndex==sdm_bins[index],]
  sliceOccs <- gridSample(sliceOccs[, c("Longitude", "Latitude")], templateForDownsamplingOccs, n = 1)

  # Format the data 
  occData <- BIOMOD_FormatingData(resp.var = rep(1, nrow(sliceOccs)), 
                                  expl.var = climate, 
                                  resp.xy = sliceOccs[, c("Longitude", "Latitude")], 
                                  resp.name = "Homo_Sapiens", 
                                  PA.nb.rep = 5, 
                                  PA.nb.absences = 5 * nrow(sliceOccs), 
                                  PA.strategy = "disk", 
                                  PA.dist.min = 150000, 
                                  PA.dist.max = 1250000)
  plot(occData)

  # Run SDM algorithms
  setwd(paste0("USER_WORKING_DIRECTORY/data/Outputs/", "SDM_Redo/"))
  models <- BIOMOD_Modeling(data = occData, 
                            models = c("MAXENT.Phillips", "GLM", "GAM", "ANN", "RF"), 
                            models.options = myBiomodOption, NbRunEval = 5, DataSplit = 80, 
                            VarImport = 3, models.eval.meth = c("TSS", "ROC"), 
                            modeling.id = paste0(occData@sp.name, sdm_bins[[index]], "k", sep=""),
                            compress = T, do.full.models = FALSE)

  # Get SDM evaluation scores and plot the results
  model_scores <- get_evaluations(models, as.data.frame=T)

  pdf(paste0("Plots_", occData@sp.name, "_", sdm_bins[index], "k.pdf", sep=""))
  models_scores_graph(models, by = "models" , metrics = c("ROC","TSS"), 
                      xlim = c(0,1), ylim = c(0.5,1), main=paste0("Model_scores_", sdm_bins[index], "k", sep=""))

  models_scores_graph(models, by="cv_run" , metrics = c("ROC","TSS"), xlim = c(0,1), ylim = c(0.5,1), 
                      main=paste0("CV_run_", sdm_bins[index], "k", sep=""))

  models_scores_graph(models, by="data_set" , metrics = c("ROC","TSS"), xlim = c(0,1), ylim = c(0.5,1), 
                      main=paste0("PA_Dataset_", sdm_bins[index], "k", sep=""))
  dev.off()
  
  models_var_import <- get_variables_importance(models,as.data.frame=T)
  var.imp <- apply(models_var_import, c(1,2), mean)
  write.csv(var.imp, paste0("VarImp_", occData@sp.name, "_", sdm_bins[index], "k.csv", sep=""))

  # Construct ensemble models
  ensemble_models <- BIOMOD_EnsembleModeling(modeling.output = models, em.by = "all", eval.metric = "ROC", 
                                             eval.metric.quality.threshold = 0.8, models.eval.meth = c("TSS", "ROC"), 
                                             prob.mean = T, prob.cv = TRUE, committee.averaging = FALSE, 
                                             prob.mean.weight = TRUE)
  ensemble_models_scores <- get_evaluations(ensemble_models)
  write.csv(ensemble_models_scores, paste0("EnsScores_", occData@sp.name, "_", sdm_bins[index], "k.csv", sep=""))
  # Do projections
  models_proj <- BIOMOD_Projection(modeling.output = models, new.env = climate,
                                   proj.name = paste0(sdm_bins[index], "k_proj_ens"),
                                   build.clamping.mask = T,
                                   output.format = ".grd", do.stack = FALSE)
  
  ensemble_models_proj <- BIOMOD_EnsembleForecasting(EM.output = ensemble_models, projection.output = models_proj, 
                             proj.name = paste0(sdm_bins[index], "k_proj_ens"),
                             output.format = ".grd", do.stack = FALSE)
  
  # Reset working directory to do it all over
  setwd(paste0(workingDirectory, "/PaleoClimateRasters/"))
}

# Get and process predictions ----
setwd(paste0("USER_WORKING_DIRECTORY/data/Outputs/", "SDM_Redo/"))
resultDirs <- dir(path = paste0("USER_WORKING_DIRECTORY/data/Outputs/", 
                                "SDM_Redo/Homo.Sapiens"), 
                  pattern = "k_proj", full.names = T)
threshFiles <- list.files(".", pattern = "EnsScores", recursive = F, full.names = T)
threshVals <- vector(mode = "list", length = length(resultDirs))
TSSvals <- vector(mode = "list", length = length(resultDirs))
for (index in 1:length(resultDirs)){
  setwd(resultDirs[[index]])
  
  #Read in layers
  clampingRegion <- raster(list.files(pattern = "ClampingMask.gri"))
  clampingRegion <- reclassify(clampingRegion, rcl = c(-Inf,0.1,1,0.9,Inf,0))
  ensProj <- raster(list.files(pattern = "EMwmeanByROC_mergedAlgo_mergedRun_mergedData.gri", 
                               full.names = T, recursive = T))
  
  #Read in TSS threshold
  setwd(paste0("USER_WORKING_DIRECTORY/data/Outputs/", "SDM_Redo/"))
  EnsScoreTable <- read.csv(file = threshFiles[index], header = T)
  threshVals[[index]] <- EnsScoreTable[1,11]/1000
  TSSvals[[index]] <- EnsScoreTable[1,10]
  
  #Process and save layers
  ensProjNoClamp <- (ensProj*clampingRegion)/1000
  writeRaster(ensProjNoClamp, filename = paste0("USER_WORKING_DIRECTORY/data/Outputs/", 
                                                "SDM_Redo/PostProcessedLayers/Continuous/ensProjNoClamp_", 
                                                sdm_bins[index], "k"), format = "ascii", overwrite = T)
  binaryEnsProj <- reclassify(ensProjNoClamp, rcl = c(-Inf, threshVals[[index]], 0, threshVals[[index]], Inf, 1))
  writeRaster(binaryEnsProj, filename = paste0("USER_WORKING_DIRECTORY/data/Outputs/", 
                                               "SDM_Redo/PostProcessedLayers/Binary/ensProjNoClampBinary_",
                                               sdm_bins[index], "k"), format = "ascii", overwrite = T)
}

# Ensemble projections for 5-11kbp ----
setwd(paste0(workingDirectory, "/PaleoClimateRasters/"))
dirs <- list.dirs()[-1]
dirs <- dirs[grepl(dirs, pattern = "X")]
sdm_sub <- sdm_bins[1:5]
occData <- vector("list", length = length(sdm_sub))

# Climate data for ensemble projections ----
projList <- 5:11
projLayers <- vector("list", length = length(projList))
for(index in 1:length(projLayers)){
  clList <- list.files(path = dirs[projList[[index]]+1], pattern = "bio", full.names = T)
  climate <- stack(clList)
  raster::crs(climate) <- raster::crs("+proj=longlat +datum=WGS84")
  projLayers[[index]] <- climate
}

setwd(paste0(workingDirectory, "/PaleoClimateRasters/"))
# Models
for (index in 1:length(sdm_sub)) {
  print(paste0("Now modeling ", sdm_dirs[[index]], "k."))
  print(Sys.time())
  clList <- list.files(path = sdm_dirs[index], pattern = "bio", full.names = T)
  climate <- stack(clList)
  raster::crs(climate) <- raster::crs("+proj=longlat +datum=WGS84")
  
  #Get the occurrences per time TimeStepIndex
  sliceOccs <- occurrences[occurrences$TimeStepIndex==sdm_sub[index],]
  
  # Format the data 
  occData <- BIOMOD_FormatingData(resp.var = rep(1, nrow(sliceOccs)), 
                                  expl.var = climate, 
                                  resp.xy = sliceOccs[, c("Longitude", "Latitude")], 
                                  resp.name = "Homo_Sapiens", 
                                  PA.nb.rep = 5, 
                                  PA.nb.absences = 5 * nrow(sliceOccs), 
                                  PA.strategy = "disk", 
                                  PA.dist.min = 150000, 
                                  PA.dist.max = 1250000)
  
  # Run SDM algorithms
  setwd(paste0("USER_WORKING_DIRECTORY/data/Outputs/", 
               "SDM_Redo/PostProcessedLayers/LateSlices"))
  models <- BIOMOD_Modeling(data = occData, 
                            models = c("MAXENT.Phillips", "GLM", "GAM", "ANN", "RF"), 
                            models.options = myBiomodOption, NbRunEval = 5, DataSplit = 80, 
                            VarImport = 3, models.eval.meth = c("TSS", "ROC"), 
                            modeling.id = paste0(occData@sp.name, sdm_sub[[index]], "k", sep=""),
                            compress = T, do.full.models = FALSE)
  ensemble_models <- BIOMOD_EnsembleModeling(modeling.output = models, em.by = "all", eval.metric = "ROC", 
                                             eval.metric.quality.threshold = 0.8, models.eval.meth = c("TSS", "ROC"), 
                                             prob.mean = T, prob.cv = TRUE, committee.averaging = FALSE, 
                                             prob.mean.weight = TRUE)
  ensemble_models_scores <- get_evaluations(ensemble_models)
  write.csv(ensemble_models_scores, paste0("EnsScores_", occData@sp.name, "_", sdm_sub[index], "k.csv", sep=""))
  
  # Do projections
  models_proj5k <- BIOMOD_Projection(modeling.output = models, new.env = projLayers[[1]],
                                   proj.name = paste0(sdm_sub[index], "k_proj_ens_5k"),
                                   build.clamping.mask = T,
                                   output.format = ".grd", do.stack = FALSE)
  
  ensemble_models_proj5k <- BIOMOD_EnsembleForecasting(EM.output = ensemble_models, projection.output = models_proj5k, 
                                                     proj.name = paste0(sdm_bins[index], "k_proj_ens_5k"),
                                                     output.format = ".grd", do.stack = FALSE)
  
  models_proj6k <- BIOMOD_Projection(modeling.output = models, new.env = projLayers[[2]],
                                     proj.name = paste0(sdm_sub[index], "k_proj_ens_6k"),
                                     build.clamping.mask = T,
                                     output.format = ".grd", do.stack = FALSE)
  
  ensemble_models_proj6k <- BIOMOD_EnsembleForecasting(EM.output = ensemble_models, projection.output = models_proj6k, 
                                                       proj.name = paste0(sdm_bins[index], "k_proj_ens_6k"),
                                                       output.format = ".grd", do.stack = FALSE)
  
  models_proj7k <- BIOMOD_Projection(modeling.output = models, new.env = projLayers[[3]],
                                     proj.name = paste0(sdm_sub[index], "k_proj_ens_7k"),
                                     build.clamping.mask = T,
                                     output.format = ".grd", do.stack = FALSE)
  
  ensemble_models_proj7k <- BIOMOD_EnsembleForecasting(EM.output = ensemble_models, projection.output = models_proj7k, 
                                                       proj.name = paste0(sdm_bins[index], "k_proj_ens_7k"),
                                                       output.format = ".grd", do.stack = FALSE)
  
  models_proj8k <- BIOMOD_Projection(modeling.output = models, new.env = projLayers[[4]],
                                     proj.name = paste0(sdm_sub[index], "k_proj_ens_8k"),
                                     build.clamping.mask = T,
                                     output.format = ".grd", do.stack = FALSE)
  
  ensemble_models_proj8k <- BIOMOD_EnsembleForecasting(EM.output = ensemble_models, projection.output = models_proj8k, 
                                                       proj.name = paste0(sdm_bins[index], "k_proj_ens_8k"),
                                                       output.format = ".grd", do.stack = FALSE)
  
  models_proj9k <- BIOMOD_Projection(modeling.output = models, new.env = projLayers[[5]],
                                     proj.name = paste0(sdm_sub[index], "k_proj_ens_9k"),
                                     build.clamping.mask = T,
                                     output.format = ".grd", do.stack = FALSE)
  
  ensemble_models_proj9k <- BIOMOD_EnsembleForecasting(EM.output = ensemble_models, projection.output = models_proj9k, 
                                                       proj.name = paste0(sdm_bins[index], "k_proj_ens_9k"),
                                                       output.format = ".grd", do.stack = FALSE)
  
  models_proj10k <- BIOMOD_Projection(modeling.output = models, new.env = projLayers[[6]],
                                     proj.name = paste0(sdm_sub[index], "k_proj_ens_10k"),
                                     build.clamping.mask = T,
                                     output.format = ".grd", do.stack = FALSE)
  
  ensemble_models_proj10k <- BIOMOD_EnsembleForecasting(EM.output = ensemble_models, projection.output = models_proj10k, 
                                                       proj.name = paste0(sdm_bins[index], "k_proj_ens_10k"),
                                                       output.format = ".grd", do.stack = FALSE)
  
  models_proj11k <- BIOMOD_Projection(modeling.output = models, new.env = projLayers[[7]],
                                      proj.name = paste0(sdm_sub[index], "k_proj_ens_11k"),
                                      build.clamping.mask = T,
                                      output.format = ".grd", do.stack = FALSE)
  
  ensemble_models_proj11k <- BIOMOD_EnsembleForecasting(EM.output = ensemble_models, projection.output = models_proj11k, 
                                                        proj.name = paste0(sdm_bins[index], "k_proj_ens_11k"),
                                                        output.format = ".grd", do.stack = FALSE)
  
  # Reset working directory to do it all over
  setwd(paste0(workingDirectory, "/PaleoClimateRasters/"))
}

# Get and process predictions for 5 to 11 kybp ensembles ----
setwd("USER_WORKING_DIRECTORY/data/Outputs/")
occurrences <- read.csv(paste0("USER_WORKING_DIRECTORY/data/Inputs/",
                               "HomoSapiensOccurrences_FinalForModeling_FINAL.csv"), 
                        header = T, stringsAsFactors = F)
resultDirs <- dir(path = paste0("USER_WORKING_DIRECTORY/data/Outputs/",
                                "SDM_Redo/PostProcessedLayers/LateSlices/Homo.Sapiens"), 
                  pattern = "k_proj", full.names = T)
clampingRemoved <- vector("list", length(resultDirs))

# Get and clamp layers
for (index in 1:length(resultDirs)){
  setwd(resultDirs[[index]])
  clampingRegion <- raster(list.files(pattern = "ClampingMask.gri"))
  clampingRegion <- reclassify(clampingRegion, rcl = c(-Inf,0.1,1,0.9,Inf,0))
  ensProj <- raster(list.files(pattern = "EMwmeanByROC_mergedAlgo_mergedRun_mergedData.gri", 
                               full.names = T, recursive = T))
  clampingRemoved[[index]] <- (ensProj*clampingRegion)/1000
}

# Combine ensembles
## Defines a function to rescale an ensemble model to input extremes
rasterRescale <- function(r, rmin, rmax){ 
  tmp <- ((r-cellStats(r, stat = "min"))*(rmax-rmin))/(cellStats(r, stat = "max")-cellStats(r, stat = "min"))
  return(tmp)
}

setwd("USER_WORKING_DIRECTORY/data/Outputs/")
ensembleMembers <- seq(12,16)
projectionSlice <- seq(5, 11)
threshFiles <- list.files("SDM_Redo/PostProcessedLayers/LateSlices", pattern = "EnsScores", recursive = F, full.names = T)
threshFiles <- lapply(threshFiles, read.csv)
threshVals <- unlist(lapply(threshFiles, function(x) x[1,11]/1000))
for(index in 1:length(projectionSlice)){
  ensMembsRast <- clampingRemoved[grep(paste0("_", projectionSlice[index], "k"), 
                                        resultDirs)]
  count <- 1
  binEnsMembs <- vector(mode = "list", length = 3)
  binEnsWeights <- vector(mode = "list", length = 3)
  sliceOccs <- occurrences[occurrences$TimeStepIndex==projectionSlice[index],]
  sliceOccs <- gridSample(sliceOccs[,c("Longitude", "Latitude")], ensMembsRast[[1]], n = 1)
  while(5 >= count){
    rcl <- c(-Inf, threshVals[count], 0, threshVals[count], Inf, 1)
    binEnsMembs[[count]] <- reclassify(ensMembsRast[[count]], rcl = rcl)
    ptOcc <- extract(binEnsMembs[[count]], sliceOccs)
    ptOcc <- ptOcc[complete.cases(ptOcc)]
    binEnsWeights[[count]] <- sum(ptOcc)/(length(ptOcc)*.95)
    ptOcc <- NULL
    count <- count + 1
  }
  # Binary ensemble assembly and writing
  binEnsCompiled <- binEnsMembs[[1]]*binEnsWeights[[1]] + 
    binEnsMembs[[2]]*binEnsWeights[[2]] + 
    binEnsMembs[[3]]*binEnsWeights[[3]] + 
    binEnsMembs[[4]]*binEnsWeights[[4]] + 
    binEnsMembs[[5]]*binEnsWeights[[5]]
  binThresh <- min(extract(binEnsCompiled, 
                           sliceOccs)[extract(binEnsCompiled, sliceOccs) != 0], na.rm = T)
  binEnsThresholded <- reclassify(binEnsCompiled, rcl = c(-Inf, binThresh, 0,
                                                          binThresh, Inf, 1))
  if(10 > projectionSlice[index]){
    writeRaster(binEnsThresholded, filename = paste0("USER_WORKING_DIRECTORY/data/Outputs/", 
                                                 "SDM_Redo/PostProcessedLayers/Binary/ensProjNoClampBinary_",
                                                 0, projectionSlice[index], "k"), 
                format = "ascii", overwrite = T)
  } else {
    writeRaster(binEnsThresholded, filename = paste0("USER_WORKING_DIRECTORY/data/Outputs/", 
                                                     "SDM_Redo/PostProcessedLayers/Binary/ensProjNoClampBinary_",
                                                     projectionSlice[index], "k"), 
                format = "ascii", overwrite = T)
  }
  # Continuous ensemble assembly and writing
  contEnsCompiled <- ensMembsRast[[1]]*binEnsWeights[[1]] + 
                      ensMembsRast[[2]]*binEnsWeights[[2]] + 
                      ensMembsRast[[3]]*binEnsWeights[[3]] + 
                      ensMembsRast[[4]]*binEnsWeights[[4]] + 
                      ensMembsRast[[5]]*binEnsWeights[[5]]
  maxFromAllMembers <- max(unlist(lapply(ensMembsRast, FUN = function(x) cellStats(x, max))))
  contEnsCompiledRescale <- rasterRescale(contEnsCompiled, rmin = 0, rmax = maxFromAllMembers)
  if(10 > projectionSlice[index]){
    writeRaster(contEnsCompiledRescale, filename = paste0("USER_WORKING_DIRECTORY/data/Outputs/", 
                                                          "SDM_Redo/PostProcessedLayers/Continuous/ensProjNoClamp_", 
                                                          0, projectionSlice[index], "k"), 
                format = "ascii", overwrite = T)
  } else {
    writeRaster(contEnsCompiledRescale, filename = paste0("USER_WORKING_DIRECTORY/data/Outputs/", 
                                                "SDM_Redo/PostProcessedLayers/Continuous/ensProjNoClamp_", 
                                                projectionSlice[index], "k"), format = "ascii", overwrite = T)
  }
}