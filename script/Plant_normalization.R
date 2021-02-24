


##########################
#This R script applies 6 different normalization methods for eliminating the effect of sequencing depth on the taxa abundance
#YAntonio Fernandez Guerra, antonio.fernandez-guerra@sund.ku.dk



#data read in 
genus = read.delim("plant_genus.csv",sep=",",stringsAsFactors = F) #a csv file with each row as a sample and each column as a taxon

#format transformation
rownames(genus) = genus$taxa
genus = genus[,-1]
genus = as.data.frame(t(genus))
genus = as.matrix(genus)


#1 Standardize abundances to the median sequencing depth
physeq_median = genus
total <- median(colSums(genus))

for (i in 1:dim(physeq_median)[2]) {
  x = physeq_median[,i]
  physeq_median[,i] = round(total * (x / sum(x)))
}
vsn::meanSdPlot(physeq_median)


#2 Standardize abundances to a common-scale
physeq_com <- genus
depths <- rowSums(physeq_com)
mindepth  <- min(depths)
physeq_com <- physeq_com * (mindepth/depths)


#3 Transform data with CSS
physeq_css = t(genus)

MGS <- newMRexperiment(
  counts = (physeq_css)
)
MGS <- cumNorm(MGS, p = cumNormStat(MGS))

physeq_css = as.matrix(MRcounts(
  MGS,
  norm = T,
  log = T,
  sl = median(unlist(normFactors(MGS)))))

myround <- function(x) { trunc(x + 0.5) }
physeq_css = myround(physeq_css)
physeq_css = t(physeq_css)


#4-6 DESeq2 count
#convert into DESeq2 format
physeq_ds = t(genus)
period = data.frame(sample=colnames(physeq_ds), 
                    period=rep(1,(dim(physeq_ds)[2])))
physeq_ds <- DESeqDataSetFromMatrix(countData = physeq_ds, colData = period, ~ 1, tidy = FALSE)

# We need a custom geometric function to deal with 0s
gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm = na.rm) / length(x))
}
geoMeans <- apply(counts(physeq_ds), 1, gm_mean)

# We want to use poscounts as its designed for metagenomic studies
diagdds <- estimateSizeFactors(physeq_ds, type = "poscounts", locfunc = genefilter::shorth, geoMeans = geoMeans)
diagdds <- estimateDispersions(diagdds, fitType = 'local')
diagvst <- varianceStabilizingTransformation(diagdds, blind = FALSE)
diagvst <- assay(diagvst)

# DESeq2 when doing the VST maps 0 to a certain value instead of to -Infinity.
# Let's find this value
get_xi <- function(X){
  ncounts <- counts(X, normalized = TRUE)
  sf <- sizeFactors(X)
  xg <- sinh(seq(asinh(0), asinh(max(ncounts)), length.out = 1000))[-1]
  xim <- mean(1/sf)
  baseVarsAtGrid <- dispersionFunction(X)(xg) * xg^2 +
    xim * xg
  integrand <- 1/sqrt(baseVarsAtGrid)
  splf <- splinefun(asinh((xg[-1] + xg[-length(xg)])/2),
                    cumsum((xg[-1] - xg[-length(xg)]) * (integrand[-1] +
                                                           integrand[-length(integrand)])/2))
  h1 <- quantile(rowMeans(ncounts), 0.95)
  h2 <- quantile(rowMeans(ncounts), 0.999)
  eta <- (log2(h2) - log2(h1))/(splf(asinh(h2)) - splf(asinh(h1)))
  xi <- log2(h1) - eta * splf(asinh(h1))
  xi
}

vst_xi <- get_xi(diagdds)

# We will do log2 transformation
diagvst_log2 <- normTransform(diagdds, f = log2, pc = 0.00001)
diagvst_log2 <- assay(diagvst_log2)

# As shown before 0s are mapped to a value to avoid -Inf
# First count number of 0s in our transformed table, no 0s are left in the table
sum(diagvst == 0)

# Let's substract this value to the transformed values.
diagvst <- diagvst - vst_xi
# We will have negative values
sum(diagvst == 0)
sum(diagvst < 0)

# As we will have negative values (probably correspond to “less than one count” after rescaling)
# we will convert them to 0 to avoid problems with distance/dissimilarity measures that have problems with
# negative values
diagvst[diagvst < 0] <- 0

# The same for the log2 transformed data
diagvst_log2[diagvst_log2 < 0] <- 0

# Those are the variance stabilized counts
physeq_vst <- t(myround(diagvst))
# Those are log2 counts
physeq_log2 <- t(myround(diagvst_log2))
# This just divides each sample by the size factor.
physeq_norm <- t(counts(diagdds, normalized = TRUE))

