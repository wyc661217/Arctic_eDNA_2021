R
source("INLA_ST_func.R")

# Build grid table
annualtemp <- tibble(read.table("climate/meanAnnualTemperature.csv",sep=",",header=TRUE))
annualtemp <- gather(annualtemp,Age,annualtemp,-Site)
annualtemp$Age <- as.numeric(gsub("X","",annualtemp$Age))*1000
annualprec <- tibble(read.table("climate/annualPrecipitation.csv",sep=",",header=TRUE))
annualprec <- gather(annualprec,Age,annualprec,-Site)
annualprec$Age <- as.numeric(gsub("X","",annualprec$Age))*1000
tempseaso <- tibble(read.table("climate/temperatureSeasonality.csv",sep=",",header=TRUE))
tempseaso <- gather(tempseaso,Age,tempseaso,-Site)
tempseaso$Age <- as.numeric(gsub("X","",tempseaso$Age))*1000
precseaso <- tibble(read.table("climate/precipitationSeasonality.csv",sep=",",header=TRUE))
precseaso <- gather(precseaso,Age,precseaso,-Site)
precseaso$Age <- as.numeric(gsub("X","",precseaso$Age))*1000
human <- tibble(read.table("human_distribution/humanPresenceAbsence.csv",sep=",",header=TRUE))
human <- gather(human,Age,human,-Site)
human$Age <- as.numeric(gsub("kbp","",gsub("X","",human$Age)))*1000
gridtab <- full_join(annualtemp,annualprec)
gridtab <- full_join(gridtab,tempseaso)
gridtab <- full_join(gridtab,precseaso)
gridtab <- full_join(gridtab,human)
allgridtimes <- sort(unique(gridtab$Age))


# Build site table
#plant <- t(read.table("plant/life_form_richness_2.csv",sep=",",header=TRUE,row.names=1))
plant <- read.table("plant/nmds_k3.csv",sep=",",header=TRUE,row.names=1)
plant <- data.frame(row.names(plant),plant)
colnames(plant)[1] <- "sample"

animal <- read.table("animal/distribution_matrix_withHerbivore.csv",sep=",",header=TRUE)
metadata <- read.table("sample_metadata/sample_metadata_v17.csv",sep=",",header=TRUE)
metadata$Age <- as.numeric(metadata$Age)
metadata <- metadata[which(metadata$Age >= 0),]
metadata <- tibble(metadata)
sitetab <- left_join(metadata,animal,by = c("Lab_ID" = "sample") )
sitetab <- left_join(sitetab,plant,by = c("Lab_ID" = "sample") )

# Plant, Animal and Environment names
allplantnames <- colnames(plant)[-1]
allanimalnames <- colnames(animal)[-1]
allenvnames <- colnames(gridtab)[-c(1,2)]

# Add covariates to site table
covtab <- tibble(do.call(rbind,apply(sitetab,1,function(vec){
    names(vec) <- colnames(sitetab)
    Lab_ID <- vec["Lab_ID"]
    Site <- vec["Site_abbreviation"]
    Age <- as.numeric(vec["Age"])
    NearestAge <- FindNearest(allgridtimes,Age)
    covars <- cbind(Lab_ID,gridtab[which(gridtab$Site == Site & gridtab$Age == NearestAge),])
    return(covars)
})))
sitetab <- left_join(sitetab,covtab,by = c("Lab_ID"))
colnames(sitetab)[which(colnames(sitetab) == "Age.x")] <- "Age"
colnames(sitetab)[which(colnames(sitetab) == "Age.y")] <- "NearestAge"
colnames(sitetab)[which(colnames(sitetab) == "Rigion")] <- "Region"


# Filter for animals with total P > 0.1
ratios <- sapply(allanimalnames,function(response){
    vec <- sitetab[,response][!is.na(sitetab[,response])]
    ratio <- sum(vec)/length(vec)
    return(ratio)
})
goodanimalnames <- allanimalnames[which(ratios > 0.1)]
goodanimalnames <- goodanimalnames[which(goodanimalnames != "Herbivore")]

# Filter out algae
goodplantnames <- allplantnames[allplantnames != "Algae"]



# Prepare spatial and temporal meshes
k <- 11
alltidx <- seq(1,k)
reducedtidx <- seq(2,k,2)
tknots <- seq(min(sitetab$Age), max(sitetab$Age), length = k)
reducedtknots <- tknots[seq(2,k,2)]
mesh.t <- inla.mesh.1d(tknots)

xytab <- cbind(sitetab$Longitude,sitetab$Latitude)
xytab <- LonLatToAzi(xytab)
n <- dim(sitetab)[1]

bound <- inla.nonconvex.hull(as.matrix(xytab))
mesh.s <- inla.mesh.2d(bound = bound,max.edge = c(200, 400),max.n=200,cutoff = 500)
plot(mesh.s)
spde.s <- inla.spde2.matern(mesh.s)



# Model with no covariates
reslist <- list()
for(response in goodanimalnames){

    # Run INLA - Binomial model with number of trials = 1 (no covariates)
    reslist[[response]] <- RunInlaBin(response,sitetab,xytab,mesh.s,mesh.t,1,namescov=c(),normcov=0)
}


# Define lists to store results
resclimnormlist <- list()
reshumnormlist <- list()
resclimhumnormlist <- list()
resaninormlist <- list()
resplanormlist <- list()
resclimhumplanormlist <- list()
resaniplanormlist <- list()
resclimhumaninormlist <- list()
rescovnormlist <- list()

for(response in goodanimalnames){    

    # Run INLA - Binomial model with number of trials = 1 with climate covariates (normalized)
    namesclim <- allenvnames[which(allenvnames != "human")]
    resclimnormlist[[response]] <- RunInlaBin(response,sitetab,xytab,mesh.s,mesh.t,1,namescov=namesclim,normcov=2)

    # Run INLA - Binomial model with number of trials = 1 with human covariates (normalized)
    reshumnormlist[[response]] <- RunInlaBin(response,sitetab,xytab,mesh.s,mesh.t,1,namescov=c("human"),normcov=2)

    # Run INLA - Binomial model with number of trials = 1 with human and climate covariates (normalized)
    resclimhumnormlist[[response]] <- RunInlaBin(response,sitetab,xytab,mesh.s,mesh.t,1,namescov=allenvnames,normcov=2)

    # Run INLA - Binomial model with number of trials = 1 with plant covariates
    resplanormlist[[response]] <- RunInlaBin(response,sitetab,xytab,mesh.s,mesh.t,1,namescov=goodplantnames,normcov=2)

    # Run INLA - Binomial model with number of trials = 1 with human, climate and plant covariates (normalized)
    namesclimhumpla <- c(allenvnames,goodplantnames)
    namesclimhumpla <- namesclimhumpla[namesclimhumpla != response]
    resclimhumplanormlist[[response]] <- RunInlaBin(response,sitetab,xytab,mesh.s,mesh.t,1,namescov=namesclimhumpla,normcov=2)

    # Run INLA - Binomial model with number of trials = 1 with animal covariates (normalized)
    animalcov <- goodanimalnames[goodanimalnames != response]
    resaninormlist[[response]] <- RunInlaBin(response,sitetab,xytab,mesh.s,mesh.t,1,namescov=animalcov,normcov=2)

    # Run INLA - Binomial model with number of trials = 1 with animal and plant covariates (normalized)
    aniplacov <- c(goodplantnames,goodanimalnames)
    aniplacov <- aniplacov[aniplacov != response]
    resaniplanormlist[[response]] <- RunInlaBin(response,sitetab,xytab,mesh.s,mesh.t,1,namescov=aniplacov,normcov=2)

    # Run INLA - Binomial model with number of trials = 1 with human, climate and animal covariates (normalized)
    namesclimhumani <- c(allenvnames,goodanimalnames)
    namesclimhumani <- namesclimhumani[namesclimhumani != response]
    resclimhumaninormlist[[response]] <- RunInlaBin(response,sitetab,xytab,mesh.s,mesh.t,1,namescov=namesclimhumani,normcov=2)

    # Run INLA - Binomial model with number of trials = 1 with all covariates (normalized)
    allcov <- c(allenvnames,goodplantnames,goodanimalnames)
    allcov <- allcov[allcov != response]
    rescovnormlist[[response]] <- RunInlaBin(response,sitetab,xytab,mesh.s,mesh.t,1,namescov=allcov,normcov=2)

}

animallists <- list()
animallists[["None"]] <- reslist
animallists[["Climate"]] <- resclimnormlist
animallists[["Humans"]] <- reshumnormlist
animallists[["Climate+Humans"]] <- resclimhumnormlist
animallists[["Plants"]] <- resplanormlist
animallists[["Climate+Humans+Plants"]] <- resclimhumplanormlist
animallists[["Animals"]] <- resaninormlist
animallists[["Animals+Plants"]] <- resaniplanormlist
animallists[["Climate+Humans+Animals"]] <- resclimhumaninormlist
animallists[["Climate+Humans+Animals+Plants"]] <- rescovnormlist


# CPO
modelcands <- names(animallists)
cpotab <- sapply(goodanimalnames, function(response){
    return(sapply( modelcands, function(modelname){ -sum(log(animallists[[modelname]][[response]]$cpo$cpo),na.rm=TRUE) }))
})
rownames(cpotab) <- modelcands
cpotabmelted <- melt(cpotab)
colnames(cpotabmelted) <- c("Model","Animal","CPO")
cpoplot <- ggplot(cpotabmelted) +
geom_point(aes(x=Model,y=CPO,group=Animal,colour=Animal)) +
geom_line(aes(x=Model,y=CPO,group=Animal,colour=Animal)) +
theme(axis.text.x = element_text(angle = 45, hjust = 1), 
axis.title.x=element_blank(),
text = element_text(size = 10)
)

# waic
modelcands <- names(animallists)
waictab <- sapply(goodanimalnames, function(response){
    return(sapply( modelcands, function(modelname){animallists[[modelname]][[response]]$waic$waic}))
})
rownames(waictab) <- modelcands
waictabmelted <- melt(waictab)
colnames(waictabmelted) <- c("Model","Animal","WAIC")
waicplot <- ggplot(waictabmelted) +
geom_point(aes(x=Model,y=WAIC,group=Animal,colour=Animal)) +
geom_line(aes(x=Model,y=WAIC,group=Animal,colour=Animal)) +
theme(axis.text.x = element_text(angle = 45, hjust = 1),
axis.title.x=element_blank(),
text = element_text(size = 10)
)

cpoplotlet <- arrangeGrob(cpoplot, top = textGrob("A", x = unit(0, "npc"), y   = unit(1, "npc"), just=c("left","top"),
         gp=gpar(col="black", fontsize=22, fontfamily="Arial")))
waicplotlet <- arrangeGrob(waicplot, top = textGrob("B", x = unit(0, "npc"), y   = unit(1, "npc"), just=c("left","top"),
         gp=gpar(col="black", fontsize=22, fontfamily="Arial")))
grid.arrange(cpoplotlet,waicplotlet, ncol=2, vp=viewport(width=0.98, height=0.98))

write.table(data.frame(Model=row.names(waictab),round(waictab,2)),file="tables/waictab_animals.tsv",col.names=TRUE,row.names=FALSE,sep="\t",quote=FALSE)
write.table(data.frame(Model=row.names(cpotab),round(cpotab,2)),file="tables/cpotab_animals.tsv",col.names=TRUE,row.names=FALSE,sep="\t",quote=FALSE)

covarplots <- CreateCovarPlots(waictab,goodanimalnames,animallists)
grid.arrange(covarplots[[1]],covarplots[[2]],covarplots[[3]],covarplots[[4]],covarplots[[5]],covarplots[[6]], vp=viewport(width=0.98, height=0.98))

covarplots <- CreateCovarPlots(waictab,goodanimalnames,animallists,forcemodel="Climate+Humans+Plants")
grid.arrange(covarplots[[1]],covarplots[[2]],covarplots[[3]],covarplots[[4]],covarplots[[5]],covarplots[[6]], vp=viewport(width=0.98, height=0.98))

covartabs <- CreateCovarTabs(waictab,goodanimalnames,animallists)
for(i in seq(1,length(goodanimalnames))){
    animal <- goodanimalnames[i]
    tab <- covartabs[[i]]
    write.table(tab,file=paste("tables/",animal,"_besttab.tsv",sep=""),col.names=TRUE,row.names=FALSE,sep="\t",quote=FALSE)
}

covartabs <- CreateCovarTabs(waictab,goodanimalnames,animallists,forcemodel="Climate+Humans+Plants")
for(i in seq(1,length(goodanimalnames))){
    animal <- goodanimalnames[i]
    tab <- covartabs[[i]]
    write.table(tab,file=paste("tables/",animal,"_climhumplatab.tsv",sep=""),col.names=TRUE,row.names=FALSE,sep="\t",quote=FALSE)
}


# Herbivore

response <- "Herbivore"    

# Run INLA - Binomial model with number of trials = 1 (no covariates)
reslist[[response]] <- RunInlaBin(response,sitetab,xytab,mesh.s,mesh.t,1,namescov=c(),normcov=0)

# Run INLA - Binomial model with number of trials = 1 with climate covariates (normalized)
namesclim <- allenvnames[which(allenvnames != "human")]
resclimnormlist[[response]] <- RunInlaBin(response,sitetab,xytab,mesh.s,mesh.t,1,namescov=namesclim,normcov=2)

# Run INLA - Binomial model with number of trials = 1 with human covariates (normalized)
reshumnormlist[[response]] <- RunInlaBin(response,sitetab,xytab,mesh.s,mesh.t,1,namescov=c("human"),normcov=2)

# Run INLA - Binomial model with number of trials = 1 with human and climate covariates (normalized)
resclimhumnormlist[[response]] <- RunInlaBin(response,sitetab,xytab,mesh.s,mesh.t,1,namescov=allenvnames,normcov=2)

# Run INLA - Binomial model with number of trials = 1 with plant covariates
resplanormlist[[response]] <- RunInlaBin(response,sitetab,xytab,mesh.s,mesh.t,1,namescov=goodplantnames,normcov=2)

# Run INLA - Binomial model with number of trials = 1 with human, climate and plant covariates (normalized)
namesclimhumpla <- c(allenvnames,goodplantnames)
namesclimhumpla <- namesclimhumpla[namesclimhumpla != response]
resclimhumplanormlist[[response]] <- RunInlaBin(response,sitetab,xytab,mesh.s,mesh.t,1,namescov=namesclimhumpla,normcov=2)


animallists <- list()
animallists[["None"]] <- reslist
animallists[["Climate"]] <- resclimnormlist
animallists[["Humans"]] <- reshumnormlist
animallists[["Climate+Humans"]] <- resclimhumnormlist
animallists[["Plants"]] <- resplanormlist
animallists[["Climate+Humans+Plants"]] <- resclimhumplanormlist
animallists[["Animals"]] <- resaninormlist
animallists[["Animals+Plants"]] <- resaniplanormlist
animallists[["Climate+Humans+Animals"]] <- resclimhumaninormlist
animallists[["Climate+Humans+Animals+Plants"]] <- rescovnormlist


# CPO
modelcands <- names(animallists)[seq(1,6)]
cpotab <- sapply(c("Herbivore"), function(response){
    return(sapply( modelcands, function(modelname){ -sum(log(animallists[[modelname]][[response]]$cpo$cpo),na.rm=TRUE) }))
})
rownames(cpotab) <- modelcands
cpotabmelted <- melt(cpotab)
colnames(cpotabmelted) <- c("Model","Animal","CPO")
cpoplot <- ggplot(cpotabmelted) +
geom_point(aes(x=Model,y=CPO,group=Animal,colour=Animal)) +
geom_line(aes(x=Model,y=CPO,group=Animal,colour=Animal)) +
theme(axis.text.x = element_text(angle = 45, hjust = 1), 
axis.title.x=element_blank(),
text = element_text(size = 10)
)

# WAIC
modelcands <- names(animallists)[seq(1,6)]
waictab <- sapply(c("Herbivore"), function(response){
    return(sapply( modelcands, function(modelname){animallists[[modelname]][[response]]$waic$waic}))
})
rownames(waictab) <- modelcands
waictabmelted <- melt(waictab)
colnames(waictabmelted) <- c("Model","Animal","WAIC")
waicplot <- ggplot(waictabmelted) +
geom_point(aes(x=Model,y=WAIC,group=Animal,colour=Animal)) +
geom_line(aes(x=Model,y=WAIC,group=Animal,colour=Animal)) +
theme(axis.text.x = element_text(angle = 45, hjust = 1),
axis.title.x=element_blank(),
text = element_text(size = 10)
)

cpoplotlet <- arrangeGrob(cpoplot, top = textGrob("A", x = unit(0, "npc"), y   = unit(1, "npc"), just=c("left","top"),
         gp=gpar(col="black", fontsize=22, fontfamily="Arial")))
waicplotlet <- arrangeGrob(waicplot, top = textGrob("B", x = unit(0, "npc"), y   = unit(1, "npc"), just=c("left","top"),
         gp=gpar(col="black", fontsize=22, fontfamily="Arial")))
grid.arrange(cpoplotlet,waicplotlet, ncol=2, vp=viewport(width=0.98, height=0.98))


write.table(data.frame(Model=row.names(waictab),round(waictab,2)),file="tables/waictab_herbivores.tsv",col.names=TRUE,row.names=FALSE,sep="\t",quote=FALSE)
write.table(data.frame(Model=row.names(cpotab),round(cpotab,2)),file="tables/cpotab_herbivores.tsv",col.names=TRUE,row.names=FALSE,sep="\t",quote=FALSE)

covarplots <- CreateCovarPlots(waictab,c("Herbivore"),animallists)
grid.arrange(covarplots[[1]], vp=viewport(width=0.98, height=0.98))


covartabs <- CreateCovarTabs(waictab,c("Herbivore"),animallists)
animal <- "Herbivore"
tab <- covartabs[[1]]
write.table(tab,file=paste("tables/",animal,"_besttab.tsv",sep=""),col.names=TRUE,row.names=FALSE,sep="\t",quote=FALSE)



# Project animal Binomial probability field
for(response in goodanimalnames){
res <- reslist[[response]]
PlotBox <- LonLatToAzi(rbind(c(135,50),c(45,50),c(-45,50),c(-135,50)))
projl <- ProjBinomMesh(res,xytab,mesh.s,mesh.t,spde.s,PlotBox=PlotBox)
minplot <- 0
maxplot <- 1
mplots <- lapply(reducedtidx,function(i){
m1 <- melt(projl[[1]][[i]])
m1 <- cbind(projl[[2]]$lattice$loc,m1[,3])
m1 <- data.frame(m1)
colnames(m1) <- c("Lon","Lat","Score")
title <- paste("t = ",round(tknots[i])," years",sep="")
p1 <- PlotPolar(m1,title=title,minplot=minplot,maxplot=maxplot,forcelegtit="Prob")
return(p1)
})
#grid.arrange(mplots[[1]], mplots[[2]], mplots[[3]], mplots[[4]], mplots[[5]], mplots[[6]], nrow = 2,top=textGrob(response,gp=gpar(fontsize=20,font=2)))
#png(paste(response,".png",sep=""),width=3000,height=1500,res=200)
#grid.arrange(mplots[[1]], mplots[[2]], mplots[[3]], mplots[[4]], mplots[[5]], mplots[[6]], nrow = 2,top=textGrob(response,gp=gpar(fontsize=20,font=2)))
#dev.off()
png(paste(response,".png",sep=""),width=4000,height=750,res=300)
grid.arrange(mplots[[1]], mplots[[2]], mplots[[3]], mplots[[4]], mplots[[5]], nrow = 1,top=textGrob(response,gp=gpar(fontsize=20,font=2)))
dev.off()
}



# Project animal field variance
for(response in goodanimalnames){
res <- reslist[[response]]
PlotBox <- LonLatToAzi(rbind(c(135,50),c(45,50),c(-45,50),c(-135,50)))
projl <- ProjBinomFieldVar(res,xytab,mesh.s,mesh.t,spde.s,PlotBox=PlotBox)
minplot <- min(unlist(projl[[1]]))
maxplot <- max(unlist(projl[[1]]))
mplots <- lapply(reducedtidx,function(i){
m1 <- melt(projl[[1]][[i]])
m1 <- cbind(projl[[2]]$lattice$loc,m1[,3])
m1 <- data.frame(m1)
colnames(m1) <- c("Lon","Lat","Score")
title <- paste("t = ",round(tknots[i])," years",sep="")
p1 <- PlotPolar(m1,title=title,minplot=minplot,maxplot=maxplot,forcelegtit="Var")
return(p1)
})
#grid.arrange(mplots[[1]], mplots[[2]], mplots[[3]], mplots[[4]], mplots[[5]], mplots[[6]], nrow = 2,top=textGrob(response,gp=gpar(fontsize=20,font=2)))
#png(paste(response,".png",sep=""),width=3000,height=1500,res=200)
#grid.arrange(mplots[[1]], mplots[[2]], mplots[[3]], mplots[[4]], mplots[[5]], mplots[[6]], nrow = 2,top=textGrob(response,gp=gpar(fontsize=20,font=2)))
#dev.off()
png(paste(response,"_Var.png",sep=""),width=4000,height=750,res=300)
grid.arrange(mplots[[1]], mplots[[2]], mplots[[3]], mplots[[4]], mplots[[5]], nrow = 1,top=textGrob(response,gp=gpar(fontsize=20,font=2)))
dev.off()
}






Z <- sitetab$Graminoid
tab <- cbind(xytab,sitetab$Age,Z)
colnames(tab) <- c("X","Y","TIME","VALS")
tab[,1] <- jitter(tab[,1],amount=0.2)
tab[,2] <- jitter(tab[,2],amount=0.2)
tab <- as.data.frame(tab)
table <- tibble(tab)
toplot <- table
#toplot <- table[which(table$VALS == 1),]
ggplot(toplot, aes(x=X, y=Y, colour=TIME)) + geom_point() + scale_colour_gradient(low="blue", high="red")






# PLANTS - USE HURDLE GAMMA MODEL

# Record presence-absence for plants
plantpres <- t(apply(sitetab[,allplantnames] > 0,1,as.numeric))
colnames(plantpres) <- paste(allplantnames,"_P",sep="")
plantscore <- ifelse(plantpres == 1, unlist(sitetab[,allplantnames]), NA)
colnames(plantscore) <- paste(allplantnames,"_S",sep="")
sitetab <- cbind(sitetab,plantpres,plantscore)



# Model with no covariates
for(response in goodplantnames){

    # Run INLA - Hurdle Gamma model (no covariates)
    reslist[[response]] <- RunInlaHurdleGamma(response,sitetab,xytab,mesh.s,mesh.t,1,namescov=c(),normcov=0)
}



for(response in goodplantnames){    

    # Run INLA - Gamma-Hurdle  model with number of trials = 1 with climate covariates (normalized)
    namesclim <- allenvnames[which(allenvnames != "human")]
    resclimnormlist[[response]] <- RunInlaHurdleGamma(response,sitetab,xytab,mesh.s,mesh.t,1,namescov=namesclim,normcov=2)

    # Run INLA - Gamma-Hurdle  model with number of trials = 1 with human covariates (normalized)
    reshumnormlist[[response]] <- RunInlaHurdleGamma(response,sitetab,xytab,mesh.s,mesh.t,1,namescov=c("human"),normcov=2)

    # Run INLA - Binomial model with number of trials = 1 with human and climate covariates (normalized)
    resclimhumnormlist[[response]] <- RunInlaHurdleGamma(response,sitetab,xytab,mesh.s,mesh.t,1,namescov=allenvnames,normcov=2)

    # Run INLA - Gamma-Hurdle  model with number of trials = 1 with plant covariates
    plantcov <- goodplantnames[goodplantnames != response]
    resplanormlist[[response]] <- RunInlaHurdleGamma(response,sitetab,xytab,mesh.s,mesh.t,1,namescov=plantcov,normcov=2)

    # Run INLA - Gamma-Hurdle  model with number of trials = 1 with human, climate and plant covariates (normalized)
    namesclimhumpla <- c(allenvnames,goodplantnames)
    namesclimhumpla <- namesclimhumpla[namesclimhumpla != response]
    resclimhumplanormlist[[response]] <- RunInlaHurdleGamma(response,sitetab,xytab,mesh.s,mesh.t,1,namescov=namesclimhumpla,normcov=2)

    # Run INLA - Gamma-Hurdle  model with number of trials = 1 with animal covariates (normalized)
    resaninormlist[[response]] <- RunInlaHurdleGamma(response,sitetab,xytab,mesh.s,mesh.t,1,namescov=goodanimalnames,normcov=2)

    # Run INLA - Gamma-Hurdle  model with number of trials = 1 with animal and plant covariates (normalized)
    aniplacov <- c(goodplantnames,goodanimalnames)
    aniplacov <- aniplacov[aniplacov != response]
    resaniplanormlist[[response]] <- RunInlaHurdleGamma(response,sitetab,xytab,mesh.s,mesh.t,1,namescov=aniplacov,normcov=2)

    # Run INLA - Gamma-Hurdle  model with number of trials = 1 with human, climate and animal covariates (normalized)
    namesclimhumani <- c(allenvnames,goodanimalnames)
    namesclimhumani <- namesclimhumani[namesclimhumani != response]
    resclimhumaninormlist[[response]] <- RunInlaHurdleGamma(response,sitetab,xytab,mesh.s,mesh.t,1,namescov=namesclimhumani,normcov=2)

    # Run INLA - Gamma-Hurdle model with number of trials = 1 with all covariates (normalized)
    allcov <- c(allenvnames,goodplantnames,goodanimalnames)
    allcov <- allcov[allcov != response]
    rescovnormlist[[response]] <- RunInlaHurdleGamma(response,sitetab,xytab,mesh.s,mesh.t,1,namescov=allcov,normcov=2)

}



plantlists <- list()
plantlists[["None"]] <- reslist
plantlists[["Climate"]] <- resclimnormlist
plantlists[["Humans"]] <- reshumnormlist
plantlists[["Climate+Humans"]] <- resclimhumnormlist
plantlists[["Plants"]] <- resplanormlist
plantlists[["Climate+Humans+Plants"]] <- resclimhumplanormlist
plantlists[["Animals"]] <- resaninormlist
plantlists[["Animals+Plants"]] <- resaniplanormlist
plantlists[["Climate+Humans+Animals"]] <- resclimhumaninormlist
plantlists[["Climate+Humans+Animals+Plants"]] <- rescovnormlist


# CPO
cpotab <- sapply(goodplantnames, function(response){
    print(response)
    vec <- sapply( names(plantlists), function(modelname){ -sum(log(plantlists[[modelname]][[response]]$cpo$cpo),na.rm=TRUE) })
    return(vec)
})
rownames(cpotab) <- names(plantlists)
cpotabmelted <- melt(cpotab)
colnames(cpotabmelted) <- c("Model","Animal","CPO")
cpoplot <- ggplot(cpotabmelted) +
geom_point(aes(x=Model,y=CPO,group=Animal,colour=Animal)) +
geom_line(aes(x=Model,y=CPO,group=Animal,colour=Animal)) +
theme(axis.text.x = element_text(angle = 45, hjust = 1), 
axis.title.x=element_blank(),
text = element_text(size = 10)
)

# waic
waictab <- sapply(goodplantnames, function(response){
    print(response)
    vec <- sapply( names(plantlists), function(modelname){plantlists[[modelname]][[response]]$waic$waic})
    return(vec)
})
rownames(waictab) <- names(plantlists)
waictabmelted <- melt(waictab)
colnames(waictabmelted) <- c("Model","Animal","WAIC")
waicplot <- ggplot(waictabmelted) +
geom_point(aes(x=Model,y=WAIC,group=Animal,colour=Animal)) +
geom_line(aes(x=Model,y=WAIC,group=Animal,colour=Animal)) +
theme(axis.text.x = element_text(angle = 45, hjust = 1),
axis.title.x=element_blank(),
text = element_text(size = 10)
)

cpoplotlet <- arrangeGrob(cpoplot, top = textGrob("A", x = unit(0, "npc"), y   = unit(1, "npc"), just=c("left","top"),
         gp=gpar(col="black", fontsize=22, fontfamily="Arial")))
waicplotlet <- arrangeGrob(waicplot, top = textGrob("B", x = unit(0, "npc"), y   = unit(1, "npc"), just=c("left","top"),
         gp=gpar(col="black", fontsize=22, fontfamily="Arial")))
grid.arrange(cpoplotlet,waicplotlet, ncol=2, vp=viewport(width=0.98, height=0.98))


covarplots <- CreateCovarPlots(waictab,goodplantnames,plantlists)
grid.arrange(covarplots[[1]],covarplots[[2]],covarplots[[3]],covarplots[[4]], vp=viewport(width=0.98, height=0.98))

covarplots <- CreateCovarPlots(waictab,goodplantnames,plantlists,forcemodel="Climate+Humans+Animals")
grid.arrange(covarplots[[1]],covarplots[[2]],covarplots[[3]],covarplots[[4]], vp=viewport(width=0.98, height=0.98))






# Project plant Gamma-Hurdle field (Bernoulli component)
for(response in goodplantnames){
res <- reslist[[response]]
PlotBox <- LonLatToAzi(rbind(c(135,50),c(45,50),c(-45,50),c(-135,50)))
projl <- ProjGammaHurdleMesh(res,xytab,mesh.s,mesh.t,spde.s,PlotBox=PlotBox)
minplot <- 0
maxplot <- 1
mplots <- lapply(reducedtidx,function(i){
m1 <- melt(projl[[1]][[i]])
m1 <- cbind(projl[[3]]$lattice$loc,m1[,3])
m1 <- data.frame(m1)
colnames(m1) <- c("Lon","Lat","Score")
title <- paste("t = ",round(tknots[i])," years",sep="")
p1 <- PlotPolar(m1,title=title,minplot=minplot,maxplot=maxplot,forcelegtit="Prob")
return(p1)
})
png(paste(response,"_Prob.png",sep=""),width=4000,height=750,res=300)
grid.arrange(mplots[[1]], mplots[[2]], mplots[[3]], mplots[[4]], mplots[[5]], nrow = 1,top=textGrob(response,gp=gpar(fontsize=20,font=2)))
dev.off()
}


# Project plant Gamma-Hurdle field (Gamma component)
for(response in goodplantnames){
res <- reslist[[response]]
PlotBox <- LonLatToAzi(rbind(c(135,50),c(45,50),c(-45,50),c(-135,50)))
projl <- ProjGammaHurdleMesh(res,xytab,mesh.s,mesh.t,spde.s,PlotBox=PlotBox)
minplot <- min(unlist(projl[[3]]))
maxplot <- max(unlist(projl[[3]]))
mplots <- lapply(reducedtidx,function(i){
m1 <- melt(projl[[3]][[i]])
m1 <- cbind(projl[[4]]$lattice$loc,m1[,3])
m1 <- data.frame(m1)
colnames(m1) <- c("Lon","Lat","Score")
title <- paste("t = ",round(tknots[i])," years",sep="")
p1 <- PlotPolar(m1,title=title,minplot=minplot,maxplot=maxplot,forcelegtit="Score")
return(p1)
})
png(paste(response,"_Div.png",sep=""),width=4000,height=750,res=300)
grid.arrange(mplots[[1]], mplots[[2]], mplots[[3]], mplots[[4]], mplots[[5]], nrow = 1,top=textGrob(response,gp=gpar(fontsize=20,font=2)))
dev.off()
}

