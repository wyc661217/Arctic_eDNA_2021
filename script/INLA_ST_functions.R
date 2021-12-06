library("tidyr")
library("INLA")
library("splancs")
library("rgeos")
library("ggplot2")
library("reshape2")
library("dplyr")
library("gridExtra")
library("deldir")
library("viridis")
library(rnaturalearth)
library(ggspatial)
source("spde-book-functions.R")

# Find nearest value in an array
FindNearest <- function(array, value){
    dists <- abs((array - value))
    idx <- which(dists == min(dists))[1]
    nearest <- array[idx]
    return(nearest)
}

MinMaxNorm <- function(x,a=0,b=1){
    numerator <- (x - min(x,na.rm=TRUE))*(b-a)
    denominator <- max(x,na.rm=TRUE) - min(x,na.rm=TRUE)
    final <- a + numerator/denominator
    return(final)
}


LonLatToAzi <- function(test){
    colnames(test) <- c("Lon","Lat")
    test <- data.frame(test)
    coordinates(test) <- ~Lon+Lat
    proj4string(test) <- CRS("+init=epsg:4326") 
    CRS.new <- CRS("+proj=aeqd +lat_0=90 +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=km +no_defs")
    test2 <- data.frame(spTransform(test, CRS.new))
    test2 <- cbind(test2[,1],test2[,2])
    return(test2)
}


AziToLonLat <- function(m1){
    colnames(m1) <- c("X","Y")
    coordinates(m1) <- ~X+Y
    proj4string(m1) <- CRS("+proj=aeqd +lat_0=90 +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=km +no_defs")
CRS.new <- CRS("+init=epsg:4326")
    m2 <- data.frame(spTransform(m1, CRS.new))
    m2 <- cbind(m2[,1],m2[,2])
    colnames(m2) <- c("Lon","Lat")
    m2 <- data.frame(m2)
    return(m2)
}

# Run binomial model in INLA
RunInlaBin <- function(var=NA,sitetab=NA,xytab=NA,mesh.s=NA,mesh.t=NA,Ntrials=1,namescov=c(),normcov=0,PlotBox=NULL,project=FALSE){

    k <- length(mesh.t$loc)

    # Projection points
    if(is.null(PlotBox)){
        limx = range(xytab[, 1])
        limy = range(xytab[, 2])
    } else {
        limx = range(PlotBox[,1])
        limy = range(PlotBox[,2])
    }
    r0 <- diff(limx) / diff(limy)
    ngrid <- 20
    prj <- inla.mesh.projector(mesh.s, xlim = limx,
        ylim = limy, dims = c(ngrid * r0, ngrid))
    prjloc <- prj$lattice$loc
    prjloc <- prjloc[rep(seq(1,dim(prjloc)[1]),k),]
    prjtime <- as.vector(sapply(mesh.t$loc, function (x) rep(x,dim(prj$lattice$loc)[1])))
    A.pred = inla.spde.make.A(mesh=mesh.s, loc=prjloc,group=prjtime,group.mesh=mesh.t)
    npred <- attributes(A.pred)$Dim[1]

    print(var)
    sitetab <- data.frame(sitetab)
    #Z <- c(sitetab[,var])
    n <- dim(sitetab)[1]

    spde <- inla.spde2.pcmatern(mesh = mesh.s, 
        prior.range = c(3000, 0.5), # P( range < 3000) = 0.5
        prior.sigma = c(10, 0.01)) # P( sigma > 10) = 0.01

    iset <- inla.spde.make.index('i', n.spde = spde$n.spde, n.group = k)

    # Make projector matrix
    A <- inla.spde.make.A(mesh = mesh.s, loc = xytab, group = sitetab$Age, group.mesh = mesh.t) 

    # Load covariates
    #namescov <- c("annualtemp","annualprec","tempseaso","precseaso","human")
    tablist <- as.list(sitetab)
    covlist <- tablist[namescov]

    # Standardize covariates
    if(normcov == 1){
        covlist <- lapply(covlist,function(x){
                return(c(base::scale(x)))
        })
    } else if(normcov == 2){
    	covlist <- lapply(covlist,function(x){
       		if( all(unique(x) %in% c(0,1,NA)) ){ return(MinMaxNorm(x,-1,1))
        	} else{ return(c(base::scale(x))) }
    	})
    }

    # Data stack
    sdat <- inla.stack( 
    	 data = tablist[var], 
    	 A = list(A,1), 
    	 effects = list( iset, c( list(b0 = rep(1, n)), covlist ) ), 
    	 tag = "stdata")

    # Projection stack
    nalist <- list(NA)
    names(nalist) <- var
    stkgrid <- inla.stack(
         data = nalist,
         A = list(A.pred, 1),
         effects = list( iset, c( list(b0 = rep(1, npred ))) ),
         tag = 'prjgr')

    stk.all <- inla.stack(sdat, stkgrid)

    # Prior for temporal auto-regressive parameter
    h.spec <- list(theta = list(prior = 'pccor1', param = c(0, 0.9)))

    # Precision prior for default data model - NOT USED
    #prec.prior <- list(prior = 'pc.prec', param = c(1, 0.01))

    # Gaussian process model formula
    formula <- as.formula( paste(var, " ~ 0 + b0 + ",paste(namescov, collapse=" + ")," + f(i, model = spde, group = i.group,control.group = list(model = 'ar1', hyper = h.spec))",sep=""))

    cres <- list(return.marginals.predictor = FALSE,
    return.marginals.random = FALSE)

    print(Sys.time())

    # Binomial data model
    resdat <- inla(formula,
    	family="binomial", Ntrials=1, data = inla.stack.data(sdat), 
    	control.family = list(),
    	control.predictor = list(link=1,compute = TRUE,A = inla.stack.A(sdat)),
    	control.compute = list(cpo = TRUE, dic = TRUE, waic = TRUE),
    	control.results = cres)

    print(Sys.time())

    if(project == TRUE){
    # Binomial data model - projection
    resall <- inla(formula,
       family="binomial", Ntrials=1, data = inla.stack.data(stk.all),
       control.family = list(),
       control.predictor = list(link=1,compute = TRUE,A = inla.stack.A(stk.all)),
       control.compute = list(cpo = TRUE, dic = TRUE, waic = TRUE),
       control.mode = list(theta=resdat$mode$theta,restart=FALSE),
       control.results = cres)
    } else{
    resall <- resdat
    }
    print(Sys.time())

    return(resall)
}



# Run Hurdle-Gamma model in INLA
RunInlaHurdleGamma <- function(var=NA,sitetab=NA,xytab=NA,mesh.s=NA,mesh.t=NA,Ntrials=1,namescov=c(),normcov=0){

    print(var)
    sitetab <- data.frame(sitetab)
    n <- dim(sitetab)[1]

    spde <- inla.spde2.pcmatern(mesh = mesh.s, 
        prior.range = c(3000, 0.5), # P( range < 3000) = 0.5
        prior.sigma = c(10, 0.01)) # P( sigma > 10) = 0.01

    A <- inla.spde.make.A(mesh = mesh.s, loc = xytab, group = sitetab$Age, group.mesh = mesh.t)

    # Define spde fields
    field.z.idx <- inla.spde.make.index(name = 'x', 
        n.spde = spde$n.spde, n.group = k)
    field.zc.idx <- inla.spde.make.index(name = 'xc', 
        n.spde = spde$n.spde, n.group = k)
    field.y.idx <- inla.spde.make.index(name = 'u', 
        n.spde = spde$n.spde, n.group = k)

    z <- sitetab[,paste(var,"_P",sep="")]
    y <- sitetab[,paste(var,"_S",sep="")]

    # Covariates
    tablist <- as.list(sitetab)
    covlist <- tablist[namescov]

    # Standardize covariates
    if(normcov == 1){
        covlist <- lapply(covlist,function(x){
                return(c(base::scale(x)))
        })
    } else if(normcov == 2){
        covlist <- lapply(covlist,function(x){
                if( all(unique(x) %in% c(0,1,NA)) ){ return(MinMaxNorm(x,-1,1))
                } else{ return(c(base::scale(x))) }
        })
    }

    # Stack for presence data
    stk.z <- inla.stack(
        data = list(Y = cbind(as.vector(z), NA), link = 1), 
        A = list(A, 1),
	effects = list( field.z.idx, c( list(z.intercept = rep(1, n)), covlist  ) ),
        tag = 'zobs') 

    # Stack for score data
    stk.y <- inla.stack(
        data = list(Y = cbind(NA, as.vector(y)), link = 2), 
        A = list(A, 1),
	effects = list( field.y.idx, c( list(y.intercept = rep(1, n)), covlist ) ),
        tag = 'yobs')

    stk.all <- inla.stack(stk.z, stk.y)

    # Prior for gamma parameter
    pcgprior <- list(prior = 'pc.gamma', param = 1)
    cff <- list(list(), list(hyper = list(theta = pcgprior)))

    cinla <- list(strategy = 'adaptive', int.strategy = 'eb') 

    cres <- list(return.marginals.predictor = FALSE, 
    return.marginals.random = FALSE)

    # Prior for temporal auto-regressive parameter
    rhoprior <- list(theta = list(prior = 'pccor1',
        param = c(0, 0.9)))
    cg <- list(model = 'ar1', hyper = rhoprior)

    # Prior for scaling parameter for space-time effect
    bprior <- list(prior = 'gaussian', param = c(0,1))

    # Define formula for model
    #formula.joint <- Y ~ -1 + z.intercept + y.intercept + 
    #  f(x, model = spde, group = x.group, control.group = cg) + 
    #  f(xc, copy = "x", fixed = FALSE, group = xc.group,
    #    hyper = list(theta = bprior)) + 
    # f(u, model = spde, group = u.group, control.group = cg)  

    formula.joint <- as.formula(paste("Y ~ -1 + z.intercept + y.intercept + ", paste(namescov, collapse=" + "),
    	"+ f(x, model = spde, group = x.group, control.group = cg) + 
        f(u, model = spde, group = u.group, control.group = cg)",sep=""))

    # Run model
    res <- inla(formula.joint, family = c("binomial", "gamma"), 
        data = inla.stack.data(stk.all), control.family = cff, 
        control.predictor = list(A = inla.stack.A(stk.all),
            link = link), 
        control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE,
            config = TRUE),
        control.results = cres, control.inla = cinla, 
        control.mode = list(restart = TRUE)) 

    return(res)
}


# Obtain projected mesh for binomial probability field
ProjBinomMesh <- function(res,xytab,mesh.s,mesh.t,spde.s,PlotBox=NULL){

    if(is.null(PlotBox)){
	limx = range(xytab[, 1])
    	limy = range(xytab[, 2])
    } else {
      	limx = range(PlotBox[,1])
	limy = range(PlotBox[,2])
    }      	 
    
    r0 <- diff(limx) / diff(limy)
    prj <- inla.mesh.projector(mesh.s, xlim = limx,
        ylim = limy, dims = c(100 * r0, 100))
    #in.pr <- inout(prj$lattice$loc, xytab)
    m.prj <- lapply(1:mesh.t$n, function(j) {
        idx <- 1:spde.s$n.spde + (j - 1) * spde.s$n.spde
	#field <- res$summary.fitted.values[,"mean"]
	field <- res$summary.fixed[1,"mean"] + res$summary.ran$i$mean[idx]
        r <- inla.mesh.project(prj,field = field)
        #r[!in.pr] <- NA
        prob <- exp(r)/(1 + exp(r))
        return(prob)
    })
    return(list(m.prj,prj))
}



# Obtain projected mesh for binomial probability field
ProjBinomFieldVar <- function(res,xytab,mesh.s,mesh.t,spde.s,PlotBox=NULL){

    if(is.null(PlotBox)){
        limx = range(xytab[, 1])
        limy = range(xytab[, 2])
    } else {
        limx = range(PlotBox[,1])
        limy = range(PlotBox[,2])
    }

    r0 <- diff(limx) / diff(limy)
    prj <- inla.mesh.projector(mesh.s, xlim = limx,
        ylim = limy, dims = c(100 * r0, 100))
    #in.pr <- inout(prj$lattice$loc, xytab)
    m.prj <- lapply(1:mesh.t$n, function(j) {
        idx <- 1:spde.s$n.spde + (j - 1) * spde.s$n.spde
        field <- (res$summary.fixed[1,"sd"])^2 + (res$summary.ran$i$sd[idx])^2
        r <- inla.mesh.project(prj,field = field)
        #r[!in.pr] <- NA
        return(r)
    })
    return(list(m.prj,prj))
}




# Obtain projected mesh for Gamma Hurdle field
ProjGammaHurdleMesh <- function(res,xytab,mesh.s,mesh.t,spde.s,PlotBox=NULL){

    if(is.null(PlotBox)){
        limx = range(xytab[, 1])
        limy = range(xytab[, 2])
    } else {
        limx = range(PlotBox[,1])
        limy = range(PlotBox[,2])
    }

    r0 <- diff(limx) / diff(limy)
    prj <- inla.mesh.projector(mesh.s, xlim = limx,
        ylim = limy, dims = c(100 * r0, 100))
    z.prj <- lapply(1:mesh.t$n, function(j) {
    	  idx <- 1:spde.s$n.spde + (j - 1) * spde.s$n.spde
          zfield <- 0 +
	  res$summary.fixed["z.intercept","mean"] +
	  res$summary.ran$x$mean[idx]
          z <- inla.mesh.project(prj,field = zfield)
	  z <- exp(z)/(1 + exp(z))
          return(z)
    })
    y.prj <- lapply(1:mesh.t$n, function(j) {
    	  idx <- 1:spde.s$n.spde + (j - 1) * spde.s$n.spde
    	  yfield <- 0 +
          res$summary.fixed["y.intercept","mean"] +
          res$summary.ran$u$mean[idx]
    	  y <- inla.mesh.project(prj,field = yfield)
	  y <- exp(y)
    	  return(y)
    })
    zy.prj <- lapply(1:mesh.t$n, function(j) {
          idx <- 1:spde.s$n.spde + (j - 1) * spde.s$n.spde
	  zfield <- 0 +
          res$summary.fixed["z.intercept","mean"] +
          res$summary.ran$x$mean[idx]
	  z <- inla.mesh.project(prj,field = zfield)
          z <- exp(z)/(1 + exp(z))
          yfield <- 0 +
          res$summary.fixed["y.intercept","mean"] +
          res$summary.ran$u$mean[idx]
          y <- inla.mesh.project(prj,field = yfield)
          y <- exp(y)
          return(z*y)
    })
    return(list(z.prj,y.prj,zy.prj,prj))
}






PlotPolar <- function(m1,title=NA,minplot=NA,maxplot=NA,forcelegtit=NA){

# custom theme without axes and annotations
theme_polar <- function(){
   list(
      theme_bw(base_size=10),
      theme(
         panel.grid.major = element_blank(),
         panel.grid.minor = element_blank(),
         axis.text = element_blank(),
         axis.ticks = element_blank()),
      labs(x='',y=''))
}

dt <- rnaturalearth::ne_coastline()
clip_boundary <- sp::SpatialPolygons(
  list(sp::Polygons(
    list(sp::Polygon(
      data.frame(lon = c(-180, 180, 180, -180), lat = c(60, 60, 90, 90)))), ID = 1)
  ), proj4string = sp::CRS(sp::proj4string(dt)))

arctic <- raster::crop(dt, clip_boundary)
arctic <- sp::spTransform(arctic, sp::CRS("+proj=aeqd +lat_0=90 +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=km +no_defs"))
arctic_sf <- sf::st_as_sf(arctic)

if(is.na(forcelegtit)){
	newlegtit="Score"
	} else{
	newlegtit = forcelegtit
}

ggplot(data=m1,aes(x = Lon, y = Lat, z = Score,fill=Score)) + 
  #geom_point(data = m1, aes(x = Lon, y = Lat, fill = Score,colour=Score)) +
  geom_tile(aes(fill = Score)) +
  ggspatial::layer_spatial(data = arctic_sf,colour="white") +
  ggtitle(title) +
  scale_fill_viridis(limits = c(minplot,maxplot),name=newlegtit) +
  theme_polar()
}



CreateCovarPlots <- function(scoretab,names,modellist,forcemodel=NA){

covarplots <- lapply(names,function(response){
    modelvec <- scoretab[,response]
    idx <- which(modelvec == min(modelvec))
    score <- modelvec[idx] 
    model <- rownames(scoretab)[idx] 
    if(!is.na(forcemodel)){model <- forcemodel}
    print(c(response,score,model))
    fitted <- modellist[[model]][[response]]
    positive <- which(fitted$summary.fixed[,"0.025quant"] > 0 & fitted$summary.fixed[,"0.975quant"] > 0 )
    positivenames <- rownames(fitted$summary.fixed[positive,])
    negative <- which(fitted$summary.fixed[,"0.025quant"] < 0 & fitted$summary.fixed[,"0.975quant"] < 0 )
    negativenames <- rownames(fitted$summary.fixed[negative,])
    zero <- which(fitted$summary.fixed[,"0.025quant"] < 0 & fitted$summary.fixed[,"0.975quant"] > 0 )
    zeronames <- rownames(fitted$summary.fixed[zero,])
    allcols <- rep("zero",dim(fitted$summary.fixed)[1])
    allcols[positive] <- "positive"
    allcols[negative] <- "negative"
    cols <- c("zero" = "dark grey", "positive" = "blue", "negative" = "red")
    validall <- cbind(rownames(fitted$summary.fixed),allcols,fitted$summary.fixed)
    colnames(validall) <- c("covariate","colour","effect","sd","quantA","quantB","quantC","mode","kld")
    # Remove intercept
    validall <- validall[rownames(validall)!="b0",]
    validall <- validall[rownames(validall)!="z.intercept",]
    validall <- validall[rownames(validall)!="y.intercept",]
    p1 <- ggplot(data=validall) + 
        geom_point(aes(y=covariate,x=effect,colour=colour)) +
        geom_errorbarh(aes(y=covariate,xmin =quantA,xmax=quantC,colour=colour)) +
        scale_color_manual(values=cols,guide=FALSE) +
        geom_vline(xintercept = 0) +
        ggtitle(response)
    return(p1)
    #barplot(validall[,"mean"],cex.names=0.5,horiz=TRUE,las=1,main=response)
    #return(validnames)
})

return(covarplots)

}



CreateCovarTabs <- function(scoretab,names,modellist,forcemodel=NA){

covartabs <- lapply(names,function(response){
    modelvec <- scoretab[,response]
    idx <- which(modelvec == min(modelvec))
    score <- modelvec[idx]
    model <- rownames(scoretab)[idx]
    if(!is.na(forcemodel)){model <- forcemodel}
    print(c(response,score,model))
    fitted <- modellist[[model]][[response]]
    positive <- which(fitted$summary.fixed[,"0.025quant"] > 0 & fitted$summary.fixed[,"0.975quant"] > 0 )
    positivenames <- rownames(fitted$summary.fixed[positive,])
    negative <- which(fitted$summary.fixed[,"0.025quant"] < 0 & fitted$summary.fixed[,"0.975quant"] < 0 )
    negativenames <- rownames(fitted$summary.fixed[negative,])
    zero <- which(fitted$summary.fixed[,"0.025quant"] < 0 & fitted$summary.fixed[,"0.975quant"] > 0 )
    zeronames <- rownames(fitted$summary.fixed[zero,])
    validall <- fitted$summary.fixed
    # Remove intercept
    validall <- validall[rownames(validall)!="b0",]
    validall <- validall[rownames(validall)!="z.intercept",]
    validall <- validall[rownames(validall)!="y.intercept",]
    validall <- data.frame(animal=row.names(validall),round(validall,3))
    #validall <- cbind(rownames(fitted$summary.fixed),fitted$summary.fixed)
    colnames(validall) <- c("covariate","mean","sd","2.5% quantile","50% quantile","97.5% quantile","mode","kld")

    return(validall)
})

return(covartabs)

}
