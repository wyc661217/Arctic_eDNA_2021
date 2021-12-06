library(fields)
library(viridisLite)
library(RColorBrewer)

book.rMatern <- function(n, coords, sigma=1, range, kappa = sqrt(8*nu)/range, variance = sigma^2, nu=1) {
    m <- as.matrix(dist(coords))
    m <- exp((1-nu)*log(2) + nu*log(kappa*m)-
             lgamma(nu))*besselK(m*kappa, nu)
    diag(m) <- 1
    return(drop(crossprod(chol(variance*m),
                          matrix(rnorm(nrow(coords)*n), ncol=n))))
}

book.rspde <- function(coords, sigma=1, range, variance=sigma^2, alpha=2, kappa = sqrt(8*(alpha-1))/range, n=1, mesh, 
                  verbose=FALSE, seed, return.attributes=FALSE) {
    t0 <- Sys.time()
    theta <- c(-0.5*log(4*pi*variance*kappa^2), log(kappa))
    if (verbose) cat('theta =', theta, '\n')
    if (missing(mesh)) {
        mesh.pars <- c(0.5, 1, 0.1, 0.5, 1)*sqrt(alpha-ncol(coords)/2)/kappa 
        if (verbose) cat('mesh.pars =', mesh.pars, '\n')
        attributes <- list(
            mesh=inla.mesh.2d(,
                coords[chull(coords), ], max.edge=mesh.pars[1:2], 
                cutoff=mesh.pars[3], offset=mesh.pars[4:5]))
        if (verbose) cat('n.mesh =', attributes$mesh$n, '\n')
    }
    else attributes <- list(mesh=mesh)
    attributes$spde <- inla.spde2.matern(attributes$mesh, alpha=alpha)
    attributes$Q <- inla.spde2.precision(attributes$spde, theta=theta)
    attributes$A <- inla.mesh.project(mesh=attributes$mesh, loc=coords)$A
    if (n==1) 
        result <- drop(attributes$A%*%inla.qsample(
            Q=attributes$Q,
            constr=attributes$spde$f$extraconstr))
    t1 <- Sys.time() 
    result <- inla.qsample(n, attributes$Q, 
                           seed=ifelse(missing(seed), 0, seed), 
                           constr=attributes$spde$f$extraconstr) 
    if (nrow(result)<nrow(attributes$A)) {
        result <- rbind(result, matrix(
            NA, nrow(attributes$A)-nrow(result), ncol(result)))
        dimnames(result)[[1]] <- paste('x', 1:nrow(result), sep='')
        for (j in 1:ncol(result)) 
            result[, j] <- drop(attributes$A%*%
                                result[1:ncol(attributes$A),j])
    }
    else {
        for (j in 1:ncol(result)) 
            result[1:nrow(attributes$A), j] <-
                drop(attributes$A%*%result[,j]) 
        result <- result[1:nrow(attributes$A), ]
    }
    t2 <- Sys.time()
    attributes$cpu <- c(prep=t1-t0, sample=t2-t1, total=t2-t0)
    if (return.attributes) 
        attributes(result) <- c(attributes(result), attributes)
    return(drop(result))
}

book.mesh.dual <- function(mesh) {
    if (mesh$manifold=='R2') {
        ce <- t(sapply(1:nrow(mesh$graph$tv), function(i)
            colMeans(mesh$loc[mesh$graph$tv[i, ], 1:2])))
        library(parallel)
        pls <- mclapply(1:mesh$n, function(i) {
            p <- unique(Reduce('rbind', lapply(1:3, function(k) {
                j <- which(mesh$graph$tv[,k]==i)
                if (length(j)>0) 
                    return(rbind(ce[j, , drop=FALSE],
                                 cbind(mesh$loc[mesh$graph$tv[j, k], 1] +
                                       mesh$loc[mesh$graph$tv[j, c(2:3,1)[k]], 1], 
                                       mesh$loc[mesh$graph$tv[j, k], 2] +
                                       mesh$loc[mesh$graph$tv[j, c(2:3,1)[k]], 2])/2))
                else return(ce[j, , drop=FALSE])
            })))
            j1 <- which(mesh$segm$bnd$idx[,1]==i)
            j2 <- which(mesh$segm$bnd$idx[,2]==i)
            if ((length(j1)>0) | (length(j2)>0)) {
                p <- unique(rbind(mesh$loc[i, 1:2], p,
                                  mesh$loc[mesh$segm$bnd$idx[j1, 1], 1:2]/2 +
                                  mesh$loc[mesh$segm$bnd$idx[j1, 2], 1:2]/2, 
                                  mesh$loc[mesh$segm$bnd$idx[j2, 1], 1:2]/2 +
                                  mesh$loc[mesh$segm$bnd$idx[j2, 2], 1:2]/2))
                yy <- p[,2]-mean(p[,2])/2-mesh$loc[i, 2]/2
                xx <- p[,1]-mean(p[,1])/2-mesh$loc[i, 1]/2
            }
            else {
                yy <- p[,2]-mesh$loc[i, 2]
                xx <- p[,1]-mesh$loc[i, 1]
            }
            Polygon(p[order(atan2(yy,xx)), ])
        })
        return(SpatialPolygons(lapply(1:mesh$n, function(i)
            Polygons(list(pls[[i]]), i))))
    }
    else stop("It only works for R2!")
}

genColor <- function(n, type=c('red', 'green', 'blue'), u=NULL) {
    cbp <- list(
        red = list(c(255, 254, 252, 252, 251, 239, 203, 165, 103), 
                   c(245, 224, 187, 146, 106, 59, 24, 15, 0), 
                   c(240, 210, 161, 114, 74, 44, 29, 21, 13)), 
        green = list(c(247, 229, 199, 161, 116, 65, 35, 0, 0), 
                     c(252, 245, 233, 217, 196, 171, 139, 109, 68), 
                     c(245, 224, 192, 155, 118, 93, 69, 44, 27)), 
        blue = list(c(247, 222, 198, 158, 107, 66, 33, 8, 8), 
                    c(251, 235, 219, 202, 174, 146, 113, 81, 48), 
                    c(255, 247, 239, 225, 214, 198, 181, 156, 107)))
    if (n<2) stop("Works for 'n>2'!")
    if (is.null(u))
        u <- 0:(n-1)/(n-1)
    u0 <- 0:8/8
    i <- findInterval(u, u0, TRUE)
    k <- pmatch(match.arg(type), c('red', 'green', 'blue'))
    w1 <- 8*(u0[i+1]-u)/255; w2 <- 8*(u-u0[i])/255
    rgb(cbp[[k]][[1]][i]*w1 + cbp[[k]][[1]][i+1]*w2, 
        cbp[[k]][[2]][i]*w1 + cbp[[k]][[2]][i+1]*w2, 
        cbp[[k]][[3]][i]*w1 + cbp[[k]][[3]][i+1]*w2)
}

plot.dgTMatrix <- function(x, y, ...) {
    cl <- match.call()
    if (is.null(cl$digits))
        digits <- 2
    z <- sort(unique(round(x@x, digits)))
    nz <- length(z)
    n1 <- sum(z<0)
    n2 <- sum(z>0)
    if (is.null(cl$colors)) 
        if (any(c(n1,n2)==0)) 
            colors <- gray(0.9*(1-(z-min(z))/diff(range(z))))
        else
            colors <- c(genColor(n1, 'red', z[z<0]/min(z)),
                        rep('white', nz-n1-n2),
                        genColor(n2, 'blue', z[z>0]/max(z)))
    z.breaks <- c(z[1]-diff(z[1:2])/2,
                  z[-nz]/2 + z[-1]/2,
                  z[nz]+diff(z[nz-1:0])/2)
    x@x <- round(x@x, digits)
    image(x, at=z.breaks, col.regions=colors, ...)
}

book.plot.field <- function(field, mesh, projector, xlim, ylim, 
			    dims=c(300,300), poly, asp = 1, 
			    axes = FALSE, xlab = '', ylab = '', 
			    col = book.color.c(), ...){
  ## you can supply field as a matrix vector or like a named list with 'x', 'y' and 'z' as for image
  ## when field is a vector, it will project it using projector, assuming projector will create a matrix 
  ## when mesh is supplied and projector not, projector will be created and used to project field
  if (missing(mesh)) {
    if (missing(projector)) {
      if (missing(xlim) | missing(ylim)) {
        image.plot(field, asp = asp, axes = axes, 
                   xlab = xlab, ylab = ylab, col = col, ...)
      } else {
        image.plot(field, xlim = xlim, ylim = ylim, asp = asp, 
                   axes = axes, xlab = xlab, ylab = ylab, col = col, ...)
      }
    } else {
      if (missing(xlim)) xlim <- range(projector$x)
      if (missing(ylim)) ylim <- range(projector$y)
      field.proj <- inla.mesh.project(projector, field)
      image.plot(x = projector$x, y = projector$y, z = field.proj, 
                 asp=asp, axes=axes, xlab = xlab, ylab = ylab, 
                 col=col, xlim=xlim, ylim=ylim, ...)
    }
  } else {
    if (missing(xlim)) xlim <- range(mesh$loc[,1])
    if (missing(ylim)) ylim <- range(mesh$loc[,2])
    projector <- inla.mesh.projector(mesh, xlim = xlim,
                                     ylim = ylim, dims=dims)
    field.proj <- inla.mesh.project(projector, field)
    image.plot(x = projector$x, y = projector$y, z = field.proj, 
               asp=asp, axes=axes, xlab = xlab, ylab = ylab, col=col, ...)
  }
  if (!missing(poly)) 
      plot(poly, add = TRUE, col = 'grey')
}

## Functions for barrier models

## Find the correlation of precision Q (defined on mesh) at location 
book.spatial.correlation <- function(Q, location, mesh) {
  ## The marginal standard deviations
  sd <- sqrt(diag(inla.qinv(Q)))

  ## Create a fake A matrix, to extract the closest mesh node index
  A.tmp <- inla.spde.make.A(mesh = mesh,
    loc = matrix(c(location[1], location[2]), 1, 2))
  id.node = which.max(A.tmp[1, ])

  ## Solve a matrix system to find just one column of the covariance matrix
  Inode <- rep(0, dim(Q)[1])
  Inode[id.node] <- 1
  covar.column <- solve(Q, Inode)
  corr <- drop(matrix(covar.column)) / (sd * sd[id.node])
  return(corr)
}

## Continuous and discrete colour scales for the book
# n=8; plot(1:n, col=brewer.pal(n = n, name = "Paired"))

# Continuous
book.color.c = function(n = 201) {
  return(viridis(n))
}

# Continuous (alternative)
book.color.c2 = function(n = 201) {
  return(magma(n))
}

# Discrete from a continuous
book.color.dc = function(n = 11) {
  return(viridis(n))
}

# Discrete (cannot be interpolated)
book.color.d = function(n=4) {
  return(brewer.pal(n = n, name = "Paired"))
}
