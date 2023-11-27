### all scripting from Petri Mikael Kemppainen

#' loads a data set from a snapp and coords file
#'
#'
#' This is a function that reads the old data format. Presented for backwards
#' compatibility
#' @param snapp.file a file in snapp data format
#' @param coords.file a file containing coordinates and other info 
#'      about the sample
#' @param n.snp the maximum number of snp to be read. if -1, all
#'      snps are read
#' @param ploidy the ploidy of the organism, 1 for haploids, 2 for diploid
#' @param ... further arguments passed to load.coord.file
#' @return an object of type origin.data with entries `genotypes`,
#'      containing the genetic data and entry `coord` containing location
#'      data
#' @example examples/example_1.r
#' @export
#' 

# raw.data <- load.data.snapp(snapp.file="GT.snapp",
#                             coords.file="coord.csv",
#                             n.snp=-1,
#                             sep=',', ploidy=2)

library(snpStats)
library(survival)
library(Matrix)

load.data.snapp <- function(snapp.file, coords.file, n.snp=-1, 
                            ploidy=2, ...){
  f <- load.snapp.file(snapp.file, n.snp=n.snp)
  coords <- load.coord.file(coords.file, sep=",")
  raw.data <- check.missing(f, coords)
  raw.data <- set.outgroups(raw.data, ploidy)
  class(raw.data) <- 'origin.data'
  return(raw.data)
}


#' Reads data from SNAPP file
#' 
#' Deprecated, use plink input format instead
#' @example examples/example_1.r
load.snapp.file <- function(snp.file=snapp.file,
                            n.snp=-1){
  data <- read.table(snp.file, sep =",", header=F, strings=F,
                     na.strings="?", nrow=n.snp, row.names=1)
  ids <- rownames(data)
  # data[1:10, 1:10]
  raw.data <- list()
  data <- apply(as.matrix(data),2,as.numeric)
  rownames(data) <- ids
  #apply(as.matrix(data)[1:10, 1:10],2,as.numeric)
  raw.data$genotypes <- as(data, "SnpMatrix")
  rownames(raw.data$genotypes) <- ids
  raw.data$fam <- rownames(data)
  return(raw.data)
}


#' Loads a coordinate data set
#'
#' Loads a file that specifies location data and other sample specific information
#' A
#'
#' @param file: the name of the file which the data are to be read from.
#' @param ... Additional arguments passed to read.table
#' @return A data.frame with columns sample, longitude, latitude, region, 
#'outgroup and country 
#' @example examples/example_1.r
#' @export

# file <-  "coord.csv"
load.coord.file <- function(file, ...){
  coords <- read.table(file, header=T, strings=F, ...)
  names(coords)[1] <- 'id'
  coords[,1] <- as.character(coords[,1])
  return(coords)
}

#' Cleans genotype and coordinate files for consistency
#'
#' 
#'
#' This is a generic function: methods can be defined for it directly
#' or via the \code{\link{Summary}} group generic. For this to work properly,
#' the arguments \code{...} should be unnamed, and dispatch is on the
#' first argument.
#'
#' @param snp.data output from load.plink.file
#' @param coord.data output from load.coord.file
#' @return An object with entries 'genotypes' and 'coords' containing the
#'    input data sets
#'
#' @example examples/example_1.r
#' @export

# snp.data <- f
# coord.data = coords
check.missing <- function(snp.data, coord.data){
  
  ids <- rownames(snp.data$genotypes)
  to.remove <- !unlist(lapply(ids,
                              function(x)x%in%coord.data$id))
  res.data <- list()
  res.data$coords <- coord.data
  res.data$genotypes <- snp.data$genotypes[!to.remove,]
  ids <- ids[!to.remove]
  
  data.ordering <- res.data$coords$id
  
  res.data$genotypes <- res.data$genotypes[data.ordering,]
  
  return(res.data)
}


#' Sets columns of outgroup individuals
#'
#' As psi requires the derived state of an allele to be known, this function
#' provides an easy way to polarize SNP with outgroup individuals. From the
#' outgroup columns in the coord file, outgroup individuals are removed
#' and SNPs swithced if necessary
#'
#' @param data output from check missing
#' @param ploidy the number of copies per individual
#' @return An object with entries 'genotypes' and 'coords' containing the
#'    input data sets
#'
#' @example examples/example_1.r
#' @export
set.outgroups <- function(data, ploidy=2
){
  if(is.null(data$coords$outgroup)){ 
    return(data)
  }
  
  outgroup.columns <- which(as.logical(data$coords$outgroup))
  n.outgroups <- length(outgroup.columns)
  if( n.outgroups > 0 ){
    data$outgroup <- data$genotypes[outgroup.columns,]
    data$outgroup <- as(data$outgroup, 'numeric')
    data$genotypes <- data$genotypes[-1* outgroup.columns,]
    
    n.alleles <- colSums(!is.na(data$outgroup))
    derived <- colSums(data$outgroup==ploidy, na.rm=T) == n.alleles
    anc <- colSums(data$outgroup==0, na.rm=T)== n.alleles
    no.outgrp <- is.na(derived+anc) | !(anc | derived)  
    
    data$genotypes <- switch.alleles(data$genotypes,derived)
    
    data$genotypes <- data$genotypes[,!no.outgrp ]
    data$coords <- data$coords[-outgroup.columns,]
  }
  
  return(data)
}




#' Sets up the popultion structure data sets
#'
#' Groups up raw data into a data structure that 
#' 
#' @param raw.data Output from set.outgroups or check.data
#' @param ploidy the ploidy of organisms
#' @return A list with entries
#' - data : an n x p matrix with derived allele counts
#' - coords : an p x 3 matrix of coords
#' - n : number of populations
#' @example examples/example_1.r
#' @export
make.pop <- function(raw.data, ploidy=2){
  if(is.null(raw.data$coord$pop)){
    raw.data<- make.pops(raw.data)
  }
  pop <- list()
  pop$data <- make.pop.data(raw.data)[,-1]
  pop$ss <- make.pop.ss(raw.data, ploidy)[,-1]
  pop$coords <- make.pop.coords(raw.data)
  pop$coords$hets <- rowMeans(pop$data/pop$ss, na.rm=T)
  pop$n <- nrow(pop$coords)
  class(pop) <- 'population'
  return(pop)
}
#hets <- (coords$hets - min(coords$hets))/(max(coords$hets) - min(coords$hets))

make.pop.data <- function(raw.data){
  data <- as(raw.data$genotypes, 'numeric')
  aggregate(data, by=list(raw.data$coords$pop), FUN=sum, na.rm=T)
}
make.pop.ss <- function(raw.data, ploidy=2){
  data <- as(raw.data$genotypes, 'numeric')
  aggregate(data, by=list(raw.data$coords$pop), 
            FUN=function(x)ploidy*sum(!is.na(x)))
}

make.pop.coords <- function(raw.data){
  lat <- aggregate(raw.data$coords$latitude, 
                   by=list(raw.data$coords$pop), 
                   FUN=mean)
  names(lat) <- c("pop", "latitude")
  long <- aggregate(raw.data$coords$longitude, 
                    by=list(raw.data$coords$pop), 
                    FUN=mean)
  names(long) <- c("pop", "longitude")
  pop.coords <- merge(lat,long)
  region.id <- match(pop.coords$pop, raw.data$coords$pop)
  pop.coords$region <- raw.data$coords$region[region.id]
  return(pop.coords)
}


#' Puts individuals into populations for further analysis
#'
#' There are two main modes of grouping individuals into populations: 
#' The first one is based on the location, all individuals with the 
#' same sample location are assigned to the same population.
#' the other one is based on a column `population` in 
#'
#'
#' @param raw.data Output from set.outgroups or check.data
#' @param mode One of 'coord' or 'custom'. If 'custom', populations are
#'   generated from the population column in the coords file. If 'coord',
#'   Individuals are grouped according to their coordinates
#' @return A list with each entry being a population
make.pops <- function(raw.data, mode='coord'){
  if(mode == 'coord'){
    s <- c('latitude', 'longitude')
    rowSums <- function(x)200*x[,1]+x[,2]
    raw.data$coords$pop <-  match(rowSums(raw.data$coords[,s]), 
                                  rowSums(raw.data$coords[,s]))
  } 
  return(raw.data)
}

# pop <- pop_sim
# bin= TRUE

flatten <- function(x) {
  if (!inherits(x, "list")) return(list(x))
  else return(unlist(c(lapply(x, flatten)), recursive = FALSE))
}


get.all.psi.mc.bin <- function(pop, n=2,cores=10){
  #this function calculates psi for columns i,j, both resampled down to
  # n samples
  
  n.pops <- pop$n
  
  ## get counts for all pops
  counts <- mclapply(1:n.pops,function(i){
    ni <- unlist(pop$ss[i,])
    fi <- unlist(pop$data[i,])
    list(ni,  fi)
  },mc.cores = cores)
  
  mat = matrix( 0, nrow=n.pops, ncol=n.pops )
  pvals = matrix( 0, nrow=n.pops, ncol=n.pops)
  cnts = matrix( 0, nrow=n.pops, ncol=n.pops)
  
  ## get psi
  # i <- 1
  # j <- 57
  out <- do.call(rbind, flatten(mclapply(1:(n.pops-1),function(i){
    
    lapply((i+1):n.pops, function(j){
      #print(c(i,j))
      c(i,j,get.psi.bin(counts[[i]][[1]], counts[[j]][[1]], counts[[i]][[2]], counts[[j]][[2]]))
      
    })
    
  },mc.cores=cores)))
  
  
  mat[out[,1:2]] <- out[,3]
  mat[lower.tri(mat)] <- -t(mat)[lower.tri(mat)]
  
  
  pvals[out[,1:2]] <- out[,4]
  pvals[lower.tri(pvals)] <- t(pvals)[lower.tri(pvals)]
  diag(pvals) <- NA
  cnts[out[,1:2]] <- out[,5]
  cnts[lower.tri(cnts)] <- t(cnts)[lower.tri(cnts)]
  
  return(list(psi=mat,pvals=pvals,cnts=cnts))
}

#list(psi=NA,pvals=NA,cnts=NA)
get.psi.bin <- function (ni, nj, fi, fj){
  n = 2
  fn <- cbind( fi, ni, fj, nj)
  
  tbl <- table(as.data.frame(fn))
  tbl <- as.data.frame( tbl )
  tbl <- tbl[tbl$Freq > 0,]
  
  tbl$fi <- as.integer(as.character( tbl$fi ))
  tbl$fj <- as.integer(as.character( tbl$fj ))
  tbl$ni <- as.integer(as.character( tbl$ni ))
  tbl$nj <- as.integer(as.character( tbl$nj ))
  
  to.exclude <- tbl$fi == 0 | tbl$fj == 0 | tbl$ni < n | tbl$nj < n
  
  tbl <- tbl[! to.exclude, ]
  
  if(nrow(tbl)==0){ return(NaN)}
  
  poly.mat <- matrix(0,nrow=n+1, ncol=n+1)
  poly.mat[2:(n+1),2:(n+1)] <- 1
  poly.mat[n+1,n+1] <- 0
  
  #psi.mat is the contribution to psi for each entry
  psi.mat <- outer(0:n,0:n,FUN=function(x,y)(y-x))
  psi.mat[1,] <- 0
  psi.mat[,1] <- 0
  
  f.contribution <- function(row, b=n){
    a <- 0:b
    f1 <- row[1]
    n1 <- row[2]
    f2 <- row[3]
    n2 <- row[4]
    cnt <- row[5]
    q1 <- choose(b, a) * choose(n1-b, f1-a)/choose(n1,f1)
    q2 <- choose(b, a) * choose(n2-b, f2-a)/choose(n2,f2)
    return( cnt * outer(q2, q1) )
  }
  
  
  resampled.mat <- matrix(rowSums(apply(tbl, 1,f.contribution)),nrow=n+1)
  sum(resampled.mat * psi.mat) / sum(resampled.mat * poly.mat)
  
  cnts <- c(as.integer(resampled.mat[2,3]),as.integer(resampled.mat[3,2]))
  
  bin <- ifelse(sum(cnts)>0, binom.test(cnts)$p.value, 1)
  
  return( c(sum(resampled.mat * psi.mat) / sum(resampled.mat * poly.mat) ,bin,sum(cnts)))
}

#' preps data for tdoa analysis
#'
#' makes a data structure of format xi, yi, xj, yj, psi
#' 
#' @param pop.coords coordination data with latitude and longitude
#' @param all.psi matrix with psi values
#' @param region region or set of regions to run analysis on
#' @param countries countries to restrict analysis to
#' @param xlen number of x points
#' @param ylen number of y points
#'
prep.tdoa.data <- function(coords, psi){
  locs <- coords[,c('longitude', 'latitude')]
  n.locs <- nrow(locs)
  
  tdoa.data <- c()
  for(i in 1:n.locs){
    for(j in (i+1):n.locs){
      if( i>=j | j >n.locs) break
      tdoa.data <- rbind( tdoa.data, c(locs[i,], locs[j,], psi[i,j]))
    }
  }
  
  tdoa.data <- matrix(unlist(tdoa.data), ncol=5)
  
  return(tdoa.data)
}

div_from_sfs <- function(pops, GTs,seq_len=1e6){
  rbindlist(lapply(unique(pops),function(pop){
    sfs <- sfs1D(pops,pop,GTs,seq_len)
    n=length(sfs)
    
    # Calculate n-1th harmonic number to use for Whatterson's theta estimation
    numChrom=n-1
    harmonicNumber = 0
    for (j in 1:(numChrom - 1)) {
      harmonicNumber = harmonicNumber + 1.0/j
    }
    
    # creates weights for each class. This assumes the spectrum is unfolded. 
    # The "weights" object stores allele frequencies for each spectrum category
    
    p = seq(0,numChrom)/numChrom
    
    # now get the terms to calculate Tajima's D
    
    a1 = sum(1/seq(1,numChrom))
    a2 = sum(1/seq(1,numChrom)^2)
    b1 = (numChrom+1)/(3*(numChrom-1))
    b2 = 2*(numChrom^2 + numChrom + 3)/(9*numChrom * (numChrom-1))
    c1 = b1 - 1/a1
    c2 = b2 - (numChrom+2)/(a1*numChrom) + a2/a1^2
    e1 = c1/a1
    e2 = c2/(a1^2+a2)
    
    ## Create a dataframe to write the statistics
    out.tab        <- list()
    #names(out.tab) <-c("pop","S","pi","ThetaW")
    out.tab$pop <- pop
    # total number of sites
    
    # Number of variable sites is the sum of the SFS excluding corners
    out.tab$S<-sum(sfs[2:(length(sfs)-1)])
    
    # Calculate Theta as K/an (where an is the harmonic numebr of n-1)
    W<-(out.tab$S/(harmonicNumber))
    out.tab$ThetaW<- W/seq_len
    # calculate pi
    P<-numChrom/(numChrom-1)*2*sum(sfs*p*(1-p))
    out.tab$pi <- P/seq_len
    
    out.tab$TajimaD<-(P - W)/sqrt(e1*out.tab$S+e2*out.tab$S*(out.tab$S-1))
    return(out.tab)
  }))
}

sfs1D <- function(pops, pop,GTs,seq_len){
  n_inds <- length(which(pops==pop))
  
  d <- table(apply(GTs[pops==pop,],2,sum))
  sfs <- d[match(1:(2*n_inds),names(d))]
  
  sfs <- as.vector(c(seq_len-sum(sfs),sfs))
  sfs[is.na(sfs)] <- 0
  return(sfs)
}

f.dist <- function(i, j){
  sqrt((i[1]-j[1])^2 + (i[2]-j[2])^2 )
}


#' finds the origin for a set of Populations
find_origin <- function(psi=tdoa.psi[,5],dist, pct=0.01,xlen=20, ylen=20){
  
  
  dists <- dist[,3]
  ij <- dist[,c(1,2)]
  
  mdlq <- list()
  for(i in 1:xlen){
    mdlq[[i]] <- list()
  }
  
  mdls <- apply(dists,1,function(d){
    l = lm( psi ~ d) 
    l$f.e = .5 * pct/ l$coefficients[2]
    l$rsq = summary(l)$r.squared
    return(l)
  })
  
  lapply(1:nrow(ij),function(r){
    i <- ij[r,1]
    j <- ij[r,2]
    mdlq[[i]][[j]] <<- mdls[[r]]
  })
  
  d0 <- matrix(NA, ncol=xlen, nrow=xlen)
  rsq <- matrix(NA, ncol=xlen, nrow=xlen)
  p <- matrix(NA, ncol=xlen, nrow=xlen)
  
  d0[ij] <- sapply(mdls, function(x)x$f.e)
  rsq[ij] <- sapply(mdls, function(x)x$rsq)
  p[ij] <- sapply(mdls, function(x)summary(x)$coefficients[2,4])
  
  res <- list(d0=d0, rsq=rsq, p, bbox=bbox, xlen=xlen, 
              ylen=ylen)
  class(res) <- 'origin.results'
  
  
  return(res)
  
}


#' gets a bounding box around the populations in samples
get.sample.bbox <- function(samples){
  s <- c('longitude', 'latitude')
  mins <- apply(samples[,s],2,min)
  maxs <- apply(samples[,s],2,max)
  return(cbind(mins, maxs))
}


#' preps data for tdoa analysis
#'
#' makes a data structure of format xi, yi, xj, yj, psi
#' 
#' @param pop.coords coordination data with latitude and longitude
#' @param all.psi matrix with psi values
#' @param region region or set of regions to run analysis on
#' @param countries countries to restrict analysis to
#' @param xlen number of x points
#' @param ylen number of y points
#'
prep.tdoa.data <- function(coords, psi){
  locs <- coords[,c('longitude', 'latitude')]
  n.locs <- nrow(locs)
  
  tdoa.data <- c()
  for(i in 1:n.locs){
    for(j in (i+1):n.locs){
      if( i>=j | j >n.locs) break
      tdoa.data <- rbind( tdoa.data, c(locs[i,], locs[j,], psi[i,j]))
    }
  }
  
  tdoa.data <- matrix(unlist(tdoa.data), ncol=5)
  
  return(tdoa.data)
}



#' Finds origin for a region
#' @param tdoa.data Object of type tdoa.data
#' @param bbox bounding box describing the location where inference should
#'   be done in
#' @param pop.coords coordinations
#' @param pct threshold parameter in model.1d
#' @param xlen, ylen parameters describing the number of points to use
#' @param exclude.ocoean boolean, whether points not on land should be
#'    excluded
#' @param exclude.land boolean, whether points on land should be
#'    excluded
single.origin <- function(tdoa.data, bbox,  pop.coords,
                          pct=0.01,
                          xlen=100, ylen=100, 
                          exclude.ocean=T,
                          exclude.land=F,
                          ...){
  
  
  #define locs for estimate
  s1<-seq(bbox[1,1],bbox[1,2], length.out=xlen)
  s2<-seq(bbox[2,1],bbox[2,2], length.out=ylen)
  coords <- expand.grid(s1,s2)
  ij <- expand.grid(1:length(s1), 1:length(s2))
  
  if(exclude.ocean){
    cc <- coords2country(coords)
    to.keep <- !is.na(cc)
  } else {
    to.keep <- rep(T, nrow(coords))
  }
  if(exclude.land){
    cc <- coords2country(coords)
    to.keep <- is.na(cc)
  }
  # init output
  d0 <- matrix(NA, ncol=ylen, nrow=xlen)
  rsq <- matrix(NA, ncol=ylen, nrow=xlen)
  mdlq <- list()
  for(i in 1:xlen){
    mdlq[[i]] <- list()
  }
  
  
  for(r in 1:nrow(coords)){
    i <- ij[r,1]
    j <- ij[r,2]
    x <- coords[r,1]
    y <- coords[r,2]
    
    if(to.keep[r]){
      mdl <- model.1d( xy=c(x,y), data=tdoa.data, pct=pct)
      d0[i, j] <- mdl$f.e
      rsq[i, j] <- mdl$rsq
      mdlq[[i]][[j]] <- mdl
    }
  }
  res <- list( d0=d0, rsq=rsq, mdlq=mdlq, bbox=bbox, xlen=xlen, 
               ylen=ylen, coords=pop.coords)
  class(res) <- 'origin.results'
  return(res)
}

# tdoa.data <- tdoa.psi
# coords <- pop.coords
single.origin_mod <- function(tdoa.data, bbox,  pop.coords,
                          pct=0.01,
                          xlen=100, ylen=100, 
                          exclude.ocean=F,
                          exclude.land=F,
                          ...){
  
  
  #define locs for estimate
  s1<-seq(bbox[1,1],bbox[1,2], length.out=xlen)
  s2<-seq(bbox[2,1],bbox[2,2], length.out=ylen)
  coords <- expand.grid(s1,s2)
  ij <- expand.grid(1:length(s1), 1:length(s2))
  
  if(exclude.ocean){
    cc <- coords2country(coords)
    to.keep <- !is.na(cc)
  } else {
    to.keep <- rep(T, nrow(coords))
  }
  if(exclude.land){
    cc <- coords2country(coords)
    to.keep <- is.na(cc)
  }
  # init output
  d0 <- matrix(NA, ncol=ylen, nrow=xlen)
  rsq <- matrix(NA, ncol=ylen, nrow=xlen)
  p_val <- matrix(NA, ncol=ylen, nrow=xlen)
  
  # for(i in 1:xlen){
  #   mdlq[[i]] <- list()
  # }
  # 
  r <- 1
  for(r in 1:nrow(coords)){
    i <- ij[r,1]
    j <- ij[r,2]
    x <- coords[r,1]
    y <- coords[r,2]
    
    if(to.keep[r]){
      mdl <- model.1d( xy=c(x,y), data=tdoa.data, pct=pct)
      #summary(mdl)
      d0[i, j] <- mdl$f.e
      rsq[i, j] <- mdl$rsq
      p_val[i, j] <- summary(mdl)$coef[2,4]
      #mdlq[[i]][[j]] <- mdl
    }
  }
  res <- list( d0=d0, rsq=rsq, p_val=p_val, bbox=bbox, xlen=xlen, 
               ylen=ylen, coords=pop.coords)
  class(res) <- 'origin.results'
  return(res)
}


#' gets a bounding box around the populations in samples
get.sample.bbox <- function(samples){
  s <- c('longitude', 'latitude')
  mins <- apply(samples[,s],2,min)
  maxs <- apply(samples[,s],2,max)
  return(cbind(mins, maxs))
}


#' calculates distance from a given point
#' @param xy coordinate of point to evalute function at
#' @param data a 5 col data frame with columns xi, yi, xj, yj, psi, with 
#'   coords from the two sample location and their psi statistic.
#'   best generated unsg prep.tdoa.data
#' @param f.dist the distance function to use. the default 'haversine'
#'   uses the haversine distance. Alternatively, 'euclidean' uses 
#'   Euclidean distance
#' @return an object of type lm describing fit
#' 

#model.1d( xy=c(x,y), data=tdoa.psi, pct=pct) 
model.1d <- function(xy, data, pct=0.01, f.dist="haversine"){
  if (f.dist=="euclidean"){
    f.dist <- function(i, j){
      sqrt((ix -jx)^2 + (iy-jy)^2 )
    }
  }else{ if(f.dist=="haversine"){
    f.dist <- distHaversine
  }}
  
  y = xy[2] 
  x = xy[1]
  ixy = data[,1:2]
  jxy = data[,3:4]
  psi = data[,5]
  
  
  d = f.dist(ixy, c(x,y)) - f.dist(jxy, c(x,y))
  l = lm( psi ~ d ) 
  l$f.e = .5 * pct/ l$coefficients[2]
  l$rsq = summary(l)$r.squared
  
  return (l)
}



#' plots the output of find.origin as a heatmap using .filled.contour
#' @param x an object of type origin.results, as obtained by 
#' @param n.levels the number of color levels
#' @param color.function a function that takes an integer argument and
#'    returns that many colors
#' @param color.negative a single color to be used for negative values
#' @param add.map boolean whether a map should be added
#' @param add.likely.origin boolean, whether origin should be marked with an
#' @param asp aspect ratio, set to 1 to keep aspect ratio with plot
#'   X
#' @export
plot.origin.results <- function(x, n.levels=100, color.function=heat.colors,
                                color.negative='grey',
                                add.map=T,
                                add.samples=T,
                                add.sample.het=T,
                                add.likely.origin=T,
                                asp=1,
                                ...){
  plot.default(NA, xlim=x$bbox[1,], ylim=x$bbox[2,], 
               xlab="", ylab="",
               xaxt="n", yaxt="n",
               xaxs='i', yaxs='i',
               asp=asp, ...)
  
  s1<-seq(x$bbox[1,1],x$bbox[1,2],length.out=x$xlen)
  s2<-seq(x$bbox[2,1],x$bbox[2,2],length.out=x$ylen)
  
  rel <- (x[[1]]>0) * (x[[2]]-min(x[[2]],na.rm=T)) /
    (max(x[[2]],na.rm=T)-min(x[[2]],na.rm=T))+0.001
  
  levels <- c(0,quantile(rel[rel>0.001], 0:n.levels/n.levels, na.rm=T) + 
                1e-6 * 0:n.levels/n.levels)
  cols <- c(color.negative, color.function(n.levels-1))
  
  
  .filled.contour(s1, s2, rel, levels, cols)
  
  
  # rect(x$bbox[1,1], x$bbox[2,1], x$bbox[1,2], x$bbox[2,2], border='black',
  #      lwd=2, col=NULL)
  
  if(add.likely.origin){
    points(summary(x)[,1:2], col='black', pch='x', cex=2)
  }
  if(add.map){
    require(rworldmap)
    m <- getMap("high")
    plot(m, add=T, lwd=1.3)
  }
  
  if(add.sample.het){
    samples <- x$coords
    hets <- (samples$hets - min(samples$hets) )/(
      max(samples$hets) - min(samples$hets))
    points( samples$longitude, samples$latitude,
            pch=16, cex=3, col=grey(hets) )
    points( samples$longitude, samples$latitude,
            pch=1, cex=3, col="black",lwd=2 )
  }
  else if(add.samples){
    points( samples$longitude, samples$latitude,
            pch=16, cex=1, col="black",lwd=1 )
  }
}
#' calculates the psi matrix for a pop object
#' 
#' @param pop population data object from make.pop
#' @param n the sample size which we downsample to
#' @param resampling mode of resampling. Currently, only
#'    hyper is supported
#' @return A n x n matrix of psi values
#' @example examples/example_1.r
#' @export
get.all.psi <- function(pop, n=2,
                        subset=NULL,
                        resampling="hyper"){
  #this function calculates psi for columns i,j, both resampled down to
  # n samples
  
  if(is.null(subset))
    subset <- 1:pop$n
  if( is.logical( subset) )
    subset <- which(subset)
  n.pops <- length(subset)
  mat = matrix( 0, nrow=n.pops, ncol=n.pops )
  for(i in 1:(n.pops-1)){
    for(j in (i+1):n.pops){
      ii <- subset[i]
      jj <- subset[j]
      ni <- unlist(pop$ss[ii,])
      nj <- unlist(pop$ss[jj,])
      fi <- unlist(pop$data[ii,])
      fj <- unlist(pop$data[jj,])
      mat[j,i] <- get.psi( ni, nj, fi, fj, 
                           resampling=resampling, n=n )
      mat[i,j] <- -mat[j,i]
      #print( c(ii, jj))
    }
  }
  
  return(mat)
}

#' the psi statistic calculation
#' the actual calculation of the psi statistic for
#' a single pair of populations. The function requires
#' 4 vectors, all of length equal to the number of snps
#' to be analyzed.
#
#'
#' @param ni the number of sampled haplotypes in population i
#' @param nj the number of sampled haplotypes in population j
#' @param fi the number of derived alleles in population i
#' @param fj the number of derived alleles in population j
#' @param n the number of samples to downsample to
#' @example examples/example_1.r
#' @return psi a matrix of pairwise psi values
get.psi <- function (ni, nj, fi, fj,
                     n=2, resampling="hyper"){
  
  fn <- cbind( fi, ni, fj, nj)
  
  tbl <- table(as.data.frame(fn))
  tbl <- as.data.frame( tbl )
  tbl <- tbl[tbl$Freq > 0,]
  
  tbl$fi <- as.integer(as.character( tbl$fi ))
  tbl$fj <- as.integer(as.character( tbl$fj ))
  tbl$ni <- as.integer(as.character( tbl$ni ))
  tbl$nj <- as.integer(as.character( tbl$nj ))
  
  to.exclude <- tbl$fi == 0 | tbl$fj == 0 | 
    tbl$ni < n | tbl$nj < n
  
  tbl <- tbl[! to.exclude, ]
  
  if(nrow(tbl)==0){ return(NaN)}
  
  
  poly.mat <- matrix(0,nrow=n+1, ncol=n+1)
  poly.mat[2:(n+1),2:(n+1)] <- 1
  poly.mat[n+1,n+1] <- 0
  
  #psi.mat is the contribution to psi for each entry
  psi.mat <- outer(0:n,0:n,FUN=function(x,y)(y-x))
  psi.mat[1,] <- 0
  psi.mat[,1] <- 0
  
  f.contribution <- function(row, b=2){
    a <- 0:b
    f1 <- row[1]
    n1 <- row[2]
    f2 <- row[3]
    n2 <- row[4]
    cnt <- row[5]
    q1 <- choose(b, a) * choose(n1-b, f1-a)/choose(n1,f1)
    q2 <- choose(b, a) * choose(n2-b, f2-a)/choose(n2,f2)
    return( cnt * outer(q2, q1) )
  }
  
  
  resampled.mat <- matrix(rowSums(apply(tbl, 1,
                                        f.contribution)),nrow=n+1)
  
  return( sum(resampled.mat * psi.mat) / sum(resampled.mat * poly.mat) )
}
