## here we are comparing the original implementation of how to calculate psi
## and compare to much faster implementation where the resulting psi-matrix is transposed

# library(viridis)
library(data.table)
library(tidyverse)
# library(cowplot)
library(parallel)
library(SNPRelate)
library(geosphere)
library(raster)
# library(ggplot2)
# need to add this package for SnpMatrix
library(snpStats)
library(rangeExpansion)

## all necessary functions are here, including original functions from 
## https://github.com/BenjaminPeter/rangeexpansion/tree/master (v.0.0.0.9000)
# note: by original, this includes functions with error
source("~/Dropbox (Personal)/research_projects-FINISHED/Singhal_etal_CentralMarginal/error/Functions2.R")

setwd("~/Dropbox (Personal)/research_projects-FINISHED/Singhal_etal_CentralMarginal/data/range_exp/")

dd = list.files(".", pattern = "coords.csv", full.names = T)
res = vector('list', length(dd))

for (i in 1:length(dd)) {
  ## read in coordinates
  coords <- fread(dd[i])
  
  ## read in gentype matrix
  gfile = gsub(".coords", "", dd[i])
  GTs <- read.csv(gfile, header = F)
  # dim(GTs)
  # GTs[, 1:10]
  
  GTs <- data.frame(apply(GTs,2, as.numeric))
  rownames(GTs) <- coords$id
  GTs <- GTs[,-1]
  
  ## only keep individuals with <80% missing data
  keep <-which(apply(GTs,1,function(x) length(which(is.na(x))))/ncol(GTs)<0.8)
  GTs <- GTs[keep,]
  coords <- coords[keep,]
  
  ## write files that can be used by functions from rangeExpansion package
  write.table(x = GTs,file = "GT.snapp",row.names = TRUE,quote = FALSE,col.names = FALSE,sep=",")
  write.table(coords,file = "coord.csv",row.names = FALSE,quote = FALSE,col.names = TRUE,sep=",")
  
  ## use functions from rangeExpansion to prepare data
  raw.data <- load.data.snapp(snapp.file="GT.snapp",
                              coords.file="coord.csv",
                              n.snp=-1,
                              sep=',', ploidy=2)
  
  
  pop_sim <- make.pop(raw.data, ploidy=2)
  
  ## now estimate psi with new function using mclapply and 
  ## calculates significance values using the binomial test
  # t1 <- Sys.time()
  PSIs <- get.all.psi.mc.bin(pop_sim, n=2, cores=4)
  # t2 <- Sys.time()
  # time_mod <- difftime(t2,t1)
  
  ## this is the original implementation 
  # WITH THE CODING ERROR
  # t1 <- Sys.time()
  # psi_org <- get.all.psi(pop_sim,n=2)
  # t2 <- Sys.time()
  # time_org <- difftime(t2,t1)
  
  ## new code is ~50 times faster 
  # time_org * 60/time_mod 
  
  ## calculate genetic diversity from sfs, not heterozygosity as calculated by 
  ## rangeExpansion package because high proportion of missing data
  # pops <- paste(pop_sim$coords$longitude,pop_sim$coords$latitude)
  # pops2 <- coords[outgroup==0,paste(longitude,latitude)]
  # div <- div_from_sfs(pop_sim$coords$pop[match(pops2,pops)],GTs,ncol(GTs))
  # pis <- div$pi
  
  # ij <- as.matrix(expand.grid(1:length(pis),1:length(pis)))
  
  # delta_pi <- matrix(NA, ncol = length(pis),nrow = length(pis))
  
  #delta_pi[as.matrix(ij)] <- apply(ij,1,function(x){
  #  pis[x[1]]-pis[x[2]]
  #})
  
  # plot(delta_pi,PSIs$psi) # negatively correlated as should be expected
  # plot(delta_pi,psi_org) # positively correlated using original psi functions
  
  ## TDOA analyses
  # tdoa.psi <- prep.tdoa.data(pop_sim$coords, PSIs$psi)
  # tdoa.psi_org <- prep.tdoa.data(pop_sim$coords, psi_org)
  # tdoa.d_pi <- prep.tdoa.data(pop_sim$coords, delta_pi)
  
  # bbox <- get.sample.bbox(pop_sim$coords)
  
  # origin_psi_org <- single.origin_mod(tdoa.psi_org, bbox,  pop.coords=coords,
  #                                 pct=0.01,
  #                                 xlen=100, ylen=100, 
  #                                 exclude.ocean=F,
  #                                 exclude.land=F)
  # origin_psi <- single.origin_mod(tdoa.psi, bbox,  pop.coords=coords,
  #                                 pct=0.01,
  #                                 xlen=100, ylen=100, 
  #                                 exclude.ocean=F,
  #                                 exclude.land=F)
  
  avgpsi = mean(lower.tri(PSIs$psi), na.rm = T)
  results1 <- run.regions(region = list("REGION_1"), pop=pop_sim, psi = PSIs$psi, xlen=10, ylen=20)
  res[[i]] = list(results1, avgpsi)
  
  # old results with old psi with typo
  # results2 <- run.regions(region = list("REGION_1"), pop=pop_sim, psi = psi_org, xlen=10, ylen=20)
  # summary.origin.results(results2$tbl[[1]])
  
  # origin_D_pi <- single.origin_mod(tdoa.d_pi, bbox,  pop.coords=coords,
  #                                 pct=0.01,
  #                                 xlen=100, ylen=100, 
  #                                 exclude.ocean=F,
  #                                exclude.land=F)
  message("done with ", dd[i])
  }


# #save.image()
# ## prepare to plot data points 
# tmp <- data.table(pop_sim$coords,div[,.(S,ThetaW,pi)],psi=rowMeans(PSIs$psi,na.rm=TRUE))
# scale01 <- function(x) (x-min(x))/max((x-min(x)))
# tmp[,y:=scale01(latitude)]
# tmp[,x:=scale01(longitude)]
# 
# col = viridis(100)[as.integer(scale01(tmp$pi)*99+1)]
# ## now print Supplementary Figure ##
# pdf("./Figures/Transposed_psi.pdf",height = 4,width = 12)
# par(mfcol=c(1,3))
# 
# mat <- (origin_psi_org$d0>0)*origin_psi_org$rsq
# mat <- t(mat)[ncol(mat):1,]
# plot(raster(mat),main="a | Original rangeexpansion R-package")
# origin <- which.max((origin_psi_org$d0>0)*origin_psi_org$rsq)
# x <- row(mat)[origin]/100
# y <- col(mat)[origin]/100
# points(tmp$x,tmp$y,col=col,pch=20,cex=2)
# points(x,y,cex=5,pch=3,col="salmon")
# points(x,y,cex=5,pch=20,col="salmon")
# 
# 
# mat <- (origin_psi$d0>0)*origin_psi$rsq
# mat <- t(mat)[ncol(mat):1,]
# plot(raster(mat),main=expression("b | "~Psi~"-matrix transposed"))
# origin <- which.max((origin_psi$d0>0)*origin_psi$rsq)
# x <- row(mat)[origin]/100
# y <- col(mat)[origin]/100
# points(tmp$x,tmp$y,col=col,pch=20,cex=2)
# points(x,y,cex=5,pch=3,col="salmon")
# points(x,y,cex=5,pch=20,col="salmon")
# 
# mat <- (origin_D_pi$d0<0)*origin_D_pi$rsq
# mat <- t(mat)[ncol(mat):1,]
# plot(raster(mat),main=expression("c | Origin based on "~Delta[HET]))
# origin <- which.max((origin_D_pi$d0<0)*origin_D_pi$rsq)
# x <- row(mat)[origin]/100
# y <- col(mat)[origin]/100
# points(tmp$x,tmp$y,col=col,pch=20,cex=2)
# points(x,y,cex=5,pch=3,col="salmon")
# points(x,y,cex=5,pch=20,col="salmon")
# dev.off()

saveRDS(res, "../../error/redone_rangeexp.Rds")

res2 = lapply(res, function(x) { summary.origin.results(x[[1]]$tbl[[1]]) })
res3 = do.call("rbind", res2)
sps = gsub("^.*\\/", "", dd)
sps = gsub(".coords.csv", "", sps)
res3$species = sps

a = read.csv("../dadi_correlation_centerdist.csv",
             stringsAsFactors = F)
a2 = full_join(res3, a, by = c("species" = "species"))


ranges = vector("list", length = nrow(a2))
for (i in 1:nrow(a2)) {
  sp_new = a2[i, "species"]
  rangedir = "~/Dropbox (Personal)/research_projects-FINISHED/Singhal_etal_SphenoGeneFlow/data/geographic_ranges/OTU_ranges/"
  rangefile = paste(rangedir, sp_new, ".shp", sep="")
  if (file.exists(rangefile)) {
    sf::sf_use_s2(FALSE)
    range = sf::st_read(rangefile) %>% sf::st_set_crs(4979)
      
    centroid = sf::st_centroid(range)
    
    rangesize = sum(sf::st_area(range)) / 1e6 # sq km
    origin = sf::st_sfc(sf::st_point(as.numeric(a2[i, c("longitude", "latitude")]))) %>%
      sf::st_set_crs(4979)
    distance1 = sf::st_distance(origin, centroid) / 1000 # km
    
    distance2 = sf::st_geometry(obj = range) %>% 
      sf::st_cast(to = 'MULTILINESTRING') %>% 
      sf::st_distance(y = origin)
    distance2 = distance2 / 1000

    ranges[[i]] = c(as.numeric(rangesize), as.numeric(distance1[1,1]),
                    as.numeric(distance2[1,1]))
  }
}

geo = as.data.frame(do.call(rbind, ranges))
names(geo) = c("range_size", "distance1", "distance2")
geo$ratio = geo$distance1 / (geo$distance1 + geo$distance2)

a3 = cbind(a2, geo)
a3$range_sig = ifelse(a3$pval.x < 0.05, TRUE, FALSE)
a3 = a3 %>% filter(nind > 9) 
write.csv(a3, "../../error/demography.csv",
          quote = F, row.names = F)

a3 = read.csv("../../error/demography.csv",
              stringsAsFactors = F)

table(a3$range_sig)
mean(a3$ratio)

### remake figure 5
a = ggplot(a3) + geom_histogram(aes(as.numeric(param1))) +
  xlab("ratio of curr. Ne to anc. Ne") + theme_classic()
d2 = a3 %>% filter(range_sig == TRUE)
b = ggplot(d2, aes(rsq)) + geom_histogram(bins = 10) +
  xlab("strength of range expansion") + theme_classic()
plt = cowplot::plot_grid(a, b, ncol = 2, align = "vh", labels = c("A", "B"))
cowplot::save_plot("../../error/Figure5.png", 
          plt, base_height = 4, 
          base_width = 5, ncol = 2)

### remake figure S6
d = ggplot(a3 %>% filter(range_sig == TRUE), aes(ratio)) + geom_histogram() + 
  xlab("relative center-edge distance of origin") + theme_classic()
cowplot::save_plot("../../error/FigureS6.png",
          d, base_width = 4)
