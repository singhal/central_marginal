library(ntbox)
library(raster)
library(rgl)
library(stringr)
library(dplyr)
library(maptools)

get_sp_points <- function(file) {
  tmp = read.csv(file, stringsAsFactors=F)
  sp = gsub("^.*\\/", "", file)
  sp = gsub(".coverage.*", "", sp)
  tmp$species = sp
  return(tmp[, c("individual", "species")])
}

######################################
# simple distance method based on PC
#####################################

# Euclidean distance
euc.dist <- function(x1, x2) sqrt(sum((x1 - x2) ^ 2))

# get gen div points
files = list.files("/Users/sonal/Dropbox/center_marginal/data/diversity/", 
                   pattern="*.csv", full.names=T)
inds = lapply(files, get_sp_points)
inds2 = do.call("rbind", inds)
loc = read.csv("~/Dropbox/center_marginal/data/spatial/gen_div_pts.csv",
               stringsAsFactors = F)
inds3 = left_join(inds2, loc, by = c("individual" = "SAMPLE_ID"))

# get rasters
envr = list.files("~/Dropbox/macroevolution/eco_IBD_oz/data/AUS_5arc/",
                  pattern = "tif", full.names = T)
envr = envr[grep("slope", envr, invert = T)]
envr = envr[grep("aspect", envr, invert = T)]
envr = envr[grep("alt", envr, invert = T)]
env = raster::stack(envr)

PCres = vector("list", length(pts))
names(PCres) = names(pts)
for (i in 1:length(pts)) {
  # first do a fairly naive estimate of distance
  rangedir = "/Users/sonal/Dropbox/publications/Sphenomorphine_Gene_Flow/data/geographic_ranges/OTU_ranges/"
  rangefile = paste(rangedir, names(pts)[i], ".shp", sep="")
  range = readShapePoly(rangefile,
                        proj4string=CRS('+proj=longlat +datum=WGS84'))
  # sample points
  randpts = spsample(range, 1000, type= "random")
  
  # actual points
  sppts = inds3 %>% 
    filter(species == names(pts)[i]) %>% 
    select(individual, LON, LAT)
  sppts2 = SpatialPoints(sppts[, c("LON", "LAT")], proj4string = crs(randpts))
  
  # get environmental data
  r_envpts = extract(env, randpts, df = T)
  r_envpts2 = r_envpts[complete.cases(r_envpts), ]
  s_envpts = raster::extract(env, sppts2, df = T)
  s_envpts$ind = sppts$individual
  s_envpts2 = s_envpts[complete.cases(s_envpts), ]
  
  # do pca
  envdata = rbind(r_envpts2, s_envpts2 %>% select(-ind))
  pcenv = prcomp(envdata, center = T, scale. = T)
  
  # get centroid
  center = apply(pcenv$x[1:nrow(r_envpts2), 1:6], 2, mean)
  
  # get euclidean distance of pts 
  s_envpts3 = pcenv$x[(nrow(r_envpts2) + 1): nrow(pcenv$x), 1:6]
  dists =  apply(s_envpts3, 1, euc.dist, center)
  PCres[[i]] = data.frame(ind = s_envpts2$ind, 
                          PC_eco_dist = dists,
                          prop6 = cumsum(pcenv$sdev^2 / sum(pcenv$sdev^2))[6])
  cat(i, "\n")
}
saveRDS(PCres, "~/Dropbox/center_marginal/data/spatial/PCres.rds")

######################################
# distance method based on MVE
#####################################

# get occ points
pts = readRDS("~/Dropbox/center_marginal/data/spatial/cleaned_occ_pts.Rds")

MVEres = vector("list", length(pts))
names(MVEres) = names(pts)
for (i in 9:length(pts)) {
  # split data by training and testing
  pts2 = data.frame(pts[[i]])
  # need to check this default split
  pts2$type = sample(c("test", "train"), nrow(pts2),
                     replace = T, prob = c(0.3, 0.7))
  pgL <- base::split(pts2, pts2$type)
  pg_train <- pgL$train
  pg_test <- pgL$test
  
  # extract data
  pg_etrain <- raster::extract(env, pg_train[,c("Longitude","Latitude")],
                               df=TRUE)
  pg_etrain = pg_etrain[, 2:25]
  pg_etest <- raster::extract(env, pg_test[,c("Longitude","Latitude")],
                               df=TRUE)
  pg_etest = pg_etest[, 2:25]
  
  # check default corr
  env_varsL <- ntbox::correlation_finder(cor(pg_etrain,method = "spearman"),
                                         threshold = 0.8,
                                         verbose = F)
  env_vars <- env_varsL$descriptors
  
  env_bg <- ntbox::sample_envbg(env, 10000)
  e_selct <- ntbox::ellipsoid_selection(env_train = pg_etrain,
                                        env_test = pg_etest,
                                        env_vars = env_vars,
                                        # what proportion of points to include
                                        # or, how much error is there in points?
                                        level = 0.99,
                                        nvarstest = c(2, 3, 4),
                                        env_bg = env_bg,
                                        omr_criteria = 0.06,
                                        proc = TRUE,
                                        parallel = T,
                                        comp_each = 100)
  
  # Best ellipsoid model for "omr_criteria" 
  bestvarcomb <- stringr::str_split(e_selct$fitted_vars,",")[[1]]
  
  # Ellipsoid model (environmental space)
  envsp1 = pg_etrain[,bestvarcomb]
  envsp2 = envsp1[complete.cases(envsp1), ]
  best_mod <- ntbox::cov_center(envsp2,
                                mve = T,
                                level = 0.99,
                                vars = 1:length(bestvarcomb))
  
  # for visualization
  # Projection model in geographic space
  # mProj <- ntbox::ellipsoidfit(env[[bestvarcomb]],
  #                             centroid = best_mod$centroid,
  #                             covar = best_mod$covariance,
  #                             level = 0.99,size = 3)
  # raster::plot(mProj$suitRaster)
  # points(pts2[,c("Longitude", "Latitude")],pch=20,cex=0.5)
  
  sppts = inds3 %>% 
    filter(species == names(pts)[i]) %>% 
    select(individual, LON, LAT)
  s_envpts = raster::extract(env, sppts %>% select(LON, LAT), df = T)
  s_envpts$ind = sppts$individual
  s_envpts2 = s_envpts[complete.cases(s_envpts), ]
  
  mhd <- stats::mahalanobis(s_envpts2[, bestvarcomb],
                            center = best_mod$centroid,
                            cov = best_mod$covariance)
  df = data.frame(MHD = mhd, ind = s_envpts2$ind, stringsAsFactors = F)
  mod_app <- e_selct %>% filter(om_rate_train<=1 & om_rate_test<=1)
  MVEres[[i]] = list(df, mod_app)
  
  cat(i, "\n")
}
saveRDS(MVEres, "~/Dropbox/center_marginal/data/spatial/MVEres.rds")
