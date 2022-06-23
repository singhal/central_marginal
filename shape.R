library(ape)
library(raster)
library(tidyr)
library(maptools)
library(sp)
library(dplyr)
library(sf)
library(rgeos)
library(lwgeom)

get_map <- function(sp_new) {
  file = paste0("~/Dropbox (Personal)/publications/center_marginal/data/diversity/", sp_new, ".coverage_filtered_10.diversity.csv")

  tmp = read.csv(file, stringsAsFactors=F)
  tmp2 = left_join(tmp, ind, by = c("individual" = "SAMPLE_ID"))
  
  # individuals with missing lat long
  missind = tmp2[!complete.cases(tmp2$LAT), "individual"]
  
  # remove individuals with too much missing data 
  # based on bootstrap results
  tmp3 = tmp2 %>% filter(denominator >= 150000)
  tmp4 = tmp3[!is.na(tmp3$LAT),]
  
  rangedir = "~/Dropbox (Personal)//publications/Sphenomorphine_Gene_Flow/data/geographic_ranges/OTU_ranges/"
  rangefile = paste(rangedir, sp_new, ".shp", sep="")
  
  if (file.exists(rangefile)) {
    range = st_read(rangefile, crs = '+proj=longlat +datum=WGS84')
  } else {
    range = NULL
  }
  return(list(range, tmp4))
}

# combine these two data frames
ind1 = read.csv('~/Dropbox (Personal)//macroevolution/eco_IBD_oz/data/metadata/individual_data_nomissing22Oct15.csv')
ind1a = ind1[, c("sample_id", "lat", "lon")]
names(ind1a) = c("SAMPLE_ID", "LAT", "LON")
ind2 = read.csv("~/Dropbox (Personal)//publications/Sphenomorphine_Gene_Flow/data/metadata/spheno_ind_data.csv", stringsAsFactors = F)
ind1b = ind1a[! ind1a$SAMPLE_ID %in% ind2$SAMPLE_ID, ]
ind = rbind(ind1b, ind2[, c("SAMPLE_ID", "LAT", "LON")])

raw = readRDS("~/Dropbox (Personal)/publications/center_marginal/data/correlation.Rds")
d = read.csv("~/Dropbox (Personal)/publications/center_marginal/data/correlation_centerdist.csv",
             stringsAsFactors = F)
names(d) = c("species", "pval", "rho", "nind")
  
res = vector("list", length = nrow(d))
for (i in 1:nrow(d)) {
  cat(d$species[i], "\n")
  x = get_map(d$species[i])
  
  range = st_make_valid(x[[1]])
  
  poly_points <- st_segmentize(range, dfMaxLength = 100) %>% 
    st_coordinates() %>% 
    as.data.frame() %>% 
    dplyr::select(X, Y) %>% 
    st_as_sf(coords = c("X", "Y"), 
             crs = '+proj=longlat +datum=WGS84')
  center = st_centroid(range)
  rand_pts = poly_points[sample(1:nrow(poly_points), size = 100),]
  dists = st_distance(center, rand_pts)
  dist_cv = sd(dists) / mean(dists)
  ex = as.vector(extent(range))
  long_ex = distHaversine(c(ex[1], center$geometry[[1]][2]), 
                c(ex[2], center$geometry[[1]][2]))
  lat_ex = distHaversine(c(center$geometry[[1]][1], ex[3]), 
                c(center$geometry[[1]][1], ex[4]))
  
  pts = x[[2]] %>% dplyr::select("LON", "LAT") %>% 
          st_as_sf(coords =c("LON", "LAT"),
           crs = '+proj=longlat +datum=WGS84')
  buff_pts_x = st_make_valid(st_union(st_buffer(pts, 100000)))
  buff_pts = st_intersection(buff_pts_x, range)
  cov = st_area(buff_pts) / st_area(range)
  
  res[[i]] = list(dist_cv, long_ex, lat_ex, cov)
}

shape = as.data.frame(apply(do.call("rbind", res), 2, as.numeric))
# eccentricity measure, longitude extent, latitude extent, 
# and sampling coverage
names(shape) =  c("dist_cv", "long_ex", "lat_ex", "cov")
shape$ex = shape$long_ex / shape$lat_ex
shape$ex2 = sapply(shape$ex, function(x) { ifelse(x < 1, x, 1/x)} )
shape$species = d$species
write.csv(shape, "~/Desktop/shape.csv",
          row.names = F, quote = F)
