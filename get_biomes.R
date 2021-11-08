library(raster)
library(sp)
library(maptools)
library(rgeos)
library(geosphere)
library(raster)
library(RColorBrewer)

biomes = readShapePoly("~/Dropbox/macroevolution/eco_IBD_oz/data/ecoregions/terr-ecoregions-aus_combined/terr_ecoregions-aus_combined.shp")
# simplify this so that it takes less time to process
x = gSimplify(biomes, 0.1, topologyPreserve=T)
x$biome = as.character(biomes$WWF_MHTNAM)

# rename the biomes so that it is easier to deal with
biomes = vector("list", length(x))
names(biomes) = x$biome
bnames1 = c("Tropical and Subtropical Moist Broadleaf Forests", "Montane Grasslands and Shrublands", "Temperate Broadleaf and Mixed Forests", "Tropical and Subtropical Grasslands, Savannas and Shrublands", "Deserts and Xeric Shrublands", "Mediterranean Forests, Woodlands and Scrub", "Temperate Grasslands, Savannas and Shrublands")
bnames2 = c("tropical forests", "montane", "temperate forests", "tropical grasslands", "desert", "Mediterranean forests", "temperature grasslands")
for (i in 1:length(bnames1)) {
  names(biomes) = gsub(bnames1[i], bnames2[i], names(biomes))
}
# helps deal with issues around illegal spatial geometries
for (i in 1:length(x)) {
  a = x[i, ]
  biomes[[i]] = gBuffer(SpatialPolygons(a@polygons,proj4string=a@proj4string), width=0)
}

d3 = read.csv("~/Dropbox/center_marginal/data/correlation_centerdist.csv", 
              stringsAsFactors = F)
bdf = data.frame(species = d3$X, biome = NA, stringsAsFactors = F)
for (i in 1:nrow(d3)) {
   sp_new = d3[i, "X"]
   rangedir = "/Users/sonal/Dropbox/publications/Sphenomorphine_Gene_Flow/data/geographic_ranges/OTU_ranges/"
   rangefile = paste(rangedir, sp_new, ".shp", sep="")
   range = readShapePoly(rangefile,
                           proj4string=CRS('+proj=longlat +datum=WGS84'))

  if (class(range) == 'SpatialPolygonsDataFrame') {
    range = SpatialPolygons(range@polygons,proj4string=range@proj4string)
  }
  
  # calculate all the overlaps in km2
  overlaps = rep(NA, length(biomes))
  names(overlaps) = names(biomes)
  for (jj in 1:length(biomes)) {
    overlap = gIntersection(gBuffer(range, width=0), biomes[[jj]])
    if (class(overlap) != 'NULL') {
      overlaps[jj] = sum(areaPolygon(overlap)) / 1e6
    }
  }	
  
  bdf[i, "biome"] = names(overlaps)[which(overlaps == max(overlaps, na.rm = T))]

  cat(i, "\n")
}
write.csv(bdf, "~/Dropbox/center_marginal/data/biomes.csv", row.names = F,
          quote = F)
