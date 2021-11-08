library(readr)
library(maptools)
data(wrld_simpl)
library(dplyr)
library(raster)

read_csv_file <- function(x) {
  x1 = read_csv(x)
  sp = gsub(".coverage.*", "", x)
  sp = gsub(".*\\/", "", sp)
  x1$species = sp
  return(x1)
}

files = list.files("~/Dropbox (Personal)/center_marginal/data/diversity/", pattern="*.csv", full.names=T)
files2 = lapply(files, read_csv_file)
files3 = do.call("rbind", files2)
spcts = table(files3$species)
keep = names(spcts)[spcts >= 10]
d = files3[files3$species %in% keep, ]

# combine these two data frames
ind1 = read.csv('~/Dropbox (Personal)/macroevolution/eco_IBD_oz/data/metadata/individual_data_nomissing22Oct15.csv')
ind1a = ind1[, c("sample_id", "lat", "lon")]
names(ind1a) = c("SAMPLE_ID", "LAT", "LON")
ind2 = read.csv("~/Dropbox (Personal)/publications/Sphenomorphine_Gene_Flow/data/metadata/spheno_ind_data.csv", stringsAsFactors = F)
ind1b = ind1a[! ind1a$SAMPLE_ID %in% ind2$SAMPLE_ID, ]
ind = rbind(ind1b, ind2[, c("SAMPLE_ID", "LAT", "LON")])

dd = left_join(d, ind, by = c("individual" = "SAMPLE_ID"))

oz_extent =  extent(wrld_simpl[wrld_simpl$NAME == 'Australia',])

# the summary raster
all_species = raster(oz_extent)
res(all_species) = c(1/12, 1/12)
all_species = setValues(all_species, 0)

for (i in 1:length(keep)) {
  spname = gsub("Ctenotus_", "C_", keep[i])
  spname = gsub("Lerista_", "L_", spname)
  
  mapname = paste0("~/Dropbox (Personal)/macroevolution/eco_IBD_oz/data/geography/new_ranges/",
                   spname, ".shp")
  map = readShapePoly(mapname)
  
  # convert the map to a raster
  map_r = raster(oz_extent)
  res(map_r) = c(1/12, 1/12)
  map_r = rasterize(map, map_r)
  
  # set empty values to 0
  map_r[is.na(map_r)] <- 0
  # add rasters
  all_species = all_species + map_r
  
  cat(i, '\n')
}


### get coast
png("~/Dropbox (Personal)/publications/Center_Marginal/figures/map_of_allranges.png", 
    width=5, height=3, res = 300, units = "in")
par(mar=c(0.5, 0.5, 0.5, 0.5))
plot(all_species, xlim=c(113, 154), ylim=c(-40, -10), axes=F, col=colorRampPalette(brewer.pal(9,"Blues"))(20), box=F)
points(dd[, c('LON','LAT')], pch=4, cex=0.2)
dev.off()