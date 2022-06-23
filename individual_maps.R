library(dplyr)
library(tidyr)
library(maptools)
library(rgeos)
library(geosphere)
library(sp)

######################
# make individual maps
######################

count_inds <- function(file) {
  f = read.csv(file, stringsAsFactors = F)
  return(nrow(f))
}


get_map <- function(file) {
  tmp = read.csv(file, stringsAsFactors=F)
  tmp2 = left_join(tmp, ind, by = c("individual" = "SAMPLE_ID"))
  
  # individuals with missing lat long
  missind = tmp2[!complete.cases(tmp2$LAT), "individual"]
  
  # remove individuals with too much missing data 
  # based on bootstrap results
  tmp3 = tmp2 %>% filter(denominator >= 150000)
  tmp4 = tmp3[!is.na(tmp3$LAT),]
  
  sp_new = gsub('.*\\/', '', file)
  sp_new = gsub(".coverage.*", "", sp_new)
  cat(sp_new, "\n")
  
  
  rangedir = "~/Dropbox (Personal)//publications/Sphenomorphine_Gene_Flow/data/geographic_ranges/OTU_ranges/"
  rangefile = paste(rangedir, sp_new, ".shp", sep="")
  
  if (file.exists(rangefile)) {
    range = readShapePoly(rangefile,
                          proj4string=CRS('+proj=longlat +datum=WGS84'))
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

files = list.files("~/Dropbox (Personal)/publications//center_marginal/data/diversity/", pattern="*.csv", full.names=T)
ninds = unlist(lapply(files, count_inds))
files2 = files[ninds >= 9]

maps = vector("list", length(files2))
for (i in 1:length(files2)) {
  maps[[i]] = get_map(files2[i])
}

library(rnaturalearthdata)
library(rnaturalearth)

world = ne_countries(scale = "medium", returnclass = "sf")
oz = world[world$admin == "Australia", ]

sps = gsub(".*/", "", files2)
sps = gsub(".coverage.*", "", sps)
sp = readxl::read_xlsx("~/Dropbox (Personal)/publications/Center_Marginal/species_names.xlsx")
sp2 = pull(sp[match(sps, sp$OTU), "sp_name"])
names(sp2) = sps

cor = read.csv("~/Dropbox (Personal)/publications/Center_Marginal/data/correlation_centerdist.csv")
cor = cor %>% arrange(rho)
cor$sig = ifelse(cor$V1 < 0.05, TRUE, FALSE)
names(maps) = sps

pdf("~/Dropbox (Personal)/publications/Center_Marginal/figures/species_ranges2.pdf", height = 10, width = 10)
par(mfrow = c(5, 5))
for (i in 1:nrow(cor)) {
  sp = cor[i, "X"]
  if (nrow(maps[[sp]][[2]]) > 9) {
    par(mar = c(0, 0, 0, 0))
    plot(1, type="n", axes=FALSE, ann=FALSE,
         xlim = c(110, 155), ylim = c(-46, -8.5), 
         xlab = "", ylab = "")
    plot(oz[1], col = "white", add = T)
    
    rho = round(cor[i, "rho"], 2)
    txt = paste(sp2[sp], ": ", rho, sep = "")
    
    mtext(txt, side = 3, font = 3, cex = 0.7)
    
    if (cor[i, "rho"] < 0) {
      colx = "#4d9221"  
      } else {
      colx = "#c51b7d" 
      }
    if (cor[i, "sig"] == TRUE) {
      alphax = 1
    } else {
      alphax = 0.2
    }
    plot(maps[[sp]][[1]], col = alpha(colx, alphax),
         border = FALSE, add = T)
    points(maps[[sp]][[2]][, c("LON", "LAT")], pch = 21, bg = "white")
  }
}
dev.off()