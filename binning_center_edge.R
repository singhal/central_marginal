library(rgeos)
library(readr)
library(dplyr)

calculate_diversity <- function(file) {
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
  
  rangedir = "~/Dropbox (Personal)/publications/Sphenomorphine_Gene_Flow/data/geographic_ranges/OTU_ranges/"
  rangefile = paste(rangedir, sp_new, ".shp", sep="")
  if (file.exists(rangefile)) {
    range = maptools::readShapePoly(rangefile,
                          proj4string=CRS('+proj=longlat +datum=WGS84'))
    points = SpatialPoints(tmp4[,c("LON", "LAT")],
                           proj4string=CRS('+proj=longlat +datum=WGS84'))
    
    rangesize = sum(geosphere::areaPolygon(range)) / 1e6
    
    # distance from center
    centroid = gCentroid(range)
    distance1 = geosphere::distHaversine(points, centroid)
    distance2 = geosphere::dist2Line(p = points, line = range)[, 1]
    
    # tot dist
    totdist = distance1 + distance2
    # % from edge
    tmp4$edgedist = distance2 / totdist
    tmp4$core_edge = ifelse(tmp4$edgedist < 0.2, "edge", "core")
    xx = aov(tmp4$pi ~ tmp4$core_edge)
  } else {
    xx = NULL
  }
   return(list(tmp4, xx))
}

count_inds <- function(file) {
  f = read_csv(file)
  return(nrow(f))
}

files = list.files("~/Dropbox (Personal)/center_marginal/data/diversity/", 
                   pattern="*.csv", full.names=T)
ninds = unlist(lapply(files, count_inds))
files2 = files[ninds >= 10]

# combine these two data frames
ind1 = read.csv('~/Dropbox (Personal)/macroevolution/eco_IBD_oz/data/metadata/individual_data_nomissing22Oct15.csv')
ind1a = ind1[, c("sample_id", "lat", "lon")]
names(ind1a) = c("SAMPLE_ID", "LAT", "LON")
ind2 = read.csv("~/Dropbox (Personal)/publications/Sphenomorphine_Gene_Flow/data/metadata/spheno_ind_data.csv", stringsAsFactors = F)
ind1b = ind1a[! ind1a$SAMPLE_ID %in% ind2$SAMPLE_ID, ]
ind = rbind(ind1b, ind2[, c("SAMPLE_ID", "LAT", "LON")])

res = lapply(files2, calculate_diversity)
res2 = lapply(res, function(x) {x[[1]]})
means = rep(NA, length(res2))
for (i in 1:length(res2)) {
  x = res2[[i]] %>% group_by(core_edge) %>% summarize(pi = mean(pi)) %>% ungroup()
  means[i] = (x$pi[2] - x$pi[1]) / x$pi[1]
}

res3 = do.call("rbind", res2)
x = res3 %>% group_by(species) %>% summarize(range = (max(pi, na.rm=T) - min(pi, na.rm=T)) / min(pi, na.rm=T))
mean(x$range)


div = list.files("~/Dropbox (Personal)/macroevolution/eco_IBD_oz/data/pop_gen/divergence_cov5_results/divergence/",
                 full.names = T)
div2 = lapply(div, read.csv)
div3 = do.call("rbind", div2)

sp = read.csv("~/Dropbox (Personal)/macroevolution/eco_IBD_oz/data/clustering/Ctenotus_Lerista.clustering.revised3.csv")
sp$LatinName = gsub("C. ", "Ctenotus_", sp$LatinName)
sp$LatinName = gsub("L. ", "Lerista_", sp$LatinName)
sp$LatinName = gsub(" ", "_", sp$LatinName)
div3$species = sp[match(div3$cl, sp$GMYC_RAxML2), "LatinName"]
div4 = div3[div3$species %in% res3$species, ]
res3$core_edge = "CORE"
res3[res3$distance3 > 0.8, "core_edge"] = "EDGE"
div4$core_edge1 = res3[match(div4$ind1, res3$individual), "core_edge"]
div4$core_edge2 = res3[match(div4$ind2, res3$individual), "core_edge"]
div5 = div4 %>% filter(complete.cases(core_edge1),
                       complete.cases(core_edge2))
div5$type = NA
for (i in 1:nrow(div5)) {
  if (div5[i, "core_edge1"] == "EDGE" &  div5[i, "core_edge2"] == "EDGE") {
    div5$type[i] = "EDGE-EDGE"
  } else if (div5[i, "core_edge1"] == "CORE" &  div5[i, "core_edge2"] == "CORE") {
    div5$type[i] = "CORE-CORE"
  } else {
    div5$type[i] = "EDGE-CORE"
  }
}
library(ggplot2)
library(lme4)
library(cowplot)
library(car)
xx = ggplot(div5, aes(type, fst)) + geom_violin() + 
  theme_classic() + xlab("") + ylab(expression(F[ST]))
save_plot("~/Desktop/FST_ranges.png", xx)

a= lmer(fst ~ type + (1 | species), data = div5)
summary(a)
Anova(a)
confint(a)

resa = readRDS("~/Dropbox (Personal)/center_marginal/data/correlation.Rds")
res2 = resa[!unlist(lapply(resa, is.null))]
cor1 = as.data.frame(do.call(rbind, lapply(res2, get_corr, "distance1")))
cor1 = cor1[which(cor1$V3 > 9 & complete.cases(cor1$rho)), ]
cor1$mean = means
cor1 %>% filter(rho <0) %>% summarize(mean(mean))
# pvals = unlist(lapply(res, function(x) {summary(x)[[1]][["Pr(>F)"]][1]}))
