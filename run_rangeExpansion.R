library(sp)
library(rangeExpansion)
library(scales)
library(maptools)
library(rgeos)
library(geosphere)
library(sp)
library(viridis)
library(dplyr)
library(ggplot2)
library(cowplot)
theme_set(theme_cowplot())


# lat long data
ind1 = read.csv('/Users/sonal/Dropbox/macroevolution/eco_IBD_oz/data/metadata/individual_data_nomissing22Oct15.csv')
ind1a = ind1[, c("sample_id", "lat", "lon")]
names(ind1a) = c("SAMPLE_ID", "LAT", "LON")
ind2 = read.csv("/Users/Sonal/Dropbox/publications/Sphenomorphine_Gene_Flow/data/metadata/spheno_ind_data.csv", stringsAsFactors = F)
ind1b = ind1a[! ind1a$SAMPLE_ID %in% ind2$SAMPLE_ID, ]
ind = rbind(ind1b, ind2[, c("SAMPLE_ID", "LAT", "LON")])
  
# ingroup vs outgroup data
cl = read.csv('/Users/Sonal/Dropbox/macroevolution/eco_IBD_oz/data/clustering/Ctenotus_Lerista.clustering.revised3.csv', stringsAsFactors = F)
cl$LatinName = gsub("C. ", "Ctenotus_", cl$LatinName)
cl$LatinName = gsub("L. ", "Lerista_", cl$LatinName)
cl$LatinName = gsub(" ", "_", cl$LatinName)

# range expansion files
re = list.files("/Users/Sonal/Dropbox/center_marginal/data/range_exp/", full.names = T)
re = re[grep("coords", re, invert = T)]

res = vector('list', length(re))
for (i in 1:length(re)) {
  # create loop to run range exp
  snp.file <- re[i]
  coord.file = gsub(".csv", ".coords.csv", snp.file)
  
  sp = gsub("^.*\\/", "", snp.file)
  sp = gsub(".csv", "", sp)
  cat(sp, "\n")
  
  # curinds = read.csv(snp.file, stringsAsFactors = F, header = F)$V1
  # latlon = ind[match(curinds, ind$SAMPLE_ID), c("LAT", "LON")]
  # names(latlon) = c("latitude", "longitude")
  # outgroup = cl[match(curinds, cl$sample), "LatinName"]
  # outgroup = ifelse(outgroup == sp, 0, 1)
  # coord = data.frame(id = curinds, latlon, region = "REGION_1", outgroup)
  # write.csv(coord, coord.file, row.names = F, quote = F)
  
  region <- list("REGION_1")
  # #diploid individuals
  ploidy <- 2 
  raw.data <- load.data.snapp(snp.file, coord.file, sep=',', 
                               ploidy=ploidy)
  pop <- make.pop(raw.data, ploidy)
  psi <- get.all.psi(pop)
  avgpsi = mean(lower.tri(psi), na.rm = T)
  results <- run.regions(region = region, pop=pop, psi = psi, xlen=10, ylen=20)
   
  res[[i]] = list(results, avgpsi)
  # # summary.origin.results(results$tbl[[1]])
  # # plot(results)
}
saveRDS(res, "~/Dropbox/center_marginal/data/rangeExpansion.rds")
  

# res = readRDS("~/Dropbox/center_marginal/data/rangeExpansion.rds") 
res2 = lapply(res, function(x) { summary.origin.results(x[[1]]$tbl[[1]]) })
res3 = do.call("rbind", res2)
sps = gsub("^.*\\/", "", re)
sps = gsub(".csv", "", sps)
res3$species = sps

a = read.csv("~/Dropbox/center_marginal/data/dadi_correlation_centerdist.csv",
             stringsAsFactors = F)
a2 = full_join(res3, a, by = c("species" = "species"))
  

ranges = vector("list", length = nrow(a2))
for (i in 1:nrow(a2)) {
  sp_new = a2[i, "species"]
  rangedir = "/Users/sonal/Dropbox/publications/Sphenomorphine_Gene_Flow/data/geographic_ranges/OTU_ranges/"
  rangefile = paste(rangedir, sp_new, ".shp", sep="")
  if (file.exists(rangefile)) {
    range = readShapePoly(rangefile,
                          proj4string=CRS('+proj=longlat +datum=WGS84'))
    centroid = gCentroid(range)
    
    rangesize = sum(areaPolygon(range)) / 1e6 # sq km
    origin = a2[i, c("longitude", "latitude")]
    distance1 = distHaversine(origin, centroid) / 1000 # km
    distance2 = geosphere::dist2Line(p = origin, line = range)[, 1] / 1000
    
    ranges[[i]] = c(rangesize, distance1, distance2)
  }
}

geo = as.data.frame(do.call(rbind, ranges))
names(geo) = c("range_size", "distance1", "distance2")
geo$ratio = geo$distance1 / (geo$distance1 + geo$distance2)

a3 = cbind(a2, geo)
a3$range_sig = ifelse(a3$pval.x < 0.05, TRUE, FALSE)
a3 = a3 %>% filter(nind > 9) 
write.csv(a3, "~/Dropbox/center_marginal/data/demography.csv",
          quote = F, row.names = F)

a3 = read.csv("~/Dropbox (Personal)/center_marginal/data/demography.csv",
              stringsAsFactors = F)

######################
# any connection between pop growth and range exp?
#####################
a = ggplot(a3, aes(param1, rsq)) +
  geom_point() + xlab("ratio of curr. Ne to anc. Ne") +
  ylab(expression(r^2 ~ " of range exp."))
save_plot("~/Dropbox/publications/Center_Marginal/figures/rangeexp_popchange.pdf",
          a, base_width = 4)
cor.test(a3$param1, a3$rsq, method = "spearman")

######################
# any connection between div corr and range exp?
#####################
b = ggplot(a3, aes(abs(rho), rsq)) +
  geom_point() + xlab(expression("abs(center distance" ~ rho ~ ")")) +
  ylab(expression(r^2 ~ " of range exp."))
save_plot("~/Dropbox/publications/Center_Marginal/figures/rangeexp_divcorr.pdf",
          b, base_width = 4)
cor.test(abs(a3$rho), a3$rsq, method = "spearman")


######################
# any connection between div corr and range exp?
#####################

a4 = a3 %>% select(cor_sig, range_sig) %>% 
  group_by(cor_sig, range_sig) %>% summarize(count = n()) %>% ungroup()
c = ggplot(a4, aes(fill=range_sig, y=count, x=cor_sig)) + 
  geom_bar(position="fill", stat="identity") +
  scale_fill_manual(values = c("#ca0020", "#0571b0"),
                    name = expression("range exp. " ~ r^2)) +
  xlab(expression("distance-diversity" ~ rho))
save_plot("~/Dropbox/publications/Center_Marginal/figures/rangeexp_divcorr2.pdf",
          c, base_width = 4)

######################
# location of origin
#####################

d = ggplot(a3 %>% filter(range_sig == TRUE), aes(ratio)) + geom_histogram() + 
  xlab("relative center-edge distance of origin") 
save_plot("~/Dropbox (Personal)/publications/Center_Marginal/figures/rangeexp_origin.png",
          d, base_width = 4)

######################
# basic range exp
#####################

e = ggplot(a3, aes(rsq)) + 
  geom_histogram(aes(fill = range_sig), bins = 10) +
  xlab(expression("range exp. " ~ r^2)) + 
  theme(legend.title=element_blank())
save_plot("~/Dropbox/publications/Center_Marginal/figures/rangeexp.pdf",
          e, base_width = 4)

