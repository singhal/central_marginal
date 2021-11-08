################################################
# test if pi changes as we approach the edges  #
# as we would expect                           #
################################################

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
  
  rangedir = "~/Dropbox (Personal)//publications/Sphenomorphine_Gene_Flow/data/geographic_ranges/OTU_ranges/"
  rangefile = paste(rangedir, sp_new, ".shp", sep="")

  
  if (file.exists(rangefile)) {
    range = readShapePoly(rangefile,
                          proj4string=CRS('+proj=longlat +datum=WGS84'))
    points = SpatialPoints(tmp4[,c("LON", "LAT")],
                           proj4string=CRS('+proj=longlat +datum=WGS84'))

    rangesize = sum(areaPolygon(range)) / 1e6

    # distance from center
    centroid = gCentroid(range)
    distance1 = distHaversine(points, centroid)

    # distance from edge
    # need to double check for species with multiple polygons
    distance2 = geosphere::dist2Line(p = points, line = range)[, 1]

    # ratio
    distance3 = distance1 / (distance1 + distance2)

    tmp5 = cbind(tmp4, distance1, distance2, distance3)  
    tmp5$species = sp_new 
      
    # combine env dist
    if (!is.null(PCenv[[sp_new]])) {
      tmp6 = left_join(tmp5, PCenv[[sp_new]], by = c("individual" = "ind"))
      tmp7 = left_join(tmp6, MVEenv[[sp_new]][[1]], by = c("individual" = "ind"))
    } else {
      tmp7 = tmp5
      tmp7$PC_eco_dist = NA
      tmp7$prop6 = NA
      tmp7$MHD = NA
    }
  } else {
    tmp7 = NULL
  }
  return(list(tmp7, sp_new))
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

count_inds <- function(file) {
  f = read.csv(file, stringsAsFactors = F)
  return(nrow(f))
}

files = list.files("~/Dropbox (Personal)//center_marginal/data/diversity/", pattern="*.csv", full.names=T)
ninds = unlist(lapply(files, count_inds))
files2 = files[ninds >= 9]

# combine these two data frames
ind1 = read.csv('~/Dropbox (Personal)//macroevolution/eco_IBD_oz/data/metadata/individual_data_nomissing22Oct15.csv')
ind1a = ind1[, c("sample_id", "lat", "lon")]
names(ind1a) = c("SAMPLE_ID", "LAT", "LON")
ind2 = read.csv("~/Dropbox (Personal)//publications/Sphenomorphine_Gene_Flow/data/metadata/spheno_ind_data.csv", stringsAsFactors = F)
ind1b = ind1a[! ind1a$SAMPLE_ID %in% ind2$SAMPLE_ID, ]
ind = rbind(ind1b, ind2[, c("SAMPLE_ID", "LAT", "LON")])

# read in environmental distance data
PCenv = readRDS("~/Dropbox (Personal)/center_marginal/data/spatial/PCres.rds")
MVEenv = readRDS("~/Dropbox (Personal)/center_marginal/data/spatial/MVEres.rds")

res = lapply(files2, calculate_diversity)

###############
# main plots
###############
get_corr <- function(df, dist) {
  df2 = df[, c("pi", dist)]
  df3 = df2[complete.cases(df2), ]
  if (nrow(df3) > 0) {
    corres = cor.test(as.numeric(df3$pi), as.numeric(df3[, dist]),
                      method = "spearman")
    pval = corres$p.value
    est = corres$estimate
    num = nrow(df3)
  } else {
    pval = NA
    est = NA
    num = 0
  }
  return(c(pval, est, num))
}

resa = lapply(res, function(x) {x[[1]]})
names(resa) = unlist(lapply(res, function(x) {x[[2]]}))
# saveRDS(resa, "~/Dropbox (Personal)/center_marginal/data/correlation.Rds")
resa = readRDS("~/Dropbox (Personal)/center_marginal/data/correlation.Rds")
res2 = resa[!unlist(lapply(resa, is.null))]

#########################
# main plot
#########################

cor1 = as.data.frame(do.call(rbind, lapply(res2, get_corr, "distance1")))
cor1 = cor1[which(cor1$V3 > 9 & complete.cases(cor1$rho)), ]
write.csv(cor1, "~/Dropbox (Personal)/center_marginal/data/correlation_centerdist.csv",
          quote = F)

cor2 = as.data.frame(do.call(rbind, lapply(res2, get_corr, "PC_eco_dist")))
cor2 = cor2[which(cor2$V3 > 9 & complete.cases(cor2$rho)), ]

plotA = ggplot(cor1, aes(rho, fill = ifelse(V1 < 0.05, "sig", "non-sig"))) + 
  geom_histogram(bins = 10, boundary = 0) + guides(fill=guide_legend(title="")) +
  scale_fill_manual(values = c("gray80", "black")) +
  xlab(expression("center distance" ~ rho))
plotB = ggplot(cor2, aes(rho, fill = ifelse(V1 < 0.05, "sig", "non-sig"))) + 
  geom_histogram(bins = 10, boundary = 0) + guides(fill=guide_legend(title="")) +
  scale_fill_manual(values = c("gray80", "black")) +
  xlab(expression("PC climate distance " ~ rho))

prow <- plot_grid(
  plotA + theme(legend.position="none"),
  plotB + theme(legend.position="none"),
  align = 'vh',
  labels = c("A", "B"),
  hjust = -1,
  nrow = 1
)
legend <- get_legend(
  # create some space to the left of the legend
  plotA + theme(legend.box.margin = margin(0, 0, 0, 4))
)

# add the legend to the row we made earlier. Give it one-third of 
# the width of one plot (via rel_widths).
ab = plot_grid(prow, legend, rel_widths = c(3, .4))
save_plot("~/Dropbox (Personal)/publications/Center_Marginal/figures/geo_eco_dist.png",
          ab, ncol = 2, base_width = 5)

############
# plots for correlations for alternate forms of distance
###########

cor3 = as.data.frame(do.call(rbind, lapply(res2, get_corr, "distance2")))
cor3 = cor3[which(cor3$V3 > 9 & complete.cases(cor3$rho)), ]

cor4 = as.data.frame(do.call(rbind, lapply(res2, get_corr, "distance3")))
cor4 = cor4[which(cor4$V3 > 9 & complete.cases(cor4$rho)), ]

cor5 = as.data.frame(do.call(rbind, lapply(res2, get_corr, "MHD")))
cor5 = cor5[which(cor5$V3 > 9 & complete.cases(cor5$rho)), ]

plotA = ggplot(cor3, aes(rho, 
                 fill = ifelse(V1 < 0.05, "sig", "non-sig"))) + 
  geom_histogram(bins = 10, boundary = 0) + guides(fill=guide_legend(title="")) +
  scale_fill_manual(values = c("gray80", "black")) +
  xlab(expression("edge distance" ~ rho))
plotB = ggplot(cor4, aes(rho, 
                         fill = ifelse(V1 < 0.05, "sig", "non-sig"))) + 
  geom_histogram(bins = 10, boundary = 0) + guides(fill=guide_legend(title="")) +
  scale_fill_manual(values = c("gray80", "black")) +
  xlab(expression("center-edge distance ratio" ~ rho))

plotC = ggplot(cor5, aes(rho, 
                         fill = ifelse(V1 < 0.05, "sig", "non-sig"))) + 
  geom_histogram(bins = 10, boundary = 0) + guides(fill=guide_legend(title="")) +
  scale_fill_manual(values = c("gray80", "black")) +
  xlab(expression("MVE distance " ~ rho))

prow <- plot_grid(
  plotA + theme(legend.position="none"),
  plotB + theme(legend.position="none"),
  plotC + theme(legend.position="none"),
  align = 'vh',
  labels = c("A", "B", "C"),
  hjust = -1,
  nrow = 1
)
legend <- get_legend(
  # create some space to the left of the legend
  plotA + theme(legend.box.margin = margin(0, 0, 0, 4))
)

# add the legend to the row we made earlier. 
abc = plot_grid(prow, legend, rel_widths = c(3, .3))
save_plot("~/Dropbox (Personal)/publications/Center_Marginal/figures/alternate_distances1.png",
          abc, ncol = 3, base_width = 5)

cortable = cbind(cor1, cor2, cor3, cor4, cor5)
for (i in 1:ncol(cortable)) {
  cortable[,i] = round(cortable[,i], 2) 
}
write.csv(cortable, "~/Dropbox (Personal)/publications/Center_Marginal/figures/cortable.csv",
          row.names = T, quote = F)

###############
# how good are the eco models
##############

# how much do the 6 PC axes summarize?
pc6 = unlist(lapply(PCenv, function(x) {return(mean(x$prop6, na.rm = T))}))
mean(pc6)
# 91% 

# what is the omission rate
get_omr <- function(x, val) {
  best = x[[2]][1, ]
  return(c(best[val]))
}

om_train = unlist(lapply(MVEenv, get_omr, 'om_rate_train'))
mean(om_train)
#  0.05027714

om_test = unlist(lapply(MVEenv, get_omr, 'om_rate_test'))
mean(om_test)
# 0.04289997

bg_prev = unlist(lapply(MVEenv, get_omr, 'bg_prevalence'))
mean(bg_prev)
# 0.3365536

###############
# correlations of alternate distance
###############

res3 = do.call("rbind", res2)
res3$distance1 = res3$distance1 / 1000
res3$distance2 = res3$distance2 / 1000

dists = c("distance1", "distance2", "distance3", "PC_eco_dist", "MHD")
for (i in 1:length(dists)) {
  for (j in (i+1):length(dists)) {
    xx = cor.test(log(res3[ , dists[i]]), 
                  log(res3[ , dists[j]]), 
                  method = "pearson")
    cat(dists[i], dists[j], round(xx$estimate, 3), "\n")
  }
}

a = ggplot(res3, aes(distance1, distance2)) + 
  geom_point(alpha = 0.5) + xlab("center distance (km)") +
  ylab("edge distance (km)")

b = ggplot(res3, aes(distance1, distance3)) + 
  geom_point(alpha = 0.5) + xlab("center distance (km)") +
  ylab("center-edge ratio")

c = ggplot(res3, aes(distance1, PC_eco_dist)) + 
  geom_point(alpha = 0.5) + xlab("center distance (km)") +
  ylab("PC clim. dist.")

d = ggplot(res3, aes(PC_eco_dist, MHD)) + 
  geom_point(alpha = 0.5) + xlab("PC clim. dist.") +
  ylab("MVE clim. dist.") + ylim(0, 20)



prow <- plot_grid(a, b, c, d, align = 'vh',
  labels = c("A", "B", "C", "D"),
  hjust = -1,
  nrow = 2
)
save_plot("~/Dropbox (Personal)/publications/Center_Marginal/figures/distance_correlations.png",
          prow, ncol = 2, nrow = 2, base_width = 5, base_height = 4)


 ##################
# combined approach - linear model
##################

lmp <- function (modelobject) {
  if (class(modelobject) != "lm") stop("Not an object of class 'lm' ")
  f <- summary(modelobject)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(p)
}

do_lm <- function(x) {
  x1 = x[ , c("distance1", "PC_eco_dist", "pi")]
  x1 = x1[complete.cases(x1), ]
  if (nrow(x1) > 9) {
    xx = lm(x1$pi ~ log(x1$distance1) + log(x1$PC_eco_dist))
    
    # coefficients for distance, PC_eco_dist, sig, sig for var, r2
    coefxx = coefficients(xx)[2:3]
    pval = lmp(xx)
    sigvar = summary(xx)$coefficients[2:3 , 4]
    r2 = summary(xx)$adj.r.squared
    
    xx1 = lm(x1$pi ~ log(x1$distance1))
    coefxx1 = coefficients(xx1)[2]
    pval1 = lmp(xx1)
    r2_1 = summary(xx1)$adj.r.squared
    
    xx2 = lm(x1$pi ~ log(x1$PC_eco_dist))
    coefxx2 = coefficients(xx2)[2]
    pval2 = lmp(xx2)
    r2_2 = summary(xx2)$adj.r.squared
    
    res = c(coefxx, pval, sigvar, r2, coefxx1, pval2, r2_1,
            coefxx1, pval2, r2_2)
  } else {
    res = rep(NA, 12)
  }
  names(res) = c("geo_coef", "eco_coef", 
                 "pval", "geo_sig",
                 "eco_sig", "r2", 
                 "geo_only_coef", "geo_only_pval", "geo_only_r2",
                 "eco_only_coef", "eco_only_pval", "eco_only_r2")
  return(res)
}

lmres = lapply(res2, do_lm)
lmres2 = do.call("rbind", lmres)
lmres2 = as.data.frame(lmres2)

plotA = ggplot(lmres2, aes(geo_coef, 
                         fill = ifelse(geo_sig < 0.05, "sig", "non-sig"))) + 
  geom_histogram(bins = 10) + guides(fill=guide_legend(title="")) +
  scale_fill_manual(values = c("gray80", "black")) +
  xlab(expression("geo. dist. coeff."))

plotB = ggplot(lmres2, aes(eco_coef, 
                           fill = ifelse(eco_sig < 0.05, "sig", "non-sig"))) + 
  geom_histogram(bins = 10) + guides(fill=guide_legend(title="")) +
  scale_fill_manual(values = c("gray80", "black")) +
  xlab(expression("clim. dist. coeff."))

prow <- plot_grid(
  plotA + theme(legend.position="none"),
  plotB + theme(legend.position="none"),
  align = 'vh',
  labels = c("A", "B"),
  hjust = -1,
  nrow = 1
)
legend <- get_legend(
  # create some space to the left of the legend
  plotA + theme(legend.box.margin = margin(0, 0, 0, 4))
)

# add the legend to the row we made earlier. Give it one-third of 
# the width of one plot (via rel_widths).
ab = plot_grid(prow, legend, rel_widths = c(3, .4))
save_plot("~/Dropbox (Personal)/publications/Center_Marginal/figures/geo_eco_lm.png",
          ab, ncol = 2, base_width = 5)


############
# plots for correlations of correlations
###########

res4 = as.data.frame(cbind(cor1$rho, cor2$rho, cor3$rho, cor4$rho, cor5$rho))
names(res4) = c("CD_dist_rho", "PC_eco_rho", 
                "ED_dist_rho", "CDED_ratio_rho", "MVE_dist_rho")

a = ggplot(res4, aes(CD_dist_rho, ED_dist_rho)) + 
  geom_point() + xlab(expression("center distance" ~ rho)) +
  ylab(expression("edge distance" ~ rho))
b = ggplot(res4, aes(CD_dist_rho, CDED_ratio_rho)) + 
  geom_point() + xlab(expression("center distance" ~ rho)) +
  ylab(expression("center-edge ratio" ~ rho))
c = ggplot(res4, aes(CD_dist_rho, PC_eco_rho)) + 
  geom_point() + xlab(expression("center distance" ~ rho)) +
  ylab(expression("PC clim distance" ~ rho))
d = ggplot(res4, aes(PC_eco_rho, MVE_dist_rho)) + 
  geom_point() + xlab(expression("PC eco." ~ rho)) +
  ylab(expression("MVE distance" ~ rho))

prow <- plot_grid(a, b, c, d, align = 'vh',
                  labels = c("A", "B", "C", "D"),
                  hjust = -1,
                  nrow = 2
)
save_plot("~/Dropbox (Personal)/publications/Center_Marginal/figures/corr_of_correlations.pdf",
          prow, ncol = 2, nrow = 2, base_width = 5, base_height = 4)

######################
# make individual maps
######################

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
sp = readxl::read_xlsx("~/Desktop/species_names.xlsx")
sp2 = pull(sp[match(sps, sp$OTU), "sp_name"])

pdf("~/Dropbox (Personal)/publications/Center_Marginal/figures/species_ranges.pdf", height = 10, width = 10)
par(mfrow = c(5, 5))
for (i in 1:length(sps)) {
  if (nrow(maps[[i]][[2]]) > 9) {
    par(mar = c(0, 0, 0, 0))
    plot(1, type="n", axes=FALSE, ann=FALSE,
         xlim = c(110, 155), ylim = c(-46, -8.5), 
         xlab = "", ylab = "")
    plot(oz[1], col = "gray95", add = T)
    mtext(sp2[i], side = 3, font = 3, cex = 0.7)
    plot(maps[[i]][[1]], col = "forestgreen", border = FALSE, add = T)
    points(maps[[i]][[2]][, c("LON", "LAT")], pch = 21, bg = "white")
  }
}
dev.off()

