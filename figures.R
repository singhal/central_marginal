library(ggplot2)
library(cowplot)
library(rgdal)
library(dplyr)
library(rworldmap)
library(patchwork)
theme_set(theme_cowplot())

setwd("~/Dropbox (Personal)/center_marginal/data/")

d = read.csv("correlation_centerdist.csv", stringsAsFactors = F)
x = readRDS("correlation.Rds")

# species
# Ctenotus_robustus_3
# Ctenotus_regius
# C euclae

worldmap <- getMap(resolution = "low")
oz = worldmap[which(worldmap[["SOVEREIGNT"]] == "Australia"), ]
oz2 = fortify(oz)

plts = vector("list", 6)
sps = c("Ctenotus_robustus_3", "Ctenotus_inornatus_1", "Ctenotus_euclae")
for (i in 1:length(sps)) {
  df = x[[sps[i]]]
  df2 = df %>% group_by(distance1) %>% summarize(pi = mean(pi))
  plts[[i + 3]] = ggplot(df2, aes(distance1 / 1000, pi)) + 
    geom_point(aes(colour = pi), pch = 16, size = 2) + 
    ylab(expression(pi)) + xlab("distance from center (km)") +
    scale_colour_viridis_c() + 
    theme(legend.position = "none")
  
  rangedir = "~/Dropbox (Personal)/publications/Sphenomorphine_Gene_Flow/data/geographic_ranges/OTU_ranges/"
  rangefile = paste(rangedir, sps[i], ".shp", sep="")
  range = readOGR(rangefile)
  range2 = fortify(range)
  
  plts[[i]] = ggplot() +
    coord_map(xlim = c(113, 154), ylim = c(-44, -12)) +
    geom_polygon(data = oz2, aes(long, lat, group = group),
                 color = "grey80", fill = "grey80", size = 0.3) +
    geom_polygon(data = range2, aes(long, lat, group = group),
                 color = "grey60", fill = "grey60", size = 0.3) +
    geom_point(data = df, aes(x = LON, y = LAT, fill = pi), 
               pch = 21, size = 2) + 
    scale_fill_viridis_c() +
    theme_void() +  theme(legend.position = "none")
}

plt = plot_grid(plts[[1]], plts[[2]], plts[[3]],
          plts[[4]], plts[[5]], plts[[6]], ncol = 3,
          rel_heights = c(1, 1), align = "vh")
save_plot("~/Desktop/test.pdf", plt, base_height = 2.5, 
          base_width = 4, ncol = 3, nrow = 2)


d = read.csv("demography.csv", stringsAsFactors = F)

a = ggplot(d, aes(param1)) + geom_histogram(bins = 10) +
  xlab("ratio of curr. Ne to anc. Ne") 
d2 = d %>% filter(range_sig == TRUE)
b = ggplot(d2, aes(rsq)) + geom_histogram(bins = 10) +
  xlab("strength of range expansion")
plt = plot_grid(a, b, ncol = 2, align = "vh", labels = c("A", "B"))
save_plot("~/Dropbox (Personal)/publications/Center_Marginal/figures/demography.png", 
          plt, base_height = 4, 
          base_width = 5, ncol = 2)
