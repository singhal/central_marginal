library(ape)
library(phytools)
library(glmulti)
library(dplyr)
library(raster)
library(tidyr)
library(maptools)
library(ggplot2)
library(phylolm)
library(sp)
library(sf)

setwd("~/Dropbox (Personal)/research_projects-FINISHED/Singhal_etal_CentralMarginal/")
raw = readRDS("data/correlation.Rds")
d = read.csv("data/correlation_centerdist.csv",
             stringsAsFactors = F)
names(d) = c("species", "pval", "rho", "nind")

pi = unlist(lapply(raw, function(x) { mean(x$pi, na.rm = T) }))
d$mean_pi = pi[d$species]

f = read.csv("~/Dropbox (Personal)/research_projects-FINISHED/Singhal_etal_SphenoGeneFlow/data/divergence/divergence.all_OTUs.csv", stringsAsFactors = F)
f1 = f %>% filter(distance_type == "geo_dist", 
             genetic_type == "fst") %>% 
             dplyr::select(species = OTU, geo_slope = slope)
f2 = f %>% filter(distance_type == "env_dist", 
                  genetic_type == "fst") %>% 
  dplyr::select(species = OTU, env_slope = slope)
d1 = left_join(d, f1)
d2 = left_join(d1, f2)

n = read.csv("~/Dropbox (Personal)/research_projects-FINISHED/Singhal_etal_SphenoGeneFlow/data/envirospatial/range_data.csv", stringsAsFactors = F)
n1 = n %>% dplyr::select(species = OTU, range_size, lat_midpoint, PC1_range, PC2_range)
d3 = left_join(d2, n1)

# not sure to do biomes bc
# mostly desert
b = read.csv("data/biomes.csv", stringsAsFactors = F)
d4 = left_join(d3, b)

# envs = vector("list", length = nrow(d3))
# r = c("aus5min_bio1.tif", "aus5min_bio12.tif")
# rdir = "~/Dropbox/macroevolution/eco_IBD_oz/data/AUS_5arc/"
# r1 = stack(paste(rdir, r, sep = ""))
# for (i in 1:nrow(d3)) {
#   sp_new = d3[i, "species"]
#   rangedir = "/Users/sonal/Dropbox/publications/Sphenomorphine_Gene_Flow/data/geographic_ranges/OTU_ranges/"
#   rangefile = paste(rangedir, sp_new, ".shp", sep="")
#   range = readShapePoly(rangefile,
#                           proj4string=CRS('+proj=longlat +datum=WGS84'))
#   sppts = spsample(range, 500, type = "random")
#   envs[[i]] = apply(raster::extract(r1, sppts), 2, mean, na.rm = 2)
# }
# d4 = cbind(d3, as.data.frame(do.call(rbind, envs)))

xx = read.csv("error/demography.csv",
              stringsAsFactors = F)
xx1 = xx %>% dplyr::select(species, popchange = param1, 
                     range_p = pval.x, founder = rsq)
d5 = left_join(d4, xx1)

# shape
s = read.csv("data/shape.csv")
d5 = left_join(d5, s %>% dplyr::select(species, dist_cv, cov, ex, ex2))

vars = c("mean_pi", "nind", "biome",
         "range_size", "geo_slope",
         "popchange", "founder",
         "cov", "dist_cv")

##########
# need to log anything?
###########
d4g = d5 %>% dplyr::select(vars)  %>% gather(metric, value) %>%
  filter(metric != "biome")
ggplot(d4g, aes(value)) + 
  geom_histogram(bins = 10) + 
  facet_wrap(~metric, scales = "free")
to_log = c("nind", "cov", 
           "PC1_range", "range_size", 
           "geo_slope", "founder")
d4l = d5 %>% mutate_at(to_log, log)
d4l$sig = ifelse(d4l$pval < 0.05, 1, 0)
d4l$posneg = ifelse(d4l$rho > 0, 1, 0)
d4l$sig_re = ifelse(d4l$range_p < 0.05, 1, 0)

phy = read.tree("~/Dropbox (Personal)/research_projects-FINISHED/Singhal_etal_SphenoGeneFlow/data/phylogeny/beast_ucln_31July17/otu.tre")
phy2 = drop.tip(phy, setdiff(phy$tip.label, d4$species))

# no phylogenetic signal
phylosig(phy2, 
         d4l[match(phy2$tip.label, d4l$species), "rho"],
         method = "lambda", test = T)
phylosig(phy2, 
         d4l[match(phy2$tip.label, d4l$species), "posneg"],
         method = "lambda", test = T)
phylosig(phy2, 
         d4l[match(phy2$tip.label, d4l$species), "sig"],
         method = "lambda", test = T)

## from https://github.com/mrhelmus/phylogeny_manipulation/blob/master/AICc.r
AICc.phylolm<-function(mod, return.K = FALSE, second.ord = TRUE, nobs = NULL, ...){
  
  if(identical(nobs, NULL)) {n <- length(mod$residuals)} else {n <- nobs}
  LL <- logLik(mod)$logLik
  K <- logLik(mod)$df  #extract correct number of parameters included in model - this includes LM
  if(second.ord == TRUE) {AICc <- -2*LL+2*K*(n/(n-K-1))}  else{AICc <- -2*LL+2*K}
  if(return.K == TRUE) AICc[1] <- K #attributes the first element of AICc to K
  return(AICc)
}




fit_models <- function(x, pi_val, variables, models, tree) {
  # fit all additive models and store results
  fits = vector("list", length=length(models))
  for (i in 1:length(models)) {
    # cat(i, "\n")
    tmp1 = x[complete.cases(x[, models[[i]]]), ]
    # rownames(tmp1) = tmp1$LatinName
    
    drop = setdiff(tree$tip.label, rownames(tmp1))
    tmp_tree = drop.tip(tree, drop)
    tmp1 = tmp1[tmp_tree$tip.label, ]
    
    fmla <- as.formula(paste(pi_val, " ~ ", paste(models[[i]], collapse= "+")))
    fit = phylolm(fmla, data=tmp1, tmp_tree, model="lambda", lower.bound=0)
    
    fits[[i]] = fit
    if (i %% 100 == 0) {
      cat("Model ", i, " of ", length(models), " done.\n", sep="")
    }
  }
  
  aics = sapply(fits, AICc.phylolm)
  best = min(aics)
  raw_weights = sapply(aics, weights <- function(x) {exp((best - x)/2)})
  weights = raw_weights / sum(raw_weights)
  
  results = vector("list", length=length(variables))
  names(results) = variables
  
  for (i in 1:length(variables)) {
    var_models = sapply(models, y <- function(x) {variables[i] %in% x})
    
    var_fits = fits[var_models]
    var_coef = sapply(var_fits, y <- function(x) {x$coefficients[variables[i]]})
    # var_pval = sapply(var_fits, y <- function(x) {summary(x)$coefficients[variables[i],'p.value']})
    var_wt = weights[var_models]
    var_raw_wt = raw_weights[var_models]
    var_rel_raw_wt = var_raw_wt / sum(var_raw_wt)
    
    coef = sum(var_coef * var_rel_raw_wt)
    # pval = sum(var_pval * var_rel_raw_wt)
    rel_imp = sum(var_wt)
    
    vals = c(coef, rel_imp)
    names(vals) = c("coefficient", "rel_importance")
    results[[i]] = vals
  }
  
  full_results = list(fits, aics, weights, results)
  names(full_results) = c("fits", "AIC", "weight", "results")
  return(full_results)
}


# create all possible models
# this is so many models
models = list()
for (i in 1:length(vars)) {
  l = combn(vars, i, simplify=FALSE)
  models = append(models, l)
}

# only use cases fit across all models
rownames(d4l) = d4l$species
all_fit1 = fit_models(d4l, "rho", vars, models, phy2)

full_fit1 = all_fit1[["results"]]
res = as.data.frame(t(as.data.frame(full_fit1)))
# best models
models[order(all_fit1[["weight"]], decreasing = T)[1:5]]

bestmod = phylolm("rho ~ cov", data=d4l, phy2, model="lambda", lower.bound=0)
summary(bestmod)

pdf("error/model_testing1.pdf",
    height = 3, width = 5)
par(mar=c(3.5, 1.1, 3.1, 4.1))
long = c("mean genetic diversity", "# of inds", "biome",
        "range size", "IBD slope", 
         "pop. size change", "range expansion",
        "sampling coverage", "range eccentricity")
names(long) = vars
midpoints <- barplot(res$rel_importance, horiz=T)
text(0.01, midpoints, labels=long[rownames(res)], adj=c(0,0.5), cex=0.8)
mtext("relative importance", side=1, line=2)
mtext(side=3, at= 0.98, text="coef", font=3, las= 1, 
      line=1, adj=c(0.5,0.5), cex=0.8)
mtext(side=4, at=midpoints, 
      text=ifelse(res$coefficient > 0, "+", "-"),
      las=2, line=2.5, cex=1.2, adj=c(0.5,0.5))
dev.off()

a = ggplot(d4l) + geom_point(aes(founder, rho)) +
  xlab("Log strength of range expansion") +
  ylab("corr. of distance-diversity")
save_plot("~/Dropbox (Personal)/publications/Center_Marginal/figures/model_testing2.pdf",
          a,
          base_height = 3, base_width = 4)
