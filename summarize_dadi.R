library(dplyr)
library(tidyr)
library(ggplot2)
library(cowplot)
theme_set(theme_cowplot())

d = list.files("/Users/Sonal/Dropbox/center_marginal/data/dadi_output/unfolded/",
               full.names = T)
sps = gsub(".*\\/", "", d)
res = vector('list', length(sps))
names(res) = sps
for (i in 1:length(d)) {
  f = paste0(d[i], "/Results_Summary_Short.txt")
  f = read.table(f, stringsAsFactors = F, header = T)
  res[[i]] = f
}

count_param <- function(col) {
  return(length(strsplit(col, ",")[[1]]))
}

summarize_tables <- function(xx) {
  xx$num_params = sapply(xx$optimized_params, count_param)
  
  xx1 = xx %>% filter(num_params < 3) %>%
    group_by(num_params) %>% 
    slice_max(order_by = log.likelihood, n = 1) %>% 
    ungroup() %>% arrange(num_params) %>% 
    mutate(stat = -2 * (lag(log.likelihood) - log.likelihood))
  xx1 = xx1 %>% mutate(pval = pchisq(stat, df = 1, lower.tail = FALSE))
  xx2 = xx1 %>% filter(pval < 0.05 | is.na(pval)) %>%
    arrange(-log.likelihood)

  # get best fitting model
  best = pull(xx2[1, "Model"])
  
  # get theta
  theta = pull(xx2[1, "theta"])
  
  # param estimates
  params = pull(xx2[1, "optimized_params"])
  params = as.numeric(strsplit(params, ",")[[1]])
  params = c(params, rep(NA, 3 - length(params)))
  
  return(c(best, theta, params))
}

res2 = lapply(res, summarize_tables)
res2 = as.data.frame(do.call(rbind, res2))
names(res2) = c("model", "theta", "param1", "param2", "param3")
res2 = res2 %>% 
  mutate_at(c("theta", "param1", "param2", "param3"), as.numeric)
# need to get L estimates
res2$species = sps

len = read.csv("~/Dropbox/center_marginal/data/dadi_output/ sequence_lengths.csv",
               stringsAsFactors = F)
res3 = inner_join(res2, len)
mu = 1e-8
res3$ne = res3$theta / (4 * mu * res3$seqlen)
res3$time = res3$ne * 2 * res3$param2
res3$pop = res3$ne * res3$param1

########################
# plot main pop diff
########################

res3$relative = ifelse(res2$param1 > 1, "greater", "less")

a = ggplot(res3, aes(param1)) + 
  geom_histogram(bins = 20) +
  geom_vline(xintercept = 1) +
  xlab("ratio of curr. to anc. Ne")
save_plot("~/Dropbox/publications/Center_Marginal/figures/dadi_popchange.pdf",
          a, base_width = 4)

cor1 = read.csv("~/Dropbox/center_marginal/data/correlation_centerdist.csv",
                stringsAsFactors = F, row.names = 1)
names(cor1) = c("pval", "rho", "nind")
cor1$species = rownames(cor1)
aa = left_join(cor1, res3)

aa$cor_sig = ifelse(aa$pval <= 0.05, "sig", "non-sig")
write.csv(aa, "~/Dropbox/center_marginal/data/dadi_correlation_centerdist.csv",
          quote = F, row.names = F)


b = ggplot(aa, aes(param1, rho)) + geom_point(aes(fill = cor_sig), 
                                          shape = 21) + 
  xlab("ratio of curr. to anc. Ne") + 
  ylab("distance-diversity correlation") + 
  scale_fill_manual(values = c("#ca0020", "#0571b0")) +
  theme(legend.title=element_blank())
save_plot("~/Dropbox/publications/Center_Marginal/figures/dadi_popchange_corr.pdf",
          b, base_width = 6)

cor.test(aa$param1, aa$rho, method = "spearman")
