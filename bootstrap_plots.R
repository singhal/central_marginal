library(ggplot2)
library(cowplot)
library(RColorBrewer)
library(patchwork)
theme_set(theme_cowplot())

boot = list.files("~/Dropbox/center_marginal/data/bootstrap/",
                  full.names = T)

# make one representative plot
# Lerista_bipes_2
d = read.csv(boot[grep("Lerista_bipes_1", boot)], stringsAsFactors = F)
cts = d %>% group_by(individual) %>% 
  filter(complete.cases(pi)) %>% 
  summarize(locict = length(unique(num_loci)))
remove = cts %>% filter(locict < 2) %>% pull(individual)
d1 = d %>% filter(!individual %in% remove) %>% filter(num_loci < 20000)
plotA = ggplot(d1, aes(num_loci, pi)) + 
  geom_point(alpha = 0.5)  +
  facet_wrap(~individual, ncol = 3) + xlab("number of sampled loci") +
  ylab(expression(pi))

# make summary plot
process_boot <- function(x) {
  d = read.csv(x, stringsAsFactors = F)
  d1 = d %>% group_by(individual, num_loci) %>% 
    summarize(var_pi = var(pi), denom = mean(denom)) %>% ungroup()
  return(d1)
}
boot2 = lapply(boot, process_boot)
boot3 = bind_rows(boot2, .id = "column_label") %>% filter(num_loci < 20e3)

plotB = ggplot(boot3, aes(num_loci, var_pi)) + 
  geom_point(alpha = 0.05) + xlab("number of sampled loci") +
  ylab(expression('variance in ' ~pi))

ab = plotA + plotB+ plot_layout(widths = c(1.5, 1))
save_plot("~/Dropbox/publications/Center_Marginal/figures/bootstrap.pdf", 
          ab, base_width = 14)
# look at denom across
boot3 %>% group_by(num_loci) %>% 
  summarize(mean(denom, na.rm = T))

