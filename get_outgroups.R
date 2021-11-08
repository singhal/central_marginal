library(ape)
library(phangorn)
library(phytools)

# first get species
res = readRDS("~/Dropbox/center_marginal/data/correlation.Rds")
sps = unlist(lapply(res, function(x) { x$species[1] }))

# then get phylogeny
t = read.tree("~/Dropbox/Sphenomorphine_Phylogeny/data/concat/ExaML_result.concat_ind0.05_loci0.05_all_n5353_2.rooted.tre")
d = read.csv("~/Dropbox/Sphenomorphine_Phylogeny/sample_data/sphenomorphine_samples_v2.csv", stringsAsFactors = F)
# some tips belong to the same lineage
dups =  d[match(t$tip.label, d$sample), "lineage"]
t = drop.tip(t, t$tip.label[duplicated(dups)])
t$tip.label = d[match(t$tip.label, d$sample), "lineage"]

# now get the reads
# ind = dbGetQuery(con, "select * from individuals")
# r = dbGetQuery(con, "select * from ddrad")

# counts
inds2 = ind %>% filter(complete.cases(OTU)) %>% pull(SAMPLE_ID)
inds3 = inds2[inds2 %in% r$SAMPLE_ID]
cts = table(ind %>% filter(SAMPLE_ID %in% inds3) %>% pull(OTU))

# then get outgroups and reads
outs = rep(NA, length(sps))
reads = vector("list", length(sps))

for (i in 1:length(sps)) {
  sp = sps[i]
  pnode = Ancestors(t, which(t$tip.label == sp), type = "parent")
  outs1 = t$tip.label[ Descendants(t, pnode, type = "tips")[[1]] ]
  outs1 = outs1[outs1 != sp]
  
  # if (sum(outs1 %in% outs) > 0) {
  #  outs2 = outs1[outs1 %in% outs]
  #  outs2 = cts[outs2]
  #  outs[i] = names(sort(outs2, decreasing = T))[1]
  #} else {
    outs2 = cts[outs1]
    outs2 = outs2[outs2 < 50]
    outs[i] = names(sort(outs2, decreasing = T))[1]
  # }
  ind2 = ind %>% filter(OTU == outs[i]) %>% select(SAMPLE_ID)
  reads[[i]] = r[r$SAMPLE_ID %in% ind2$SAMPLE_ID, c("SAMPLE_ID", "READ1", "READ2")]
  reads[[i]]$OTU = outs[i]
}

length(unique(outs))
reads2 = do.call("rbind", reads)
reads3 = reads2[!duplicated(reads2), ]
dim(reads3)

allreads = data.frame(reads = c(reads3$READ1, reads3$READ2))
write.csv(allreads, "~/Desktop/reads.txt", row.names = F, quote = FALSE)

# add target PRG
cl = read.csv("~/Dropbox/macroevolution/eco_IBD_oz/data/clustering/Ctenotus_Lerista.clustering.revised3.csv", stringsAsFactors = F)
cl$LatinName = gsub("C\\. ", "Ctenotus_", cl$LatinName)
cl$LatinName = gsub("L\\. ", "Lerista_", cl$LatinName)
cl$LatinName = gsub(" ", "_", cl$LatinName)

# add target species
for (i in 1:length(reads)) {
  reads[[i]]$REF_SPECIES = sps[i] 
  lin = unique(cl[which(cl$LatinName == sps[i]), "GMYC_RAxML2"])
  if (length(lin) > 1) {
    print(sps[i])
  } else {
    reads[[i]]$LINEAGE = lin
  }
}
reads4 = do.call("rbind", reads)
reads4$READ1 = gsub("/nfs/turbo/lsa-rabosky/Lab/skink_ddrad/individual_reads/",
                "/home/user/Desktop/center_marginal/raw_reads/", reads4$READ1)
reads4$READ2 = gsub("/nfs/turbo/lsa-rabosky/Lab/skink_ddrad/individual_reads/",
                    "/home/user/Desktop/center_marginal/raw_reads/", reads4$READ2)
write.csv(reads4, "~/Desktop/cm_samples.csv", quote = FALSE, row.names = FALSE)

