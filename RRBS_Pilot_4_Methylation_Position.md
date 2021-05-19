---
title: "RRBS Pilot Methylation by Position"
author: "Annat Haber"
date: '2021-05-19'
output: 
  html_document:
    toc: true
    toc_depth: 3
    code_folding: hide
---


```r
knitr::opts_knit$set(root.dir = rprojroot::find_rstudio_root_file())
knitr::opts_chunk$set(message=FALSE, warning=FALSE, cache = FALSE, cache.path = "output/RRBS_Pilot/4_Methylation_Position", eval=FALSE)
```


```r
rm(list=ls())

suppressPackageStartupMessages({
  library(tidyverse)
  library(GenomicRanges)
  library(biomaRt)
  })

source("scripts/AnnaTyler/center.on.feature.R")
source("scripts/AnnaTyler/plot.centered.vals.R")
source("scripts/AnnaTyler/plot.poly.xy.R")
```

Load data

```r
# a matrix of methyl percent per site by sample; only sites that are shared by all samples
meth <- readRDS("output/2_Analysis/precentMethylation.RDS")

# metadata
meta <- read_csv("output/2_Analysis/sampleInfo.globalmeth.csv")

# reference annotations from http://genome.ucsc.edu/cgi-bin/hgTables
refseq <- read_tsv("data/mm10.refseq.genes.txt")
```

Gene info for all known mouse genes  
rows are genes; columns are external_gene_name, start_position, end_position, strand

```r
mart <- useMart(biomart="ensembl", dataset="mmusculus_gene_ensembl")
gene.info <- getBM(attributes = c("external_gene_name", "chromosome_name", "start_position", "end_position", "strand"), filters = c("external_gene_name"), values = unique(refseq$name2), mart = mart)
del <- which(is.na(as.numeric(gene.info$chromosome_name)))
gene.info <- gene.info[-del,]
```

Map snps to genes, including downstream and upstream buffer zones

```r
buffer <- 2000
ext.start <- gene.info$start_position - buffer
ext.end <- gene.info$end_position + buffer

snps <- rownames(meth) %>%
  strsplit(split=":") %>% unlist() %>% 
  matrix(nr=nrow(meth), nc=2, byrow = TRUE) %>%
  apply(2, as.numeric)
snps <- snps[!is.na(snps[,1]),] # removing X and Y

snps.genes <- c()
for (i in 1:nrow(snps)) {
  chr <- snps[i,1]
  pos <- snps[i,2]
  genes.i <- which(gene.info$chromosome_name == chr &
                 ext.start < pos &
                 ext.end > pos)
  if (length(genes.i)>0) {
    snps.genes <- rbind(snps.genes, 
                      cbind(pos, gene=gene.info$external_gene_name[genes.i]))
  }
}
```
Get methyl values by gene, averaged across replicates  
Coordinates are normalized within gene   
One vector for CC and one for TT for each gene

```r
gene.methyl.norm <- list()
for (tr in c("CC", "TT")) {
  for (g in unique(snps.genes[,"gene"])) {
    chr <- gene.info$chromosome_name[gene.info$external_gene_name==g]
    pos.g <-  paste(chr, snps.genes[snps.genes[,"gene"]==g, "pos"], sep=":")
    if (length(pos.g)>1) {
      # extracting and averaging methyl values across replicates
      # mg <- meth[pos.g, meta$genotype==tr & meta$age=="13"] %>% rowMeans()
      mg <- meth[pos.g, meta$genotype==tr] %>% rowMeans()
     names(mg) <- matrix(unlist(strsplit(names(mg), split=":")), nc=2, byrow=TRUE)[,2]
      # centering on TSS and normalizing by gene length
      gene.methyl.norm[[tr]][[g]] <- center.on.feature(
                                      gene.name = g, # mgi_symbol
                                      gene.info = gene.info, 
                                      vals = mg, # methyl vals of one gene with coordinates as names
                                      feature = "full")
    }
  }
}
```

Aligning and ploting

```r
# aligning the centered normallized methylation values for CC averaged across genes
cc.aligned.means <- plot.centered.vals(val.list = gene.methyl.norm[["CC"]], 
plot.label = "CC Genes", ylim = c(0, 500), min.representation = 1, seq.by = 0.01, merge.by = 0,
plot.individual = FALSE, return.means=TRUE, min.upstream = -2000, max.downstream = 2000)

# centered normallized methylation values for TT averaged across genes
tt.aligned.means <- plot.centered.vals(val.list = gene.methyl.norm[["TT"]], 
plot.label = "TT Genes", ylim = c(0, 500), min.representation =1, seq.by = 0.01, merge.by = 0,
plot.individual = FALSE, return.means=TRUE, min.upstream = -2000, max.downstream = 2000)

# position values
cc.pos <- as.numeric(names(cc.aligned.means))
tt.pos <- as.numeric(names(tt.aligned.means)) 

line.col <- c(rgb(216/256, 179/256, 101/256), rgb(1/256, 133/256, 113/256, alpha = 0.5))

# plot
#quartz(width = 9, height = 5)
plot.new()
plot.window(ylim = c(0, 100), xlim = c(-2, 2))
abline(v = c(0, 1), lwd = 1, col = "black")
text(x = 0.5, y = 100, labels = "gene body")
arrows(x0 = 0, x1 = 1, y0 = 95)
points(cc.pos, cc.aligned.means, type = "l", col = line.col[1])
points(tt.pos, tt.aligned.means, type = "l", col = line.col[2])
mtext("Relative Position", side = 1, line = 2)
mtext("Percent Methylation", side = 2, line = 2)
legend("topright", col = line.col, lty = 1, 
legend = c("CC", "TT"), lwd = 3)
axis(1);axis(2)


plot(cc.aligned.means, tt.aligned.means)
segments(0,0,100,100, col="red")

### BINNED
# aligning the centered normallized methylation values for CC averaged across genes
cc.aligned.means <- plot.centered.vals(val.list = gene.methyl.norm[["CC"]], 
plot.label = "CC Genes", ylim = c(0, 500), min.representation = 1, seq.by = 0.05, merge.by = 0,
plot.individual = FALSE, return.means=TRUE, min.upstream = -2000, max.downstream = 2000)

# centered normallized methylation values for TT averaged across genes
tt.aligned.means <- plot.centered.vals(val.list = gene.methyl.norm[["TT"]], 
plot.label = "TT Genes", ylim = c(0, 500), min.representation =1, seq.by = 0.05, merge.by = 0,
plot.individual = FALSE, return.means=TRUE, min.upstream = -2000, max.downstream = 2000)

# position values
cc.pos <- as.numeric(names(cc.aligned.means))
tt.pos <- as.numeric(names(tt.aligned.means)) 

line.col <- c(rgb(216/256, 179/256, 101/256), rgb(1/256, 133/256, 113/256, alpha = 0.5))

# plot
#quartz(width = 9, height = 5)
plot.new()
plot.window(ylim = c(0, 100), xlim = c(-2, 2))
abline(v = c(0, 1), lwd = 1, col = "black")
text(x = 0.5, y = 100, labels = "gene body")
arrows(x0 = 0, x1 = 1, y0 = 95)
points(cc.pos, cc.aligned.means, type = "l", col = line.col[1])
points(tt.pos, tt.aligned.means, type = "l", col = line.col[2])
mtext("Relative Position", side = 1, line = 2)
mtext("Percent Methylation", side = 2, line = 2)
legend("topright", col = line.col, lty = 1, 
legend = c("CC", "TT"), lwd = 3)
axis(1);axis(2)
```

