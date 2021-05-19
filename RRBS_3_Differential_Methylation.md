---
title: "RRBS Differential Methylation Annotation"
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
knitr::opts_chunk$set(cache = TRUE, 
                      #cache.path = "output/RRBS/3_Differential_Methylation",
                      message=FALSE, warning=FALSE)
```


```r
rm(list=ls())

suppressPackageStartupMessages({
  library(methylKit)
  library(tidyverse)

  library(GenomicRanges)
  library(genomation)
  library(ChIPpeakAnno)
  library(clusterProfiler)
  library(GOsummaries)
  library(org.Mm.eg.db)
  library(biomaRt)
  library(vcd)
})
```

## Input

```r
united_meth <- readRDS(file="output/RRBS/1_Processing_QC/united_meth.RDS")
meta <- read_csv(file="output/RRBS/2_Exploratory_Analysis/RRBS_meta_globalmeth.csv",
                 col_types = cols(.default = "c"))
```

## differentially Methylated CpG sites  {.tabset .tabset-fade .tabset-pills}
Finding differentially methylated CpG sites using age as covariate.  
Mapping to known genes by finding overlap between the snp's and mm10 reference.  
The reference annotation was downloaded directly from http://genome.ucsc.edu/cgi-bin/hgTables 


```r
# read reference and convert to GRanges
refseq <- read_tsv("data/mm10.refseq.genes.txt") 

ref.gr <- GRanges(
        seqnames = substr(refseq$chrom, 4,5),
        ranges=IRanges(start = refseq$txStart, end = refseq$txEnd),
        mgi_name=refseq$name2
        )
names(ref.gr) <- refseq$name
rm(refseq)
```

### TT vs CC 

```r
mm <- united_meth
cov <- meta %>%
  filter(genotype!="CT") %>%
  mutate(genotype = ifelse(genotype=="TT", 1, 0))
  
diffM_TT <- reorganize(mm, 
                       sample.ids = cov$sample_id, 
                       treatment= cov$genotype) %>%
            calculateDiffMeth(covariates = as.data.frame(cov$age))

rm(mm, cov)
```



```r
diffM <- diffM_TT
# extract significant sites and convert to GRanges
snp.gr <- getMethylDiff(diffM, difference = 10, qvalue=0.01, type="all") %>%
  as("GRanges")
names(snp.gr) <- paste(snp.gr@seqnames@values, snp.gr@ranges@start, sep=":")

# Find overlap with reference
overlap <- findOverlaps(snp.gr, ref.gr)
hits <- ref.gr$mgi_name[slot(overlap, "to")]
names(hits) <- names(snp.gr)[slot(overlap, "from")]
hits <- hits[!duplicated(hits)]
snpDMgenes <- snp.gr[names(hits),]
mcols(snpDMgenes)$mgi_name <- hits

snpDMgenes <- data.frame(snpDMgenes) %>%
  dplyr::select(chromosome=seqnames, everything()) %>%
  arrange(chromosome, meth.diff)

# separate hypo and hyper methylated gene sets
snpDMgenes_TT <- list("hypo"= snpDMgenes %>%
                        filter(qvalue < 0.01 & meth.diff < 0) %>%
                        arrange(meth.diff),
                      "hyper"= snpDMgenes %>%
                        filter(qvalue < 0.01 & meth.diff > 0) %>%
                        arrange(-meth.diff)
                      )

rm(diffM, snp.gr, snpDMgenes, hits, overlap)
```


```r
diffMethPerChr(diffM_TT, plot=TRUE, qvalue.cutoff=0.01, meth.cutoff=10)
```

![](RRBS_3_Differential_Methylation_files/figure-html/diffM_TT_plot-1.png)<!-- -->

```r
diffMethPerChr(diffM_TT, plot=TRUE, qvalue.cutoff=0.01, meth.cutoff=0)
```

![](RRBS_3_Differential_Methylation_files/figure-html/diffM_TT_plot-2.png)<!-- -->

### CT vs CC 

```r
mm <- united_meth
cov <- meta %>%
  filter(genotype!="TT") %>%
  mutate(genotype = ifelse(genotype=="CT", 1, 0))
  
diffM_CT <- reorganize(mm, 
                       sample.ids = cov$sample_id, 
                       treatment= cov$genotype) %>%
            calculateDiffMeth(covariates = as.data.frame(cov$age))

rm(mm, cov)
```


```r
diffM <- diffM_CT
# extract significant sites and convert to GRanges
snp.gr <- getMethylDiff(diffM, difference = 0, qvalue=0.01, type="all") %>%
  as("GRanges")
names(snp.gr) <- paste(snp.gr@seqnames@values, snp.gr@ranges@start, sep=":")

# Find overlap with reference
overlap <- findOverlaps(snp.gr, ref.gr)
hits <- ref.gr$mgi_name[slot(overlap, "to")]
names(hits) <- names(snp.gr)[slot(overlap, "from")]
hits <- hits[!duplicated(hits)]
snpDMgenes <- snp.gr[names(hits),]
mcols(snpDMgenes)$mgi_name <- hits

snpDMgenes <- data.frame(snpDMgenes) %>%
  dplyr::select(chromosome=seqnames, everything()) %>%
  arrange(chromosome, meth.diff)

# separate hypo and hyper methylated gene sets
snpDMgenes_CT <- list("hypo"= snpDMgenes %>%
                        filter(qvalue < 0.01 & meth.diff < 0) %>%
                        arrange(meth.diff),
                      "hyper"= snpDMgenes %>%
                        filter(qvalue < 0.01 & meth.diff > 0) %>%
                        arrange(-meth.diff)
                      )

rm(diffM, snp.gr, snpDMgenes, hits, overlap)
```


```r
diffMethPerChr(diffM_CT, plot=TRUE, qvalue.cutoff=0.01, meth.cutoff=10)
```

![](RRBS_3_Differential_Methylation_files/figure-html/diffM_CT_plot-1.png)<!-- -->

```r
diffMethPerChr(diffM_CT, plot=TRUE, qvalue.cutoff=0.01, meth.cutoff=0)
```

![](RRBS_3_Differential_Methylation_files/figure-html/diffM_CT_plot-2.png)<!-- -->

## Gene parts annotation {.tabset .tabset-fade .tabset-pills}
Gene parts distribution for site that are differentially-methylated 

### TT vs CC

```r
diffM <- diffM_TT

# Get hypermethylated sites with difference > 10% and convert to GRanges
diffM10.hyper <- getMethylDiff(diffM, difference=10, qvalue=0.01, type="hyper") %>%
  as("GRanges")

# Get hypomethylated sites with difference > 10% and convert to GRanges
diffM10.hypo <- getMethylDiff(diffM, difference=10, qvalue=0.01, type="hypo") %>%
  as("GRanges")


# Get gene parts annotation
gene.obj=readTranscriptFeatures("data/mm10.refseq.genes.bed")

hypo10.ann <- annotateWithGeneParts(diffM10.hypo,  gene.obj)

hyper10.ann <- annotateWithGeneParts(diffM10.hyper, gene.obj)

# Get flank annotation
cpg.obj=readFeatureFlank("data/mm10.refseq.genes.bed",
                           feature.flank.name=c("CpGi","shores"))

diffCpGann.hypo10 <- annotateWithFeatureFlank(diffM10.hypo,
                                    cpg.obj$CpGi,cpg.obj$shores,
                         feature.name="CpGi",flank.name="shores")

diffCpGann.hyper10 <- annotateWithFeatureFlank(diffM10.hyper,
                                    cpg.obj$CpGi,cpg.obj$shores,
                         feature.name="CpGi",flank.name="shores")
rm(diffM)
```


```r
plotTargetAnnotation(hypo10.ann, precedence=TRUE,
    main="Hypo-methylated gene parts; > 10% diff.")

plotTargetAnnotation(hyper10.ann, precedence=TRUE,
    main="Hyper-methylated gene parts; > 10% diff.")


plotTargetAnnotation(diffCpGann.hypo10,col=c("green","gray","white"),
       main="Hypo-methylated flank annotation; > 10% diff.")

plotTargetAnnotation(diffCpGann.hyper10,col=c("green","gray","white"),
       main="Hyper-methylated flank annotation; > 10% diff.")
```

<img src="RRBS_3_Differential_Methylation_files/figure-html/gene_parts_plots_TT-1.png" width="50%" /><img src="RRBS_3_Differential_Methylation_files/figure-html/gene_parts_plots_TT-2.png" width="50%" /><img src="RRBS_3_Differential_Methylation_files/figure-html/gene_parts_plots_TT-3.png" width="50%" /><img src="RRBS_3_Differential_Methylation_files/figure-html/gene_parts_plots_TT-4.png" width="50%" />

### CT vs CC

```r
diffM <- diffM_CT

# Get hypermethylated sites with difference > 10% and convert to GRanges
diffM10.hyper <- getMethylDiff(diffM, difference=10, qvalue=0.01, type="hyper") %>%
  as("GRanges")

# Get hypomethylated sites with difference > 10% and convert to GRanges
diffM10.hypo <- getMethylDiff(diffM, difference=10, qvalue=0.01, type="hypo") %>%
  as("GRanges")


# Get gene parts annotation
gene.obj=readTranscriptFeatures("data/mm10.refseq.genes.bed")

hypo10.ann <- annotateWithGeneParts(diffM10.hypo,  gene.obj)

hyper10.ann <- annotateWithGeneParts(diffM10.hyper, gene.obj)

# Get flank annotation
cpg.obj=readFeatureFlank("data/mm10.refseq.genes.bed",
                           feature.flank.name=c("CpGi","shores"))

diffCpGann.hypo10 <- annotateWithFeatureFlank(diffM10.hypo,
                                    cpg.obj$CpGi,cpg.obj$shores,
                         feature.name="CpGi",flank.name="shores")

diffCpGann.hyper10 <- annotateWithFeatureFlank(diffM10.hyper,
                                    cpg.obj$CpGi,cpg.obj$shores,
                         feature.name="CpGi",flank.name="shores")
rm(diffM)
```


```r
plotTargetAnnotation(hypo10.ann, precedence=TRUE,
    main="Hypo-methylated gene parts")

plotTargetAnnotation(hyper10.ann, precedence=TRUE,
    main="Hyper-methylated gene parts")


plotTargetAnnotation(diffCpGann.hypo10,col=c("green","gray","white"),
       main="Hypo-methylated flank annotation")

plotTargetAnnotation(diffCpGann.hyper10,col=c("green","gray","white"),
       main="Hyper-methylated flank annotation")
```

<img src="RRBS_3_Differential_Methylation_files/figure-html/gene_parts_plots_CT-1.png" width="50%" /><img src="RRBS_3_Differential_Methylation_files/figure-html/gene_parts_plots_CT-2.png" width="50%" /><img src="RRBS_3_Differential_Methylation_files/figure-html/gene_parts_plots_CT-3.png" width="50%" /><img src="RRBS_3_Differential_Methylation_files/figure-html/gene_parts_plots_CT-4.png" width="50%" />

## KEGG Annotation {.tabset .tabset-fade .tabset-pills}

### TT vs CC


```r
sig.gs <- snpDMgenes_TT

# get entrez id
sig.gs.ez <- list()
for (gs in names(sig.gs)) {
  sig.gs.ez[[gs]] <- bitr(sig.gs[[gs]]$mgi_name, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Mm.eg.db")$ENTREZID
}

kegg_anno_TT <- list()
for (gs in names(sig.gs.ez)) {
  cat("\n Annotating:", gs,"-methylated genes\n")
  
  anno <- enrichKEGG(gene = sig.gs.ez[[gs]],
                organism = "mmu",
                keyType = "ncbi-geneid",
                pvalueCutoff = 0.5,
                pAdjustMethod = "BH",
                minGSSize = 10,
                maxGSSize = 500,
                qvalueCutoff = 0.5,
                use_internal_data = FALSE)
  
  print(barplot(anno, showCategory=20, title=paste0("KEGG pathways for ", gs, "-methylated sites")))
  anno <- setReadable(anno, OrgDb = org.Mm.eg.db, keyType="ENTREZID")
  anno <- as.data.frame(anno@result)
  write_csv(anno, file=paste0("output/RRBS/3_Differential_Methylation/kegg_anno_TT_", gs, ".csv"))

  kegg_anno_TT[[gs]] <- anno
  rm(anno)
  }
```

![](RRBS_3_Differential_Methylation_files/figure-html/kegg_annotation_TT-1.png)<!-- -->![](RRBS_3_Differential_Methylation_files/figure-html/kegg_annotation_TT-2.png)<!-- -->

```r
rm(sig.gs, sig.gs.ez)
```

### CT vs CC

```r
sig.gs <- snpDMgenes_CT

# get entrez id
sig.gs.ez <- list()
for (gs in names(sig.gs)) {
  sig.gs.ez[[gs]] <- bitr(sig.gs[[gs]]$mgi_name, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Mm.eg.db")$ENTREZID
}
kegg_anno_CT <- list()
for (gs in names(sig.gs.ez)) {
  cat("\n Annotating:", gs,"-methylated genes\n")
  
  anno <- enrichKEGG(gene = sig.gs.ez[[gs]],
                organism = "mmu",
                keyType = "ncbi-geneid",
                pvalueCutoff = 0.5,
                pAdjustMethod = "BH",
                minGSSize = 10,
                maxGSSize = 500,
                qvalueCutoff = 0.5,
                use_internal_data = FALSE)
  
  print(barplot(anno, showCategory=20, title=paste0("KEGG pathways for ", gs, "-methylated sites")))
  anno <- setReadable(anno, OrgDb = org.Mm.eg.db, keyType="ENTREZID")
  anno <- as.data.frame(anno@result)
  write_csv(anno, file=paste0("output/RRBS/3_Differential_Methylation/kegg_anno_CT_", gs, ".csv"))

  kegg_anno_CT[[gs]] <- anno
  rm(anno)
  }
```

![](RRBS_3_Differential_Methylation_files/figure-html/kegg_annotation_CT-1.png)<!-- -->![](RRBS_3_Differential_Methylation_files/figure-html/kegg_annotation_CT-2.png)<!-- -->

```r
rm(sig.gs, sig.gs.ez)
```

## Cell-specificity {.tabset .tabset-fade .tabset-pills}
Cell type markers in the brain are taken from McKenzie et al. (Supplemental File 1, DOI:10.1038/s41598-018-27293-5).
Cell specificity for a given gene is the minimum log fold change when pairwise comparing the expression of the gene in the given cell type with the expression of the gene in all other cell types.
Here I'm using cell specificity score to annotate the differentially-methylated genes, with chi-square to test for biased representations of any cell type genes in either the hypo- or hyper-methylated genes, and a plot of deviations from expectation for visualization.  


```r
cell.types <- c("Astroctyes", "Endothelial", "Microglia", "Neurons", "Oligodendrocytes")
sets <- c("hypo","hyper")

# read cell-specificity reference
hum.specific <- read.csv("data/McKenzie_2017_Brain_Cell_Specific_Markers_Human_Specificity.csv", header=T, stringsAsFactors=F)
hum.specific <- split(hum.specific, hum.specific$Celltype, drop=TRUE) 
names(hum.specific) <- cell.types

# convert HGNC.ID in the specificity reference to MGI names
human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")

hum.specific.mg <- list()
for (ct in names(hum.specific)) {
  hum.specific.mg[[ct]] = getLDS(attributes = c("hgnc_symbol"), filters = "hgnc_symbol", 
                              values = hum.specific[[ct]]$gene, mart = human,
                              attributesL = c("mgi_symbol"), martL = mouse, uniqueRows=T)$MGI.symbol
  }
```

### TT vs CT
There is in general little overlap between the cell markers and our differentially-methylated genes, and no significant association between them. Hypomethylated genes are somewhat enriched for astrocytes and endothelials (BBB?) marker, and possibly neurons. The hypermethylated genes are more microglia and oligodendrocytes, but none of it is significant.
Note, the chi-square test compares hypomethylated sites relative to hypermethylated sites, not within each set. 

```r
sig.gs <- snpDMgenes_TT

# get gene sets
gene.sets <- lapply(sig.gs, function(l){dplyr::select(l, mgi_name)})

gene.sets.hg <- list()
for (s in names(sig.gs)) {
    gene.sets.hg[[s]] = getLDS(attributes = c("mgi_symbol"), filters = "mgi_symbol", 
                              values = sig.gs[[s]]$mgi_name, mart = mouse,
                              attributesL = c("hgnc_symbol"), martL = human, uniqueRows=T)$HGNC.symbol

}

# Number of marker genes for each cell type in each gene set
specific.set <- sapply(sets, function(s) {
  sapply(cell.types, function(ct) length(intersect(hum.specific[[ct]]$gene, gene.sets.hg[[s]])), USE.NAMES=T)
})

print(specific.set)
```

```
##                  hypo hyper
## Astroctyes          5    13
## Endothelial         6     2
## Microglia           4     2
## Neurons            12     9
## Oligodendrocytes    4     3
```

```r
chisq.test(specific.set)
```

```
## 
## 	Pearson's Chi-squared test
## 
## data:  specific.set
## X-squared = 6.7345, df = 4, p-value = 0.1506
```

```r
assoc(t(specific.set), shade=TRUE,labeling_args = list(abbreviate_labs = c(6, 9)))
```

![](RRBS_3_Differential_Methylation_files/figure-html/cell_spec_TT-1.png)<!-- -->

```r
rm(sig.gs, gene.sets, gene.sets.hg, specific.set)
```


### CT vs CC

```r
sig.gs <- snpDMgenes_CT

# get gene sets
gene.sets <- lapply(sig.gs, function(l){dplyr::select(l, mgi_name)})

gene.sets.hg <- list()
for (s in names(sig.gs)) {
    gene.sets.hg[[s]] = getLDS(attributes = c("mgi_symbol"), filters = "mgi_symbol", 
                              values = sig.gs[[s]]$mgi_name, mart = mouse,
                              attributesL = c("hgnc_symbol"), martL = human, uniqueRows=T)$HGNC.symbol

}

# Number of marker genes for each cell type in each gene set
specific.set <- sapply(sets, function(s) {
  sapply(cell.types, function(ct) length(intersect(hum.specific[[ct]]$gene, gene.sets.hg[[s]])), USE.NAMES=T)
})

print(specific.set)
```

```
##                  hypo hyper
## Astroctyes        106   113
## Endothelial        97    87
## Microglia          62    86
## Neurons           119   133
## Oligodendrocytes   77   106
```

```r
chisq.test(specific.set)
```

```
## 
## 	Pearson's Chi-squared test
## 
## data:  specific.set
## X-squared = 5.9032, df = 4, p-value = 0.2065
```

```r
assoc(t(specific.set), shade=TRUE,labeling_args = list(abbreviate_labs = c(6, 9)))
```

![](RRBS_3_Differential_Methylation_files/figure-html/cell_spec_CT-1.png)<!-- -->

```r
rm(sig.gs, gene.sets, gene.sets.hg, specific.set)
```

## Output
Saving differential methylation objects and kegg annotations as RDS for subsequent analyses.   Saving lists of significantly methylated sites and their genes as csv.

```r
saveRDS(diffM_TT, file="output/RRBS/3_Differential_Methylation/diffM_TT.RDS")
saveRDS(diffM_CT, file="output/RRBS/3_Differential_Methylation/diffM_CT.RDS")

write_csv(snpDMgenes_TT[[1]], file = "output/RRBS/3_Differential_Methylation/snpDMgenes_TT.csv")
write_csv(snpDMgenes_TT[[2]], file = "output/RRBS/3_Differential_Methylation/snpDMgenes_TT.csv", append=TRUE)
write_csv(snpDMgenes_CT[[1]], file = "output/RRBS/3_Differential_Methylation/snpDMgenes_CT.csv")
write_csv(snpDMgenes_CT[[2]], file = "output/RRBS/3_Differential_Methylation/snpDMgenes_CT.csv", append=TRUE)

saveRDS(kegg_anno_TT, file="output/RRBS/3_Differential_Methylation/kegg_anno_TT.RDS")
saveRDS(kegg_anno_CT, file="output/RRBS/3_Differential_Methylation/kegg_anno_CT.RDS")
```
