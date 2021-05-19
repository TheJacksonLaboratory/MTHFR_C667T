---
title: "MTHFR RRBS_Pilot Exploratory Analysis"
author: "Annat Haber"
date: '2021-05-19'
output: 
  html_document:
    toc: true
    toc_depth: 3
    code_folding: hide
---




```r
rm(list=ls())

suppressPackageStartupMessages({
  library(synapser)
  #library(methylKit)
  library(tidyverse)
  library(betareg)
  library(emmeans)
  library(ComplexHeatmap)
})
```

## Input
Input is percent methylation in each site for each sample (samples x sites matrix), including only sites within CpG context that are covered in all samples with at least 10x coverage. 

```r
site_meth <- readRDS(file = "output/RRBS_Pilot/1_Processing_QC/persite_meth.RDS")
# meta
meta <- readRDS(file = "output/RRBS_Pilot/sample_meta.RDS")
```

## Pattern of Methylation
### PCA

```r
pc <- prcomp(site_meth, center = TRUE)
plot(pc)
```

![](RRBS_Pilot_2_Exploratory_Analysis_files/figure-html/pattern_withOutliers-1.png)<!-- -->

```r
df <- left_join(dplyr::select(meta, age, genotype, sample_id),
                as_tibble(pc$x[,1:2], rownames = "sample_id"),
                by="sample_id")

#outliers <- df %>% filter(abs(PC1) > 2000 | abs(PC2) > 2000)
  
ggplot(df, aes(x=PC1,y=PC2)) + 
  geom_point(aes(color=genotype, size=age)) +
  scale_size_discrete(range=c(3,1.5)) 
```

![](RRBS_Pilot_2_Exploratory_Analysis_files/figure-html/pattern_withOutliers-2.png)<!-- -->

```r
#+
#  geom_text(data=outliers, aes(label=sample_id), vjust = "inward",  hjust = -0.5)
```
  
### Methylation by site  

```r
df = t(site_meth) #- colMeans(site_meth) #relative to the mean across samples

cols = list(age = c("4" = "blue", "13" = "darkblue"),
           genotype = c("CC" = "yellow", "CT" = "orange", "TT" = "red"))
# Create the heatmap annotation
ha <- HeatmapAnnotation(age = meta$age, genotype = meta$genotype, col=cols)

Heatmap(df,
        show_row_dend = FALSE,
        show_row_names = FALSE,
        heatmap_legend_param = list(title = "Methlation %"),
        #cluster_rows = FALSE, 
        #clustering_method_columns="ward.D2", 
        col = hcl.colors(12, "YlOrRd", rev = TRUE),
        top_annotation = ha
        )
```

![](RRBS_Pilot_2_Exploratory_Analysis_files/figure-html/bySite_withOutliers-1.png)<!-- -->



## Global methylation {.tabset .tabset-fade .tabset-pills}
### Counts
Calculating global methylation as the fraction of sites determined as methylated by bismark for each sample.

```r
global_meth_counts <- read_tsv("output/RRBS_Pilot/1_Processing_QC/nf_methylseq_output/bismark_summary/bismark_summary_report.txt") %>%
  tibble(.name_repair = "universal") %>%
  mutate(sample_id = substr(File, 1, 5),
         global_meth_counts = Methylated.CpGs / (Methylated.CpGs + Unmethylated.CpGs) ) %>%
  select(global_meth_counts, sample_id)
```

Using Beta regression to test differences in global methylation as a function of sex, age, and genotype. None of the factors and interactions is significant.

```r
meta.gm <-  left_join(meta, global_meth_counts, by="sample_id")

br.df <- meta.gm %>%
  mutate(genotype = recode(genotype, "CC"=0, "CT"=1, "TT"=2),
         age = recode(age, "4"=0, "13"=1))
   
model = betareg(global_meth_counts ~ genotype + age, data = br.df)
emmeans::joint_tests(model)
```

```
##  model term df1 df2 F.ratio p.value
##  genotype     1 Inf  18.793 <.0001 
##  age          1 Inf   0.289 0.5907
```

```r
ggplot(meta.gm, aes(x=genotype, y=global_meth_counts, color=age)) +
  geom_point() +
  facet_wrap(~ age)
```

![](RRBS_Pilot_2_Exploratory_Analysis_files/figure-html/global_meth_counts_betareg-1.png)<!-- -->

### Mean Percent
Calculating global methylation as the mean percent methylation for each sample, averaged across sites.

```r
global_meth_perc <- data.frame(global_meth_perc=rowMeans(site_meth)/100, sample_id=rownames(site_meth))
```

Using Beta regression to test differences in global methylation as a function of sex, age, and genotype. None of the factors and interactions is significant.

```r
meta.gm <-  left_join(meta.gm, global_meth_perc, by="sample_id")

br.df <- meta.gm %>%
  mutate(genotype = recode(genotype, "CC"=0, "CT"=1, "TT"=2),
         age = recode(age, "4"=0, "13"=1))
   
model = betareg(global_meth_perc ~ genotype + age, data = br.df)
emmeans::joint_tests(model)
```

```
##  model term df1 df2 F.ratio p.value
##  genotype     1 Inf  29.405 <.0001 
##  age          1 Inf   0.068 0.7946
```

```r
ggplot(meta.gm, aes(x=genotype, y=global_meth_perc, color=age)) +
  geom_point() +
  facet_wrap(~ age)
```

![](RRBS_Pilot_2_Exploratory_Analysis_files/figure-html/global_meth_perc_betareg-1.png)<!-- -->

## Output
Uploading sample metadata with their global methylation values to [Synapse](https://www.synapse.org/#!Synapse:syn23573590/wiki/607402)

```r
write_csv(meta.gm, file="output/RRBS_Pilot/2_Exploratory_Analysis/RRBS_Pilot_meta_globalmeth.csv")

# upload meta_global_meth to synapse
synLogin() 

file <- File("output/RRBS_Pilot/2_Exploratory_Analysis/RRBS_Pilot_meta_globalmeth.csv",
             contentType = "text/plain",
             description = "RRBS Pilot sample metadata with their global methylation values",
             parent = "syn25174488")

provenance <- Activity(used = c("syn25175205", "syn25175218"),
                       executed = "https://github.com/TheJacksonLaboratory/MTHFR_C667T/blob/master/RRBS_Pilot_2_Exploratory_Analysis.Rmd")

file <- synStore(file, activity=provenance)
```
