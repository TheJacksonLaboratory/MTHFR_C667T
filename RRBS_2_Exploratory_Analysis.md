---
title: "MTHFR RRBS Exploratory Analysis"
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
site_meth <- readRDS(file = "output/RRBS/1_Processing_QC/persite_meth.RDS")
# meta
meta <- readRDS(file = "output/RRBS/sample_meta.RDS")

hmc <- read_tsv("data/homocysteine.txt", col_types = cols("c","n"))

meta <- left_join(meta, hmc, by="tag_id")
```

## Pattern of Methylation {.tabset .tabset-fade .tabset-pills}

### With Outliers
#### PCA

```r
pc <- prcomp(site_meth, center = TRUE)
plot(pc)
```

![](RRBS_2_Exploratory_Analysis_files/figure-html/pattern_withOutliers-1.png)<!-- -->

```r
df <- left_join(dplyr::select(meta, age, sex, genotype, sample_id),
                as_tibble(pc$x[,1:2], rownames = "sample_id"),
                by="sample_id")

outliers <- df %>% filter(abs(PC1) > 2000 | abs(PC2) > 2000)
  
ggplot(df, aes(x=PC1,y=PC2)) + 
  geom_point(aes(shape=sex, color=genotype, size=age)) +
  scale_size_discrete(range=c(3,1.5)) +
  geom_text(data=outliers, aes(label=sample_id), vjust = "inward",  hjust = "inward")
```

![](RRBS_2_Exploratory_Analysis_files/figure-html/pattern_withOutliers-2.png)<!-- -->

```r
rm(df)
```



### Without Outliers

```r
site_meth.wou <- site_meth[-which(rownames(site_meth) %in% outliers$sample_id),]
meta.wou <- filter(meta, !sample_id %in% outliers$sample_id)
```

#### PCA

```r
pc <- prcomp(site_meth.wou, center = TRUE)
plot(pc)
```

![](RRBS_2_Exploratory_Analysis_files/figure-html/pattern__noOutliers-1.png)<!-- -->

```r
left_join(select(meta.wou, age, sex, genotype, sample_id),
                as_tibble(pc$x[,1:2], rownames = "sample_id"),
                by="sample_id") %>%
  ggplot(aes(x=PC1,y=PC2)) + 
  geom_point(aes(shape=sex, color=genotype, size=age)) +
  scale_size_discrete(range=c(3,1.5)) 
```

![](RRBS_2_Exploratory_Analysis_files/figure-html/pattern__noOutliers-2.png)<!-- -->

```r
#+
#  geom_text(data=outliers, aes(label=sample_id), vjust = -0.5,  hjust = -0.5)
```
  
 


## Global methylation {.tabset .tabset-fade .tabset-pills}
### Counts
Calculating global methylation as the fraction of sites determined as methylated by bismark for each sample.

```r
global_meth_counts <- read_tsv("output/RRBS/1_Processing_QC/nf_methylseq_output/bismark_summary/bismark_summary_report.txt") %>%
  tibble(.name_repair = "universal") %>%
  mutate(sample_id = substr(File, 1, 5),
         global_meth_counts = Methylated.CpGs / (Methylated.CpGs + Unmethylated.CpGs) ) %>%
  select(global_meth_counts, sample_id)

meta.gm <-  left_join(meta, global_meth_counts, by="sample_id")
```

Using Beta regression to test differences in global methylation as a function of sex, age, and genotype. None of the factors and interactions is significant.

```r
br.df <- meta.gm %>%
  mutate(genotype = recode(genotype, "CC"=0, "CT"=1, "TT"=2),
         sex = recode(sex, "F"=0, "M"=1),
         age = recode(age, "6"=0, "18"=1))
   
model = betareg(global_meth_counts ~ genotype + age + sex, data = br.df)
emmeans::joint_tests(model)
```

```
##  model term df1 df2 F.ratio p.value
##  genotype     1 Inf   0.235 0.6278 
##  age          1 Inf   0.550 0.4582 
##  sex          1 Inf   0.315 0.5744
```

```r
ggplot(meta.gm, aes(x=genotype, y=global_meth_counts, color=age)) +
  geom_point() +
  facet_wrap(~ sex)
```

![](RRBS_2_Exploratory_Analysis_files/figure-html/global_meth_counts_betareg-1.png)<!-- -->

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
         sex = recode(sex, "F"=0, "M"=1),
         age = recode(age, "6"=0, "18"=1))
   
model = betareg(global_meth_perc ~ genotype + age + sex, data = br.df)
emmeans::joint_tests(model)
```

```
##  model term df1 df2 F.ratio p.value
##  genotype     1 Inf   0.221 0.6379 
##  age          1 Inf   1.159 0.2818 
##  sex          1 Inf   0.027 0.8694
```

```r
ggplot(meta.gm, aes(x=genotype, y=global_meth_perc, color=age)) +
  geom_point() +
  facet_wrap(~ sex)
```

![](RRBS_2_Exploratory_Analysis_files/figure-html/global_meth_percent_betareg-1.png)<!-- -->

### Correlation with homocysteine level


```r
ggplot(meta.gm, aes(x = hmc, y = global_meth_counts)) + 
  geom_point(aes(shape=sex, color=genotype, size=age)) +
  scale_size_discrete(range=c(3,1.5)) +
  ggtitle("Counts")
```

<img src="RRBS_2_Exploratory_Analysis_files/figure-html/hmc-1.png" width="50%" />

```r
ggplot(meta.gm, aes(x = hmc, y = global_meth_perc)) + 
  geom_point(aes(shape=sex, color=genotype, size=age)) +
  scale_size_discrete(range=c(3,1.5)) +
  ggtitle("Percent")
```

<img src="RRBS_2_Exploratory_Analysis_files/figure-html/hmc-2.png" width="50%" />


## Output
Uploading sample metadata with their global methylation values to [Synapse](https://www.synapse.org/#!Synapse:syn23573590/wiki/607402)

```r
write_csv(meta.gm, file="output/RRBS/2_Exploratory_Analysis/RRBS_meta_globalmeth.csv")

# upload meta_global_meth to synapse
synLogin() 

file <- File("output/RRBS/2_Exploratory_Analysis/RRBS_meta_globalmeth.csv",
             contentType = "text/plain",
             description = "RRBS sample metadata with their global methylation values",
             parent = "syn25174486")

provenance <- Activity(used = c("syn25174542", "syn25174517"),
                       executed = "https://github.com/TheJacksonLaboratory/MTHFR_C667T/blob/master/RRBS_2_Exploratory_Analysis.Rmd")

file <- synStore(file, activity=provenance)
```

