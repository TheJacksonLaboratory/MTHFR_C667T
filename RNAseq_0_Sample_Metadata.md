---
title: "RNAseq Sample Metadata"
author: "Annat Haber"
date: '2021-05-19'
output:
  html_document:
    toc: true
    code_folding: hide
    df_print: "kable"
---


```r
knitr::opts_knit$set(root.dir = rprojroot::find_rstudio_root_file())
knitr::opts_chunk$set(message=FALSE, warning=FALSE, cache = FALSE)
```


```r
suppressPackageStartupMessages({
  library(tidyverse)
})
```

Sample size included in the RNAseq analysis by sex, age, and genotype.  

```r
meta <- read_csv("data/RNAseq/RNAseq_sample_metadata.csv", col_types = cols(.default = "c")) %>%
  mutate(genotype = recode(genotype, "-/-"="TT", "-/+"="CT", "+/+"="CC"),
         group = paste(genotype, age, sex, sep="_"),
         sample_id = str_pad(tag_id, 5, pad = "0")) %>%
  arrange(sample_id)

meta %>% 
  group_by(sex, age, genotype) %>%
  summarise(n=n())
```

<div class="kable-table">

|sex |age |genotype |  n|
|:---|:---|:--------|--:|
|F   |18  |CC       |  6|
|F   |18  |CT       |  6|
|F   |18  |TT       |  6|
|F   |6   |CC       |  6|
|F   |6   |CT       |  6|
|F   |6   |TT       |  7|
|M   |18  |CC       |  6|
|M   |18  |CT       |  6|
|M   |18  |TT       |  6|
|M   |6   |CC       |  6|
|M   |6   |CT       |  6|
|M   |6   |TT       |  6|

</div>

```r
write_csv(meta, file="output/RNAseq/RNAseq_sample_meta.csv")
```
