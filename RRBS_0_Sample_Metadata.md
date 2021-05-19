---
title: "RRBS Sample Metadata"
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
  library(synapser)
})
```

Sample size included in the RRBS analysis by sex, age, and genotype.  

```r
meta <- read_csv("data/RRBS/RRBS_sample_metadata.csv", col_types = cols(.default = "c")) %>%
  mutate(genotype = recode(genotype, "-/-"="TT", "-/+"="CT", "+/+"="CC"),
         group = paste(genotype, age, sex, sep="_"),
         sample_id = str_pad(tag_id, 5, pad = "0")) %>%
  arrange(sample_id)

# add treatment vector for methylKit operations
meta <- left_join(meta,
              data.frame(treatment=0:11, group=unique(meta$group)),
              by="group")

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
saveRDS(meta, file="output/RRBS/sample_meta.RDS")
write_csv(meta, file="output/RRBS/RRBS_sample_meta.csv")

synLogin() 

file <- File("output/RRBS/RRBS_sample_meta.csv",
             contentType = "text/plain",
             description = "RRBS sample metadata",
             parent = "syn25174486")
provenance <- Activity(executed = "https://github.com/TheJacksonLaboratory/MTHFR_C667T")
file <- synStore(file, activity=provenance)
```
