---
title: "RRBS Pilot Sample Metadata"
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

Sample size included in the RRBS Pilot analysis by age and genotype.  
All samples are female.

```r
meta <- read_csv("data/RRBS_Pilot/RRBS_Pilot_sample_metadata.csv", 
                 col_types = cols(.default = "c")) %>%
  mutate(group = paste(genotype, age, sep="_")) %>%
  arrange(sample_id)

# add treatment vector for methylKit operations
meta <- left_join(meta,
              data.frame(treatment=0:5, group=unique(meta$group)),
              by="group")

meta %>% 
  group_by(age, genotype) %>%
  summarise(n=n())
```

<div class="kable-table">

|age |genotype |  n|
|:---|:--------|--:|
|13  |CC       |  3|
|13  |CT       |  4|
|13  |TT       |  4|
|4   |CC       |  4|
|4   |CT       |  4|
|4   |TT       |  4|

</div>


```r
saveRDS(meta, file="output/RRBS_Pilot/sample_meta.RDS")
write_csv(meta, file="output/RRBS_Pilot/RRBS_Pilot_sample_meta.csv")

synLogin() 

file <- File("output/RRBS_Pilot/RRBS_Pilot_sample_meta.csv",
             contentType = "text/plain",
             description = "RRBS_Pilot sample metadata",
             parent = "syn25174488")
provenance <- Activity(executed = "https://github.com/TheJacksonLaboratory/MTHFR_C667T")
file <- synStore(file, activity=provenance)
```
```
