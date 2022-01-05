# Microbial genomics extensions / extra functions in R


### Dependencies

* dplyr
* ggplot2
* stringr
* gplots
* RColorBrewer  
* ggtree
* ape
* castor
* purrr
* glue
* tibble
* seqinr
* genoplotR
* vegan
* UpSetR
* cowplot


# Install dependencies

```
R
install.packages(c("dplyr", "ggplot2", "stringr", "gplots", "RColorBrewer", "ape", "castor", "purrr", "glue", "tibble", "seqinr","genoplotR","vegan", "UpSetR", "cowplot"))

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(c("ggtree"))
```


# Install micro.gen.extra

```R
install.packages("devtools")
library(devtools)
install_github("djw533/micro.gen.extra")
```
