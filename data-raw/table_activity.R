## code to prepare `table_activity` dataset goes here

library(tidyverse)
library(magrittr)
library(hitac2)

# read data from 20200207
design <- read_csv("~/Seafile/ToshibaWorkCrypt/Rfolder/Project8/20200207/20200207_design_001.csv")
write_csv(design, "inst/extdata/design.csv")
treatment <- read_csv("~/Seafile/ToshibaWorkCrypt/Rfolder/Project8/20200207/20200207_treatment.csv")
write_csv(treatment, "inst/extdata/treatment.csv")

ha_data <- left_join(ha_data, design, by = "Hitachi_ID")
table_activity <- hitachi_tibble(data = ha_data, treatment = treatment, V_RuBP = 4e-06, C_RuBP = 6e-04, convert = TRUE)

usethis::use_data(table_activity, overwrite = TRUE)
