## code to prepare `michaelis_menten` dataset goes here

library(tidyverse)
library(magrittr)
library(hitac2)

michaelis_menten <- read.csv(file = "~/Seafile/ToshibaWorkCrypt/Rfolder/NewPackages/Examples/michaelismenten_example.csv", sep = ",") %>%
  as_tibble() %>%
  transmute(ExpID = X, Concentration = Concentration_RuBP, Time = Time, Activity = Activity_nmol)

write_csv(michaelis_menten, "data-raw/michaelis_menten.csv")
usethis::use_data(michaelis_menten, overwrite = TRUE)
