library(tidyverse)

results <-
    read.delim("data/raw/analysis_result.txt",
               header = TRUE) %>%
    janitor::clean_names() %>%
    select(-c(na:na_29))  %>%
    filter(str_detect(well, "([0-9]{1,2})")) %>%
    separate(sample_name, c('conidtion', 'sample'))

melt <- read.delim("data/raw/meltcurve_result.txt",
                   header = TRUE) %>%
    janitor::clean_names() %>%
    select(-c(na:na_44)) %>%
    filter(str_detect(well, "([0-9]{1,2})")) %>%
    separate(sample_name, c('conidtion', 'sample'))
