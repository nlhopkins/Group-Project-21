library(tidyverse)
library(ggpubr)
library(rstatix)

theme <- theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    text = element_text(size = 12),
    axis.text.x = element_text(hjust = 0.5, size = 12),
    axis.title.x = element_text(face = "bold", size = 14),
    axis.text.y = element_text(size = 12),
    axis.title.y = element_text(face = "bold", size = 14),
    axis.ticks = element_blank(),
    strip.placement = "outside",
    strip.background = element_rect(colour = "white", fill = "white"),
    strip.text.x = element_text(size = 12, face = "bold"),
    panel.spacing = unit(0, "lines"),
    plot.background = element_rect(fill = "transparent", color = NA),
    panel.grid.major.x = element_blank()
)



## ----✷----✷---- location ----✷----✷----✷----✷----✷----✷----✷----✷----
location <-
    read.csv("../data/raw/periphery.csv",
             header = TRUE) %>%
    janitor::clean_names() %>%
    mutate(name = paste(location, condition, sep = "_")) %>%
    filter(location != "Static")

## ----✷----✷---- Statistical Tests ----✷----✷----
ks.test(location$fold_change, "pnorm")


location %>%
    group_by(primer) %>%
    t_test(fold_change ~ location, p.adjust.method = "bonferroni") %>%
    add_xy_position(x = "name") %>%
    add_significance("p")

location <- location %>%
    filter(condition == "Control", location != "Periphery") %>%  distinct(avg_ct, .keep_all = TRUE) %>%
    mutate(across(primer, factor, levels = c("AXIN2", "ANGPT2", "THBS1")))

## ----✷----✷---- ANGPT2 plot ----✷----✷----✷----✷----✷----✷----
location_plot <- location %>%
    ggplot(aes(y = fold_change, x = condition)) +
    geom_bar(
        aes(fill = condition),
        show.legend = FALSE,
        position = position_dodge(width = 1),
        stat = "identity",
        colour = "black"
    ) +
    geom_errorbar(
        aes(ymin = fold_change,
            ymax = fold_change + fold_sem),
        position = position_dodge(1),
        size = 0.4,
        width = 0.3
    ) +
    facet_wrap(. ~ primer, strip.position = "bottom") +
    scale_x_discrete(labels = c("", "")) +
    labs(title = "",
         x = "",
         y = "Relative mRNA Expression",
         fill = "") +
    theme_classic() +
    theme +
    scale_fill_grey(start = 0,
                    end = 0) +
    scale_y_continuous(
        expand = c(0, 0),
        breaks = seq(0, 0.6, 0.2),
        limits = c(0, 0.6)
    )










## ----✷----✷---- drug ----✷----✷----✷----✷----✷----✷----✷----✷----
drug <-
    read.csv("../data/raw/dmso.csv",
             header = TRUE) %>%
    janitor::clean_names() %>%
    mutate(name = paste(location, condition, sep = "_")) %>%
    filter(location != "Static")

## ----✷----✷---- Statistical Tests ----✷----✷----
ks.test(drug$fold_change, "pnorm")

drug %>%
    group_by(primer) %>%
    t_test(fold_change ~ condition, p.adjust.method = "bonferroni") %>%
    add_xy_position(x = "name") %>%
    add_significance("p")


drug <- drug %>%
    filter(condition != "Control") %>%  distinct(avg_ct, .keep_all = TRUE) %>%
    mutate(across(primer, factor, levels = c("AXIN2", "ANGPT2", "THBS1")))


## ----✷----✷---- plot ----✷----✷----✷----✷----✷----✷----
drug_plot <- drug %>%
    ggplot(aes(y = fold_change, x = location)) +
    geom_bar(
        aes(fill = location),
        show.legend = FALSE,
        position = position_dodge(width = 1),
        stat = "identity",
        colour = "black"
    ) +
    geom_errorbar(
        aes(ymin = fold_change,
            ymax = fold_change + fold_sem),
        position = position_dodge(1),
        size = 0.4,
        width = 0.3
    ) +
    facet_wrap(. ~ primer, strip.position = "bottom") +
    scale_x_discrete(labels = c("Low", "High")) +
    labs(title = "",
         x = "",
         y = "Relative mRNA Expression",
         fill = "") +
    theme_classic() +
    theme +
    scale_fill_grey(start = 0,
                    end = 1) +
    scale_y_continuous(
        expand = c(0, 0),
        breaks = seq(0, 50, 10),
        limits = c(0, 50)
    )








