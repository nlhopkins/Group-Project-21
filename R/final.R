library(tidyverse)
library(ggpubr)
library(rstatix)

theme <- theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    text = element_text(size = 14),
    axis.text.x = element_text(hjust = 0.5, size = 14),
    axis.text.y = element_text(size = 14),
    axis.title.y = element_text(face = "bold", size = 14),
    axis.ticks = element_blank(),
    strip.placement = "outside",
    strip.background = element_rect(colour = "white", fill = "white"),
    strip.text.x = element_text(size = 16, face = "bold"),
    panel.spacing = unit(0, "lines"),
    plot.background = element_rect(fill = "transparent", color = NA),
    panel.grid.major.x = element_blank()
)
## ----✷----✷---- dmso ----✷----✷----✷----✷----✷----✷----✷----✷----
dmso <-
    read.csv("../data/raw/dmso.csv",
             header = TRUE) %>%
    janitor::clean_names() %>%
    mutate(name = paste(location, condition, sep = "_")) %>%
    filter(location != "Static")

## ----✷----✷---- Statistical Tests ----✷----✷----
ks.test(dmso$fold_change, "pnorm")


dmso %>%
    group_by(primer) %>%
    t_test(fold_change ~ location, p.adjust.method = "bonferroni") %>%
    add_xy_position(x = "name") %>%
    add_significance("p")

dmso <- dmso %>%
    filter(condition != "Control") %>%  distinct(fold_change, .keep_all = TRUE) %>%
    mutate(across(primer, factor, levels = c("AXIN2", "ANGPT2", "THBS1")))

## ----✷----✷---- ANGPT2 plot ----✷----✷----✷----✷----✷----✷----
dmso_plot <- dmso %>%
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
    scale_fill_grey(start = 1,
                    end = 0) +
    scale_y_continuous(
        expand = c(0, 0),
        breaks = seq(0, 50, 10),
        limits = c(0, 50)
    ) +
    geom_hline(yintercept = 1, linetype = 3)













## ----✷----✷---- periphery ----✷----✷----✷----✷----✷----✷----✷----✷----
periphery <-
    read.csv("../data/raw/periphery.csv",
             header = TRUE) %>%
    janitor::clean_names() %>%
    mutate(name = paste(location, condition, sep = "_")) %>%
    filter(location != "Static")

## ----✷----✷---- Statistical Tests ----✷----✷----
ks.test(periphery$fold_change, "pnorm")

periphery %>%
    group_by(primer) %>%
    t_test(fold_change ~ condition, p.adjust.method = "bonferroni") %>%
    add_xy_position(x = "name") %>%
    add_significance("p")


periphery <- periphery %>%
    filter(location == "Centre") %>%  distinct(fold_change, .keep_all = TRUE) %>%
    mutate(across(primer, factor, levels = c("AXIN2", "ANGPT2", "THBS1")))


## ----✷----✷---- plot ----✷----✷----✷----✷----✷----✷----
periphery_plot <- periphery %>%
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
    scale_x_discrete(labels = c("-", "+")) +
    labs(title = "",
         x = "",
         y = "Relative mRNA Expression",
         fill = "") +
    theme_classic() +
    theme +
    scale_fill_grey(start = 1,
                    end = 0) +
    scale_y_continuous(
        expand = c(0, 0),
        breaks = seq(0, 20, 5),
        limits = c(0, 20)
    ) +
    geom_hline(yintercept = 1, linetype = 3)
