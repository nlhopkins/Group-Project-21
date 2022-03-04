library(tidyverse)
library(ggpubr)
library(rstatix)


# Results ####
results <-
    read.csv("../data/raw/results.csv",
             header = TRUE) %>%
    janitor::clean_names() %>%
    mutate(fold_sem = ((2 ^ -(dd_ct - sem_dd_ct)) - fold_change) - (fold_change - (2 ^ -(dd_ct + sem_dd_ct)))) %>%
    mutate(name = paste(location, condition, sep = "_"))


# Tidy Results ####
ang <- results %>%
    filter(primer == "ANGPT2" & !str_detect(name, "Static"))

axin <- results %>%
    filter(primer == "AXIN2" & !str_detect(name, "Static"))

thsb1 <- results %>%
    filter(primer == "THSB1" & !str_detect(name, "Static"))



# Statistical Tests ####
ang_test <-
    ang %>% dunn_test(fold_change ~ name, p.adjust.method = "none") %>%
    add_xy_position(x = "name") %>%
    filter(group1 != "Centre_Control" |
               group2 != "Periphery_Drug") %>%
    filter(group1 != "Centre_Drug" | group2 != "Periphery_Control")


axin_test <-
    axin %>% dunn_test(fold_change ~ name, p.adjust.method = "none") %>%
    add_xy_position(x = "name") %>%
    filter(group1 != "Centre_Control" |
               group2 != "Periphery_Drug") %>%
    filter(group1 != "Centre_Drug" | group2 != "Periphery_Control")
thsb1_test <-
    thsb1 %>% dunn_test(fold_change ~ name, p.adjust.method = "none") %>%
    add_xy_position(x = "name") %>%
    filter(group1 != "Centre_Control" |
               group2 != "Periphery_Drug") %>%
    filter(group1 != "Centre_Drug" | group2 != "Periphery_Control")



# ANGPT2 plot ####
ang_plot <- ang %>%
    ggplot(aes(y = fold_change, x = name)) +
    geom_bar(
        aes(fill = condition),
        position = position_dodge(width = 1),
        stat = "identity",
        colour = "black"
    ) +
    geom_errorbar(
        aes(ymin = fold_change,
            ymax = fold_change + fold_sem),
        position = position_dodge(1),
        size = 0.5,
        width = 0.2
    ) +
    scale_x_discrete(labels = c("Low", "", "High", "")) +
    labs(title = "ANGPT2",
         x = "",
         y = "Relative mRNA Expression",
         fill = "") +
    theme_classic() +
    theme(
        plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
        text = element_text(size = 14),
        axis.text.x = element_text(face = "bold", hjust = -0.5),
        axis.title.y = element_text(face = "bold"),
        axis.ticks = element_blank()
    ) +
    scale_fill_grey(
        start = 1,
        end = 0,
        labels = c("Control", "XAV939")
    ) +
    scale_y_continuous(
        expand = c(0, 0),
        breaks = seq(0, 0.15, 0.04),
        limits = c(0.0, 0.15)
    ) + stat_pvalue_manual(
        ang_test,
        label = "p.adj.signif",
        tip.length = 0,
        hide.ns = TRUE,
        y.position = c(0.14, 0.13, 0.12)
    )

ang_plot


# Axin Plot ####
axin_plot <- axin %>%
    ggplot(aes(y = fold_change, x = name)) +
    geom_bar(
        aes(fill = condition),
        position = position_dodge(width = 1),
        stat = "identity",
        colour = "black"
    ) +
    geom_errorbar(
        aes(ymin = fold_change,
            ymax = fold_change + fold_sem),
        position = position_dodge(1),
        size = 0.5,
        width = 0.2
    ) +
    scale_x_discrete(labels = c("Low","","High","")) +
    labs(title = "AXIN2",
         x = "",
         y = "Relative mRNA Expression",
         fill = "") +
    theme_classic() +
    theme(
        plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
        text = element_text(size = 14),
        axis.text.x = element_text(face = "bold", hjust = -0.9),
        axis.title.y = element_text(face = "bold"),
        axis.ticks = element_blank()) +
    scale_fill_grey(
        start = 1,
        end = 0,
        labels = c("Control", "XAV939")
    ) +
    scale_y_continuous(
        expand = c(0, 0),
        breaks = seq(0, 0.7, 0.2),
        limits = c(0.0, 0.7)
    ) + stat_pvalue_manual(
        axin_test,
        label = "p.adj.signif",
        tip.length = 0,
        hide.ns = TRUE,
        y.position = c(0.65, 0.55)
    )
axin_plot

# THSB1 Plot ####
thsb1_plot <- thsb1 %>%
    ggplot(aes(y = fold_change, x = name)) +
    geom_bar(
        aes(fill = condition),
        position = position_dodge(width = 1),
        stat = "identity",
        colour = "black"
    ) +
    geom_errorbar(
        aes(ymin = fold_change,
            ymax = fold_change + fold_sem),
        position = position_dodge(1),
        size = 0.5,
        width = 0.2
    ) +
    scale_x_discrete(labels = c("Low", "", "High", "")) +
    labs(title = "THSB1",
         x = "",
         y = "Relative mRNA Expression",
         fill = "") +
    theme_classic() +
    theme(
        plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
        text = element_text(size = 14),
        axis.text.x = element_text(face = "bold", hjust = -0.9),
        axis.title.y = element_text(face = "bold"),
        axis.ticks = element_blank()
    ) +
    scale_fill_grey(
        start = 1,
        end = 0,
        labels = c("Control", "XAV939")
    ) +
    scale_y_continuous(
        expand = c(0, 0),
        breaks = seq(0, 0.8, 0.2),
        limits = c(0.0, 0.8)
    ) +
    stat_pvalue_manual(
        thsb1_test,
        label = "p.adj.signif",
        tip.length = 0,
        hide.ns = TRUE,
        y.position = c(0.7)
    )

thsb1_plot

# plot ####

plot <- ggarrange(
    ang_plot,
    axin_plot,
    thsb1_plot,
    labels = c("A", "B", "C"),
    font.label = list(face = "bold", size = 20),
    label.x = -0.03,
    label.y = 1.01,
    ncol = 3,
    nrow = 1,
    common.legend = TRUE,
    legend = "bottom"
)
plot
