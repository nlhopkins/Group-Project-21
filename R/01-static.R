library(tidyverse)
library(ggpubr)
library(rstatix)


## ----✷----✷---- Results ----✷----✷----✷----✷----✷----✷----✷----✷----
results <-
    read.csv("../data/raw/static.csv",
             header = TRUE) %>%
    janitor::clean_names() %>%
    mutate(name = paste(location, condition, sep = "_"))


## ----✷----✷---- Tidy Results  ----✷----✷----✷----✷----✷----✷----
ang <- results %>%
    filter(primer == "ANGPT2" & !str_detect(name, "Static"))

ang2 <- results %>%
    filter(primer == "ANGPT2")


axin <- results %>%
    filter(primer == "AXIN2" & !str_detect(name, "Static"))
axin2 <- results %>%
    filter(primer == "AXIN2")


thbs1 <- results %>%
    filter(primer == "THSB1" & !str_detect(name, "Static"))
thbs12 <- results %>%
    filter(primer == "THSB1")


## ----✷----✷---- Statistical Tests ----✷----✷----
ks_ang <- ang$dd_ct[!duplicated(ang$fold_change)]
ks.test(ks_ang, "pnorm")

ks_axin <- axin$dd_ct[!duplicated(axin$fold_change)]
ks.test(ks_axin, "pnorm")

ks_thbs1 <- thbs1$dd_ct[!duplicated(thbs1$fold_change)]
ks.test(ks_thbs1, "pnorm")


ang %>% kruskal_test(fold_change ~ name)
axin %>% kruskal_test(fold_change ~ name)
thbs1 %>% kruskal_test(fold_change ~ name)


ang_test <-
    ang %>% dunn_test(dd_ct ~ name, p.adjust.method = "none") %>%
    add_xy_position(x = "name") %>%
    filter(group1 != "Centre_Control" |
               group2 != "Periphery_Drug") %>%
    filter(group1 != "Centre_Drug" | group2 != "Periphery_Control")

ang2_test <-
    ang2 %>% tukey_hsd(dd_ct ~ name) %>%
    add_xy_position(x = "name") %>%
    filter(
        group1 == "Centre_Drug" & group2 == "Static_Drug" |
            group1 == "Periphery_Drug" & group2 == "Static_Drug" |
            group1 == "Centre_Control" &
            group2 == "Static_Control" |
            group1 == "Periphery_Control" &
            group2 == "Static_Control"
    )



axin_test <-
    axin %>% dunn_test(dd_ct ~ name, p.adjust.method = "none") %>%
    add_xy_position(x = "name") %>%
    filter(group1 != "Centre_Control" |
               group2 != "Periphery_Drug") %>%
    filter(group1 != "Centre_Drug" | group2 != "Periphery_Control")

axin2_test <-
    axin2 %>% tukey_hsd(dd_ct ~ name) %>%
    add_xy_position(x = "name") %>%
    filter(
        group1 == "Centre_Drug" & group2 == "Static_Drug" |
            group1 == "Periphery_Drug" & group2 == "Static_Drug" |
            group1 == "Centre_Control" &
            group2 == "Static_Control" |
            group1 == "Periphery_Control" &
            group2 == "Static_Control"
    )



thbs1_test <-
    thbs1 %>% dunn_test(dd_ct ~ name, p.adjust.method = "none") %>%
    add_xy_position(x = "name") %>%
    filter(group1 != "Centre_Control" |
               group2 != "Periphery_Drug") %>%
    filter(group1 != "Centre_Drug" | group2 != "Periphery_Control")

thbs12_test <-
    thbs12 %>% tukey_hsd(dd_ct ~ name) %>%
    add_xy_position(x = "name") %>%
    filter(
        group1 == "Centre_Drug" & group2 == "Static_Drug" |
            group1 == "Periphery_Drug" & group2 == "Static_Drug" |
            group1 == "Centre_Control" &
            group2 == "Static_Control" |
            group1 == "Periphery_Control" &
            group2 == "Static_Control"
    )


## ----✷----✷---- ANGPT2 plot ----✷----✷----✷----✷----✷----✷----
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
        y.position = c(0.14, 0.13, 0.12),
        size = 6
    ) 

# + stat_pvalue_manual(
#         ang2_test,
#         y.position = c((0.0553 +
#                             0.033933692 + 0.0005), (0.0067
#                                                     + 0.002183904 + 0.0005)),
#         label = "p.adj.signif",
#         remove.bracket = TRUE,
#         hide.ns = TRUE,
#         color = "red",
#         x =  "xmin",
#         size = 6
#     )
ang_plot


## ----✷----✷---- Axin Plot ----✷----✷----✷----✷----✷----✷----
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
    scale_x_discrete(labels = c("Low", "", "High", "")) +
    labs(title = "AXIN2",
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
        breaks = seq(0, 0.7, 0.1),
        limits = c(0.0, 0.75)
    ) + stat_pvalue_manual(
        axin_test,
        label = "p.adj.signif",
        tip.length = 0,
        hide.ns = TRUE,
        y.position = c(0.65, 0.6, 0.55),
        size = 6
    ) 
# + stat_pvalue_manual(
#         axin2_test,
#         y.position = c((0.2323 + 0.03291077 + 0.005), (0.0297 + 0.02158658 + 0.005)),
#         label = "p.adj.signif",
#         remove.bracket = TRUE,
#         hide.ns = TRUE,
#         color = "red",
#         x =  "group1",
#         size = 6
#     )

axin_plot

## ----✷----✷---- THSB1 Plot ----✷----✷----✷----✷----✷----✷----
thbs1_plot <- thbs1 %>%
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
    labs(title = "THBS1",
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
        breaks = seq(0, 0.7, 0.2),
        limits = c(0.0, 0.75)
    ) +
    stat_pvalue_manual(
        thbs1_test,
        label = "p.adj.signif",
        tip.length = 0,
        hide.ns = TRUE,
        y.position = c(0.7, 0.65, 0.6),
        size = 6
    ) 
# + stat_pvalue_manual(
#         thbs12_test,
#         y.position = c((0.1246 + 0.015852085 + 0.005), (0.0263 + 0.009946131 + 0.005)),
#         label = "p.adj.signif",
#         remove.bracket = TRUE,
#         hide.ns = TRUE,
#         color = "red",
#         x =  "group1",
#         size = 6
#     )

thbs1_plot

## ----✷----✷---- Plot ----✷----✷----✷----✷----✷----✷----

plot1 <- ggarrange(
    ang_plot,
    axin_plot,
    thbs1_plot,
    labels = c("A", "B", "C"),
    font.label = list(face = "bold", size = 20),
    label.x = -0.03,
    label.y = 1.01,
    ncol = 3,
    nrow = 1,
    common.legend = TRUE,
    legend = "bottom"
)
plot1
