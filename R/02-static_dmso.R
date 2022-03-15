library(tidyverse)
library(ggpubr)
library(rstatix)
library(grid)
library(gridExtra)

ymax <-
    function(x) {
        (x$fold_change + (x$fold_change > 0) * x$fold_sem)
    }

ymin <-
    function(x) {
        (x$fold_change - (x$fold_change < 0) * x$fold_sem)
    }


theme <- theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
    text = element_text(size = 14),
    axis.text.x = element_text(hjust = -0.5),
    axis.title.y = element_text(face = "bold"),
    axis.ticks = element_blank(),
    strip.placement = "outside",
    strip.background = element_rect(colour = "white", fill = "white"),
    strip.text.x = element_text(size = 12, face = "bold"),
    panel.spacing = unit(0, "lines"),
    plot.background = element_rect(fill = "transparent", color = NA),
    panel.grid.major.x = element_blank()
)

## ----✷----✷---- results2 ----✷----✷----✷----✷----✷----✷----✷----✷----
results2 <-
    read.csv("../data/raw/static_dmso.csv",
             header = TRUE) %>%
    janitor::clean_names() %>%
    mutate(name = paste(location, condition, sep = "_"))



## ----✷----✷---- Tidy results2  ----✷----✷----✷----✷----✷----✷----
ang <- results2 %>%
    filter(primer == "ANGPT2" & location != "Static")

ang2 <- results2 %>%
    filter(primer == "ANGPT2")


axin <- results2 %>%
    filter(primer == "AXIN2" & location != "Static")
axin2 <- results2 %>%
    filter(primer == "AXIN2")


thbs1 <- results2 %>%
    filter(primer == "THBS1" & location != "Static")
thbs12 <- results2 %>%
    filter(primer == "THBS1")


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
    ang %>%
    dunn_test(fold_change ~ name, p.adjust.method = "bonferroni") %>%
    add_xy_position(x = "name") %>%
    filter(group1 != "Centre_Control" |
               group2 != "Periphery_Drug") %>%
    filter(group1 != "Centre_Drug" | group2 != "Periphery_Control")

ang2_test <-
    ang2 %>% tukey_hsd(fold_change ~ name) %>%
    add_xy_position(x = "name") %>%
    filter(
        group1 == "Centre_Drug" & group2 == "Static_Control" |
            group1 == "Periphery_Drug" &
            group2 == "Static_Control" |
            group1 == "Centre_Control" &
            group2 == "Static_Control" |
            group1 == "Periphery_Control" &
            group2 == "Static_Control"
    )



axin_test <-
    axin %>%
    group_by(location) %>%
    dunn_test(fold_change ~ name, p.adjust.method = "bonferroni") %>%
    add_xy_position(x = "name") %>%
    filter(group1 != "Centre_Control" |
               group2 != "Periphery_Drug") %>%
    filter(group1 != "Centre_Drug" | group2 != "Periphery_Control")

axin2_test <-
    axin2 %>% tukey_hsd(fold_change ~ name) %>%
    add_xy_position(x = "name") %>%
    filter(
        group1 == "Centre_Drug" & group2 == "Static_Control" |
            group1 == "Periphery_Drug" &
            group2 == "Static_Control" |
            group1 == "Centre_Control" &
            group2 == "Static_Control" |
            group1 == "Periphery_Control" &
            group2 == "Static_Control"
    )



thbs1_test <-
    thbs1 %>% dunn_test(fold_change ~ name, p.adjust.method = "bonferroni") %>%
    add_xy_position(x = "name") %>%
    filter(group1 != "Centre_Control" |
               group2 != "Periphery_Drug") %>%
    filter(group1 != "Centre_Drug" | group2 != "Periphery_Control")

thbs12_test <-
    thbs12 %>% tukey_hsd(fold_change ~ name) %>%
    add_xy_position(x = "name") %>%
    filter(
        group1 == "Centre_Drug" & group2 == "Static_Control" |
            group1 == "Periphery_Drug" &
            group2 == "Static_Control" |
            group1 == "Centre_Control" &
            group2 == "Static_Control" |
            group1 == "Periphery_Control" &
            group2 == "Static_Control"
    )


## ----✷----✷---- ANGPT2 plot ----✷----✷----✷----✷----✷----✷----
ang_plot_2 <- ang %>%
    ggplot(aes(y = fold_change, x = name)) +
    geom_bar(
        aes(fill = condition),
        position = position_dodge(width = 1),
        stat = "identity",
        colour = "black"
    ) +
    geom_errorbar(
        aes(ymin = ymin(ang),
            ymax = ymax(ang)),
        position = position_dodge(1),
        size = 0.4,
        width = 0.3
    ) +
    scale_x_discrete(labels = c("Low", "", "High", ""), ) +
    labs(title = "ANGPT2",
         x = "",
         y = "Relative mRNA Expression",
         fill = "") +
    theme_classic() +
    theme +
    scale_fill_grey(start = 1,
                    end = 0,
                    labels = c("DMSO", "XAV939")) +
    scale_y_continuous(
        expand = c(0, 0),
        breaks = seq(0, .75, 0.25),
        limits = c(0, .75)
    ) +
    geom_hline(yintercept = 1, linetype = 3) +
    stat_pvalue_manual(
        ang_test,
        label = "p.adj.signif",
        tip.length = 0,
        hide.ns = TRUE,
        y.position = .65,
        size = 6
    )

ang_plot_2

# + stat_pvalue_manual(
#         ang_test,
#         label = "p.adj.signif",
#         tip.length = 0,
#         hide.ns = TRUE,
#         y.position = c(0.14, 0.13, 0.12),
#         size = 6
#     )

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


## ----✷----✷---- Axin Plot ----✷----✷----✷----✷----✷----✷----
axin_plot_2 <- axin %>%
    ggplot(aes(y = fold_change, x = name)) +
    geom_bar(
        aes(fill = condition),
        
        position = position_dodge(width = 1),
        stat = "identity",
        colour = "black"
    ) +
    geom_errorbar(
        aes(ymin = ymin(axin),
            ymax = ymax(axin)),
        position = position_dodge(1),
        size = 0.4,
        width = 0.3
    ) +
    scale_x_discrete(labels = c("Low", "", "High", "")) +
    labs(title = "AXIN2",
         x = "",
         y = "Relative mRNA Expression",
         fill = "") +
    theme_classic() +
    theme +
    scale_fill_grey(start = 1,
                    end = 0,
                    labels = c("DMSO", "XAV939")) +
    scale_y_continuous(
        expand = c(0, 0),
        breaks = seq(0, 2.5, 0.5),
        limits = c(0, 2.5)
    ) +
    geom_hline(yintercept = 1, linetype = 3) +
    stat_pvalue_manual(
        ang_test,
        label = "p.adj.signif",
        tip.length = 0,
        hide.ns = TRUE,
        y.position = 2.2,
        size = 6
    )

axin_plot_2


# + stat_pvalue_manual(
#         axin_test,
#         label = "p.adj.signif",
#         tip.length = 0,
#         hide.ns = TRUE,
#         y.position = c(0.65, 0.55),
#         size = 6
#     )

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

axin_plot_2

## ----✷----✷---- THBS1 Plot ----✷----✷----✷----✷----✷----✷----
thbs1_plot_2 <- thbs1 %>%
    ggplot(aes(y = fold_change, x = name)) +
    geom_bar(
        aes(fill = condition),
        
        position = position_dodge(width = 1),
        stat = "identity",
        colour = "black"
    ) +
    geom_errorbar(
        aes(ymin = ymin(thbs1),
            ymax = ymax(thbs1)),
        position = position_dodge(1),
        size = 0.4,
        width = 0.4
    ) +
    scale_x_discrete(labels = c("Low", "", "High", "")) +
    labs(title = "THBS1",
         x = "",
         y = "Relative mRNA Expression",
         fill = "") +
    theme_classic() +
    theme +
    scale_fill_grey(start = 1,
                    end = 0,
                    labels = c("DMSO", "XAV939")) +
    scale_y_continuous(
        expand = c(0, 0),
        breaks = seq(0, 2.5, 0.5),
        limits = c(0, 2.5)
    ) +
    geom_hline(yintercept = 1, linetype = 3) +
    stat_pvalue_manual(
        ang_test,
        label = "p.adj.signif",
        tip.length = 0,
        hide.ns = TRUE,
        y.position = 2.4,
        size = 6
    )



# +
#     stat_pvalue_manual(
#         thbs1_test,
#         label = "p.adj.signif",
#         tip.length = 0,
#         hide.ns = TRUE,
#         y.position = c(0.7),
#         size = 6
#     )
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

thbs1_plot_2

## ----✷----✷---- Plot ----✷----✷----✷----✷----✷----✷----

plot2 <- ggarrange(
    axin_plot_2,
    ang_plot_2,
    thbs1_plot_2,
    labels = c("A", "B", "C"),
    font.label = list(face = "bold", size = 20),
    label.x = -0.03,
    label.y = 1.01,
    ncol = 3,
    nrow = 1,
    common.legend = TRUE,
    legend = "bottom"
)
plot2
