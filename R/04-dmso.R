library(tidyverse)
library(ggpubr)
library(rstatix)


ymax <-
    function(x) {
        (x$fold_change + (x$fold_change > 0) * x$fold_sem)
    }

ymin <-
    function(x) {
        (x$fold_change - (x$fold_change < 0) * x$fold_sem)
    }



## ----✷----✷---- results4 ----✷----✷----✷----✷----✷----✷----✷----✷----
results4 <-
    read.csv("../data/raw/dmso.csv",
             header = TRUE) %>%
    janitor::clean_names() %>%
    mutate(name = paste(location, condition, sep = "_"))



## ----✷----✷---- Tidy results4  ----✷----✷----✷----✷----✷----✷----
ang <- results4 %>%
    filter(primer == "ANGPT2" &
               condition == "Drug", location != "Static")

ang2 <- results4 %>%
    filter(primer == "ANGPT2")


axin <- results4 %>%
    filter(primer == "AXIN2" &
               condition == "Drug", location != "Static")
axin2 <- results4 %>%
    filter(primer == "AXIN2")


THBS1 <- results4 %>%
    filter(primer == "THBS1" &
               condition == "Drug", location != "Static")
THBS12 <- results4 %>%
    filter(primer == "THBS1")


## ----✷----✷---- Statistical Tests ----✷----✷----
ks_ang <- ang$dd_ct[!duplicated(ang$fold_change)]
ks.test(ks_ang, "pnorm")

ks_axin <- axin$dd_ct[!duplicated(axin$fold_change)]
ks.test(ks_axin, "pnorm")

ks_THBS1 <- THBS1$dd_ct[!duplicated(THBS1$fold_change)]
ks.test(ks_THBS1, "pnorm")


ang %>% kruskal_test(fold_change ~ name)
axin %>% kruskal_test(fold_change ~ name)
THBS1 %>% kruskal_test(fold_change ~ name)

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
        group1 == "Centre_Control" & group2 == "Periphery_Control" |
            group1 == "Centre_Drug" &
            group2 == "Periphery_Drug"
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
        group1 == "Centre_Control" & group2 == "Periphery_Control" |
            group1 == "Centre_Drug" &
            group2 == "Periphery_Drug"
    )



THBS1_test <-
    THBS1 %>% dunn_test(dd_ct ~ name, p.adjust.method = "none") %>%
    add_xy_position(x = "name") %>%
    filter(group1 != "Centre_Control" |
               group2 != "Periphery_Drug") %>%
    filter(group1 != "Centre_Drug" | group2 != "Periphery_Control")

THBS12_test <-
    THBS12 %>% tukey_hsd(dd_ct ~ name) %>%
    add_xy_position(x = "name") %>%
    filter(
        group1 == "Centre_Control" & group2 == "Periphery_Control" |
            group1 == "Centre_Drug" &
            group2 == "Periphery_Drug"
    )


## ----✷----✷---- ANGPT2 plot ----✷----✷----✷----✷----✷----✷----
ang_plot_3 <- ang %>%
    ggplot(aes(y = fold_change, x = name)) +
    geom_bar(
        aes(fill = condition),
        show.legend = FALSE,
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
    scale_x_discrete(labels = c("Centre", "Periphery")) +
    labs(title = "ANGPT2",
         x = "",
         y = "Relative mRNA Expression",
         fill = "") +
    theme_classic() +
    theme(
        plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
        text = element_text(size = 14),
        axis.text.x = element_text(face = "bold"),
        axis.title.y = element_text(face = "bold"),
        axis.ticks = element_blank()
    ) +
    scale_fill_grey(start = 1,
                    end = 0) +
    scale_y_continuous(
        expand = c(0, 0),
        breaks = seq(0, 50, 10),
        limits = c(0, 50)
    ) +
    geom_hline(yintercept = 1, linetype = 3)

ang_plot_3

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
axin_plot_3 <- axin %>%
    ggplot(aes(y = fold_change, x = name)) +
    geom_bar(
        aes(fill = condition),
        show.legend = FALSE,
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
    scale_x_discrete(labels = c("Centre", "Periphery")) +
    labs(title = "AXIN2",
         x = "",
         y = "Relative mRNA Expression",
         fill = "") +
    theme_classic() +
    theme(
        plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
        text = element_text(size = 14),
        axis.text.x = element_text(face = "bold"),
        axis.title.y = element_text(face = "bold"),
        axis.ticks = element_blank()
    ) +
    scale_fill_grey(start = 1,
                    end = 0) +
    scale_y_continuous(
        expand = c(0, 0),
        breaks = seq(0, 50, 10),
        limits = c(0, 50)
    ) +
    geom_hline(yintercept = 1, linetype = 3)

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

axin_plot_3

## ----✷----✷---- THSB1 Plot ----✷----✷----✷----✷----✷----✷----
THBS1_plot_3 <- THBS1 %>%
    ggplot(aes(y = fold_change, x = name)) +
    geom_bar(
        aes(fill = condition),
        show.legend = FALSE,
        position = position_dodge(width = 1),
        stat = "identity",
        colour = "black"
    ) +
    geom_errorbar(
        aes(ymin = ymin(THBS1),
            ymax = ymax(THBS1)),
        position = position_dodge(1),
        size = 0.4,
        width = 0.3
    ) +
    scale_x_discrete(labels = c("Centre", "Periphery")) +
    labs(title = "THBS1",
         x = "",
         y = "Relative mRNA Expression",
         fill = "") +
    theme_classic() +
    theme(
        plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
        text = element_text(size = 14),
        axis.text.x = element_text(face = "bold"),
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
        breaks = seq(0, 50, 10),
        limits = c(0, 50)
    ) +
    geom_hline(yintercept = 1, linetype = 3)

# +
#     stat_pvalue_manual(
#         THBS1_test,
#         label = "p.adj.signif",
#         tip.length = 0,
#         hide.ns = TRUE,
#         y.position = c(0.7),
#         size = 6
#     )
# + stat_pvalue_manual(
#         THBS12_test,
#         y.position = c((0.1246 + 0.015852085 + 0.005), (0.0263 + 0.009946131 + 0.005)),
#         label = "p.adj.signif",
#         remove.bracket = TRUE,
#         hide.ns = TRUE,
#         color = "red",
#         x =  "group1",
#         size = 6
#     )

THBS1_plot_3

## ----✷----✷---- Plot ----✷----✷----✷----✷----✷----✷----

plot4 <- ggarrange(
    ang_plot_3,
    axin_plot_3,
    THBS1_plot_3,
    labels = c("D", "E", "F"),
    font.label = list(face = "bold", size = 20),
    label.x = -0.03,
    label.y = 1.01,
    ncol = 3,
    nrow = 1
)
plot4
