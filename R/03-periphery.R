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



## ----✷----✷---- results3 ----✷----✷----✷----✷----✷----✷----✷----✷----
results3 <-
    read.csv("../data/raw/periphery.csv",
             header = TRUE) %>%
    janitor::clean_names() %>%
    mutate(name = paste(location, condition, sep = "_"))



## ----✷----✷---- Tidy results3  ----✷----✷----✷----✷----✷----✷----
ang <- results3 %>%
    filter(primer == "ANGPT2" & location != "Periphery",
               primer == "ANGPT2" & location != "Static")

ang2 <- results3 %>%
    filter(primer == "ANGPT2" & location != "Static")


axin <- results3 %>%
    filter(primer == "AXIN2" & location != "Periphery",
               primer == "AXIN2" & location != "Static")

axin2 <- results3 %>%
    filter(primer == "AXIN2" & location != "Static")


thbs1 <- results3 %>%
    filter(primer == "THBS1" & location != "Periphery",
               primer == "THBS1" & location != "Static")
thbs12 <- results3 %>%
    filter(primer == "THBS1" & location != "Static")





## ----✷----✷---- Statistical Tests ----✷----✷----
ks_ang <- ang$dd_ct[!duplicated(ang$name)]
ks.test(ks_ang, "pnorm")

ks_axin <- axin$dd_ct[!duplicated(axin$name)]
ks.test(ks_axin, "pnorm")

ks_thbs1 <- thbs1$dd_ct[!duplicated(thbs1$name)]
ks.test(ks_thbs1, "pnorm")


ang2 %>% wilcox_test(fold_change ~ name, paired = TRUE)
axin2 %>% wilcox_test(fold_change ~ name, paired = TRUE)
thbs12 %>% wilcox_test(fold_change ~ name, paired = TRUE)

## ----✷----✷---- ANGPT2 plot ----✷----✷----✷----✷----✷----✷----
ang_plot_3 <- ang %>%
    ggplot(aes(y = fold_change, x = condition)) +
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
    scale_x_discrete(labels = c("-", "+")) +
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
        breaks = seq(0, 20, 5),
        limits = c(0, 20)
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
    ggplot(aes(y = fold_change, x = condition)) +
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
    scale_x_discrete(labels = c("-", "+")) +
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
        breaks = seq(0, 10, 2),
        limits = c(0, 10)
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
thbs1_plot_3 <- thbs1 %>%
    ggplot(aes(y = fold_change, x = condition)) +
    geom_bar(
        aes(fill = condition),
        show.legend = FALSE,
        position = position_dodge(width = 1),
        stat = "identity",
        colour = "black"
    ) +
    geom_errorbar(
        aes(ymin = ymin(thbs1),
            ymax = ymax(thbs1)),
        position = position_dodge(1),
        size = 0.4,
        width = 0.3
    ) +
    scale_x_discrete(labels = c("-", "+")) +
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
        breaks = seq(0, 6, 1),
        limits = c(0, 6)
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

thbs1_plot_3

## ----✷----✷---- Plot ----✷----✷----✷----✷----✷----✷----

plot3 <- ggarrange(
    axin_plot_3,
    ang_plot_3,
    thbs1_plot_3,
    labels = c("D", "E", "F"),
    font.label = list(face = "bold", size = 20),
    label.x = -0.03,
    label.y = 1.01,
    ncol = 3,
    nrow = 1
)
plot3
