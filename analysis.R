library(tidyverse)
library(ggpubr)

results <-
    read.csv("../data/raw/results.csv",
             header = TRUE) %>%
    janitor::clean_names() %>%
    mutate(sem = sd / sqrt(3))


ang_plot <- results %>%
    filter(primer == "Ang-2" & !str_detect(location, "Static")) %>%
    ggplot(aes(fill = condition, y = dd_ct, x = location)) +
    geom_bar(position = position_dodge(width = 1),
             stat = "identity",
             colour = "black") +
    geom_errorbar(
        aes(ymin = dd_ct,
            ymax = dd_ct + sem),
        position = position_dodge(1),
        size = 0.5,
        width = 0.2
    ) +
    scale_x_discrete(labels = c("Low", "High")) +
    labs(
        title = "Ang-2",
        x = "",
        y = expression(~ Delta * ~ Delta * bold("Cq")),
        fill = ""
    ) +
    theme_classic() +
    theme(
        plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
        text = element_text(size = 16),
        axis.text.x = element_text(face = "bold"),
        axis.title.y = element_text(face = "bold")
    ) +
    scale_fill_grey(start = 0,
                    end = 1) +
    scale_y_continuous(
        expand = c(0, 0),
        breaks = seq(0, 10, 2),
        limits = c(0.0, 10.0)
    )

axin_plot <- results %>%
    filter(primer == "Axin" & !str_detect(location, "Static")) %>%
    ggplot(aes(fill = condition, y = dd_ct, x = location)) +
    geom_bar(position = position_dodge(width = 1),
             stat = "identity",
             colour = "black") +
    geom_errorbar(
        aes(ymin = dd_ct,
            ymax = dd_ct + sem),
        position = position_dodge(1),
        size = 0.5,
        width = 0.2
    ) +
    scale_x_discrete(labels = c("Low", "High")) +
    labs(
        title = "Axin",
        x = "",
        y = expression(~ Delta * ~ Delta * bold("Cq")),
        fill = ""
    ) +
    theme_classic() +
    theme(
        plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
        text = element_text(size = 16),
        axis.text.x = element_text(face = "bold"),
        axis.title.y = element_text(face = "bold")
    ) +
    scale_fill_grey(start = 0,
                    end = 1) +
    scale_y_continuous(
        expand = c(0, 0),
        breaks = seq(0, 8, 2),
        limits = c(0.0, 8)
    )


thsb1_plot <- results %>%
    filter(primer == "THSB-1" & !str_detect(location, "Static")) %>%
    ggplot(aes(fill = condition, y = dd_ct, x = location)) +
    geom_bar(position = position_dodge(width = 1),
             stat = "identity",
             colour = "black") +
    geom_errorbar(
        aes(ymin = dd_ct,
            ymax = dd_ct + sem),
        position = position_dodge(1),
        size = 0.5,
        width = 0.2
    ) +
    scale_x_discrete(labels = c("Low", "High")) +
    labs(
        title = "THSB-1",
        x = "",
        y = expression(~ Delta * ~ Delta * bold("Cq")),
        fill = ""
    ) +
    theme_classic() +
    theme(
        plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
        text = element_text(size = 16),
        axis.text.x = element_text(face = "bold"),
        axis.title.y = element_text(face = "bold")
    ) +
    scale_fill_grey(start = 0,
                    end = 1) +
    scale_y_continuous(
        expand = c(0, 0),
        breaks = seq(0, 8, 2),
        limits = c(0.0, 8)
    )


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