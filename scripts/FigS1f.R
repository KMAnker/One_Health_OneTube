library(tidyverse)
library(dplyr)
library(ggplot2)
library(viridis)

setwd("PATH/TO/PROJECT/DIRECTORY")
gel_data <- read.table("data/OT1_OT2_OT5_gel_scores.txt", header=TRUE)


plotting_data <- gel_data %>%
  mutate(OT_primers = paste(OT, primers, sep= "-"))


gel_scores <- ggplot(plotting_data, aes(x = EX, y = score, fill = OT_primers)) +
  geom_boxplot(position = position_dodge(width = 1), alpha = 0.9, size = 0.3, color = "grey40", outlier.size = 0.5) +
  scale_fill_manual(values = c("#e76f51", "#f4a261", "#2a9d8f", "#264653", "#b8af23", "#808000")) +
  labs(x = "", y = "Score") +
  theme_classic() +
  theme(legend.position = "top") + 
  facet_grid(rows = vars(host)) +
  theme_classic() +
  theme(plot.title = element_blank(),
        legend.position = "top",
        legend.text = element_text(size = 8),
        legend.title = element_blank(),
        legend.key.size = unit(8, "points"),
        axis.text.x = element_text(size = 8, vjust = 0.5),
        axis.text.y = element_text(size = 8),
        axis.title.y = element_text(size = 8),
        strip.text.y = element_text(size = 8),
        strip.text.x = element_text(size = 0),
        axis.line = element_line(size = 0.5),
        axis.ticks = element_line(size = 0.5),
        strip.background = element_rect(size = 0.6), 
        panel.spacing.y = unit(0.5, "lines"),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line(size = 0.5, linetype = 2),
        legend.box.spacing = unit(0, "lines"),
        plot.margin = unit(c(0.1,0.2,0.1,0.2), "lines"))

ggsave("figures/figureS1f.png", plot = gel_scores, dpi = 300, width = 9, height = 7, units = "cm" )



