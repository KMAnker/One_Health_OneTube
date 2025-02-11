library(tidyverse)
library(dplyr)
library(ggplot2)
library(ggpubr)

setwd("PATH/TO/PROJECT/DIRECTORY")
setwd("/Users/klara_marie/Desktop/One_health_one_tube")

################################################################################

### Read in coverage_all.txt ###
read_coverage_table <- read.table("data/coverage_all.txt", header=TRUE, sep = "\t", na.strings="")

read_coverage <- read_coverage_table %>%
  group_by(host, extraction, OT, primers, segment2, Position) %>%
  mutate(OT_primers = paste(OT, primers, sep= "-"),
         average_coverage_depth = mean(Coverage.Depth)) %>%
  select(sample, host, extraction, OT, primers, OT_primers, segment2, segment, Position, Coverage.Depth, average_coverage_depth) %>%
  ungroup()

# order columns
read_coverage$segment <- factor(read_coverage$segment, levels = c("PB2", "PB1", "PA", "HA", "NP", "NA", "MP", "NS"))
read_coverage$OT_primers <- factor(read_coverage$OT_primers, levels = c("OT1-MBT", "OT1-UNI", "OT2-MBT", "OT2-UNI", "OT1-MBT2"))
read_coverage$host <- factor(read_coverage$host, levels = c("HU", "SW", "AV"))

################################################################################

### Read in ct data ###
sample_ct <- read.table("data/sample_cts.txt", header = TRUE)

################################################################################

### prepare data for plotting ###

#################################

### combine coverage table with ct data ###
summarized_read_coverage <- read_coverage %>%
  left_join(sample_ct, by = c("host", "sample", "extraction"))
summarized_read_coverage$host <- factor(summarized_read_coverage$host, levels = c("HU", "SW", "AV"))

### get all samples tested with extra primer set
extra_primers_samples <- summarized_read_coverage %>%
  group_by(sample) %>%
  filter(any(primers == "MBT2")) %>%
  distinct(sample)

### filter table to represent those samples only
extra_primers_coverage <- summarized_read_coverage %>%
  filter(sample %in% extra_primers_samples$sample)

### filter away samples with ct < 36
extra_primers_coverage_filtered <- extra_primers_coverage %>%
  group_by(sample) %>%
  filter(any(ct < 36)) %>%
  ungroup()

# order categories
extra_primers_coverage_filtered$OT_primers <- 
  factor(extra_primers_coverage_filtered$OT_primers, 
         levels = c("OT1-MBT", "OT1-UNI", "OT2-MBT", "OT2-UNI", "OT1-MBT2"))

### order segments
extra_primers_coverage_filtered$segment2 <- factor(extra_primers_coverage_filtered$segment2, levels = c("PB2", "PB1", "PA", 
                                                                                            "H1", "H2", "H3", "H4", "H5", "H6", "H7", "H8", "H9", "H10", "H11", "H12", "H13", "H14", "H15", "H16", 
                                                                                            "NP", 
                                                                                            "N1", "N2", "N3", "N4", "N5", "N6", "N7", "N8", "N9",
                                                                                            "MP", "NS"))

### seperate by host 
extra_primers_coverage_filtered_hu <- extra_primers_coverage_filtered %>%
  filter(host == "HU")
extra_primers_coverage_filtered_sw <- extra_primers_coverage_filtered %>%
  filter(host == "SW")
extra_primers_coverage_filtered_av <- extra_primers_coverage_filtered %>%
  filter(host == "AV")

### Plotting ###

theme1 <- theme(plot.title = element_text(size = 8, hjust = 0.5, margin = margin(0, 0, 0, 0)),
                legend.position = "top",
                legend.text = element_text(size = 8),
                legend.title = element_blank(),
                legend.key.size = unit(8, "points"),
                axis.text.x=element_blank(),
                axis.ticks.x=element_blank(), 
                axis.text.y = element_text(size = 8),
                axis.title.y = element_text(size = 8),
                strip.text = element_text(size = 8),
                strip.text.y = element_blank(),
                axis.line = element_line(size = 0),
                strip.background = element_rect(size = 0),
                panel.spacing.x = unit(0.2, "lines"),
                panel.spacing.y = unit(0.1, "lines"),
                panel.grid.minor.x = element_blank(),
                panel.grid.major.x = element_blank(),
                legend.box.spacing = unit(-1, "lines"),
                plot.margin = unit(c(0,0.2,0,0.4), "lines"))

theme2 <- theme(plot.title = element_text(size = 8, hjust = 0.5, margin = margin(0, 0, 0, 0)),
                legend.position = "top",
                legend.text = element_text(size = 8),
                legend.title = element_blank(),
                legend.key.size = unit(8, "points"),
                axis.text.x=element_blank(), 
                axis.ticks.x=element_blank(),  
                axis.text.y = element_text(size = 8),
                axis.title.y = element_text(size = 8),
                strip.text = element_text(size = 8),
                axis.line = element_line(size = 0),
                strip.background = element_rect(size = 0),
                panel.spacing.x = unit(0.2, "lines"),
                panel.spacing.y = unit(0.1, "lines"),
                panel.grid.minor.x = element_blank(),
                panel.grid.major.x = element_blank(),
                legend.box.spacing = unit(-1, "lines"),
                plot.margin = unit(c(0,0.2,0,0.4), "lines"))


### FigS4: Coverage per segment extra primers ###

coverage_human_extra <- ggplot() +
  geom_line(data = extra_primers_coverage_filtered_hu, aes(x = Position, y = Coverage.Depth, color = OT_primers, group = interaction(sample, OT_primers)), linewidth = 0.05, alpha = 0.5) +
  geom_line(data = extra_primers_coverage_filtered_hu, aes(x = Position, y = average_coverage_depth, color = OT_primers), linewidth = 0.3) +
  geom_hline(yintercept = 1000, linetype = 2, size = 0.3, alpha = 0.5) +
  labs(title = "Human Samples", x = "", y = "Read depth", color = "") +
  theme_minimal() +
  scale_color_manual(values = c("#e76f51", "#f4a261", "#2a9d8f", "#264653", "#613653")) +
  theme2 +
  facet_grid(extraction~segment2, scales = "free_x", switch = "x", labeller = labeller(host = c(HU = "Human Samples", SW = "Swine Samples", AV = "Avian Samples"))) +
  scale_y_log10()

coverage_swine_extra <- ggplot() +
  geom_line(data = extra_primers_coverage_filtered_sw, aes(x = Position, y = Coverage.Depth, color = OT_primers, group = interaction(sample, OT_primers)), linewidth = 0.05, alpha = 0.5) +
  geom_line(data = extra_primers_coverage_filtered_sw, aes(x = Position, y = average_coverage_depth, color = OT_primers), linewidth = 0.3) +
  geom_hline(yintercept = 1000, linetype = 2, size = 0.3, alpha = 0.5) +
  labs(title = "Swine Samples", x = "", y = "Read depth", color = "") +
  theme_minimal() +
  scale_color_manual(values = c("#e76f51", "#f4a261", "#2a9d8f", "#264653", "#613653")) +
  theme2 +
  facet_grid(extraction~segment2, scales = "free_x", switch = "x", labeller = labeller(host = c(HU = "Human Samples", SW = "Swine Samples", AV = "Avian Samples"))) +
  scale_y_log10()

coverage_avian_extra <- ggplot() +
  geom_line(data = extra_primers_coverage_filtered_av, aes(x = Position, y = Coverage.Depth, color = OT_primers, group = interaction(sample, OT_primers)), linewidth = 0.05, alpha = 0.5) +
  geom_line(data = extra_primers_coverage_filtered_av, aes(x = Position, y = average_coverage_depth, color = OT_primers), linewidth = 0.3) +
  geom_hline(yintercept = 1000, linetype = 2, size = 0.3, alpha = 0.5) +
  labs(title = "Avian Samples", x = "", y = "Read depth", color = "") +
  theme_minimal() +
  scale_color_manual(values = c("#e76f51", "#f4a261", "#2a9d8f", "#264653", "#613653")) +
  theme2 +
  facet_grid(extraction~segment2, scales = "free_x", switch = "x", labeller = labeller(host = c(HU = "Human Samples", SW = "Swine Samples", AV = "Avian Samples"))) +
  scale_y_log10()

figS4 <- ggarrange(coverage_human_extra, coverage_swine_extra, coverage_avian_extra,
                  ncol = 1,
                  nrow = 3,
                  common.legend = TRUE) +
  theme(plot.margin = unit(c(0,0,0,0), "null"))


ggsave("figures/FigureS4.png", plot = figS4, dpi = 600, width = 16 , height = 20, units = "cm", bg = "white")


