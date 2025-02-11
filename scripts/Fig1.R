library(tidyverse)
library(dplyr)
library(ggplot2)
library(ggpubr)

setwd("PATH/TO/PROJECT/DIRECTORY")

################################################################################

### Read in READ_COUNTS_all.txt ###
read_counts_table <- read.table("data/READ_COUNTS_all.txt", header=TRUE, sep = "\t")

# Function to extract sample, method, OT, and primers from the Record column
extract_sample_method <- function(Record) {
  record <- str_extract(Record, "(?<=_).*")
  sample_method <- str_extract(Record, "^[^_]+")
  parts_sample_method <- unlist(strsplit(sample_method, "-", fixed = TRUE))
  # Assign values or NA for missing parts
  host <- paste(parts_sample_method[1])
  sample <- paste(parts_sample_method[1], parts_sample_method[2], sep = "-")
  extraction <- ifelse(length(parts_sample_method) >= 3, parts_sample_method[3], NA)
  OT <- ifelse(length(parts_sample_method) >= 4, parts_sample_method[4], NA)
  primers <- ifelse(length(parts_sample_method) >= 5, parts_sample_method[5], NA)
  return(data.frame(host = host, sample = sample, extraction = extraction, OT = OT, primers = primers, record = record))
}

# Apply the function to create the new columns
read_counts <- read_counts_table %>%
  rowwise() %>%
  do({
    extracted_cols <- extract_sample_method(.$Record)
    data.frame(extracted_cols, .)
  }) %>%
  ungroup() %>%
  select(-Record)
read_counts$Reads <- as.numeric(read_counts$Reads)

# filter the table to get data from first extraction test only
read_counts_extraction <- read_counts %>%
  filter(sample %in% c("AV-w", "AV-m", "AV-s", "SW-w", "SW-m1", "SW-m2", "SW-s", "HU-w", "HU-m", "HU-s" ))

# filter the table to get total read no. per sample, no. of reads passing filters and no. of matched reads
filtered_read_counts <- read_counts_extraction %>%
  filter(record %in% c("1-initial", "3-match")) %>%
  mutate(record = case_when(
    record == "1-initial" ~ "Total",
    record == "3-match" ~ "Matched")) %>%
  group_by(host, sample, extraction) %>%
  summarize(
    total_reads = sum(ifelse(record == "Total", as.numeric(Reads), 0)),
    reads_matched = sum(ifelse(record == "Matched", as.numeric(Reads), 0)),
    percent_matched = (reads_matched / total_reads) * 100
  ) %>%
  ungroup()

################################################################################

### Read in passing segments data ###
passing_segments <- read.table("data/passing_segments_all.txt", header=TRUE)

# filter the table to get data from first extraction test only
passing_segments_extraction <- passing_segments %>%
  filter(Sample %in% c("AV-w", "AV-m", "AV-s", "SW-w", "SW-m1", "SW-m2", "SW-s", "HU-w", "HU-m", "HU-s" )) %>%
  rename(sample = Sample, extraction = Method)

  
################################################################################

### combine filtered read counts table with passing segments data ###
filtered_read_counts_segments <- filtered_read_counts %>%
  left_join(passing_segments_extraction, by = c("sample", "extraction"))
filtered_read_counts_segments$host <- factor(filtered_read_counts_segments$host, levels = c("HU", "SW", "AV"))


################################################################################

### Read in coverage_all.txt for plotting ###
read_coverage_table <- read.table("data/coverage_all.txt", header=TRUE, sep = "\t", na.strings="")

# filter the table to get data from first extraction test only and add average depth column
read_coverage_extraction <- read_coverage_table %>%
  filter(sample %in% c("AV-w", "AV-m", "AV-s", "SW-w", "SW-m1", "SW-m2", "SW-s", "HU-w", "HU-m", "HU-s" )) %>%
  group_by(host, extraction, segment2, Position) %>%
  mutate(average_coverage_depth = mean(Coverage.Depth)) %>%
  select(sample, host, extraction, segment2, segment, Position, Coverage.Depth, average_coverage_depth)
read_coverage_extraction$segment <- factor(read_coverage_extraction$segment, levels = c("PB2", "PB1", "PA", "HA", "NP", "NA", "MP", "NS"))
read_coverage_extraction$sample <- factor(read_coverage_extraction$sample, levels = c("AV-w", "AV-m", "AV-s", "SW-w", "SW-m1", "SW-m2", "SW-s", "HU-w", "HU-m", "HU-s"))



################################################################################

### Make summaries of means and standard error

summary_total <- filtered_read_counts_segments %>%
  group_by(extraction) %>%
  summarize(mean_total_reads = mean(total_reads),
            sd_total_reads = sd(total_reads),
            n = n(),
            se_total_reads = sd(total_reads) / sqrt(n()), .groups = 'drop')

summary_matched <- filtered_read_counts_segments %>%
  group_by(extraction) %>%
  summarize(mean_percent_matched = mean(percent_matched),
            sd_percent_matched = sd(percent_matched),
            n = n(),
            se_percent_matched = sd(percent_matched) / sqrt(n()), .groups = 'drop')

summary_passing_segments <- filtered_read_counts_segments %>%
  group_by(extraction) %>%
  summarize(mean_passing_segments = mean(passing_segments),
            sd_passing_segments = sd(passing_segments),
            n = n(),
            se_passing_segments = sd(passing_segments) / sqrt(n()), .groups = 'drop')


### Plotting ###

theme1 <-  theme(legend.position = "bottom",
                legend.text = element_text(size = 8),
                legend.key.size = unit(8, "points"),
                legend.title = element_blank(),
                axis.text.x = element_text(size = 8, vjust = 0.5),
                axis.ticks.x=element_blank(),
                axis.text.y = element_text(size = 8),
                axis.title.y = element_text(size = 8),
                strip.text.y = element_text(size = 8),
                strip.text.x = element_text(size = 8),
                axis.line = element_line(linewidth = 0.5),
                axis.ticks = element_line(linewidth = 0.5),
                strip.background = element_rect(size = 0.6),
                panel.spacing.x = unit(0.2, "lines"), 
                panel.spacing.y = unit(0.1, "lines"),
                panel.grid.minor.x = element_blank(),
                panel.grid.major.x = element_blank(),
                legend.box.spacing = unit(-1, "lines"),
                plot.margin = unit(c(0,0.2,0,0.2), "lines"))

theme2 <- theme(legend.position = "bottom",
                legend.text = element_text(size = 8),
                legend.key.size = unit(8, "points"),
                legend.title = element_blank(),
                axis.text.x=element_blank(),
                axis.ticks.x=element_blank(),
                axis.text.y = element_text(size = 8),
                axis.title.y = element_text(size = 8),
                strip.text.y = element_text(size = 8),
                strip.text.x = element_text(size = 8),
                axis.line = element_line(linewidth = 0),
                axis.ticks = element_line(linewidth = 0),
                strip.background = element_rect(size = 0),
                panel.spacing.x = unit(0.3, "lines"), 
                panel.spacing.y = unit(0.8, "lines"),
                panel.grid.minor.x = element_blank(),
                panel.grid.major.x = element_blank(),
                legend.box.spacing = unit(-1, "lines"),
                plot.margin = unit(c(0,0.2,0,0.2), "lines"))



total_reads_bar_plot <- ggplot(summary_total, aes(x = extraction, y = mean_total_reads, fill = extraction)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.9), alpha = 0.9, color = "grey40", size = 0.4) +
  geom_errorbar(aes(ymin = mean_total_reads - se_total_reads, ymax = mean_total_reads + se_total_reads), 
                position = position_dodge(width = 0.9), width = 0.1, size = 0.4, color = "grey20") +
  geom_point(data = filtered_read_counts_segments, aes(x = extraction, y = total_reads, fill = extraction), 
             position = position_dodge(width = 0.9), shape = 21, size = 1, alpha = 0.1) +
  labs(title = "", 
       y = "Mean Read Count",
       x = "",
       fill = "") +
  theme_classic() +
  theme1 +
  scale_fill_manual(values = c("#287271", "#e76f51"))+
  scale_y_log10() +
  geom_text(aes(label = n, y = max(filtered_read_counts_segments$total_reads)+500000), 
            position = position_dodge(width = 0.9), size = 2)


mapped_reads_bar_plot <- ggplot(summary_matched, aes(x = extraction, y = mean_percent_matched, fill = extraction)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.9), alpha = 0.9, color = "grey40", size = 0.4) +
  geom_errorbar(aes(ymin = mean_percent_matched - se_percent_matched, ymax = mean_percent_matched + se_percent_matched), 
                position = position_dodge(width = 0.9), width = 0.1, size = 0.4, color = "grey20") +
  geom_point(data = filtered_read_counts_segments, aes(x = extraction, y = percent_matched, fill = extraction), 
             position = position_dodge(width = 0.9), shape = 21, size = 1, alpha = 0.1) +
  labs(title = "", 
       y = "% Influenza Reads",
       x = "",
       fill = "") +
  theme_classic() +
  theme1 +
  scale_fill_manual(values = c("#287271", "#e76f51"))+
  geom_text(aes(label = n, y = 105), 
            position = position_dodge(width = 0.9), size = 2)

passing_segments_bar_plot <- ggplot(summary_passing_segments, aes(x = extraction, y = mean_passing_segments, fill = extraction)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.9), alpha = 0.9, color = "grey40", size = 0.4) +
  geom_errorbar(aes(ymin = mean_passing_segments - se_passing_segments, ymax = mean_passing_segments + se_passing_segments), 
                position = position_dodge(width = 0.9), width = 0.1, size = 0.4, color = "grey20") +
  geom_point(data = filtered_read_counts_segments, aes(x = extraction, y = passing_segments, fill = extraction), 
             position = position_dodge(width = 0.9), shape = 21, size = 1, alpha = 0.1) +
  labs(title = "", 
       y = "Passing Segments",
       x = "",
       fill = "") +
  theme_classic() +
  theme1 +
  scale_fill_manual(values = c("#287271", "#e76f51"))+
  scale_y_continuous(limits = c(0,8.8), breaks = c(0,1,2,3,4,5,6,7,8)) +
  geom_text(aes(label = n, y = 8.8), 
            position = position_dodge(width = 0.9), size = 2)


combined_barplots <- ggarrange(total_reads_bar_plot, mapped_reads_bar_plot, passing_segments_bar_plot, ncol = 1, common.legend = TRUE, legend = "bottom", labels = c("a", "b", "c"), font.label = list(size = 10)) +
  theme(plot.margin = unit(c(0,0,0,0), "null"))


### Coverage plots

coverage_extraction <- ggplot() +
  geom_line(data = read_coverage_extraction, aes(x = Position, y = Coverage.Depth, color = extraction, group = interaction(sample, extraction)), linewidth = 0.4) +
  labs(title = "", x = "", y = "Read depth") +
  theme_minimal() +
  theme2 +
  scale_color_manual(values = c("#287271", "#e76f51")) +
  facet_grid(sample~segment, scales = "free_x", switch = "x") +
  scale_y_continuous(breaks = c(1, 100, 10000), trans = "log10")


### Final figure 1 ###

figure1 <- ggarrange(combined_barplots, coverage_extraction, 
                                      ncol = 2, 
                                      labels = c("a", "d"), 
                                      font.label = list(size = 10),
                                      widths = c(0.8, 2)) +
  theme(plot.margin = unit(c(0,0,0,0), "null"))

ggsave("figures/Figure1.png", plot = figure1, dpi = 600, width = 16 , height = 15, units = "cm", bg = "white")


