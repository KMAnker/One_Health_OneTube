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
  OT_primers <- paste(OT, primers, sep="-")
  segment2 = str_extract(record, "(?<=_)[^_]+$")
  segment = case_when(
    str_detect(segment2, "H\\d+") ~ "HA",
    str_detect(segment2, "N\\d+") ~ "NA",
    TRUE ~ segment2)
  return(data.frame(host = host, sample = sample, extraction = extraction, OT = OT, 
                    primers = primers, OT_primers = OT_primers, segment = segment, 
                    segment2 = segment2, record = record))
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


# summarize the table to get total read no. per sample, no. of reads passing filters and no. of matched reads
summarized_read_counts <- read_counts %>%
  filter(record %in% c("1-initial", "3-match")) %>%
  mutate(record = case_when(
    record == "1-initial" ~ "Total",
    record == "3-match" ~ "Matched")) %>%
  group_by(host, sample, extraction, OT, primers) %>%
  summarize(
    total_reads = sum(ifelse(record == "Total", as.numeric(Reads), 0)),
    reads_matched = sum(ifelse(record == "Matched", as.numeric(Reads), 0)),
    percent_matched = (reads_matched / total_reads) * 100
  ) %>%
  ungroup()

# summarize the table to get segment specific reads and total matched reads, and calculate percentage of matched reads to each segment for each sample
segment_read_counts <- read_counts %>%
  filter(record == "3-match" | str_detect(record, "^4")) %>%
  select(-c(Patterns, PairsAndWidows)) %>%
  mutate(Reads = as.numeric(Reads),
         matched_reads = ifelse(record == "3-match", as.numeric(Reads), NA))%>%
  group_by(sample, extraction, OT, primers, OT_primers) %>%
  fill(matched_reads, .direction = "down") %>%
  fill(matched_reads, .direction = "up") %>%
  filter(record != "3-match") %>%
  mutate(percent_of_matched_reads = ifelse(!is.na(matched_reads), (Reads / matched_reads) * 100, NA)) %>%
  ungroup()

# Order categories
segment_read_counts$OT_primers <- factor(segment_read_counts$OT_primers, levels = c("OT1-MBT", "OT1-UNI", "OT2-MBT", "OT2-UNI", "OT1-MBT2"))
segment_read_counts$segment <- factor(segment_read_counts$segment, levels = c("PB2", "PB1", "PA", "HA", "NP","NA", "MP", "NS"))
segment_read_counts$host <- factor(segment_read_counts$host, levels = c("HU", "SW", "AV"))

################################################################################

### Read in passing segments data ###
passing_segments <- read.table("data/passing_segments_all.txt", header=TRUE)

filtered_passing_segments <- passing_segments %>%
  mutate(
    parts = str_split_fixed(Method, "-", 3),
    extraction = na_if(parts[, 1], ""),
    OT = na_if(parts[, 2], ""),
    primers = na_if(parts[, 3], ""),
    OT_primers = paste(OT,primers, sep = "-")) %>%
  select(-parts, -Method) %>%
  rename(sample = Sample)

################################################################################

### Get missing sites percentages for each segment from mapping_stats file ###

mapping_stats_table <- read.table("data/mapping_stats_all.txt", header=TRUE)
colnames(mapping_stats_table)[8] ="percent_missing_sites"

# Function to extract sample, method, OT, primers and segment from the Sample column
extract_sample_method_2 <- function(Sample) {
  sample_method <- str_extract(Sample, "^[^_]+")
  parts_sample_method <- unlist(strsplit(sample_method, "-", fixed = TRUE))
  # Assign values or NA for missing parts
  host <- paste(parts_sample_method[1])
  sample <- paste(parts_sample_method[1], parts_sample_method[2], sep = "-")
  extraction <- ifelse(length(parts_sample_method) >= 3, parts_sample_method[3], NA)
  OT <- ifelse(length(parts_sample_method) >= 4, parts_sample_method[4], NA)
  primers <- ifelse(length(parts_sample_method) >= 5, parts_sample_method[5], NA)
  OT_primers <- paste(OT, primers, sep="-")
  segment2 = str_extract(Sample, "(?<=_)[^_]+$")
  segment = case_when(
    str_detect(segment2, "H\\d+") ~ "HA",
    str_detect(segment2, "N\\d+") ~ "NA",
    TRUE ~ segment2)
  return(data.frame(host = host, sample = sample, extraction = extraction, OT = OT, 
                    primers = primers, OT_primers = OT_primers, segment = segment, 
                    segment2 = segment2))
}

# Apply the function to create the new columns
mapping_stats <- mapping_stats_table %>%
  rowwise() %>%
  do({
    extracted_cols <- extract_sample_method_2(.$Sample)
    data.frame(extracted_cols, .)
  }) %>%
  ungroup() %>%
  select(-Sample)

# Order categories
mapping_stats$OT_primers <- factor(mapping_stats$OT_primers, levels = c("OT1-MBT", "OT1-UNI", "OT2-MBT", "OT2-UNI", "OT1-MBT2"))
mapping_stats$segment <- factor(mapping_stats$segment, levels = c("PB2", "PB1", "PA", "HA", "NP","NA", "MP", "NS"))
mapping_stats$host <- factor(mapping_stats$host, levels = c("HU", "SW", "AV"))

################################################################################

### Read in ct data ###
sample_ct <- read.table("data/sample_cts.txt", header = TRUE)

################################################################################

# prepare data for plotting

#################################

# Read counts #

### combine summarized read counts table with passing segments and ct data ###
summarized_data <- summarized_read_counts %>%
  left_join(filtered_passing_segments, by = c("sample", "extraction", "OT", "primers")) %>%
  left_join(sample_ct, by = c("host", "sample", "extraction"))
summarized_data$host <- factor(summarized_data$host, levels = c("HU", "SW", "AV"))

### get all samples tested with extra primer set
extra_primers_samples <- summarized_data %>%
  group_by(sample) %>%
  filter(any(primers == "MBT2")) %>%
  distinct(sample)

### filter table to represent those samples only
extra_primers_data <- summarized_data %>%
  filter(sample %in% extra_primers_samples$sample)

### filter away samples with ct < 36
extra_primers_filtered <- extra_primers_data %>%
  group_by(sample) %>%
  filter(any(ct < 36)) %>%
  ungroup()

# order categories
extra_primers_filtered$OT_primers <- 
  factor(extra_primers_filtered$OT_primers, 
         levels = c("OT1-MBT", "OT1-UNI", "OT2-MBT", "OT2-UNI", "OT1-MBT2"))


### Make summaries of means and standard error

summary_total_extra_primers <- extra_primers_filtered %>%
  group_by(host, extraction, OT_primers) %>%
  summarize(mean_total_reads = mean(total_reads),
            sd_total_reads = sd(total_reads),
            n = n(),
            se_total_reads = sd(total_reads) / sqrt(n()), .groups = 'drop')

summary_matched_extra_primers <- extra_primers_filtered %>%
  group_by(host, extraction, OT_primers) %>%
  summarize(mean_percent_matched = mean(percent_matched),
            sd_percent_matched = sd(percent_matched),
            n = n(),
            se_percent_matched = sd(percent_matched) / sqrt(n()), .groups = 'drop')

summary_passing_segments_extra_primers <- extra_primers_filtered %>%
  group_by(host, extraction, OT_primers) %>%
  summarize(mean_passing_segments = mean(passing_segments),
            sd_passing_segments = sd(passing_segments),
            n = n(),
            se_passing_segments = sd(passing_segments) / sqrt(n()), .groups = 'drop')

#################################

# Segment read percentages #

### combine segment read counts table with ct data ###
summarized_segment_reads <- segment_read_counts %>%
  left_join(sample_ct, by = c("host", "sample", "extraction"))
summarized_segment_reads$host <- factor(summarized_segment_reads$host, levels = c("HU", "SW", "AV"))

### get all samples tested with extra primer set
extra_primers_samples <- summarized_segment_reads %>%
  group_by(sample) %>%
  filter(any(primers == "MBT2")) %>%
  distinct(sample)

### filter table to represent those samples only
extra_primers_segment_reads <- summarized_segment_reads %>%
  filter(sample %in% extra_primers_samples$sample)

### filter away samples with ct < 36
extra_primers_segment_reads_filtered <- extra_primers_segment_reads %>%
  group_by(sample) %>%
  filter(any(ct < 36)) %>%
  ungroup()

# order categories
extra_primers_segment_reads_filtered$OT_primers <- 
  factor(extra_primers_segment_reads_filtered$OT_primers, 
         levels = c("OT1-MBT", "OT1-UNI", "OT2-MBT", "OT2-UNI", "OT1-MBT2"))


### summarize data to get mean and standard error of segment reads
extra_primers_segment_reads_summarized <- extra_primers_segment_reads_filtered %>%
  group_by(host, extraction, OT_primers, segment) %>%
  summarize(mean_segment_reads = mean(percent_of_matched_reads),
            n = n(),
            se_segment_reads = sd(percent_of_matched_reads) / sqrt(n()), .groups = 'drop')

n_segment_reads_extra <- extra_primers_segment_reads_summarized %>%
  group_by(host, extraction) %>%
  summarize(min_n = min(n), max_n = max(n), .groups = 'drop')

#################################

# Percent missing sites #

### combine mapping stats table with ct data ###
summarized_mapping_stats <- mapping_stats %>%
  left_join(sample_ct, by = c("host", "sample", "extraction"))
summarized_mapping_stats$host <- factor(summarized_mapping_stats$host, levels = c("HU", "SW", "AV"))

### get all samples tested with extra primer set
extra_primers_samples <- summarized_mapping_stats %>%
  group_by(sample) %>%
  filter(any(primers == "MBT2")) %>%
  distinct(sample)

### filter table to represent those samples only
extra_primers_mapping_stats <- summarized_mapping_stats %>%
  filter(sample %in% extra_primers_samples$sample)

### filter away samples with ct < 36
extra_primers_mapping_stats_filtered <- extra_primers_mapping_stats %>%
  group_by(sample) %>%
  filter(any(ct < 36)) %>%
  ungroup()

# order categories
extra_primers_mapping_stats_filtered$OT_primers <- 
  factor(extra_primers_mapping_stats_filtered$OT_primers, 
         levels = c("OT1-MBT", "OT1-UNI", "OT2-MBT", "OT2-UNI", "OT1-MBT2"))

### seperate by host 
extra_primers_mapping_stats_filtered_hu <- extra_primers_mapping_stats_filtered %>%
  filter(host == "HU")
extra_primers_mapping_stats_filtered_sw <- extra_primers_mapping_stats_filtered %>%
  filter(host == "SW")
extra_primers_mapping_stats_filtered_av <- extra_primers_mapping_stats_filtered %>%
  filter(host == "AV")

################################################################################

### Plotting ###

theme1 <- theme(legend.position = "top",
                legend.text = element_text(size = 7),
                legend.title = element_blank(),
                legend.key.size = unit(7, "points"),
                axis.text.x = element_text(size = 7, vjust = 0.5),
                axis.ticks.x=element_blank(),
                axis.text.y = element_text(size = 7),
                axis.title.y = element_text(size = 7),
                strip.text = element_text(size = 7),
                axis.line = element_line(size = 0.5),
                axis.ticks = element_line(size = 0.5),
                strip.background = element_rect(size = 0.6),
                panel.spacing.x = unit(0.2, "lines"), 
                panel.spacing.y = unit(0.1, "lines"),
                panel.grid.minor.x = element_blank(),
                panel.grid.major.x = element_blank(),
                legend.box.spacing = unit(0, "lines"),
                plot.margin = unit(c(0,0.2,0,0.2), "lines"))

theme2 <- theme(legend.position = "top",
                legend.text = element_text(size = 7),
                legend.title = element_blank(),
                legend.key.size = unit(7, "points"),
                axis.text.x = element_text(size = 7, angle = 45, vjust = 0.7),
                axis.text.y = element_text(size = 7),
                axis.title.y = element_text(size = 7),
                strip.text.y = element_text(size = 7),
                strip.text.x = element_text(size = 7),
                axis.line = element_line(size = 0.5),
                axis.ticks = element_line(size = 0.5),
                strip.background = element_rect(size = 0.6),
                panel.spacing.x = unit(0.2, "lines"), 
                panel.spacing.y = unit(0.2, "lines"),
                panel.grid.minor.x = element_blank(),
                panel.grid.major.x = element_blank(),
                legend.box.spacing = unit(0, "lines"),
                plot.margin = unit(c(2,0.4,2,0.4), "lines"))

theme3 <- theme(legend.position = "top", 
                legend.text = element_text(size = 7),
                legend.title = element_blank(),
                legend.key.size = unit(7, "points"),
                plot.title = element_text(hjust = 0.5, size = 7),
                axis.text.x = element_text(size = 7, angle = 45, vjust = 0.7),
                axis.text.y = element_text(size = 7),
                axis.title.y = element_text(size = 7),
                strip.text.y = element_text(size = 0),
                strip.text.x = element_text(size = 7),
                axis.line = element_line(size = 0.5),
                axis.ticks = element_line(size = 0.5),
                strip.background = element_rect(size = 0.6),
                strip.background.y = element_blank(),
                panel.spacing.x = unit(0.2, "lines"),
                panel.spacing.y = unit(0.1, "lines"),
                panel.grid.minor.x = element_blank(),
                panel.grid.major.x = element_blank(),
                legend.box.spacing = unit(0, "lines"),
                plot.margin = unit(c(0,0,0,1), "lines"))

theme4 <- theme(legend.position = "top", 
                legend.text = element_text(size = 7),
                legend.title = element_blank(),
                legend.key.size = unit(7, "points"),
                plot.title = element_text(hjust = 0.5, size = 7),
                axis.text.x = element_text(size = 7, angle = 45, vjust = 0.7),
                axis.text.y = element_blank(),
                axis.title.y = element_blank(),
                strip.text.y = element_text(size = 0),
                strip.text.x = element_text(size = 7),
                axis.line = element_line(size = 0.5),
                axis.ticks = element_line(size = 0.5),
                strip.background = element_rect(size = 0.6),
                strip.background.y = element_blank(),
                panel.spacing.x = unit(0.2, "lines"),
                panel.spacing.y = unit(0.1, "lines"),
                panel.grid.minor.x = element_blank(),
                panel.grid.major.x = element_blank(),
                legend.box.spacing = unit(0, "lines"),
                plot.margin = unit(c(0,0,0,0), "lines"))

################################################################################

total_reads_extra_primers_bar_plot <- ggplot(summary_total_extra_primers, aes(x = extraction, y = mean_total_reads, fill = OT_primers)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.9), alpha = 0.9, color = "grey40", size = 0.4) +
  geom_errorbar(aes(ymin = mean_total_reads - se_total_reads, ymax = mean_total_reads + se_total_reads), 
                position = position_dodge(width = 0.9), width = 0.25, size = 0.4, color = "grey20") +
  geom_point(data = extra_primers_filtered, aes(x = extraction, y = total_reads, fill = OT_primers), 
             position = position_dodge(width = 0.9), shape = 21, size = 1, alpha = 0.1) +
  labs(title = "", 
       y = "Mean Read Count",
       x = "",
       fill = "") +
  theme_classic() +
  theme1 +
  scale_fill_manual(values = c("#e76f51", "#f4a261", "#2a9d8f", "#264653", "#613653")) +
  facet_wrap(~host, ncol = 3) +
  scale_y_log10() +
  geom_text(aes(label = n, y = max(extra_primers_filtered$total_reads)+500000), 
            position = position_dodge(width = 0.9), size = 2)


mapped_reads_extra_primers_bar_plot <- ggplot(summary_matched_extra_primers, aes(x = extraction, y = mean_percent_matched, fill = OT_primers)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.9), alpha = 0.9, color = "grey40", size = 0.4) +
  geom_errorbar(aes(ymin = mean_percent_matched - se_percent_matched, ymax = mean_percent_matched + se_percent_matched), 
                position = position_dodge(width = 0.9), width = 0.25, size = 0.4, color = "grey20") +
  geom_point(data = extra_primers_filtered, aes(x = extraction, y = percent_matched, fill = OT_primers), 
             position = position_dodge(width = 0.9), shape = 21, size = 1, alpha = 0.1) +
  labs(title = "", 
       y = "% Influenza Reads",
       x = "",
       fill = "") +
  theme_classic() +
  theme1 +
  scale_fill_manual(values = c("#e76f51", "#f4a261", "#2a9d8f","#264653", "#613653")) +
  facet_wrap(~host, ncol = 3) +
  geom_text(aes(label = n, y = 105), 
            position = position_dodge(width = 0.9), size = 2)

passing_segments_extra_primers_bar_plot <- ggplot(summary_passing_segments_extra_primers, aes(x = extraction, y = mean_passing_segments, fill = OT_primers)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.9), alpha = 0.9, color = "grey40", size = 0.4) +
  geom_errorbar(aes(ymin = mean_passing_segments - se_passing_segments, ymax = mean_passing_segments + se_passing_segments), 
                position = position_dodge(width = 0.9), width = 0.25, size = 0.4, color = "grey20") +
  geom_point(data = extra_primers_filtered, aes(x = extraction, y = passing_segments, fill = OT_primers), 
             position = position_dodge(width = 0.9), shape = 21, size = 1, alpha = 0.1) +
  labs(title = "", 
       y = "Passing Segments",
       x = "",
       fill = "") +
  theme_classic() +
  theme1 +
  scale_fill_manual(values = c("#e76f51", "#f4a261", "#2a9d8f", "#264653", "#613653")) +
  scale_y_continuous(limits = c(0,8.5), breaks = c(0,1,2,3,4,5,6,7,8)) +
  facet_wrap(~host, ncol = 3) +
  geom_text(aes(label = n, y = 8.5), 
            position = position_dodge(width = 0.9), size = 2)


read_stats_total <- ggarrange(total_reads_extra_primers_bar_plot, mapped_reads_extra_primers_bar_plot, passing_segments_extra_primers_bar_plot, 
                              ncol = 1, 
                              common.legend = TRUE, 
                              legend = "top")



### FigS3d: Read percentage per segment ###
segment_read_percentages_plot_extra <- ggplot(extra_primers_segment_reads_summarized, aes(x = segment, y = mean_segment_reads, color = OT_primers)) +
  geom_point(size = 1) +
  geom_line(aes(group = OT_primers), size = 0.7) +
  geom_errorbar(aes(ymin = mean_segment_reads - se_segment_reads, ymax = mean_segment_reads + se_segment_reads), 
                , width = 0.2, size = 0.2) +
  geom_hline(yintercept = 12.5, linetype = 2, size = 0.3, alpha = 0.6) +
  labs(title = "", x = "", y = "% Reads") +
  theme_classic() +
  facet_grid(extraction~host, labeller = labeller(host = c(HU = "Human Samples", SW = "Swine Samples", AV = "Avian Samples"))) +
  scale_color_manual(values = c("#e76f51", "#f4a261", "#2a9d8f", "#264653", "#613653")) +
  scale_y_continuous(limits = c(0,45), breaks = seq(0, 40, by = 10)) +
  theme2 +
  geom_text(data = n_segment_reads_extra, aes(x = 4.5, y = 45, 
                                              label = ifelse(min_n == max_n, paste0("n = ", min_n), 
                                                             paste0("n = ", min_n, "-", max_n))), 
            hjust = 0.5, vjust = 2, size = 2, inherit.aes = FALSE)


### FigS3e: Missing sites per segment ###

percent_missing_sites_extra_hu <- ggplot(extra_primers_mapping_stats_filtered_hu, aes(x = segment, y = percent_missing_sites, fill = OT_primers)) +
  geom_boxplot(position = position_dodge(width = 1), alpha = 1, size = 0.3, outlier.size = 0.3, color = "grey40") +
  labs(title = "Human samples", x = "", y = "% Missing sites") +
  theme_classic() +
  facet_grid(OT_primers~extraction, labeller = labeller(host = c(HU = "Human Samples", SW = "Swine Samples", AV = "Avian Samples"))) +
  scale_fill_manual(values = c("#e76f51", "#f4a261", "#2a9d8f", "#264653", "#613653")) +
  geom_vline(xintercept = c(0.4,8.6), size = 0.4) +
  theme3

percent_missing_sites_extra_sw <- ggplot(extra_primers_mapping_stats_filtered_sw, aes(x = segment, y = percent_missing_sites, fill = OT_primers)) +
  geom_boxplot(position = position_dodge(width = 1), alpha = 1, size = 0.3, outlier.size = 0.3, color = "grey40") +
  labs(title = "Swine samples", x = "", y = "") +
  theme_classic() +
  facet_grid(OT_primers~extraction, labeller = labeller(host = c(HU = "Human Samples", SW = "Swine Samples", AV = "Avian Samples"))) +
  scale_fill_manual(values = c("#e76f51", "#f4a261", "#2a9d8f", "#264653", "#613653")) +
  geom_vline(xintercept = c(0.4,8.6), size = 0.4) +
  theme4

percent_missing_sites_extra_av <- ggplot(extra_primers_mapping_stats_filtered_av, aes(x = segment, y = percent_missing_sites, fill = OT_primers)) +
  geom_boxplot(position = position_dodge(width = 1), alpha = 1, size = 0.3, outlier.size = 0.3, color = "grey40") +
  labs(title = "Avian samples", x = "", y = "") +
  theme_classic() +
  facet_grid(OT_primers~extraction, labeller = labeller(host = c(HU = "Human Samples", SW = "Swine Samples", AV = "Avian Samples"))) +
  scale_fill_manual(values = c("#e76f51", "#f4a261", "#2a9d8f", "#264653", "#613653")) +
  geom_vline(xintercept = c(0.4,8.6), size = 0.4) +
  theme4

missing_sites_plot_extra <- ggarrange(percent_missing_sites_extra_hu, percent_missing_sites_extra_sw, percent_missing_sites_extra_av,
                                      ncol = 3, 
                                      common.legend = TRUE,
                                      widths = c(1, 0.85, 0.85))

read_stats_all <- ggarrange(read_stats_total, segment_read_percentages_plot_extra, 
                            ncol = 2, 
                            labels = c("a", "b"), 
                            font.label = list(size = 10))

FigS3 <- ggarrange(read_stats_all, missing_sites_plot_extra,
                   ncol = 1, 
                   labels = c("a", "c"), 
                   font.label = list(size = 10),
                   heights = c(1,0.8))

ggsave("figures/FigureS3.png", plot = FigS3, dpi = 600, width = 20 , height = 25, units = "cm", bg = "white")

