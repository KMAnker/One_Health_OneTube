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

### Read in ct data ###
sample_ct <- read.table("data/sample_cts.txt", header = TRUE)

################################################################################

# prepare data for plotting

### combine summarized read counts table with passing segments and ct data ###
summarized_data <- summarized_read_counts %>%
  left_join(filtered_passing_segments, by = c("sample", "extraction", "OT", "primers")) %>%
  left_join(sample_ct, by = c("host", "sample", "extraction"))
summarized_data$host <- factor(summarized_data$host, levels = c("HU", "SW", "AV"))


### get all samples tested with all 4 PCR methods (not necessarily both EX)
all_pcr_samples <- summarized_data %>%
  filter(OT_primers %in% c("OT1-MBT", "OT1-UNI", "OT2-MBT", "OT2-UNI")) %>%
  distinct(sample)


### filter table to represent those samples only (and leave out extra PCR method for now)
all_pcr_data <- summarized_data %>%
  filter(sample %in% all_pcr_samples$sample) %>%
  filter(!primers == "MBT2")

### filter away samples with ct > 36 for all extractions and human samples from another test of pcr methods 
all_pcr_filtered <- all_pcr_data %>%
  filter(!sample %in% c("HU-22", "HU-23", "HU-24", "HU-25", "HU-26", "HU-27", "HU-28", "HU-29", "HU-30")) %>%
  group_by(sample) %>%
  filter(any(ct < 36)) %>%
  ungroup()

# order categories
all_pcr_filtered$OT_primers <- 
  factor(all_pcr_filtered$OT_primers, 
         levels = c("OT1-MBT", "OT1-UNI", "OT2-MBT", "OT2-UNI"))


### Make summaries of means and standard error

summary_total <- all_pcr_filtered %>%
  group_by(host, extraction, OT_primers) %>%
  summarize(mean_total_reads = mean(total_reads),
            sd_total_reads = sd(total_reads),
            n = n(),
            se_total_reads = sd(total_reads) / sqrt(n()), .groups = 'drop')

summary_matched <- all_pcr_filtered %>%
  group_by(host, extraction, OT_primers) %>%
  summarize(mean_percent_matched = mean(percent_matched),
            sd_percent_matched = sd(percent_matched),
            n = n(),
            se_percent_matched = sd(percent_matched) / sqrt(n()), .groups = 'drop')

summary_passing_segments <- all_pcr_filtered %>%
  group_by(host, extraction, OT_primers) %>%
  summarize(mean_passing_segments = mean(passing_segments),
            sd_passing_segments = sd(passing_segments),
            n = n(),
            se_passing_segments = sd(passing_segments) / sqrt(n()), .groups = 'drop')

################################################################################

### Plotting ###

theme1 <- theme(legend.position = "top",
                legend.text = element_text(size = 8),
                legend.title = element_blank(),
                legend.key.size = unit(8, "points"),
                axis.text.x = element_text(size = 8, vjust = 0.5),
                axis.ticks.x=element_blank(),
                axis.text.y = element_text(size = 8),
                axis.title.y = element_text(size = 8),
                strip.text = element_text(size = 8),
                axis.line = element_line(size = 0.5),
                axis.ticks = element_line(size = 0.5),
                strip.background = element_rect(size = 0.6),
                panel.spacing.x = unit(0.2, "lines"), 
                panel.spacing.y = unit(0.1, "lines"),
                panel.grid.minor.x = element_blank(),
                panel.grid.major.x = element_blank(),
                legend.box.spacing = unit(0, "lines"),
                plot.margin = unit(c(0,0.4,0,0.4), "lines"))



total_reads_bar_plot <- ggplot(summary_total, aes(x = extraction, y = mean_total_reads, fill = OT_primers)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.9), alpha = 0.9, color = "grey40", size = 0.4) +
  geom_errorbar(aes(ymin = mean_total_reads - se_total_reads, ymax = mean_total_reads + se_total_reads), 
                position = position_dodge(width = 0.9), width = 0.25, size = 0.4, color = "grey20") +
  geom_point(data = all_pcr_filtered, aes(x = extraction, y = total_reads, fill = OT_primers), 
             position = position_dodge(width = 0.9), shape = 21, size = 1, alpha = 0.1) +
  labs(title = "", 
       y = "Mean Read Count",
       x = "",
       fill = "") +
  theme_classic() +
  theme1 +
  scale_fill_manual(values = c("#e76f51", "#f4a261", "#2a9d8f", "#264653")) +
  facet_wrap(~host, ncol = 3) +
  scale_y_log10() +
  geom_text(aes(label = n, y = max(all_pcr_filtered$total_reads)+500000), 
            position = position_dodge(width = 0.9), size = 1.5)


mapped_reads_bar_plot <- ggplot(summary_matched, aes(x = extraction, y = mean_percent_matched, fill = OT_primers)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.9), alpha = 0.9, color = "grey40", size = 0.4) +
  geom_errorbar(aes(ymin = mean_percent_matched - se_percent_matched, ymax = mean_percent_matched + se_percent_matched), 
                position = position_dodge(width = 0.9), width = 0.25, size = 0.4, color = "grey20") +
  geom_point(data = all_pcr_filtered, aes(x = extraction, y = percent_matched, fill = OT_primers), 
             position = position_dodge(width = 0.9), shape = 21, size = 1, alpha = 0.1) +
  labs(title = "", 
       y = "% Influenza Reads",
       x = "",
       fill = "") +
  theme_classic() +
  theme1 +
  theme(axis.title.y = element_text(margin = unit(c(0,0.8,0,0), "lines"))) +
  scale_fill_manual(values = c("#e76f51", "#f4a261", "#2a9d8f","#264653")) +
  facet_wrap(~host, ncol = 3) +
  geom_text(aes(label = n, y = 105), 
            position = position_dodge(width = 0.9), size = 1.5)

passing_segments_bar_plot <- ggplot(summary_passing_segments, aes(x = extraction, y = mean_passing_segments, fill = OT_primers)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.9), alpha = 0.9, color = "grey40", size = 0.4) +
  geom_errorbar(aes(ymin = mean_passing_segments - se_passing_segments, ymax = mean_passing_segments + se_passing_segments), 
                position = position_dodge(width = 0.9), width = 0.25, size = 0.4, color = "grey20") +
  geom_point(data = all_pcr_filtered, aes(x = extraction, y = passing_segments, fill = OT_primers), 
             position = position_dodge(width = 0.9), shape = 21, size = 1, alpha = 0.1) +
  labs(title = "", 
       y = "Passing Segments",
       x = "",
       fill = "") +
  theme_classic() +
  theme1 +
  theme(axis.title.y = element_text(margin = unit(c(0,1.5,0,0), "lines"))) +
  scale_fill_manual(values = c("#e76f51", "#f4a261", "#2a9d8f", "#264653")) +
  scale_y_continuous(limits = c(0,8.5), breaks = c(0,1,2,3,4,5,6,7,8)) +
  facet_wrap(~host, ncol = 3) +
  geom_text(aes(label = n, y = 8.5), 
            position = position_dodge(width = 0.9), size = 1.5)


figure2 <- ggarrange(total_reads_bar_plot, mapped_reads_bar_plot, passing_segments_bar_plot, 
                     ncol = 1,
                     common.legend = TRUE, 
                     legend = "top", 
                     labels = c("a", "b", "c"), 
                     font.label = list(size = 10)) +
  theme(plot.margin = unit(c(0,0,0,0), "null"))

ggsave("figures/Figure2.png", plot = figure2, dpi = 600, width = 12 , height = 15, units = "cm", bg = "white")

