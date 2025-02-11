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
  group_by(host, sample, extraction, OT, primers, OT_primers) %>%
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

################################################################################

### prepare data for plotting ###

#################################


### get all human samples tested in the three OT-PCR methods
hu_other_OT_samples <- summarized_read_counts %>%
  group_by(sample) %>%
  filter(all(c("OT1-MBTR", "OT2-UNI", "OT3-MBTA") %in% OT_primers)) %>%
  distinct(sample) %>%
  filter(!sample == "HU-22")

### filter tables to represent those samples only
hu_other_OT_total_reads <- summarized_read_counts %>%
  filter(sample %in% hu_other_OT_samples$sample)

### filter table to represent those samples only
hu_other_OT_segment_reads <- segment_read_counts %>%
  filter(sample %in% hu_other_OT_samples$sample)
hu_other_OT_segment_reads$segment <- factor(hu_other_OT_segment_reads$segment, levels = c("PB2", "PB1", "PA", "HA", "NP","NA", "MP", "NS"))


### filter table to represent those samples only
hu_other_OT_mapping_stats <- mapping_stats %>%
  filter(sample %in% hu_other_OT_samples$sample)
hu_other_OT_mapping_stats$segment <- factor(hu_other_OT_mapping_stats$segment, levels = c("PB2", "PB1", "PA", "HA", "NP","NA", "MP", "NS"))


##################################################################


### Plotting ###

theme1 <- theme(legend.position = "top",
                legend.text = element_text(size = 8),
                legend.title = element_blank(),
                legend.key.size = unit(8, "points"),
                axis.text.x = element_text(size = 6, vjust = 0.5),
                axis.ticks.x = element_blank(),
                axis.text.y = element_text(size = 6),
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
                plot.margin = unit(c(0,0.2,0,0.2), "lines"))

theme2 <- theme(legend.position = "top",
                legend.text = element_text(size = 8),
                legend.title = element_blank(),
                legend.key.size = unit(8, "points"),
                axis.text.x = element_text(size = 6, vjust = 0.5, angle=45),
                axis.text.y = element_text(size = 6),
                axis.title.y = element_text(size = 8),
                strip.text.y = element_text(size = 8),
                strip.text.x = element_text(size = 0),
                axis.line = element_line(size = 0.5),
                axis.ticks = element_line(size = 0.5),
                strip.background = element_rect(size = 0.6),
                panel.spacing.x = unit(0.2, "lines"), 
                panel.spacing.y = unit(0.1, "lines"),
                panel.grid.minor.x = element_blank(),
                panel.grid.major.x = element_blank(),
                legend.box.spacing = unit(0, "lines"),
                plot.margin = unit(c(0.1,0.2,0.1,0.2), "lines"))



read_counts_human_plot <- ggplot(hu_other_OT_total_reads, aes(x = OT_primers, y = total_reads, fill = OT_primers )) +
  geom_boxplot(position = position_dodge(width = 1), alpha = 0.9, size = 0.3, color = "grey40", outlier.size = 0.5) +
  labs(title = "", 
       y = "Read Count",
       x = "",
       legend = "") +
  theme_classic() +
  theme1 +
  scale_fill_manual(values = c("#e76f51", "#264653", "#fadfb3")) +
  scale_y_log10() 


mapped_reads_plot <- ggplot(hu_other_OT_total_reads, aes(x = OT_primers, y = percent_matched, fill = OT_primers)) +
  geom_boxplot(position = position_dodge(width = 1), alpha = 0.9, size = 0.3, color = "grey40", outlier.size = 0.5) +
  labs(title = "", 
       y = "% Influenza Reads",
       x = "",
       legend = "") +
  theme_classic() +
  theme1 +
  scale_fill_manual(values = c("#e76f51", "#264653", "#fadfb3"))


segment_read_percentages_plot <- ggplot(hu_other_OT_segment_reads, aes(x = segment, y = percent_of_matched_reads, fill = OT_primers)) +
  geom_boxplot(position = position_dodge(width = 1), alpha = 0.9, size = 0.3, color = "grey40", outlier.size = 0.5) +
  geom_hline(yintercept = 12.5, linetype = 2, size = 0.5, alpha = 0.6) +
  labs(title = "", x = "", y = "% Reads") +
  theme_classic() +
  facet_grid(~OT_primers)+
  scale_fill_manual(values = c("#e76f51", "#264653", "#fadfb3")) +
  theme2 +
  scale_y_continuous(limits = c(0, 100))


missing_sites_plot <- ggplot(hu_other_OT_mapping_stats, aes(x = segment, y = percent_missing_sites, fill = OT_primers)) +
  geom_boxplot(position = position_dodge(width = 1), alpha = 0.9, size = 0.3, color = "grey40", outlier.size = 0.5) +
  labs(title = "", x = "", y = "% Missing sites") +
  theme_classic() +
  facet_grid(~OT_primers) +
  scale_fill_manual(values = c("#e76f51", "#264653", "#fadfb3")) +
  theme2


combined_plot <- ggarrange(read_counts_human_plot, mapped_reads_plot, segment_read_percentages_plot, missing_sites_plot,
                                             ncol = 2, nrow = 2,
                                             common.legend = TRUE,
                                             labels = c("a", "b", "c", "d"), 
                                             font.label = list(size= 10)) +
  theme(plot.margin = unit(c(0.05,0,0,0), "null"))


ggsave("figures/FigureS2.png", plot = combined_plot, dpi = 300, width = 18, height = 20, units = "cm" , bg = "white")
