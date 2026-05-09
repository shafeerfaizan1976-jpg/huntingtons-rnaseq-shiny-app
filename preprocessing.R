library(tidyverse)

# paths to raw GEO files
counts_path <- "data/GSE64810_mlhd_DESeq2_norm_counts_adjust.txt.gz"
de_path     <- "data/GSE64810_mlhd_DESeq2_diffexp_DESeq2_outlier_trimmed_adjust.txt.gz"
meta_path   <- "data/GSE64810_series_matrix.txt.gz"

# read normalized counts and rename first column
counts_raw <- read_tsv(counts_path, show_col_types = FALSE)
counts_clean <- counts_raw %>% rename(gene = `...1`)

# read DE results and rename first column
de_raw <- read_tsv(de_path, show_col_types = FALSE)
de_clean <- de_raw %>% rename(gene = `...1`)

# parse sample metadata from the series matrix file
meta_lines <- read_lines(meta_path)

# split a tab-separated line, drop the leading tag, strip surrounding quotes
parse_row <- function(line) {
  parts <- strsplit(line, "\t")[[1]]
  parts <- parts[-1]
  gsub('^"|"$', '', parts)
}

sample_ids <- parse_row(meta_lines[grep("^!Sample_title",         meta_lines)])
geo_ids    <- parse_row(meta_lines[grep("^!Sample_geo_accession", meta_lines)])

# each characteristics line has a "key: value" format per sample
char_lines <- meta_lines[grep("^!Sample_characteristics_ch1", meta_lines)]

extract_field <- function(line) {
  values   <- parse_row(line)
  nonempty <- values[values != ""]
  key      <- sub(":\\s*.*$", "", nonempty[1])
  vals     <- sub("^[^:]+:\\s*", "", values)
  vals[values == ""] <- NA
  list(key = key, values = vals)
}

parsed <- lapply(char_lines, extract_field)

sample_info <- tibble(sample_id = sample_ids, geo_accession = geo_ids)

for (p in parsed) {
  col_name <- tolower(p$key)
  col_name <- gsub("[^a-z0-9]+", "_", col_name)
  col_name <- gsub("^_|_$", "", col_name)
  sample_info[[col_name]] <- p$values
}

# convert numeric fields (HD-only clinical fields will have NAs for controls)
numeric_cols <- c("pmi", "age_of_death", "rin", "mrna_seq_reads",
                  "age_of_onset", "duration", "cag",
                  "vonsattel_grade", "h_v_striatal_score", "h_v_cortical_score")
numeric_cols <- intersect(numeric_cols, colnames(sample_info))

sample_info <- sample_info %>%
  mutate(across(all_of(numeric_cols), as.numeric))

# sanity checks
stopifnot(nrow(sample_info) == 69)
stopifnot(all(sample_info$sample_id %in% colnames(counts_clean)))

# save cleaned files
write_csv(sample_info, "data/sample_info.csv")
write_csv(counts_clean, "data/norm_counts.csv")
write_csv(de_clean, "data/de_results.csv")