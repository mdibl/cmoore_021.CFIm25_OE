#Rscript combine_all_sample_PAS_count_tables.R --csv your_input_dir --out merged.csv --pattern ".*\\.csv$"

if (!require('tidyverse')) install.packages('tidyverse'); library('tidyverse')

args <- commandArgs(TRUE)
args_v <- NA
if (length(args) > 1) {
  args_v <- args[c(FALSE, TRUE)]
  names(args_v) <- sub("\\-+", "", args[c(TRUE, FALSE)])
  print(args_v)
}
args <- NULL

csv_dir <- ifelse(is.na(args_v["csv"]), "_csv_", args_v["csv"])
out_dir <- ifelse(is.na(args_v["out"]), "_out_", args_v["out"])
file_pattern <- ifelse(is.na(args_v["pattern"]), "*.all.reads.csv$", args_v["pattern"]) # sy

filenames <- list.files(csv_dir, pattern = file_pattern, full.names = TRUE)
print(filenames)
ldf <- lapply(filenames, read.csv)

cname <- colnames(ldf[[1]])[1] # sy
cluster.all.reads <- Reduce(function(x, y) merge(x, y, by = cname, all = TRUE,sort=F), ldf) # sy
cluster.all.reads[is.na(cluster.all.reads)] <- 0

# cluster.all.reads <- cluster.all.reads %>% separate(PasID, sep = ":", into = c("chromosome", "strand", "position"), remove = T) #sy

write.csv(cluster.all.reads, out_dir, row.names = FALSE, quote = FALSE)
