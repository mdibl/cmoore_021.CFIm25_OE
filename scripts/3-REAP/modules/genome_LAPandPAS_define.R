if (!require('tidyverse')) install.packages('tidyverse'); suppressMessages(library('tidyverse'))
if (!require('plyranges')) install.packages('plyranges'); suppressMessages(library('plyranges'))
if (!require('GenomicRanges')) install.packages('GenomicRanges'); suppressMessages(library('GenomicRanges'))
if (!require('GenomicAlignments')) install.packages('GenomicAlignments'); suppressMessages(library('GenomicAlignments'))
if (!require('dplyr')) install.packages('dplyr'); suppressMessages(library('dplyr'))
if (!require('ggplot2')) install.packages('ggplot2'); suppressMessages(library('ggplot2'))
if (!require('patchwork')) install.packages('patchwork'); suppressMessages(library('patchwork'))

# no . in sample prefix like K562.ctrl

args<- commandArgs(TRUE)
args_v <- setNames(rep(NA, length(args) / 2), character(0))

if (length(args) > 1) {
    keys <- args[c(TRUE, FALSE)]
    values <- args[c(FALSE, TRUE)]
    names(values) <- sub("^\\-+", "", keys)
    args_v <- values
}

bed_dir <- ifelse(!is.na(args_v["bedLAP"]), args_v["bedLAP"], "NULL")
out_dir <- ifelse(!is.na(args_v["out"]), args_v["out"], "NULL")
refPAS_file <- ifelse(!is.na(args_v["refPAS"]), args_v["refPAS"], "NULL")
dist <- ifelse(!is.na(args_v["dist"]), args_v["dist"], 24)
misM <- ifelse(!is.na(args_v["misM"]), args_v["misM"], 2)
dist <- as.numeric(dist)
misM <- as.numeric(misM)

stats = numeric()

#sample <- gsub(".sorted.bed","",basename(bed_dir)) # HepG2_Ctrl_R2.dedup
sample <- unlist(strsplit(basename(bed_dir),split="[.]"))[1]
print(paste0("sample name, ", sample))

# read in bed files
bed <- read.table(bed_dir, header = F, sep = "\t", quote = "", stringsAsFactors = FALSE)
names(bed) <- c("mapped_chr", "mapped_start", "mapped_end", "readID", "MAPQ", "mapped_strand", "CIGAR")
stats[1] <- nrow(bed)
## because the reverse reads are used, the strand should be reversed:
bed$strand_rev <- ifelse(bed$mapped_strand == "-", "+", "-")
bed <- bed %>% separate(readID, sep = "_", into = c("readID", "UMI"), remove = F)


df2 <- bed

# define LAP position
## for PAS in + strand, the LAP should be the end of the mapped region. end-1 due to right open of bed files[);
## for PAS in - strand, the LAP should be the start of the mapped region.
df2$LAP <- ifelse(df2$strand_rev == "+",df2$mapped_end-1, df2$mapped_start) # 0-based

# annotate if the LAP is for real by soft clip
## extract the first and last alignment information
df2$CIGAR_head_symbol <- sapply(explodeCigarOps(df2$CIGAR), head, 1)
df2$CIGAR_head_value <- sapply(explodeCigarOpLengths(df2$CIGAR), head, 1)
df2$CIGAR_tail_symbol <- sapply(explodeCigarOps(df2$CIGAR), tail, 1)
df2$CIGAR_tail_value <- sapply(explodeCigarOpLengths(df2$CIGAR), tail, 1)

## for PAS in + strand, the LAP should be not with end soft clip, with misM nt mismatch allowed;
## for PAS in - strand, the LAP should be not with head soft clip, with misM nt mismatch allowed
df2$LAP_anno <- "LAP"
df2$LAP_anno[df2$strand_rev == "+" & df2$CIGAR_tail_symbol=="S" & df2$CIGAR_tail_value > misM] <- "NO"
df2$LAP_anno[df2$strand_rev == "-" & df2$CIGAR_head_symbol=="S" & df2$CIGAR_head_value > misM] <- "NO"

# LAP_sClip: LAP + [soft clip length]
#df2$LAP_sClip <- df2$LAP
#plus_mapped <- which(df2$strand_rev == "+" & df2$CIGAR_tail_symbol=="S")
#minus_mapped <- which(df2$strand_rev == "-" & df2$CIGAR_head_symbol=="S")
#df2$LAP_sClip[plus_mapped] <- df2$LAP_sClip[plus_mapped]+df2$CIGAR_tail_value[plus_mapped]
#df2$LAP_sClip[minus_mapped] <- df2$LAP_sClip[minus_mapped]+df2$CIGAR_head_value[minus_mapped]

saveRDS(df2, paste0(out_dir,"_LAP_bed.rds"))

## plot the head and tail soft clip length
df2$CIGAR_head_value <- as.numeric(as.character(df2$CIGAR_head_value))
df2$CIGAR_tail_value <- as.numeric(as.character(df2$CIGAR_tail_value))
print(class(df2$CIGAR_head_value));print(class(df2$CIGAR_tail_value))
reso=0.2;wi=.1;st=9
input1=subset(df2,CIGAR_head_symbol=="S")
input2=subset(df2,CIGAR_tail_symbol=="S")
print(class(input1$CIGAR_head_value));print(class(input1$CIGAR_tail_value))
head.p <- ggplot(input1,aes(x=CIGAR_head_value))+geom_density(adjust = reso,linewidth=wi)+
  theme_bw()+theme(text = element_text(size=st))+geom_vline(xintercept=misM,linetype=5,col="grey",linewidth=wi)+
  facet_grid(strand_rev ~ .,scales = "free")+
  labs(title=paste0(sample,": head soft clip of LAP"))
tail.p <- ggplot(input2,aes(x=CIGAR_tail_value))+geom_density(adjust = reso,linewidth=wi)+
  theme_bw()+theme(text = element_text(size=st))+geom_vline(xintercept=misM,linetype=5,col="grey",linewidth=wi)+
  facet_grid(strand_rev ~ .,scales = "free")+
  labs(title=paste0(sample,": tail soft clip of LAP"))
#ggsave(plot=head.p+tail.p,paste0(out_dir,"_CIGAR_distrib_",sample,".pdf"),w=6.21,h=4.26)
pdf(paste0(out_dir, "_CIGAR_distrib.pdf"), w = 6.21, h = 4.26)
tryCatch({
  print(head.p)
}, error = function(e) {
  message("Warning: Failed to plot head.p due to error: ", e$message)
})
tryCatch({
  print(tail.p)
}, error = function(e) {
  message("Warning: Failed to plot tail.p due to error: ", e$message)
})
dev.off()

#-----------------------------------------
# assign LAPs to PAS no matter which PASref they mapped to
#max.dist = 500000000
#refPAS_file <- "/wistar/tian_publicdata/PolyA_DBv3.2/human/human.PAS.hg38.txt"
## polyAdb
anno <- read.delim(file = refPAS_file, header = TRUE)
# Support PolyA_DB v4.1: PAS_ID encodes "chr:strand:pos" with no separate columns.
# v3.2 had explicit Chromosome/Strand/Position columns; parse from PAS_ID if absent.
if (!all(c("Chromosome", "Position", "Strand") %in% colnames(anno))) {
  pas_parts <- do.call(rbind, strsplit(as.character(anno$PAS_ID), ":"))
  anno$Chromosome <- pas_parts[, 1]
  anno$Strand     <- pas_parts[, 2]
  anno$Position   <- as.numeric(pas_parts[, 3])
}
anno$Position_end <- anno$Position+1
## GenomicRanges object is 1-based!
GR.polyA.db = GenomicRanges::makeGRangesFromDataFrame(anno, keep.extra.columns = TRUE,starts.in.df.are.0based = TRUE,
    seqnames.field = "Chromosome", start.field = "Position",
    end.field = "Position_end", strand.field = "Strand")

## LAP
df2$LAP_end <- df2$LAP+1
ranges = GenomicRanges::makeGRangesFromDataFrame(df2, keep.extra.columns = TRUE,starts.in.df.are.0based = TRUE,
    seqnames.field = "mapped_chr", start.field = "LAP",
    end.field = "LAP_end", strand.field = "strand_rev")

OL = suppressWarnings(GenomicRanges::distanceToNearest(x = ranges,
    subject = GR.polyA.db, ignore.strand = FALSE))

## extract the nearest PAS position
query_end <- end(ranges[queryHits(OL)])
subject_end <- end(GR.polyA.db[subjectHits(OL)])
# annotate upstrem or downstream
distance <- subject_end - query_end
strand_info <- as.character(strand(ranges[queryHits(OL)]))
distance_signed <- ifelse(strand_info == "+", distance, -distance)
# add to mcols
mcols(OL)$distance_signed <- distance_signed

# distance：
ol <- as.data.frame(OL)
#  queryHits subjectHits distance
#1         1       45851        4
#2         2       45853       95
#5         7      228531     2179

# Create a row number column for matching
df2$index <- seq_len(nrow(df2))
# matching
meta.new <- merge(df2, ol, by.x = "index", by.y = "queryHits", all.x = TRUE,sort=F)
colnames(anno) <- paste0("hit_",colnames(anno))
anno$index <- seq_len(nrow(anno))
meta.new <- merge(meta.new, anno, by.x = "subjectHits", by.y="index", all.x = TRUE,sort=F)
meta.new <- meta.new[, !names(meta.new) %in% c("subjectHits", "index")]
saveRDS(meta.new,paste0(out_dir,"_LAP_polyAdb.rds"))

LAPs <- subset(meta.new, LAP_anno=="LAP")
stats[length(stats)+1] <- nrow(LAPs)
stats[length(stats)+1] <- stats[length(stats)]/stats[1]*100

stats[length(stats)+1] <- nrow(subset(LAPs, mapped_chr != "chrM"))
#stats[length(stats)+1] <- nrow(subset(LAPs, mapped_chr %in% unique(anno$Chromosome)))
stats[length(stats)+1] <- nrow(subset(LAPs, ! is.na(hit_PAS_ID)))

sum(is.na(LAPs$hit_PAS_ID))

## plot distance distribution
get_mode <- function(x) {
  freq_table <- table(x)
  mode_values <- as.numeric(names(freq_table[freq_table == max(freq_table)]))
  return(mode_values)
}
mode = get_mode(LAPs$distance) # 0
print(paste0("distance mode: ",mode," and the freq is ",nrow(subset(LAPs,distance==mode))))
print(summary(LAPs$distance))

#ggplot(LAPs,aes(x=distance))+geom_density(adjust = 0.1)+geom_vline(xintercept=c(mode, dist),linetype=5,col="grey")+theme_bw()+#xlim(1,1000)+
#labs(title=paste0(sample,": Distance between LAP and PAS"))
LAPs$distance_signed <- as.numeric(LAPs$distance_signed)
p1=ggplot(LAPs,aes(x=distance_signed))+geom_density(adjust = 0.1)+geom_vline(xintercept=c(mode,-dist,dist),linetype=5,col="grey")+theme_bw()+#xlim(1,1000)+
labs(title=paste0(sample,": Distance of LAP to PAS"))
p2=ggplot(LAPs,aes(x=distance_signed))+geom_density(adjust = 0.1)+geom_vline(xintercept=c(mode,-dist,dist),linetype=5,col="grey")+theme_bw()+xlim(-100,100)+
labs(title=paste0(sample,": Distance of LAP to PAS"))
pdf(paste0(out_dir,"_LAP_polyAdb_distance.pdf"))
tryCatch({
  print(p1)
}, error = function(e) {
  message("Warning: Failed to plot p1 due to error: ", e$message)
})
tryCatch({
  print(p2)
}, error = function(e) {
  message("Warning: Failed to plot p2 due to error: ", e$message)
})
dev.off()

## filter by dist nt
LAPs_dist <- subset(LAPs, distance < dist)
stats[length(stats)+1] <- nrow(LAPs_dist)
stats[length(stats)+1] <- stats[length(stats)]/stats[1]*100

saveRDS(LAPs_dist,paste0(out_dir,"_LAPs_",dist,".PASS.rds"))

cluster.all.reads <- data.frame(table(LAPs_dist$hit_PAS_ID))
names(cluster.all.reads) <- c("hit_PAS_ID",sample)
write.csv(cluster.all.reads, paste0(out_dir,"_cluster.all.reads.csv"), row.names = F, quote = F)

names(stats) <- c("mapped_reads","LAP_reads","LAP_pct","non-chrM","PAS_matched",paste0("LAP_f",dist),paste0("LAP_f",dist,"_pct"))
stats.df <- data.frame(stats)
colnames(stats.df) <- paste0(sample,"_stats")
write.csv(stats.df, paste0(out_dir,"_stats.csv"), row.names = T, quote = F)

PASS_bw_bed <- LAPs_dist[c("mapped_chr", "mapped_start", "mapped_end", "readID", "MAPQ", "strand_rev")]
options(scipen = 999)
write.table(PASS_bw_bed, paste0(out_dir,"_PASS_bw_LAP",dist,".bed"), col.names = F, row.names = F, sep = "\t", quote = F)

PASS_bw_bed_point <- LAPs_dist[c("mapped_chr", "LAP", "LAP_end", "readID", "MAPQ", "strand_rev")]
options(scipen = 999)
write.table(PASS_bw_bed_point, paste0(out_dir,"_PASS_bw_LAP",dist,"_positon.bed"), col.names = F, row.names = F, sep = "\t", quote = F)

LAPs_dist$PAS_end <- LAPs_dist$hit_Position+1
PASS_bw_bed_PAS <- LAPs_dist[c("hit_Chromosome", "hit_Position", "PAS_end", "hit_PAS_ID", "MAPQ", "hit_Strand")]
options(scipen = 999)
write.table(PASS_bw_bed_PAS, paste0(out_dir,"_PASS_bw_LAP",dist,"_PAS.bed"), col.names = F, row.names = F, sep = "\t", quote = F)

mapped_read_bw_bed <- df2[c("mapped_chr", "mapped_start", "mapped_end", "readID", "MAPQ", "strand_rev")]
options(scipen = 999)
write.table(mapped_read_bw_bed, paste0(out_dir,"_mapped_read.bed"), col.names = F, row.names = F, sep = "\t", quote = F) 

print("done!")