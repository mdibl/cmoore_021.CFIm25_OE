# ============================================================
# DEXSeq-based Alternative Polyadenylation (APA) Analysis
# Generic script — configure the USER CONFIGURATION section only.
#
# Usage:
#   Set BASE_DIR, input file paths, and AXES. Each entry in AXES runs one
#   DEXSeq comparison (comparison_var: ref_label vs alt_label) independently
#   within every level of stratify_var -- e.g. by_condition compares condition
#   (control vs TGFB) within each plasmid, while by_plasmid compares plasmid
#   (control vs CFIm25) within each condition.
#   Add/remove entries in AXES to add more cross-sections.
#
# Input:
#   - PAS count matrix (rows = PAS IDs, columns = samples)
#   - PAS annotation file (PolyA_DB v3.2 or v4.1 auto-detected)
#   - Design file (tab-delimited; rows = samples, columns = metadata)
#   - (optional) sample name list
# ============================================================

# -- Libraries ---------------------------------------------------------------
library(DEXSeq)
library(SummarizedExperiment)
library(GenomicRanges)
library(DESeq2)
library(ggplot2)
library(ggrepel)
library(scales)
library(limma)
library(dplyr)
library(tidyr)
library(tibble)
library(matrixStats)
library(stringr)

# ============================================================
#  USER CONFIGURATION — EDIT THESE
# ============================================================

BASE_DIR <- "/compbio/analysis/ClaireMoore/cmoore_021.CFIm25_OE"   # project root

# Input files (relative to BASE_DIR)
COUNTS_CSV  <- "data/cluster.all.reads.csv"
SAMPLES_TXT <- "data/sample_name.txt"       # set to NULL to use all count_mat columns
ANNO_FILE   <- "data/PolyA_DB_v4.1/hg38.PAS.main.tsv"
DESIGN_FILE <- "data/design.txt"

# RefSeq GTF (same one used to build the STAR index for this project) --
# provides real exon/intron structure for the APA genome-map figures below.
# Set to NULL to skip those figures entirely (e.g. if no GTF is available).
ANNO_GTF <- "data/ref/hg38.refGene.gtf"

# Output base directory
OUT_BASE <- file.path(BASE_DIR, "output/DEXseq/")

# Comparison axes -- each axis runs one DEXSeq comparison (comparison_var: ref_label
# vs alt_label) independently within every level of stratify_var (or once overall if
# stratify_var's level is left unspecified). Add/remove axes here to add more
# cross-sections without touching any code below.
#   comparison_var/stratify_var — colData column names
#   ref_label/alt_label         — the two comparison_var levels to test (ref = fold-
#                                  change denominator)
#   comparison_label            — display label for plot legends/titles
#   stratify_levels              — NULL = auto-detect from colData (alphabetical order)
#   stratify_dir                 — output subdirectory name per stratify level;
#                                  NULL = use the level name as-is
AXES <- list(
  by_condition = list(                     
    comparison_var   = "condition",
    ref_label         = "control",
    alt_label         = "TGFB",
    comparison_label  = "condition",
    stratify_var      = "plasmid",
    stratify_levels   = c("control", "CFIm25"),
    stratify_dir      = NULL
  ),
  by_plasmid = list(                     
    comparison_var    = "plasmid",
    ref_label         = "control",
    alt_label         = "CFIm25",
    comparison_label  = "plasmid",
    stratify_var      = "condition",
    stratify_levels   = c("control", "TGFB"),
    stratify_dir      = NULL
  )
)

# PAS type filter — keep rows whose PAS_type matches this regex
PAS_TYPE_REGEX <- "^3'UTR"

# Analysis thresholds
MIN_TOTAL_READS   <- 10    # min "this"-read row-sum per group to keep a PAS (hard filter)
MIN_PER_CONDITION <- 1     # min samples with count > 0 per condition; 0 = disable
PADJ_CUT          <- 0.05  # FDR threshold for significance calls
LFC_CUT           <- 0.5   # |log2FC| threshold for "big" hits in gene summary

# thin_counts is a warning flag, not a filter -- PAS below MIN_TOTAL_READS are
# already excluded above; this just flags surviving PAS that are still thin
# (e.g. cleared 10 total reads but sit well under a stricter 40-read bar).
THIN_COUNTS_CUT <- 40

# single_dominant_site (see build_gene_apa_summary()): a gene's most-abundant
# PAS must carry at least DOMINANT_USAGE_CUT of usage in BOTH conditions, and
# move by less than DOMINANT_STABLE_CUT between them, to count as "the" site.
DOMINANT_USAGE_CUT  <- 0.90
DOMINANT_STABLE_CUT <- 0.05

# RED (Relative Expression Difference) -- our N-site-generalized version of the RED
# metric used from MAAPER (which is limited to exactly 2 sites/gene). A small
# pseudocount added to both n(i) and m(i) (Haldane-Anscombe style) so a PAS at 0% or
# 100% usage in a condition still yields a large, finite, rankable value instead of
# +-Inf/NA -- exactly the dramatic on/off shifts a ranking rubric most cares about.
RED_PSEUDOCOUNT <- 0.5

# Candidate-gene rubric (see build_candidate_genes()) -- our replacement for the old
# MAAPER/RED-based candidate-selection rubric. A gene's single biggest-moving PAS
# must clear ALL of: usage change, DEXSeq significance (reuses PADJ_CUT below), read
# depth (reuses thin_counts), 3' UTR length shift, and RED magnitude -- and the gene
# must not be a single_dominant_site gene (a shift confined to a minor site there
# isn't the kind of finding this rubric is meant to surface).
USAGE_CHANGE_CUT <- 0.05  # min |meanUsage_Treatment - meanUsage_Control| for the top-moving PAS
WUTR_CHANGE_CUT  <- 50    # min |delta_wUTR| (bp) for the gene overall
RED_MAG_CUT      <- 1     # min |delta_RED| (log2 scale) for the top-moving PAS

# PCA settings
NTOP_PCA      <- 2000        # top-N variable PAS for PCA (by row variance)
PCA_LABEL_COL <- "replicate" # colData column to use as point labels (NULL = none)

# Genes to plot APA genome-map figures for; set to character(0) to skip
# Cherry picked genes (by me):
GENES_OF_INTEREST <- c("PCF11", "TAB2", "ICAM1", "MCAM", "RCOR3", "OTUD7B", "PCNX4", "LATS2", "CENPA")

# RUVseq batch correction
#   USE_RUV <- TRUE  when replicates cluster by batch in PCA rather than condition.
#   RUVs estimates W factors from within-condition replicate variation; W_1:exon is
#   then added to the DEXSeq design to correct for batch effects on relative PAS usage.
#   RUV_K: number of factors to estimate (1 is almost always sufficient; raise to 2
#   only if the first factor does not account for the outlier clustering).
#   RUV_BATCH_VAR: colData column standing in for batch/processing-day structure --
#   used only by the covariate report below (checks whether each W factor is
#   explaining this batch proxy rather than real biology).
USE_RUV       <- TRUE
RUV_K         <- 1
RUV_BATCH_VAR <- "replicate"

# Cache the DEXSeq/DESeq2 model fits (the slow part of this script) to disk,
# keyed on input-file mtimes + the modeling config above. Re-running the
# script to change a plot color or GENES_OF_INTEREST reuses the cached fit
# instead of refitting; changing the count/annotation/design files or any of
# the modeling settings above automatically invalidates it. Set to FALSE to
# always refit from scratch (e.g. while debugging the modeling code itself).
USE_CACHE <- TRUE

# ============================================================
#  DERIVED PATHS — do not edit below this line
# ============================================================

path_counts  <- file.path(BASE_DIR, COUNTS_CSV)
path_anno    <- file.path(BASE_DIR, ANNO_FILE)
path_design  <- file.path(BASE_DIR, DESIGN_FILE)
path_samples <- if (!is.null(SAMPLES_TXT)) file.path(BASE_DIR, SAMPLES_TXT) else NULL
path_gtf     <- if (!is.null(ANNO_GTF)) file.path(BASE_DIR, ANNO_GTF) else NULL

# Sequential color ramp for usage magnitude (validated blue single-hue ramp;
# see dataviz skill references/palette.md steps 100/700)
SEQ_LOW  <- "#cde2fb"
SEQ_HIGH <- "#0d366b"

# Validated categorical pair for Control/Experimental overlay lines (dataviz
# skill references/palette.md slots 1 & 2; passes CVD + normal-vision checks)
CAT_CTRL <- "#2a78d6"
CAT_TRT  <- "#eb6834"

# Diverging pair for gene-level DGE direction (dataviz skill: blue<->red poles,
# distinct from the categorical Control/Experimental colors above so "which
# condition" and "which direction" never collide in the same figure)
DIV_UP   <- "#e34948"
DIV_DOWN <- "#0d366b"
DIV_NS   <- "grey60"

dir_results     <- OUT_BASE
dir_plots       <- file.path(OUT_BASE, "plots")
dir_qc          <- file.path(OUT_BASE, "plots/qc")
dir_apa_genome  <- file.path(OUT_BASE, "plots/apa_genome")
dir_apa_zoom    <- file.path(OUT_BASE, "plots/apa_zoom")
dir_cache       <- file.path(OUT_BASE, "cache")
for (d in c(dir_results, dir_plots, dir_qc, dir_apa_genome, dir_apa_zoom, dir_cache)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

# Per-axis derived column names (fold-change, usage, RED, wUTR, etc.) are computed
# inline from each axis's ref_label/alt_label wherever needed (see run_axis()),
# rather than as fixed globals -- there are now two axes with different labels.

# ============================================================
#  BUILD COUNT MATRIX
# ============================================================

message("Reading count matrix...")
counts_df <- read.csv(path_counts, check.names = FALSE)
stopifnot("hit_PAS_ID" %in% names(counts_df))

if (!is.null(path_samples)) {
  samples <- scan(path_samples, what = character(), quiet = TRUE)
} else {
  samples <- setdiff(names(counts_df), "hit_PAS_ID")
}
stopifnot(all(samples %in% names(counts_df)))

count_mat_raw <- as.matrix(counts_df[, samples])            # unfiltered — used for library sizes
rownames(count_mat_raw) <- counts_df$hit_PAS_ID
storage.mode(count_mat_raw) <- "integer"

count_mat <- count_mat_raw   # will be row-filtered after annotation alignment below

# ============================================================
#  BUILD PAS ANNOTATION
# ============================================================

message("Reading PAS annotation...")
pas_anno <- read.delim(path_anno, stringsAsFactors = FALSE, check.names = FALSE)

# Auto-detect PolyA_DB format:
#   v3.2: explicit Chromosome / Position / Strand columns
#   v4.1: all three encoded in PAS_ID as "chr:strand:pos"
has_coord_cols <- all(c("Chromosome", "Position", "Strand") %in% names(pas_anno))

if (!has_coord_cols) {
  if ("PAS_ID" %in% names(pas_anno)) {
    message("PolyA_DB v4.1 format detected: parsing coordinates from PAS_ID.")
    parsed              <- strsplit(as.character(pas_anno$PAS_ID), ":", fixed = TRUE)
    pas_anno$Chromosome <- vapply(parsed, `[`, character(1), 1)
    pas_anno$Strand     <- vapply(parsed, `[`, character(1), 2)
    pas_anno$Position   <- as.integer(vapply(parsed, `[`, character(1), 3))
  } else {
    message("No coordinate columns found — GRanges will not be built; APA direction will be NA.")
  }
}

pas_anno <- pas_anno[grepl(PAS_TYPE_REGEX, pas_anno$PAS_type), ]

blank_to_na <- function(x) {
  x <- stringr::str_squish(x)
  x[tolower(x) %in% c("", "na")] <- NA_character_
  x
}

id_cols <- intersect(c("Gene Symbol", "GeneSymbol", "Ensembl ID", "RefSeq Gene ID", "FANTOM ID"),
                     names(pas_anno))
if (length(id_cols) == 0) stop("No recognized gene-ID columns in annotation file.")
for (col in id_cols) pas_anno[[col]] <- blank_to_na(pas_anno[[col]])
pas_anno$groupID <- do.call(dplyr::coalesce, as.list(pas_anno[, id_cols, drop = FALSE]))

# Rescue PAS IDs that are in the count matrix but have no gene name
unmatched_ids <- setdiff(counts_df$hit_PAS_ID, pas_anno$PAS_ID[!is.na(pas_anno$groupID)])
rescue <- pas_anno$PAS_ID %in% unmatched_ids & is.na(pas_anno$groupID)
pas_anno$groupID[rescue] <- pas_anno$PAS_ID[rescue]
pas_anno <- pas_anno[!is.na(pas_anno$groupID), ]

if (anyDuplicated(pas_anno[c("PAS_ID", "groupID")])) {
  pas_anno <- pas_anno[!duplicated(pas_anno[c("PAS_ID", "groupID")]), ]
}

# ============================================================
#  ALIGN MATRIX & ANNOTATION; SET groupID / featureID
# ============================================================

common_pas <- intersect(rownames(count_mat), pas_anno$PAS_ID)
count_mat  <- count_mat[common_pas, , drop = FALSE]
pas_anno   <- pas_anno[match(common_pas, pas_anno$PAS_ID), ]

groupID   <- pas_anno$groupID
# PolyA_DB v4.1 PAS_IDs contain ':' (chr:strand:pos). DEXSeq strips ':' from
# featureIDs, which would desync GRanges names from internal rownames and crash.
# Sanitize here so DEXSeq sees IDs it doesn't need to modify.
featureID           <- gsub("[: ]", "_", pas_anno$PAS_ID)
rownames(count_mat) <- featureID          # count_mat rows are in pas_anno order here
pas_anno$featureID_safe <- featureID      # retained for annotation merge below

# ============================================================
#  BUILD GRanges (enables APA direction labelling)
# ============================================================

if (all(c("Chromosome", "Position", "Strand") %in% names(pas_anno))) {
  gr <- GRanges(
    seqnames = pas_anno$Chromosome,
    ranges   = IRanges(start = as.integer(pas_anno$Position), width = 1L),
    strand   = pas_anno$Strand
  )
  # DEXSeq composite key — allows gr[rrn] to subset correctly after per-group filtering
  names(gr) <- paste(groupID, featureID, sep = ":")
} else {
  gr <- NULL
}

# ============================================================
#  LOAD TRANSCRIPT ANNOTATION (GTF) — optional, for APA genome-map figures
# ============================================================
# Provides real exon/intron structure (PolyA_DB alone has no gene structure,
# only flat PAS positions). Skipped gracefully -- with the two figures below
# simply not produced -- if ANNO_GTF is NULL, the file is missing, or
# rtracklayer isn't installed; nothing else in the pipeline depends on it.

gtf_exons <- NULL
if (!is.null(path_gtf)) {
  if (!file.exists(path_gtf)) {
    warning("ANNO_GTF set but file not found at ", path_gtf,
            " — APA genome-map figures will be skipped.")
  } else if (!requireNamespace("rtracklayer", quietly = TRUE)) {
    warning("rtracklayer not installed — APA genome-map figures will be skipped. ",
            "Install with: BiocManager::install('rtracklayer')")
  } else {
    message("Reading transcript annotation (GTF) for APA genome-map figures...")
    gtf_all   <- rtracklayer::import(path_gtf)
    gtf_exons <- as.data.frame(gtf_all[gtf_all$type == "exon"])
    gtf_exons$transcript_id <- as.character(gtf_exons$transcript_id)
    gtf_exons$strand        <- as.character(gtf_exons$strand)
    rm(gtf_all)
  }
}

# ============================================================
#  LOAD DESIGN / colData
# ============================================================

colData <- read.delim(path_design, row.names = 1, check.names = FALSE)
stopifnot(all(colnames(count_mat) == rownames(colData)))

# Validates every AXES entry's comparison_var/stratify_var exist in colData, levels
# comparison_var as c(ref_label, alt_label) (the reference level fold-changes are
# measured against), and auto-detects+backfills stratify_levels (alphabetical) when
# an axis leaves it NULL. A colData column referenced by two axes (e.g. "condition"
# is by_condition's comparison_var AND by_plasmid's stratify_var) gets leveled once
# per axis that references it -- harmless as long as every axis agrees on that
# column's level order, which AXES above does.
prepare_colData_axes <- function(colData, axes) {
  for (nm in names(axes)) {
    ax <- axes[[nm]]
    for (v in c(ax$comparison_var, ax$stratify_var)) {
      if (!v %in% names(colData))
        stop(sprintf("AXES$%s references colData column '%s', which is not in the design file.", nm, v))
    }
    colData[[ax$comparison_var]] <- factor(colData[[ax$comparison_var]],
                                            levels = c(ax$ref_label, ax$alt_label))
    strat_levels <- ax$stratify_levels
    if (is.null(strat_levels)) strat_levels <- sort(unique(as.character(colData[[ax$stratify_var]])))
    axes[[nm]]$stratify_levels <- strat_levels
    colData[[ax$stratify_var]] <- factor(colData[[ax$stratify_var]], levels = strat_levels)
  }
  list(colData = colData, axes = axes)
}
prepared <- prepare_colData_axes(colData, AXES)
colData  <- prepared$colData
AXES     <- prepared$axes
rm(prepared)

# Second design factor to show as point shape on the global pre-analysis PCA
# (which colors by "condition"). Derived from AXES rather than hardcoded to a
# specific column name -- picks the stratify_var of whichever axis uses
# "condition" as its comparison_var, so a project whose second factor isn't
# called "timepoint" (e.g. "plasmid") still gets that contrast shown, instead
# of the shape aesthetic silently matching nothing.
condition_axis <- Find(function(ax) ax$comparison_var == "condition", AXES)
PCA_SHAPE_VAR  <- if (!is.null(condition_axis)) condition_axis$stratify_var else NULL

# ============================================================
#  RUVseq BATCH CORRECTION (optional — set USE_RUV <- TRUE)
# ============================================================

if (USE_RUV) {
  if (!requireNamespace("RUVSeq", quietly = TRUE))
    stop("RUVSeq is required when USE_RUV = TRUE.  Install with: BiocManager::install('RUVSeq')")
  library(RUVSeq)

  message("Estimating unwanted variation with RUVs (k = ", RUV_K, ")...")

  # scIdx: each row is a set of within-condition replicates — any variation among
  # them is treated as unwanted.  Pad uneven groups with -1.
  # (Depends on some axis in AXES using "condition" as its comparison_var --
  # prepare_colData_axes() above is what levels colData$condition ref-label-first.)
  condition_levels <- levels(colData$condition)
  ctrl_idx <- which(colData$condition == condition_levels[1])
  trt_idx  <- which(colData$condition == condition_levels[2])
  nc       <- max(length(ctrl_idx), length(trt_idx))
  # scIdx must be a numeric (double) matrix; -1 marks empty padding slots
  scIdx         <- matrix(-1, nrow = 2, ncol = nc)
  scIdx[1, seq_along(ctrl_idx)] <- ctrl_idx
  scIdx[2, seq_along(trt_idx)]  <- trt_idx

  ruv_set <- EDASeq::newSeqExpressionSet(
    counts    = count_mat,
    phenoData = data.frame(condition = colData$condition, row.names = rownames(colData))
  )
  # SeqExpressionSet method requires cIdx = "character" (feature names, not indices)
  ruv_fit <- RUVSeq::RUVs(ruv_set,
                           cIdx  = rownames(count_mat),
                           k     = as.numeric(RUV_K),
                           scIdx = scIdx)

  # Attach W factors to colData so run_dexseq_group picks them up automatically
  w_cols <- paste0("W_", seq_len(RUV_K))
  for (wc in w_cols) {
    colData[[wc]] <- pData(ruv_fit)[rownames(colData), wc]
  }
  message("  W factors attached to colData: ", paste(w_cols, collapse = ", "))

  # ---- Covariate reporting: is W absorbing real biology, or real batch? ----
  # For each W_k: a low p-value against any axis's biological comparison_var
  # (condition, plasmid, ...) is a red flag -- W is eating a real treatment
  # effect, not unwanted variation. A high R^2 against RUV_BATCH_VAR is the
  # reassuring case -- W is explaining known batch/processing structure.
  stopifnot("RUV_BATCH_VAR column not found in design file" = RUV_BATCH_VAR %in% names(colData))

  # One-way lm() of a W factor against a grouping column -- reused for both checks
  # below (a biological comparison_var, or the batch proxy). Works for 2-level and
  # multi-level factors alike (anova() F-test generalizes; a 2-level factor's F-test
  # p-value is equivalent to the usual two-sample t-test p-value).
  covariate_assoc <- function(w, grp) {
    grp <- factor(grp)
    if (nlevels(grp) < 2 || any(table(grp) < 2)) return(c(pvalue = NA_real_, r_squared = NA_real_))
    fit <- lm(w ~ grp)
    c(pvalue = anova(fit)[["Pr(>F)"]][1], r_squared = summary(fit)$r.squared)
  }

  # Every distinct biological comparison_var across AXES (e.g. "condition", "plasmid")
  # -- W ~ these should come back NS; W ~ RUV_BATCH_VAR should come back high-R^2.
  biological_vars <- unique(vapply(AXES, function(ax) ax$comparison_var, character(1)))

  covariate_report <- do.call(rbind, lapply(w_cols, function(wc) {
    w   <- colData[[wc]]
    row <- list(W = wc)
    for (bv in biological_vars) {
      row[[paste0(bv, "_pvalue")]] <- covariate_assoc(w, colData[[bv]])[["pvalue"]]
    }
    row[[paste0(RUV_BATCH_VAR, "_R2")]] <- covariate_assoc(w, colData[[RUV_BATCH_VAR]])[["r_squared"]]
    as.data.frame(row, stringsAsFactors = FALSE)
  }))

  message("  Covariate report (NS pvalue = W isn't absorbing biology; high R2 = W is correcting batch):")
  for (i in seq_len(nrow(covariate_report))) {
    message(sprintf("    %s", paste(names(covariate_report), unlist(covariate_report[i, ]),
                                     sep = "=", collapse = "  ")))
  }
  write.csv(covariate_report, file.path(dir_qc, "RUV_covariate_report.csv"), row.names = FALSE, quote = FALSE)

  # Raw W factors (samples x k), for anyone who wants to inspect/reuse them directly
  w_export <- tibble::rownames_to_column(as.data.frame(colData[, w_cols, drop = FALSE]), "sample")
  write.csv(w_export, file.path(dir_qc, "RUV_W_factors.csv"), row.names = FALSE, quote = FALSE)
}

# ============================================================
#  FUNCTION: global PCA on raw VST PAS counts (pre-DEXSeq)
# ============================================================

make_pca_global <- function(count_mat, col_data,
                             color_col      = "condition",
                             shape_col      = NULL,
                             label_col      = PCA_LABEL_COL,
                             ntop           = NTOP_PCA,
                             batch_covars   = NULL,
                             outfile        = NULL,
                             matrix_outfile = NULL) {
  dds     <- DESeq2::DESeqDataSetFromMatrix(countData = count_mat,
                                             colData  = col_data,
                                             design   = ~ 1)
  vst_mat <- SummarizedExperiment::assay(DESeq2::vst(dds, blind = TRUE))

  # Export the full (pre-ntop-filter, pre-batch-correction) VST matrix on request --
  # this is the "real" normalized data; batch correction below is visualization-only.
  if (!is.null(matrix_outfile)) {
    vst_out <- tibble::rownames_to_column(as.data.frame(vst_mat), "featureID")
    write.csv(vst_out, matrix_outfile, row.names = FALSE, quote = FALSE)
  }

  # Optionally remove batch covariates from VST matrix for visualization only.
  # This does not affect the statistical model — it is purely for the PCA plot.
  if (!is.null(batch_covars) && length(batch_covars) > 0) {
    stopifnot(all(batch_covars %in% colnames(col_data)))
    vst_mat <- limma::removeBatchEffect(
      vst_mat,
      covariates = as.matrix(col_data[, batch_covars, drop = FALSE])
    )
  }

  if (!is.null(ntop) && ntop > 0 && nrow(vst_mat) > ntop) {
    rv  <- matrixStats::rowVars(vst_mat)
    vst_mat <- vst_mat[order(rv, decreasing = TRUE)[seq_len(ntop)], , drop = FALSE]
  }

  p   <- prcomp(t(vst_mat), center = TRUE, scale. = FALSE)
  pct <- round(100 * (p$sdev^2 / sum(p$sdev^2))[1:2], 1)

  pcs <- as.data.frame(p$x[, 1:2]) |>
    tibble::rownames_to_column("sample_id") |>
    dplyr::left_join(tibble::rownames_to_column(as.data.frame(col_data), "sample_id"),
                     by = "sample_id")

  aes_base <- if (!is.null(color_col) && color_col %in% names(pcs))
    ggplot2::aes(x = PC1, y = PC2, color = .data[[color_col]])
  else
    ggplot2::aes(x = PC1, y = PC2)

  pca_title <- if (!is.null(batch_covars))
    "Global PCA — VST 3'UTR PAS counts (RUV-corrected, visual only)"
  else
    "Global PCA — VST 3'UTR PAS counts"

  g <- ggplot2::ggplot(pcs, aes_base) +
    ggplot2::geom_point(size = 3) +
    ggplot2::labs(x = paste0("PC1 (", pct[1], "%)"),
                  y = paste0("PC2 (", pct[2], "%)"),
                  title = pca_title) +
    ggplot2::theme_minimal(base_size = 12)

  if (!is.null(shape_col) && shape_col %in% names(pcs)) {
    g <- g + ggplot2::aes(shape = .data[[shape_col]])
  }
  if (!is.null(label_col) && label_col %in% names(pcs)) {
    g <- g + ggrepel::geom_text_repel(ggplot2::aes(label = .data[[label_col]]), size = 3)
  }

  if (!is.null(outfile)) ggplot2::ggsave(outfile, g, width = 7, height = 5, dpi = 300)
  invisible(list(plot = g, pcs = pcs, prcomp = p))
}

# ============================================================
#  FUNCTION: library size bar chart
# ============================================================

# Intentionally hardcoded to "condition" rather than axis-parameterized -- this is a
# single global pre-analysis QC plot (raw library depth), not specific to any axis.
make_library_sizes <- function(count_mat_raw, col_data, outfile = NULL) {
  df <- data.frame(
    sample    = colnames(count_mat_raw),
    lib_size  = colSums(count_mat_raw),
    condition = as.character(col_data[colnames(count_mat_raw), "condition"]),
    stringsAsFactors = FALSE
  )
  df <- df[order(df$condition, df$sample), ]
  df$sample <- factor(df$sample, levels = df$sample)

  g <- ggplot2::ggplot(df, ggplot2::aes(x = sample, y = lib_size / 1e6, fill = condition)) +
    ggplot2::geom_col() +
    ggplot2::labs(x = "Sample", y = "Total PAS reads (M)",
                  title = "Library sizes (total PAS read counts)") +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))

  if (!is.null(outfile)) {
    w <- max(4, 0.5 * nrow(df) + 2)
    ggplot2::ggsave(outfile, g, width = w, height = 4, dpi = 300)
  }
  invisible(list(plot = g, data = df))
}

# ============================================================
#  FUNCTION: per-group PCA using DEXSeqDataSet (post-normalization)
# ============================================================

make_pca_dxd <- function(dxd, grp_label = "",
                          color_col      = "condition",
                          label_col      = PCA_LABEL_COL,
                          ntop           = NTOP_PCA,
                          outfile        = NULL,
                          matrix_outfile = NULL) {
  cd_full <- as.data.frame(SummarizedExperiment::colData(dxd))
  this_idx <- which(cd_full$exon == "this")
  cd  <- cd_full[this_idx, , drop = FALSE]
  mat <- SummarizedExperiment::assay(dxd)[, this_idx, drop = FALSE]

  dds     <- DESeq2::DESeqDataSetFromMatrix(countData = mat, colData = cd, design = ~ 1)
  vst_mat <- SummarizedExperiment::assay(DESeq2::vst(dds, blind = TRUE))

  if (!is.null(matrix_outfile)) {
    vst_out <- tibble::rownames_to_column(as.data.frame(vst_mat), "featureID")
    write.csv(vst_out, matrix_outfile, row.names = FALSE, quote = FALSE)
  }

  if (!is.null(ntop) && ntop > 0 && nrow(vst_mat) > ntop) {
    rv  <- matrixStats::rowVars(vst_mat)
    vst_mat <- vst_mat[order(rv, decreasing = TRUE)[seq_len(ntop)], , drop = FALSE]
  }

  p   <- prcomp(t(vst_mat), center = TRUE, scale. = FALSE)
  pct <- round(100 * (p$sdev^2 / sum(p$sdev^2))[1:2], 1)

  pcs <- as.data.frame(p$x[, 1:2]) |>
    tibble::rownames_to_column("sample_id") |>
    dplyr::left_join(tibble::rownames_to_column(cd, "sample_id"), by = "sample_id")

  aes_base <- if (!is.null(color_col) && color_col %in% names(pcs))
    ggplot2::aes(x = PC1, y = PC2, color = .data[[color_col]])
  else
    ggplot2::aes(x = PC1, y = PC2)

  ttl <- paste0("PCA — normalized PAS counts",
                if (nchar(grp_label) > 0) paste0(" (", grp_label, ")") else "")

  g <- ggplot2::ggplot(pcs, aes_base) +
    ggplot2::geom_point(size = 3) +
    ggplot2::labs(x = paste0("PC1 (", pct[1], "%)"),
                  y = paste0("PC2 (", pct[2], "%)"),
                  title = ttl) +
    ggplot2::theme_minimal(base_size = 12)

  if (!is.null(label_col) && label_col %in% names(pcs)) {
    g <- g + ggrepel::geom_text_repel(ggplot2::aes(label = .data[[label_col]]), size = 3)
  }

  if (!is.null(outfile)) ggplot2::ggsave(outfile, g, width = 7, height = 5, dpi = 300)
  invisible(list(plot = g, pcs = pcs))
}

# ============================================================
#  FUNCTION: size factor bar chart
# ============================================================

make_size_factors <- function(dxd, grp_label = "", comparison_var = "condition", outfile = NULL) {
  cd_full  <- as.data.frame(SummarizedExperiment::colData(dxd))
  this_idx <- which(cd_full$exon == "this")
  cd       <- cd_full[this_idx, , drop = FALSE]
  sf       <- sizeFactors(dxd)[this_idx]

  df <- data.frame(
    sample      = as.character(cd$sample),
    size_factor = as.numeric(sf),
    condition   = as.character(cd[[comparison_var]]),
    stringsAsFactors = FALSE
  )
  df <- df[order(df$condition, df$sample), ]
  df$sample <- factor(df$sample, levels = df$sample)

  g <- ggplot2::ggplot(df, ggplot2::aes(x = sample, y = size_factor, fill = condition)) +
    ggplot2::geom_col() +
    ggplot2::geom_hline(yintercept = 1, linetype = "dashed", color = "gray50") +
    ggplot2::labs(
      x     = "Sample",
      y     = "Size factor",
      title = paste0("DEXSeq size factors",
                     if (nchar(grp_label) > 0) paste0(" (", grp_label, ")") else "")
    ) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))

  if (!is.null(outfile)) {
    w <- max(4, 0.5 * nrow(df) + 2)
    ggplot2::ggsave(outfile, g, width = w, height = 4, dpi = 300)
  }
  invisible(list(plot = g, data = df))
}

# ============================================================
#  FUNCTION: dispersion plot (base R via DEXSeq)
# ============================================================

make_dispersion_plot <- function(dxd, grp_label = "", outfile = NULL) {
  ttl <- paste0("Dispersion estimates",
                if (nchar(grp_label) > 0) paste0(" (", grp_label, ")") else "")
  if (!is.null(outfile)) png(outfile, width = 800, height = 600, res = 100)
  DESeq2::plotDispEsts(dxd, main = ttl)
  if (!is.null(outfile)) { dev.off(); invisible(NULL) }
}

# ============================================================
#  FUNCTION: p-value histogram
# ============================================================

make_pvalue_hist <- function(res_df, grp_label = "", outfile = NULL) {
  pv <- res_df$pvalue[!is.na(res_df$pvalue)]
  df <- data.frame(pvalue = pv)
  n_bins <- 20
  expected_height <- length(pv) / n_bins

  g <- ggplot2::ggplot(df, ggplot2::aes(x = pvalue)) +
    ggplot2::geom_histogram(bins = n_bins, fill = "steelblue", color = "white", boundary = 0) +
    ggplot2::geom_hline(yintercept = expected_height, linetype = "dashed",
                        color = "firebrick", linewidth = 0.8) +
    ggplot2::scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) +
    ggplot2::labs(
      x     = "p-value",
      y     = "Count",
      title = paste0("p-value distribution",
                     if (nchar(grp_label) > 0) paste0(" (", grp_label, ")") else ""),
      caption = "Dashed red line = uniform expectation"
    ) +
    ggplot2::theme_minimal(base_size = 11)

  if (!is.null(outfile)) ggplot2::ggsave(outfile, g, width = 5, height = 4, dpi = 300)
  invisible(g)
}

# ============================================================
#  FUNCTION: MA plot
# ============================================================

make_ma_plot <- function(res_df, grp_label = "",
                          padj_cut = PADJ_CUT, lfc_cut = LFC_CUT,
                          fc_col, outfile = NULL) {
  if (!all(c("exonBaseMean", fc_col, "padj") %in% names(res_df))) {
    warning("MA plot skipped — required columns missing.")
    return(invisible(NULL))
  }

  df <- res_df %>%
    dplyr::filter(!is.na(exonBaseMean), !is.na(.data[[fc_col]]), exonBaseMean > 0) %>%
    dplyr::mutate(
      logMean = log10(exonBaseMean),
      sig     = !is.na(padj) & padj <= padj_cut
    )

  n_up   <- sum(df$sig & df[[fc_col]] > 0, na.rm = TRUE)
  n_down <- sum(df$sig & df[[fc_col]] < 0, na.rm = TRUE)

  g <- ggplot2::ggplot(df, ggplot2::aes(x = logMean, y = .data[[fc_col]], color = sig)) +
    ggplot2::geom_point(size = 0.8, alpha = 0.6) +
    ggplot2::scale_color_manual(values = c("FALSE" = "gray60", "TRUE" = "firebrick"),
                                labels = c("Not sig", paste0("padj<", padj_cut)),
                                name = NULL) +
    ggplot2::geom_hline(yintercept = 0, linetype = "solid",  color = "black",   linewidth = 0.4) +
    ggplot2::geom_hline(yintercept = c(-lfc_cut, lfc_cut),
                        linetype = "dashed", color = "gray40", linewidth = 0.4) +
    ggplot2::annotate("text", x = -Inf, y = Inf,  hjust = -0.1, vjust = 1.5,
                      label = paste0("Up: ", n_up), color = "firebrick", size = 3.5) +
    ggplot2::annotate("text", x = -Inf, y = -Inf, hjust = -0.1, vjust = -0.5,
                      label = paste0("Down: ", n_down), color = "steelblue", size = 3.5) +
    ggplot2::labs(
      x     = expression(log[10](mean~PAS~count)),
      y     = expression(log[2]~fold~change),
      title = paste0("MA plot",
                     if (nchar(grp_label) > 0) paste0(" (", grp_label, ")") else "")
    ) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(legend.position = "top")

  if (!is.null(outfile)) ggplot2::ggsave(outfile, g, width = 6, height = 5, dpi = 300)
  invisible(g)
}

# ============================================================
#  FUNCTION: Volcano plot
# ============================================================

make_volcano_plot <- function(res_df, grp_label = "",
                               padj_cut = PADJ_CUT, lfc_cut = LFC_CUT,
                               fc_col, n_label = 10, outfile = NULL) {
  if (!all(c(fc_col, "padj", "gene") %in% names(res_df))) {
    warning("Volcano plot skipped — required columns missing.")
    return(invisible(NULL))
  }

  df <- res_df %>%
    dplyr::filter(!is.na(.data[[fc_col]]), !is.na(padj)) %>%
    dplyr::mutate(
      neg_log10_padj = -log10(pmax(padj, 1e-300)),
      category = dplyr::case_when(
        padj <= padj_cut & abs(.data[[fc_col]]) >= lfc_cut & .data[[fc_col]] > 0 ~ "Up (sig + |L2FC|)",
        padj <= padj_cut & abs(.data[[fc_col]]) >= lfc_cut & .data[[fc_col]] < 0 ~ "Down (sig + |L2FC|)",
        padj <= padj_cut ~ "Sig (small FC)",
        TRUE ~ "Not sig"
      )
    )

  top_labs <- df %>%
    dplyr::filter(padj <= padj_cut) %>%
    dplyr::slice_min(padj, n = n_label, with_ties = FALSE)

  cols <- c("Up (sig + |L2FC|)"   = "firebrick",
             "Down (sig + |L2FC|)" = "steelblue",
             "Sig (small FC)"      = "goldenrod3",
             "Not sig"             = "gray70")

  g <- ggplot2::ggplot(df, ggplot2::aes(x = .data[[fc_col]],
                                         y = neg_log10_padj,
                                         color = category)) +
    ggplot2::geom_point(size = 0.8, alpha = 0.6) +
    ggplot2::scale_color_manual(values = cols, name = NULL) +
    ggplot2::geom_vline(xintercept = c(-lfc_cut, lfc_cut),
                        linetype = "dashed", color = "gray40", linewidth = 0.4) +
    ggplot2::geom_hline(yintercept = -log10(padj_cut),
                        linetype = "dashed", color = "gray40", linewidth = 0.4) +
    ggrepel::geom_text_repel(data = top_labs,
                              ggplot2::aes(label = gene),
                              size = 2.5, max.overlaps = 20) +
    ggplot2::labs(
      x     = expression(log[2]~fold~change),
      y     = expression(-log[10](padj)),
      title = paste0("Volcano plot",
                     if (nchar(grp_label) > 0) paste0(" (", grp_label, ")") else "")
    ) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(legend.position = "top")

  if (!is.null(outfile)) ggplot2::ggsave(outfile, g, width = 6, height = 6, dpi = 300)
  invisible(g)
}

# ============================================================
#  FUNCTION: run DEXSeq for one group of samples
# ============================================================

run_dexseq_group <- function(grp_label, sub_samples, count_mat, colData,
                              featureID, groupID,
                              comparison_var, ref_label, alt_label,
                              gr = NULL, pas_anno = NULL,
                              min_total = 10, min_per_condition = 0,
                              use_ruv = FALSE, ruv_k = 1) {

  message(sprintf("[%s] Subsetting samples...", grp_label))
  sub_cd  <- droplevels(colData[sub_samples, , drop = FALSE])
  sub_cnt <- count_mat[, sub_samples, drop = FALSE]
  stopifnot(all(colnames(sub_cnt) == rownames(sub_cd)))

  sub_cd$sample <- factor(rownames(sub_cd))
  sub_cd[[comparison_var]] <- droplevels(factor(sub_cd[[comparison_var]], levels = c(ref_label, alt_label)))
  char_cols <- vapply(sub_cd, is.character, logical(1))
  sub_cd[char_cols] <- lapply(sub_cd[char_cols], factor)

  if (length(levels(sub_cd[[comparison_var]])) < 2 || any(table(sub_cd[[comparison_var]]) == 0)) {
    stop(sprintf("Group '%s' is missing one or both %s levels.", grp_label, comparison_var))
  }

  # W_k:exon corrects for batch effects on *relative* PAS usage (PSI); it is not
  # collinear with the `sample` term, which only absorbs gene-level expression totals.
  dex_design <- if (use_ruv && ruv_k >= 1) {
    w_terms <- paste(paste0("W_", seq_len(ruv_k)), "exon", sep = ":")
    as.formula(paste("~ sample + exon +", paste(c(w_terms, paste0(comparison_var, ":exon")), collapse = " + ")))
  } else {
    as.formula(paste0("~ sample + exon + ", comparison_var, ":exon"))
  }
  message(sprintf("[%s] Creating DEXSeqDataSet (design: %s)...",
                  grp_label, deparse(dex_design)))
  dxd <- DEXSeqDataSet(
    countData     = sub_cnt,
    sampleData    = sub_cd,
    design        = dex_design,
    featureID     = featureID,
    groupID       = groupID,
    featureRanges = gr
  )

  # Filter on "this" columns only — the doubled matrix inflates apparent counts
  message(sprintf("[%s] Filtering low-count PAS...", grp_label))
  this_idx <- which(colData(dxd)$exon == "this")
  mat_this <- counts(dxd)[, this_idx, drop = FALSE]
  dxd <- dxd[rowSums(mat_this) >= min_total, , drop = FALSE]

  if (min_per_condition > 0) {
    this_idx  <- which(colData(dxd)$exon == "this")
    mat_this  <- counts(dxd)[, this_idx, drop = FALSE]
    cond_this <- colData(dxd)[[comparison_var]][this_idx]
    keep2 <- rowSums(mat_this[, cond_this == ref_label, drop = FALSE] > 0) >= min_per_condition &
              rowSums(mat_this[, cond_this == alt_label, drop = FALSE] > 0) >= min_per_condition
    dxd <- dxd[keep2, , drop = FALSE]
  }

  # Drop genes with < 2 PAS remaining (no within-gene contrast)
  tab <- table(rowData(dxd)$groupID)
  dxd <- dxd[rowData(dxd)$groupID %in% names(tab[tab >= 2]), ]

  message(sprintf("[%s] Normalizing and testing (%d PAS, %d genes)...",
                  grp_label, nrow(dxd), length(unique(rowData(dxd)$groupID))))
  dxd <- estimateSizeFactors(dxd)
  dxd <- estimateDispersions(dxd)
  dxd <- testForDEU(dxd)
  dxr <- DEXSeqResults(dxd)
  dxd <- estimateExonFoldChanges(dxd, fitExpToVar = comparison_var,
                                  denominator = ref_label)

  message(sprintf("[%s] Building results table...", grp_label))
  res     <- S4Vectors::as.data.frame(dxr)
  rd      <- S4Vectors::mcols(dxd)
  fc_cols <- grep("^log2fold_", colnames(rd), value = TRUE)
  if (length(fc_cols) == 0) stop("No fold-change columns found in rowData(dxd).")

  fc_df <- base::as.data.frame(
    lapply(S4Vectors::as.list(rd[, fc_cols, drop = FALSE]), as.numeric),
    check.names = FALSE, stringsAsFactors = FALSE
  )
  key       <- rownames(dxd)
  res$key   <- rownames(res)
  fc_df$key <- key
  res_tp    <- merge(res, fc_df, by = "key", all.x = TRUE, sort = FALSE)
  rownames(res_tp) <- res_tp$key
  res_tp$key       <- NULL

  res_tp$gene    <- rowData(dxd)$groupID[match(rownames(res_tp), rownames(dxd))]
  res_tp$feature <- rowData(dxd)$featureID[match(rownames(res_tp), rownames(dxd))]

  if (!is.null(pas_anno)) {
    # Join on featureID_safe (sanitized); bring original PAS_ID through for readability
    ann_cols <- intersect(c("featureID_safe", "PAS_ID", "Intron_exon_location", "PAS_type"),
                          names(pas_anno))
    if (length(ann_cols) > 1) {
      ann <- pas_anno[, ann_cols, drop = FALSE]
      names(ann)[names(ann) == "featureID_safe"] <- "featureID"
      ann <- ann[!duplicated(ann$featureID), , drop = FALSE]
      ann$featureID    <- as.character(ann$featureID)
      res_tp$featureID <- as.character(res_tp$featureID)
      # merge() always resets rownames to "1","2",... regardless of `sort` --
      # carry the real per-PAS key through as an explicit column (same pattern
      # as the fc_df merge above) so every later step keyed on rownames(res_tp)
      # (the pos_meta join below, the normalized-counts merge further down)
      # still lines up with rownames(dxd) instead of matching against nothing.
      res_tp$key <- rownames(res_tp)
      res_tp     <- merge(res_tp, ann, by = "featureID", all.x = TRUE, sort = FALSE)
      rownames(res_tp) <- res_tp$key
      res_tp$key       <- NULL
    }
  }

  # APA_direction reports only the sign of THIS PAS's own usage change
  # (Increased/Decreased) -- it does not assert an overall UTR-lengthening/
  # shortening call, because with 3+ PAS per gene a shift onto a middle site
  # doesn't map cleanly onto "shorter" vs "longer" 3' UTR. `is_distal` below
  # is retained as positional context only; the actual lengthening/
  # shortening narrative is summarized properly across all PAS at once by
  # delta_wUTR in gene_summary.csv (see build_gene_apa_summary()).
  message(sprintf("[%s] Assigning usage direction...", grp_label))
  res_tp$APA_direction <- NA_character_

  if (!is.null(gr)) {
    rrn <- rownames(dxd)
    if (!is.null(names(gr)) && all(rrn %in% names(gr))) {
      pos      <- as.integer(GenomicRanges::start(gr[rrn]))
      strand   <- as.character(GenomicRanges::strand(gr[rrn]))
      pos_meta <- data.frame(row = rrn, pos = pos, strand = strand,
                             stringsAsFactors = FALSE)
      res_tp$row <- rownames(res_tp)
      res_tp     <- dplyr::left_join(res_tp, pos_meta, by = "row")

      res_tp <- res_tp |>
        dplyr::group_by(gene) |>
        dplyr::mutate(
          g_max     = suppressWarnings(max(pos, na.rm = TRUE)),
          g_min     = suppressWarnings(min(pos, na.rm = TRUE)),
          g_max     = ifelse(is.finite(g_max), g_max, NA_real_),
          g_min     = ifelse(is.finite(g_min), g_min, NA_real_),
          is_distal = dplyr::case_when(
            strand == "+" & !is.na(pos) & !is.na(g_max) & pos == g_max ~ TRUE,
            strand == "-" & !is.na(pos) & !is.na(g_min) & pos == g_min ~ TRUE,
            TRUE ~ FALSE
          )
        ) |>
        dplyr::ungroup()

      fc_col_dir <- grep("^log2fold_", names(res_tp), value = TRUE)[1]
      if (!is.na(fc_col_dir)) {
        res_tp <- dplyr::mutate(res_tp,
          APA_direction = dplyr::case_when(
            .data[[fc_col_dir]] > 0 ~ "Increased",
            .data[[fc_col_dir]] < 0 ~ "Decreased",
            TRUE ~ "Ambiguous"
          )
        )
        res_tp$APA_direction[is.na(res_tp[[fc_col_dir]])] <- "Ambiguous"
      }
      res_tp <- res_tp[, setdiff(names(res_tp), c("row", "g_max", "g_min")), drop = FALSE]
    } else {
      warning(sprintf("[%s] GRanges names don't cover all dxd features — usage direction skipped.", grp_label))
    }
  }

  # Normalized counts alongside the raw countData.* columns already in res_tp --
  # DEXSeq's own size-factor-corrected values (dxd has had estimateSizeFactors()
  # applied above), so a PAS's normalized and raw depth sit side by side per sample.
  message(sprintf("[%s] Adding normalized counts...", grp_label))
  norm_this_idx <- which(SummarizedExperiment::colData(dxd)$exon == "this")
  norm_mat      <- DESeq2::counts(dxd, normalized = TRUE)[, norm_this_idx, drop = FALSE]
  colnames(norm_mat) <- paste0("normCountData.",
                                as.character(SummarizedExperiment::colData(dxd)$sample[norm_this_idx]))
  norm_df     <- as.data.frame(norm_mat)
  norm_df$key <- rownames(norm_mat)
  # rownames(res_tp) is NOT trustworthy here -- the dplyr pipe above (like every
  # merge()/dplyr call in this function) silently resets data.frame rownames to
  # "1","2",... on its way through. Rebuild the true per-PAS key from `gene`/
  # `feature`, which are ordinary columns and always survive those calls intact,
  # using the same "groupID:featureID" format DEXSeq itself uses for rownames(dxd).
  res_tp$key <- paste(res_tp$gene, res_tp$feature, sep = ":")
  res_tp     <- merge(res_tp, norm_df, by = "key", all.x = TRUE, sort = FALSE)
  rownames(res_tp) <- res_tp$key
  res_tp$key       <- NULL

  fc_out_cols <- grep("^log2fold_", colnames(res_tp), value = TRUE)
  count_cols  <- grep("^countData\\.", colnames(res_tp), value = TRUE)
  # Group countData.<sample> columns by actual colData membership, not by matching
  # ref_label/alt_label as a substring of the sample name -- sample names aren't
  # guaranteed to be prefixed with the comparison label (e.g. the by_treatment axis's
  # ref_label "3day" is not a prefix of any sample name in a Control_3day_repN scheme).
  sample_of_col <- sub("^countData\\.", "", count_cols)
  grp_of_col    <- as.character(sub_cd[sample_of_col, comparison_var])
  ctrl_cols     <- count_cols[grp_of_col == ref_label]
  treat_cols    <- count_cols[grp_of_col == alt_label]

  # Interleave each raw count column with its normalized counterpart so the two
  # sit next to each other per sample, rather than in two separate blocks.
  interleave_counts <- function(raw_cols) {
    as.vector(rbind(raw_cols, paste0("normCountData.", sub("^countData\\.", "", raw_cols))))
  }

  desired_order <- c(
    "groupID", "featureID", "PAS_ID", "gene", "feature",
    "Intron_exon_location", "PAS_type", "APA_direction",
    "exonBaseMean", "dispersion", "stat", "pvalue", "padj",
    fc_out_cols, interleave_counts(ctrl_cols), interleave_counts(treat_cols),
    "genomicData.seqnames", "genomicData.start", "genomicData.end",
    "genomicData.width", "genomicData.strand"
  )
  desired_order <- desired_order[desired_order %in% colnames(res_tp)]
  res_tp <- res_tp[, desired_order, drop = FALSE]

  list(dxd = dxd, dxr = dxr, results = res_tp)
}

# ============================================================
#  FUNCTION: collapse per-PAS results to gene-level summary
# ============================================================

collapse_gene <- function(res_df, gene_q = NULL,
                           padj_cut = PADJ_CUT, lfc_cut = LFC_CUT,
                           fc_col) {
  safe_min  <- function(x) { x <- x[is.finite(x)]; if (length(x)) min(x)  else NA_real_ }
  safe_max  <- function(x) { x <- x[is.finite(x)]; if (length(x)) max(x)  else NA_real_ }
  safe_mean <- function(x) { x <- x[is.finite(x)]; if (length(x)) mean(x) else NA_real_ }
  has_dir   <- "APA_direction" %in% names(res_df)

  res_df %>%
    dplyr::mutate(
      padj_clean = ifelse(is.infinite(padj), NA_real_, padj),
      lfc_clean  = ifelse(is.infinite(.data[[fc_col]]), NA_real_, .data[[fc_col]])
    ) %>%
    dplyr::mutate(
      sig = !is.na(padj_clean) & padj_clean <= padj_cut,
      big = sig & abs(lfc_clean) >= lfc_cut
    ) %>%
    dplyr::group_by(groupID) %>%
    dplyr::summarise(
      n_PAS        = dplyr::n(),
      n_sig        = sum(sig, na.rm = TRUE),
      n_big        = sum(big, na.rm = TRUE),
      min_padj     = safe_min(padj_clean),
      max_absL2FC  = safe_max(abs(lfc_clean[sig])),
      mean_absL2FC = safe_mean(abs(lfc_clean[sig])),
      dir_consensus = if (has_dir) {
        d <- APA_direction[sig & !is.na(APA_direction)]
        if (!length(d))                 NA_character_
        else if (all(d == "Increased")) "Increased_only"
        else if (all(d == "Decreased")) "Decreased_only"
        else                            "Mixed"
      } else NA_character_,
      .groups = "drop"
    ) %>%
    dplyr::mutate(perGeneQ = if (!is.null(gene_q)) unname(gene_q[groupID]) else NA_real_) %>%
    dplyr::arrange(is.na(perGeneQ), perGeneQ, min_padj, dplyr::desc(max_absL2FC))
}

# ============================================================
#  FUNCTION: per-condition mean PSI + SEM + RED
# ============================================================

usage_by_condition <- function(res_df, sub_cd,
                                comparison_var, ref_label, alt_label,
                                usage_ctrl_col, usage_trtmt_col,
                                se_ctrl_col, se_trtmt_col,
                                red_ctrl_col, red_trtmt_col,
                                thin_cut = THIN_COUNTS_CUT) {
  cnt_cols <- grep("^countData\\.", names(res_df), value = TRUE)
  counts   <- as.matrix(res_df[, cnt_cols, drop = FALSE])
  colnames(counts) <- sub("^countData\\.", "", colnames(counts))
  stopifnot(all(colnames(counts) %in% rownames(sub_cd)))

  gene_totals <- rowsum(counts, res_df$groupID)
  totals_row  <- gene_totals[res_df$groupID, colnames(counts), drop = FALSE]
  # A sample with zero reads across all of a gene's retained PAS has undefined
  # usage there, not 0% usage -- mark it NA so it is excluded (not averaged in
  # as a real observation) below, instead of silently dragging the mean down.
  usage <- counts / totals_row
  usage[totals_row == 0] <- NA_real_

  cond     <- sub_cd[colnames(counts), comparison_var, drop = TRUE]
  ctrl_idx <- which(cond == ref_label)
  trt_idx  <- which(cond == alt_label)

  se_fn <- function(x) {
    x <- x[!is.na(x)]
    if (length(x) > 1) sd(x) / sqrt(length(x)) else NA_real_
  }

  out <- res_df
  out[["meanUsage_All"]] <- rowMeans(usage, na.rm = TRUE)
  out[[usage_ctrl_col]]  <- if (length(ctrl_idx)) rowMeans(usage[, ctrl_idx, drop = FALSE], na.rm = TRUE) else NA_real_
  out[[usage_trtmt_col]] <- if (length(trt_idx))  rowMeans(usage[, trt_idx,  drop = FALSE], na.rm = TRUE) else NA_real_
  out[[se_ctrl_col]]     <- if (length(ctrl_idx) > 1) apply(usage[, ctrl_idx, drop = FALSE], 1, se_fn) else rep(NA_real_, nrow(usage))
  out[[se_trtmt_col]]    <- if (length(trt_idx)  > 1) apply(usage[, trt_idx,  drop = FALSE], 1, se_fn) else rep(NA_real_, nrow(usage))

  # RED (Relative Expression Difference), our N-site-generalized version of the
  # MAAPER RED metric -- see build_gene_apa_summary()'s max_abs_delta_RED for the
  # gene-level rollup. Computed from summed *normalized* counts per condition, not
  # the meanUsage_* fractions above, because the pseudocount correction below only
  # makes sense added to an absolute count scale, not a 0-1 fraction.
  norm_cnt_cols <- grep("^normCountData\\.", names(res_df), value = TRUE)
  stopifnot(length(norm_cnt_cols) > 0)
  norm_counts <- as.matrix(res_df[, norm_cnt_cols, drop = FALSE])
  colnames(norm_counts) <- sub("^normCountData\\.", "", colnames(norm_counts))
  norm_counts <- norm_counts[, colnames(counts), drop = FALSE]  # align to counts'/ctrl_idx's/trt_idx's sample order

  n_ctrl <- if (length(ctrl_idx)) rowSums(norm_counts[, ctrl_idx, drop = FALSE], na.rm = TRUE) else rep(NA_real_, nrow(norm_counts))
  n_trt  <- if (length(trt_idx))  rowSums(norm_counts[, trt_idx,  drop = FALSE], na.rm = TRUE) else rep(NA_real_, nrow(norm_counts))

  gene_norm_totals <- rowsum(norm_counts, res_df$groupID)
  gene_tot_ctrl    <- rowSums(gene_norm_totals[, ctrl_idx, drop = FALSE], na.rm = TRUE)[res_df$groupID]
  gene_tot_trt     <- rowSums(gene_norm_totals[, trt_idx,  drop = FALSE], na.rm = TRUE)[res_df$groupID]

  m_ctrl         <- gene_tot_ctrl - n_ctrl
  m_trt          <- gene_tot_trt  - n_trt
  n_pas_per_gene <- as.integer(table(res_df$groupID)[res_df$groupID])

  red_fn <- function(n, m, N) log2(N * (n + RED_PSEUDOCOUNT) / (m + RED_PSEUDOCOUNT))
  out[[red_ctrl_col]]  <- red_fn(n_ctrl, m_ctrl, n_pas_per_gene)
  out[[red_trtmt_col]] <- red_fn(n_trt,  m_trt,  n_pas_per_gene)
  out[["delta_RED"]]   <- out[[red_trtmt_col]] - out[[red_ctrl_col]]

  # thin_counts is a soft warning, not a filter: MIN_TOTAL_READS already excludes
  # PAS below its (gene-total) bar before DEXSeq ever sees them, so this flags
  # PAS that cleared that hard cutoff but still have < thin_cut raw reads
  # summed across replicates in either condition on their own.
  ctrl_total <- if (length(ctrl_idx)) rowSums(counts[, ctrl_idx, drop = FALSE], na.rm = TRUE) else rep(0, nrow(counts))
  trt_total  <- if (length(trt_idx))  rowSums(counts[, trt_idx,  drop = FALSE], na.rm = TRUE) else rep(0, nrow(counts))
  out[["thin_counts"]] <- ctrl_total < thin_cut | trt_total < thin_cut
  out
}

# ============================================================
#  FUNCTION: gene-level differential expression (real DESeq2 test,
#            not derivable from DEXSeq's per-PAS exon-usage results)
# ============================================================
#
# Answers a different question than everything else in this script: is the
# gene as a whole up or down, independent of how reads are distributed across
# its PAS. Sums each gene's retained-PAS raw counts (the same "this"-only
# counts already used everywhere else in this pipeline) to gene level, then
# runs a standard single-factor DESeq2 test -- ~condition, or ~W_1+...+
# condition when USE_RUV is on, matching how batch correction is already
# handled for the exon-usage test. Reuses the already-built, already-filtered
# dxd rather than rebuilding gene totals from scratch, so this stays exactly
# consistent with what's tested elsewhere for the same group.

build_gene_dge <- function(dxd, comparison_var, ctrl_label, trt_label) {
  cd_full  <- as.data.frame(SummarizedExperiment::colData(dxd))
  this_idx <- which(cd_full$exon == "this")
  cd       <- cd_full[this_idx, , drop = FALSE]
  mat      <- SummarizedExperiment::assay(dxd)[, this_idx, drop = FALSE]
  colnames(mat) <- as.character(cd$sample)
  rownames(cd)  <- as.character(cd$sample)

  gene_counts <- rowsum(mat, S4Vectors::mcols(dxd)$groupID)
  storage.mode(gene_counts) <- "integer"

  w_cols <- grep("^W_", names(cd), value = TRUE)
  design <- if (length(w_cols) > 0) {
    as.formula(paste("~", paste(c(w_cols, comparison_var), collapse = " + ")))
  } else {
    as.formula(paste0("~ ", comparison_var))
  }

  dds <- DESeq2::DESeqDataSetFromMatrix(countData = gene_counts, colData = cd, design = design)
  dds <- DESeq2::DESeq(dds, quiet = TRUE)
  res <- DESeq2::results(dds, contrast = c(comparison_var, trt_label, ctrl_label))

  # Mean size-factor-normalized gene-level count per condition, reusing this same
  # DESeqDataSet (no extra model fit) -- feeds the gene-count bar panel in
  # plot_gene_apa_zoom(), placed next to the log2FC bar it's derived from.
  norm_counts    <- DESeq2::counts(dds, normalized = TRUE)
  cond           <- as.character(cd[[comparison_var]])
  ctrl_idx       <- which(cond == ctrl_label)
  trt_idx        <- which(cond == trt_label)
  norm_ctrl_mean <- if (length(ctrl_idx)) rowMeans(norm_counts[, ctrl_idx, drop = FALSE]) else rep(NA_real_, nrow(norm_counts))
  norm_trt_mean  <- if (length(trt_idx))  rowMeans(norm_counts[, trt_idx,  drop = FALSE]) else rep(NA_real_, nrow(norm_counts))

  # Per-sample normalized counts too (not just the per-condition mean) --
  # feeds the individual sample dots + error bars on the gene-count bar panel.
  norm_counts_df <- as.data.frame(norm_counts)
  names(norm_counts_df) <- paste0("normCount.", names(norm_counts_df))
  norm_counts_df$groupID <- rownames(norm_counts_df)

  normcount_ctrl_col <- paste0("normCount_", ctrl_label)
  normcount_trtmt_col <- paste0("normCount_", trt_label)

  out <- as.data.frame(res)
  out$groupID <- rownames(out)
  out[[normcount_ctrl_col]]  <- norm_ctrl_mean[out$groupID]
  out[[normcount_trtmt_col]] <- norm_trt_mean[out$groupID]
  out <- merge(out, norm_counts_df, by = "groupID", all.x = TRUE, sort = FALSE)
  rownames(out) <- NULL
  out
}

# ============================================================
#  FUNCTION: gene-level APA summary — weighted UTR-length shift,
#            dominant-site tracking, RED rollup, chi-squared cross-check
# ============================================================
#
# delta_wUTR ("weighted UTR length") replaces the old binary
# Lengthened/Shortened call with a continuous metric that generalizes to any
# number of PAS per gene: for each condition, wUTR = sum(usage_i * distance_i),
# where distance_i is each PAS's genomic distance from the gene's own
# most-proximal surviving PAS. (No stop-codon coordinate is available in the
# PolyA_DB annotation, so the proximal-most PAS is used as an internally
# consistent reference point instead.) A positive delta_wUTR means usage
# shifted toward more distal sites on average (longer 3' UTR); negative means
# a shift toward proximal sites -- and unlike APA_direction, this is not
# limited to a distal-vs-everything-else binary.
#
# single_dominant_site marks genes that functionally have one PAS: the most-
# abundant PAS carries >= dominant_usage_cut (default 90%) of the gene's
# usage in BOTH conditions, and that share doesn't shift between them. For
# these genes, whatever the remaining (minor) PAS are doing individually is
# unlikely to be biologically meaningful -- one site is running the show
# either way. It is a flag, not a filter -- rows are never dropped, so
# excluding them from downstream use is a deliberate choice made on read,
# not silent data loss.
#
# chisq_stat/chisq_pvalue is an independent cross-check computed directly on
# summed raw counts (a PAS x condition contingency table per gene) -- it does
# not depend on estimateSizeFactors/DEXSeq's NB model, so it does not inherit
# any bias from the usage-averaging step above. Expect it to run hotter
# (more sensitive, no overdispersion correction) than the DEXSeq padj values.
#
# max_abs_delta_RED/max_delta_RED_featureID roll the per-PAS RED metric (see
# usage_by_condition()) up to one number per gene -- the single biggest
# relative-expression shift among that gene's PAS -- so genes can be ranked
# by RED the way they used to be under MAAPER, but correctly accounting for
# every tested PAS instead of just a proximal/distal pair.

build_gene_apa_summary <- function(res_u, group_of_sample,
                                    ctrl_col, trt_col,
                                    ref_label, alt_label,
                                    wutr_ctrl_col, wutr_trtmt_col,
                                    dom_usage_ctrl_col, dom_usage_trtmt_col,
                                    dominant_stable_cut = DOMINANT_STABLE_CUT,
                                    dominant_usage_cut  = DOMINANT_USAGE_CUT) {
  cnt_cols <- grep("^countData\\.", names(res_u), value = TRUE)
  # Group countData.<sample> columns by actual colData membership (group_of_sample),
  # not by matching ref_label/alt_label as a substring of the sample name -- see the
  # identical fix in run_dexseq_group() for why sample-name matching isn't safe here.
  sample_of_col <- sub("^countData\\.", "", cnt_cols)
  grp_of_col    <- unname(group_of_sample[sample_of_col])
  ctrl_cnt      <- cnt_cols[grp_of_col == ref_label]
  trt_cnt       <- cnt_cols[grp_of_col == alt_label]

  df <- res_u
  df$pos        <- dplyr::coalesce(df$genomicData.start, NA_integer_)
  df$strand     <- as.character(dplyr::coalesce(df$genomicData.strand, NA))
  df$ctrl_total <- rowSums(as.matrix(res_u[, ctrl_cnt, drop = FALSE]), na.rm = TRUE)
  df$trt_total  <- rowSums(as.matrix(res_u[, trt_cnt,  drop = FALSE]), na.rm = TRUE)

  df <- df %>%
    dplyr::group_by(groupID) %>%
    dplyr::mutate(
      g_max    = suppressWarnings(max(pos, na.rm = TRUE)),
      g_min    = suppressWarnings(min(pos, na.rm = TRUE)),
      proximal = if (all(na.omit(strand) == "-")) g_max else g_min,
      distance = abs(pos - proximal)
    ) %>%
    dplyr::ungroup()

  chisq_for_gene <- function(sub) {
    tab <- rbind(sub$ctrl_total, sub$trt_total)
    tab <- tab[, colSums(tab) > 0, drop = FALSE]
    if (ncol(tab) < 2 || any(rowSums(tab) == 0)) return(c(stat = NA_real_, pvalue = NA_real_))
    res <- tryCatch(suppressWarnings(chisq.test(tab)), error = function(e) NULL)
    if (is.null(res)) return(c(stat = NA_real_, pvalue = NA_real_))
    c(stat = unname(res$statistic), pvalue = unname(res$p.value))
  }

  out <- df %>%
    dplyr::group_by(groupID) %>%
    dplyr::group_modify(function(sub, key) {
      n_pas <- nrow(sub)
      dom_i <- which.max(sub$meanUsage_All)
      if (length(dom_i) == 0) dom_i <- 1L

      wutr_ctrl       <- sum(sub[[ctrl_col]] * sub$distance, na.rm = TRUE)
      wutr_trt        <- sum(sub[[trt_col]]  * sub$distance, na.rm = TRUE)
      dom_usage_ctrl  <- sub[[ctrl_col]][dom_i]
      dom_usage_trt   <- sub[[trt_col]][dom_i]
      dom_delta_usage <- dom_usage_trt - dom_usage_ctrl
      dom_stable      <- is.na(dom_delta_usage) || abs(dom_delta_usage) < dominant_stable_cut
      dom_dominant    <- !is.na(dom_usage_ctrl) & !is.na(dom_usage_trt) &
                          dom_usage_ctrl >= dominant_usage_cut & dom_usage_trt >= dominant_usage_cut
      cs <- chisq_for_gene(sub)

      # Gene-level RED rollup: the PAS with the single largest |delta_RED| --
      # the direct replacement for "rank genes by their biggest RED shift"
      # (see usage_by_condition() for the per-PAS RED_Control/RED_Treatment/
      # delta_RED this is built from). Independent of dominant_featureID above:
      # the PAS driving the biggest relative-expression shift need not be the
      # gene's most-abundant PAS.
      red_i   <- which.max(abs(sub$delta_RED))
      has_red <- length(red_i) > 0

      data.frame(
        n_PAS                 = n_pas,
        dominant_featureID    = sub$featureID[dom_i],
        dominant_abundance    = sub$exonBaseMean[dom_i],
        dom_usage_ctrl_ph     = dom_usage_ctrl,
        dom_usage_trt_ph      = dom_usage_trt,
        dominant_delta_usage  = dom_delta_usage,
        dominant_padj         = sub$padj[dom_i],
        wUTR_Control_ph       = wutr_ctrl,
        wUTR_Treatment_ph     = wutr_trt,
        delta_wUTR            = wutr_trt - wutr_ctrl,
        chisq_stat            = cs["stat"],
        chisq_pvalue          = cs["pvalue"],
        single_dominant_site  = isTRUE(dom_dominant) & dom_stable,
        max_abs_delta_RED       = if (has_red) abs(sub$delta_RED[red_i]) else NA_real_,
        max_delta_RED_featureID = if (has_red) sub$featureID[red_i]      else NA_character_,
        stringsAsFactors      = FALSE
      )
    }) %>%
    dplyr::ungroup() %>%
    as.data.frame()

  out$chisq_padj <- p.adjust(out$chisq_pvalue, method = "BH")

  names(out)[names(out) == "wUTR_Control_ph"]   <- wutr_ctrl_col
  names(out)[names(out) == "wUTR_Treatment_ph"] <- wutr_trtmt_col
  names(out)[names(out) == "dom_usage_ctrl_ph"] <- dom_usage_ctrl_col
  names(out)[names(out) == "dom_usage_trt_ph"]  <- dom_usage_trtmt_col
  out
}

# ============================================================
#  FUNCTION: candidate-gene rubric
# ============================================================
#
# Our replacement for the old MAAPER/RED-based candidate-selection rubric
# (Kathleen's process): for each gene, take its single biggest-moving PAS (by
# |usage change|, same "greatest change in percent usage of any PAS" idea the
# old rubric used) and require that PAS/gene to clear every one of:
#   - usage change      > USAGE_CHANGE_CUT   (meaningful absolute shift, not
#                                              just a big relative one)
#   - padj              < padj_cut           (DEXSeq's own significance call
#                                              on that PAS, replacing the old
#                                              separate MAAPER p-value)
#   - not thin_counts                        (adequate read depth, replacing
#                                              the old ">40 reads" bar)
#   - |delta_wUTR|      > WUTR_CHANGE_CUT     (a real gene-wide 3' UTR length
#                                              shift, not just one site moving)
#   - not single_dominant_site                (excludes genes where a shift is
#                                              confined to a site that carries
#                                              minor content either way)
#   - |delta_RED|       >= RED_MAG_CUT        (our generalized analog of the
#                                              old |REDu1| >= 1 magnitude bar)
# This is a filter, not a flag -- genes that don't clear all six simply don't
# appear in the output, which is the point of this table (a curated
# candidate list, not a full annotated universe like gene_summary.csv).

build_candidate_genes <- function(res_u, gene_summary,
                                   ctrl_col, trt_col,
                                   usage_change_cut = USAGE_CHANGE_CUT,
                                   padj_cut          = PADJ_CUT,
                                   wutr_change_cut   = WUTR_CHANGE_CUT,
                                   red_mag_cut       = RED_MAG_CUT) {
  top_per_gene <- res_u %>%
    dplyr::mutate(meanUsage_change = abs(.data[[trt_col]] - .data[[ctrl_col]])) %>%
    dplyr::filter(!is.na(meanUsage_change)) %>%
    dplyr::group_by(groupID) %>%
    dplyr::slice_max(meanUsage_change, n = 1, with_ties = FALSE) %>%
    dplyr::ungroup()

  top_per_gene %>%
    dplyr::left_join(gene_summary, by = "groupID", suffix = c("", "_gene")) %>%
    dplyr::filter(
      meanUsage_change > usage_change_cut,
      !is.na(padj), padj < padj_cut,
      !thin_counts,
      !is.na(delta_wUTR), abs(delta_wUTR) > wutr_change_cut,
      !dplyr::coalesce(single_dominant_site, FALSE),
      !is.na(delta_RED), abs(delta_RED) >= red_mag_cut
    ) %>%
    dplyr::arrange(padj) %>%
    dplyr::select(groupID, gene, featureID, meanUsage_change, padj, delta_wUTR,
                  delta_RED, max_abs_delta_RED, max_delta_RED_featureID,
                  perGeneQ, dir_consensus, dominant_featureID) %>%
    as.data.frame()
}

# ============================================================
#  FUNCTION: whole-gene structure map (real coordinates, all
#            annotated isoforms, highlighting AS/APA differences)
# ============================================================
#
# Requires gtf_exons (see "LOAD TRANSCRIPT ANNOTATION" above) -- real exon/
# intron structure from the RefSeq GTF used to build this project's STAR
# index. Purpose: show how many isoforms are annotated and where they
# structurally differ (alternative splicing, alternative terminal exons) --
# usage/expression per PAS live in plot_gene_apa_zoom() instead, since
# cramming both concerns into one figure made neither easy to read.
#
# Only an isoform's TERMINAL exon can be colored Alternative -- this project
# cares about AS+APA interaction (a different terminal exon means a
# different splice choice AND a different poly(A) site at once), not
# alternative splicing on its own. Internal cassette exons/alternative splice
# sites elsewhere in the gene body are real AS but unrelated to which PAS a
# transcript ends at, so they're drawn as plain (grey) regardless of how much
# they vary between isoforms -- highlighting them would just be noise for
# this question. Constitutive/Alternative among terminal exons still means
# the same thing as before: identical start/end shared by every OTHER
# isoform that overlaps this genomic region, isoforms compared only against
# others whose span actually reaches that region.
#
# ggforce::facet_zoom() draws a second, cropped panel below the whole-gene
# view (zoomed to the terminal exon(s) near the detected PAS) with an
# automatic connector between the two -- tightly clustered PAS are illegible
# at whole-gene scale (e.g. 5 PAS within ~1kb of a ~200kb gene), so this
# isn't drawn to the same scale as the panel above it. An isoform with no
# data in that cropped region (e.g. a transcript from a distant alternative
# promoter) simply shows as an empty row there, which is accurate -- it
# isn't present in that region, not excluded by any filtering step.
#
# NOTE: facet_zoom's connector is drawn via custom "zoom.x"/"zoom.y" theme
# elements it registers; theme_minimal() (or any other *complete* theme) is
# a full replacement that doesn't carry those over, which silently makes the
# connector disappear. Re-declaring them explicitly in theme() below is what
# keeps it working.

plot_gene_apa_genome <- function(res_u, gene_symbol, gtf_exons, title_suffix = NULL,
                                  near_bp = 5000, pad_frac = 0.15, zoom_size = 1) {
  tx_exons <- gtf_exons[gtf_exons$gene_name == gene_symbol, ]
  if (nrow(tx_exons) == 0) stop(sprintf("Gene '%s' not found in GTF (gene_name mismatch?).", gene_symbol))
  strand_sym <- tx_exons$strand[1]

  pas_df <- res_u %>%
    dplyr::filter(gene == gene_symbol) %>%
    dplyr::mutate(pos = dplyr::coalesce(.data$genomicData.start, NA_integer_)) %>%
    dplyr::filter(!is.na(pos)) %>%
    dplyr::arrange(pos)
  if (nrow(pas_df) == 0) stop(sprintf("No positioned PAS for gene '%s'.", gene_symbol))
  pas_range <- range(pas_df$pos)

  tx_span <- tx_exons %>%
    dplyr::group_by(transcript_id) %>%
    dplyr::summarise(tx_start = min(start), tx_end = max(end), .groups = "drop")
  tx_order <- tx_span$transcript_id[order(tx_span$tx_end - tx_span$tx_start)]

  tx_exons <- tx_exons %>%
    dplyr::group_by(transcript_id) %>%
    dplyr::mutate(is_terminal_exon = if (strand_sym == "-") start == min(start) else end == max(end)) %>%
    dplyr::ungroup()

  key_summary <- tx_exons %>%
    dplyr::distinct(start, end, transcript_id) %>%
    dplyr::group_by(start, end) %>%
    dplyr::summarise(having = list(unique(transcript_id)), .groups = "drop") %>%
    dplyr::rowwise() %>%
    dplyr::mutate(
      n_relevant = sum(tx_span$tx_start <= end & tx_span$tx_end >= start),
      n_having   = length(having),
      exon_class = ifelse(n_having == n_relevant, "Constitutive", "Alternative")
    ) %>%
    dplyr::ungroup() %>%
    dplyr::select(start, end, exon_class)
  tx_exons <- tx_exons %>% dplyr::left_join(key_summary, by = c("start", "end"))
  tx_exons$exon_class <- ifelse(tx_exons$is_terminal_exon, tx_exons$exon_class, "Constitutive")
  n_alt <- sum(tx_exons$exon_class == "Alternative")
  fill_levels <- if (n_alt > 0) c("Constitutive", "Alternative") else "Constitutive"
  tx_exons$exon_class <- factor(tx_exons$exon_class, levels = fill_levels)
  tx_exons$y <- factor(tx_exons$transcript_id, levels = tx_order)
  introns    <- ggtranscript::to_intron(tx_exons, "transcript_id")
  introns$y  <- factor(introns$transcript_id, levels = tx_order)

  # terminal exon extended into intron, wherever a PAS lies beyond it (same as plot_gene_apa_zoom)
  term_exon <- tx_exons %>%
    dplyr::group_by(transcript_id) %>%
    dplyr::slice(if (strand_sym == "-") which.min(start) else which.max(end)) %>%
    dplyr::ungroup() %>%
    dplyr::filter(start <= pas_range[2] + near_bp, end >= pas_range[1] - near_bp)
  ext <- term_exon %>%
    dplyr::rowwise() %>%
    dplyr::mutate(
      ext_start = if (strand_sym == "-") min(pas_df$pos[pas_df$pos < start], start) else end,
      ext_end   = if (strand_sym == "-") start else max(pas_df$pos[pas_df$pos > end], end)
    ) %>% dplyr::ungroup() %>%
    dplyr::filter(ext_end > ext_start) %>%
    dplyr::transmute(transcript_id, start = ext_start, end = ext_end)
  if (nrow(ext) > 0) ext$y <- factor(ext$transcript_id, levels = tx_order)

  win      <- range(c(term_exon$start, term_exon$end, ext$start, ext$end, pas_df$pos))
  zoom_pad <- diff(win) * pad_frac + 300
  zoom_xr  <- c(win[1] - zoom_pad, win[2] + zoom_pad)

  ttl <- paste0(gene_symbol, if (!is.null(title_suffix)) paste0(" — ", title_suffix) else "")
  fill_scale <- ggplot2::scale_fill_manual(values = c(Constitutive = "grey35", Alternative = CAT_TRT), name = NULL)

  g <- ggplot2::ggplot() +
    ggplot2::geom_vline(data = pas_df, ggplot2::aes(xintercept = pos),
                         linetype = "dotted", color = "grey30", linewidth = 0.75) +
    ggtranscript::geom_range(data = tx_exons, ggplot2::aes(xstart = start, xend = end, y = y, fill = exon_class),
                              height = 0.4) +
    ggtranscript::geom_intron(data = introns, ggplot2::aes(xstart = start, xend = end, y = y),
                               strand = strand_sym, arrow.min.intron.length = 200, color = "grey50")
  if (nrow(ext) > 0)
    g <- g + ggtranscript::geom_range(data = ext, ggplot2::aes(xstart = start, xend = end, y = y),
                                       fill = NA, color = "grey35", linetype = "dashed", height = 0.4)

  g +
    fill_scale +
    ggplot2::scale_x_continuous(labels = scales::comma) +
    ggplot2::labs(x = paste0(tx_exons$seqnames[1], " (", strand_sym, " strand), ",
                             length(unique(pas_df$pos)), " PAS"), y = NULL, title = ttl) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(panel.grid.minor = ggplot2::element_blank(), legend.position = "top",
                   zoom.x = ggplot2::element_rect(fill = "grey88", colour = NA),
                   zoom.y = ggplot2::element_rect(fill = "grey88", colour = NA)) +
    ggforce::facet_zoom(xlim = zoom_xr, horizontal = FALSE, zoom.size = zoom_size, show.area = TRUE)
}

# ============================================================
#  FUNCTION: terminal-exon zoom (real coordinates, expanded into
#            intron wherever a detected PAS lies beyond the
#            annotated exon boundary)
# ============================================================
#
# Three panels stacked on a shared genomic x-axis (patchwork): gene-level DGE
# summary on top (is the gene as a whole up or down -- see build_gene_dge()),
# PAS usage (0-100%) in the middle, gene model on the bottom. An earlier
# version plotted per-PAS normalized counts in the top panel instead, but
# that's misleading on its own: if a gene's total expression differs between
# conditions, every one of its PAS shifts together in raw count terms even
# when relative usage hasn't changed at all, so a consistent per-PAS offset
# there was often just restating one gene-level fact N times rather than
# adding information. The one-line gene DGE summary answers that question
# directly instead. Usage shows every replicate as a jittered point (not
# just the mean) with +-1 SEM error bars, so you can see whether a crossing
# pattern between conditions is a robust shift or within replicate noise.
#
# Isoforms whose terminal exon is nowhere near the PAS cluster (e.g. a
# retained-intron transcript at a different locus) are excluded from the
# gene-model panel here, not because they're wrong, but because they're not
# part of this APA story. Where a PAS falls beyond a transcript's annotated
# terminal-exon boundary, a dashed-outline extension is drawn out to that
# PAS -- exon "extended into intron" by the detected site.

plot_gene_apa_zoom <- function(res_u, gene_symbol, gtf_exons, sub_cd, gene_dge, title_suffix = NULL,
                                ctrl_col, trt_col, se_ctrl_col, se_trt_col,
                                comparison_var, ref_label, alt_label, comparison_label,
                                near_bp = 5000, pad_frac = 0.15) {
  tx_exons <- gtf_exons[gtf_exons$gene_name == gene_symbol, ]
  if (nrow(tx_exons) == 0) stop(sprintf("Gene '%s' not found in GTF.", gene_symbol))
  strand_sym <- tx_exons$strand[1]
  # Minus-strand genes are transcribed right-to-left in real genomic coordinates;
  # reverse the x-axis so every gene reads 5'->3' left to right regardless of
  # strand. This is a display-only flip (scale_x_reverse) -- none of the
  # position/window/extension math above or below needs to change, since it all
  # still operates on real ascending genomic coordinates.
  flip_x <- strand_sym == "-"

  pas_df <- res_u %>%
    dplyr::filter(gene == gene_symbol) %>%
    dplyr::mutate(pos = dplyr::coalesce(.data$genomicData.start, NA_integer_)) %>%
    dplyr::filter(!is.na(pos)) %>%
    dplyr::arrange(pos)
  if (nrow(pas_df) == 0) stop(sprintf("No positioned PAS for gene '%s'.", gene_symbol))
  pas_range <- range(pas_df$pos)

  term_exon <- tx_exons %>%
    dplyr::group_by(transcript_id) %>%
    dplyr::slice(if (strand_sym == "-") which.min(start) else which.max(end)) %>%
    dplyr::ungroup() %>%
    dplyr::filter(start <= pas_range[2] + near_bp, end >= pas_range[1] - near_bp)
  if (nrow(term_exon) == 0)
    stop(sprintf("No terminal exon within %d bp of gene '%s's PAS sites.", near_bp, gene_symbol))
  keep_tx <- term_exon$transcript_id

  ext <- term_exon %>%
    dplyr::rowwise() %>%
    dplyr::mutate(
      ext_start = if (strand_sym == "-") min(pas_df$pos[pas_df$pos < start], start) else end,
      ext_end   = if (strand_sym == "-") start else max(pas_df$pos[pas_df$pos > end], end)
    ) %>%
    dplyr::ungroup() %>%
    dplyr::filter(ext_end > ext_start) %>%
    dplyr::transmute(transcript_id, start = ext_start, end = ext_end)

  tx_order <- term_exon$transcript_id[order(term_exon$end - term_exon$start)]

  win <- range(c(term_exon$start, term_exon$end, ext$start, ext$end, pas_df$pos))
  pad <- diff(win) * pad_frac + 300
  xr  <- c(win[1] - pad, win[2] + pad)

  # Keep every exon of each retained transcript (not just those within xr) so that
  # to_intron() below can still compute the intron leading into the terminal exon --
  # the view is cropped to xr via coord_cartesian(), not by discarding data, so a
  # partially-visible upstream intron (and its direction arrow) still renders, the
  # same way ggforce::facet_zoom() crops plot_gene_apa_genome()'s zoomed inset.
  ctx_exons <- tx_exons %>%
    dplyr::filter(transcript_id %in% keep_tx)

  cond_colors <- c(setNames(CAT_CTRL, ref_label), setNames(CAT_TRT, alt_label))
  vlines <- ggplot2::geom_vline(xintercept = pas_df$pos, linetype = "dotted", color = "grey30", linewidth = 0.75)

  # --- panel 1: gene-level DGE summary ---
  # No error bar here (the lfcSE addition was more confusing than informative for a
  # single-number summary). Bar is vertical -- up for a positive log2FC, down for a
  # negative one -- with the log2FC/padj label anchored beside the bar at a fixed
  # y (0, the zero line) so its position never depends on the bar's sign or length.
  dge_row <- gene_dge[gene_dge$groupID == gene_symbol, ]
  lfc    <- if (nrow(dge_row)) dge_row$log2FoldChange[1] else NA_real_
  padj_g <- if (nrow(dge_row)) dge_row$padj[1] else NA_real_
  direction <- if (is.na(lfc) || is.na(padj_g) || padj_g >= PADJ_CUT) "NS" else if (lfc > 0) "Up" else "Down"
  lfc0 <- ifelse(is.na(lfc), 0, lfc)
  dge_buf  <- max(abs(lfc0) * 0.3, 0.5)  # always reserve a minimum gap above/below the bar
  dge_ylim <- c(min(0, lfc0) - dge_buf, max(0, lfc0) + dge_buf)
  dge_df   <- data.frame(x = 1, lfc = lfc0, direction = direction)
  dge_label <- sprintf("log2FC=%s\npadj=%s",
                        ifelse(is.na(lfc), "NA", sprintf("%.2f", lfc)),
                        ifelse(is.na(padj_g), "NA", formatC(padj_g, digits = 2, format = "g")))

  p_dge <- ggplot2::ggplot(dge_df, ggplot2::aes(x = x, y = lfc, fill = direction)) +
    ggplot2::geom_col(width = 0.5) +
    ggplot2::geom_hline(yintercept = 0, linewidth = 0.4) +
    ggplot2::annotate("text", x = 1.55, y = 0, label = dge_label,
                       hjust = 0, vjust = 0.5, size = 3, lineheight = 0.95) +
    ggplot2::scale_fill_manual(values = c(Up = DIV_UP, Down = DIV_DOWN, NS = DIV_NS), guide = "none") +
    ggplot2::scale_x_continuous(breaks = 1, labels = "Gene", limits = c(0.5, 2.3)) +
    ggplot2::scale_y_continuous(limits = dge_ylim) +
    ggplot2::labs(x = NULL, y = "Gene log2FC") +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(panel.grid.minor = ggplot2::element_blank(), panel.grid.major = ggplot2::element_blank())

  # --- panel 1b: gene-level normalized counts, next to the log2FC bar it's
  # derived from -- same DESeq2 fit as build_gene_dge(), no extra model run.
  # Built from the per-sample normCount.<sample> columns (not just the two
  # per-condition means) so replicate spread is visible, the same way the
  # PAS-usage panel below shows jittered per-replicate points.
  se_fn2 <- function(x) { x <- x[!is.na(x)]; if (length(x) > 1) sd(x) / sqrt(length(x)) else NA_real_ }
  norm_cols <- grep("^normCount\\.", names(dge_row), value = TRUE)
  count_rep <- if (nrow(dge_row) && length(norm_cols) > 0) {
    data.frame(sample    = sub("^normCount\\.", "", norm_cols),
               count     = as.numeric(dge_row[1, norm_cols]),
               stringsAsFactors = FALSE)
  } else {
    data.frame(sample = character(0), count = numeric(0))
  }
  count_rep$condition <- factor(sub_cd[count_rep$sample, comparison_var], levels = c(ref_label, alt_label))
  count_summary <- data.frame(condition = factor(c(ref_label, alt_label), levels = c(ref_label, alt_label))) %>%
    dplyr::left_join(
      count_rep %>% dplyr::group_by(condition) %>%
        dplyr::summarise(mean = mean(count, na.rm = TRUE), se = se_fn2(count), .groups = "drop"),
      by = "condition"
    )

  p_counts <- ggplot2::ggplot(count_summary, ggplot2::aes(x = condition, y = mean, fill = condition)) +
    ggplot2::geom_col(width = 0.6, alpha = 0.5, show.legend = FALSE) +
    ggplot2::geom_point(data = count_rep, ggplot2::aes(x = condition, y = count, color = condition),
                         inherit.aes = FALSE, position = ggplot2::position_jitter(width = 0.08, height = 0),
                         size = 1.4, alpha = 0.8, show.legend = FALSE) +
    ggplot2::geom_errorbar(ggplot2::aes(ymin = pmax(mean - se, 0), ymax = mean + se, color = condition),
                            width = 0.2, linewidth = 0.6, na.rm = TRUE, show.legend = FALSE) +
    ggplot2::scale_fill_manual(values = cond_colors) +
    ggplot2::scale_color_manual(values = cond_colors) +
    ggplot2::labs(x = NULL, y = "Expression") +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(panel.grid.minor = ggplot2::element_blank(), panel.grid.major.x = ggplot2::element_blank())

  # --- panel 2: PAS usage -- replicate jitter + mean line + SEM error bars ---
  usage_long <- dplyr::bind_rows(
    dplyr::transmute(pas_df, pos, condition = ref_label, val = .data[[ctrl_col]], se = .data[[se_ctrl_col]]),
    dplyr::transmute(pas_df, pos, condition = alt_label, val = .data[[trt_col]],  se = .data[[se_trt_col]])
  ) %>% dplyr::mutate(val = ifelse(is.na(val), 0, val))

  cnt_cols   <- grep("^countData\\.", names(pas_df), value = TRUE)
  counts_mat <- as.matrix(pas_df[, cnt_cols, drop = FALSE])
  colnames(counts_mat) <- sub("^countData\\.", "", colnames(counts_mat))
  gene_tot   <- colSums(counts_mat)
  usage_mat  <- sweep(counts_mat, 2, pmax(gene_tot, 1), "/")
  usage_mat[, gene_tot == 0] <- NA
  rep_long <- as.data.frame(usage_mat) %>%
    dplyr::mutate(pos = pas_df$pos) %>%
    tidyr::pivot_longer(-pos, names_to = "sample", values_to = "usage") %>%
    dplyr::mutate(condition = sub_cd[sample, comparison_var]) %>%
    dplyr::filter(!is.na(usage))

  jitter_w <- diff(xr) * 0.004
  dodge_w  <- diff(xr) * 0.01

  # Height to anchor each PAS's significance label at: just above whichever
  # condition's point/error-bar reaches highest at that position (not a fixed
  # spot), so the label always sits right next to what it's actually marking
  # instead of floating somewhere it might collide with -- or be far from --
  # the data.
  star_height <- usage_long %>%
    dplyr::mutate(top = ifelse(is.na(se), val, pmin(val + se, 1))) %>%
    dplyr::group_by(pos) %>%
    dplyr::summarise(top = max(top, na.rm = TRUE), .groups = "drop")

  # Per-PAS significance stars (1-3, DEXSeq padj, same thresholds used
  # elsewhere in this script). Given as an opaque label (not plain text) so
  # it stays legible against the dotted guide line and data points behind it
  # regardless of color or position -- plain colored text kept blending into
  # both.
  star_df <- pas_df %>%
    dplyr::mutate(n_star = dplyr::case_when(
      is.na(padj)  ~ 0L,
      padj < 0.001 ~ 3L,
      padj < 0.01  ~ 2L,
      padj < 0.05  ~ 1L,
      TRUE         ~ 0L
    )) %>%
    dplyr::filter(n_star > 0) %>%
    dplyr::mutate(label = strrep("★", n_star)) %>%
    dplyr::left_join(star_height, by = "pos") %>%
    dplyr::mutate(y_star = pmin(top + 0.05, 0.97))

  x_scale_usage <- if (flip_x) {
    ggplot2::scale_x_reverse(limits = xr, labels = scales::comma)
  } else {
    ggplot2::scale_x_continuous(limits = xr, labels = scales::comma)
  }

  p_usage <- ggplot2::ggplot() +
    vlines +
    ggplot2::geom_point(data = rep_long, ggplot2::aes(x = pos, y = usage, color = condition),
                         position = ggplot2::position_jitterdodge(dodge.width = dodge_w, jitter.width = jitter_w, seed = 1),
                         size = 1, alpha = 0.45, show.legend = FALSE) +
    ggplot2::geom_errorbar(data = usage_long,
                            ggplot2::aes(x = pos, ymin = pmax(val - se, 0), ymax = pmin(val + se, 1), color = condition),
                            width = dodge_w, linewidth = 0.6, position = ggplot2::position_dodge(width = dodge_w * 2)) +
    ggplot2::geom_line(data = usage_long, ggplot2::aes(x = pos, y = val, color = condition), linewidth = 0.7) +
    ggplot2::geom_point(data = usage_long, ggplot2::aes(x = pos, y = val, color = condition), size = 1.8) +
    ggplot2::geom_label(data = star_df, ggplot2::aes(x = pos, y = y_star, label = label), inherit.aes = FALSE,
                         vjust = 0.1, size = 3, fontface = "bold", color = "red", fill = "white",
                         linewidth = 0.3, label.padding = ggplot2::unit(0.12, "lines"),
                         label.r = ggplot2::unit(0.08, "lines")) +
    ggplot2::scale_color_manual(values = cond_colors, name = comparison_label) +
    ggplot2::scale_y_continuous(labels = scales::percent_format(accuracy = 1), limits = c(0, 1)) +
    x_scale_usage +
    ggplot2::labs(y = "PAS usage", x = NULL) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(panel.grid.minor = ggplot2::element_blank(), panel.grid.major.x = ggplot2::element_blank(),
                   axis.text.x = ggplot2::element_blank(), axis.ticks.x = ggplot2::element_blank(),
                   legend.position = "inside", legend.position.inside = c(0.92, 0.88),
                   legend.background = ggplot2::element_rect(fill = scales::alpha("white", 0.75), color = NA),
                   legend.title = ggplot2::element_text(size = 9), legend.text = ggplot2::element_text(size = 8),
                   legend.key.size = grid::unit(0.8, "lines"))

  # --- panel 3: gene model ---
  # ggtranscript's geom_range()/geom_intron() don't render correctly under
  # scale_x_reverse() -- inspecting ggplot_build() output shows their internal
  # xstart/xend -> xmin/xmax conversion only routes ONE of the two endpoints
  # through the scale's transform, producing an inverted (xmin > xmax)
  # rectangle that draws as one solid block spanning the whole panel instead
  # of a proper exon shape. So minus-strand genes mirror the DATA itself here
  # (negate + swap start/end) and plot it on a perfectly ordinary,
  # untransformed continuous scale, with tick labels that un-negate back to
  # the real coordinate -- visually identical to a reversed axis, but
  # something ggtranscript's custom geoms can actually render. After this
  # mirroring, transcription always proceeds in the increasing-x direction
  # regardless of the gene's real strand, so geom_intron() is always told
  # strand = "+" here.
  mirror <- function(d) {
    if (!flip_x || nrow(d) == 0) return(d)
    tmp     <- d$start
    d$start <- -d$end
    d$end   <- -tmp
    d
  }
  ctx_exons_disp <- mirror(ctx_exons)
  ext_disp       <- mirror(ext)
  introns_disp   <- ggtranscript::to_intron(ctx_exons_disp, "transcript_id")
  xr_disp        <- if (flip_x) sort(-xr) else xr
  pas_disp       <- if (flip_x) -pas_df$pos else pas_df$pos
  vlines_gene    <- ggplot2::geom_vline(xintercept = pas_disp, linetype = "dotted", color = "grey30", linewidth = 0.75)

  y_pos <- setNames(seq_along(tx_order), tx_order)
  ctx_exons_disp$y <- y_pos[ctx_exons_disp$transcript_id]
  if (nrow(introns_disp) > 0) introns_disp$y <- y_pos[introns_disp$transcript_id]
  if (nrow(ext_disp) > 0) ext_disp$y <- y_pos[ext_disp$transcript_id]

  p_gene <- ggplot2::ggplot() +
    vlines_gene +
    ggtranscript::geom_range(data = ctx_exons_disp, ggplot2::aes(xstart = start, xend = end, y = y),
                              fill = "grey35", height = 0.4)
  if (nrow(introns_disp) > 0)
    p_gene <- p_gene + ggtranscript::geom_intron(data = introns_disp, ggplot2::aes(xstart = start, xend = end, y = y),
                                                  strand = "+", arrow.min.intron.length = 200, color = "grey50")
  if (nrow(ext_disp) > 0)
    p_gene <- p_gene + ggtranscript::geom_range(data = ext_disp, ggplot2::aes(xstart = start, xend = end, y = y),
                                                 fill = NA, color = "grey35", linetype = "dashed", height = 0.4)

  p_gene <- p_gene +
    ggplot2::scale_y_continuous(breaks = y_pos, labels = names(y_pos), limits = c(0.5, max(y_pos) + 0.5)) +
    ggplot2::scale_x_continuous(labels = function(x) scales::comma(if (flip_x) -x else x)) +
    ggplot2::coord_cartesian(xlim = xr_disp) +
    ggplot2::labs(x = paste0("Genomic position (", tx_exons$seqnames[1], strand_sym, ")"), y = NULL) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(panel.grid.minor = ggplot2::element_blank(), panel.grid.major = ggplot2::element_blank())

  ttl <- paste0(gene_symbol, if (!is.null(title_suffix)) paste0(" — ", title_suffix) else "")

  # Gene-model panel height is capped rather than growing without bound --
  # for genes with many isoforms it was crowding out the DGE/counts/usage
  # panels above, which carry the more decision-relevant information.
  gene_panel_height <- min(0.18 + 0.12 * length(tx_order), 1.3)

  ((p_dge | p_counts) / p_usage / p_gene) +
    patchwork::plot_layout(heights = c(0.45, 1.3, gene_panel_height)) +
    patchwork::plot_annotation(title = ttl)
}

# ============================================================
#  PRE-ANALYSIS QC
# ============================================================

message("--- Pre-analysis QC ---")

make_library_sizes(
  count_mat_raw, colData,
  outfile = file.path(dir_qc, "library_sizes.png")
)

if (USE_RUV) {
  # W factor bar chart — shows how much each sample is pulled by the batch factor
  w_cols <- paste0("W_", seq_len(RUV_K))
  w_df   <- tibble::rownames_to_column(
               as.data.frame(colData[, c("condition", w_cols), drop = FALSE]), "sample")
  w_long <- tidyr::pivot_longer(w_df, cols = dplyr::all_of(w_cols),
                                 names_to = "factor", values_to = "value")
  gw <- ggplot2::ggplot(w_long, ggplot2::aes(x = sample, y = value, fill = condition)) +
    ggplot2::geom_col() +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
    ggplot2::facet_wrap(~ factor, ncol = 1, scales = "free_y") +
    ggplot2::labs(x = "Sample", y = "W (RUV factor)",
                  title = "RUVs estimated batch factors") +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
  ggplot2::ggsave(file.path(dir_qc, "RUV_W_factors.png"), gw,
                  width = max(4, 0.5 * nrow(colData) + 2),
                  height = 3 * RUV_K + 1, dpi = 300)

  # Uncorrected PCA — shows the raw batch structure you are correcting for
  make_pca_global(count_mat, colData, shape_col = PCA_SHAPE_VAR,
                  outfile        = file.path(dir_qc, "PCA_global_uncorrected.png"),
                  matrix_outfile = file.path(dir_qc, "vst_global.csv"))
  # RUV-corrected PCA — W factors regressed from VST matrix for visualization only
  make_pca_global(count_mat, colData, shape_col = PCA_SHAPE_VAR,
                  batch_covars = w_cols,
                  outfile = file.path(dir_qc, "PCA_global_RUV_corrected.png"))
} else {
  make_pca_global(count_mat, colData, shape_col = PCA_SHAPE_VAR,
                  outfile        = file.path(dir_qc, "PCA_global.png"),
                  matrix_outfile = file.path(dir_qc, "vst_global.csv"))
}

# ============================================================
#  FUNCTION: cache wrapper for expensive per-group model fits
# ============================================================
#
# run_dexseq_group()/build_gene_dge() are the dominant cost in this script --
# full DEXSeq/DESeq2 model fits across every tested PAS. Most re-runs of this
# script are iterating on a plot or GENES_OF_INTEREST, not changing the input
# data or modeling config, so refitting from scratch every time is wasted
# work. with_cache() saves a compute_fn()'s result plus a small manifest
# (input-file mtimes + the modeling config that actually affects the fit) to
# dir_cache; on the next run, if the manifest is unchanged, the cached result
# is returned instead of recomputing. USE_CACHE <- FALSE always recomputes.

file_mtime_safe <- function(path) {
  if (is.null(path) || !file.exists(path)) return(NA_real_)
  as.numeric(file.info(path)$mtime)
}

with_cache <- function(cache_key, manifest, compute_fn,
                        cache_dir = dir_cache, use_cache = USE_CACHE) {
  safe_key      <- gsub("[^A-Za-z0-9_.-]", "_", cache_key)
  rds_file      <- file.path(cache_dir, paste0(safe_key, ".rds"))
  manifest_file <- file.path(cache_dir, paste0(safe_key, "_manifest.rds"))

  if (use_cache && file.exists(rds_file) && file.exists(manifest_file) &&
      identical(readRDS(manifest_file), manifest)) {
    message(sprintf("[cache] Reusing cached result for '%s' (inputs/config unchanged).", cache_key))
    return(readRDS(rds_file))
  }

  if (use_cache) message(sprintf("[cache] No valid cache for '%s' -- computing.", cache_key))
  result <- compute_fn()
  saveRDS(result, rds_file)
  saveRDS(manifest, manifest_file)
  result
}

# Manifest shared by both cached steps below -- build_gene_dge() operates
# entirely on the dxd that run_dexseq_group() produces, so anything that
# would invalidate one fit invalidates the other identically.
group_cache_manifest <- function(axis_name, axis_cfg, grp, group_list) {
  list(
    inputs = list(
      counts  = file_mtime_safe(path_counts),
      anno    = file_mtime_safe(path_anno),
      design  = file_mtime_safe(path_design),
      samples = file_mtime_safe(path_samples)
    ),
    axis_name    = axis_name,
    samples_used = sort(group_list[[grp]]),
    config = list(
      comparison_var = axis_cfg$comparison_var,
      ref_label = axis_cfg$ref_label, alt_label = axis_cfg$alt_label,
      PAS_TYPE_REGEX = PAS_TYPE_REGEX, MIN_TOTAL_READS = MIN_TOTAL_READS,
      MIN_PER_CONDITION = MIN_PER_CONDITION, USE_RUV = USE_RUV, RUV_K = RUV_K
    )
  )
}

# ============================================================
#  FUNCTION: run one comparison axis end-to-end -- grouping, cached DEXSeq
#  fits, gene-level DGE, QC plots, results export, and APA figures -- all
#  written into a nested <axis_name>/<stratify level>/ output tree.
# ============================================================

run_axis <- function(axis_name, axis_cfg) {
  message(sprintf("=== Axis '%s' (%s: %s vs %s, stratified by %s) ===",
                  axis_name, axis_cfg$comparison_var, axis_cfg$alt_label,
                  axis_cfg$ref_label, axis_cfg$stratify_var))

  fc_col    <- paste0("log2fold_", axis_cfg$alt_label, "_", axis_cfg$ref_label)
  fc_export <- paste0("log2fold_", axis_cfg$alt_label, "_v_", axis_cfg$ref_label)
  # Plot-title comparison label: alt_label first, ref_label second, matching the
  # fold-change sign convention used everywhere in this script (estimateExonFoldChanges'
  # denominator = ref_label, and build_gene_dge()'s contrast = c(var, alt_label,
  # ref_label)) -- a positive log2FC/log2FoldChange always means "up in alt_label",
  # so naming it "<alt> v <ref>" keeps title order and sign meaning consistent.
  comparison_title <- paste0(axis_cfg$alt_label, " v ", axis_cfg$ref_label)
  usage_ctrl_col     <- paste0("meanUsage_", axis_cfg$ref_label)
  usage_trtmt_col    <- paste0("meanUsage_", axis_cfg$alt_label)
  se_ctrl_col        <- paste0("seUsage_", axis_cfg$ref_label)
  se_trtmt_col       <- paste0("seUsage_", axis_cfg$alt_label)
  red_ctrl_col       <- paste0("RED_", axis_cfg$ref_label)
  red_trtmt_col      <- paste0("RED_", axis_cfg$alt_label)
  wutr_ctrl_col      <- paste0("wUTR_", axis_cfg$ref_label)
  wutr_trtmt_col     <- paste0("wUTR_", axis_cfg$alt_label)
  dom_usage_ctrl_col <- paste0("dominant_meanUsage_", axis_cfg$ref_label)
  dom_usage_trtmt_col<- paste0("dominant_meanUsage_", axis_cfg$alt_label)

  # --- grouping setup (stratify samples by axis_cfg$stratify_var) ---
  stratify_levels <- axis_cfg$stratify_levels
  group_list <- setNames(
    lapply(stratify_levels, function(lv) {
      rownames(colData)[as.character(colData[[axis_cfg$stratify_var]]) == lv]
    }),
    stratify_levels
  )
  level_dir <- setNames(vapply(stratify_levels, function(lv) {
    d <- axis_cfg$stratify_dir[lv]
    if (is.null(axis_cfg$stratify_dir) || is.na(d)) lv else unname(d)
  }, character(1)), stratify_levels)

  for (lv in stratify_levels) {
    for (base in c(dir_results, dir_qc, dir_apa_genome, dir_apa_zoom)) {
      dir.create(file.path(base, axis_name, level_dir[[lv]]), recursive = TRUE, showWarnings = FALSE)
    }
  }

  # ---- cached DEXSeq fits ----
  message(sprintf("[%s] Running DEXSeq analyses...", axis_name))
  runs <- lapply(names(group_list), function(grp) {
    with_cache(
      cache_key = paste0("dexseq_", axis_name, "_", grp),
      manifest  = group_cache_manifest(axis_name, axis_cfg, grp, group_list),
      compute_fn = function() {
        run_dexseq_group(
          grp_label         = paste(axis_name, grp, sep = "/"),
          sub_samples       = group_list[[grp]],
          count_mat         = count_mat,
          colData           = colData,
          featureID         = featureID,
          groupID           = groupID,
          comparison_var    = axis_cfg$comparison_var,
          ref_label         = axis_cfg$ref_label,
          alt_label         = axis_cfg$alt_label,
          gr                = gr,
          pas_anno          = pas_anno,
          min_total         = MIN_TOTAL_READS,
          min_per_condition = MIN_PER_CONDITION,
          use_ruv           = USE_RUV,
          ruv_k             = RUV_K
        )
      }
    )
  })
  names(runs) <- names(group_list)

  # ---- gene-level DGE (separate from DEXSeq's per-PAS exon-usage test) ----
  message(sprintf("[%s] Running gene-level DGE...", axis_name))
  gene_dge_results <- lapply(names(runs), function(grp) {
    with_cache(
      cache_key  = paste0("genedge_", axis_name, "_", grp),
      manifest   = group_cache_manifest(axis_name, axis_cfg, grp, group_list),
      compute_fn = function() build_gene_dge(runs[[grp]]$dxd,
                     comparison_var = axis_cfg$comparison_var,
                     ctrl_label = axis_cfg$ref_label, trt_label = axis_cfg$alt_label)
    )
  })
  names(gene_dge_results) <- names(runs)

  for (grp in names(gene_dge_results)) {
    write.csv(gene_dge_results[[grp]],
              file.path(dir_results, axis_name, level_dir[[grp]], "gene_dge.csv"),
              row.names = FALSE, quote = FALSE)
  }

  # ---- post-analysis QC plots (per stratify level) ----
  message(sprintf("[%s] Post-analysis QC plots...", axis_name))
  for (grp in names(runs)) {
    qdir <- file.path(dir_qc, axis_name, level_dir[[grp]])
    dxd  <- runs[[grp]]$dxd
    res  <- runs[[grp]]$results
    grp_label <- paste0(comparison_title, ", ", grp)

    make_pca_dxd(
      dxd, grp_label = grp_label, color_col = axis_cfg$comparison_var,
      outfile        = file.path(qdir, "PCA_normalized.png"),
      matrix_outfile = file.path(qdir, "vst_normalized.csv")
    )
    make_size_factors(
      dxd, grp_label = grp_label, comparison_var = axis_cfg$comparison_var,
      outfile = file.path(qdir, "size_factors.png")
    )
    make_dispersion_plot(
      dxd, grp_label = grp_label,
      outfile = file.path(qdir, "dispersion.png")
    )
    make_pvalue_hist(
      res, grp_label = grp_label,
      outfile = file.path(qdir, "pvalue_hist.png")
    )
    make_ma_plot(
      res, grp_label = grp_label, fc_col = fc_col,
      outfile = file.path(qdir, "MA.png")
    )
    make_volcano_plot(
      res, grp_label = grp_label, fc_col = fc_col,
      outfile = file.path(qdir, "volcano.png")
    )
  }

  message(sprintf("[%s] Exporting results...", axis_name))

  # ---- per-PAS results + usage (DEXSeq results, annotation, raw/normalized
  #      counts, and mean PSI + SEM per condition) ----
  usage_results <- lapply(names(group_list), function(grp) {
    cd_sub <- colData[group_list[[grp]], , drop = FALSE]
    usage_by_condition(runs[[grp]]$results, cd_sub,
                        comparison_var = axis_cfg$comparison_var,
                        ref_label = axis_cfg$ref_label, alt_label = axis_cfg$alt_label,
                        usage_ctrl_col = usage_ctrl_col, usage_trtmt_col = usage_trtmt_col,
                        se_ctrl_col = se_ctrl_col, se_trtmt_col = se_trtmt_col,
                        red_ctrl_col = red_ctrl_col, red_trtmt_col = red_trtmt_col)
  })
  names(usage_results) <- names(group_list)

  for (grp in names(usage_results)) {
    out <- usage_results[[grp]]
    names(out) <- sub(fc_col, fc_export, names(out), fixed = TRUE)
    write.csv(out,
              file.path(dir_results, axis_name, level_dir[[grp]], "pas_usage.csv"),
              row.names = FALSE, quote = FALSE)
  }

  # ---- gene-level summary (collapse_gene() significance rollup, joined with
  #      weighted UTR-length shift, dominant-site tracking, RED rollup, and the
  #      chi-squared cross-check from build_gene_apa_summary()) + candidate genes ----
  message(sprintf("[%s] Building gene-level summary and candidate-gene list...", axis_name))

  for (grp in names(runs)) {
    gene_q <- perGeneQValue(runs[[grp]]$dxr)
    group_of_sample <- setNames(as.character(colData[group_list[[grp]], axis_cfg$comparison_var]),
                                 group_list[[grp]])

    sig_summary <- collapse_gene(runs[[grp]]$results, gene_q,
                                  padj_cut = PADJ_CUT, lfc_cut = LFC_CUT, fc_col = fc_col)
    apa_summary <- build_gene_apa_summary(usage_results[[grp]], group_of_sample,
                     ctrl_col = usage_ctrl_col, trt_col = usage_trtmt_col,
                     ref_label = axis_cfg$ref_label, alt_label = axis_cfg$alt_label,
                     wutr_ctrl_col = wutr_ctrl_col, wutr_trtmt_col = wutr_trtmt_col,
                     dom_usage_ctrl_col = dom_usage_ctrl_col, dom_usage_trtmt_col = dom_usage_trtmt_col)
    apa_summary <- apa_summary[, setdiff(names(apa_summary), "n_PAS"), drop = FALSE]  # already in sig_summary

    gene_summary <- merge(sig_summary, apa_summary, by = "groupID", all = TRUE, sort = FALSE)
    write.csv(gene_summary,
              file.path(dir_results, axis_name, level_dir[[grp]], "gene_summary.csv"),
              row.names = FALSE, quote = FALSE)

    candidate_genes <- build_candidate_genes(usage_results[[grp]], gene_summary,
                          ctrl_col = usage_ctrl_col, trt_col = usage_trtmt_col,
                          usage_change_cut = USAGE_CHANGE_CUT, padj_cut = PADJ_CUT,
                          wutr_change_cut = WUTR_CHANGE_CUT, red_mag_cut = RED_MAG_CUT)
    write.csv(candidate_genes,
              file.path(dir_results, axis_name, level_dir[[grp]], "candidate_genes.csv"),
              row.names = FALSE, quote = FALSE)
  }

  # ---- APA genome-map + terminal-exon zoom figures for GENES_OF_INTEREST ----
  if (length(GENES_OF_INTEREST) > 0) {
    message(sprintf("[%s] Plotting APA genome-map figures...", axis_name))

    for (g in GENES_OF_INTEREST) {
      for (grp in names(usage_results)) {
        ttl_sfx <- paste0(comparison_title, ", ", grp)

        if (!is.null(gtf_exons) && requireNamespace("ggtranscript", quietly = TRUE)) {
          n_tx <- length(unique(gtf_exons$transcript_id[gtf_exons$gene_name == g]))

          if (requireNamespace("ggforce", quietly = TRUE)) {
            genome_outfile <- file.path(dir_apa_genome, axis_name, level_dir[[grp]], sprintf("%s.png", g))
            tryCatch({
              p <- plot_gene_apa_genome(usage_results[[grp]], g, gtf_exons, title_suffix = ttl_sfx)
              ggplot2::ggsave(genome_outfile, p, width = 10, height = 3 + 1.1 * n_tx, dpi = 300)
            }, error = function(e) {
              warning(sprintf("[%s/%s] Could not plot APA genome map for '%s': %s", axis_name, grp, g, conditionMessage(e)))
            })
          } else {
            warning("ggforce not installed — APA genome map figure skipped for '", g, "'.")
          }

          if (requireNamespace("patchwork", quietly = TRUE)) {
            zoom_outfile <- file.path(dir_apa_zoom, axis_name, level_dir[[grp]], sprintf("%s.png", g))
            tryCatch({
              cd_sub <- colData[group_list[[grp]], , drop = FALSE]
              p <- plot_gene_apa_zoom(usage_results[[grp]], g, gtf_exons, cd_sub, gene_dge_results[[grp]],
                     title_suffix = ttl_sfx,
                     ctrl_col = usage_ctrl_col, trt_col = usage_trtmt_col,
                     se_ctrl_col = se_ctrl_col, se_trt_col = se_trtmt_col,
                     comparison_var = axis_cfg$comparison_var,
                     ref_label = axis_cfg$ref_label, alt_label = axis_cfg$alt_label,
                     comparison_label = axis_cfg$comparison_label)
              ggplot2::ggsave(zoom_outfile, p, width = 9, height = 4.5 + 0.3 * n_tx, dpi = 300)
            }, error = function(e) {
              warning(sprintf("[%s/%s] Could not plot APA terminal-exon zoom for '%s': %s", axis_name, grp, g, conditionMessage(e)))
            })
          } else {
            warning("patchwork not installed — APA terminal-exon zoom figure skipped for '", g, "'.")
          }
        }
      }
    }
  }

  invisible(list(runs = runs, gene_dge_results = gene_dge_results, usage_results = usage_results))
}

# ============================================================
#  MAIN LOOP -- run every configured axis
# ============================================================

for (axis_name in names(AXES)) run_axis(axis_name, AXES[[axis_name]])

message("Done. Outputs written to: ", OUT_BASE)
