# DEXseq.R — Differential Alternative Polyadenylation Analysis

**Script:** `scripts/DEXseq.R`
**Last updated:** 2026-08-19

---

## What this script does — and why

### The biological question

After the 3-REAP pipeline runs, we know how many RNA molecules from each sample ended at each poly-A site (PAS). But raw counts alone don't answer the biological question. What we want to know is: **does the relative usage of each PAS site change between the control and treatment conditions?**

Consider a gene that has two PAS: a proximal one (shorter 3' UTR) and a distal one (longer 3' UTR). If treatment causes cells to preferentially use the distal site, that gene is said to have undergone **3' UTR lengthening**. This shifts regulatory control — the longer 3' UTR exposes binding sites for microRNAs and RNA-binding proteins that weren't accessible before, potentially changing how much protein is made, how stable the mRNA is, or where in the cell it goes.

The fraction of reads from a gene that end at a specific PAS is called **PSI (Percent Site Inclusion)** or **fractional usage**. A PSI of 0.7 at the distal site means 70% of that gene's poly-A reads end at the distal PAS. We test whether PSI changes significantly between conditions.

### The statistical method: DEXSeq

This script uses **DEXSeq**, a Bioconductor package originally developed to detect differential exon usage in RNA-seq data. We repurpose it for APA by treating each PAS as if it were an exon within its gene. DEXSeq's core insight is that comparing *relative* usage (each PAS as a fraction of its gene's total counts) is statistically more powerful than comparing absolute counts, because it controls for variation in overall gene expression.

The DEXSeq model tests each PAS independently with the design `~ sample + exon + condition:exon`:
- `sample` absorbs per-sample baseline differences (the blocking factor)
- `condition:exon` is the actual test: does the Treatment condition change the fractional usage of this PAS?

A significant result (adjusted p-value < threshold) means the relative usage of that PAS is statistically different between control and treatment, after accounting for how many reads the gene received overall.

### Usage direction: Increased vs. Decreased, and the gene-level APA summary

For each PAS, `APA_direction` reports only the sign of *that PAS's own* usage change — `Increased` or `Decreased` — not an overall UTR-lengthening/shortening call. Earlier versions of this script tried to infer "Lengthened"/"Shortened" from a distal-vs-everything-else binary, but that breaks down once a gene has 3+ PAS: a usage gain at a middle site doesn't map cleanly onto "shorter" or "longer" 3' UTR, only onto "distal vs. not distal." `is_distal` (still in the output) tells you whether a given PAS is the gene's most 3' site, as positional context.

The actual lengthening/shortening narrative is captured properly, across all of a gene's PAS at once, in `gene_summary.csv` via `delta_wUTR` (weighted UTR length) — see that section below. `dir_consensus` in the same table (per-gene rollup of `APA_direction`) is now `Increased_only` / `Decreased_only` / `Mixed`.

---

## Quick-start: how to run the script

### 1. Edit the configuration block

Open `DEXseq.R`. Change only the `USER CONFIGURATION — EDIT THESE` section at the top:

```r
BASE_DIR    <- "/path/to/your/project"          # project root
COUNTS_CSV  <- "path/to/cluster.all.reads.csv"          # PAS count matrix from 3-REAP
SAMPLES_TXT <- "path/to/sample_name.txt"                # sample list (or set to NULL)
ANNO_FILE   <- "data/PolyA_DB_v4.1/hg38.PAS.main.tsv"  # same file used by 3-REAP
DESIGN_FILE <- "path/to/design.txt"                     # experimental design table

ANNO_GTF <- "data/ref/hg38.refGene.gtf"   # RefSeq GTF for the APA genome-map figures;
                                          # set to NULL to skip those figures entirely

CTRL_LABEL  <- "Control"     # must match values in design file's 'condition' column
TRTMT_LABEL <- "Treatment"   # must match values in design file's 'condition' column

GROUPING_VAR <- NULL   # NULL = run one analysis on all samples
                       # set to "timepoint" (or any colData column) to run per level

GENES_OF_INTEREST <- c("GENE1", "GENE2")  # genes to plot APA genome-map figures for

USE_RUV <- FALSE   # set to TRUE to enable RUVseq batch correction (see below)
RUV_K   <- 1       # number of unwanted variation factors to estimate (1 is usually enough)

USE_CACHE <- TRUE  # reuse cached DEXSeq/DESeq2 fits when inputs/config are unchanged (see below)
```

### 2. Prepare input files

**Count matrix** (`cluster.all.reads.csv`): produced by 3-REAP step 6. A CSV where the first column is `hit_PAS_ID` and the remaining columns are samples.

**PAS annotation file** (`data/PolyA_DB_v4.1/hg38.PAS.main.tsv`): use the **same file** that was provided to the 3-REAP pipeline. The `hit_PAS_ID` values in the count matrix are taken directly from this database, so the IDs will match. The script auto-detects the PolyA_DB format:
- v3.2: explicit `Chromosome`, `Position`, `Strand` columns
- v4.1: coordinates encoded in `PAS_ID` as `chr:strand:position` (auto-parsed; default for this project)

**Design file** (`design.txt`): a tab-delimited file, first column = sample names (matching column names in the count matrix), remaining columns = metadata. Must include at minimum:

```
sample_name    condition
Control_rep1   Control
Control_rep2   Control
Treatment_rep1 Treatment
Treatment_rep2 Treatment
```

If you have a timepoint or batch column, include it here. Set `GROUPING_VAR <- "timepoint"` to run a separate DEXSeq analysis per timepoint level.

**Sample name file** (`sample_name.txt`, optional): a plain text file listing sample names, one per line. If set to NULL, all non-`hit_PAS_ID` columns in the count matrix are used.

**GTF file** (`data/ref/hg38.refGene.gtf`, optional — controlled by `ANNO_GTF`): a RefSeq-format GTF, ideally the **same one used to build your STAR index** (check `data/ref/star_index_hg38/genomeParameters.txt` for the `--sjdbGTFfile` path it was built with, if you're not sure which one that is). This is what gives the APA genome-map figures (below) real exon/intron structure — PolyA_DB alone only has flat PAS positions, no gene structure. Set `ANNO_GTF <- NULL` to skip those two figures entirely; nothing else in the pipeline depends on this file.

### 3. Run

```r
Rscript scripts/DEXseq.R
```

Or open in RStudio and run interactively. The script prints progress messages for each step.

---

## Software requirements (R packages)

Install via Bioconductor/CRAN if missing:

```r
BiocManager::install(c("DEXSeq", "DESeq2", "GenomicRanges", "SummarizedExperiment"))
install.packages(c("ggplot2", "ggrepel", "scales", "limma", "dplyr",
                   "tidyr", "tibble", "matrixStats", "stringr"))

# Optional — only needed for the APA genome-map figures (plots/apa_genome/, plots/apa_zoom/).
# Both are skipped with a warning if ANNO_GTF is unset/missing or any package below
# isn't installed; nothing else in the pipeline depends on them.
BiocManager::install("rtracklayer")          # parses the GTF into exon coordinates
install.packages("patchwork")                # stacks the zoom figure's 3 aligned panels
install.packages("ggforce")                  # main+zoom-inset connector on the genome map figure

# Not on CRAN/Bioconductor:
remotes::install_github("dzhang32/ggtranscript")
```

---

## Understanding the analysis: step by step

### Filtering
Before testing, PAS sites are filtered to retain only:
1. **3' UTR PAS** — sites in the 3' untranslated region (set by `PAS_TYPE_REGEX`). Intronic and other PAS are excluded because APA in 3' UTRs is the biologically relevant type for most questions.
2. **Sufficient read coverage** — PAS with fewer than `MIN_TOTAL_READS` supporting reads in a group are removed. This prevents testing sites where we have too little evidence.
3. **Genes with ≥2 PAS** — DEXSeq tests relative usage within genes; a gene with only one observable PAS cannot have differential isoform usage.

### Normalization
DEXSeq estimates a **size factor** for each sample — a scaling factor that accounts for differences in sequencing depth (some samples have more reads than others). All counts are internally divided by their sample's size factor before testing.

### Testing
DEXSeq fits a statistical model to each PAS and tests whether the Treatment:exon interaction term is significantly non-zero — i.e., whether the fractional usage of that PAS changes in treatment vs. control. This produces a p-value per PAS. Multiple testing correction (Benjamini-Hochberg FDR) is applied to give an adjusted p-value (`padj`).

Additionally, the Lancaster procedure (`perGeneQValue`) combines per-PAS p-values into a per-gene q-value that summarizes whether *any* PAS in the gene shows significant differential usage.

### Fold changes
DEXSeq reports a **log2 fold change** for each PAS: the log2 ratio of normalized fractional usage (Treatment / Control). However, for biological interpretation we report **ΔPSI** (Treatment PSI − Control PSI) in the plots, because it is directly interpretable as a percentage-point change in isoform usage.

A ΔPSI of +0.20 at the distal PAS means the distal site gained 20 percentage points of relative usage in treatment.

---

## Output files

All outputs are written under `OUT_BASE` (set in configuration). Structure:

```
{OUT_BASE}/
├── pas_usage[.group].csv        ← per-PAS DEXSeq results + raw/normalized counts + usage (PSI) + SEM + RED
├── gene_summary[.group].csv     ← per-gene summary (significance, wUTR, dominant site, RED rollup, chi-squared)
├── candidate_genes[.group].csv  ← curated shortlist: genes clearing every rubric criterion at once
├── gene_dge[.group].csv         ← per-gene differential expression (DESeq2, independent of DEXSeq)
├── cache/                   ← cached DEXSeq/DESeq2 fits (see "Caching" below); safe to delete
└── plots/
    ├── qc/
    │   ├── library_sizes.png
    │   ├── vst_global.csv                 ← full VST matrix, all samples (pre-DEXSeq)
    │   ├── PCA_global.png                  ← (USE_RUV=FALSE only)
    │   ├── PCA_global_uncorrected.png      ← (USE_RUV=TRUE only)
    │   ├── PCA_global_RUV_corrected.png    ← (USE_RUV=TRUE only; visual only)
    │   ├── RUV_W_factors.png               ← (USE_RUV=TRUE only)
    │   ├── vst_normalized[.group].csv      ← full VST matrix, post-DEXSeq-normalization, per group
    │   ├── PCA_normalized[.group].png
    │   ├── size_factors[.group].png
    │   ├── dispersion[.group].png
    │   ├── pvalue_hist[.group].png
    │   ├── MA[.group].png
    │   └── volcano[.group].png
    ├── apa_genome/
    │   └── {GENE}[.group].png   ← whole-gene genome map, all isoforms (needs ANNO_GTF + ggtranscript)
    └── apa_zoom/
        └── {GENE}[.group].png   ← terminal-exon zoom, same requirements
```

`[.group]` is appended only when `GROUPING_VAR` is set (e.g., `.3d` for a "3day" timepoint level). With `GROUPING_VAR <- NULL`, files have no suffix: `pas_usage.csv`, `gene_summary.csv`, etc.

---

### `gene_summary.csv` — per-gene summary

One row per gene. Joins two complementary views computed by `collapse_gene()` (a rollup of the per-PAS DEXSeq significance calls) and `build_gene_apa_summary()` (weighted UTR-length shift, dominant-site tracking, RED rollup, chi-squared cross-check — metrics that generalize cleanly to genes with 3+ PAS) into a single table, merged on `groupID`.

| Column | Meaning |
|---|---|
| `groupID` | Gene identifier |
| `n_PAS` | Total PAS tested for this gene |
| `n_sig` | Number of PAS with `padj < PADJ_CUT` |
| `n_big` | Number of significant PAS with `\|log2FC\| >= LFC_CUT` |
| `min_padj` | Most significant adjusted p-value among this gene's PAS |
| `max_absL2FC` | Largest absolute log2 fold change among significant PAS |
| `dir_consensus` | `Increased_only`, `Decreased_only`, or `Mixed` |
| `perGeneQ` | Lancaster-combined gene-level q-value (primary ranking metric) |
| `dominant_featureID` | The PAS with the highest overall mean usage (`meanUsage_All`) — i.e., the gene's "default" isoform |
| `dominant_abundance` | That PAS's `exonBaseMean` |
| `dominant_meanUsage_Control` / `_Treatment` | That PAS's mean PSI per condition |
| `dominant_delta_usage` | Treatment − Control PSI for the dominant PAS only |
| `dominant_padj` | DEXSeq `padj` for the dominant PAS only |
| `wUTR_Control` / `wUTR_Treatment` | Weighted UTR length (bp): usage-weighted average distance of the gene's PAS from its most-proximal surviving PAS, per condition: `sum(usage_i * distance_i)` |
| `delta_wUTR` | `wUTR_Treatment − wUTR_Control`. Positive = usage shifted toward more distal sites on average (net lengthening); negative = net shortening. Generalizes to any number of PAS, unlike the old distal-vs-everything-else `APA_direction` binary — there's no stop-codon coordinate in the PolyA_DB annotation, so the gene's own most-proximal PAS is used as the reference point instead. |
| `chisq_stat` / `chisq_pvalue` / `chisq_padj` | An independent cross-check: chi-squared test of association between PAS identity and condition, on summed raw counts (a PAS × condition contingency table). Doesn't depend on DEXSeq's negative-binomial model or size-factor normalization, so it isn't affected by anything upstream in that model — but it also doesn't correct for overdispersion, so expect it to run hotter (more significant) than `padj` on noisy/low-count genes. Use as a corroborating signal, not a replacement. |
| `single_dominant_site` | `TRUE` when this gene functionally has one PAS: the dominant PAS carries ≥90% of usage in **both** conditions, and that share doesn't shift between them. This is a **flag, not a filter** — rows are never dropped from this table; whether to exclude them from a given analysis is left to you. |
| `max_abs_delta_RED` / `max_delta_RED_featureID` | The single largest `\|delta_RED\|` among this gene's PAS (see `pas_usage.csv`'s "How RED is computed"), and which PAS it belongs to. Not necessarily `dominant_featureID` — the PAS driving the biggest relative-expression shift doesn't have to be the gene's most-abundant one. This is the direct replacement for ranking genes by RED the way the old MAAPER-based rubric did, generalized to any number of PAS. |

**How to read this table:** sort by `perGeneQ` ascending for the strongest overall statistical evidence of differential APA. `dir_consensus` tells you whether the change is directional (one dominant direction) or complex (multiple PAS moving in different directions). `delta_wUTR` and `max_abs_delta_RED`/`max_delta_RED_featureID` give you two complementary, N-PAS-aware effect-size views (net UTR-length shift, and the single biggest per-site relative-expression shift) to weigh alongside the significance columns — none of `perGeneQ`, `delta_wUTR`, or RED alone tells the whole story.

#### How `dominant_featureID` and `single_dominant_site` are determined

All of this happens once per gene inside `build_gene_apa_summary()`, on the full set of that gene's tested PAS.

**Step 1 — pick the dominant PAS.** `meanUsage_All` is each PAS's mean PSI across *every* sample in both conditions pooled together (not per-condition). Whichever PAS has the single highest `meanUsage_All` becomes `dominant_featureID` — the PAS carrying the plurality of that gene's poly-A reads overall, i.e. its "default" cleavage site. This pick is made once from the pooled data and is not re-evaluated separately per condition.

**Step 2 — check whether the dominant PAS itself moved.**
```
dominant_delta_usage = meanUsage_Treatment[dominant] − meanUsage_Control[dominant]
dom_stable            = is.na(dominant_delta_usage)  OR  abs(dominant_delta_usage) < DOMINANT_STABLE_CUT   (default 0.05)
```
An `NA` delta (e.g. the dominant PAS had zero reads in one condition) counts as "stable" by convention — there's no evidence it moved, so it isn't allowed to itself trigger the flag.

**Step 3 — check whether the dominant PAS is actually dominant, in both conditions.**
```
dom_dominant = dominant_meanUsage_Control >= DOMINANT_USAGE_CUT (default 0.90)
               AND dominant_meanUsage_Treatment >= DOMINANT_USAGE_CUT
```
This is a usage-magnitude check, not a significance check — it asks nothing about `padj`/`log2FC` on any PAS, only whether one site carries the overwhelming majority of the gene's reads in *both* conditions.

**Step 4 — combine.**
```
single_dominant_site = dom_dominant  AND  dom_stable
```
Both halves must hold simultaneously. A gene where the top PAS holds 92% of usage in Control but only 70% in Treatment is *not* flagged — usage moved too much for "single dominant site" to still describe the Treatment condition, even though a single site was clearly dominant in Control alone.

**Why it matters:** for a gene like this, whatever the minor PAS are individually doing (statistically significant or not) is unlikely to be biologically meaningful — one site is running the show, in both conditions, the whole time. It's a quick way to separate "this gene functionally has one PAS" from genes where usage is genuinely split across two or more real, biologically comparable options. It is purely informational: no rows are dropped or reordered because of it.

**Worked example:** a gene has 3 PAS. PAS_B has the highest pooled `meanUsage_All`, and holds 95% of usage in Control and 96% in Treatment (`dom_dominant = TRUE`). `dominant_delta_usage = 0.96 − 0.95 = 0.01` → `dom_stable = TRUE` (< 0.05). Result: `single_dominant_site = TRUE` — this gene functionally behaves as if it only has one PAS; PAS_A and PAS_C are minor enough, in both conditions, that whatever they're individually doing is unlikely to matter.

---

### `pas_usage.csv` — per-PAS results, counts, and usage (PSI)

The main results table. One row per PAS site that passed filtering and was tested — DEXSeq's own results, raw and normalized counts, and mean PSI + SEM per condition, all in one place.

| Column | Meaning |
|---|---|
| `groupID` | Gene identifier (Gene Symbol if available, else Ensembl/RefSeq) |
| `featureID` | PAS site ID (from PolyA_DB) |
| `gene` | Same as groupID |
| `PAS_type` | Annotation category (e.g., `3'UTR`, `intron`) |
| `APA_direction` | `Increased`, `Decreased`, `Ambiguous`, or `NA` — the sign of *this PAS's own* usage change only (see "Usage direction" above) |
| `exonBaseMean` | Mean normalized read count at this PAS across all samples |
| `pvalue` | DEXSeq per-PAS p-value |
| `padj` | FDR-adjusted p-value (Benjamini-Hochberg) |
| `log2fold_Treatment_v_Control` | log2(Treatment fractional usage / Control fractional usage) |
| `countData.Control_rep1`, ... | Raw read counts per sample |
| `normCountData.Control_rep1`, ... | DEXSeq's size-factor-normalized counts, same sample, sitting right next to its raw counterpart |
| `genomicData.seqnames/start/strand` | Genomic coordinates of the PAS |
| `meanUsage_Control` | Mean PSI across control replicates (0 to 1; 0.5 = 50% of gene's reads) |
| `meanUsage_Treatment` | Mean PSI across treatment replicates |
| `seUsage_Control` | Standard error of PSI across control replicates |
| `seUsage_Treatment` | Standard error of PSI across treatment replicates |
| `RED_Control` / `RED_Treatment` | Relative Expression Difference for this PAS, per condition — see "How RED is computed" below |
| `delta_RED` | `RED_Treatment − RED_Control`. The actual treatment-effect comparison; large in magnitude (either sign) means this PAS's relative usage shifted a lot between conditions. |
| `thin_counts` | `TRUE` when this PAS's raw reads, summed across replicates, are under `THIN_COUNTS_CUT` (default 40) in **either** condition. A soft warning, not a filter: `MIN_TOTAL_READS` already excludes PAS below its (lower, gene-total) bar before DEXSeq ever tests them — this just flags PAS that cleared that hard cutoff but are still thin on their own. |

PSI (fractional usage) = reads at this PAS / total reads at all PAS of this gene, per sample, then averaged. This is the most biologically interpretable metric: a PAS with `meanUsage_Control = 0.25` accounts for 25% of that gene's poly-A reads in controls.

**How to identify significant APA events:** filter for `padj < 0.05` (or your threshold). Further filter by `APA_direction == "Increased"` or `"Decreased"` to focus on which PAS moved, or use `gene_summary.csv`'s `delta_wUTR` for the gene-wide lengthening/shortening call.

#### How RED is computed

RED ("Relative Expression Difference") is our own generalization of the RED metric from MAAPER — computed here because MAAPER's version only works for exactly two PAS per gene (a "proximal" and "distal" site); genes with 3+ PAS get silently collapsed down to two representative sites and lose information. This version works for any number of tested PAS.

For PAS `i` in a gene with `N` tested PAS, per condition:
```
n(i)   = sum of that PAS's normalized counts (normCountData.*) across that condition's replicates
m(i)   = sum of every OTHER PAS's normalized counts in that condition (gene total − n(i))
RED(i) = log2( N * (n(i) + RED_PSEUDOCOUNT) / (m(i) + RED_PSEUDOCOUNT) )
```
`RED_PSEUDOCOUNT` (default 0.5, a Haldane-Anscombe-style correction) is added to both `n(i)` and `m(i)` before dividing. Without it, a PAS at exactly 0% or 100% usage in a condition would produce `log2(0)` or a division by zero — but those on/off shifts are exactly the dramatic cases a gene-ranking rubric cares most about, so rather than lose them to `NA`, the pseudocount keeps them as large, finite, rankable values. This is computed from summed *normalized* counts (`normCountData.*`), not the `meanUsage_*` fractions above — the pseudocount only makes sense added on an absolute count-like scale, not a 0–1 fraction.

`RED(i)` on its own answers "how far is this PAS's usage from the null expectation of equal use of all `N` sites" (positive = used more than its 1/N share; negative = less). `delta_RED = RED_Treatment − RED_Control` is the actual between-condition comparison — the `1/N` term cancels exactly in that subtraction (same `N` for both conditions of the same gene), which is why only the single formula above is needed for both the per-condition value and the delta.

See `gene_summary.csv`'s `max_abs_delta_RED` for the gene-level rollup (the PAS with the single biggest RED shift) — the direct replacement for ranking genes by RED under the old MAAPER-based approach.

---

### `candidate_genes.csv` — the curated shortlist

Computed by `build_candidate_genes()`. This is a **filter, not a flag**: unlike every other table in this pipeline, genes that don't clear every criterion below simply don't appear here — this file is meant to be a small, curated shortlist, not a fully annotated universe.

It's the direct successor to the old MAAPER/RED-based candidate-selection rubric (Kathleen's process): for each gene, take its single biggest-moving PAS (the same "greatest change in percent usage of any PAS" idea the old rubric used) and require that PAS/gene to clear **all** of:

| Criterion | Threshold | Replaces (old rubric) |
|---|---|---|
| Usage change | `\|meanUsage_change\| > USAGE_CHANGE_CUT` (default 0.05) | Same idea, same default (0.05) |
| Statistical significance | `padj < PADJ_CUT` (default 0.05) | MAAPER's separate `pval.fisher.adj < 0.001` — now DEXSeq's own per-PAS `padj`, on the same PAS the usage-change criterion already picked out |
| Read depth | not `thin_counts` (≥40 reads) | MAAPER's "read count at pPAS/dPAS > 40" |
| 3′ UTR length shift | `\|delta_wUTR\| > WUTR_CHANGE_CUT` (default 50 bp) | Joel's manually-computed "mean change in 3' UTR length" |
| Not a single-dominant-site gene | `single_dominant_site == FALSE` | The old rubric had no equivalent — this is new, and catches the same "big relative change, tiny real impact" false positives that criterion was already trying to filter, generalized past 2 PAS |
| RED magnitude | `\|delta_RED\| >= RED_MAG_CUT` (default 1) | MAAPER's `\|REDu1\| >= 1` |

Every threshold is a top-level constant (see "Key thresholds" below) — none are hardcoded inside `build_candidate_genes()`. Columns in the output are a mix of the winning PAS's own values (`featureID`, `meanUsage_change`, `padj`, `delta_RED`) and gene-level context (`delta_wUTR`, `max_abs_delta_RED`/`max_delta_RED_featureID`, `perGeneQ`, `dir_consensus`, `dominant_featureID`), sorted by `padj` ascending.

**A note on thresholds:** the significance bar here reuses `PADJ_CUT` (0.05) — noticeably looser than the old rubric's `pval.fisher.adj < 0.001`. At `PADJ_CUT` left at its default, expect a few dozen genes on a typical run; tightening `PADJ_CUT` to `0.001` to match the original rubric's actual stringency typically drops that to a small handful. Similarly, `RED_MAG_CUT` (our analog of `\|REDu1\| >= 1`) is the single most consequential of these six thresholds — check how the candidate count changes with and without it before treating a given run's list as final.

---

## Quality control plots

All QC plots (and the VST matrices behind the PCA plots) are saved to `plots/qc/`.

### `library_sizes.png`
Bar chart of total PAS read counts per sample (before any filtering). Bars should be roughly similar height across samples. Large differences (>3×) may indicate a failed sample or a sample that needs closer attention during normalization. This is assessed before the analysis runs.

### `vst_global.csv` / `vst_normalized[.group].csv`
The full variance-stabilized (VST) count matrix — one row per PAS, one column per sample — exported alongside the PCA plots that are already computed from it (`vst_global.csv` is pre-DEXSeq, all samples; `vst_normalized[.group].csv` is post-DEXSeq-normalization, per analysis group). Useful for spot-checking a specific gene's overall expression change against its PAS-level usage change directly, or for any downstream analysis (clustering, custom plots) that wants normalized values rather than raw counts.

### `PCA_global.png` / `PCA_global_uncorrected.png`
Principal component analysis on all samples simultaneously, using variance-stabilized 3' UTR PAS counts. Points that cluster by condition (e.g., all Controls together) indicate that the experimental effect is larger than technical variation — a good sign. Points that cluster by batch or unexpected metadata suggest confounding. When `USE_RUV=FALSE` this file is named `PCA_global.png`; when `USE_RUV=TRUE` it is named `PCA_global_uncorrected.png` so the pre-correction structure is preserved for comparison.

### `PCA_global_RUV_corrected.png` *(USE_RUV=TRUE only)*
The same PCA after regressing the estimated W factors out of the VST matrix using `limma::removeBatchEffect`. This is **for visualization only** — it shows what the data look like after removing the batch signal, but does not affect the statistical model. Compare this to the uncorrected PCA to see how much the W factors explained.

### `RUV_W_factors.png` *(USE_RUV=TRUE only)*
Bar chart of the estimated RUV W factor(s) per sample, colored by condition. Samples that drove the batch factor will have large positive or negative bars. Ideally this factor should not be correlated with condition — if it is, the correction may over-remove true signal.

### `PCA_normalized[.group].png`
PCA on the post-normalization counts within each analysis group. Similar interpretation as above, but after DEXSeq's size factor normalization.

### `size_factors[.group].png`
Bar chart of the normalization factors assigned by DEXSeq to each sample. All bars should be close to 1.0 (the dashed reference line). A sample with a very high or low size factor (>2× away from 1) may be an outlier.

### `dispersion[.group].png`
This plot shows the estimated variance of each PAS's counts across samples. DEXSeq fits a curve (shown in red) through the cloud of per-PAS dispersion estimates (black dots) and uses that curve to smooth individual estimates toward the genome-wide trend — a technique called shrinkage. Well-behaved data has a tight cloud around the fitted curve with a mean that decreases at higher mean counts.

### `pvalue_hist[.group].png`
Histogram of raw (unadjusted) p-values. The dashed red line shows what a uniform distribution would look like. A well-powered experiment with true effects shows a spike of small p-values near 0 (real signal) and a roughly flat distribution from 0.1–1.0 (null PAS). A completely flat or anti-conservative (U-shaped) histogram suggests model misspecification.

---

## Diagnostic plots

Saved to `plots/qc/`.

### `MA[.group].png`
Plots log2 fold change (Y) versus mean count (X) for every tested PAS. Significant hits (padj < threshold) are shown in red. A symmetric cloud centered at 0 with random scatter indicates no systematic bias. Up/down counts are annotated in the corners.

### `volcano[.group].png`
Plots log2 fold change (X) versus −log10(padj) (Y). Points in the upper corners are the most significant and most changed PAS. The top hits by significance are labeled with gene names. Color categories:
- **Red** — significant AND large fold change (`padj < threshold` AND `|log2FC| >= LFC_CUT`)
- **Gold** — significant but small fold change
- **Gray** — not significant

---

## APA genome-map figures (`plots/apa_genome/`, `plots/apa_zoom/`)

Two PNGs per gene in `GENES_OF_INTEREST`, per analysis group — both require `ANNO_GTF` to be set to a real GTF and `rtracklayer`/`ggtranscript` installed (skipped with a warning otherwise; see Software requirements). These use **real genomic coordinates and real exon/intron structure** parsed from that GTF (`gtf_exons`, built once near the top of the script) — PolyA_DB alone has no gene structure, only flat PAS positions, so nothing before this could show actual isoform geometry. The two figures answer different questions and are deliberately not the same design:

- **`apa_genome/{GENE}[.group].png`** — structure only, no usage/expression. The whole gene at real scale on top, and a `ggforce::facet_zoom()` inset of just the terminal-exon/PAS region below, connected by an automatic shaded funnel — tightly clustered PAS (e.g. 5 within ~1kb) are illegible at whole-gene scale, so the inset exists purely to make that count legible. Every annotated RefSeq transcript is its own row. Only an isoform's **terminal exon** can be colored **Alternative** (orange — this isoform's terminal exon has a different start/end than every other isoform overlapping that region); everything else, including genuinely alternatively-spliced internal exons elsewhere in the gene body, is drawn plain **Constitutive** grey regardless of how much it varies between isoforms. This is deliberate: this project cares about AS+APA *interaction* (a different terminal exon is simultaneously a splicing choice and a polyadenylation-site choice), not alternative splicing on its own — highlighting internal cassette-exon variation unrelated to which PAS a transcript ends at would just be noise for this question. Isoforms at a completely different, unrelated locus (e.g. a retained-intron transcript far upstream) are still shown as a row in both panels — in the zoomed one it's simply empty, since that isoform has no data there, which is accurate rather than a filtering choice.
- **`apa_zoom/{GENE}[.group].png`** — gene-level DGE and PAS usage, cropped tightly to just the terminal exon(s) near the detected PAS (unrelated distant isoforms excluded here — not because they're wrong, but because they're not part of this specific APA story and would blow up the crop window). A top row of two side-by-side panels, then two full-width panels stacked below on a shared x-axis (via `patchwork`):
  - **Gene log2FC** (top left) — a single vertical bar showing this gene's log2 fold change from `build_gene_dge()` (real DESeq2 test, gene-summed counts, independent of the DEXSeq exon-usage test below): up for a positive log2FC, down for a negative one, colored by direction and significance (red = up, blue = down, grey = not significant at `PADJ_CUT`). The log2FC/padj label sits beside the bar at a fixed position (the zero line) so it never depends on the bar's sign or length.
  - **Expression** (top right) — bars showing mean DESeq2 size-factor-normalized gene-level count per condition (the same normalized values as `normCount_Control`/`normCount_Treatment` in `gene_dge.csv`), colored by condition, with ±1 SEM error bars and every individual replicate plotted as a jittered dot on top — the same "don't just show the mean" treatment the PAS-usage panel already gets below. Sits next to the log2FC bar it's derived from, so "how much did overall gene expression change" and "in what direction, how confidently" read together in one glance. Together these two panels answer "is the gene as a whole up or down, and by how much" — a question the per-PAS panels below can't answer on their own, since if a gene's total expression differs between conditions, every one of its PAS shifts together in raw count terms even when relative usage hasn't changed at all.
  - **PAS usage** (middle, full width, 0–100%) — Control and Experimental overlaid as two colored connected lines (direct visual comparison, no per-point labels to parse), with every replicate plotted as a small jittered point behind the mean and ±1 SEM error bars, so you can see whether an apparent crossing pattern between conditions is a robust shift or within replicate noise. Each PAS's own DEXSeq `padj` is marked with 1–3 stars stacked vertically along its guide line, anchored to the bottom of the panel at a fixed position regardless of where that PAS's data sits (`*` < 0.05, `**` < 0.01, `***` < 0.001; nothing for non-significant PAS) — the same significance you'd read off `pas_usage.csv` directly, without leaving the plot. The Condition legend sits *inside* the panel (top-right, on a translucent background) rather than in a separate margin, to avoid burning a whole column of whitespace on it.
  - **Gene model** (bottom, full width) — same as before, but sized more modestly: its height grows with isoform count but is capped, so a gene with many isoforms doesn't crowd out the DGE/expression/usage panels above, which carry the more decision-relevant information. **Where a PAS falls beyond a transcript's annotated terminal-exon boundary, the exon is drawn extended into the intron with a dashed outline** — RefSeq's model stops at the annotated end, but a real detected cleavage site downstream means the transcript actually extends further than annotated; the dashed region marks exactly that gap.

**Minus-strand genes read left to right too.** A gene on the `-` strand is transcribed from a higher genomic coordinate toward a lower one, so plotting raw ascending coordinates left-to-right would make it read backwards. `apa_zoom` reverses the x-axis for `-` strand genes (a display-only flip — none of the underlying position math changes) so every gene, regardless of strand, reads 5'→3' left to right; the x-axis label says so explicitly, and the intron arrows still point the correct (now-rightward) transcription direction after the flip.

Across both figures, the dotted vertical guide lines marking each detected PAS are drawn dark and reasonably thick (`grey30`, linewidth 0.5) so they stay legible against the gene model and usage lines rather than blending into the plot background.

**Not implemented: full IsoformSwitchAnalyzeR integration.** IsoformSwitchAnalyzeR was considered for the ORF/NMD-consequence analysis needed to properly answer "does this APA event also change the coding sequence" (see the stop-codon section below) — it's the right tool for that question, but it expects transcript-level abundance quantification (e.g. from Salmon/kallisto), not just a GTF, which this PAS-count-based pipeline doesn't produce or consume. Integrating it would mean adding a parallel transcript-quantification workflow, which is out of scope for this script. If ORF/NMD consequences matter for a specific gene, run IsoformSwitchAnalyzeR as a separate analysis alongside this one rather than through `DEXseq.R`.

---

## Interpreting the results: a guide for collaborators

**What does it mean when a gene shows differential APA?**

It means cells under the two conditions preferentially choose different poly-A sites for that gene. This changes the length of the mRNA's 3' UTR without necessarily changing how much of the protein is made (though it can do that too).

A longer 3' UTR exposes more regulatory sequences — microRNA binding sites, RNA-binding protein sites, adenylation signals — that can reduce mRNA stability or protein output. A shorter 3' UTR "hides" those sites. APA is therefore a mechanism for fine-tuning gene regulation without changing the gene's transcription rate.

**Which genes should I pay attention to?**

The gene-summary table (`gene_summary.csv`) ranked by `perGeneQ` is the best place to start. Genes with:
- Low `perGeneQ` (close to 0) — strong statistical evidence for differential APA
- `dir_consensus = "Increased_only"` or `"Decreased_only"` — clean, directional changes
- High `n_sig` — multiple PAS in the gene are changing, not just one

**What does a significant padj mean vs. a large ΔPSI?**

`padj` tells you how confident we are that the change is real (statistical significance). ΔPSI tells you how large the change is (biological effect size). A small ΔPSI can be highly significant if we have many replicates. A large ΔPSI may not reach significance if the data are noisy or replicates are few. Both matter: ideally you want both low padj and large |ΔPSI|.

**What is "perGeneQ"?**

The per-gene q-value (Lancaster method) combines the evidence from all PAS sites within a gene into a single gene-level significance score. It is more conservative than just taking the minimum padj, but it accounts for the fact that testing multiple PAS per gene inflates the chance of at least one false positive. Use `perGeneQ` for gene-level prioritization.

---

## Running for different experimental designs

### Simple two-condition (no grouping)
Leave `GROUPING_VAR <- NULL`. All samples in the design file are analyzed together.

### Multiple timepoints (or batches as separate analyses)
Set `GROUPING_VAR <- "timepoint"` (or whatever the column is called in your design file). The script automatically detects all levels in that column and runs a separate DEXSeq analysis for each — for example, one for day 3 samples only, and one for day 5 samples only. Output files are suffixed accordingly (`.3d`, `.5d`).

### Adding a batch covariate (blocking factor)
DEXSeq's `~ sample + exon + condition:exon` model absorbs all per-sample gene-level variation through the `sample` term, which handles differences in overall expression between samples. However, it does **not** correct for batch effects on *relative* PAS usage (PSI). If replicates cluster by batch rather than condition in the PCA, use RUVseq.

### RUVseq batch correction (`USE_RUV <- TRUE`)
Set `USE_RUV <- TRUE` when your PCA shows replicates clustering by something other than condition — for example, if one replicate from each condition groups separately from the others, suggesting a technical batch effect.

RUVseq uses the within-condition replicate structure to estimate "unwanted variation" W factors — sources of variability that are shared within replicates of the same condition and therefore cannot be biological. `RUV_K` controls how many such factors to estimate (1 is almost always sufficient; try 2 if the first factor does not resolve the outlier clustering).

When `USE_RUV=TRUE` the DEXSeq design becomes:
```
~ sample + exon + W_1:exon + condition:exon
```
The `W_1:exon` term tests whether the batch factor changes relative PAS usage fractions within genes — this is the dimension that `sample` alone cannot capture. The `condition:exon` term retains the same interpretation as before.

**Important:** RUVseq corrects for unwanted variation that exists consistently across replicates. It cannot rescue experiments where replicates are fundamentally incomparable. Always compare the uncorrected and corrected PCA plots before trusting the RUV-corrected results.

---

## Caching (`USE_CACHE <- TRUE`)

The DEXSeq and DESeq2 model fits (`run_dexseq_group()`, `build_gene_dge()`) are by far the slowest part of this script — full GLM fits across every tested PAS. Most re-runs are really about iterating on a plot or `GENES_OF_INTEREST`, not changing the underlying data or modeling settings, so refitting from scratch every time is wasted time.

With `USE_CACHE <- TRUE` (the default), each group's fit is saved to `{OUT_BASE}/cache/` alongside a small manifest (mtimes of the count/annotation/design files, the exact sample list for that group, and the modeling settings that affect the fit: `CTRL_LABEL`, `TRTMT_LABEL`, `PAS_TYPE_REGEX`, `MIN_TOTAL_READS`, `MIN_PER_CONDITION`, `USE_RUV`, `RUV_K`). On the next run, if that manifest is unchanged, the cached fit is loaded instead of recomputed — the console prints `[cache] Reusing cached result for '...'` when this happens, or `[cache] No valid cache for '...' -- computing.` when it doesn't. Editing anything downstream of the fit (plot code, `GENES_OF_INTEREST`, color palettes, thresholds like `PADJ_CUT`/`LFC_CUT` that only affect summarization/plotting) leaves the cache valid. Changing the count matrix, annotation, design file, or any of the modeling settings above automatically invalidates it — no manual cleanup needed.

Set `USE_CACHE <- FALSE` to always refit from scratch (e.g. while actively debugging `run_dexseq_group()`/`build_gene_dge()` themselves, where you don't want a stale-but-manifest-matching cache masking a code change). The `cache/` directory is otherwise safe to delete any time — it will simply be rebuilt on the next run.

---

## Key thresholds and what they control

| Parameter | Default | Effect |
|---|---|---|
| `PAS_TYPE_REGEX` | `^3'UTR` | Which PAS types to include. Change to `""` to include all types. |
| `MIN_TOTAL_READS` | 10 | Minimum total reads at a PAS in the group. Lower = more PAS tested, more noise. |
| `MIN_PER_CONDITION` | 1 | Minimum samples with >0 reads per condition. Ensures both conditions have detectable signal. |
| `PADJ_CUT` | 0.05 | FDR threshold. Lower = fewer but more confident hits. |
| `LFC_CUT` | 0.5 | log2 fold change threshold for "big" hits in the gene summary and volcano plot coloring. |
| `THIN_COUNTS_CUT` | 40 | Per-PAS, per-condition warning threshold for `thin_counts` in `pas_usage.csv`. A flag, not a filter — `MIN_TOTAL_READS` above is what actually excludes PAS. |
| `DOMINANT_USAGE_CUT` | 0.90 | Minimum usage share (in **both** conditions) a gene's top PAS needs to count as "the" dominant site for `single_dominant_site` in `gene_summary.csv`. |
| `DOMINANT_STABLE_CUT` | 0.05 | Maximum allowed shift in the dominant PAS's own usage between conditions to still call it "stable" for `single_dominant_site`. |
| `RED_PSEUDOCOUNT` | 0.5 | Haldane-Anscombe-style constant added to both `n(i)` and `m(i)` when computing `RED_Control`/`RED_Treatment`/`delta_RED`, so a PAS at 0%/100% usage still yields a large finite value instead of `NA`/`-Inf`. |
| `USAGE_CHANGE_CUT` | 0.05 | Minimum `\|meanUsage_change\|` for a gene's top-moving PAS to qualify for `candidate_genes.csv`. |
| `WUTR_CHANGE_CUT` | 50 | Minimum `\|delta_wUTR\|` (bp) for `candidate_genes.csv`. |
| `RED_MAG_CUT` | 1 | Minimum `\|delta_RED\|` (log2 scale) for `candidate_genes.csv` — our analog of the old MAAPER rubric's `\|REDu1\| >= 1`. |
| `NTOP_PCA` | 2000 | Number of most-variable PAS to use for PCA. 500–5000 is typical. |

---

## Notes on the method

**Why DEXSeq and not a simpler test?**

DEXSeq's key advantage is that it explicitly models the within-gene relative usage using a negative binomial model, which handles the overdispersion typical of count data. The `sample` term absorbs all inter-sample variability, making the test highly specific to isoform usage changes rather than overall expression changes. Simpler approaches (e.g., Fisher's exact test on counts) do not account for between-replicate variance and over-report false positives.

**Why PSI and not log2FC for the bar plots?**

The DEXSeq log2 fold change is a ratio-of-ratios on the count scale (internally: Treatment PSI / Control PSI, log2-transformed). While it is useful for volcano/MA plots, it is not intuitively interpretable. ΔPSI (simply Treatment PSI − Control PSI) is expressed in the same units as the bars themselves and is directly relatable to the plot — a ΔPSI of +0.15 means the bar got 15 percentage points taller.
