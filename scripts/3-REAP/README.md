# 3-REAP: Paired-Read Alternative Polyadenylation Pipeline

**Script:** `scripts/3-REAP_workflow.sh`
**Last updated:** 2026-06-25

---

## What this pipeline does — and why

### The biology

Every gene has a "stop signal" at the end of its messenger RNA (mRNA) — a short sequence called a **poly-adenylation site (PAS)** where the cell stops copying the gene and adds a protective poly-A tail. Most genes have **more than one PAS**, which means the same gene can produce mRNA molecules of different lengths. This is called **alternative polyadenylation (APA)**.

The choice of which PAS to use changes the length of the 3' untranslated region (3' UTR) — the tail end of the mRNA that doesn't code for protein. The 3' UTR contains regulatory elements: binding sites for microRNAs and RNA-binding proteins that control how stable the mRNA is, how much protein gets made, and where in the cell the mRNA goes. **APA is therefore a major post-transcriptional regulatory mechanism**, and dysregulation of APA is associated with cancer, neurodegeneration, and development.

This pipeline analyzes sequencing data designed specifically to capture APA. The sequencing strategy uses **paired-end reads** where one read (R2) is anchored near the poly-A tail of mRNAs, allowing us to pinpoint exactly which PAS was used in each sequenced molecule.

### What we measure

From sequencing data → we count how many reads ended at each known PAS site, for each sample. The output is a **PAS count matrix**: rows are PAS sites, columns are samples, values are read counts. This is the direct equivalent of a gene count matrix in standard RNA-seq, but at the resolution of individual poly-A sites within each gene.

That count matrix is the input to the downstream statistical analysis (see `README.DEXseq.md`).

---

## Quick-start: how to run the pipeline

### 1. Edit the configuration block

Open `3-REAP_workflow.sh` in a text editor. The only section you need to change is at the top, clearly marked `USER CONFIGURATION — EDIT THESE`:

```
fastq_path   — directory containing raw FASTQ files from the sequencer
work_path    — where all output files will be written
sample_file  — path to the sample manifest (see format below)
threads_per_job / trim_threads / max_jobs  — parallelism (scale to your server)
lap_dist     — assignment window in nucleotides (default 24; don't change without reason)
```

Everything else (output directory names, reference file paths) is derived automatically.

### 2. Prepare the sample manifest

`data/sample.txt` is a tab-separated file with three columns, one row per sample, no header:

```
sample_name    /absolute/path/to/sample_R1.fastq.gz    /absolute/path/to/sample_R2.fastq.gz
```

Sample names must not contain spaces. R2 path can be left blank if it follows the standard naming convention `_R1_` → `_R2_` (the script derives it automatically).

### 3. Check that reference files exist

The pipeline uses:
- `data/ref/star_index_hg38/` — STAR genome index (human hg38)
- `data/PolyA_DB_v4.1/hg38.PAS.main.tsv` — PolyA_DB v4.1 reference PAS database
- `data/ref/hg38.chrom.sizes` — chromosome size file

### 4. Activate conda env with appropriate dependancies

Required tools: 
- umi_tools 1.1.5
- cutadapt 4.9
- bbmap 39.28
- STAR 2.7.11b
- samtools 1.22
- bedtools 2.31.1

### 5. Run

```bash
bash 3-REAP_workflow.sh 2>&1 | tee pipeline_run_$(date +%Y%m%d_%H%M).log
```

The `tee` command saves a timestamped log alongside the script. The pipeline prints progress at each step and skips any sample/step where output already exists, so it is safe to re-run after a partial failure.

---

## Software dependencies

| Tool | Purpose |
|---|---|
| `cutadapt` | Adapter and poly-A trimming |
| `repair.sh` (BBTools) | Re-pair FASTQ read pairs after independent trimming |
| `STAR` | Splice-aware genome alignment |
| `samtools` | BAM indexing, filtering |
| `bedtools` | BAM to BED conversion |
| `umi_tools` | UMI extraction/deduplication (optional, currently disabled) |
| `genomeCoverageBed` (bedtools) | Bedgraph generation |
| `bedGraphToBigWig` (UCSC) | Convert bedgraph to bigWig |
| `Rscript` | LAP/PAS assignment and count table merging |
| Python 3 | Custom 5'-T trimming script |

All tools must be on `$PATH` before running.

---

## Pipeline steps in detail

### Step 1 — UMI extraction (disabled)
If the library was prepared with unique molecular identifiers (UMIs), this step extracts them from the reads and encodes them in the read name. It is currently commented out because the cmoore_020 library does not use UMIs.

### Step 2 — Adapter and poly-A trimming

Sequencing adapters (synthetic sequences from the library preparation kit) and internal poly-A/poly-T stretches must be removed before alignment. R1 and R2 are treated differently because they serve different roles in this protocol.

**R1** goes through four passes of `cutadapt`:
1. Strip poly-A runs ≥20 nt and poly-G quality artifacts (NextSeq two-color chemistry artifact)
2. Remove Illumina TruSeq Read 1 adapter
3. Remove any remaining internal poly-A
4. Discard reads that still contain the adapter after the above (internal priming artifacts)

**R2** goes through two passes, followed by a custom Python script:
1. Remove Illumina TruSeq Read 2 adapter
2. Discard R2 reads that still contain the adapter (artifact filter)
3. `modules/trim_N-5T_fq.py` — strips 5'-terminal T runs from R2. R2 is anchored to the poly-A tail; the 5' end of R2 often carries trailing T nucleotides from the tail itself, which must be removed before alignment.

Minimum read length after trimming: 20 nt (R1), 18 nt (R2). Shorter reads are discarded.

Samples run in parallel (up to `max_jobs` at a time). Already-trimmed samples are skipped.

**Output:** `01_trim/{sample}/{sample}_R1.fastq.gz` and `{sample}_R2_trimAdapt.5Ttrimmed.fastq.gz`

### Step 3 — Read-pair repair and genome alignment

After independent trimming, some read pairs may have become unpaired (one read was discarded, the other wasn't). `repair.sh` (BBTools) re-synchronizes pairs before alignment.

Alignment is performed with **STAR** — a splice-aware aligner that can handle reads spanning exon-exon junctions. Key parameters:
- `--outFilterMultimapNmax 1` — discard reads that map equally well to multiple locations (multimappers)
- `--outFilterScoreMinOverLread 0.2` / `--outFilterMatchNminOverLread 0.2` — allow short R2 reads (after 5'-T trimming) to align successfully
- Output is a coordinate-sorted BAM file

After all samples align, alignment statistics are compiled into a single summary CSV (`02_align/stats_STAR_summary.csv`).

**Output:** `02_align/{sample}/star/Aligned.sortedByCoord.out.bam` and `stats_STAR_summary.csv`

### Step 4 — Deduplication

PCR amplification during library preparation can create duplicate read pairs from the same original molecule, inflating read counts. This step is designed to remove them.

**Currently:** UMI deduplication is disabled (this library does not have UMIs). Instead, the aligned BAM is symlinked to `03_dedup/` as a placeholder. To enable UMI deduplication with `umi_tools dedup`, replace the `ln -s` command in the script with the commented-out `umi_tools` command.

**Output:** `03_dedup/{sample}/{sample}.dedup.bam` (symlink or deduplicated BAM)

### Step 5 — BAM to BED (R2 reads only)

R2 reads are the ones whose 5' end marks the poly-A cleavage site. We extract only R2 reads from the BAM (`-f 128` = read2 flag) and convert them to BED format, which is what the LAP/PAS assignment script expects.

A per-sample and total read count is printed. If all BED files are empty, the pipeline exits with an error because that indicates R2 trimming or alignment failed (commonly caused by missing R2 paths in the sample file).

**Output:** `04_bed/{sample}.dedup_R2.sorted.bed`, `stats_deduped_reads.txt`

### Step 6 — LAP definition and PAS assignment

This is the core analytical step. The R script `modules/genome_LAPandPAS_define.R` does two things:

1. **Define LAPs (Last Aligned Positions):** The 5' end of each R2 read marks where on the genome the poly-A cleavage occurred. Groups of R2 reads whose 5' ends cluster together are summarized as a **LAP** — a candidate cleavage site.

2. **Assign LAPs to known PAS:** Each LAP is matched to the nearest entry in the **PolyA_DB v4.1** reference database (a curated catalog of validated human poly-A sites). The assignment window is ±`lap_dist` nucleotides (default: ±24 nt). LAPs with no nearby reference PAS are discarded; LAPs within range are assigned to that PAS and their read counts are recorded.

The result per sample is a table: `{sample}_cluster.all.reads.csv`, where each row is a reference PAS and the value is the number of supporting R2 reads in that sample.

These per-sample tables are then merged by `modules/combine_all_sample_PAS_count_tables.R` into a single matrix with all PAS as rows and all samples as columns.

**Key parameter:** `lap_dist` (default 24). Increasing it assigns more reads but risks false assignments; decreasing it is more conservative. This value is set in the USER CONFIGURATION block.

**Output:**
- `05_LAP_PAS/{sample}_cluster.all.reads.csv` — per-sample PAS count table
- `05_LAP_PAS/PAS_quant/cluster.all.reads.csv` — **combined PAS count matrix (primary output)**
- `05_LAP_PAS/PAS_quant/stats_PAS.csv` — combined per-sample detection statistics

### Step 7 — bigWig track generation

Genome browser tracks are generated so reads can be visualized in tools like **IGV** or **UCSC Genome Browser**. Tracks are CPM-normalized (counts per million) so samples are comparable regardless of sequencing depth.

Four track types are produced per sample, each useful for a different visual question:

| Track type | What it shows |
|---|---|
| `PASS_bw_LAP{dist}` | All R2 reads successfully assigned to a PAS within the window |
| `PASS_bw_LAP{dist}_position` | LAP positions only (peak centers, no width) |
| `PASS_bw_LAP{dist}_PAS` | Reads collapsed to the assigned PAS coordinate |
| `mapped_read` | All mapped R2 reads (no assignment filtering) |

Each track type has three files per sample:
- `plus.{sample}.*.bw` — plus-strand signal (positive values)
- `minus.{sample}.*.bw` — minus-strand signal (positive values, for tools that display both strands separately)
- `mi.{sample}.*.bw` — minus-strand signal (negative values, for single-track IGV display)

**Output:** `06_bigwig/{track_type}/`

---

## Output directory summary

```
{work_path}/
├── 01_trim/            Trimmed FASTQs (intermediate; large; can delete after step 3)
├── 02_align/
│   ├── {sample}/star/  Aligned BAMs
│   └── stats_STAR_summary.csv
├── 03_dedup/           Deduplicated BAMs (or symlinks)
├── 04_bed/             R2 read BED files
├── 05_LAP_PAS/
│   ├── {sample}_cluster.all.reads.csv   Per-sample PAS tables
│   └── PAS_quant/
│       ├── cluster.all.reads.csv        ← PRIMARY OUTPUT: PAS count matrix
│       └── stats_PAS.csv                Per-sample PAS detection statistics
└── 06_bigwig/
    ├── PASS_bw_LAP{dist}/               Genome browser tracks
    ├── PASS_bw_LAP{dist}_positon/
    ├── PASS_bw_LAP{dist}_PAS/
    └── mapped_read/
```

### The primary output: `cluster.all.reads.csv`

This is what matters most. It is a matrix where:
- **Rows** = poly-A sites (PAS IDs from PolyA_DB v4.1 in `chr:strand:position` format)
- **Columns** = samples
- **Values** = integer read counts (how many R2 reads from that sample ended at that PAS)

This file feeds directly into the DEXSeq differential APA analysis. See `README.DEXseq.md`.

---

## Reference database: PolyA_DB v4.1

PolyA_DB is a curated database of human poly-adenylation sites supported by experimental evidence from multiple sequencing datasets. Version 4.1 (used here) encodes the genomic coordinates of each PAS in its ID field as `chr:strand:position`.

The database file is at: `data/PolyA_DB_v4.1/hg38.PAS.main.tsv`

Each row is a validated PAS site with annotations including which gene it belongs to, whether it is intronic or 3' UTR, and its FANTOM/RefSeq/Ensembl identifiers.

---

## Troubleshooting

**All R2 BED files are empty (pipeline exits at step 5)**
The most common cause is that `sample.txt` is missing the R2 path column, or the R2 file paths are wrong. Check that all three columns are present and that R2 FASTQ files exist.

**STAR alignment rate is very low (<40% uniquely mapped)**
Check that `star_index_path` points to the correct genome build. Also check trimming: very aggressive trimming can make reads too short to align.

**LAP step produces very few PAS assignments**
Try increasing `lap_dist` (e.g., 30–50). Also confirm that the chromosome naming in your BAM (chr1 vs 1) matches the reference PAS file and STAR index.

**A sample is missing from the combined matrix**
The combine script expects a specific filename pattern (`_cluster.all.reads.csv`). Check that the sample's step 6 output exists and was moved to `PAS_quant/`.

---

## Module scripts (do not call directly)

| Script | Called by | Purpose |
|---|---|---|
| `modules/trim_N-5T_fq.py` | Step 2 | Strip 5'-terminal T runs from R2 |
| `modules/genome_LAPandPAS_define.R` | Step 6 | Define LAPs and assign to PolyA_DB |
| `modules/combine_all_sample_PAS_count_tables.R` | Step 6 | Merge per-sample count tables |
| `modules/norm_bedgraph.pl` | Step 7 | CPM normalize bedgraph files |
