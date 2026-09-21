#!/usr/bin/env bash
### 3-REAP: Paired-read Alternative Polyadenylation pipeline

#####################################################
###  USER CONFIGURATION — EDIT THESE            ###
#####################################################

# Auto-detect project root (one directory above scripts/3-REAP/)
PROJECT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"

# Paths to raw FASTQ files and pipeline outputs
fastq_path=/compbio/data/ClaireMoore/cmoore_019.admera/raw/       # directory containing raw FASTQ files
work_path=/compbio/data/ClaireMoore/cmoore_019.admera/processed/     # all pipeline outputs go here

# Sample manifest: tab-separated, one sample per line
#   <sample_name>  <R1_absolute_path>  <R2_absolute_path>
sample_file=${PROJECT_DIR}/data/sample.txt

# Parallelism
threads_per_job=8    # STAR alignment threads per sample
trim_threads=2       # cutadapt threads per invocation (piped chains multiply this)
max_jobs=4           # maximum samples processed in parallel

# Reference files
lap_dist=24          # LAP-to-PAS assignment window (±nt); propagates to Step 6 and Step 7 track names
star_index_path=${PROJECT_DIR}/data/ref/star_index_hg38
refPAS_file=${PROJECT_DIR}/data/PolyA_DB_v4.1/hg38.PAS.main.tsv
chromsizes=${PROJECT_DIR}/data/ref/hg38.chrom.sizes

#####################################################
###  DERIVED PATHS — DO NOT EDIT               ###
#####################################################
modules_path=${PROJECT_DIR}/scripts/3-REAP/modules

echo "=== 3-REAP pipeline ==="
echo "    Project : ${PROJECT_DIR}"
echo "    Output  : ${work_path}"
echo "    Samples : ${sample_file}"
echo "    PolyADB : ${refPAS_file}"
echo ""

# ─────────────────────────────────────────────────────────────────────────────
# Step 1. UMI extraction  (disabled — enable if UMI library prep was used)
# ─────────────────────────────────────────────────────────────────────────────
# while IFS=$'\t' read -r sample fastq_1 fastq_2; do
#     outpath=${work_path}/00_umi/${sample}; mkdir -p "${outpath}"
#     umi_tools extract --extract-method=string --bc-pattern X --bc-pattern2 NNNNNNNNNN \
#         -I "${fastq_1}" --read2-in="${fastq_2}" \
#         -L "${outpath}/extraction_log" \
#         -S "${outpath}/${sample}_R1.fastq.gz" \
#         --read2-out="${outpath}/${sample}_R2.fastq.gz"
# done < "${sample_file}"

# ─────────────────────────────────────────────────────────────────────────────
# Step 2. Adapter and polyA/T trimming  (R1 and R2 per sample, parallel)
# ─────────────────────────────────────────────────────────────────────────────
echo "[Step 2] Adapter / polyA trimming..."

job_count=0
while IFS=$'\t' read -r sample fastq_1 fastq_2; do
    # Strip Windows carriage returns; derive R2 from R1 if the column was left empty
    fastq_2="${fastq_2%$'\r'}"
    [[ -z "${fastq_2}" ]] && fastq_2="${fastq_1/_R1_/_R2_}"
    {
        outpath=${work_path}/01_trim/${sample}
        mkdir -p "${outpath}"

        # ── R1 ────────────────────────────────────────────────────────────────
        outfile_r1=${outpath}/${sample}_R1.fastq.gz
        if [[ -f "${outfile_r1}" ]]; then
            echo "  [R1 skip] ${sample}"
        else
            echo "  [R1] ${sample}"
            cutadapt --quiet -j "${trim_threads}" -m 20 -O 20 \
                -a "polyA=A{20}" -a "QUALITY=G{20}" -n 2 "${fastq_1}" | \
            cutadapt --quiet -j "${trim_threads}" -m 20 -O 3 --nextseq-trim=10 \
                -a "r1adapter=AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC;min_overlap=3;max_error_rate=0.100000" - | \
            cutadapt --quiet -j "${trim_threads}" -m 20 -O 3 \
                -a "r1polyA=A{18}" - | \
            cutadapt --quiet -j "${trim_threads}" -m 20 -O 20 \
                -g "r1adapter=AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC;min_overlap=20" \
                --discard-trimmed -o "${outfile_r1}" -
        fi

        # ── R2 ────────────────────────────────────────────────────────────────
        outfile_r2=${outpath}/${sample}_R2_trimAdapt.5Ttrimmed.fastq.gz
        if [[ -f "${outfile_r2}" ]]; then
            echo "  [R2 skip] ${sample}"
        else
            echo "  [R2] ${sample}"
            r2_adapt=${outpath}/${sample}_R2_trimAdapt.fastq.gz
            cutadapt --quiet -j "${trim_threads}" -m 18 -O 3 --nextseq-trim=10 \
                -a "rd2Adapter=AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT;min_overlap=3;max_error_rate=0.100000" "${fastq_2}" | \
            cutadapt --quiet -j "${trim_threads}" -m 18 -O 20 \
                -g "rd2Adapter=AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT;min_overlap=20" \
                --discard-trimmed -o "${r2_adapt}" -
            python "${modules_path}/trim_N-5T_fq.py" \
                --rawfastq "${r2_adapt}" --out_dir "${outpath}" --random_NT_len 0
            rm -f "${r2_adapt}"
        fi
    } &
    ((job_count++))
    if ((job_count >= max_jobs)); then
        wait -n
        ((job_count--))
    fi
done < "${sample_file}"
wait

# ─────────────────────────────────────────────────────────────────────────────
# Step 3. Repair read pairs and align with STAR
# ─────────────────────────────────────────────────────────────────────────────
echo "[Step 3] Repair pairs + STAR alignment..."

job_count=0
while IFS=$'\t' read -r sample fastq_1 fastq_2; do
    {
        trim_dir=${work_path}/01_trim/${sample}
        pair_dir=${work_path}/02_align/${sample}/repaired
        star_dir=${work_path}/02_align/${sample}/star
        star_bam=${star_dir}/Aligned.sortedByCoord.out.bam

        if [[ -f "${star_bam}" ]]; then
            echo "  [skip] ${sample}: already aligned"
        else
            echo "  [align] ${sample}"
            mkdir -p "${pair_dir}" "${star_dir}"

            repair.sh overwrite=t \
                in1="${trim_dir}/${sample}_R1.fastq.gz" \
                in2="${trim_dir}/${sample}_R2_trimAdapt.5Ttrimmed.fastq.gz" \
                out1="${pair_dir}/${sample}_R1.fastq.gz" \
                out2="${pair_dir}/${sample}_R2.fastq.gz" \
                outs="${pair_dir}/${sample}_unpaired.fastq.gz" \
                repair

            STAR --runThreadN "${threads_per_job}" \
                --readFilesCommand zcat \
                --genomeDir "${star_index_path}" \
                --readFilesIn "${pair_dir}/${sample}_R1.fastq.gz" "${pair_dir}/${sample}_R2.fastq.gz" \
                --outReadsUnmapped None \
                --outFilterType BySJout \
                --limitBAMsortRAM 4000000000 \
                --outSAMattributes NH HI NM MD \
                --outSAMtype BAM SortedByCoordinate \
                --outFileNamePrefix "${star_dir}/" \
                --outFilterMultimapNmax 1 \
                --outFilterScoreMinOverLread 0.2 \
                --outFilterMatchNminOverLread 0.2

            # Remove repaired FASTQs after confirmed successful alignment
            [[ -f "${star_bam}" ]] && rm -rf "${pair_dir}"
        fi
    } &
    ((job_count++))
    if ((job_count >= max_jobs)); then
        wait -n
        ((job_count--))
    fi
done < "${sample_file}"
wait

echo "[Step 3] Compiling alignment statistics..."
stats_out=${work_path}/02_align/stats_STAR_summary.csv
echo "Sample,Total Reads,Uniquely Mapped (%),Multi-Loci (%),Unmapped:Mismatches (%),Unmapped:TooShort (%),Unmapped:Other (%)" > "${stats_out}"
for logfile in "${work_path}"/02_align/*/star/*Log.final.out; do
    sample=$(basename "$(dirname "$(dirname "${logfile}")")")
    total_reads=$(grep    "Number of input reads"                       "${logfile}" | awk '{print $NF}')
    uniquely_mapped=$(grep "Uniquely mapped reads %"                    "${logfile}" | awk '{print $NF}')
    multiple_mapped=$(grep "% of reads mapped to multiple loci"         "${logfile}" | awk '{print $NF}')
    unmapped_mm=$(grep     "% of reads unmapped: too many mismatches"   "${logfile}" | awk '{print $NF}')
    unmapped_short=$(grep  "% of reads unmapped: too short"             "${logfile}" | awk '{print $NF}')
    unmapped_other=$(grep  "% of reads unmapped: other"                 "${logfile}" | awk '{print $NF}')
    echo "${sample},${total_reads},${uniquely_mapped},${multiple_mapped},${unmapped_mm},${unmapped_short},${unmapped_other}" >> "${stats_out}"
done

# ─────────────────────────────────────────────────────────────────────────────
# Step 4. PCR deduplication  (symlink placeholder — swap ln for umi_tools dedup to enable)
# ─────────────────────────────────────────────────────────────────────────────
echo "[Step 4] Deduplication..."

job_count=0
while IFS=$'\t' read -r sample fastq_1 fastq_2; do
    {
        star_bam=${work_path}/02_align/${sample}/star/Aligned.sortedByCoord.out.bam
        outpath=${work_path}/03_dedup/${sample}
        bam_out=${outpath}/${sample}.dedup.bam
        mkdir -p "${outpath}"

        if [[ ! -f "${star_bam}.bai" ]] || [[ "${star_bam}" -nt "${star_bam}.bai" ]]; then
            samtools index "${star_bam}"
        fi

        if [[ -e "${bam_out}" ]]; then
            echo "  [skip] ${sample}: dedup link exists"
        else
            ln -s "${star_bam}" "${bam_out}"
            # To enable UMI deduplication, replace the ln above with:
            # umi_tools dedup --paired -I "${star_bam}" -S "${bam_out}" \
            #     --multimapping-detection-method=NH \
            #     --output-stats="${outpath}/deduplicated.txt" \
            #     --log="${outpath}/deduplication.log"
        fi
    } &
    ((job_count++))
    if ((job_count >= max_jobs)); then
        wait -n
        ((job_count--))
    fi
done < "${sample_file}"
wait

# ─────────────────────────────────────────────────────────────────────────────
# Step 5. BAM → BED  (R2 reads only, streamed — no intermediate R2 BAM saved)
# ─────────────────────────────────────────────────────────────────────────────
# samtools view flags: -f 128 select read2; -F 4 exclude unmapped mates
echo "[Step 5] BAM → BED (R2 reads)..."

bed_outdir=${work_path}/04_bed
mkdir -p "${bed_outdir}"

job_count=0
while IFS=$'\t' read -r sample fastq_1 fastq_2; do
    {
        final_bed=${bed_outdir}/${sample}.dedup_R2.sorted.bed
        if [[ -f "${final_bed}" ]]; then
            echo "  [skip] ${sample}: BED exists"
        else
            echo "  [bed] ${sample}"
            samtools view -f 128 -F 4 -b \
                "${work_path}/03_dedup/${sample}/${sample}.dedup.bam" | \
            bedtools bamtobed -cigar -i stdin | \
            sort -k1,1 > "${final_bed}"
        fi
    } &
    ((job_count++))
    if ((job_count >= max_jobs)); then
        wait -n
        ((job_count--))
    fi
done < "${sample_file}"
wait

echo "  Read counts per sample:"
wc -l "${bed_outdir}"/*.sorted.bed | tee "${bed_outdir}/stats_deduped_reads.txt"

total_bed_reads=$(awk '/total/{print $1}' "${bed_outdir}/stats_deduped_reads.txt")
if [[ "${total_bed_reads:-0}" -eq 0 ]]; then
    echo "ERROR: All R2 BED files are empty. R2 trimming likely failed — check that sample.txt has valid R2 paths and that the raw R2 FASTQs exist." >&2
    exit 1
fi

# ─────────────────────────────────────────────────────────────────────────────
# Step 6. Define LAPs and assign to PAS  (PolyA_DB v4.1; parallel per sample)
# ─────────────────────────────────────────────────────────────────────────────
# -misM 2  : allow up to 2 nt soft-clipping at read ends
# -dist N  : assign LAPs to PAS within ±N nt (set by lap_dist above)
echo "[Step 6] LAP/PAS analysis (dist=${lap_dist} nt)..."

quant_dir=${work_path}/05_LAP_PAS/PAS_quant
final_combined=${quant_dir}/cluster.all.reads.csv

if [[ -f "${final_combined}" ]]; then
    echo "  Output already exists, skipping"
else
    lap_dir=${work_path}/05_LAP_PAS
    mkdir -p "${lap_dir}"

    job_count=0
    while IFS=$'\t' read -r sample fastq_1 fastq_2; do
        {
            echo "  [LAP] ${sample}"
            Rscript "${modules_path}/genome_LAPandPAS_define.R" \
                -bedLAP "${bed_outdir}/${sample}.dedup_R2.sorted.bed" \
                -out    "${lap_dir}/${sample}" \
                -refPAS "${refPAS_file}" \
                -misM 2 -dist "${lap_dist}"
        } &
        ((job_count++))
        if ((job_count >= max_jobs)); then
            wait -n
            ((job_count--))
        fi
    done < "${sample_file}"
    wait

    mkdir -p "${quant_dir}"
    mv "${lap_dir}"/*_cluster.all.reads.csv "${quant_dir}/"
    Rscript "${modules_path}/combine_all_sample_PAS_count_tables.R" \
        -csv "${quant_dir}" -out "${quant_dir}/cluster.all.reads.csv" \
        --pattern "_cluster\\.all\\.reads\\.csv$"
    Rscript "${modules_path}/combine_all_sample_PAS_count_tables.R" \
        -csv "${lap_dir}" -out "${quant_dir}/stats_PAS.csv" \
        --pattern ".*_stats\\.csv$"
fi

# ─────────────────────────────────────────────────────────────────────────────
# Step 7. bigWig generation  (parallelized across sample × track type)
# ─────────────────────────────────────────────────────────────────────────────
echo "[Step 7] Generating bigWig tracks..."

lap_dir=${work_path}/05_LAP_PAS

job_count=0
for type in "PASS_bw_LAP${lap_dist}" "PASS_bw_LAP${lap_dist}_positon" "PASS_bw_LAP${lap_dist}_PAS" "mapped_read"; do
    bw_outdir=${work_path}/06_bigwig/${type}
    mkdir -p "${bw_outdir}"
    while IFS=$'\t' read -r sample fastq_1 fastq_2; do
        {
            echo "  [bw] ${type} / ${sample}"
            input_bed=${lap_dir}/${sample}_${type}.bed
            total_reads=$(wc -l < "${input_bed}")

            # Sort to a temp file; genomeCoverageBed requires sorted BED and is called twice
            sorted_bed=$(mktemp)
            sort -k1,1 "${input_bed}" > "${sorted_bed}"

            genomeCoverageBed -bg -split -i "${sorted_bed}" -strand '+' -g "${chromsizes}" \
                > "${bw_outdir}/${sample}.plus.bedgraph"
            genomeCoverageBed -bg -split -i "${sorted_bed}" -strand '-' -g "${chromsizes}" \
                > "${bw_outdir}/${sample}.minus.bedgraph"
            rm -f "${sorted_bed}"

            "${modules_path}/norm_bedgraph.pl" -t "${total_reads}" -i "${bw_outdir}/${sample}.plus.bedgraph"
            "${modules_path}/norm_bedgraph.pl" -t "${total_reads}" -i "${bw_outdir}/${sample}.minus.bedgraph"

            # Give the minus strand negative values for IGV display
            awk -v FS="\t" -v OFS="\t" 'NR>1 {print $1,$2,$3,-$4}' \
                "${bw_outdir}/${sample}.minus.bedgraph.normalized" \
                > "${bw_outdir}/${sample}.minus.bedgraph.normalized1"

            bedGraphToBigWig "${bw_outdir}/${sample}.plus.bedgraph.normalized"   "${chromsizes}" \
                "${bw_outdir}/plus.${sample}.${type}.bw"
            bedGraphToBigWig "${bw_outdir}/${sample}.minus.bedgraph.normalized"  "${chromsizes}" \
                "${bw_outdir}/minus.${sample}.${type}.bw"
            bedGraphToBigWig "${bw_outdir}/${sample}.minus.bedgraph.normalized1" "${chromsizes}" \
                "${bw_outdir}/mi.${sample}.${type}.bw"

            rm -f "${bw_outdir}/${sample}".{plus,minus}.bedgraph \
                  "${bw_outdir}/${sample}".{plus,minus}.bedgraph.normalized \
                  "${bw_outdir}/${sample}.minus.bedgraph.normalized1"
        } &
        ((job_count++))
        if ((job_count >= max_jobs)); then
            wait -n
            ((job_count--))
        fi
    done < "${sample_file}"
done
wait
chmod 775 "${work_path}"/06_bigwig/*/*.bw

echo ""
echo "=== 3-REAP pipeline complete! ==="
