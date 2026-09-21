#!/usr/bin/env bash
# Convert PolyA_DB v4.1 TSV to bigBed for UCSC track hub.
# Run this on the compute server where the TSV lives.
# Output: output/trackhub/hg38.PAS.v4.1.bb  (upload to EC2 hg38/ directory)
#
# Usage: bash make_polyaDB_bigbed.sh
#
# Requires: python3, sort, bedToBigBed (UCSC tools)

set -euo pipefail

PROJECT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
PAS_TSV="${PROJECT_DIR}/data/PolyA_DB_v4.1/hg38.PAS.main.tsv"
CHROM_SIZES="${PROJECT_DIR}/data/ref/hg38.chrom.sizes"
TMP_BED=$(mktemp /tmp/polyaDB_XXXXXX.bed)
OUT_BB="${PROJECT_DIR}/output/trackhub/hg38.PAS.v4.1.bb"

echo "Input TSV   : ${PAS_TSV}"
echo "Chrom sizes : ${CHROM_SIZES}"
echo "Output      : ${OUT_BB}"
echo ""

# ── Step 1: print the header columns so we can verify the format ────────────
echo "=== TSV columns ==="
head -1 "${PAS_TSV}" | tr '\t' '\n' | nl -ba
echo ""

# ── Step 2: TSV → sorted BED6 ────────────────────────────────────────────────
# PAS_ID format: chr:strand:position  (position is 0-based per PolyA_DB v4.1)
# name field   : PAS_ID|GeneSymbol|PAS_type
# score field  : PolyaStrength_percentile × 10  (0–1000; shades sites by strength in dense view)
#
# NOTE: pipe goes on the opening <<'PYEOF' line; closing PYEOF must be alone on its line.
python3 - "${PAS_TSV}" <<'PYEOF' | sort -k1,1 -k2,2n > "${TMP_BED}"
import sys

infile = sys.argv[1]

GENE_CANDIDATES = [
    "GeneSymbol", "Gene_Symbol", "Gene_symbol", "GeneName",
    "gene_symbol", "Gene", "gene", "Gene_Name", "GENE",
]
TYPE_CANDIDATES = [
    "PAS_type_RefSeq_label", "PAS_type", "Type", "type",
]
STRENGTH_CANDIDATES = [
    "PolyaStrength_percentile", "PolyaStrength", "strength",
]

with open(infile) as fh:
    header = fh.readline().rstrip("\n").split("\t")

    if "PAS_ID" not in header:
        sys.exit(f"ERROR: PAS_ID column not found. Columns: {header}")
    pas_col = header.index("PAS_ID")

    gene_col = next((header.index(c) for c in GENE_CANDIDATES if c in header), None)
    type_col = next((header.index(c) for c in TYPE_CANDIDATES if c in header), None)
    strength_col = next((header.index(c) for c in STRENGTH_CANDIDATES if c in header), None)

    print(f"# gene col     : {header[gene_col] if gene_col is not None else 'not found'}", file=sys.stderr)
    print(f"# type col     : {header[type_col] if type_col is not None else 'not found'}", file=sys.stderr)
    print(f"# strength col : {header[strength_col] if strength_col is not None else 'not found'}", file=sys.stderr)

    def clean(val):
        """Strip leading/trailing whitespace and remove embedded tabs/newlines
        that would create extra BED fields and break bedToBigBed parsing."""
        return val.replace('\t', ' ').replace('\r', '').replace('\n', '').strip()

    skipped = 0
    for line in fh:
        fields = line.rstrip("\n").split("\t")
        pas_id = fields[pas_col]

        try:
            parts  = pas_id.split(":")
            chrom  = parts[0]
            strand = parts[1]
            pos    = int(parts[2])   # 0-based
        except (IndexError, ValueError):
            skipped += 1
            continue

        # name: PAS_ID|Gene|Type  (tabs sanitized so they don't split into extra fields)
        parts_name = [pas_id]
        if gene_col is not None and gene_col < len(fields):
            g = clean(fields[gene_col])
            if g and g not in (".", "NA", ""):
                parts_name.append(g)
        if type_col is not None and type_col < len(fields):
            t = clean(fields[type_col])
            if t and t not in (".", "NA", ""):
                parts_name.append(t)
        name = "|".join(parts_name).replace(' ', '_')

        # score: PolyaStrength_percentile * 10 → 0–1000
        score = 0
        if strength_col is not None and strength_col < len(fields):
            try:
                score = min(1000, int(float(fields[strength_col]) * 10))
            except (ValueError, TypeError):
                score = 0

        print(f"{chrom}\t{pos}\t{pos + 1}\t{name}\t{score}\t{strand}")

    if skipped:
        print(f"# Skipped {skipped} malformed PAS_ID lines", file=sys.stderr)
PYEOF

echo "BED entries : $(wc -l < "${TMP_BED}")"
echo "Chroms found: $(cut -f1 "${TMP_BED}" | sort -u | tr '\n' ' ')"
echo ""

# ── Step 3: BED → bigBed ─────────────────────────────────────────────────────
echo "=== First line of BED (field check) ==="
head -1 "${TMP_BED}" | cat -A
echo ""

if ! bedToBigBed "${TMP_BED}" "${CHROM_SIZES}" "${OUT_BB}"; then
    echo ""
    echo "ERROR: bedToBigBed failed — BED file preserved at ${TMP_BED}"
    echo "Inspect with: cut -f5 ${TMP_BED} | sort -u | head -20"
    exit 1
fi
rm -f "${TMP_BED}"

echo ""
echo "Done: ${OUT_BB}"
echo "Upload this file to the EC2 under the hub's hg38/ directory."
