#!/usr/bin/env python3
"""Generate the UCSC track hub for cmoore_021.CFIm25_OE (3-REAP output).

Layout mirrors the cmoore_020 hub (test_data/trackHub_example/trackhub):
  <outdir>/hub.txt
  <outdir>/genomes.txt
  <outdir>/hg38/trackDb.txt
bigWigs are expected on the server at <base-url>/<set>/<strand>.<sample>.<set>.bw,
i.e. the 06_bigwig/* subfolders copied into the hub's hg38/ directory.

Usage:
  python3 scripts/make_trackhub.py [--base-url URL] [--outdir DIR] [--tree TREE]
"""
import argparse
import os
import re
import sys

PROJECT_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

HUB_NAME = "cmoore_021"
EMAIL = "cnobrega@mdibl.org"
DEFAULT_BASE_URL = "http://3.83.106.170/cmoore_021/hg38"
PAS_BB = "hg38.PAS.v4.1.bb"

# condition key -> (genotype, treatment, short tag, long label, reps, rep colours)
CONDITIONS = [
    ("control_control", "control", "control", "Ctl",       "Control / untreated",
     [1, 2, 3], ["60,60,60", "120,120,120", "180,180,180"]),
    ("control_TGFB",    "control", "TGFB",    "CtlTGFB",   "Control / TGFB",
     [1, 2, 3], ["200,90,0", "235,130,30", "250,175,90"]),
    ("CFIm25_control",  "CFIm25",  "control", "CF25",      "CFIm25 OE / untreated",
     [4, 5, 6], ["30,50,170", "60,100,220", "110,150,250"]),
    ("CFIm25_TGFB",     "CFIm25",  "TGFB",    "CF25TGFB",  "CFIm25 OE / TGFB",
     [4, 5, 6], ["170,0,60", "215,40,100", "245,100,150"]),
]

STRANDS = [("plus", "+", "plus strand"), ("mi", "mi", "mirror strand")]

SETS = {
    "PAS": "PASS_bw_LAP24_PAS",
    "PASSreads": "PASS_bw_LAP24",
    "mapped": "mapped_read",
    "LAPs": "PASS_bw_LAP24_positon",
}


def samples():
    for key, geno, treat, tag, label, reps, colours in CONDITIONS:
        for rep, col in zip(reps, colours):
            yield dict(cond=key, geno=geno, treat=treat, tag=tag, label=label,
                       rep=rep, color=col, name=f"{key}_rep_{rep}")


def bw_url(base, set_key, strand, sample):
    s = SETS[set_key]
    return f"{base}/{s}/{strand}.{sample}.{s}.bw"


def block(**kv):
    lines = []
    for k, v in kv.items():
        lines.append(f"{k} {v}")
    return "\n".join(lines) + "\n\n"


def header(title):
    bar = "#" * 46
    return f"{bar}\n# {title}\n{bar}\n"


def build_trackdb(base):
    out = []

    # 1. PolyA_DB reference -------------------------------------------------
    out.append(header("Reference: PolyA_DB v4.1 PAS annotations"))
    out.append(block(track="polyaDB_v41", bigDataUrl=f"{base}/{PAS_BB}",
                     shortLabel="PolyA_DB v4.1",
                     longLabel="PolyA_DB v4.1 — hg38 polyadenylation sites (PAS_ID|Gene)",
                     type="bigBed 6", visibility="pack", priority="0.5",
                     color="160,32,240", itemRgb="off"))

    # 2. PAS composites: condition x strand ---------------------------------
    out.append(header("PAS composites  [8 tracks: condition x strand]"))
    prio = 1.0
    for strand, sym, sdesc in STRANDS:
        for key, _g, _t, tag, label, _r, _c in CONDITIONS:
            parent = f"PAS_{key}_{strand}"
            out.append(block(track=parent, compositeTrack="on", priority=f"{prio:.2f}",
                             shortLabel=f"{tag} PAS {sym}",
                             longLabel=f"Detected PASs — {label}, {sdesc}",
                             type="bigWig",
                             visibility="full" if strand == "plus" else "hide",
                             autoScale="group", alwaysZero="on", smoothingWindow="off",
                             windowingFunction="maximum", maxHeightPixels="120:80:40"))
            prio += 0.05
            for s in samples():
                if s["cond"] != key:
                    continue
                out.append(block(track=f"PAS_{s['name']}_{strand}", parent=parent,
                                 bigDataUrl=bw_url(base, "PAS", strand, s["name"]),
                                 shortLabel=f"PAS{sym} {tag} r{s['rep']}",
                                 longLabel=f"Detected PASs — {label} rep {s['rep']}, {sdesc}",
                                 type="bigWig", visibility="full",
                                 graphTypeDefault="bar", color=s["color"]))

    # 3. PASS-read overlays (multiWig) --------------------------------------
    out.append(header("PASSreads overlays  [8 tracks: condition x strand, points]"))
    prio = 2.0
    for strand, sym, sdesc in STRANDS:
        for key, _g, _t, tag, label, _r, _c in CONDITIONS:
            parent = f"PASSreads_{key}_{strand}"
            out.append(block(track=parent, type="bigWig", container="multiWig",
                             aggregate="transparentOverlay", showSubtrackColorOnUi="on",
                             shortLabel=f"{tag} PASS {sym}",
                             longLabel=f"PASS reads — {label}, {sdesc} (reps overlay)",
                             visibility="full" if strand == "plus" else "hide",
                             autoScale="on", alwaysZero="on", graphTypeDefault="points",
                             maxHeightPixels="120:80:40", priority=f"{prio:.2f}"))
            prio += 0.05
            for s in samples():
                if s["cond"] != key:
                    continue
                out.append(block(track=f"PASSreads_ov_{s['name']}_{strand}", parent=parent,
                                 bigDataUrl=bw_url(base, "PASSreads", strand, s["name"]),
                                 shortLabel=f"rep {s['rep']}",
                                 longLabel=f"{label} rep {s['rep']}, {sdesc}",
                                 type="bigWig", smoothingWindow="2", color=s["color"]))

    # 4-6. Per-replicate matrix composites ----------------------------------
    matrix = [
        ("PASSreads", "PASSreads", "PASS (per-rep)", "PASS reads — per-replicate detail",
         "hide", "points", 2.5, "Composite: PASS reads  [hide — per-rep detail view]"),
        ("mapped", "mapped", "mapped", "All aligned reads",
         "dense", "bar", 3.0, "Composite: All aligned reads  [dense — QC overview]"),
        ("LAPs", "LAPs", "LAPs", "PASS LAPs (last aligned positions)",
         "hide", "bar", 4.0, "Composite: PASS LAPs  [hide — enable when needed]"),
    ]
    rep_ids = sorted({s["rep"] for s in samples()})
    rep_sg = " ".join(f"rep{r}=rep{r}" for r in rep_ids)
    for comp, set_key, slab, llab, vis, graph, prio, title in matrix:
        out.append(header(title))
        out.append(block(track=comp, compositeTrack="on", priority=f"{prio}",
                         shortLabel=slab, longLabel=llab, type="bigWig", visibility=vis,
                         subGroup1="genotype genotype control=Control CFIm25=CFIm25_OE",
                         subGroup2="treatment treatment control=Untreated TGFB=TGFB",
                         subGroup3=f"rep rep {rep_sg}",
                         subGroup4="strand strand mi=Mirror plus=Plus",
                         sortOrder="strand=+ genotype=+ treatment=+ rep=+",
                         dimensions="dimX=genotype dimY=treatment dimA=rep",
                         dragAndDrop="subtracks", filterComposite="dimA",
                         autoScale="group"))
        for strand, sym, sdesc in STRANDS:
            for s in samples():
                kv = dict(track=f"{comp}_{s['name']}_{strand}", parent=comp,
                          bigDataUrl=bw_url(base, set_key, strand, s["name"]),
                          shortLabel=f"{s['tag']} r{s['rep']} {sym}",
                          longLabel=f"{s['label']} rep {s['rep']}, {sdesc}",
                          type="bigWig",
                          subGroups=(f"genotype={s['geno']} treatment={s['treat']} "
                                     f"rep=rep{s['rep']} strand={strand}"),
                          visibility="full", alwaysZero="on", autoScale="on",
                          graphTypeDefault=graph, maxHeightPixels="120:80:40",
                          color=s["color"])
                out.append(block(**kv))

    return "".join(out).rstrip() + "\n"


def validate(trackdb, tree_path):
    """Sanity checks: unique names, parents exist, label lengths, files in tree.txt."""
    errors = []
    stanzas = [s for s in trackdb.split("\n\n") if "track " in s]
    names, parents, files = [], [], []
    for st in stanzas:
        d = {}
        for line in st.splitlines():
            if line.startswith("#") or not line.strip():
                continue
            k, _, v = line.partition(" ")
            d[k] = v
        names.append(d["track"])
        if "parent" in d:
            parents.append(d["parent"])
        if len(d.get("shortLabel", "")) > 17:
            errors.append(f"shortLabel >17 chars: {d['track']}: {d['shortLabel']!r}")
        if len(d.get("longLabel", "")) > 76:
            errors.append(f"longLabel >76 chars: {d['track']}")
        if "bigDataUrl" in d:
            files.append(d["bigDataUrl"].rsplit("/", 1)[1])
    dups = {n for n in names if names.count(n) > 1}
    if dups:
        errors.append(f"duplicate track names: {sorted(dups)}")
    missing_parents = set(parents) - set(names)
    if missing_parents:
        errors.append(f"undefined parents: {sorted(missing_parents)}")

    if tree_path and os.path.exists(tree_path):
        with open(tree_path, encoding="utf-8") as fh:
            tree_files = set(re.findall(r"[─ ]+(\S+\.(?:bw|bb))\s*$", fh.read(), re.M))
        bw_missing = sorted({f for f in files if f.endswith(".bw")} - tree_files)
        if bw_missing:
            errors.append(f"{len(bw_missing)} bigWigs not in tree.txt, e.g. {bw_missing[:5]}")
    print(f"tracks: {len(names)}  bigDataUrls: {len(files)} "
          f"(unique files: {len(set(files))})")
    return errors


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--base-url", default=DEFAULT_BASE_URL)
    ap.add_argument("--outdir", default=os.path.join(PROJECT_DIR, "output", "3REAP", "trackhub"))
    ap.add_argument("--tree", default=os.path.join(PROJECT_DIR, "output", "3REAP", "tree.txt"))
    args = ap.parse_args()
    base = args.base_url.rstrip("/")

    os.makedirs(os.path.join(args.outdir, "hg38"), exist_ok=True)
    with open(os.path.join(args.outdir, "hub.txt"), "w") as fh:
        fh.write(f"hub {HUB_NAME}\n"
                 "shortLabel Moore CFIm25 OE\n"
                 "longLabel 3-REAP APA results: CFIm25 OE x TGFB (Claire Moore Lab, Tufts)\n"
                 "genomesFile genomes.txt\n"
                 f"email {EMAIL}\n")
    with open(os.path.join(args.outdir, "genomes.txt"), "w") as fh:
        fh.write("genome hg38\ntrackDb hg38/trackDb.txt\n")

    trackdb = build_trackdb(base)
    with open(os.path.join(args.outdir, "hg38", "trackDb.txt"), "w") as fh:
        fh.write(trackdb)

    errors = validate(trackdb, args.tree)
    print(f"wrote hub to {args.outdir}")
    if errors:
        print("VALIDATION ERRORS:", *errors, sep="\n  ", file=sys.stderr)
        sys.exit(1)
    print("validation OK")


if __name__ == "__main__":
    main()
