#!/usr/bin/env python3
import os
import pathlib
import glob
import pandas as pd
import argparse
import zipfile, io, re

# Notes
"""
Generates TSVs of tiered fusion candidates filtered by support/readthrough logic,
enriched with STAR-Fusion + AGFusion, then merged with pVACfuse epitopes.
Originally conceptualized by Malachi Griffith, Obi Griffith, and Kelsy Cotto.

Author: Kelsy Cotto (updated with pVACfuse fallback merge)
Date: August 2024 (rev)
"""

# ---- PARSE ARGUMENTS -------------------------------------------------------
def parse_arguments():
    parser = argparse.ArgumentParser(
        description="Filter fusions identified by FusionInspector and consolidate information from multiple fusion tools."
    )
    parser.add_argument(
        "-WB",
        required=True,
        help="Path to the gcp_immuno folder of the trial (WORKING_BASE). Should include trailing slash (e.g., /.../gcp_immuno/).",
    )
    parser.add_argument(
        "-f", "--fin_results",
        required=True,
        help="Final results folder name (under -WB). Output fusion_review/ will be created here."
    )
    return parser.parse_args()

# ---- HELPERS ---------------------------------------------------------------
def _safe_glob_one(pattern: str):
    """Return a single matching path or raise a helpful error."""
    paths = glob.glob(pattern)
    if not paths:
        raise FileNotFoundError(f"No file matched pattern: {pattern}")
    if len(paths) > 1:
        # If multiple, pick the first but warn (could be refined as needed)
        print(f"[WARN] Multiple files matched pattern; using first:\n  {pattern}\n  -> {paths[0]}")
    return paths[0]

def _parse_fasta_seqs(text: str):
    """Return list of sequences from FASTA text (no Biopython)."""
    if not text:
        return []
    seqs, cur = [], []
    for line in text.splitlines():
        if line.startswith(">"):
            if cur:
                seqs.append("".join(cur))
                cur = []
        else:
            cur.append(line.strip())
    if cur:
        seqs.append("".join(cur))
    return seqs

def _first_nonnull(df: pd.DataFrame, cols_regex):
    """Pick first existing column matching any regex; return first non-null value (as str)."""
    for pat in cols_regex:
        for c in df.columns:
            if re.search(pat, c, flags=re.I):
                col = df[c]
                if col.notna().any():
                    return str(col.dropna().iloc[0])
                return None
    return None

def _collect_agfusion_view(ag_dir=None, ag_zip=None):
    """
    Build a table with:
      FusionPartners, AGFUS_LEFT_POS, AGFUS_RIGHT_POS, AGFUS_INFRAME,
      AGFUS_CDS_SEQUENCE, AGFUS_TRANSL
    by scanning either agfusion_results/ (unzipped) or agfusion_results.zip.
    Folder format expected: GENE1-pos1_GENE2-pos2
    """
    rows = []

    def handle_folder(folder_name, file_reader, names_lister, raw_reader=None):
        # Expect: HIC2-21442857_PI4KA-20734553  (we then map to "HIC2_PI4KA")
        m = re.match(r"^([A-Za-z0-9]+)-(-?\d+)_([A-Za-z0-9]+)-(-?\d+)$", folder_name)
        if not m:
            return
        lg, lp, rg, rp = m.groups()
        base = f"{lg}_{rg}"

        # Prefer fusion_transcripts.csv; fallback to any transcript/summary csv
        preferred = f"{folder_name}/{base}.fusion_transcripts.csv"
        target = preferred if preferred in names_lister() else None
        if not target:
            for n in names_lister():
                low = n.lower()
                if n.startswith(folder_name + "/") and low.endswith(".csv") and ("transcript" in low or "summary" in low):
                    target = n
                    break

        ag_inframe = None
        if target:
            try:
                df_csv = file_reader(target)
                ag_inframe = _first_nonnull(df_csv, [r"fusion[_-]?effect", r"\bin[_-]?frame\b", r"\bframe\b", r"\borf"])
            except Exception:
                pass

        cds_member = f"{folder_name}/{base}_cds.fa"
        prot_member = f"{folder_name}/{base}_protein.fa"

        cds_text = raw_reader(cds_member) if raw_reader else None
        prot_text = raw_reader(prot_member) if raw_reader else None

        if ag_dir:
            cds_path  = os.path.join(ag_dir, cds_member)
            prot_path = os.path.join(ag_dir, prot_member)
            if os.path.exists(cds_path):
                cds_text = open(cds_path, "r").read()
            if os.path.exists(prot_path):
                prot_text = open(prot_path, "r").read()

        cds_seqs  = _parse_fasta_seqs(cds_text)
        prot_seqs = _parse_fasta_seqs(prot_text)

        rows.append({
            "FusionPartners": base,
            "AGFUS_LEFT_POS":  int(lp),
            "AGFUS_RIGHT_POS": int(rp),
            "AGFUS_INFRAME": ag_inframe,
            "AGFUS_CDS_SEQUENCE":  ";".join(cds_seqs)  if cds_seqs  else None,
            "AGFUS_TRANSL":        ";".join(prot_seqs) if prot_seqs else None,
        })

    if ag_dir and os.path.isdir(ag_dir):
        def list_names():
            out = []
            for d in os.listdir(ag_dir):
                p = os.path.join(ag_dir, d)
                if os.path.isdir(p):
                    for f in os.listdir(p):
                        out.append(f"{d}/{f}")
            return out

        def file_reader(path_rel):
            return pd.read_csv(os.path.join(ag_dir, path_rel))

        for d in os.listdir(ag_dir):
            if os.path.isdir(os.path.join(ag_dir, d)):
                handle_folder(d, file_reader, list_names, raw_reader=None)

    elif ag_zip and os.path.isfile(ag_zip):
        with zipfile.ZipFile(ag_zip) as zf:
            names = zf.namelist()
            folders = sorted({n.split("/")[0] for n in names if "/" in n})

            def list_names():
                return names

            def file_reader(path_rel):
                with zf.open(path_rel) as fh:
                    return pd.read_csv(io.TextIOWrapper(fh, encoding="utf-8"))

            def raw_reader(member):
                try:
                    with zf.open(member) as fh:
                        return io.TextIOWrapper(fh, encoding="utf-8").read()
                except Exception:
                    return None

            for d in folders:
                handle_folder(d, file_reader, list_names, raw_reader)
    else:
        return None

    return pd.DataFrame(rows).drop_duplicates(subset=["FusionPartners"]) if rows else None

def _build_pvac_maps(pattern):
    """
    Returns two dicts:
      - map_nocat:   'Partners.Transcripts' -> 'frameshift_fusion'|'inframe_fusion'
      - map_partners:'Partners' -> list of (Transcripts, category)

    Parses pVACfuse aggregated 'ID' which looks like:
      <num>.<FusionPartners.Transcripts.Category>.<num>
    """
    files = glob.glob(pattern)
    if not files:
        return {}, {}
    df = pd.read_csv(files[0], sep="\t")
    fk_full = df["ID"].str.replace(
        r"^\d+\.(.*?\.(?:frameshift_fusion|inframe_fusion))\.\d+$", r"\1", regex=True
    )
    fk_nocat = fk_full.str.rsplit(".", n=1).str[0]  # Partners.Transcripts
    fk_cat   = fk_full.str.rsplit(".", n=1).str[1]  # inframe_fusion|frameshift_fusion
    partners = fk_nocat.str.split(".").str[0]       # Partners
    transcripts = fk_nocat.str.split(".").str[1]    # ENSTx_ENSTy

    map_nocat, map_partners = {}, {}
    for p, t, c, k in zip(partners, transcripts, fk_cat, fk_nocat):
        if c not in {"inframe_fusion", "frameshift_fusion"}:
            continue
        map_nocat.setdefault(k, c)
        map_partners.setdefault(p, []).append((t, c))
    return map_nocat, map_partners

def _norm_cds(x):
    """Normalize transcript IDs: convert '.', 'nan', 'None' -> None."""
    x = str(x)
    return None if x in {".", "nan", "None"} else x

# ---- MAIN ------------------------------------------------------------------
def main():
    args = parse_arguments()
    final_result = args.fin_results

    fusioninspector_path = (
        f"{args.WB}{final_result}/rnaseq/fusioninspector_evidence/finspector.FusionInspector.fusions.tsv"
    )
    fusion_candidates_df = pd.read_csv(fusioninspector_path, sep="\t")

    fusion_dir = f"{args.WB}{final_result}/fusion_review"
    pathlib.Path(fusion_dir).mkdir(parents=True, exist_ok=True)

    # --- Support filters
    filt1 = (fusion_candidates_df["JunctionReadCount"] + fusion_candidates_df["SpanningFragCount"]) > 5
    filt2 = (fusion_candidates_df["JunctionReadCount"] > 0)
    fusion_candidates_df_filt = fusion_candidates_df.loc[filt1 & filt2].copy()

    # --- Readthrough logic (KEEP if any of the following holds):
    # (a) gene partner in Cancer Gene Census (Role contains 'fusion'), OR
    # (b) partners on different chromosomes, OR
    # (c) same chromosome AND (different strands OR abs(posL - posR) > 1,000,000)
    # NOTE: The code keeps same-chr/same-strand fusions only if they are FAR apart (>1Mb).
    # (This matches your code; original comment said '< 1,000,000', which would be opposite.)
    cancer_genes = pd.read_csv("~/workspace/fusion_review/scripts/CancerGeneCensus-Mar2023.tsv", sep="\t")
    fusion_genes = cancer_genes[cancer_genes["Role in Cancer"].str.contains("fusion", case=False, na=False)]["Gene Symbol"].tolist()

    # Build boolean condition
    left_chr  = fusion_candidates_df_filt["LeftBreakpoint"].str.split(":").str[0]
    right_chr = fusion_candidates_df_filt["RightBreakpoint"].str.split(":").str[0]
    left_pos  = fusion_candidates_df_filt["LeftBreakpoint"].str.split(":").str[1].astype(int)
    right_pos = fusion_candidates_df_filt["RightBreakpoint"].str.split(":").str[1].astype(int)
    left_str  = fusion_candidates_df_filt["LeftBreakpoint"].str.split(":").str[2]
    right_str = fusion_candidates_df_filt["RightBreakpoint"].str.split(":").str[2]

    partners_a = fusion_candidates_df_filt["#FusionName"].str.split("--").str[0].isin(fusion_genes)
    partners_b = fusion_candidates_df_filt["#FusionName"].str.split("--").str[1].isin(fusion_genes)
    diff_chr   = (left_chr != right_chr)
    same_chr   = ~diff_chr
    diff_str   = (left_str != right_str)
    far_apart  = (left_pos - right_pos).abs() > 1_000_000

    condition = partners_a | partners_b | diff_chr | (same_chr & (diff_str | far_apart))
    filtered_df = fusion_candidates_df_filt.loc[condition].copy()

    # --- Merge STAR-Fusion coding effects
    star_fusion_path = f"{args.WB}{final_result}/rnaseq/star_fusion/results/star-fusion.fusion_predictions.abridged.coding_effect.tsv"
    star_fusion_df = pd.read_csv(star_fusion_path, sep="\t")

    merged_filtered_df = pd.merge(
        filtered_df,
        star_fusion_df,
        on=["#FusionName", "LeftBreakpoint", "RightBreakpoint"],
        suffixes=('_filtered_fusioninspector', '_starfusion'),
        how="left"
    )

    # Normalize partner label & numeric positions for FI
    merged_filtered_df["FusionPartners"] = merged_filtered_df["#FusionName"].str.replace("--", "_", regex=False)
    merged_filtered_df["FI_LEFT_POS"]  = merged_filtered_df["LeftBreakpoint"].str.split(":").str[1].astype(int)
    merged_filtered_df["FI_RIGHT_POS"] = merged_filtered_df["RightBreakpoint"].str.split(":").str[1].astype(int)

    # --- AGFusion enrichment (optional)
    BP_TOLERANCE = 25
    ag_dir = f"{args.WB}{final_result}/rnaseq/star_fusion/results/agfusion_results"
    ag_zip = f"{args.WB}{final_result}/rnaseq/star_fusion/results/agfusion_results.zip"
    agfusion_view = _collect_agfusion_view(ag_dir=ag_dir if os.path.isdir(ag_dir) else None,
                                           ag_zip=ag_zip if os.path.isfile(ag_zip) else None)

    if agfusion_view is not None and len(agfusion_view) > 0:
        merged_filtered_df = merged_filtered_df.merge(agfusion_view, on="FusionPartners", how="left")
        merged_filtered_df["AGFUS_BP_MATCH"] = (
            merged_filtered_df["AGFUS_LEFT_POS"].notna()
            & merged_filtered_df["AGFUS_RIGHT_POS"].notna()
            & ((merged_filtered_df["FI_LEFT_POS"]  - merged_filtered_df["AGFUS_LEFT_POS"]).abs()  <= BP_TOLERANCE)
            & ((merged_filtered_df["FI_RIGHT_POS"] - merged_filtered_df["AGFUS_RIGHT_POS"]).abs() <= BP_TOLERANCE)
        )

    for col in ["AGFUS_LEFT_CDS","AGFUS_RIGHT_CDS","AGFUS_LEFT_TX","AGFUS_RIGHT_TX"]:
        if col in merged_filtered_df.columns:
            merged_filtered_df.drop(columns=[col], inplace=True)

    # --- Build FusionCategory (STAR/AGFusion) and TranscriptKey
    star_cat = merged_filtered_df["PROT_FUSION_TYPE"].map({
        "FRAMESHIFT": "frameshift_fusion",
        "INFRAME":    "inframe_fusion"
    })

    def _ag_to_cat(v):
        if v is None:
            return None
        s = str(v).lower()
        if "in-frame" in s or s.strip() == "inframe":
            return "inframe_fusion"
        if "out-of-frame" in s or "frameshift" in s:
            return "frameshift_fusion"
        return None

    ag_cat = merged_filtered_df.get("AGFUS_INFRAME", pd.Series([None]*len(merged_filtered_df))).map(_ag_to_cat)
    merged_filtered_df["FusionCategory"] = star_cat.fillna(ag_cat).fillna("no_cds_prediction")

    # Normalize transcript IDs and construct TranscriptKey
    merged_filtered_df["CDS_LEFT_ID_N"]  = merged_filtered_df["CDS_LEFT_ID"].map(_norm_cds)
    merged_filtered_df["CDS_RIGHT_ID_N"] = merged_filtered_df["CDS_RIGHT_ID"].map(_norm_cds)
    merged_filtered_df["TranscriptKey"] = (
        merged_filtered_df["CDS_LEFT_ID_N"].fillna("NA").astype(str) + "_" +
        merged_filtered_df["CDS_RIGHT_ID_N"].fillna("NA").astype(str)
    )

    # --- pVACfuse-driven upgrade for no_cds_prediction and/or NA transcripts
    map_i_nocat,  map_i_part  = _build_pvac_maps(f"{args.WB}{final_result}/pVACfuse/mhc_i/*.all_epitopes.aggregated.tsv")
    map_ii_nocat, map_ii_part = _build_pvac_maps(f"{args.WB}{final_result}/pVACfuse/mhc_ii/*.all_epitopes.aggregated.tsv")

    merged_filtered_df["FusionKey_nocat"] = merged_filtered_df["FusionPartners"] + "." + merged_filtered_df["TranscriptKey"]

    def _upgrade_cat_and_transcripts(row):
        # If we already have a good category and non-NA transcripts, keep them
        if row["FusionCategory"] != "no_cds_prediction" and "NA" not in row["TranscriptKey"]:
            return pd.Series([row["FusionCategory"], row["TranscriptKey"]])

        partners = row["FusionPartners"]
        nocat    = row["FusionKey_nocat"]

        # 1) Exact Partners.Transcripts match: use its category
        cat = map_i_nocat.get(nocat) or map_ii_nocat.get(nocat)
        if cat in {"inframe_fusion","frameshift_fusion"}:
            return pd.Series([cat, row["TranscriptKey"]])

        # 2) Fallback: partners-only — adopt first pVACfuse transcript pair & category
        cand_list = (map_i_part.get(partners) or []) + (map_ii_part.get(partners) or [])
        if cand_list:
            t, c = cand_list[0]  # pick the first; refine if needed
            return pd.Series([c, t])

        # 3) No upgrade available
        return pd.Series([row["FusionCategory"], row["TranscriptKey"]])

    merged_filtered_df[["FusionCategory","TranscriptKey"]] = (
        merged_filtered_df.apply(_upgrade_cat_and_transcripts, axis=1)
    )

    # Final FusionKey
    merged_filtered_df["FusionKey"] = (
        merged_filtered_df["FusionPartners"] + "." +
        merged_filtered_df["TranscriptKey"] + "." +
        merged_filtered_df["FusionCategory"]
    )

    # --- Merge with pVACfuse aggregated outputs (MHC I & II)
    try:
        pvacfuse_file_mhc_i = _safe_glob_one(f"{args.WB}{final_result}/pVACfuse/mhc_i/*.all_epitopes.aggregated.tsv")
        pvacfuse_df_mhc_i = pd.read_csv(pvacfuse_file_mhc_i, sep="\t")
        pvacfuse_df_mhc_i["FusionKey"] = pvacfuse_df_mhc_i["ID"].str.replace(
            r"^\d+\.(.*?\.(?:frameshift_fusion|inframe_fusion))\.\d+$", r"\1", regex=True
        )
        merged_class_i_df = pd.merge(
            pvacfuse_df_mhc_i,
            merged_filtered_df,
            on="FusionKey",
            how="inner"
        )
        for c in ["JunctionReads","SpanningFrags","CounterFusionLeftReads","CounterFusionRightReads"]:
            if c in merged_class_i_df.columns:
                merged_class_i_df.drop(columns=[c], inplace=True)
        merged_class_i_df.to_csv(
            f"{fusion_dir}/tumor-exome.all_epitopes.aggregated.mhc_i.review.tsv",
            sep="\t",
            index=False
        )
        print(f"[OK] Wrote MHC I review: {fusion_dir}/tumor-exome.all_epitopes.aggregated.mhc_i.review.tsv  "
              f"(rows: {len(merged_class_i_df)})")
    except FileNotFoundError as e:
        print(f"[WARN] Skipping MHC I merge: {e}")

    try:
        pvacfuse_file_mhc_ii = _safe_glob_one(f"{args.WB}{final_result}/pVACfuse/mhc_ii/*.all_epitopes.aggregated.tsv")
        pvacfuse_df_mhc_ii = pd.read_csv(pvacfuse_file_mhc_ii, sep="\t")
        pvacfuse_df_mhc_ii["FusionKey"] = pvacfuse_df_mhc_ii["ID"].str.replace(
            r"^\d+\.(.*?\.(?:frameshift_fusion|inframe_fusion))\.\d+$", r"\1", regex=True
        )
        merged_class_ii_df = pd.merge(
            pvacfuse_df_mhc_ii,
            merged_filtered_df,
            on="FusionKey",
            how="inner"
        )
        for c in ["JunctionReads","SpanningFrags","CounterFusionLeftReads","CounterFusionRightReads"]:
            if c in merged_class_ii_df.columns:
                merged_class_ii_df.drop(columns=[c], inplace=True)
        merged_class_ii_df.to_csv(
            f"{fusion_dir}/tumor-exome.all_epitopes.aggregated.mhc_ii.review.tsv",
            sep="\t",
            index=False
        )
        print(f"[OK] Wrote MHC II review: {fusion_dir}/tumor-exome.all_epitopes.aggregated.mhc_ii.review.tsv "
              f"(rows: {len(merged_class_ii_df)})")
    except FileNotFoundError as e:
        print(f"[WARN] Skipping MHC II merge: {e}")

    # Small summary
    print(f"[Summary] FusionInspector input: {len(fusion_candidates_df)}")
    print(f"[Summary] After support filters: {len(fusion_candidates_df_filt)}")
    print(f"[Summary] After readthrough filters: {len(filtered_df)}")
    print(f"[Summary] Unique fusions after STAR/AGFusion merge: {merged_filtered_df['#FusionName'].nunique()}")

if __name__ == "__main__":
    main()
