import os
import subprocess
import pathlib
import glob
import pandas as pd
import argparse

# Notes
"""
Generates a tsv file containing tiered fusion candidates that have been filtered based on criteria set by Malachi Griffith, Obi Griffith, and Kelsy Cotto.

Author: Kelsy Cotto
Date: August 2024
"""


# ---- PARSE ARGUMENTS -------------------------------------------------------
# Parses command line arguments
# Enables user help
# Future impovements: require the user to enter either the WB OR the list of files
def parse_arguments():
    # Parse command line arugments
    parser = argparse.ArgumentParser(
        description="Filter fusions identified by FusionInspector and consolidate information from multiple fusion tools."
    )

    parser.add_argument(
        "-WB",
        help="the path to the gcp_immuno folder of the trial you wish to run script on, defined as WORKING_BASE in envs.txt",
    )

    # The name of the final results folder
    parser.add_argument(
        "-f", "--fin_results", help="Name of the folder you want the results to be written in. Note, will create a new folder called fusion_review in this directory to write output files to"
    )

    return parser.parse_args()


# home_dir = os.path.expanduser("~")
# sample = "JLF-100-021"

args = parse_arguments()
print(args)
final_result = args.fin_results
fusion_candidates_df = pd.read_csv(
    args.WB
    + final_result
    + "/rnaseq/fusioninspector_evidence/finspector.FusionInspector.fusions.tsv",
    delimiter="\t",
)

fusion_dir = f"{args.WB}{final_result}/fusion_review"
pathlib.Path(fusion_dir).mkdir(parents=True, exist_ok=True)

##filter based on support

filt_condition = (
    fusion_candidates_df["JunctionReadCount"]
    + fusion_candidates_df["SpanningFragCount"]
) > 5
fusion_candidates_df_filt = fusion_candidates_df.loc[filt_condition]


filt_condition = fusion_candidates_df["JunctionReadCount"] > 0
fusion_candidates_df_filt = fusion_candidates_df_filt.loc[filt_condition]
##filter readthroughs
# NOT a readthrough. Defined as:
# [Left Chr and Right Chr are different] OR
# [chromosome are the same BUT Left Strand and Right Strand are different] OR
# [chromosome and strand are the same BUT ABS(Left Pos - Right Pos) < 1,000,000] OR
# [Fusion GeneA Name OR Fusion GeneB Name matches a known fusion driver gene]

cancer_genes = pd.read_csv(f"/opt/scripts/CancerGeneCensus-Mar2023.tsv", delimiter="\t")
fusion_genes_df = cancer_genes[
    cancer_genes["Role in Cancer"].str.contains("fusion", case=False, na=False)
]
fusion_genes = fusion_genes_df["Gene Symbol"].tolist()

# Create the OR condition
condition = (
    #is the gene a known fusion gene partner according to cancer gene census
    (fusion_candidates_df_filt["#FusionName"].str.split("--").str[0].isin(fusion_genes))
    | (
        fusion_candidates_df_filt["#FusionName"]
        .str.split("--")
        .str[1]
        .isin(fusion_genes)
    )
    | (
        #are the fusion partners on different chromosomes
        fusion_candidates_df_filt["LeftBreakpoint"].str.split(":").str[0]
        != fusion_candidates_df_filt["RightBreakpoint"].str.split(":").str[0]
    )
    | (
        (
            # if the fusion partners are on the same chromosome, are they on different strands
            fusion_candidates_df_filt["LeftBreakpoint"].str.split(":").str[0]
            == fusion_candidates_df_filt["RightBreakpoint"].str.split(":").str[0]
        )
        & (
            (
                fusion_candidates_df_filt["LeftBreakpoint"].str.split(":").str[2]
                != fusion_candidates_df_filt["RightBreakpoint"].str.split(":").str[2]
            )
            | (
                # if the fusion partners are on the same chromosome and strand, are they 1,000,000 bp apart? i.e. not a readthrough event (maybe we want to change this in the future?)
                abs(
                    fusion_candidates_df_filt["LeftBreakpoint"]
                    .str.split(":")
                    .str[1]
                    .astype(int)
                    - fusion_candidates_df_filt["RightBreakpoint"]
                    .str.split(":")
                    .str[1]
                    .astype(int)
                )
                > 1000000
            )
        )
    )
)

# Filter the DataFrame for Class I
filtered_df = fusion_candidates_df_filt.loc[condition]
star_fusion_df = pd.read_csv(
    args.WB
    + final_result
    + "/rnaseq/star_fusion/results/star-fusion.fusion_predictions.abridged.coding_effect.tsv",
    delimiter="\t",
)

merged_filtered_df = pd.merge(
    filtered_df,
    star_fusion_df,
    on=["#FusionName", "LeftBreakpoint", "RightBreakpoint"],
    suffixes=('_filtered_fusioninspector', '_starfusion'), # Optional: distinguish overlapping columns
    how="left" #only keep things that passed filtering
    )

# Step 1: Normalize your fusion name
merged_filtered_df["FusionPartners"] = merged_filtered_df["#FusionName"].str.replace("--", "_", regex=False)

# # Step 2: Create a concatenated transcript field for matching
# merged_filtered_df["TranscriptKey"] = merged_filtered_df["CDS_LEFT_ID"].astype(str) + "_" + merged_filtered_df["CDS_RIGHT_ID"].astype(str)

# merged_filtered_df["FusionCategory"] = merged_filtered_df["PROT_FUSION_TYPE"].map({
#     "FRAMESHIFT": "frameshift_fusion",
#     "INFRAME": "inframe_fusion"
# })

# merged_filtered_df["FusionKey"] = merged_filtered_df["FusionPartners"] + '.' + merged_filtered_df["TranscriptKey"] + '.' + merged_filtered_df["FusionCategory"]

# ---------- AGFusion sequences + breakpoint cross-check ----------
import os, zipfile, io, re

BP_TOLERANCE = 25  # bp window to accept FI↔AGFusion breakpoint agreement

# Parse numeric positions from FusionInspector breakpoints for cross-check
merged_filtered_df["FI_LEFT_POS"]  = merged_filtered_df["LeftBreakpoint"].str.split(":").str[1].astype(int)
merged_filtered_df["FI_RIGHT_POS"] = merged_filtered_df["RightBreakpoint"].str.split(":").str[1].astype(int)

def _read_text_from_zip(zip_path, member):
    try:
        with zipfile.ZipFile(zip_path) as zf:
            if member in zf.namelist():
                with zf.open(member) as fh:
                    return io.TextIOWrapper(fh, encoding="utf-8").read()
    except Exception:
        pass
    return None

def _parse_fasta_seqs(text):
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

def _first_nonnull(df, cols_regex):
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
    Build a small table with:
      FusionPartners, AGFUS_LEFT_POS, AGFUS_RIGHT_POS,
      AGFUS_INFRAME, AGFUS_CDS_SEQUENCE, AGFUS_TRANSL
    by scanning either the unzipped agfusion_results/ or agfusion_results.zip.
    """
    rows = []

    def handle_folder(folder_name, file_reader, names_lister, raw_reader=None):
        # Expect folder like: HIC2-21442857_PI4KA-20734553
        m = re.match(r"^([A-Za-z0-9]+)-(-?\d+)_([A-Za-z0-9]+)-(-?\d+)$", folder_name)
        if not m:
            return
        lg, lp, rg, rp = m.groups()
        base = f"{lg}_{rg}"

        # Locate fusion_transcripts.csv
        preferred = f"{folder_name}/{base}.fusion_transcripts.csv"
        target = preferred if preferred in names_lister() else None
        if not target:
            # fallback: any CSV mentioning transcript/summary
            for n in names_lister():
                low = n.lower()
                if n.startswith(folder_name + "/") and low.endswith(".csv") and ("transcript" in low or "summary" in low):
                    target = n
                    break

        # Pull AGFUS_INFRAME from CSV (Fusion_effect or similar)
        ag_inframe = None
        if target:
            try:
                import pandas as _pd
                df_csv = file_reader(target)
                ag_inframe = _first_nonnull(df_csv, [r"fusion[_-]?effect", r"\bin[_-]?frame\b", r"\bframe\b", r"\borf"])
            except Exception:
                pass

        # Get sequences from *_cds.fa and *_protein.fa
        cds_member = f"{folder_name}/{base}_cds.fa"
        prot_member = f"{folder_name}/{base}_protein.fa"

        cds_text = raw_reader(cds_member) if raw_reader else None
        prot_text = raw_reader(prot_member) if raw_reader else None

        if ag_dir:
            # unzipped dir path
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
        # unzipped
        def list_names():
            out = []
            for d in os.listdir(ag_dir):
                p = os.path.join(ag_dir, d)
                if os.path.isdir(p):
                    for f in os.listdir(p):
                        out.append(f"{d}/{f}")
            return out
        def file_reader(path_rel):
            import pandas as _pd
            return _pd.read_csv(os.path.join(ag_dir, path_rel))
        for d in os.listdir(ag_dir):
            if os.path.isdir(os.path.join(ag_dir, d)):
                handle_folder(d, file_reader, list_names, raw_reader=None)
    elif ag_zip and os.path.isfile(ag_zip):
        # zip
        with zipfile.ZipFile(ag_zip) as zf:
            names = zf.namelist()
            folders = sorted({n.split("/")[0] for n in names if "/" in n})
            def list_names():
                return names
            def file_reader(path_rel):
                import pandas as _pd
                with zf.open(path_rel) as fh:
                    return _pd.read_csv(io.TextIOWrapper(fh, encoding="utf-8"))
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

ag_dir = args.WB + final_result + "/rnaseq/star_fusion/results/agfusion_results"
ag_zip = args.WB + final_result + "/rnaseq/star_fusion/results/agfusion_results.zip"
agfusion_view = _collect_agfusion_view(ag_dir=ag_dir if os.path.isdir(ag_dir) else None,
                                       ag_zip=ag_zip if os.path.isfile(ag_zip) else None)

if agfusion_view is not None:
    merged_filtered_df = merged_filtered_df.merge(agfusion_view, on="FusionPartners", how="left")
    # Breakpoint agreement check
    merged_filtered_df["AGFUS_BP_MATCH"] = (
        merged_filtered_df["AGFUS_LEFT_POS"].notna()
        & merged_filtered_df["AGFUS_RIGHT_POS"].notna()
        & ((merged_filtered_df["FI_LEFT_POS"]  - merged_filtered_df["AGFUS_LEFT_POS"]).abs()  <= BP_TOLERANCE)
        & ((merged_filtered_df["FI_RIGHT_POS"] - merged_filtered_df["AGFUS_RIGHT_POS"]).abs() <= BP_TOLERANCE)
    )

# OPTIONAL: drop legacy AGFUS_* ID columns if they exist
for col in ["AGFUS_LEFT_CDS","AGFUS_RIGHT_CDS","AGFUS_LEFT_TX","AGFUS_RIGHT_TX"]:
    if col in merged_filtered_df.columns:
        merged_filtered_df.drop(columns=[col], inplace=True)

# Now (re)build keys AFTER enrichment
star_cat = merged_filtered_df["PROT_FUSION_TYPE"].map({"FRAMESHIFT": "frameshift_fusion", "INFRAME": "inframe_fusion"})
def _ag_to_cat(v):
    if v is None: return None
    s = str(v).lower()
    if "in-frame" in s or s.strip() == "inframe": return "inframe_fusion"
    if "out-of-frame" in s or "frameshift" in s:  return "frameshift_fusion"
    return None
ag_cat = merged_filtered_df.get("AGFUS_INFRAME", pd.Series([None]*len(merged_filtered_df))).map(_ag_to_cat)
merged_filtered_df["FusionCategory"] = star_cat.fillna(ag_cat).fillna("no_cds_prediction")

merged_filtered_df["TranscriptKey"] = (
    merged_filtered_df["CDS_LEFT_ID"].fillna("NA").astype(str) + "_" + merged_filtered_df["CDS_RIGHT_ID"].fillna("NA").astype(str)
)
merged_filtered_df["FusionKey"] = (
    merged_filtered_df["FusionPartners"] + "." + merged_filtered_df["TranscriptKey"] + "." + merged_filtered_df["FusionCategory"]
)
# ---------- end AGFusion sequences block ----------



#read in pvacfuse classi file to merge
pvacfuse_file_mhc_i = glob.glob(
    f"{args.WB}{final_result}/pVACfuse/mhc_i/*.all_epitopes.aggregated.tsv"
)[0]
pvacfuse_df_mhc_i = pd.read_csv(
    pvacfuse_file_mhc_i,
    delimiter="\t",
)

# Pull out the fusion key from the pvacfuse file
pvacfuse_df_mhc_i["FusionKey"] = pvacfuse_df_mhc_i["ID"].str.replace(r"^\d+\.(.*?\.(?:frameshift_fusion|inframe_fusion))\.\d+$", r"\1", regex=True)


# Step 4: Merge on both keys
merged_class_i_df = pd.merge(
    pvacfuse_df_mhc_i,
    merged_filtered_df,
    on='FusionKey',
    how="inner"  # or "inner" if you want only matches
)

cols_to_drop = [
    "JunctionReads",
    "SpanningFrags",
    "CounterFusionLeftReads",
    "CounterFusionRightReads",
]

for c in cols_to_drop:
    if c in merged_class_i_df.columns:
        merged_class_i_df.drop(columns=[c], inplace=True)

merged_class_i_df.to_csv(
    f"{fusion_dir}/tumor-exome.all_epitopes.aggregated.mhc_i.review.tsv",
    sep='\t',
    index=False
)

#read in pvacfuse classii file to merge
pvacfuse_file_mhc_ii = glob.glob(
    f"{args.WB}{final_result}/pVACfuse/mhc_ii/*.all_epitopes.aggregated.tsv"
)[0]

pvacfuse_df_mhc_ii = pd.read_csv(
    pvacfuse_file_mhc_ii,
    delimiter="\t",
)

# Pull out the fusion key from the pvacfuse file
pvacfuse_df_mhc_ii["FusionKey"] = pvacfuse_df_mhc_ii["ID"].str.replace(r"^\d+\.(.*?\.(?:frameshift_fusion|inframe_fusion))\.\d+$", r"\1", regex=True)


# Step 4: Merge on both keys
merged_class_ii_df = pd.merge(
    pvacfuse_df_mhc_ii,
    merged_filtered_df,
    on='FusionKey',
    how="inner"  # or "inner" if you want only matches
)

for c in cols_to_drop:
    if c in merged_class_ii_df.columns:
        merged_class_ii_df.drop(columns=[c], inplace=True)

merged_class_ii_df.to_csv(
    f"{fusion_dir}/tumor-exome.all_epitopes.aggregated.mhc_ii.review.tsv",
    sep='\t',
    index=False
)