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
        help="the path to the gcp_immuno folder of the trial you wish to tun script on, defined as WORKING_BASE in envs.txt",
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

# Step 2: Create a concatenated transcript field for matching
merged_filtered_df["TranscriptKey"] = merged_filtered_df["CDS_LEFT_ID"].astype(str) + "_" + merged_filtered_df["CDS_RIGHT_ID"].astype(str)

merged_filtered_df["FusionCategory"] = merged_filtered_df["PROT_FUSION_TYPE"].map({
    "FRAMESHIFT": "frameshift_fusion",
    "INFRAME": "inframe_fusion"
})

merged_filtered_df["FusionKey"] = merged_filtered_df["FusionPartners"] + '.' + merged_filtered_df["TranscriptKey"] + '.' + merged_filtered_df["FusionCategory"]

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

merged_class_i_df.to_csv(
    f"{fusion_dir}/tumor-exome.all_epitopes.aggregated.mhc_i.review.csv", index=False
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

merged_class_ii_df.to_csv(
    f"{fusion_dir}/tumor-exome.all_epitopes.aggregated.mhc_ii.review.csv", index=False
)

# condition = (
#     pvacfuse_df_mhc_i["Gene"].str.replace("_", "--").isin(passed_candidates)
# ) & (pvacfuse_df_mhc_i["Num Passing Peptides"] > 0)
# pvacfuse_df_mhc_i["Tier"] = pvacfuse_df_mhc_i["Tier"].astype("string")

# # Set 'Status' column to 'Fusion Detected' where condition is True
# pvacfuse_df_mhc_i.loc[condition, "Tier"] = "Review"

# # Set 'Status' column to 'No Fusion' where condition is False
# pvacfuse_df_mhc_i.loc[~condition, "Tier"] = "Poor"

# pvacfuse_df_mhc_i.to_csv(
#     f"{fusion_dir}/tumor-exome.all_epitopes.aggregated.mhc_i.review.csv", index=False
# )


# # Filter the DataFrame for Class II
# pvacfuse_file_mhc_ii = glob.glob(
#     f"{args.WB}{final_result}/pVACfuse/mhc_ii/*.all_epitopes.aggregated.tsv"
# )[0]
# pvacfuse_df_mhc_ii = fusion_candidates_df = pd.read_csv(
#     pvacfuse_file_mhc_ii,
#     delimiter="\t",
# )

# condition = (
#     pvacfuse_df_mhc_ii["Gene"].str.replace("_", "--").isin(passed_candidates)
# ) & (pvacfuse_df_mhc_ii["Num Passing Peptides"] > 0)
# pvacfuse_df_mhc_ii["Tier"] = pvacfuse_df_mhc_ii["Tier"].astype("string")

# # Set 'Status' column to 'Fusion Detected' where condition is True
# pvacfuse_df_mhc_ii.loc[condition, "Tier"] = "Review"

# # Set 'Status' column to 'No Fusion' where condition is False
# pvacfuse_df_mhc_ii.loc[~condition, "Tier"] = "Poor"

# pvacfuse_df_mhc_ii.to_csv(
#     f"{fusion_dir}/tumor-exome.all_epitopes.aggregated.mhc_ii.review.csv", index=False
# )
