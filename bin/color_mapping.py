import pandas as pd
import os

# load monocle metadata
monocle = pd.read_csv("/data/pam/team230/sm71/scratch/gps_project/metadata/results_combined.csv", low_memory=False)

# report and drop rows with no GPSC
nan_gpsc = monocle["GPSC"].isna().sum()
monocle = monocle.dropna(subset=["GPSC"])

# for samples with multiple GPSC assignments (e.g. "1215;5"), pick the smallest numeric value.
# smaller = current canonical label under v11 (higher labels were merged into lower ones).
multi_gpsc = monocle["GPSC"].astype(str).str.contains(";")
multi_gpsc_df = monocle.loc[multi_gpsc, ["Sample_ID", "GPSC"]].copy()
monocle["GPSC"] = (
    monocle["GPSC"]
    .astype(str)
    .str.split(";")
    .apply(lambda parts: min(int(p.strip()) for p in parts))
)

# sort by GPSC number, then Sample_ID lexicographically within each group
monocle = monocle.sort_values(["GPSC", "Sample_ID"])

# get assembly paths from disk
assembly_dir = "/data/pam/team230/sm71/scratch/gps_project/assemblies/"
assembly_files = set(os.listdir(assembly_dir))

# order assembly paths based on sorted monocle dataframe, only include samples with an assembly on disk
monocle["filename"] = monocle["Sample_ID"] + ".contigs.fasta"
no_assembly = ~monocle["filename"].isin(assembly_files)
on_disk_no_metadata = assembly_files - set(monocle["filename"])
monocle = monocle[monocle["filename"].isin(assembly_files)].reset_index(drop=True)
monocle["file_path"] = monocle["filename"].apply(lambda f: os.path.join(assembly_dir, f))

# write file_colors.txt — sorted paths for --file-colors flag
monocle["file_path"].to_csv(
    "/data/pam/team230/sm71/scratch/gps_project/themisto2/file_colors_input.txt",
    index=False, header=False
)

# write gpsc_mapping.tsv — one row per Sample_ID, multiple GPSCs semicolon-separated
gpsc_mapping = (
    monocle.groupby("Sample_ID", sort=False)["GPSC"]
    .apply(lambda x: ";".join(x.astype(str)))
    .reset_index()
)
gpsc_mapping.to_csv(
    "/data/pam/team230/sm71/scratch/gps_project/themisto2/gpsc_mapping.tsv",
    index=False, sep="\t"
)

lines = [
    f"Samples in metadata with no assembly on disk (dropped): {no_assembly.sum()}",
    f"Assemblies on disk with no metadata entry (excluded): {len(on_disk_no_metadata)}",
    f"Total assemblies written: {len(monocle)}",
    f"Samples with GPSC labels: {len(gpsc_mapping)}",
    f"Samples with no GPSC assignment (dropped): {nan_gpsc}",
    f"Samples with multiple GPSC assignments (resolved to smallest): {len(multi_gpsc_df)}",
]
if not multi_gpsc_df.empty:
    resolved = multi_gpsc_df.copy()
    resolved["GPSC_resolved"] = resolved["GPSC"].str.split(";").apply(
        lambda parts: min(int(p.strip()) for p in parts)
    )
    lines.append(resolved.to_string(index=False))

stats_path = "/data/pam/team230/sm71/scratch/gps_project/themisto2/stats.txt"
with open(stats_path, "w") as f:
    f.write("\n".join(lines) + "\n")

for line in lines:
    print(line)
