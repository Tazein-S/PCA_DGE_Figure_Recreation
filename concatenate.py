files = glob.glob("data/*.readCounts.txt")
dfs = []

for f in files:
    # sample name from filename
    name = os.path.basename(f)
    sample = re.sub(r"\.readCounts\.txt$", "", name)
    sample = sample.split("_", 1)[1]   # optional: remove GSM prefix
    
    # load file, robust to spaces or tabs
    df = pd.read_csv(
        f, 
        sep=r"\s+",          # handles tab OR spaces
        header=None, 
        names=["gene", sample],
        engine="python"      # needed for regex separators
    )
    
    dfs.append(df)

# Merge all files by gene
merged = dfs[0]
for df in dfs[1:]:
    merged = merged.merge(df, on="gene", how="outer")

# Replace NaNs with zeros
merged = merged.fillna(0)

merged.to_csv("merged_counts_matrix.txt", sep="\t", index=False)
