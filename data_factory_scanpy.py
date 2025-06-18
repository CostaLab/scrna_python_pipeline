# We load the libraries.
import numpy as np
import pandas as pd
import scanpy as sc
import scanpy.external as sce
import seaborn as sns
import harmonypy as hm
import anndata as ad
import re
import gc
import scrublet as scr
import skmisc
import matplotlib.pyplot as plt
import sys
import os
sc.logging.print_header()
sc.settings.set_figure_params(dpi = 300, facecolor = "white")


# Defining custom functions.
def check_numeric(value, var_name = "value"):
    if not isinstance(value, (int, float)):
        sys.exit(f"Error: {var_name} must be a number. Got {type(value).__name__} instead.")


def check_numeric_or_inf(value, var_name = "value"):
    if isinstance(value, str) and value == "Inf":
        return float("inf")
    elif isinstance(value, (int, float)):
        return value
    else:
        sys.exit(f"Error: {var_name} must be a number or 'Inf' (string). Got {type(value).__name__} instead.")


if len(sys.argv) != 2:
    raise ValueError("Usage: python data_factory_scanpy.py /path/to/config.py")


# Where is the config file located?
config_path = sys.argv[1]


if not os.path.isfile(config_path):
    raise FileNotFoundError(f"Config file not found: {config_path}")


print("Loading the config file.")
exec(open(config_path).read())
sc.settings.verbosity = SETTING_VERBOSITY


print("Generating the output folder if needed.")
os.makedirs(output, exist_ok = True)


# Checking the parameters.
check_numeric(MINCELLS, "MINCELLS")
check_numeric(MINGENES, "MINGENES")
check_numeric(MINREADS, "MINREADS")

MAXGENES         = check_numeric_or_inf(MAXGENES, "MAXGENES")
MAXREADS         = check_numeric_or_inf(MAXREADS, "MAXREADS")
PCT_MITO_CEILING = check_numeric_or_inf(PCT_MITO_CEILING, "PCT_MITO_CEILING")
PCT_MITO_FLOOR   = check_numeric_or_inf(PCT_MITO_FLOOR, "PCT_MITO_FLOOR")
PCT_RIBO_CEILING = check_numeric_or_inf(PCT_RIBO_CEILING, "PCT_RIBO_CEILING")
PCT_RIBO_FLOOR   = check_numeric_or_inf(PCT_RIBO_FLOOR, "PCT_RIBO_FLOOR")


print("We get the list of sample names.")
sample_names = list(stage_lst.keys())


print("Loading the data.")
scrna = sc.read_10x_mtx(data_src[sample_names[0]], cache = True)
if(len(sample_names) > 1):
    scrna = [scrna]
    for x in sample_names[1:]:
        print(x)
        scrna.append(sc.read_10x_mtx(data_src[x], cache = True))


scrna_dict = {name: adata for name, adata in zip(sample_names, scrna)}
scrna = ad.concat(scrna, label = "batch", keys = sample_names, index_unique = "-")
# Changing the names to have the form [SAMPLE]_[CELL_ID]-1.
# This is done, because of tradition.
scrna.obs_names = [
    f"{batch}_{cell_id.rsplit('-', 1)[0]}"
    for batch, cell_id in zip(scrna.obs["batch"], scrna.obs_names)
]


scrna.obs["stage"] = [stage_lst[x] for x in list(scrna.obs["batch"])]
scrna.write_h5ad(filename = output+"/scrna_raw_data.h5ad")
scrna = sc.read_h5ad(output+"/scrna_raw_data.h5ad")#, backed = "r+") # Backed is less memory intensive, but does not work. The matrix cannot be read by the QC functions.


print("Looking into the total counts per condition and the number of genes per counts.")
scrna.var["mito"] = scrna.var_names.str.upper().str.startswith("MT-")
scrna.var["ribo"] = scrna.var_names.str.upper().str.match(r"^RP[SL]")
sc.pp.calculate_qc_metrics(scrna, qc_vars = ["mito", "ribo"], percent_top = None, log1p = False, inplace = True)


qc_dir = os.path.join(output, "qc_plots")
os.makedirs(qc_dir, exist_ok = True)
qc_vars = ["total_counts", "n_genes_by_counts", "pct_counts_mito", "pct_counts_ribo"]
y_labels = {
    "total_counts": "Number of Reads",
    "n_genes_by_counts": "Number of Genes",
    "pct_counts_mito": "Mitochondrial Reads (%)",
    "pct_counts_ribo": "Ribosomal Reads (%)"
}


print("Generating QC figures before quality control.")
n_samples = len(scrna.obs["batch"].unique())
fig_width = max(6, n_samples)
for var in qc_vars:
    plt.figure(figsize = (fig_width, 5))
    ax = sns.violinplot(x = "batch", y = var, data = scrna.obs, density_norm = "width", inner = "box", cut = 0)
    ax.set_xlabel("Sample")
    ax.set_ylabel(y_labels.get(var, var))
    ax.set_title(f"{y_labels.get(var, var)} by Sample")
    plt.xticks(rotation = 45, ha = "right")
    plt.tight_layout()
    pdf_path = os.path.join(qc_dir, f"{var}_violin_raw.pdf")
    plt.savefig(pdf_path)
    png_path = os.path.join(qc_dir, f"{var}_violin_raw.png")
    plt.savefig(png_path, dpi = 300)
    plt.close()


print("Filtering genes and cells.")
sc.pp.filter_cells(scrna, min_genes = MINGENES)
sc.pp.filter_cells(scrna, max_genes = MAXGENES)
sc.pp.filter_genes(scrna, min_cells = MINCELLS)
scrna = scrna[scrna.obs.total_counts < MAXREADS, :]
scrna = scrna[scrna.obs.total_counts > MINREADS, :]
scrna = scrna[scrna.obs.pct_counts_mito < PCT_MITO_CEILING, :]
scrna = scrna[scrna.obs.pct_counts_mito > PCT_MITO_FLOOR, :]
scrna = scrna[scrna.obs.pct_counts_ribo < PCT_RIBO_CEILING, :]
scrna = scrna[scrna.obs.pct_counts_ribo > PCT_RIBO_FLOOR, :]


print("Detecting doublets.")
scrna.layers["counts"] = scrna.X.copy()
doublet_scores = []
predicted_doublets = []
for batch in scrna.obs["batch"].unique():
    adata_batch = scrna[scrna.obs["batch"] == batch]
    counts_matrix = adata_batch.layers["counts"].toarray() if hasattr(adata_batch.layers["counts"], "toarray") else adata_batch.layers["counts"]
    scrub = scr.Scrublet(counts_matrix)
    scores, preds = scrub.scrub_doublets()
    scrna.obs.loc[adata_batch.obs_names, "scrublet_score"] = scores
    scrna.obs.loc[adata_batch.obs_names, "scrublet_predicted_doublet"] = preds


print("Plotting the number of doublets per sample.")
doublet_counts = (
    scrna.obs
    .groupby("batch")["scrublet_predicted_doublet"]
    .sum()
    .sort_index()
)
plt.figure(figsize=(fig_width, 5))
sns.barplot(x = doublet_counts.index, y = doublet_counts.values, palette = "Blues_d")
plt.ylabel("Number of Suspected Doublets")
plt.xlabel("Sample")
plt.title("Suspected Doublets per Sample")
plt.xticks(rotation = 45, ha = "right")
plt.tight_layout()
plt.savefig(os.path.join(qc_dir, "doublet_barplot_filtered.pdf"))
plt.savefig(os.path.join(qc_dir, "doublet_barplot_filtered.png"), dpi = 300)
plt.close()


print("Removing doublets if so desired.")
if doublet_switch:
    print("Removing doublets.")
    scrna = scrna[scrna.obs["scrublet_predicted_doublet"] == False, :]


print("Generating QC figures after quality control.")
qc_vars = ["total_counts", "n_genes_by_counts", "pct_counts_mito", "pct_counts_ribo", "scrublet_score"]
for var in qc_vars:
    plt.figure(figsize = (fig_width, 5))
    ax = sns.violinplot( x = "batch", y = var, a = scrna.obs, density_norm = "width", inner = "box", cut = 0)
    ax.set_xlabel("Sample")
    ax.set_ylabel(y_labels.get(var, var))
    ax.set_title(f"{y_labels.get(var, var)} by Sample")
    plt.xticks(rotation = 45, ha = "right")
    plt.tight_layout()
    pdf_path = os.path.join(qc_dir, f"{var}_violin_filtered.pdf")
    plt.savefig(pdf_path)
    png_path = os.path.join(qc_dir, f"{var}_violin_filtered.png")
    plt.savefig(png_path, dpi = 300)
    plt.close()


print("Log normalizing the data.")
sc.pp.normalize_total(scrna, target_sum = 1e4)
sc.pp.log1p(scrna)
scrna.layers["lognorm"] = scrna.X.copy()


print("Cell Cycle Scoring.")
# To reduce the number of files required to run the analysis, the cell cycle genes are saved here.
# They are from Regev Lab (regev_lab_cell_cycle_genes.txt)
# Right now, they are only for humans. An implementation for mouse is necessary.
cell_cycle_genes = ["MCM5", "PCNA", "TYMS", "FEN1", "MCM2", "MCM4", "RRM1", "UNG", "GINS2", "MCM6", "CDCA7", "DTL", "PRIM1", "UHRF1", "MLF1IP", "HELLS", "RFC2", "RPA2", "NASP", "RAD51AP1", "GMNN", "WDR76", "SLBP", "CCNE2", "UBR7",
                    "POLD3", "MSH2", "ATAD2", "RAD51", "RRM2", "CDC45", "CDC6", "EXO1", "TIPIN", "DSCC1", "BLM", "CASP8AP2", "USP1", "CLSPN", "POLA1", "CHAF1B", "BRIP1", "E2F8", "HMGB2", "CDK1", "NUSAP1", "UBE2C", "BIRC5", "TPX2",
                    "TOP2A", "NDC80", "CKS2", "NUF2", "CKS1B", "MKI67", "TMPO", "CENPF", "TACC3", "FAM64A", "SMC4", "CCNB2", "CKAP2L", "CKAP2", "AURKB", "BUB1", "KIF11", "ANP32E", "TUBB4B", "GTSE1", "KIF20B", "HJURP", "CDCA3", "HN1",
                    "CDC20", "TTK", "CDC25C", "KIF2C", "RANGAP1", "NCAPD2", "DLGAP5", "CDCA2", "CDCA8", "ECT2", "KIF23", "HMMR", "AURKA", "PSRC1", "ANLN", "LBR", "CKAP5", "CENPE", "CTCF", "NEK2", "G2E3", "GAS2L3", "CBX5", "CENPA"]
s_genes = ["MCM5", "PCNA", "TYMS", "FEN1", "MCM2", "MCM4", "RRM1", "UNG", "GINS2", "MCM6", "CDCA7", "DTL", "PRIM1", "UHRF1", "MLF1IP", "HELLS", "RFC2", "RPA2", "NASP", "RAD51AP1", "GMNN", "WDR76", "SLBP", "CCNE2", "UBR7", "POLD3",
           "MSH2", "ATAD2", "RAD51", "RRM2", "CDC45", "CDC6", "EXO1", "TIPIN", "DSCC1", "BLM", "CASP8AP2", "USP1", "CLSPN", "POLA1", "CHAF1B", "BRIP1", "E2F8"]
g2m_genes = ["HMGB2", "CDK1", "NUSAP1", "UBE2C", "BIRC5", "TPX2", "TOP2A", "NDC80", "CKS2", "NUF2", "CKS1B", "MKI67", "TMPO", "CENPF", "TACC3", "FAM64A", "SMC4", "CCNB2", "CKAP2L", "CKAP2", "AURKB", "BUB1", "KIF11", "ANP32E", "TUBB4B",
             "GTSE1", "KIF20B", "HJURP", "CDCA3", "HN1", "CDC20", "TTK", "CDC25C", "KIF2C", "RANGAP1", "NCAPD2", "DLGAP5", "CDCA2", "CDCA8", "ECT2", "KIF23", "HMMR", "AURKA", "PSRC1", "ANLN", "LBR", "CKAP5", "CENPE", "CTCF", "NEK2",
             "G2E3", "GAS2L3", "CBX5", "CENPA"]
cell_cycle_genes = [x for x in cell_cycle_genes if x in scrna.var_names]
s_genes = [x for x in s_genes if x in scrna.var_names]
g2m_genes = [x for x in g2m_genes if x in scrna.var_names]
sc.tl.score_genes_cell_cycle(scrna, s_genes = s_genes, g2m_genes = g2m_genes)


print("Determining highly variable genes.")
# We only scale on the highly variable genes to reduce the amount of memory needed.
sc.pp.highly_variable_genes(scrna, n_top_genes = 2000, flavor = "seurat_v3", subset = False, layer = "counts")
# Generating a second object to reduce the memory needed in the integration.
scrna_hvgs = scrna[:, scrna.var["highly_variable"]].copy()


print("Regressing out confounders.")
sc.pp.regress_out(scrna_hvgs, keys = ["total_counts", "pct_counts_mito", "S_score", "G2M_score"], n_jobs = WORKER_NUM)
sc.pp.scale(scrna_hvgs)


print("Dimension Reduction.")
sc.tl.pca(scrna_hvgs, svd_solver = "arpack")


print("Looking into the PC variance contribution.")
sc.pl.pca_variance_ratio(scrna_hvgs, log = True, n_pcs = 50)


print("PCA.")
sc.pp.neighbors(scrna_hvgs, n_neighbors = 10, n_pcs = 40, use_rep = "X_pca")
sc.tl.umap(scrna_hvgs)
scrna_hvgs.obsm["X_pca_umap"] = scrna_hvgs.obsm["X_umap"] 


print("Sample Integration.")
sce.pp.harmony_integrate(scrna_hvgs, key = "batch")
sc.pp.neighbors(scrna_hvgs, n_neighbors = 10, n_pcs = 40, use_rep = "X_pca_harmony")
sc.tl.umap(scrna_hvgs)
scrna_hvgs.obsm["X_umap_harmony"] = scrna_hvgs.obsm["X_umap"] 


print("Adding the integration information to the full object.")
scrna.obsm["X_pca"]          = scrna_hvgs.obsm["X_pca"]
scrna.obsm["X_umap"]         = scrna_hvgs.obsm["X_pca_umap"]
scrna.obsm["X_umap_harmony"] = scrna_hvgs.obsm["X_umap"]
scrna.obsm["X_pca_harmony"]  = scrna_hvgs.obsm["X_pca_harmony"]
scrna.obsp["distances"]      = scrna_hvgs.obsp["distances"]
scrna.obsp["connectivities"] = scrna_hvgs.obsp["connectivities"]
scrna.uns["neighbors"]       = scrna_hvgs.uns["neighbors"]
scrna.uns["umap"]            = scrna_hvgs.uns["umap"]
scrna.layers["scaled_hvg"]   = scrna_hvgs.X.copy()
scrna.varm["PCs"]            = scrna_hvgs.varm["PCs"]


print("Saving the preprocessed data.")
scrna = sc.read_h5ad(output+"/scrna_preprocessed_data.h5ad")


print("Data Clustering.")
clustering_dir = os.path.join(output, "clustering_plots")
os.makedirs(clustering_dir, exist_ok = True)
for res in np.arange(0.1, 0.9, 0.1):
    key_added = f"leiden_{res:.1f}"
    print(key_added)
    sc.tl.leiden(scrna, resolution = res, key_added = key_added)
    fig = sc.pl.umap(scrna, color = key_added, basis = "X_umap_harmony", show = False, return_fig = True)
    fig.savefig(os.path.join(clustering_dir, f"umap_{key_added}.pdf"))
    fig.savefig(os.path.join(clustering_dir, f"umap_{key_added}.png"), dpi = 300)
    plt.close(fig)
    deg_dir = os.path.join(output, f"degs_{key_added}")
    os.makedirs(deg_dir, exist_ok = True)
    writer = pd.ExcelWriter(os.path.join(deg_dir, f"{key_added}_degs.xlsx"), engine = "xlsxwriter")
    for cluster in sorted(scrna.obs[key_added].unique()):
        print("Cluster: "+str(cluster))
        sc.tl.rank_genes_groups(scrna, groupby = key_added, groups = [cluster], reference = "rest", method = DEG_METHOD, key_added = "degs_genes_temp")
        result = sc.get.rank_genes_groups_df(scrna, group = cluster, key = "degs_genes_temp")
        result.to_excel(writer, sheet_name = f"Cluster_{cluster}", index = False)
        top_up = result[result["logfoldchanges"] > 0].nlargest(10, "logfoldchanges")
        top_down = result[result["logfoldchanges"] < 0].nsmallest(10, "logfoldchanges")
        top_genes = pd.concat([top_down[::-1], top_up])
        plt.figure(figsize = (6, 6))
        bar_colors = ["red"] * len(top_down) + ["blue"] * len(top_up)
        sns.barplot(x = "logfoldchanges", y = "names", data = top_genes, palette = bar_colors)
        plt.axvline(0, color = "gray", linestyle = "--")
        plt.title(f"Top DEGs for Cluster {cluster} ({key_added})")
        plt.tight_layout()
        plt.savefig(os.path.join(deg_dir, f"degs_{key_added}_cluster_{cluster}.pdf"))
        plt.savefig(os.path.join(deg_dir, f"degs_{key_added}_cluster_{cluster}.png"), dpi = 300)
        plt.close()
    writer.close()


print("We compare the stages.")
stage_dir = os.path.join(output, "degs_stage")
os.makedirs(stage_dir, exist_ok = True)
sc.tl.rank_genes_groups(scrna, groupby = "stage", method = DEG_METHOD, key_added = "deg_genes_stage")


stages = scrna.obs["stage"].unique().tolist()
writer = pd.ExcelWriter(os.path.join(stage_dir, "stage_degs.xlsx"), engine = "xlsxwriter")
for stage in stages:
    result = sc.get.rank_genes_groups_df(scrna, group = stage, key = "deg_genes_stage")
    result.to_excel(writer, sheet_name = f"Stage_{stage}", index = False)
    top_up = result[result["logfoldchanges"] > 0].nlargest(10, "logfoldchanges")
    top_down = result[result["logfoldchanges"] < 0].nsmallest(10, "logfoldchanges")
    top_genes = pd.concat([top_down[::-1], top_up])
    plt.figure(figsize = (6, 6))
    bar_colors = ["red"] * len(top_down) + ["blue"] * len(top_up)
    sns.barplot(x = "logfoldchanges", y = "names", data = top_genes, palette = bar_colors)
    plt.axvline(0, color = "gray", linestyle = "--")
    plt.title(f"Top DEGs for Stage {stage}")
    plt.tight_layout()
    plt.savefig(os.path.join(stage_dir, f"degs_stage_{stage}.pdf"))
    plt.savefig(os.path.join(stage_dir, f"degs_stage_{stage}.png"), dpi = 300)
    plt.close()


writer.close()

scrna.write_h5ad(os.path.join(output, "scrna_final_data.h5ad"))
