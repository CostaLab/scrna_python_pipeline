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
# import multiprocessing as mp
from collections import defaultdict
import glob
sc.logging.print_header()
sc.settings.set_figure_params(dpi = 300, facecolor = "white")


# Defining custom functions.
def check_numeric(value, var_name = "value"):
    if not isinstance(value, (int, float)):
        sys.exit(f"Error: {var_name} must be a number. Got {type(value).__name__} instead.")


def check_strings(value, allowed_values, var_name = "value"):
    if not isinstance(value, str):
        sys.exit(f"Error: {var_name} must be a string. Got {type(value).__name__} instead.")
    if not value in allowed_values:
        sys.exit(f"Error: {var_name} must be one of {allowed_values}.")


def check_boolean(value, var_name = "value"):
    if not isinstance(value, (bool)):
        sys.exit(f"Error: phase {var_name} must be logical. Got {type(value).__name__} instead.")


def check_numeric_or_inf(value, var_name = "value"):
    if isinstance(value, str) and value == "Inf":
        return float("inf")
    elif isinstance(value, (int, float)):
        return value
    else:
        sys.exit(f"Error: {var_name} must be a number or 'Inf' (string). Got {type(value).__name__} instead.")


def process_cluster_task(args):
    res, cluster, key_added, output, deg_method, max_iter = args
    deg_dir = os.path.join(output, f"degs_{key_added}")
    os.makedirs(deg_dir, exist_ok = True)
    
    print(f"Running: {key_added}, Cluster {cluster}", flush = True)
    
    adata_cluster = scrna.copy()
    
    # Run DEG
    sc.tl.rank_genes_groups(adata_cluster, groupby = key_added, groups = [cluster],
                            reference = "rest", method = deg_method,
                            max_iter = max_iter,
                            key_added = f"degs_{key_added}_cluster_{cluster}")
    result = sc.get.rank_genes_groups_df(adata_cluster, group = cluster,
                                         key = f"degs_{key_added}_cluster_{cluster}")
    
    # Save CSV
    csv_path = os.path.join(deg_dir, f"degs_{key_added}_cluster_{cluster}.csv")
    result.to_csv(csv_path, index=False)
    
    # Plot top genes
    top_up = result[result["logfoldchanges"] > 0].nlargest(10, "logfoldchanges")
    top_down = result[result["logfoldchanges"] < 0].nsmallest(10, "logfoldchanges")
    top_genes = pd.concat([top_down[::-1], top_up])
    plt.figure(figsize=(6, 6))
    bar_colors = ["red"] * len(top_down) + ["blue"] * len(top_up)
    sns.barplot(x = "logfoldchanges", y="names", data = top_genes, palette = bar_colors)
    plt.axvline(0, color = "gray", linestyle = "--")
    plt.title(f"Top DEGs for Cluster {cluster} ({key_added})")
    plt.tight_layout()
    plt.savefig(os.path.join(deg_dir, f"degs_{key_added}_cluster_{cluster}.pdf"), bbox_inches = "tight")
    plt.savefig(os.path.join(deg_dir, f"degs_{key_added}_cluster_{cluster}.png"), bbox_inches = "tight", dpi=300)
    plt.close()
    
    return key_added, cluster, csv_path


if len(sys.argv) != 2:
    raise ValueError("Usage: python data_factory_scanpy.py /path/to/config.py")


# Where is the config file located?
config_path = sys.argv[1]


if not os.path.isfile(config_path):
    raise FileNotFoundError(f"Config file not found: {config_path}.")


print("Loading the config file.", flush = True)
exec(open(config_path).read())
sc.settings.verbosity = SETTING_VERBOSITY


print("Generating the output folder if needed.", flush = True)
os.makedirs(output, exist_ok = True)
qc_dir = os.path.join(output, "qc_plots")
os.makedirs(qc_dir, exist_ok = True)
clustering_dir = os.path.join(output, "clustering_plots")
os.makedirs(clustering_dir, exist_ok = True)


# Checking the parameters.
check_numeric(MINCELLS, "MINCELLS")
check_numeric(MINGENES, "MINGENES")
check_numeric(MINREADS, "MINREADS")
check_numeric(LOGREG_MAXITER, "LOGREG_MAXITER")
check_numeric(N_PCS_VARIANCE_CONTRIBUTION, "N_PCS_VARIANCE_CONTRIBUTION")
check_numeric(N_NEIGHBORS_PCA, "N_NEIGHBORS_PCA")
check_numeric(N_PCS_PCA, "N_PCS_PCA")
check_numeric(N_TOP_GENES_INTEGRATION, "N_TOP_GENES_INTEGRATION")
check_numeric(MAX_ITER_HARMONY, "MAX_ITER_HARMONY")
check_numeric(N_NEIGHBORS_UMAP_HARMONY, "N_NEIGHBORS_UMAP_HARMONY")
check_numeric(N_PCS_UMAP_HARMONY, "N_PCS_UMAP_HARMONY")

MAXGENES         = check_numeric_or_inf(MAXGENES, "MAXGENES")
MAXREADS         = check_numeric_or_inf(MAXREADS, "MAXREADS")
PCT_MITO_CEILING = check_numeric_or_inf(PCT_MITO_CEILING, "PCT_MITO_CEILING")
PCT_MITO_FLOOR   = check_numeric_or_inf(PCT_MITO_FLOOR, "PCT_MITO_FLOOR")
PCT_RIBO_CEILING = check_numeric_or_inf(PCT_RIBO_CEILING, "PCT_RIBO_CEILING")
PCT_RIBO_FLOOR   = check_numeric_or_inf(PCT_RIBO_FLOOR, "PCT_RIBO_FLOOR")


if(len(phases) == 0):
    raise ValueError("You are not supplying any phases.")


for key, val in phases.items():
    check_strings(val, allowed_values = ("execute", "load", "skip"), var_name = key)


# We check if a phase has a required phase.
required_phases_for_phases = {
    "filter": ["integration"],
    "integration": ["cluster"],
    "cluster": ["comparison"]
}

for phase, dependents in required_phases_for_phases.items():
    if phases.get(phase) == "skip":
        for dependent_use in dependents:
            if phases.get(dependent_use) == "execute":
                raise ValueError(f"Cannot execute '{dependent_use}' phase without loading or executing '{phase}' phase.")


print("We get the list of sample names.", flush = True)
sample_names = list(stage_lst.keys())




if phases.get("raw") == "skip":
    print("Skipping data loading.", flush = True)
elif phases.get("raw") == "load":
    print("Loading the raw data.", flush = True)
    if not os.path.isfile(os.path.join(output, "scrna_raw_data.h5ad")):
        raise FileNotFoundError(f"scrna_raw_data not found: {output}scrna_raw_data.h5ad.")
    scrna = sc.read_h5ad(os.path.join(output, "scrna_raw_data.h5ad"))#, backed = "r+") # Backed is less memory intensive, but does not work. The matrix cannot be read by the QC functions.
elif phases.get("raw") == "execute":
    print("Loading the data.", flush = True)
    scrna = sc.read_10x_mtx(data_src[sample_names[0]], cache = True)
    if(len(sample_names) > 1):
        scrna = [scrna]
        for x in sample_names[1:]:
            print(x, flush = True)
            scrna.append(sc.read_10x_mtx(data_src[x], cache = True))
        scrna_dict = {name: adata for name, adata in zip(sample_names, scrna)}
        scrna = ad.concat(scrna, label = "batch", keys = sample_names, index_unique = "-")
    else:
        scrna.obs["batch"] = sample_names[0]
    
    
    # Changing the names to have the form [SAMPLE]_[CELL_ID]-1.
    # This is done, because of tradition.
    scrna.obs_names = [
        f"{batch}_{cell_id.rsplit('-', 1)[0]}"
        for batch, cell_id in zip(scrna.obs["batch"], scrna.obs_names)
    ]
    
    
    scrna.obs["stage"] = [stage_lst[x] for x in list(scrna.obs["batch"])]
    print("Saving the raw data.", flush = True)
    scrna.write_h5ad(filename = os.path.join(output, "scrna_raw_data.h5ad"))




if phases.get("filter") == "skip":
    print("Skipping data filtering.", flush = True)
elif phases.get("filter") == "load":
    print("Loading the filtered data.", flush = True)
    if not os.path.isfile(os.path.join(output, "scrna_filtered_data.h5ad")):
        raise FileNotFoundError(f"scrna_filtered_data not found: {output}scrna_filtered_data.h5ad.")
    scrna = sc.read_h5ad(os.path.join(output, "scrna_filtered_data.h5ad"))
    scrna.obs["scrublet_predicted_doublet"] = scrna.obs["scrublet_predicted_doublet_str"] == "True"
    del scrna.obs["scrublet_predicted_doublet_str"]
elif phases.get("filter") == "execute":
    print("Looking into the total counts per condition and the number of genes per counts.", flush = True)
    scrna.var["mito"] = scrna.var_names.str.upper().str.startswith("MT-")
    scrna.var["ribo"] = scrna.var_names.str.upper().str.match(r"^RP[SL]")
    sc.pp.calculate_qc_metrics(scrna, qc_vars = ["mito", "ribo"], percent_top = None, log1p = False, inplace = True)
    
    
    qc_vars = ["total_counts", "n_genes_by_counts", "pct_counts_mito", "pct_counts_ribo"]
    y_labels = {
        "total_counts": "Number of Reads",
        "n_genes_by_counts": "Number of Genes",
        "pct_counts_mito": "Mitochondrial Reads (%)",
        "pct_counts_ribo": "Ribosomal Reads (%)"
    }
    
    
    print("Generating QC figures before quality control.", flush = True)
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
        plt.savefig(pdf_path, bbox_inches = "tight")
        png_path = os.path.join(qc_dir, f"{var}_violin_raw.png")
        plt.savefig(png_path, bbox_inches = "tight", dpi = 300)
        plt.close()
    
    
    print("Filtering genes and cells.", flush = True)
    sc.pp.filter_cells(scrna, min_genes = MINGENES)
    sc.pp.filter_cells(scrna, max_genes = MAXGENES)
    sc.pp.filter_genes(scrna, min_cells = MINCELLS)
    scrna = scrna[scrna.obs.total_counts < MAXREADS, :]
    scrna = scrna[scrna.obs.total_counts > MINREADS, :]
    scrna = scrna[scrna.obs.pct_counts_mito < PCT_MITO_CEILING, :]
    scrna = scrna[scrna.obs.pct_counts_mito > PCT_MITO_FLOOR, :]
    scrna = scrna[scrna.obs.pct_counts_ribo < PCT_RIBO_CEILING, :]
    scrna = scrna[scrna.obs.pct_counts_ribo > PCT_RIBO_FLOOR, :]
    
    
    print("Detecting doublets.", flush = True)
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
    
    
    print("Plotting the number of doublets per sample.", flush = True)
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
    plt.savefig(os.path.join(qc_dir, "doublet_barplot_filtered.pdf"), bbox_inches = "tight")
    plt.savefig(os.path.join(qc_dir, "doublet_barplot_filtered.png"), bbox_inches = "tight", dpi = 300)
    plt.close()
    
    
    print("Removing doublets if so desired.", flush = True)
    if doublet_switch:
        print("Removing doublets.", flush = True)
        scrna = scrna[scrna.obs["scrublet_predicted_doublet"] == False, :].copy()
    
    
    print("Generating QC figures after quality control.", flush = True)
    qc_vars = ["total_counts", "n_genes_by_counts", "pct_counts_mito", "pct_counts_ribo", "scrublet_score"]
    for var in qc_vars:
        plt.figure(figsize = (fig_width, 5))
        ax = sns.violinplot(x = "batch", y = var, data = scrna.obs, density_norm = "width", inner = "box", cut = 0)
        ax.set_xlabel("Sample")
        ax.set_ylabel(y_labels.get(var, var))
        ax.set_title(f"{y_labels.get(var, var)} by Sample")
        plt.xticks(rotation = 45, ha = "right")
        plt.tight_layout()
        pdf_path = os.path.join(qc_dir, f"{var}_violin_filtered.pdf")
        plt.savefig(pdf_path, bbox_inches = "tight")
        png_path = os.path.join(qc_dir, f"{var}_violin_filtered.png")
        plt.savefig(png_path, bbox_inches = "tight", dpi = 300)
        plt.close()
    
    
    print("Log normalizing the data.", flush = True)
    sc.pp.normalize_total(scrna, target_sum = 1e4)
    sc.pp.log1p(scrna)
    scrna.layers["lognorm"] = scrna.X.copy()
    
    
    print("Cell Cycle Scoring.", flush = True)
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
    
    print("Saving the filtered data.", flush = True)
    scrna.obs["scrublet_predicted_doublet_str"] = scrna.obs["scrublet_predicted_doublet"].astype(str)
    del scrna.obs["scrublet_predicted_doublet"]
    scrna.write_h5ad(filename = os.path.join(output, "scrna_filtered_data.h5ad"))
    scrna.obs["scrublet_predicted_doublet"] = scrna.obs["scrublet_predicted_doublet_str"] == "True"
    del scrna.obs["scrublet_predicted_doublet_str"]




if phases.get("integration") == "skip":
    print("Skipping data integration.", flush = True)
elif phases.get("integration") == "load":
    print("Loading the integrated data.", flush = True)
    if not os.path.isfile(os.path.join(output, "scrna_integrated_data.h5ad")):
        raise FileNotFoundError(f"scrna_integrated_data not found: {output}scrna_integrated_data.h5ad.")
    scrna = sc.read_h5ad(filename = os.path.join(output, "scrna_integrated_data.h5ad"))
    scrna.obs["scrublet_predicted_doublet"] = scrna.obs["scrublet_predicted_doublet_str"] == "True"
    del scrna.obs["scrublet_predicted_doublet_str"]
elif phases.get("integration") == "execute":
    print("Determining highly variable genes.", flush = True)
    # We only scale on the highly variable genes to reduce the amount of memory needed.
    sc.pp.highly_variable_genes(scrna, n_top_genes = N_TOP_GENES_INTEGRATION, flavor = "seurat_v3", subset = False, layer = "counts")
    # Generating a second object to reduce the memory needed in the integration.
    scrna_hvgs = scrna[:, scrna.var["highly_variable"]].copy()
    
    
    print("Regressing out confounders.", flush = True)
    sc.pp.regress_out(scrna_hvgs, keys = ["total_counts", "pct_counts_mito", "S_score", "G2M_score"], n_jobs = WORKER_NUM)
    sc.pp.scale(scrna_hvgs)
    
    
    print("Dimension Reduction.", flush = True)
    sc.tl.pca(scrna_hvgs, svd_solver = "arpack")
    
    
    print("Looking into the PC variance contribution.", flush = True)
    plt.figure(figsize = (5, 5))
    sc.pl.pca_variance_ratio(scrna_hvgs, log = True, n_pcs = N_PCS_VARIANCE_CONTRIBUTION, show = False)
    plt.savefig(os.path.join(qc_dir, "pca_variance_ratio.png"), dpi = 300, bbox_inches = "tight")
    plt.savefig(os.path.join(qc_dir, "pca_variance_ratio.pdf"), bbox_inches = "tight")
    plt.close()
    
    
    print("PCA.", flush = True)
    sc.pp.neighbors(scrna_hvgs, n_neighbors = N_NEIGHBORS_PCA, n_pcs = N_PCS_PCA, use_rep = "X_pca")
    sc.tl.umap(scrna_hvgs)
    scrna_hvgs.obsm["X_pca_umap"] = scrna_hvgs.obsm["X_umap"] 
    
    
    print("Sample Integration.", flush = True)
    # In case of only one sample (batch), harmony_integrate simply returns the original PCA.
    # This does not cause a crash. To keep everything simple, we will ignore this fact.
    sce.pp.harmony_integrate(scrna_hvgs, key = "batch", max_iter_harmony = MAX_ITER_HARMONY)
    sc.pp.neighbors(scrna_hvgs, n_neighbors = N_NEIGHBORS_UMAP_HARMONY, n_pcs = N_PCS_UMAP_HARMONY, use_rep = "X_pca_harmony")
    sc.tl.umap(scrna_hvgs)
    scrna_hvgs.obsm["X_umap_harmony"] = scrna_hvgs.obsm["X_umap"] 
    
    
    print("Adding the integration information to the full object.", flush = True)
    scrna.obsm["X_pca"]          = scrna_hvgs.obsm["X_pca"]
    scrna.obsm["X_umap"]         = scrna_hvgs.obsm["X_pca_umap"]
    scrna.obsm["X_umap_harmony"] = scrna_hvgs.obsm["X_umap_harmony"]
    scrna.obsm["X_pca_harmony"]  = scrna_hvgs.obsm["X_pca_harmony"]
    scrna.obsp["distances"]      = scrna_hvgs.obsp["distances"]
    scrna.obsp["connectivities"] = scrna_hvgs.obsp["connectivities"]
    scrna.uns["neighbors"]       = scrna_hvgs.uns["neighbors"]
    scrna.uns["umap"]            = scrna_hvgs.uns["umap"]
    scrna.uns["scaled_hvg"]      = scrna_hvgs.X.copy()
    scrna.uns["PCs"]             = scrna_hvgs.varm["PCs"]
    del(scrna_hvgs)
    gc.collect()
    
    
    print("Generating additional plots.", flush = True)
    print("Plot: samples", flush = True)
    fig = sc.pl.embedding(scrna, basis = "X_umap_harmony", color = "batch", 
                          show = False, return_fig = True, legend_fontsize = 8)
    fig.savefig(os.path.join(clustering_dir, "umap_sample.pdf"), bbox_inches = "tight")
    fig.savefig(os.path.join(clustering_dir, "umap_sample.png"), bbox_inches = "tight", dpi = 300)
    plt.close(fig)
    
    
    print("Plot: stage", flush = True)
    fig = sc.pl.embedding(scrna, basis = "X_umap_harmony", color = "stage", 
                          show = False, return_fig = True, legend_fontsize = 8)
    fig.savefig(os.path.join(clustering_dir, "umap_stage.pdf"), bbox_inches = "tight")
    fig.savefig(os.path.join(clustering_dir, "umap_stage.png"), bbox_inches = "tight", dpi = 300)
    plt.close(fig)
    
    
    print("Plot: cell cycle", flush = True)
    fig = sc.pl.embedding(scrna, basis = "X_umap_harmony", color = "phase", 
                          show = False, return_fig = True, legend_fontsize = 8)
    fig.savefig(os.path.join(clustering_dir, "umap_phase.pdf"), bbox_inches = "tight")
    fig.savefig(os.path.join(clustering_dir, "umap_phase.png"), bbox_inches = "tight", dpi = 300)
    plt.close(fig)
    
    
    print("Saving the integrated data.", flush = True)
    scrna.obs["scrublet_predicted_doublet_str"] = scrna.obs["scrublet_predicted_doublet"].astype(str)
    del scrna.obs["scrublet_predicted_doublet"]
    scrna.write_h5ad(filename = os.path.join(output, "scrna_integrated_data.h5ad"))
    scrna.obs["scrublet_predicted_doublet"] = scrna.obs["scrublet_predicted_doublet_str"] == "True"
    del scrna.obs["scrublet_predicted_doublet_str"]




if phases.get("cluster") == "skip":
    print("Skipping data clustering.", flush = True)
elif phases.get("cluster") == "load":
    print("Loading the clustered data.", flush = True)
    if not os.path.isfile(os.path.join(output, "scrna_clustered_data.h5ad")):
        raise FileNotFoundError(f"scrna_clustered_data not found: {output}scrna_clustered_data.h5ad.")
    scrna = sc.read_h5ad(filename = os.path.join(output, "scrna_clustered_data.h5ad"))
    scrna.obs["scrublet_predicted_doublet"] = scrna.obs["scrublet_predicted_doublet_str"] == "True"
    del scrna.obs["scrublet_predicted_doublet_str"]
elif phases.get("cluster") == "execute":
    print("Data Clustering.", flush = True)
    for res in np.arange(0.1, 0.9, 0.1):
        key_added = f"leiden_{res:.1f}"
        print(key_added)
        sc.tl.leiden(scrna, resolution = res, key_added = key_added)
        # fig = sc.pl.umap(scrna, color = key_added, basis = "X_umap_harmony", show = False, return_fig = True)
        fig = sc.pl.embedding(scrna, basis = "X_umap_harmony", color = key_added, 
                              show = False, return_fig = True, legend_fontsize = 8)
        fig.savefig(os.path.join(clustering_dir, f"umap_{key_added}.pdf"), bbox_inches = "tight")
        fig.savefig(os.path.join(clustering_dir, f"umap_{key_added}.png"), bbox_inches = "tight", dpi = 300)
        plt.close(fig)
    
    
    print("Save the clustered object.", flush = True)
    scrna.obs["scrublet_predicted_doublet_str"] = scrna.obs["scrublet_predicted_doublet"].astype(str)
    del scrna.obs["scrublet_predicted_doublet"]
    scrna.write_h5ad(filename = os.path.join(output, "scrna_clustered_data.h5ad"))
    scrna.obs["scrublet_predicted_doublet"] = scrna.obs["scrublet_predicted_doublet_str"] == "True"
    del scrna.obs["scrublet_predicted_doublet_str"]




if phases.get("comparison") == "skip":
    print("Skipping data comparison.")
elif phases.get("comparison") == "load":
    print("Loading the data with comparison results.")
    if not os.path.isfile(os.path.join(output, "scrna_comparison_data.h5ad")):
        raise FileNotFoundError(f"scrna_comparison_data not found: {output}scrna_comparison_data.h5ad.")
    scrna = sc.read_h5ad(filename = os.path.join(output, "scrna_comparison_data.h5ad"))
    scrna.obs["scrublet_predicted_doublet"] = scrna.obs["scrublet_predicted_doublet_str"] == "True"
    del scrna.obs["scrublet_predicted_doublet_str"]
elif phases.get("comparison") == "execute":
    print("Comparing clusters for marker genes.", flush = True)
    for res in np.arange(0.1, 0.9, 0.1):
        key_added = f"leiden_{res:.1f}"
        print(key_added)
        deg_dir = os.path.join(output, f"degs_{key_added}")
        os.makedirs(deg_dir, exist_ok = True)
        writer = pd.ExcelWriter(os.path.join(deg_dir, f"{key_added}_degs.xlsx"), engine = "xlsxwriter")
        sc.tl.rank_genes_groups(scrna, groupby = key_added, reference = "rest", method = DEG_METHOD, max_iter = LOGREG_MAXITER, key_added = f"degs_{key_added}")
        for cluster in sorted(scrna.obs[key_added].unique()):
            print("Cluster: "+str(cluster))
            result = sc.get.rank_genes_groups_df(scrna, group = cluster, key = f"degs_{key_added}")
            result.to_excel(writer, sheet_name = f"Cluster_{cluster}", index = False)
            top_up = result[result["logfoldchanges"] > 0].nlargest(10, "logfoldchanges")
            top_down = result[result["logfoldchanges"] < 0].nsmallest(10, "logfoldchanges")
            top_genes = pd.concat([top_up, top_down[::-1]])
            plt.figure(figsize = (6, 6))
            bar_colors = ["red"] * len(top_up) + ["blue"] * len(top_down)
            sns.barplot(x = "logfoldchanges", y = "names", data = top_genes, palette = bar_colors)
            plt.axvline(0, color = "gray", linestyle = "--")
            plt.title(f"Top DEGs for Cluster {cluster} ({key_added})")
            plt.tight_layout()
            plt.savefig(os.path.join(deg_dir, f"degs_{key_added}_cluster_{cluster}.pdf"), bbox_inches = "tight")
            plt.savefig(os.path.join(deg_dir, f"degs_{key_added}_cluster_{cluster}.png"), bbox_inches = "tight", dpi = 300)
            plt.close()
        writer.close()
    
    
    print("Comparing stages.", flush = True)
    if(scrna.obs["stage"].nunique() > 1):
        print("We compare the stages.", flush = True)
        stage_dir = os.path.join(output, "degs_stage")
        os.makedirs(stage_dir, exist_ok = True)
        sc.tl.rank_genes_groups(scrna, groupby = "stage", method = DEG_METHOD, max_iter = LOGREG_MAXITER, key_added = "deg_genes_stage")
        stages = scrna.obs["stage"].unique().tolist()
        writer = pd.ExcelWriter(os.path.join(stage_dir, "stage_degs.xlsx"), engine = "xlsxwriter")
        for stage in stages:
            result = sc.get.rank_genes_groups_df(scrna, group = stage, key = "deg_genes_stage")
            result.to_excel(writer, sheet_name = f"Stage_{stage}", index = False)
            top_up = result[result["logfoldchanges"] > 0].nlargest(10, "logfoldchanges")
            top_down = result[result["logfoldchanges"] < 0].nsmallest(10, "logfoldchanges")
            top_genes = pd.concat([top_up, top_down[::-1]])
            plt.figure(figsize = (6, 6))
            bar_colors = ["red"] * len(top_up) + ["blue"] * len(top_down)
            sns.barplot(x = "logfoldchanges", y = "names", data = top_genes, palette = bar_colors)
            plt.axvline(0, color = "gray", linestyle = "--")
            plt.title(f"Top DEGs for Stage {stage}")
            plt.tight_layout()
            plt.savefig(os.path.join(stage_dir, f"degs_stage_{stage}.pdf"), bbox_inches = "tight")
            plt.savefig(os.path.join(stage_dir, f"degs_stage_{stage}.png"), bbox_inches = "tight", dpi = 300)
            plt.close()
        writer.close()
    
    print("Save the comparison object.", flush = True)
    scrna.obs["scrublet_predicted_doublet_str"] = scrna.obs["scrublet_predicted_doublet"].astype(str)
    del scrna.obs["scrublet_predicted_doublet"]
    scrna.write_h5ad(filename = os.path.join(output, "scrna_comparison_data.h5ad"))

