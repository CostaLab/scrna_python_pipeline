# We load the libraries.
from __future__ import annotations
from pathlib import Path
from typing import Dict, Literal # , Optional
from pydantic import BaseModel, field_validator, ValidationError
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
import math
# import multiprocessing as mp
from collections import defaultdict
import glob
# In python 3.11+ we could use the built-in tomllib.
# For python 3.10-, we have to use tomli.
import tomli


# Defining some additional parameters.
PhaseAction = Literal["execute", "load", "skip"]
sc.logging.print_header()
sc.settings.set_figure_params(dpi = 300, facecolor = "white")


# Defining custom functions and classes.
def _str_inf_to_math_inf(v):
    """Convert 'Inf'/'inf'/'+inf'/'-inf' to +/- math.inf; leave others unchanged."""
    if isinstance(v, str):
        s = v.strip().lower()
        if s in {"inf", "+inf"}:
            return math.inf
        if s == "-inf":
            return -math.inf
    return v


class Options(BaseModel):
    output: Path
    
    # QC
    min_cells: int | float
    min_genes: int | float
    max_genes: float | str
    max_reads: int | float | str
    min_reads: int | float
    pct_mito_ceiling: int | float | str
    pct_mito_floor: int | float
    pct_ribo_ceiling: int | float | str
    pct_ribo_floor: int | float
    doublet_switch: bool
    
    # PCA / neighbors / UMAP
    n_pcs_variance_contribution: int
    n_neighbors_pca: int
    n_pcs_pca: int
    
    # Harmony
    n_top_genes_integration: int
    max_iter_harmony: int
    n_neighbors_umap_harmony: int
    n_pcs_umap_harmony: int
    
    # Verbosity / DEG / workers
    setting_verbosity: int
    deg_method: Literal["logreg", "t-test", "wilcoxon", "t-test_overestim_var"]
    logreg_maxiter: int
    p_value_cutoff: float
    worker_num: int
    
    @field_validator("output", mode = "after")
    def _expand_output(cls, v: Path) -> Path:
        return v.expanduser()
    
    # Convert any "Inf" fields to floats
    @field_validator(
        "max_genes",
        "max_reads",
        "pct_mito_ceiling",
        "pct_ribo_ceiling",
        mode = "before",
    )
    def _convert_inf_strings(cls, v):
        return _str_inf_to_math_inf(v)


class Config(BaseModel):
    phases: Dict[str, PhaseAction]
    data_src: Dict[str, Path]
    stage_lst: Dict[str, str]
    options: Options
    
    # Enforce known phase names.
    @field_validator("phases")
    def _known_phases(cls, d: Dict[str, PhaseAction]) -> Dict[str, PhaseAction]:
        allowed = {"raw", "filter", "integration", "cluster", "comparison"}
        unknown = set(d) - allowed
        if unknown:
            raise ValueError(f"Unknown phase(s): {sorted(unknown)}; allowed: {sorted(allowed)}")
        
        
        return d
    
    
    # Expand user paths in data_src
    @field_validator("data_src", mode = "after")
    def _expand_paths(cls, d: Dict[str, Path]) -> Dict[str, Path]:
        return {k: Path(p).expanduser() for k, p in d.items()}


def load_config(path: str | Path) -> Config:
    with open(path, "rb") as f:
        data = tomli.load(f)
    try:
        cfg = Config.model_validate(data)
    except ValidationError as e:
        # Friendly message for non-programmers
        print("\nConfiguration error:\n", e, file=sys.stderr)
        raise
    
    
    return cfg


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
# exec(open(config_path).read())
cfg = load_config(config_path)
sc.settings.verbosity = cfg.options.setting_verbosity


print("Generating the output folder if needed.", flush = True)
os.makedirs(cfg.options.output, exist_ok = True)
qc_dir = os.path.join(cfg.options.output, "qc_plots")
os.makedirs(qc_dir, exist_ok = True)
clustering_dir = os.path.join(cfg.options.output, "clustering_plots")
os.makedirs(clustering_dir, exist_ok = True)


# Checking the parameters.
check_numeric(cfg.options.min_cells, "min_cells")
check_numeric(cfg.options.min_genes, "min_genes")
check_numeric(cfg.options.min_reads, "min_reads")
check_numeric(cfg.options.logreg_maxiter, "logreg_maxiter")
check_numeric(cfg.options.n_pcs_variance_contribution, "n_pcs_variance_contribution")
check_numeric(cfg.options.n_neighbors_pca, "n_neighbors_pca")
check_numeric(cfg.options.n_pcs_pca, "n_pcs_pca")
check_numeric(cfg.options.n_top_genes_integration, "n_top_genes_integration")
check_numeric(cfg.options.max_iter_harmony, "max_iter_harmony")
check_numeric(cfg.options.n_neighbors_umap_harmony, "n_neighbors_umap_harmony")
check_numeric(cfg.options.n_pcs_umap_harmony, "n_pcs_umap_harmony")
check_numeric(cfg.options.p_value_cutoff, "p_value_cutoff")

cfg.options.max_genes        = check_numeric_or_inf(cfg.options.max_genes, "max_genes")
cfg.options.max_reads        = check_numeric_or_inf(cfg.options.max_reads, "max_reads")
cfg.options.pct_mito_ceiling = check_numeric_or_inf(cfg.options.pct_mito_ceiling, "pct_mito_ceiling")
cfg.options.pct_mito_floor   = check_numeric_or_inf(cfg.options.pct_mito_floor, "pct_mito_floor")
cfg.options.pct_ribo_ceiling = check_numeric_or_inf(cfg.options.pct_ribo_ceiling, "pct_ribo_ceiling")
cfg.options.pct_ribo_floor   = check_numeric_or_inf(cfg.options.pct_ribo_floor, "pct_ribo_floor")


if(len(cfg.phases) == 0):
    raise ValueError("You are not supplying any phases.")


for key, val in cfg.phases.items():
    check_strings(val, allowed_values = ("execute", "load", "skip"), var_name = key)


# We check if a phase has a required phase.
required_phases_for_phases = {
    "filter": ["integration"],
    "integration": ["cluster"],
    "cluster": ["comparison"]
}

for phase, dependents in required_phases_for_phases.items():
    if cfg.phases.get(phase) == "skip":
        for dependent_use in dependents:
            if cfg.phases.get(dependent_use) == "execute":
                raise ValueError(f"Cannot execute '{dependent_use}' phase without loading or executing '{phase}' phase.")


print("We get the list of sample names.", flush = True)
sample_names = list(cfg.stage_lst.keys())




if cfg.phases.get("raw") == "skip":
    print("Skipping data loading.", flush = True)
elif cfg.phases.get("raw") == "load":
    print("Loading the raw data.", flush = True)
    if not os.path.isfile(os.path.join(cfg.options.output, "scrna_raw_data.h5ad")):
        raise FileNotFoundError(f"scrna_raw_data not found: {cfg.options.output}scrna_raw_data.h5ad.")
    
    
    scrna = sc.read_h5ad(os.path.join(cfg.options.output, "scrna_raw_data.h5ad"))#, backed = "r+") # Backed is less memory intensive, but does not work. The matrix cannot be read by the QC functions.
elif cfg.phases.get("raw") == "execute":
    print("Loading the data.", flush = True)
    scrna = sc.read_10x_mtx(cfg.data_src[sample_names[0]], cache = True)
    if(len(sample_names) > 1):
        scrna = [scrna]
        for x in sample_names[1:]:
            print(x, flush = True)
            scrna.append(sc.read_10x_mtx(cfg.data_src[x], cache = True))
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
    
    
    scrna.obs["stage"] = [cfg.stage_lst[x] for x in list(scrna.obs["batch"])]
    print("Saving the raw data.", flush = True)
    scrna.write_h5ad(filename = os.path.join(cfg.options.output, "scrna_raw_data.h5ad"))




if cfg.phases.get("filter") == "skip":
    print("Skipping data filtering.", flush = True)
elif cfg.phases.get("filter") == "load":
    print("Loading the filtered data.", flush = True)
    if not os.path.isfile(os.path.join(cfg.options.output, "scrna_filtered_data.h5ad")):
        raise FileNotFoundError(f"scrna_filtered_data not found: {cfg.options.output}scrna_filtered_data.h5ad.")
    
    
    scrna = sc.read_h5ad(os.path.join(cfg.options.output, "scrna_filtered_data.h5ad"))
    scrna.obs["scrublet_predicted_doublet"] = scrna.obs["scrublet_predicted_doublet_str"] == "True"
    del scrna.obs["scrublet_predicted_doublet_str"]
elif cfg.phases.get("filter") == "execute":
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
    sc.pp.filter_cells(scrna, min_genes = cfg.options.min_cells)
    sc.pp.filter_genes(scrna, min_cells = cfg.options.min_genes)
    scrna = scrna[scrna.obs.total_counts > cfg.options.min_reads, :]
    scrna = scrna[scrna.obs.pct_counts_mito < cfg.options.pct_mito_ceiling, :]
    scrna = scrna[scrna.obs.pct_counts_mito > cfg.options.pct_mito_floor, :]
    scrna = scrna[scrna.obs.pct_counts_ribo < cfg.options.pct_ribo_ceiling, :]
    scrna = scrna[scrna.obs.pct_counts_ribo > cfg.options.pct_ribo_floor, :]
    
    
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
        scrub.set_embedding("UMAP", scr.get_umap(scrub.manifold_obs_, n_neighbors = 15, min_dist = 0.3))
        x = scrub._embeddings["UMAP"][:,0]
        y = scrub._embeddings["UMAP"][:,1]
        predicted_doublets_index = np.argsort(preds)
        plt.figure(figsize=(7, 6))
        plt.scatter(x, y, c = preds[predicted_doublets_index], cmap = scr.custom_cmap([[.7,.7,.7], [1,0,0]]), s = 2)
        plt.xlabel("UMAP 1")
        plt.ylabel("UMAP 2")
        plt.title(batch)
        plt.tight_layout()
        pdf_path = os.path.join(qc_dir, f"{batch}_doublet_UMAP.pdf")
        plt.savefig(pdf_path, bbox_inches = "tight")
        png_path = os.path.join(qc_dir, f"{batch}_doublet_UMAP.png")
        plt.savefig(png_path, bbox_inches = "tight", dpi = 300)
        plt.close()
    
    
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
    if cfg.options.doublet_switch:
        print("Removing doublets.", flush = True)
        scrna = scrna[scrna.obs["scrublet_predicted_doublet"] == False, :].copy()
    
    
    print("Filtering genes and cells after doublet removal.", flush = True)
    sc.pp.filter_cells(scrna, max_genes = cfg.options.max_genes)
    scrna = scrna[scrna.obs.total_counts < cfg.options.max_reads, :]
    
    
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
    scrna.write_h5ad(filename = os.path.join(cfg.options.output, "scrna_filtered_data.h5ad"))
    scrna.obs["scrublet_predicted_doublet"] = scrna.obs["scrublet_predicted_doublet_str"] == "True"
    del scrna.obs["scrublet_predicted_doublet_str"]




if cfg.phases.get("integration") == "skip":
    print("Skipping data integration.", flush = True)
elif cfg.phases.get("integration") == "load":
    print("Loading the integrated data.", flush = True)
    if not os.path.isfile(os.path.join(cfg.options.output, "scrna_integrated_data.h5ad")):
        raise FileNotFoundError(f"scrna_integrated_data not found: {cfg.options.output}scrna_integrated_data.h5ad.")
    
    
    scrna = sc.read_h5ad(filename = os.path.join(cfg.options.output, "scrna_integrated_data.h5ad"))
    scrna.obs["scrublet_predicted_doublet"] = scrna.obs["scrublet_predicted_doublet_str"] == "True"
    del scrna.obs["scrublet_predicted_doublet_str"]
elif cfg.phases.get("integration") == "execute":
    print("Determining highly variable genes.", flush = True)
    # We only scale on the highly variable genes to reduce the amount of memory needed.
    sc.pp.highly_variable_genes(scrna, n_top_genes = cfg.options.n_top_genes_integration, flavor = "seurat_v3", subset = False, layer = "counts")
    # Generating a second object to reduce the memory needed in the integration.
    scrna_hvgs = scrna[:, scrna.var["highly_variable"]].copy()
    
    
    print("Regressing out confounders.", flush = True)
    sc.pp.regress_out(scrna_hvgs, keys = ["total_counts", "pct_counts_mito", "S_score", "G2M_score"], n_jobs = cfg.options.worker_num)
    sc.pp.scale(scrna_hvgs)
    
    
    print("Dimension Reduction.", flush = True)
    sc.tl.pca(scrna_hvgs, svd_solver = "arpack")
    
    
    print("Looking into the PC variance contribution.", flush = True)
    plt.figure(figsize = (5, 5))
    sc.pl.pca_variance_ratio(scrna_hvgs, log = True, n_pcs = cfg.options.n_pcs_variance_contribution, show = False)
    plt.savefig(os.path.join(qc_dir, "pca_variance_ratio.png"), dpi = 300, bbox_inches = "tight")
    plt.savefig(os.path.join(qc_dir, "pca_variance_ratio.pdf"), bbox_inches = "tight")
    plt.close()
    
    
    print("PCA.", flush = True)
    sc.pp.neighbors(scrna_hvgs, n_neighbors = cfg.options.n_neighbors_pca, n_pcs = cfg.options.n_pcs_pca, use_rep = "X_pca")
    sc.tl.umap(scrna_hvgs)
    scrna_hvgs.obsm["X_pca_umap"] = scrna_hvgs.obsm["X_umap"] 
    
    
    print("Sample Integration.", flush = True)
    # In case of only one sample (batch), harmony_integrate simply returns the original PCA.
    # This does not cause a crash. To keep everything simple, we will ignore this fact.
    sce.pp.harmony_integrate(scrna_hvgs, key = "batch", max_iter_harmony = cfg.options.max_iter_harmony)
    sc.pp.neighbors(scrna_hvgs, n_neighbors = cfg.options.n_neighbors_umap_harmony, n_pcs = cfg.options.n_pcs_umap_harmony, use_rep = "X_pca_harmony")
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
    scrna.write_h5ad(filename = os.path.join(cfg.options.output, "scrna_integrated_data.h5ad"))
    scrna.obs["scrublet_predicted_doublet"] = scrna.obs["scrublet_predicted_doublet_str"] == "True"
    del scrna.obs["scrublet_predicted_doublet_str"]




if cfg.phases.get("cluster") == "skip":
    print("Skipping data clustering.", flush = True)
elif cfg.phases.get("cluster") == "load":
    print("Loading the clustered data.", flush = True)
    if not os.path.isfile(os.path.join(cfg.options.output, "scrna_clustered_data.h5ad")):
        raise FileNotFoundError(f"scrna_clustered_data not found: {cfg.options.output}scrna_clustered_data.h5ad.")
    scrna = sc.read_h5ad(filename = os.path.join(cfg.options.output, "scrna_clustered_data.h5ad"))
    scrna.obs["scrublet_predicted_doublet"] = scrna.obs["scrublet_predicted_doublet_str"] == "True"
    del scrna.obs["scrublet_predicted_doublet_str"]
elif cfg.phases.get("cluster") == "execute":
    print("Data Clustering.", flush = True)
    for res in np.arange(0.1, 0.9, 0.1):
        key_added = f"leiden_{res:.1f}"
        print(key_added, flush = True)
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
    scrna.write_h5ad(filename = os.path.join(cfg.options.output, "scrna_clustered_data.h5ad"))
    scrna.obs["scrublet_predicted_doublet"] = scrna.obs["scrublet_predicted_doublet_str"] == "True"
    del scrna.obs["scrublet_predicted_doublet_str"]




if cfg.phases.get("comparison") == "skip":
    print("Skipping data comparison.", flush = True)
elif cfg.phases.get("comparison") == "load":
    print("Loading the data with comparison results.", flush = True)
    if not os.path.isfile(os.path.join(cfg.options.output, "scrna_comparison_data.h5ad")):
        raise FileNotFoundError(f"scrna_comparison_data not found: {cfg.options.output}scrna_comparison_data.h5ad.")
    scrna = sc.read_h5ad(filename = os.path.join(cfg.options.output, "scrna_comparison_data.h5ad"))
    scrna.obs["scrublet_predicted_doublet"] = scrna.obs["scrublet_predicted_doublet_str"] == "True"
    del scrna.obs["scrublet_predicted_doublet_str"]
elif cfg.phases.get("comparison") == "execute":
    print("Comparing clusters for marker genes.", flush = True)
    for res in np.arange(0.1, 0.9, 0.1):
        key_added = f"leiden_{res:.1f}"
        print(key_added, flush = True)
        deg_dir = os.path.join(cfg.options.output, f"degs_{key_added}")
        os.makedirs(deg_dir, exist_ok = True)
        writer = pd.ExcelWriter(os.path.join(deg_dir, f"{key_added}_degs.xlsx"), engine = "xlsxwriter")
        sc.tl.rank_genes_groups(scrna, groupby = key_added, reference = "rest", method = cfg.options.deg_method, max_iter = cfg.options.logreg_maxiter, key_added = f"degs_{key_added}")
        for cluster in sorted(scrna.obs[key_added].unique()):
            print(f"Cluster: {cluster}", flush = True)
            result = sc.get.rank_genes_groups_df(scrna, group = cluster, key = f"degs_{key_added}")
            result.to_excel(writer, sheet_name = f"Cluster_{cluster}", index = False)
            
            # Selecting the top 10 up regulated genes that are also statistically significant.
            sig_up = result[(result["logfoldchanges"] > 0) & (result["pvals_adj"] < cfg.options.p_value_cutoff)].nlargest(10, "logfoldchanges")
            # If there are less than 10 sig_up genes, we fill up with non significant genes.
            non_sig_up = (result[(result["logfoldchanges"] > 0) & (result["pvals_adj"] >= cfg.options.p_value_cutoff)].head(10 - len(sig_up))).nlargest(10, "logfoldchanges")
            up_selected = pd.concat([sig_up, non_sig_up], ignore_index = True)
            
            # Selecting the top 10 down regulated genes that are also statistically significant.
            sig_down = result[(result["logfoldchanges"] <= 0) & (result["pvals_adj"] < cfg.options.p_value_cutoff)].nsmallest(10, "logfoldchanges")
            sig_down = sig_down.sort_values("logfoldchanges", ascending = False)
            # If there are less than 10 sig_up genes, we fill up with non significant genes.
            non_sig_down = (result[(result["logfoldchanges"] <= 0) & (result["pvals_adj"] >= cfg.options.p_value_cutoff)].head(10 - len(sig_down))).nsmallest(10, "logfoldchanges")
            non_sig_down = non_sig_down.sort_values("logfoldchanges", ascending = True)
            down_selected = pd.concat([sig_down, non_sig_down], ignore_index = True)
            
            top_genes = pd.concat([up_selected, down_selected], ignore_index = True)
            plt.figure(figsize = (6, 6))
            bar_colors = ["red"] * len(up_selected) + ["blue"] * len(down_selected)
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
        stage_dir = os.path.join(cfg.options.output, "degs_stage")
        os.makedirs(stage_dir, exist_ok = True)
        sc.tl.rank_genes_groups(scrna, groupby = "stage", method = cfg.options.deg_method, max_iter = cfg.options.logreg_maxiter, key_added = "deg_genes_stage")
        stages = scrna.obs["stage"].unique().tolist()
        writer = pd.ExcelWriter(os.path.join(stage_dir, "stage_degs.xlsx"), engine = "xlsxwriter")
        for stage in stages:
            result = sc.get.rank_genes_groups_df(scrna, group = stage, key = "deg_genes_stage")
            result.to_excel(writer, sheet_name = f"Stage_{stage}", index = False)
            
            # Selecting the top 10 up regulated genes that are also statistically significant.
            sig_up = result[(result["logfoldchanges"] > 0) & (result["pvals_adj"] < cfg.options.p_value_cutoff)].nlargest(10, "logfoldchanges")
            # If there are less than 10 sig_up genes, we fill up with non significant genes.
            non_sig_up = (result[(result["logfoldchanges"] > 0) & (result["pvals_adj"] >= cfg.options.p_value_cutoff)].head(10 - len(sig_up))).nlargest(10, "logfoldchanges")
            up_selected = pd.concat([sig_up, non_sig_up], ignore_index = True)
            
            # Selecting the top 10 down regulated genes that are also statistically significant.
            sig_down = result[(result["logfoldchanges"] <= 0) & (result["pvals_adj"] < cfg.options.p_value_cutoff)].nsmallest(10, "logfoldchanges")
            sig_down = sig_down.sort_values("logfoldchanges", ascending = False)
            # If there are less than 10 sig_up genes, we fill up with non significant genes.
            non_sig_down = (result[(result["logfoldchanges"] <= 0) & (result["pvals_adj"] >= cfg.options.p_value_cutoff)].head(10 - len(sig_down))).nsmallest(10, "logfoldchanges")
            non_sig_down = non_sig_down.sort_values("logfoldchanges", ascending = True)
            down_selected = pd.concat([sig_down, non_sig_down], ignore_index = True)
            
            top_genes = pd.concat([up_selected, down_selected], ignore_index = True)
            plt.figure(figsize = (6, 6))
            bar_colors = ["red"] * len(up_selected) + ["blue"] * len(down_selected)
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
    scrna.write_h5ad(filename = os.path.join(cfg.options.output, "scrna_comparison_data.h5ad"))

