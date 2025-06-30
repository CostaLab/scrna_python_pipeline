# Phases: the different parts of the pipeline
# Which parts of the pipeline do you want to execute? load, execute or skip
# The comparison phase takes the longest time.
phases = {
    "raw"         : "execute", # Load the data.
    "filter"      : "execute", # Perform QC
    "integration" : "execute", # Integrate the data using harmony
    "cluster"     : "execute", # Perform the clustering using Leiden
    "comparison"  : "execute"  # Determine marker genes for each cluster and compare the stages.
}

# Where is your output supposed to go?
output = "/path/to/your/output/folder/"

# QC parameters
# Minimum number of cells per gene. int OR float
MINCELLS         = 5
# Minimum number of genes per cell. int OR float
MINGENES         = 50
# Maximum number of genes per cell. int OR float OR "Inf" (string)
MAXGENES         = "Inf"
# Maximum number of reads per cell. int OR float OR "Inf" (string)
MAXREADS         = 20000
# Maximum number of reads per cell. int OR float OR "Inf" (string)
MINREADS         = 500
# Maximum percentage of mitochondrial reads per cell. int OR float OR "Inf" (string)
PCT_MITO_CEILING = 15
# Minimum percentage of mitochondrial reads per cell. int OR float
PCT_MITO_FLOOR   = 0
# Maximum percentage of ribosome reads per cell. int OR float OR "Inf" (string)
PCT_RIBO_CEILING = 100
# Minimum percentage of ribosome reads per cell. int OR float
PCT_RIBO_FLOOR   = 0
# Do you want to remove doublets (True) or do you just want to determine them (False)?
doublet_switch   = False

# How many PCs do you want to use to determine the the variance contribution?
N_PCS_VARIANCE_CONTRIBUTION = 50
# How many neighbors do you want to use for PCA UMAP?
N_NEIGHBORS_PCA = 10
# How many PCs do you want to use for the PCA UMAP?
N_PCS_PCA = 40

# How many highly variable genes do you want to use for the Harmony integration?
N_TOP_GENES_INTEGRATION = 2000
# How many iterations should harmony try before terminating? Default is 10.
MAX_ITER_HARMONY = 10
# How many neighbors do you want to use for the UMAP?
N_NEIGHBORS_UMAP_HARMONY = 10
# How many PCs do you want to use to determine the neighbors for the UAMP?
N_PCS_UMAP_HARMONY = 40

# How verbose should the script be?
SETTING_VERBOSITY = 3 # verbosity: errors (0), warnings (1), info (2), hints (3)

# What is the test used for DEGs? Select from: "logreg", "t-test", "wilcoxon", "t-test_overestim_var"
DEG_METHOD = "wilcoxon"
# The maximum number of iterations for the logreg method.
LOGREG_MAXITER = 100

# Number of workers
WORKER_NUM = 1

### -------------- Data SRC-----------------------------
data_src = {
    "A_MxCre" : "data/A_MxCre/",
    "B_MxCre" : "data/B_MxCre/",
    "C_Csnk"  : "data/C_Csnk/",
    "D_Csnk"  : "data/D_Csnk/"
}


##------------------ SET REPLICATE GROUP --------------
stage_lst = {
    "A_MxCre" : "MxCre",
    "B_MxCre" : "MxCre",
    "C_Csnk"  : "Csnk",
    "D_Csnk"  : "Csnk"
}
