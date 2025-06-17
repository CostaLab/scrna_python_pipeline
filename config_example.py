output = "/path/to/your/output/folder/"

# filtering params 
MINCELLS  = 5
MINGENES  = 50
MAXGENES = 20000
PCT_MITO_CEILING = 15
PCT_MITO_FLOOR = 0
PCT_RIBO_CEILING = 100
PCT_RIBO_FLOOR = 0

doublet_switch = False

# How verbose should the script be?
SETTING_VERBOSITY = 3 # verbosity: errors (0), warnings (1), info (2), hints (3)

# Number of workers
WORKER_NUM = 1

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
