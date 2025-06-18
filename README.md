# scrna_python_pipeline
Respository for the analysis pipeline using scanpy.

Usage:
run_example.sh
-- call data_factory.py 

config.py
-- the config file similar to the config file from the Seurat pipeline
-- this file should contain all parameters --> no parameters should be supplied directly to data_factory.py

data_factory_scanpy.py
-- the actual executing script
-- the script performs initial QC and generates some plots.
-- doublets are estimated using scrublet
-- samples are integrated using harmony 
-- clusters are determined for a resolution between 0.1 to 0.8 in steps of 0.1
-- DEGs are determined for each cluster for the different resolutions 
