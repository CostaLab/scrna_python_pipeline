#!/bin/bash
#SBATCH -c 1
#SBATCH --mem=400G
#SBATCH --output=/path/to/your/logs/folder/output.%J.%x.txt
#SBATCH --job-name=your_job_name
#SBATCH --mail-type=END
#SBATCH --mail-user=your_email@provider.topleveldomain
#SBATCH --time=10-10:00:00

echo "Running script:"
echo "run_example.sh"
echo "Starting time:"
date

# Loading the scRNA module.
module unload R
module unload scRNA
. /activating/your/miniconda3/etc/profile.d/conda.sh
conda activate single_cell_pipeline_python
unset PYTHONPATH

Rscript "data_factory_scanpy.py" \
	"/path/to/config.toml"
