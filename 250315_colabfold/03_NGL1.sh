#!/bin/bash
#SBATCH --job-name=03_NGL1          # Job name
#SBATCH --nodes=1                         # Number of nodes
#SBATCH --ntasks=1                        # Number of tasks (processes)
#SBATCH --cpus-per-task=4                 # Number of CPU cores per task
#SBATCH --mem=16G                         # Total memory per node
#SBATCH --gres=gpu:v100:1                 # Request 1 A100 GPU
#SBATCH --time=24:00:00                   # Time limit (48 hours)
#SBATCH --output=colabfold_%j.out         # Standard output log
#SBATCH --error=colabfold_%j.err          # Standard error log

# Load necessary modules
module load anaconda3/2023.09-0
module load cuda/12.3.0

# Activate the ColabFold environment
source activate colabfold_env5

# Navigate to the working directory
cd /project/ealexov/compbio/shamrat

# make input file executable
chmod u+r cf_input2/03_NGL1.fasta

# Run ColabFold
colabfold_batch cf_input2/03_NGL1.fasta cf_output2/03_NGL1/
