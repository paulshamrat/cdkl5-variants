#!/bin/bash
#SBATCH --job-name=HADDOCK_Docking
#SBATCH --output=logs/haddock_%j.out
#SBATCH --error=logs/haddock_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=40
#SBATCH --mem=128G
#SBATCH --time=4:00:00

# activate conda
module load anaconda3/2023.09-0
source activate

# Load any necessary modules (if required)
module load apptainer  # Only if needed; otherwise remove this line

# Create logs directory if it doesn't exist
mkdir -p logs

# Define paths
HADDOCKER_PATH="/project/ealexov/compbio/shamrat/250324_HADDOCKer/haddock3.sif"
CFG_FOLDER="05_cfg_files"

# Process each .cfg file sequentially
for cfg_file in "$CFG_FOLDER"/*.cfg; do
    echo "Running HADDOCK for $cfg_file"
    apptainer exec "$HADDOCKER_PATH" haddock3 "$cfg_file"

    # Check if HADDOCK failed
    if [ $? -ne 0 ]; then
        echo "Error in processing $cfg_file. Check logs." >> logs/errors.log
    fi
done

echo "All docking jobs completed."
