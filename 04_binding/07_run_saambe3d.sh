#!/bin/bash
#SBATCH --job-name=SAAMBE3D_Batch
#SBATCH --output=logs/saambe3d_%j.out
#SBATCH --error=logs/saambe3d_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=64G
#SBATCH --time=08:00:00

# — PREPARE ENVIRONMENT —
module load anaconda3/2023.09-0
source activate /project/ealexov/compbio/software/py37_webservers

# Move to your binding directory
cd /project/ealexov/compbio/shamrat/250519_energy/04_binding

# Ensure output folders exist
mkdir -p logs outputs/saambe_3d

# Paths
SAAMBE_SCRIPT="/project/ealexov/compbio/shamrat/250415_binding/250416_saambe3d/SAAMBE-3D/saambe-3d.py"
MUTATION_LIST="/project/ealexov/compbio/shamrat/250519_energy/04_binding/02_saambe3d/mutations_list.txt"
MODEL_TYPE=1

# — LOOP OVER ALL PDBS —
for pdb_file in 06_cluster1_models/*.pdb; do
    echo "=== Processing $pdb_file ==="
    # Strip off the suffix to get e.g. P48436_SOX9_197-202
    base_name=$(basename "$pdb_file" "_cluster_1_model_1.pdb")

    # Run SAAMBE-3D with the mutation list
    python "$SAAMBE_SCRIPT" -i "$pdb_file" -f "$MUTATION_LIST" -d "$MODEL_TYPE"

    # Move the output into outputs/saambe_3d/
    mv output.out "outputs/saambe_3d/${base_name}.out"
    echo "Saved → outputs/saambe_3d/${base_name}.out"
done

echo "✅ All SAAMBE‑3D predictions completed."
