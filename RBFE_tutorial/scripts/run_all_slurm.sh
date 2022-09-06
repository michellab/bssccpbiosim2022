#! /bin/bash 
#
# Run all the lig prep, FEP prep, and the production runs
# Prior to this, the execution model should have been set up using the overall_network.ipynb
# This should have also tested if all the required parts are installed.

# TODO set the main directory filepath, make sure BSS branch and engines are sourced correctly
# export MAINDIRECTORY="RBFE_tutorial" # Set file path
# module load cuda/11.6
# module load amber/22
# module load gromacs/20.4
# somd="/export/users/XXX/sire.app/bin/somd-freenrg"
# export BSS="/export/users/XXX/anaconda3/bin/activate biosimspace-dev"

# export amber="/usr/local/amber22/amber.sh" # sourced in each script
# export gromacs="/usr/local/gromacs/bin/GMXRC" # sourced in each script
# source $BSS
# source $amber
# source $gromacs

# export all execution model files for later scripts
export lig_file="$MAINDIRECTORY/ligands.dat"
export net_file="$MAINDIRECTORY/network.dat"
export prot_file="$MAINDIRECTORY/protocol.dat"

# sourcing - as needed in the othe sh scripts
source extract_execution_model_bash.sh

###########
echo "The folder for all these runs is $MAINDIRECTORY"
echo ${lig_array[@]}
echo ${trans_array[@]}
echo ${eng_array[@]}
echo ${win_array[@]}

# make output dir for slurm out and err files
if [[ ! -d ../slurm_logs ]]; then
    mkdir ../slurm_logs
fi

# chmod all files so can be executed by sbatch.
chmod u+x run_ligprep_slurm.sh
chmod u+x run_fepprep_slurm.sh
chmod u+x run_production_slurm.sh

# Run the runs
# ligand prep
jidlig=$(sbatch --parsable --array=0-$((${#lig_array[@]}-1)) run_ligprep_slurm.sh)
echo "ligand prep jobid is $jidlig"

# FEP prep --dependency=afterany:${jidlig}
jidfep=$(sbatch  --parsable --array=0-$((${#trans_array[@]}-1)) run_fepprep_slurm.sh)
echo "FEP prep jobid is $jidfep"

# Production runs for the transformation --dependency=afterany:${jidfep}
for i in "${!trans_array[@]}"; do
jidprod=$(sbatch  --parsable --array=0-$((${win_array[i]}-1)) run_production_slurm.sh ${trans_array[i]} ${eng_array[i]} ${win_array[i]} $repeats)
echo "Production jobid for ${trans_array[i]}, ${eng_array[i]} is $jidprod"
done
