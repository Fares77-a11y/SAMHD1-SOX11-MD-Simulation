#!/bin/bash
#SBATCH --job-name=Center_SAMHD1xSOX11
#SBATCH --partition=tetralith
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --gpus-per-task=1
#SBATCH --time=6-00:00:00
#SBATCH --account=naiss2025-22-427  # Added your account

# Load the GPU version of GROMACS
module purge --force
module load GROMACS/2023.4-gpu-hpc1-g9

# Set the number of OpenMP threads for optimal GPU performance
export OMP_NUM_THREADS=16

# Define input and output file prefixes
init=step3_input
mini_prefix=step4.0_minimization
equi_prefix=step4.1_equilibration
equi_prefix_2=step4.2_equilibration
prod_prefix=step5_production

# for CONTINUATION
prev_step=step5_production
prod_step=step6_production

# Check if step4.2_equilibration.mdp exists, if not, create it
#if [ ! -f "${equi_prefix_2}.mdp" ]; then
    #echo "Creating ${equi_prefix_2}.mdp from ${equi_prefix}.mdp with NPT settings"
    #cp ${equi_prefix}.mdp ${equi_prefix_2}.mdp
    #sed -i 's/pcoupl\s*=\s*No/pcoupl = Parrinello-Rahman/g' ${equi_prefix_2}.mdp
    #grep -q "pcoupltype" ${equi_prefix_2}.mdp || echo "pcoupltype = isotropic" >> ${equi_prefix_2}.mdp
    #grep -q "tau_p" ${equi_prefix_2}.mdp || echo "tau_p = 2.0" >> ${equi_prefix_2}.mdp
    #grep -q "compressibility" ${equi_prefix_2}.mdp || echo "compressibility = 4.5e-5" >> ${equi_prefix_2}.mdp
    #grep -q "ref_p" ${equi_prefix_2}.mdp || echo "ref_p = 1.0" >> ${equi_prefix_2}.mdp
#fi

# minimization
#gmx grompp -f ${mini_prefix}.mdp -o ${mini_prefix}.tpr -c ${init}.gro -r ${init}.gro -p topol.top -n index.ndx -maxwarn 1
#gmx mdrun -ntmpi 1 -v -deffnm step4.0_minimization
#rm step4.0_minimization.trr

# Equilibration 4.1 NVT
#gmx grompp -f ${equi_prefix}.mdp -o ${equi_prefix}.tpr -c ${mini_prefix}.gro -r ${init}.gro -p topol.top -n index.ndx
#gmx mdrun -ntmpi 1 -v -deffnm step4.1_equilibration
#rm step4.1_equilibration.trr

# Equilibration 4.2 NPT
#gmx grompp -f ${equi_prefix_2}.mdp -o ${equi_prefix_2}.tpr -c ${equi_prefix}.gro -r ${init}.gro -t step4.1_equilibration.cpt -p topol.top -n index.ndx
#gmx mdrun -ntmpi 1 -v -deffnm step4.2_equilibration
#rm step4.2_equilibration.trr

# Production
#gmx grompp -f ${prod_prefix}.mdp -o ${prod_prefix}.tpr -c ${equi_prefix_2}.gro -r ${equi_prefix_2}.gro -t ${equi_prefix_2}.cpt -p topol.top -n index.ndx
#gmx mdrun -ntmpi 1 -v -deffnm ${prod_prefix}
#rm ${prod_prefix}.trr

# Continue run (commented out as in original)
gmx mdrun -ntmpi 1 -v -deffnm step5_production -noappend -cpi step5_production.cpt