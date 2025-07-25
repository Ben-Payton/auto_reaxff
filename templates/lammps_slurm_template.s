#!/bin/bash

#SBATCH --job-name={job_name}
#SBATCH --account=2101080209
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH --time=71:59:59
#SBATCH -o stdout
#SBATCH -e stderr
#SBATCH --export=ALL

# Go to the directoy from which our job was launched
cd $SLURM_SUBMIT_DIR

module purge
module load libs/mkl/2021.1
module load mpi/impi/2021.1
module load compilers/intel/2021.1
module load apps/lammps/intel-impi-plumed/23Jun2022
ulimit -s unlimited

echo "running job"


srun -n 4 lmp -in {infile}


echo "job has finished"
