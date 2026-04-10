#!/bin/bash

#SBATCH --account=project_2003847
#SBATCH --partition=small
#SBATCH --time=7:00:00
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=30G

module load gsl
srun /projappl/project_2003847/nemo/bin/nemo2.4.0 chr10_V0.5_c1_rep1.ini
mv chr10_V0.5_c1_rep1.ini chr10_V0.5_c1_rep1/chr10_V0.5_c1_rep1.ini
mv run_sim_chr10_V0.5_c1_rep1.sh chr10_V0.5_c1_rep1/run_sim_chr10_V0.5_c1_rep1.sh
tar -zcvf chr10_V0.5_c1_rep1.tar.gz chr10_V0.5_c1_rep1
rm -r chr10_V0.5_c1_rep1