#!/bin/bash

#SBATCH --account=project_2003847
#SBATCH --partition=small
#SBATCH --time=7:00:00
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=30G

module load gsl
srun /projappl/project_2003847/nemo/bin/nemo2.4.0 name.ini
mv name.ini name/name.ini
mv run_sim_name.sh name/run_sim_name.sh
tar -zcvf name.tar.gz name
rm -r name