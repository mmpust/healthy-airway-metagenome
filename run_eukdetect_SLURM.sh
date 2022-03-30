#!/bin/bash

# Last updated: 02 Marc 2022
# Author: Marie-Madlen Pust

# SLURM parameters
# set partition
#SBATCH -p normal

# set run on x MB node only
#SBATCH --mem 40000

# set run x cpus
#SBATCH --cpus-per-task 4

# set name of job
#SBATCH --job-name=eukdetect

# Add miniconda3 to PATH
. /mariep/miniconda3/etc/profile.d/conda.sh

# Activate env on cluster node
conda activate eukdetect >> /dev/null

eukdetect --mode runall --configfile /mariep/programs/EukDetect/metagenome_sputum_configfile.yml
