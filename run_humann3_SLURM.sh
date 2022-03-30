#!/bin/bash

# Last updated: 05 March 2022
# Author: Marie-Madlen Pust

# SLURM parameters
# set partition
#SBATCH -p normal

# set run on x MB node only
#SBATCH --mem 40000

# set run x cpus
#SBATCH --cpus-per-task 4

# set name of job
#SBATCH --job-name=HUMAnN3

# Add miniconda3 to PATH
. /mariep/miniconda3/etc/profile.d/conda.sh

# Activate env on cluster node
conda activate humann3_env >> /dev/null

input_fastq=$1
output_dir=${input_fastq%_kneaddata.fastq}

humann --input $input_fastq --input-format fastq --remove-stratified-output  --output $output_dir
