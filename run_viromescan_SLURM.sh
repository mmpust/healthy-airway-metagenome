#!/bin/bash

# Last updated: 02 March 2022
# Author: Marie-Madlen Pust

# SLURM parameters
# set partition
#SBATCH -p normal

# set run on x MB node only
#SBATCH --mem 40000

# set run x cpus
#SBATCH --cpus-per-task 4

# set name of job
#SBATCH --job-name=viromescan

# Add miniconda3 to PATH
. /mariep/miniconda3/etc/profile.d/conda.sh

# Activate env on cluster node
conda activate viromescan >> /dev/null

fastq_in_1=$1
fastq_in_2=${fastq_in_1%_R1.fastq}_R2.fastq
out_dir=${fastq_in_1%_R1.fastq}_viromescan

srun /mnt/sfb900nfs/groups/tuemmler/mariep/miniconda3/envs/viromescan/viromescan/viromescan.sh \
 -m /mnt/sfb900nfs/groups/tuemmler/mariep/miniconda3/envs/viromescan/ \
 -d customised_phages_viruses \
 -1 $fastq_in_1 \
 -2 $fastq_in_2 \
 -o $out_dir
