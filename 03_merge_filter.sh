#!/bin/sh
#SBATCH --chdir=./
#SBATCH --job-name=merge
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --partition=nocona
#SBATCH --time=48:00:00
#SBATCH --mem-per-cpu=8G
#SBATCH --array=1-28

source activate bcftools

# define main working directory
workdir=/lustre/scratch/jmanthey/03_pygmy_owl

# define variables
region_array=$( head -n${SLURM_ARRAY_TASK_ID} ${workdir}/scaffolds.txt | tail -n1 )

# run bcftools to merge the vcf files
bcftools merge -m id --regions ${region_array} ${workdir}/03_vcf/*vcf.gz > \
${workdir}/04_vcf/${region_array}.vcf

# filter based on missing data for each of the subsets of data
vcftools --vcf ${workdir}/04_vcf/${region_array}.vcf --max-missing 0.5 \
--max-alleles 2 --max-maf 0.49 --recode --recode-INFO-all \
--out ${workdir}/05_phylo/${region_array}

bgzip ${workdir}/05_phylo/${region_array}.recode.vcf

tabix ${workdir}/05_phylo/${region_array}.recode.vcf.gz
