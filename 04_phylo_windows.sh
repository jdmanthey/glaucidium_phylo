#!/bin/sh
#SBATCH --chdir=./
#SBATCH --job-name=stats
#SBATCH --nodes=1 --ntasks=2
#SBATCH --partition nocona
#SBATCH --time=48:00:00
#SBATCH --mem-per-cpu=8G
#SBATCH --array=1-491

source activate bcftools

# Set the number of runs that each SLURM task should do
PER_TASK=43

# Calculate the starting and ending values for this task based
# on the SLURM task and the number of runs per task.
START_NUM=$(( ($SLURM_ARRAY_TASK_ID - 1) * $PER_TASK + 1 ))
END_NUM=$(( $SLURM_ARRAY_TASK_ID * $PER_TASK ))

# Print the task and run range
echo This is task $SLURM_ARRAY_TASK_ID, which will do runs $START_NUM to $END_NUM

# Run the loop of runs for this task.
for (( run=$START_NUM; run<=$END_NUM; run++ )); do
	echo This is SLURM task $SLURM_ARRAY_TASK_ID, run number $run

	chrom_array=$( head -n${run} stat_helper.txt | cut -f1 | tail -n1 )

	start_array=$( head -n${run} stat_helper.txt | cut -f2 | tail -n1 )

	end_array=$( head -n${run} stat_helper.txt | cut -f3 | tail -n1 )

	gunzip -cd /lustre/scratch/jmanthey/03_pygmy_owl/05_phylo/QEEU01000028.1.recode.vcf.gz | grep "#" > /lustre/scratch/jmanthey/03_pygmy_owl/05_phylo/windows/${chrom_array}__${start_array}__${end_array}.recode.vcf

	tabix /lustre/scratch/jmanthey/03_pygmy_owl/05_phylo/${chrom_array}.recode.vcf.gz ${chrom_array}:${start_array}-${end_array} >> /lustre/scratch/jmanthey/03_pygmy_owl/05_phylo/windows/${chrom_array}__${start_array}__${end_array}.recode.vcf

	Rscript _window_phylo.r /lustre/scratch/jmanthey/03_pygmy_owl/05_phylo/windows/${chrom_array}__${start_array}__${end_array}.recode.vcf

	rm /lustre/scratch/jmanthey/03_pygmy_owl/05_phylo/windows/${chrom_array}__${start_array}__${end_array}.recode.vcf
	raxmlHPC-PTHREADS-SSE3 -T 2 -f a -x 50 -m GTRCAT -p 253 -N 10 -s /lustre/scratch/jmanthey/03_pygmy_owl/05_phylo/windows/${chrom_array}__${start_array}__${end_array}.fasta -n ${chrom_array}__${start_array}__${end_array}.tre -w /lustre/scratch/jmanthey/03_pygmy_owl/05_phylo/windows/

	rm /lustre/scratch/jmanthey/03_pygmy_owl/05_phylo/windows/${chrom_array}__${start_array}__${end_array}.fasta
	rm /lustre/scratch/jmanthey/03_pygmy_owl/05_phylo/windows/RAxML_bestTree.${chrom_array}__${start_array}__${end_array}.tre
	rm /lustre/scratch/jmanthey/03_pygmy_owl/05_phylo/windows/RAxML_bipartitionsBranchLabels.${chrom_array}__${start_array}__${end_array}.tre
	rm /lustre/scratch/jmanthey/03_pygmy_owl/05_phylo/windows/RAxML_info.${chrom_array}__${start_array}__${end_array}.tre

done
