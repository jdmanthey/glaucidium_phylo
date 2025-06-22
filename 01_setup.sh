# define main working directory
workdir=/lustre/scratch/jmanthey/03_pygmy_owl

# make output directories
cd ${workdir}

mkdir 01_bam_files
mkdir 02_vcf
mkdir 03_vcf
mkdir 04_vcf
mkdir 05_phylo
mkdir 05_phylo/windows
mkdir 10_genotype
mkdir 11_merge_filter

# copy all bam files to 01_bam_files folder

cd  01_bam_files

while read -r name1 name2; do
	mv $name1 $name2
done < rename.txt
