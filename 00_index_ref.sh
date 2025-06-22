interactive -p nocona

cd references

source activate bcftools

samtools faidx GCF_003259725.1_athCun1_genomic.fna

java -jar picard.jar CreateSequenceDictionary \
R=/home/jmanthey/references/GCF_003259725.1_athCun1_genomic.fna \
O=/home/jmanthey/references/GCF_003259725.1_athCun1_genomic.dict
