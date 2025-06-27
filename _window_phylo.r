options(scipen=999)

# input args <- input file, popmap
args <- commandArgs(trailingOnly = TRUE)

# what do you want to calculate or use?
make_FASTA <- TRUE				# make a FASTA file (e.g., for RAxML)


# files in windows or whole chromosomes
vcf_type <- "windows" # two options, "windows" or "chromosomes"

# base file name for parsing chromosome and window bounds if applicable
# modify this section as needed to work with your file names and filtering strategy
# the chromosome_info object should end as a vector with 3 strings -> chromosome, start, end
# if any of this info is missing, feel free to just make those values 0 so they are still a number, etc.
# modify the below lines as needed 
filename_simple <- args[1]
filename_simple <- strsplit(filename_simple, "windows/")[[1]][2]
filename_simple <- strsplit(filename_simple, ".recode")[[1]][1]
if(vcf_type == "chromosomes") {
	chromosome_info <- c(strsplit(filename_simple, ".recode.vcf")[[1]][1], 0, 0)
} else if( vcf_type == "windows") {
	chromosome_info <- c(strsplit(filename_simple, "__")[[1]][1], 
                    as.numeric(strsplit(strsplit(filename_simple, "__")[[1]][2], "__")[[1]][1]),
                    as.numeric(strsplit(filename_simple, "__")[[1]][3]))
}

# determine base outname (edit as needed)
outname_base <- strsplit(filename_simple, ".recode.vcf")[[1]][1]


# define minimum number of sites to keep a fasta file 
min_sites <- 10000


############################################################################################################
############################################################################################################
############################################################################################################
############################################################################################################
# do not modify below this point
############################################################################################################
############################################################################################################
############################################################################################################
############################################################################################################



####################################
####################################
####################################
####################################
# functions used in below sections
####################################
####################################
####################################
####################################


####################################
# read in a small vcf (don't use 
# for large vcf files)
####################################
read_vcf <- function(input_file) {
  header <- readLines(input_file)
  header <- header[grep('^#C', header)]
  header <- strsplit(header, "\t")[[1]]
  vcf <- read.table(input_file, header=F)
  colnames(vcf) <- header
  return(vcf)
}


####################################
# make FASTA from VCF
# heterozygous sites choose one allele
# at random
####################################
create_fasta_from_vcf <- function(xxx, outname, num_sites_needed) {
	# keep going only if enough sites
	if(num_sites_needed <= nrow(xxx)) {
		
		# subset the genotypes from the allele info
		allele_info <- xxx[,c(4,5)]
		genotypes <- xxx[,10:ncol(xxx)]
		# define names of individuals in output fasta
		output_names_fasta <- paste(">", colnames(genotypes), sep="")
		
		# define heterozygous ambiguities
    		het <- rep("N", nrow(allele_info))
    		het[allele_info[,1] == "A" & allele_info[,2] == "C"] <- "M"
     	het[allele_info[,1] == "A" & allele_info[,2] == "G"] <- "R"
      	het[allele_info[,1] == "A" & allele_info[,2] == "T"] <- "W"
       	het[allele_info[,1] == "C" & allele_info[,2] == "A"] <- "M"
       	het[allele_info[,1] == "C" & allele_info[,2] == "G"] <- "S"
       	het[allele_info[,1] == "C" & allele_info[,2] == "T"] <- "Y"
       	het[allele_info[,1] == "G" & allele_info[,2] == "A"] <- "R"
       	het[allele_info[,1] == "G" & allele_info[,2] == "C"] <- "S"
       	het[allele_info[,1] == "G" & allele_info[,2] == "T"] <- "K"
       	het[allele_info[,1] == "T" & allele_info[,2] == "A"] <- "W"
       	het[allele_info[,1] == "T" & allele_info[,2] == "C"] <- "Y"
       	het[allele_info[,1] == "T" & allele_info[,2] == "G"] <- "K"
       	
       	# convert all numbers in genotypes to actual bases and ambiguities		
		for(a in 1:ncol(genotypes)) {
			# keep only genotype
			genotypes[,a] <- substr(genotypes[,a], 1, 3)
			# convert to bases
			genotypes[,a][genotypes[,a] == "0/0"] <- allele_info[genotypes[,a] == "0/0",1]
			genotypes[,a][genotypes[,a] == "1/1"] <- allele_info[genotypes[,a] == "1/1",2]
			genotypes[,a][genotypes[,a] == "./."] <- "?"
			genotypes[,a][genotypes[,a] == "0/1"] <- het[genotypes[,a] == "0/1"]
			# choose random allele for heterozygotes
			genotypes[,a][genotypes[,a] == "R"] <- sample(c("A", "G"), length(genotypes[,a][genotypes[,a] == "R"]), replace=T)
			genotypes[,a][genotypes[,a] == "K"] <- sample(c("T", "G"), length(genotypes[,a][genotypes[,a] == "K"]), replace=T)
			genotypes[,a][genotypes[,a] == "S"] <- sample(c("C", "G"), length(genotypes[,a][genotypes[,a] == "S"]), replace=T)
			genotypes[,a][genotypes[,a] == "Y"] <- sample(c("C", "T"), length(genotypes[,a][genotypes[,a] == "Y"]), replace=T)
			genotypes[,a][genotypes[,a] == "M"] <- sample(c("A", "C"), length(genotypes[,a][genotypes[,a] == "M"]), replace=T)
			genotypes[,a][genotypes[,a] == "W"] <- sample(c("A", "T"), length(genotypes[,a][genotypes[,a] == "W"]), replace=T)

		}
    
		# write output
		for(a in 1:ncol(genotypes)) {
			if(a == 1) {
				write(output_names_fasta[a], file=outname, ncolumns=1)
				write(paste(genotypes[,a], collapse=""), file=outname, ncolumns=1, append=T)
			} else {
				write(output_names_fasta[a], file=outname, ncolumns=1, append=T)
				write(paste(genotypes[,a], collapse=""), file=outname, ncolumns=1, append=T)
			}
		}
	}
}



####################################
# make FASTA file
####################################
if(make_FASTA == TRUE) {
	# define output fasta name
	output_fasta <- paste(strsplit(args[1], ".recode.vcf")[[1]][1], ".fasta", sep="")
	# read vcf
	input_file <- read_vcf(args[1])
	# create fasta sequence alignments for all files with a minimum number of sites
	create_fasta_from_vcf(input_file, output_fasta, min_sites)
}













