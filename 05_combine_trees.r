source activate bcftools

R

options(scipen=999)

# list all the files in the trees directory
x_files <- list.files(pattern="*tre")

# find the chromosome, start, and end for each tree
x_names <- sapply(strsplit(sapply(strsplit(x_files, "bipartitions."), "[[", 2), ".tre"), "[[", 1)
x_chrom <- sapply(strsplit(x_names, "__"), "[[", 1)
x_start <- sapply(strsplit(x_names, "__"), "[[", 2)
x_end <- sapply(strsplit(x_names, "__"), "[[", 3)

# write tree info
write.table(cbind(x_chrom, x_start, x_end), file="_glauc_tree_info.txt", sep="\t", quote=F, row.names=F)

# trees into one file
tree_list <- list()
for(a in 1:length(x_files)) {
	tree_list[[a]] <- scan(x_files[a], what="character")
}
tree_list <- unlist(tree_list)
write(tree_list, file="_glauc_50kbp_windows.trees", ncolumns=1)


