options(scipen=999)

library(ape)

x <- read.tree("_glauc_50kbp_windows.trees")

rename <- read.table("rename_tips.txt", header=F, sep="\t")

# rename tips
for(a in 1:length(x)) {
	for(b in 1:nrow(rename)) {
		x[[a]]$tip.label[x[[a]]$tip.label == rename[b,1]] <- rename[b,2]
	}
}

write.tree(x, "_glauc_50kbp_windows_renamed.trees")





