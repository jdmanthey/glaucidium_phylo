# maximum clade credibility tree
sumtrees.py --output=_glauc_50kbp_consensus.tre --min-clade-freq=0.01 _glauc_50kbp_windows_renamed.trees

# ASTRAL 
java -jar ~/Astral/astral.5.6.3.jar -i_glauc_50kbp_windows_renamed.trees -o glauc_50kbp_astral.tre 2> glauc_50kbp_astral.log
