#file:config.R

# folder name in which all output will be written
run_name = "validation_monsternrs_1"
# binary variable: run function, yes(1) or no(0)
algorithm = 1
ratios = 1
violin = 1
# path: to DIMS excel file
#path_DIMSfile = "/Volumes/LAB/metab/Metabolomics/Research Metabolic Diagnostics/Metabolomics Projects/Projects 2018/Project 2018_012 Algoritme_TestSet_DBS/RES_DBS_20181001_Run1_Algoritmetest/Bioinformatics_master_ppm5/RES_DBS_20181001_Run1_Algoritmetest_master_ppm5_filter.xlsx"
path_DIMSfile = "~/Documents/dIEM/test_Helixcodes.xlsx"
# integer: are the sample names headers on row 1 or row 2 in the DIMS excel? (default 1)
header_row = 1
# column name where the data starts (default B)
col_start = "B"

### if violin = 1 
# directory: folder in which all metabolite lists are (.txt)
path_txtfiles = "~/Documents/dIEM/stofgroups"
zscore_cutoff = 5
xaxis_cutoff = 20


### if algorithm = 1
path_expected = "~/Documents/dIEM/Expected_version20191108_CSV.csv"

### if ratios = 1
path_ratios = "~/Documents/dIEM/Ratios_20191108_CSV.csv"


### if algoritm = 1 & violin = 1
#top_diseases = 5
#top_metab = 20
