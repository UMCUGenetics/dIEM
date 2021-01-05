#file:config.R


# do not remove the quotes in this file, everything that is not a number
# in this file should be in double quotes
# every line with one or more #'s at the beginning is a comment and will not be 
# read by R.



# name of the run, will be the name of output folder and should not already exist 
run_name = "2020-008-22q11-zscore5"

# path: to DIMS excel file
path_DIMSfile = "/Volumes/metab/Metabolomics/Research Metabolic Diagnostics/Metabolomics Projects/Projects 2020/Project 2020_008 22q11/Bioinformatics/Project_2020_008_22q11.xlsx"




# binary variable: run function, yes(1) or no(0)
algorithm = 1
ratios = 1
violin = 1

# path: where should the outputfolder be made? Use your own output folder
path_output = "~/Documents/dIEM/output"


# standard input data files:


### if violin = 1 
# directory: folder in which all metabolite lists are (.txt)
path_txtfiles = "~/github/dIEM/stofgroups"

### if algorithm = 1
path_expected = "~/github/dIEM/data/Expected_version20191108_CSV.csv"

### if ratios = 1
path_ratios = "~/github/dIEM/data/Ratios_20191108_CSV.csv"








#### All settings below are not used yet.

#top_diseases = 5
#top_metab = 20
# integer: are the sample names headers on row 1 or row 2 in the DIMS excel? (default 1)
header_row = 1
# column name where the data starts (default B)
col_start = "B"
zscore_cutoff = 5
xaxis_cutoff = 20