#file:config.R


# do not remove the quotes in this file, everything that is not a number
# in this file should be in double quotes
# every line with one or more #'s at the beginning is a comment and will not be 
# read by R.



# name of the run, will be the name of output folder and should not already exist 
run_name = "test_slicemaxmin"

# path: to DIMS excel file
#path_DIMSfile = "Y:/Untargeted\ Metabolomics/Algoritme/dIEM/voorbeeld_PL_Coco_Nodding_RUN1.xlsx"

path_DIMSfile = "/Users/ihoek3/Documents/dIEM/voorbeeld_PL_Diagn17_RUN10.xlsx"


# path: where should the outputfolder be made? Use your own output folder
#path_output = "Y:/Untargeted\ Metabolomics/Algoritme/dIEM/test-output"
path_output = "/Users/ihoek3/Documents/dIEM/output"

# standard input data files:


### if violin = 1 
# directory: folder in which all metabolite lists are (.txt)
#path_txtfiles = "Y:/Untargeted\ Metabolomics/Algoritme/dIEM/dIEM-main/stofgroups"
path_txtfiles = "/Users/ihoek3/github/dIEM/stofgroups"
### if algorithm = 1
#path_expected = "Y:/Untargeted\ Metabolomics/Algoritme/dIEM/dIEM-main/data/Expected_version20191108_CSV.csv"
path_expected = "/Users/ihoek3/github/dIEM/data/Expected_version2020_CSV.csv"
### if ratios = 1
#path_ratios = "Y:/Untargeted\ Metabolomics/Algoritme/dIEM/dIEM-main/data/Ratios_20191108_CSV.csv"
path_ratios = "/Users/ihoek3/github/dIEM/data/Ratios_20210216_CSV.csv"









####### All settings below are not used yet.




# binary variable: run function, yes(1) or no(0)
algorithm = 1
ratios = 1
violin = 1
#top_diseases = 5
#top_metab = 20
# integer: are the sample names headers on row 1 or row 2 in the DIMS excel? (default 1)
header_row = 1
# column name where the data starts (default B)
col_start = "B"
zscore_cutoff = 5
xaxis_cutoff = 20