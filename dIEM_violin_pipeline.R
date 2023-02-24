# For untargeted metabolomics, this tool calculates probability scores for 
# metabolic disorders. In addition, it provides visual support with violin plots 
# of the DIMS measurements for the lab specialists.
# Input needed: 
# 1. Excel file in which metabolites are listed with their intensities for
#    controls (with C in samplename) and patients (with P in samplename) and their
#    corresponding Z-scores. 
# 2. All files from github: https://github.com/UMCUGenetics/dIEM

# Feb 2023: code refactored to run as last step of DIMS pipeline on HPC
# setwd("/Users/mraves2/Metabolomics/DIMS_Test_run"); outdir <- getwd() # test env

library(dplyr)
library(reshape2)
library(data.table)
library(openxlsx) # for opening Excel file
library(ggplot2) # for plotting
library(gghighlight)
library(sys)
library(tidyr)
library(gridExtra)
library(grid)
library(testthat) # for unit testing

# define parameters - check after addition to run.sh
cmd_args <- commandArgs(trailingOnly = TRUE)
for (arg in cmd_args) cat("  ", arg, "\n", sep = "")

outdir <- cmd_args[1] # "/Users/nunen/Documents/Metab/test_set"
run_name <- cmd_args[2] # run_name <- "Test_run_5"
z_score <- as.numeric(cmd_args[3]) # calculate Z-scores or not? If no Z-scores, then no violin plots
dims_matrix <- cmd_args[4] # "DBS"

w <- 0.01 # seconds that system waits in between steps. Remove on HPC?

# The list of parameters can be shortened for HPC. Leave for now.
check.lists <- FALSE # check the lists of metabolites
rest <- FALSE
top <- 5 # number of diseases that score highest in algorithm to plot
threshold_IEM = 5 # probability score cut-off for plotting the top diseases
ratios_cutoff = -5 # z-score cutoff of axis on the left for top diseases

# Settings from config.R
# binary variable: run function, yes(1) or no(0). Can be removed at later stage
algorithm = 1
ratios = 1
violin = 1
# replace with: if (z_score == 1) { algorithm <- ratios <- violin <- 1 } ? Or omit?
# integer: are the sample names headers on row 1 or row 2 in the DIMS excel? (default 1)
header_row = 1
# column name where the data starts (default B)
col_start = "B"
zscore_cutoff = 5
xaxis_cutoff = 20

# source functions; check command later.
# source("/Users/mraves2/Development/dIEM/dIEM/functions/*.R")
source("/hpc/dbg_mz/production/DIMS/pipeline/scripts/AddOnFunctions/other_isobaric.R")
source("/hpc/dbg_mz/production/DIMS/pipeline/scripts/AddOnFunctions/make_plots.R")


#### STEP 0  load config settings ####  
# sink(file="dIEM_log.txt") # create log file in run.sh

# path to DIMS excel file
# path_DIMSfile = "/Users/ihoek3/Documents/dIEM/voorbeeld_PL_Diagn17_RUN10.xlsx"
path_DIMSfile = paste0(outdir, "/", run_name, ".xlsx") # ${outdir} in run.sh

# path: output folder 
# path_output = "/Users/ihoek3/Documents/dIEM/output"
path_output = output_dir = paste0(outdir, "/dIEM_refactored") # change to dIEM later
dir.create(path_output, showWarnings = F) # remove either path_output or output_dir later

# directory: folder in which all metabolite lists are (.txt)
# path_txtfiles = "/Volumes/LAB/metab/Untargeted Metabolomics/Algoritme/dIEM/dIEM-main/stofgroups"
# path_txtfiles = "/Users/mraves2/Development/dIEM/dIEM/metabolite_groups"
path_txtfiles = "/hpc/dbg_mz/tools/db/metabolite_groups"
### for algorithm step 4
# path_expected = "/Volumes/LAB/metab/Untargeted Metabolomics/Algoritme/dIEM/dIEM-main/data/Expected_version20191108_CSV.csv"
# path_expected = "/Users/mraves2/Metabolomics/dIEM-main/data/Expected_version20191108_CSV.csv"
path_expected = "/hpc/dbg_mz/tools/db/dIEM/Expected_version2020_CSV.csv"
### for ratios step 3
# path_ratios = "/Volumes/LAB/metab/Untargeted Metabolomics/Algoritme/dIEM/dIEM-main/data/Ratios_20191108_CSV.csv"
# path_ratios = "/Users/mraves2/Metabolomics/dIEM-main/data/Ratios_20230220_CSV.csv"
path_ratios = "/hpc/dbg_mz/tools/db/ratios/Ratios_20210216_CSV.csv"

# where do print and cat statements go on HPC?
if (exists("run_name")) {
  cat("\nThe config file is succesfully loaded from working directory.")
} else {
  cat("\n**** Error: Could not find a config file. please check if Working Directory is set in 'Session'. \n")
}


#### STEP 1: Preparation #### 
# in: run_name, path_DIMSfile, header_row ||| out: output_dir, DIMS

# Load the excel file.
dimsxls <- readWorkbook(xlsxFile = path_DIMSfile, sheet = 1, startRow = header_row)
if (exists("dimsxls")) {
  cat(paste0("\nThe excel file is succesfully loaded:\n -> ",path_DIMSfile))
} else {
  cat(paste0("\n\n**** Error: Could not find an excel file. Please check if path to excel file is correct in config.R:\n -> ",path_DIMSfile,"\n"))
}

Sys.sleep(w)

#### STEP 2: Edit DIMS data #####      
# in: DIMSxls ||| out: Data, nrcontr, nrpat
# Input: the xlsx file that comes out of the pipeline with format:
# [plots] [C] [P] [summary columns] [C_Zscore] [P_Zscore]
# Output: "_CSV.csv" file that is suited for the algorithm in shiny.

dims2 <- dimsxls
# Calculate the number of Cs and Ps in column names to extract the following numbers:
nrcontr <- length(grep("C",names(dims2)))/2   # Number of control samples
nrpat <- length(grep("P",names(dims2)))/2     # Number of patient samples
if (nrcontr + nrpat != length(grep("_Zscore", names(dims2)))) {
  cat("\n**** Error: there aren't as many intensities listed as Zscores")
}
cat(paste0("\n\n------------\n",nrcontr, " controls \n",nrpat," patients\n------------\n\n"))

# Get the columns HMDB_code and HMDB_name to the beginning. 
dims2 <- select(dims2, c(HMDB_code, HMDB_name), everything())
# Remove the columns from 'name' to 'pathway'
if (!is.na(dims2[1,3])){ # in case the excel had no empty "plots" column
  dims2 <- subset( dims2, select = -c(name : pathway ))
} else {
  dims2 <- subset( dims2, select = -c( name : pathway ))[-3] 
}
# Rename the columns 
names(dims2) <- gsub("avg.ctrls", "Mean_controls", names(dims2))
names(dims2) <- gsub("sd.ctrls",  "SD_controls", names(dims2))
names(dims2) <- gsub("HMDB_code", "HMDB.code", names(dims2))
names(dims2) <- gsub("HMDB_name", "HMDB.name", names(dims2))

#first, select the intensity columns by all cols minus nrsamples
nrsamples = nrcontr + nrpat
# start of the section with Z-scores
beginZscores = ncol(dims2) - nrsamples
# intensity columns and mean and standard deviation of controls
i <- c(3:(nrsamples+2))
#change the intensities to numeric values
dims2[, i] <- sapply(dims2[, i], as.numeric)

if (exists("dims2") & (length(dims2)<length(dimsxls))) {
  cat("\n### Step 2 # Edit dims data is done.\n")
} else {
  cat("\n**** Error: Could not execute step 2 \n")
}
#cat("### Step 2 # Edit dims data is done.\n")

Sys.sleep(w)

#### STEP 3: Calculate ratios ####      
# in: ratios, path_ratios, dims2, nrcontr, nrpat ||| out: Zscore (+file)
# This script loads the file with Ratios (path_ratios) and calculates 
# the ratios of the intensities of the given metabolites. It also calculates
# Zs-cores based on the avg and sd of the ratios of the controls.

# Input: dataframe with intenstities and Zscores of controls and patients:
# [HMDB.code] [HMDB.name] [C] [P] [Mean_controls] [SD_controls] [C_Zscore] [P_Zscore]

# Output: "_CSV.csv" file that is suited for the algorithm, with format:
# "_Ratios_CSV.csv" file, same file as above, but with ratio rows added.

dims3 <- dims2
if (ratios == 1) { # should be default
  cat(paste0("\nloading ratios file:\n ->  ", path_ratios, "\n"))
  RatioInput <- read.csv(path_ratios, sep=';', stringsAsFactors=FALSE)
  
  # Prepare empty data frame to fill with ratios
  Ratios <- setNames(data.frame(matrix(ncol=ncol(dims3), 
                                       nrow=nrow(RatioInput))), colnames(dims3))

  # put HMDB info into first two columns
  Ratios[ ,1:2] <- RatioInput[ ,1:2]

  # use ratios with or without log2? For now, keep both.
  Ratios_log <- Ratios
  # look for intensity columns (exclude Zscore columns)
  control_cols <- grep("C", colnames(Ratios)[1:which(colnames(Ratios) == "Mean_controls")])
  patient_cols <- grep("P", colnames(Ratios)[1:which(colnames(Ratios) == "Mean_controls")])
  intensity_cols <- c(control_cols, patient_cols)
  # calculate each of the ratios
  for (ratio_index in 1:nrow(RatioInput)) {
    ratio_numerator   <- RatioInput[ratio_index, "HMDB_numerator"] 
    ratio_numerator   <- strsplit(ratio_numerator, "plus")[[1]]
    ratio_denominator <- RatioInput[ratio_index, "HMDB_denominator"] 
    ratio_denominator <- strsplit(ratio_denominator, "plus")[[1]]
    # find these HMDB IDs in dataset. Could be a sum of multiple metabolites
    sel_denominator <- sel_numerator <- c()
    for (len in 1:length(ratio_numerator)) { 
      sel_numerator <- c(sel_numerator, which(dims3[ , "HMDB.code"] == ratio_numerator[len])) 
    }
    for (len in 1:length(ratio_denominator)) { 
      sel_denominator <- c(sel_denominator, which(dims3[ , "HMDB.code"] == ratio_denominator[len])) 
    }
    # calculate ratio
    Ratios[ratio_index, intensity_cols] <- apply(dims3[sel_numerator, intensity_cols], 2, sum) /
      apply(dims3[sel_denominator, intensity_cols], 2, sum)
    # calculate log of ratio
    Ratios_log[ratio_index, intensity_cols]<- log2(Ratios[ratio_index, intensity_cols])
  }
  
  # Calculate means and SD's of the calculated ratios, add them in 2 columns in ratio df.
  Ratios[ , "Mean_controls"] <- apply(Ratios[ , control_cols], 1, mean)
  Ratios[ , "SD_controls"]   <- apply(Ratios[ , control_cols], 1, sd)
  Ratios_log[ , "Mean_controls"] <- apply(Ratios_log[ , control_cols], 1, mean)
  Ratios_log[ , "SD_controls"]   <- apply(Ratios_log[ , control_cols], 1, sd)
  
  # Calc z-scores with the means and SD's 
  zscore_cols <- grep("Zscore", colnames(Ratios))
  for (zscore_col in zscore_cols) {
    # matching intensity column; add the two columns for mean and sd when calculating Z-scores
    int_col <- zscore_col - (nrsamples + 2) 
    # add unit test 
    test_that("column names are the same", {
      expect_equal(paste0(colnames(Ratios)[int_col], "_Zscore"), colnames(Ratios)[zscore_col])
    })
    # calculate Z-scores
    Ratios[ , zscore_col] <- (Ratios[ , int_col] - Ratios[ , "Mean_controls"]) / Ratios[ , "SD_controls"]
    Ratios_log[ , zscore_col] <- (Ratios_log[ , int_col] - Ratios_log[ , "Mean_controls"]) / Ratios_log[ , "SD_controls"]
  }
  
  # Add rows of the ratio hmdb codes to the data of zscores from the pipeline.
  dims4 <- rbind(Ratios, dims3)
  dims4_log <- rbind(Ratios_log, dims3)
  
  # Select only the cols with zscores of only the patients (Zscore)
  Zscore <- dims4[, c(1, 2, zscore_cols[grep("P", colnames(dims4)[zscore_cols])])]
  Zscore_log <- dims4_log[, c(1, 2, zscore_cols[grep("P", colnames(dims4)[zscore_cols])])]
  # And with the controls
  Zscore_all <- dims4[ , c(1, 2, zscore_cols)]
  Zscore_all_log <- dims4_log[ , c(1, 2, zscore_cols)]

}

Sys.sleep(w)

#### STEP 4: Run the algorithm #########      
# in: algoritm, path_expected, Zscore ||| out: ProbScore0 (+file)

# Zscore <- ObsZscore
if (algorithm == 1) {
  # Load data
  cat(paste0("\nloading expected file:\n ->  ", path_expected, "\n"))
  Expected <- read.csv(path_expected,sep=';',stringsAsFactors=FALSE)
  
  # prepare dataframe scaffold Rank
  Rank <- Zscore
  r <- nrow(Rank)
  # Fill df Rank with the ranks for each patient
  for (patient in 3:ncol(Zscore)) {
    # number of positive zscores in patient
    pos <- length(Zscore[patient][(Zscore[patient] > 0) == TRUE])
    # sort the column on zscore
    Rank <- Rank[order(-Rank[patient]), ]
    # Rank all positive zscores highest to lowest
    Rank[1:pos, patient] <- as.numeric(ordered(-Rank[1:pos, patient]))
    # Rank all negative zscores lowest to highest
    Rank[(pos+1):r, patient] <- as.numeric(ordered(Rank[(pos+1):r, patient]))
  }
  
  # Calculate metabolite score, using the dataframes with only values, and later add the cols without values (1&2).
  # new expected: -c(29,30) delete Blood column and hmdb name column
  Exp_Zscores <- merge(x=Expected, y=Zscore, by.x = c("HMDB.code"), by.y = c("HMDB.code"))[-c(29,30)]
  Exp_Zscores0 <- Exp_Zscores
  
  cat("setting some zscores to zero.\t")
  Exp_Zscores0[which(Exp_Zscores0$Change=="Increase" & Exp_Zscores0$Dispensability=="Indispensable"),29:ncol(Exp_Zscores0)] <- lapply(Exp_Zscores0[which(Exp_Zscores0$Change=="Increase" & Exp_Zscores0$Dispensability=="Indispensable"),29:ncol(Exp_Zscores0)], function(x) ifelse(x<=1.6 , 0, x))
  Exp_Zscores0[which(Exp_Zscores0$Change=="Decrease" & Exp_Zscores0$Dispensability=="Indispensable"),29:ncol(Exp_Zscores0)] <- lapply(Exp_Zscores0[which(Exp_Zscores0$Change=="Decrease" & Exp_Zscores0$Dispensability=="Indispensable"),29:ncol(Exp_Zscores0)], function(x) ifelse(x>=-1.2 , 0, x))
  cat("done.\n")
  # new expected: -c(29,30) delete Blood column and hmdb name column
  Exp_Rank <- merge(x=Expected, y=Rank, by.x = c("HMDB.code"), by.y = c("HMDB.code"))[-c(29,30)]
  cat("calculate rank score.\t")
  Exp_Metabscore <- cbind(Exp_Rank[order(Exp_Zscores0$HMDB.code),][,1:28],((Exp_Zscores0[order(Exp_Zscores0$HMDB.code),][-c(1:28)])/((Exp_Rank[order(Exp_Rank$HMDB.code),][-c(1:28)])*0.9)))
  
  cat("multiplying weight score and rank score.\t")
  Wscore <- Exp_Zscores0
  Wscore[29:ncol(Exp_Metabscore)] <- Exp_Metabscore$Total_Weight*Exp_Metabscore[29:ncol(Exp_Metabscore)]
  
  #De expected weight scores moeten (kolom 27) op sort staan. De weights horen van 9 tot -9 te gaan per m/z waarde.
  Wscore <- Wscore[order(Wscore$Disease,Wscore$Absolute_Weight,decreasing = TRUE),]
  cat("done.\n")
  
  dup = Wscore[,c('Disease','M.z')] # select columns to check duplicates
  uni <- Wscore[!duplicated(dup) | !duplicated(dup, fromLast=FALSE),]
  
  ProbScore <-  aggregate(uni[29:ncol(uni)],uni["Disease"],sum)
  ProbScore0 <- ProbScore
  
  for (pt in 2:ncol(ProbScore0)) {
    # for each patient (each column in df)
    # list of all diseases that have at least one metabolite zscore at 0
    dis_setto0 <- unique(Exp_Zscores0[which(Exp_Zscores0[pt+27]==0),][["Disease"]])
    # set the prob score of these diseases to 0 
    ProbScore0[which(ProbScore0$Disease %in% dis_setto0),pt]<- 0
  }
  
  disRank <- ProbScore0
  disRank[2:ncol(disRank)] <- lapply(2:ncol(disRank), function(x) as.numeric(ordered(-disRank[1:nrow(disRank),x])))
  # col names aanpassen van _Zscore naar _ProbScore, omdat het geen Zscore meer is.
  names(ProbScore0) <- gsub("_Zscore","_ProbScore",names(ProbScore0))
  
  # Create conditional formatting for output excel sheet. Colors according to values.
  wb <- createWorkbook()
  addWorksheet(wb, "Probability Scores")
  writeData(wb, "Probability Scores", ProbScore0)
  conditionalFormatting(wb, "Probability Scores", cols = 2:ncol(ProbScore0), rows = 1:nrow(ProbScore0), type = "colourScale", style = c("white","#FFFDA2","red"), rule = c(1, 10, 100)) # middle color is yellow
  saveWorkbook(wb, file = paste0(output_dir,"/algoritme_output_",run_name,".xlsx"), overwrite = TRUE)
  # check whether ProbScore df exists and is in expected dimensions.
  if (exists("Expected") & (length(disRank)==length(ProbScore0))) {
    cat("\n### Step 4 # Run the algorithm is done.\n\n")
  } else {
    cat("\n**** Error: Could not run algorithm. Check if path to Expected csv-file is correct in config.R. \n")
  }
  
  rm(wb)
}
Sys.sleep(w)


#### STEP 5: Make violin plots #####      
# in: algorithm / Zscore, violin, nrcontr, nrpat, Data, path_textfiles, zscore_cutoff, xaxis_cutoff, top_diseases, top_metab, output_dir ||| out: pdf file

# object Data is nowhere to be found.
# if ( algorithm == 0 ) {
#   Zscore <- Data[,c(1:2,(2*nrcontr+nrpat+5):(2*(nrcontr+nrpat)+4))]                                                         
# }

if (violin == 1) { # make violin plots
  isobarics_txt <- c()
  # Edit the DIMS output Zscores of all patients in format:
  # HMDB_code patientname1  patientname2
  names(Zscore) <- gsub("HMDB.code","HMDB_code", names(Zscore))
  names(Zscore) <- gsub("HMDB.name", "HMDB_name", names(Zscore))
  summed_nonames <- Zscore[-2] # remove the HMDB.name column
  summed <- Zscore
  if (check.lists) {
    summed.c <- Zscore
  }
  # remove the string "_Zscore" from columnnames
  names(summed) <- gsub("_Zscore", "", names(summed)) # remove the _Zscore from columnnames
  
  # Make a patient list so it can be looped over later in lapply.
  patient_list <- names(summed)[-2]
  patient_list[1] <- "all"
  # add list for IEM plots 
  patient_list2 <- gsub("all_IEM", "all_volle_x-as", paste0(patient_list, "_IEM"))
  
  # keep these two patient lists separate
  # patient_list <- c(patient_list, patient_list2) # otherwise the result is too complicated:
  # [1]  "all"               "P1002.1"           "P1003.1"           "P1005.2"           "P2021M04435.1"    
  # [6]  "P2021M04781.1"     "P2021M05527.1"     "P2021M05533.1"     "P2021M05536.1"     "P2021M05538.1"    
  # [11] "P2021M05543.1"     "P2021M05546.1"     "P2021M05552.1"     "P2021M05570.1"     "P2021M05622.1"    
  # [16] "P2021M05629.1"     "P2021M05633.1"     "P2021M05637.1"     "P2021M05640.1"     "P2021M05644.1"    
  # [21] "P2021M05647.1"     "P2021M05658.1"     "P2021M05663.1"     "alle_volle_x-as"   "P1002.1_IEM"      
  # [26] "P1003.1_IEM"       "P1005.2_IEM"       "P2021M04435.1_IEM" "P2021M04781.1_IEM" "P2021M05527.1_IEM"
  # [31] "P2021M05533.1_IEM" "P2021M05536.1_IEM" "P2021M05538.1_IEM" "P2021M05543.1_IEM" "P2021M05546.1_IEM"
  # [36] "P2021M05552.1_IEM" "P2021M05570.1_IEM" "P2021M05622.1_IEM" "P2021M05629.1_IEM" "P2021M05633.1_IEM"
  # [41] "P2021M05637.1_IEM" "P2021M05640.1_IEM" "P2021M05644.1_IEM" "P2021M05647.1_IEM" "P2021M05658.1_IEM"
  # [46] "P2021M05663.1_IEM"
  
  # from table Expected, choose columns "Disease", "HMDB.code", "Metabolite" and "Isobaric_Compounds_Names"
  Expected_red <- Expected[ ,c(5,14,13,18)]
  Expected_red <- Expected_red[!duplicated(Expected_red[,c(1,2)]),]
  names(Expected_red) <- gsub("HMDB.code", "HMDB_code",  names(Expected_red))  
  names(Expected_red) <- gsub("Metabolite", "HMDB_name", names(Expected_red))
  
  
  # Find all text files in the given folder, which contain metabolite lists of which
  # each file will be a page in the pdf with violin plots.
  # Make a PDF file for each of the categories in metabolite_dirs
  metabolite_dirs <- list.files(path=path_txtfiles, full.names=FALSE, recursive=FALSE)
  for (metabolite_dir in metabolite_dirs) {
    # create a directory for the output PDFs
    pdf_dir <- paste(output_dir, metabolite_dir, sep="/")
    dir.create(pdf_dir)
    metabolite_files <- list.files(path=paste(path_txtfiles, metabolite_dir, sep="/"), pattern="*.txt", full.names=FALSE, recursive=FALSE)

    # put all metabolites into one list
    metab_list_all <- list()
    index = 0
    cat("making plots from the input files:")
    # open the text files and add each to a list of dataframes (metab.list0)
    for (infile in metabolite_files) {
      index = index + 1
      metab.list1 <- read.table(paste(path_txtfiles, metabolite_dir, infile, sep="/"), sep = "\t", header = TRUE, quote="")
      isos <- other_isobaric(dims4, metab.list1, infile, check.lists) # puts qqqzzz before every HMDB_name
      metab.list2 <- rbind(metab.list1, isos)
      metab_list_all[[index]] <- metab.list2
      cat(paste0("\n", infile))
    } 
    # cat("\n")
      
    # if (rest) { # MPR: not sure what this does; check later.
    #   # reduce the columns from the expected df for further use and remove double 
    #   # metabolites per disease.
    #   Expected_rest <- Expected_red[!duplicated(Expected_red[,2]),][,c(2,3)]
    #   # Remove the metab-list from the expected df.
    # 
    #   cc = 0
    #   for (metab.list in metab.list0){
    #     if (ncol(metab.list) > 2) {
    #       metab.list4 <- metab.list
    #       metab.list <- metab.list[c(1,2)]
    #     }
    #     Expected_rest <- anti_join(Expected_rest, metab.list, by='HMDB_code')
    #     # when checking metabolite lists is needed 
    #     stoftest <- gsub(".txt","",stofgroup_files[cc]) # (string) stoftest name
    #     if (check.lists) {
    #       cc = cc + 1
    #       joined.c <- inner_join(metab.list, summed.c, by = "HMDB_code")
    #       check <- joined.c[,c(1:4)]
    #       write.xlsx(check, paste0(output_dir,"/", stoftest, "_check.xlsx"))
    #     }
    #   }
    #   
    #   index = index + 1 
    #   metab.list0[[index]] <- Expected_rest
    #   stofgroup_files[index] <- "rest_automatic"
    # }
    
    # voeg additionele lijsten toe, rest, de top 20 hoogst scorende en top 10 laagst.
    # leave out for now. List summed is not sorted, so why take the top 20?
    # metab.list0[[index+1]] <- summed[c(1:20),c(1,2)]
    # stofgroup_files[index+1] <- "top20_hoogst"
    # metab.list0[[index+2]] <- summed[c(1:10),c(1,2)]
    # stofgroup_files[index+2] <- "top10_laagst"
    
    # set parameters for plots
    n <- 0
    i_tot <- 24
    plot_height <- 0.40 * i_tot
    plot_width <- 6
    
    # build list of patients, replaced with clearer construction
    # First loop: no IEM plots
    for (pt_nr in 1:length(patient_list)) {
      pt <- patient_list[pt_nr]
      # normal violin plots:
      ThisProbScore = 0 # means that this is no IEM plot
      cat(paste0("\n","For ", pt, ", done with plot nr. "))
      
      if (startsWith(pt,"all")){
        # overview plots
        pdf(paste0(pdf_dir,"/", pt, "_patients_overview.pdf"), onefile = TRUE,
            width = plot_width, height = plot_height) # create the PDF device
      } else {
        # patient plots
        pdf(paste0(pdf_dir, "/", pt, ".pdf"), onefile = TRUE,
            width = plot_width, height = plot_height) # create the PDF device
      }
      
      for (metab_list_nr in 1:length(metab_list_all)) {
        # c = c + 1
        metabolite_class <- gsub(".txt","", metabolite_files[metab_list_nr]) # (string) stoftest name
        metab_list <- metab_list_all[[metab_list_nr]]
        # send to function
        make_plots(metab_list, metabolite_class, pt, zscore_cutoff, xaxis_cutoff, ThisProbScore, ratios_cutoff, pt_nr)
        # if (c%%1==0){
        cat(paste0(metab_list_nr, " "))
        # }
      }
      print(paste0("patientnr: ", pt_nr)) # remove later
      
      k <- dev.off()
      
    } # end for patient in first for-loop
    
    # Second loop: IEM plots
    for (pt_nr in 1:length(patient_list2)) {
      pt <- patient_list2[pt_nr]
      top_IEM <- c()
      ProbScore_top_IEM <- c()
      integer_list <- c(1:top)
      # Select the metabolites that are associated with the top highest scoring IEM, for each patient
      # disRank is from step 4: the algorithm. The lower the value, the more likely.
      IEMs <- disRank[disRank[[(pt_nr+1)]] %in% integer_list, ][[1]] # why +1? First column is Disease
      
      #### NB: this goes wrong, disRank has 22 columns for 21 patients, while patient_list has 23 patients ####
      # this object is generated in step 4 (algorithm)
      # add check on column name and patient name, they should match. Is that a unit test?
      
      for (IEM in IEMs) {
        ProbScore_IEM <- ProbScore0[which(ProbScore0$Disease==IEM), (pt_nr+1)]
        if (ProbScore_IEM >= threshold_IEM){
          top_IEM <- c(top_IEM, IEM)
          ProbScore_top_IEM <- c(ProbScore_top_IEM, ProbScore_IEM)
        }
      }
      l <- length(top_IEM)
      
      #If ProbScore_top_IEM is an empty list, don't continue to make_plots.
      if (length(top_IEM)==0){
        # If no Prob scores were above set threshold (5), send to log file
        cat(paste0("\n\n**** Note that this patient had no ProbScores higher than ", threshold_IEM,". 
                   Therefore, this pdf was not made:\t ",pt,"_top0 \n"))
      } else {
        # make plots of top 5 diseases
        cat(paste0("\n", "For ", pt,", done with plot nr. "))
        pdf(paste0(pdf_dir,"/", pt, "_top" , l , ".pdf"), onefile = TRUE,
            width = plot_width, height = plot_height) # create the PDF device
        # Sorting from high to low, both ProbScore_top_IEM as well as top_IEM.
        ind <- order(-ProbScore_top_IEM)
        ProbScore_top_IEM_sorted <- ProbScore_top_IEM[ind]
        top_IEM_sorted <- top_IEM[ind]
        # getting metabolites for each top_IEM disease exactly like in metab.list0
        dis_list <- Expected_red[Expected_red$Disease %in% top_IEM_sorted,]
        dis_list <- setDT(dis_list, key = "Disease")[top_IEM_sorted]
        dis_list$Disease <- factor(dis_list$Disease, levels=unique(dis_list$Disease))
        #disnames <- unique(dis_list$Disease)
        dis_list_all <- split(dis_list, f = dis_list$Disease)
        # d = 0
        # forloop over metab.list0 homolog
        for (dis_nr in 1:length(dis_list_all)) {
          # d = d + 1
          dis <- dis_list_all[[dis_nr]] # lots of duplicate info in dis
          disease <- as.character(dis[[1]][1]) # same as stoftest for normal plots
          dis <- dis[ , -c(1,4)] # same as metab.list in normal plots
          ThisProbScore <- ProbScore_top_IEM_s[dis_nr]
          #cat(paste0("d  ",d,"\n",disease,"\t",ThisProbScore))
          make_plots(dis, disease, pt, zscore_cutoff, xaxis_cutoff, ThisProbScore, ratios_cutoff, 0)
        }
      k <- dev.off()
      } # end make plots 
      
    } # end for patient in second for-loop
    
  } # end for metabolite_dir
  
} # end if violin = 1
      
# original code for violin plots:
      # ptcount <- 0
      # lapply(patient_list, function(pt) {
      #   # for each patient, go into the for loop for as many text files there are.
      #   c = 0
      #   if (ptcount <= nrpat) {
      #     ptcount <<- ptcount + 1
      #   } else {
      #     ptcount <<- 0
      #   }
      #   
      #   # IEM violin plots:
      #   if (endsWith(pt, "_IEM")) {
      #     n <<- n + 1
      #     cat(paste0("n",n))
      #     top_IEM <- c()
      #     ProbScore_top_IEM <- c()
      #     integer_list <- c(1:top)
      #     # Select the metabolites that are associated with the top 5 highest scoring IEM, for each patient
      #     IEMs <- disRank[disRank[[n+1]] %in% integer_list,][[1]]
      #     for (IEM in IEMs) {
      #       ProbScore_IEM <- ProbScore0[which(ProbScore0$Disease==IEM),(n+1)]
      #       if (ProbScore_IEM>=threshold_IEM){
      #         top_IEM <- c(top_IEM, IEM)
      #         ProbScore_top_IEM <- c(ProbScore_top_IEM, ProbScore_IEM)
      #       }
      #     }
      #     l <- length(top_IEM)
      #     
      #     #If ProbScore_top_IEM is an empty list, don't continue to make_plots.
      #     if (length(top_IEM)==0){
      #       # If no Prob scores were above set threshold (5), send to log file
      #       cat(paste0("\n\n**** Note that this patient had no ProbScores higher than ",threshold_IEM,". Therefore, this pdf was not made:\t ",pt,"_top0 \n"))
      #     } else {
      #       # make plots of top 5 diseases
      #       cat(paste0("\n","For ",pt,", done with plot nr. "))
      #       pdf(paste0(pdf_dir,"/", pt, "_top" , l , ".pdf"), onefile = TRUE,
      #           width = plot_width, height = plot_height) # create the PDF device
      #       # Sorting from high to low, both ProbScore_top_IEM as well as top_IEM.
      #       ind <- order(-ProbScore_top_IEM)
      #       ProbScore_top_IEM_s <- ProbScore_top_IEM[ind]
      #       top_IEM_s <- top_IEM[ind]
      #       # getting metabolites for each top_IEM disease exactly like in metab.list0
      #       dis_list <- Expected_red[Expected_red$Disease %in% top_IEM_s,]
      #       dis_list <- setDT(dis_list, key = "Disease")[top_IEM_s]
      #       dis_list$Disease <- factor(dis_list$Disease, levels=unique(dis_list$Disease))
      #       #disnames <- unique(dis_list$Disease)
      #       dis_list0 <- split(dis_list, f = dis_list$Disease)
      #       d = 0
      #       # forloop over metab.list0 homolog
      #       for (dis in dis_list0){
      #         d = d + 1
      #         
      #         disease <- as.character(dis[[1]][1]) # same as stoftest for normal plots
      #         dis <- dis[,-c(1,4)] # same as metab.list in normal plots
      #         #names(dis) <- gsub("HMDB.code", "HMDB_code", gsub("Metabolite", "HMDB_name", names(dis) ) )
      #         ThisProbScore <- ProbScore_top_IEM_s[d]
      #         #cat(paste0("d  ",d,"\n",disease,"\t",ThisProbScore))
      #         make_plots(dis,disease,pt,zscore_cutoff,xaxis_cutoff,ThisProbScore,ratios_cutoff,0)
      #       }
      #       check.lists <- FALSE
      #       k <- dev.off()
      #     }
      #     
      #   } else {
      #     
      #     # normal violin plots:
      #     cat(paste0("\n","For ",pt,", done with plot nr. "))
      #     
      #     if (startsWith(pt,"all")){
      #       # overview plots
      #       pdf(paste0(pdf_dir,"/", pt, "_patienten_overview.pdf"),onefile = TRUE,
      #           width = plot_width, height = plot_height) # create the PDF device
      #     } else {
      #       # patient plots
      #       pdf(paste0(pdf_dir,"/", pt, ".pdf"),onefile = TRUE,
      #           width = plot_width, height = plot_height) # create the PDF device
      #     }
      #     
      #     for (metab.list in metab.list0){
      #       ThisProbScore = 0 # means that this is no IEM plot
      #       c = c + 1
      #       stoftest <- gsub(".txt","",stofgroup_files[c]) # (string) stoftest name
      #       
      #       # send to function
      #       make_plots(metab.list, stoftest, pt, zscore_cutoff, xaxis_cutoff,ThisProbScore,ratios_cutoff,ptcount)
      #       if (c%%1==0){
      #         cat(paste0(c," "))
      #       }
      #     }
      #     print(paste0("patientnr: ",ptcount))
      #     
      #     k <- dev.off()
      #   }
      # }
      # ) # end of lapply
      

  # wrap up:
  writeLines(isobarics_txt, paste0(path_output, "/", run_name, "/isobarics.txt"), sep = "\n")
  outputfiles <- list.files(path=paste0(path_output, "/", run_name), pattern="*.pdf", full.names=FALSE, recursive=FALSE)
  if (exists("stofgroup_files") & exists("metab.list1") & (length(outputfiles) >= (nrpat+2))) {
    cat("\n\n### Step 5 # Make the violin plots is done. \n")
  } else {
    cat("\n\n**** Error: Could not make all violin plots or output folder already existed. pdf's made: \n")
  }
  #cat("### Step 5 # Make the violin plots is done.\n")
  cat("\nviolin plots pdf files made: ")
  cat(paste0("\n",outputfiles))
}


Sys.sleep(w)
#############################
cat(paste0("\n\nAll steps are executed, find output files here:\n -> ",output_dir))
end_time <- Sys.time()
cat("\n\nRun ended, end time: \t")
end_time-start_time
sink()
able_to_copy <- file.copy("log.txt",paste0(output_dir, "/log_",run_name,".txt"))
if (able_to_copy) {
  cat(paste0("\nlog file successfully copied to:\n -> ",output_dir,"\n\n"))
} else {
  cat("\n---- Warning: Could not copy log file. Look in working directory folder for a log.txt file. \n")
  file.copy("log.txt",paste0(output_dir, "/log_",run_name,"_____read_warning.txt"))
}
beep(1)
