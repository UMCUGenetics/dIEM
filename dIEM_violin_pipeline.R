#For untargeted metabolomics, this tool calculates probability scores for 
# metabolic disorders. In addition, it provides visual support with violin plots 
# of the mass spectrometry (DI-HRMS) measurements for the lab specialists.
# Input needed: 
# 1. excel file in which metabolites are listed with their intensities among
#    controls (with C in samplename) and patients (with P in samplename) and their
#    corresponding zscores. 
# 2. All files from github: https://github.com/UMCUGenetics/dIEM

library(beepr)
library(dplyr)
library(reshape2)
library(data.table)
library(openxlsx)
library(ggplot2)
library(gghighlight)
library(sys)

rm(list = ls())
w <- 0.5 # seconds that system waits in between steps
shorter <- 0 # shorter list of patients
low_memory <- 0 # if RStudio crashes during script
top <- 5 # number of diseases that score highest in algorithm to plot
threshold_IEM = 5 # probability score cut-off for plotting the top diseases
ratios_cutoff = -5 # z-score cutoff of axis on the left for top diseases


#############################
########## STEP # 0 #########      load config settings
#############################
sink(file="log.txt")
source("config.R")
start_time <- Sys.time()
cat(paste0("Running dIEM violin pipeline, start time: \t",start_time,"\n"))
if (exists("run_name")) {
  cat("\nThe config file is succesfully loaded from working directory.")
} else {
  cat("\n**** Error: Could not find a config file. please check if working directory is set in 'Session'. \n")
}


#############################
########## STEP # 1 #########      Preparation
############################# in: run_name, path_DIMSfile, header_row ||| out: output_dir, DIMS
#############################
# create new output folder and set as output directory.
dir.create(file.path(path_output, run_name))
output_dir <- paste0(path_output,"/",run_name)
able_to_copy <- file.copy("config.R",output_dir)
if (able_to_copy) {
  cat(paste0("\nconfig file successfully copied to:\n -> ",output_dir))
} else {
  # If the config file could not be copied, the run name probably already exists.
  cat("\n---- Warning: please use a new run name for every run. Now, a time-stamp is added to the runname. \n")
  # Thus, add timestamp to run name.
  run_name <- paste0(run_name,"_",gsub("CET","",gsub(" |-|:", "",Sys.time())))
  dir.create(file.path(path_output, run_name))
  output_dir <- paste0(path_output,"/",run_name)
  able_to_copy <- file.copy("config.R",output_dir)
  cat(paste0("\nconfig file successfully copied to ",output_dir," = ",able_to_copy))
}

# Load the excel file.
dimsxls <- readWorkbook(xlsxFile = path_DIMSfile, sheet = 1, startRow = header_row)
if (exists("dimsxls")) {
  cat(paste0("\nThe excel file is succesfully loaded:\n -> ",path_DIMSfile))
} else {
  cat(paste0("\n\n**** Error: Could not find an excel file. Please check if path to excel file is correct in config.R:\n -> ",path_DIMSfile,"\n"))
}
beep("coin")

Sys.sleep(w)
#############################
########## STEP # 2 #########      Edit DIMS data 
############################# in: DIMSxls ||| out: Data, nrcontr, nrpat
#############################
# It edits a few column names, removes irrelevant columns.
# Input: 
# The xlsx file that comes out of the pipeline (v.2.0.0) with format:
# [plots] [C] [P] [summary columns] [C_Zscore] [P_Zscore]

# Output: 
# dims2 dataframe 
# "_CSV.csv" file that is suited for the algorithm in shiny.

dims2 <- dimsxls
# Calculate the number of C's and P's in column names to extract the following numbers:
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
# Rename the columns from..to..
names(dims2) <- gsub("avg.ctrls", "Mean_controls", gsub("sd.ctrls", "SD_controls", names(dims2) ) )
names(dims2) <- gsub("HMDB_code", "HMDB.code", gsub("HMDB_name", "HMDB.name", names(dims2) ) )

#first, select the intensity columns by all cols minus nrsamples
nrsamples = nrcontr + nrpat
beginZscores = ncol(dims2) - nrsamples
i <- c(3:(nrsamples+2))
#change the intensities to numeric values
dims2[, i] <- sapply(dims2[, i], as.numeric)



if (shorter==1){
  pos_ctr <- 13 # number of positive ctrs that are between controls and patients in df
  number_patients <- 15 # number of patients that are sampled
  dims2 <- dims2[,c(1:2,
                         3:(nrcontr+2), 
                         (3+nrcontr+pos_ctr):(2+ nrcontr +pos_ctr+number_patients),
                 (3+ nrcontr +nrpat):(4+ 2*nrcontr +nrpat),
                 (5+ 2*nrcontr +nrpat+pos_ctr):(4+ 2*nrcontr +nrpat+ pos_ctr +number_patients)
                 
                         )] #(5+ 2*nrcontr + 2*pos_ctr +number_patients):(4+ 2*nrcontr + 2*pos_ctr + 2*number_patients)
  nrcontr <- length(grep("C",names(dims2)))/2   # Number of control samples
  nrpat <- length(grep("P",names(dims2)))/2     # Number of patient samples
  i <- c(3:(nrsamples+2))
}
if (exists("dims2") & (length(dims2)<length(dimsxls))) {
  cat("\n### Step 2 # Edit dims data is done.\n")
} else {
  cat("\n**** Error: Could not execute step 2 \n")
}
#cat("### Step 2 # Edit dims data is done.\n")


Sys.sleep(w)
#############################
########## STEP # 3 #########      Calculate ratios
############################# in: ratios, path_ratios, dims2, nrcontr, nrpat ||| out: Zscore (+file)
############################# 
# This script loads the file with Ratios (path_ratios) and calculates 
# the ratios of the intensities of the given metabolites. It also calculates
# Zscores based on the avg and sd of the ratios of the controls.

# Input:
# The dataframe with intenstities and Zscores of controls and patients:
# [HMDB.code] [HMDB.name] [C] [P] [Mean_controls] [SD_controls] [C_Zscore] [P_Zscore]

# Output:
# "_CSV.csv" file that is suited for the algorithm, with format:
# "_Ratios_CSV.csv" file, same file as above, but with ratio rows added.
dims3 <- dims2
if (ratios == 1) { # ratios in settings is 1
cat(paste0("\nloading ratios file:\n ->  ",path_ratios,"\n"))
RatioInput<-read.csv(path_ratios,sep=';',stringsAsFactors=FALSE)

# Prepare empty data frame to fill with ratios
Ratios<-setNames(data.frame(matrix(ncol=ncol(dims3),nrow=nrow(RatioInput))),colnames(dims3))
Ratios[,1:2]<-RatioInput[,1:2]
### idea: test without log10, look into expected for ratios
for (controls in c(3:(nrcontr+2),(nrcontr+3):(nrcontr+nrpat+2))) {
  Ratios[1,controls]<-log10(dims3[which(dims3[,1]=='HMDB00159'),controls]/dims3[which(dims3[,1]=='HMDB00158'),controls])
  Ratios[2,controls]<-log10(dims3[which(dims3[,1]=='HMDB00161'),controls]/dims3[which(dims3[,1]=='HMDB00182'),controls])
  Ratios[3,controls]<-log10(dims3[which(dims3[,1]=='HMDB00161'),controls]/
                              (dims3[which(dims3[,1]=='HMDB00159'),controls]+dims3[which(dims3[,1]=='HMDB00158'),controls]))
  Ratios[4,controls]<-log10(dims3[which(dims3[,1]=='HMDB00062'),controls]/
                              (dims3[which(dims3[,1]=='HMDB00222'),controls]+dims3[which(dims3[,1]=='HMDB00848'),controls]))
  Ratios[5,controls]<-log10((dims3[which(dims3[,1]=='HMDB00222'),controls]+dims3[which(dims3[,1]=='HMDB05065'),controls])
                            /dims3[which(dims3[,1]=='HMDB00201'),controls])
  Ratios[6,controls]<-log10(dims3[which(dims3[,1]=='HMDB02014'),controls]/dims3[which(dims3[,1]=='HMDB00201'),controls])
  Ratios[7,controls]<-log10(dims3[which(dims3[,1]=='HMDB00791'),controls]/dims3[which(dims3[,1]=='HMDB00201'),controls])
  Ratios[8,controls]<-log10(dims3[which(dims3[,1]=='HMDB13127'),controls]/dims3[which(dims3[,1]=='HMDB00201'),controls])
  Ratios[9,controls]<-log10(dims3[which(dims3[,1]=='HMDB00064'),controls]/dims3[which(dims3[,1]=='HMDB00562'),controls])
  Ratios[10,controls]<-log10(dims3[which(dims3[,1]=='HMDB00824'),controls]/dims3[which(dims3[,1]=='HMDB00696'),controls])
  Ratios[11,controls]<-log10(dims3[which(dims3[,1]=='HMDB00159'),controls]/
                               (dims3[which(dims3[,1]=='HMDB00824'),controls]+dims3[which(dims3[,1]=='HMDB00222'),controls]))
  Ratios[12,controls]<-log10(dims3[which(dims3[,1]=='HMDB00118'),controls]/dims3[which(dims3[,1]=='HMDB00763'),controls])
  Ratios[13,controls]<-log10(dims3[which(dims3[,1]=='HMDB01325'),controls]/dims3[which(dims3[,1]=='HMDB06831'),controls])
  Ratios[14,controls]<-log10(dims3[which(dims3[,1]=='HMDB00791'),controls]/dims3[which(dims3[,1]=='HMDB00651'),controls])
}

# Calc means and SD's of the calculated ratios, add them in 2 columns in ratio df.
for (calc in 1:nrow(Ratios)) {
  Ratios[calc,(nrcontr+nrpat+3)]<-mean(as.numeric(Ratios[calc,3:(nrcontr+2)]))
  Ratios[calc,(nrcontr+nrpat+4)]<-sd(as.numeric(Ratios[calc,3:(nrcontr+2)]))
}

# Calc z-scores with the means and SD's 
for (Zscores in (nrcontr+nrpat+5):(2*(nrcontr+nrpat)+4)) {
  for (rows in 1:nrow(Ratios)) {
    Ratios[rows,Zscores]<-(Ratios[rows,(Zscores-nrcontr-nrpat-2)]-Ratios[rows,(nrcontr+nrpat+3)])/Ratios[rows,(nrcontr+nrpat+4)]
    }
}

# Add rows of the ratio hmdb codes to the data of zscores from the pipeline.
Combined<-rbind.data.frame(Ratios,dims3)

# Select only the cols with zscores of only the patients (Zscore)
Zscore <- Combined[,c(1:2,(2*nrcontr+nrpat+5):(2*(nrcontr+nrpat)+4))]
# And with the controls
Zscore_all <- Combined[,c(1:2,(nrcontr+nrpat+5):(2*(nrcontr+nrpat)+4))]
write.table(Zscore,file=paste(output_dir,"/inputshiny_",run_name,"_CSV.csv",sep=""),quote=FALSE,sep=";",row.names=FALSE)

# check whether Zscore and Zscore_all are as expected: implies success of calc ratios
if (exists("Combined") & (length(Zscore)<length(Zscore_all))) {
  cat("\n### Step 3 # Calculate ratios is done.\n")
} else {
  cat("\n**** Error: Could not calculate ratios. Check if path to ratios-file is correct in config.R. \n")
}

if (low_memory == 1) {
  rm(Zscore_all,dims2,dims3,dimsxls,Combined)
}
}

Sys.sleep(w)
#############################
########## STEP # 4 #########      Run the algorithm
############################# in: algoritm, path_expected, Zscore ||| out: ProbScore0 (+file)
#############################


# Zscore <- ObsZscore
if (algorithm == 1){
# Load data
cat(paste0("\nloading expected file:\n ->  ",path_expected,"\n"))
Expected<-read.csv(path_expected,sep=';',stringsAsFactors=FALSE)

# prepare dataframe scaffold Rank
Rank <- Zscore
r <- nrow(Rank)
# Fill df Rank with the ranks for each patient
for (patients in 3:ncol(Zscore)) {
  # number of positive zscores in patient
  pos <- length(Zscore[patients][(Zscore[patients]>0)==TRUE])
  # sort the column on zscore
  Rank <- Rank[order(-Rank[patients]),]
  # Rank all positive zscores highest to lowest
  Rank[1:pos,patients]<-as.numeric(ordered(-Rank[1:pos,patients]))
  # Rank all negative zscores lowest to highest
  Rank[(pos+1):r,patients]<-as.numeric(ordered(Rank[(pos+1):r,patients]))
}

# Calculate metabolite score, using the dataframes with only values, and later add the cols without values (1&2).
Exp_Zscores <- merge(x=Expected, y=Zscore, by.x = c("HMDB.code"), by.y = c("HMDB.code"))[-29]
Exp_Zscores0 <- Exp_Zscores

cat("setting some zscores to zero.\t")
Exp_Zscores0[which(Exp_Zscores0$Change=="Increase" & Exp_Zscores0$Dispensability=="Indispensable"),29:ncol(Exp_Zscores0)] <- lapply(Exp_Zscores0[which(Exp_Zscores0$Change=="Increase" & Exp_Zscores0$Dispensability=="Indispensable"),29:ncol(Exp_Zscores0)], function(x) ifelse(x<=1.6 , 0, x))
Exp_Zscores0[which(Exp_Zscores0$Change=="Decrease" & Exp_Zscores0$Dispensability=="Indispensable"),29:ncol(Exp_Zscores0)] <- lapply(Exp_Zscores0[which(Exp_Zscores0$Change=="Decrease" & Exp_Zscores0$Dispensability=="Indispensable"),29:ncol(Exp_Zscores0)], function(x) ifelse(x>=-1.2 , 0, x))
cat("done.\n")
Exp_Rank <- merge(x=Expected, y=Rank, by.x = c("HMDB.code"), by.y = c("HMDB.code"))[-29]
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

if (low_memory == 1) {
  rm(Rank, Exp_Metabscore,Exp_Rank,Exp_Zscores,Exp_Zscores0,ProbScore,dup,uni,Wscore,Ratios)
}
rm(wb)
}
Sys.sleep(w)





#############################
########## STEP # 5 #########      Make violin plots
############################# in: algorithm / Zscore, violin, nrcontr, nrpat, Data, path_textfiles, zscore_cutoff, xaxis_cutoff, top_diseases, top_metab, output_dir ||| out: pdf file
#############################

if ( algorithm == 0 ){
  #nrcontr <- Data
  Zscore <- Data[,c(1:2,(2*nrcontr+nrpat+5):(2*(nrcontr+nrpat)+4))]                                                         
}
if (violin == 1) {
#Edit the DIMS output Zscores of all patients in format:
# HMDB_code patientname1  patientname2
names(Zscore) <- gsub("HMDB.code","HMDB_code", names(Zscore) )
summed <- Zscore[-2] # remove the HMDB.name column
names(summed) <- gsub("_Zscore", "", names(summed)) # remove the _Zscore from columnnames

# Make a patient list so it can be looped over later in lapply.
patient_list <- names(summed)
patient_list[1] <- "alle"
# add list for IEM plots 
patient_list2 <- gsub("alle_IEM","alle_uitgerekt",paste0(patient_list, "_IEM"))
patient_list <- c(patient_list, patient_list2)
# Find all text files in the given folder, which contain metabolite lists of which
# each file will be a page in the pdf with violin plots.
stofgroup_files <- list.files(path=path_txtfiles, pattern="*.txt", full.names=FALSE, recursive=FALSE)
metab.list0 <- list()
index = 0
cat("making plots from the input files:")
# open the text files and add each to a list of dataframes (metab.list0)
for (infile in stofgroup_files) {
  index = index + 1
  metab.list1 <- unique(read.table(paste0(path_txtfiles,"/",infile), sep = "\t", header = TRUE, quote=""))
  metab.list0[[index]] <- metab.list1
  cat(paste0("\n",infile))
}
cat("\n")
# reduce the columns from the expected df for further use and remove double 
# metabolites per disease.
Expected_red <- Expected[,c(5,14,13,18)]
Expected_red <- Expected_red[!duplicated(Expected_red[,c(1,2)]),]

make_plots <- function(metab.list,stoftest,pt,zscore_cutoff,xaxis_cutoff,ThisProbScore,ratios_cutoff) {
  i_tot <- nrow(metab.list)
  # Filter summed on the metabolites of interest (moi)
  joined <- inner_join(metab.list, summed, by = "HMDB_code")
  moi <- joined[,-2]
  moi_names <- joined[,-1]
  # remove "_Zscore" from col names
  names(moi_names) <- gsub("_Zscore", "", names(moi_names))
  names(moi) <- gsub("_Zscore", "", names(moi))
  # melt the dataframe because that is easily plotted.
  moi_m <- melt(moi_names, id.vars = "HMDB_name")
  # collapse all compounds with the same zscore, as they are isobaric compounds
  # (= different molecular formula, same nominal mass)
  moi_m<- aggregate(HMDB_name ~ variable+value, moi_m, paste0, collapse = "\n")
  # make selection of scores higher than the cut-off that will be colored according
  # to their values. They will be plotted according to their values in moi_m
  group_highZ <- moi_m %>%
    group_by(value) %>% 
    filter(value > zscore_cutoff) %>%
    ungroup()
  # change all values above the xaxis cutoff to the cut-off value. So that
  # all plotted data is within a shorter range on the x-axis.
  moi_m_max20 <- moi_m
  moi_m_max20$value <- as.numeric(lapply(moi_m$value, function(x) ifelse(x > xaxis_cutoff, as.numeric(xaxis_cutoff), x)))
  # Because the range of diagnostic ratios is often far below zero, change all
  # values below the xaxis cutoff to the cut-off value to prevent a very long
  # x-axis in the plots with other metabolites (in IEM plots).
  if (endsWith(pt,"_IEM")){
    moi_m_max20$value <- as.numeric(lapply(moi_m_max20$value, function(x) ifelse(x < ratios_cutoff, as.numeric(ratios_cutoff), x)))
  }
  # 
  group_highZ_max20 <- moi_m_max20 %>%
    group_by(value) %>% 
    filter(value > zscore_cutoff) %>%
    ungroup()

  if (!startsWith(pt,"all")){
    # patient one by one
    pt <- gsub("_IEM", "", pt)
    pt_colname <- pt
    pt_data <- moi_m[which(moi_m$variable==pt_colname),]
    pt_data_max20 <- moi_m_max20[which(moi_m_max20$variable==pt_colname),]
    pt_values <- pt_data$value
    #colors <- c("#b7355d", "#7435b7", "#3592b7","#35b740","#f8f32b","#e4a010")
    colors <- c("#22E4AC", "#00B0F0", "#504FFF","#A704FD","#F36265","#DA0641")
    #             green     blue      blue/purple purple    orange    red
    plot_height <- 80 * i_tot
    file_png <- paste0(output_dir,"/", pt, "_",stoftest,"_",i,".png")
    if (ThisProbScore==0){
      
      g <- ggplot(moi_m_max20, aes(x=value, y=HMDB_name, color = value))+
        geom_violin(scale="width")+
        geom_point(data = pt_data_max20, aes(color=pt_data$value),size = 3.5,shape=21, fill="white")+
        geom_jitter(data = group_highZ_max20, aes(color=group_highZ$value), size = 1.5, position = position_dodge(1.5))+ #,colour = "#3592b7" 
        scale_fill_gradientn(colors = colors,values = NULL,space = "Lab",na.value = "grey50",guide = "colourbar",aesthetics = "colour")+
        labs(x = "Z-scores",y = "Metabolites",title = paste0("Results for patient ",pt), subtitle = paste0(stoftest,"\nZ-score > ",zscore_cutoff))+
        geom_vline(xintercept = 2, col = "grey", lwd = 0.5,lty=2)
    } else {
      
      g <- ggplot(moi_m_max20, aes(x=value, y=HMDB_name, color = value))+
        geom_violin(scale="width")+
        geom_point(data = pt_data_max20, aes(color=pt_data$value),size = 3.5,shape=21, fill="white")+
        geom_jitter(data = group_highZ_max20, aes(color=group_highZ$value), size = 1.5, position = position_dodge(1.5))+ #,colour = "#3592b7" 
        scale_fill_gradientn(colors = colors,values = NULL,space = "Lab",na.value = "grey50",guide = "colourbar",aesthetics = "colour")+
        labs(x = "Z-scores",y = "Metabolites",title = paste0("Algorithm results for patient ",pt), subtitle = paste0("Disease: ",stoftest,"\nProbability Score = ",format(round(ThisProbScore, 2), nsmall = 2)))+
        geom_vline(xintercept = 2, col = "grey", lwd = 0.5,lty=2)
    }
    print(g)
  }
  
  if (startsWith(pt,"all")){
    # overview plots
    colors <- c("#b7355d", "#7435b7", "#3592b7","#35b740","#f8f32b","#e4a010")
    plot_height <- 80 * i_tot + (i_tot/(nrow(group_highZ) *0.25))*5
    plot_width <- 800+(max(group_highZ$value)/2)
    file_png <- paste0(output_dir,"/", pt, "_",stoftest,"_fullaxis_",i,".png")
    
    g <- ggplot(moi_m, mapping = aes(x=value, y=HMDB_name))+
      geom_violin(scale="width")+
      geom_jitter(data = group_highZ,aes(color = variable), size = 2.5, position = position_dodge(0.8))+ 
      labs(x = "Z-scores",y = "Metabolites",title = "Overview plot", subtitle = paste0(stoftest,"\nZ-score > ",zscore_cutoff), color = "patients")+
      geom_vline(xintercept = 2, col = "grey", lwd = 0.5,lty=2)
    print(g)

    
    plot_height <- 80 * i_tot + (i_tot/(nrow(group_highZ) *0.25))*5
    plot_width <- 800
    file_png <- paste0(output_dir,"/", pt, "_",stoftest,"_max20_",i,".png")
    
    g <- ggplot(moi_m_max20, mapping = aes(x=value, y=HMDB_name))+
      geom_violin(scale="width")+
      geom_jitter(data = group_highZ_max20,aes(color = group_highZ$variable), size = 2.5, position = position_dodge(0.8))+ 
      labs(x = "Z-scores",y = "Metabolites",title = "Overview plot", subtitle = paste0(stoftest,"\nZ-score > ",zscore_cutoff), color = "patients")+
      geom_vline(xintercept = 2, col = "grey", lwd = 0.5,lty=2)
    print(g)

  }
}
n <- 0
i_tot <- 24
plot_height <- 0.40 * i_tot
plot_width <- 6
cat(paste0("\n(plot dimensions: ",plot_width, " x ", plot_height,")\n" ))

lapply(patient_list, function(pt) {
  # for each patient, go into the for loop for as many text files there are.
  cat(paste0("\n","For ",pt,", done with plot nr. "))
  c = 0
  # IEM violin plots:
  if (endsWith(pt,"_IEM")) {
    n <<- n + 1
    cat(paste0("n",n))
    top_IEM <- c()
    ProbScore_top_IEM <- c()
    integer_list <- c(1:top)
    # Select the metabolites that are associated with the top 5 highest scoring IEM, for each patient
    IEMs <- disRank[disRank[[n+1]] %in% integer_list,][[1]]
    for (IEM in IEMs) {
      ProbScore_IEM <- ProbScore0[which(ProbScore0$Disease==IEM),(n+1)]
      if (ProbScore_IEM>=threshold_IEM){
        top_IEM <- c(top_IEM, IEM)
        ProbScore_top_IEM <- c(ProbScore_top_IEM, ProbScore_IEM)
      }
    }
    l <- length(top_IEM)
    pdf(paste0(output_dir,"/", pt, "_top" , l , ".pdf"),onefile = TRUE,
        width = plot_width, height = plot_height) # create the PDF device
    #If ProbScore_top_IEM is an empty list, don't continue to make_plots.
    if (length(top_IEM)==0){
      # If no Prob scores were above set threshold (5), send to log file
      cat(paste0(pt," had no ProbScores higher than ",threshold_IEM))
    } else {
      
      # Sorting from high to low, both ProbScore_top_IEM as well as top_IEM.
      ind <- order(-ProbScore_top_IEM)
      ProbScore_top_IEM_s <- ProbScore_top_IEM[ind]
      top_IEM_s <- top_IEM[ind]
      # getting metabolites for each top_IEM disease exactly like in metab.list0
      dis_list <- Expected_red[Expected_red$Disease %in% top_IEM_s,]
      dis_list <- setDT(dis_list, key = "Disease")[top_IEM_s]
      dis_list0 <- split(dis_list, f = dis_list$Disease)
      d = 0
      # forloop over metab.list0 homolog
      for (dis in dis_list0){
        d = d + 1
        cat(paste0("d",d))
        disease <- dis[[1]][1] # same as stoftest for normal plots
        dis <- dis[,-c(1,4)] # same as metab.list in normal plots
        names(dis) <- gsub("HMDB.code", "HMDB_code", gsub("Metabolite", "HMDB_name", names(dis) ) )
        ThisProbScore <- ProbScore_top_IEM_s[d]
        make_plots(dis,disease,pt,zscore_cutoff,xaxis_cutoff,ThisProbScore,ratios_cutoff)
      }
    }
  } else {
    # normal violin plots:
    if (startsWith(pt,"all")){
      # overview plots
      pdf(paste0(output_dir,"/", pt, "_patienten_overview.pdf"),onefile = TRUE,
          width = plot_width, height = plot_height) # create the PDF device
    } else {
      # patient plots
      pdf(paste0(output_dir,"/", pt, ".pdf"),onefile = TRUE,
          width = plot_width, height = plot_height) # create the PDF device
    }
    for (metab.list in metab.list0){
      ThisProbScore = 0 # means that this is no IEM plot
      c = c + 1
      stoftest <- gsub(".txt","",stofgroup_files[c]) # (string) stoftest name
      # generate a new line when ; is in metabolite name.
      metab.list$HMDB_name <- gsub(';','\n',metab.list$HMDB_name)
      # send to function
      make_plots(metab.list, stoftest, pt, zscore_cutoff, xaxis_cutoff,ThisProbScore,ratios_cutoff)
      if (c%%5==0){
        cat(paste0(c," "))
      }
    }
  }
  k <- dev.off()
  
})
outputfiles <- list.files(path=paste0(path_output,"/",run_name), pattern="*.pdf", full.names=FALSE, recursive=FALSE)
if (exists("stofgroup_files") & exists("metab.list1") & (length(outputfiles)>=(nrpat+2))) {
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
