# With the output of the DIMS pipeline, this script generates probability scores
# for each patient. In addition, it provides visual support with violin plots 
# of the DIMS measurements for the lab specialists.


#packused <- list.functions.in.file("dIEM_violin_pipeline.R", alphabetic = TRUE)




library(beepr)
library(dplyr)
library(reshape2)
library(data.table)
library(openxlsx)
library(ggplot2)
library(gghighlight)

rm(list = ls())

shorter <- 0
low_memory <- 1
sink(file="log.txt")
#############################
########## STEP # 0 #########      config settings
#############################


source("config.R")


#Use example config settings if no config.R file could be loaded
if (!exists("run_name")) {
  # folder name in which all output will be written
  run_name = "00_default_PLRUN10"
  # binary variable: run function, yes(1) or no(0)
  algorithm = 1
  ratios = 1
  violin = 1
  path_output = "~/Documents/dIEM/"
  # path: to DIMS excel file
  #path_DIMSfile = "/Volumes/LAB/metab/Metabolomics/Research Metabolic Diagnostics/Metabolomics Projects/Projects 2019/Project 2019_006 Saoudi Arabie/Bioinformatics_RUN3/2020-02-22_SA_RUN3.xlsx"
  path_DIMSfile = "~/Documents/dIEM/RES_PL_20200907_Diagnosis2017_RUN10_5ppm.xlsx"
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
  top_diseases = 5
  top_metab = 20
  
  cat("loading config.R not successful, now sample config settings are used.")
}
s <- 1 #suffix for plot filenames
#############################
########## STEP # 1 #########      Preparation
############################# in: run_name, path_DIMSfile, header_row ||| out: output_dir, DIMS
#############################
# create new output folder and set as output directory.
dir.create(file.path(path_output, run_name))
output_dir <- paste0(path_output,"/",run_name)
cat (paste("directory made:",output_dir, sep=""))
able_to_copy <- file.copy("config.R",output_dir)
cat(paste0("\n config file successfully copied to ",output_dir," = ",able_to_copy))
# Load the excel file.
dimsxls <- readWorkbook(xlsxFile = path_DIMSfile, sheet = 1, startRow = header_row)
cat (paste("\n loaded",path_DIMSfile, sep=""))
cat("\n ### Step 1 # Preparation is done.\n")
beep("coin")

#############################
########## STEP # 2 #########      Edit DIMS data 
############################# in: DIMSxls ||| out: Data, nrcontr, nrpat
#############################
# It edits a few column names, removes irrelevant columns.
# Input: 
# The xlsx file that comes out of the pipeline (v.2.0.0) with format:
# [plots] [C] [P] [summary columns] [C_Zscore] [P_Zscore]

# Output: 
# Data dataframe 
# "_CSV.csv" file that is suited for the algorithm in shiny.



dims2 <- dimsxls
# Calculate the number of C's and P's in column names to extract the following numbers:
nrcontr <- length(grep("C",names(dims2)))/2   # Number of control samples
nrpat <- length(grep("P",names(dims2)))/2     # Number of patient samples
if (nrcontr + nrpat != length(grep("_Zscore", names(dims2)))) {
  cat("Error: there aren't as many intensities listed as Zscores")
}
cat(paste0(nrcontr, " controls \n",nrpat," patients \n"))

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
cat("### Step 2 # Edit dims data is done.\n")


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
# 
# "_Ratios_CSV.csv" file, same file as above, but with ratio rows added.
dims3 <- dims2
#rm(dims2)

if (ratios == 1) { # ratios in settings is 1
cat(paste0("\n loading ratios file: \n",path_ratios))
RatioInput<-read.csv(path_ratios,sep=';',stringsAsFactors=FALSE)

# Prepare empty data frame to fill with ratios
Ratios<-setNames(data.frame(matrix(ncol=ncol(dims3),nrow=nrow(RatioInput))),colnames(dims3))
Ratios[,1:2]<-RatioInput[,1:2]
### idea: test without log10 , look into expected for ratios
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
    #Ratios[rows,Zscores]<-(Ratios[rows,(Zscores-nrcontr-nrpat-2)]-Ratios[rows,(nrcontr+nrpat+3)])/Ratios[rows,(nrcontr+nrpat+4)]
  }
}

# Add rows of the ratio hmdb codes to the data of zscores from the pipeline.
Combined<-rbind.data.frame(Ratios,dims3)

# explain.
Zscore <- Combined[,c(1:2,(2*nrcontr+nrpat+5):(2*(nrcontr+nrpat)+4))]
#Zscoretest <- Combined[,c(1:2,grep("P",colnames(Data))[-c(1:nrpat)] )]
Zscore_all <- Combined[,c(1:2,(nrcontr+nrpat+5):(2*(nrcontr+nrpat)+4))]
#Zscore_alltest <- Combined[,c(1:2,grep("P|C",colnames(Data))[-c(1:(nrpat+nrcontr))] )]
# _Ratios: full dataframe, with intensities, intensity ratios & zscores and zscores ratios
#write.csv(Combined, file=paste(output_dir,"/",run_name,"_Ratios_CSV.csv",sep=""))
# _inputshiny: only zscores and zscores of ratio hmdb's, to be used as input for algorithm shiny app
write.table(Zscore,file=paste(output_dir,"/",run_name,"inputshiny_CSV.csv",sep=""),quote=FALSE,sep=";",row.names=FALSE)
cat("### Step 3 # Calculate ratios is done.\n")
if (low_memory == 1) {
  rm(Zscore_all,dims2,dims3,dimsxls,Combined)
}
}


#############################
########## STEP # 4 #########      Run the algorithm
############################# in: algoritm, path_expected, Zscore ||| out: ProbScore0 (+file)
#############################
# This script....

# Zscore <- ObsZscore
if (algorithm == 1){
# Load data
cat(paste0("\n loading expected file: \n",path_expected))
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

#Exp_Metabscore <- merge(x=Expected, y=Metabscore, by = "HMDB.code")
#Exp_Metabscore[30:ncol(Exp_Metabscore)] <- lapply(Exp_Metabscore[30:ncol(Metabscore)], function(x) ifelse((Exp_Metabscore$Change=="Increase")&(x<=0), 0, x)) 
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
# col names aanpassen! van _Zscore naar niks


write.xlsx(ProbScore0, paste0(output_dir,"/",run_name,"algoritme_output.xlsx"))

cat("### Step 4 # Run the algorithm is done.\n")
if (low_memory == 1) {
  rm(Rank, disRank, Exp_Metabscore,Exp_Rank,Exp_Zscores,Expected,Exp_Zscores0,ProbScore,dup,uni,Wscore,Ratios)
}
}

beep("coin")
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


patient_list <- names(summed)
patient_list[1] <- "all"
stofgroup_files <- list.files(path=path_txtfiles, pattern="*.txt", full.names=FALSE, recursive=FALSE)
metab.list0 <- list()
index = 0
cat("making plots from the input files:")
for (infile in stofgroup_files) {
  index = index + 1
  metab.list1 <- unique(read.table(paste0(path_txtfiles,"/",infile), sep = "\t", header = TRUE, quote=""))
  metab.list0[[index]] <- metab.list1
  cat(paste0("\n",infile))
}

make_plots <- function(metab.list,stoftest,pt,zscore_cutoff,xaxis_cutoff) {
  i_tot <- nrow(metab.list)
  #moi <- summed[summed$HMDB_code %in% metab.list[[1]], ]
  #moi <- summed[summed$HMDB_code %in% metab.list$HMDB_code, ]
  #metab.list <- metab.list0[[1]]
  joined <- inner_join(metab.list, summed, by = "HMDB_code")
  
  moi <- joined[-2]
  moi_names <- joined[-1]
  #moi_names <- merge(metab.list, summed, by = "HMDB_code")[-1]
  names(moi_names) <- gsub("_Zscore", "", names(moi_names))
  names(moi) <- gsub("_Zscore", "", names(moi))
  moi_m <- melt(moi_names, id.vars = "HMDB_name")
  moi_m<- aggregate(HMDB_name ~ variable+value, moi_m, paste0, collapse = "\n")
  group_highZ <- moi_m %>%
    group_by(value) %>% 
    filter(value > zscore_cutoff) %>%
    ungroup()
  
  moi_m_max20 <- moi_m
  moi_m_max20$value <- as.numeric(lapply(moi_m$value, function(x) ifelse(x > xaxis_cutoff, as.numeric(xaxis_cutoff), x)))
  
  group_highZ_max20 <- moi_m_max20 %>%
    group_by(value) %>% 
    filter(value > zscore_cutoff) %>%
    ungroup()
  
  if (pt!="all"){
    pt_colname <- pt
    pt_data <- moi_m[which(moi_m$variable==pt_colname),]
    pt_data_max20 <- moi_m_max20[which(moi_m_max20$variable==pt_colname),]
    pt_values <- pt_data$value
    colors <- c("#b7355d", "#7435b7", "#3592b7","#35b740","#e4a010","#ffff00")
    
    plot_height <- 80 * i_tot
    file_png <- paste0(output_dir,"/", pt, "_",stoftest,"_",i,".png")
    
    g <- ggplot(moi_m_max20, aes(x=value, y=HMDB_name, color = value))+
      geom_violin(scale="width")+
      geom_point(data = pt_data_max20, aes(color=pt_data$value),size = 3.5,shape=21, fill="white")+
      geom_jitter(data = group_highZ_max20, aes(color=group_highZ$value), size = 1.5, position = position_dodge(1.5))+ #,colour = "#3592b7" 
      scale_fill_gradientn(colors = colors,values = NULL,space = "Lab",na.value = "grey50",guide = "colourbar",aesthetics = "colour")+
      labs(x = "Z-scores",y = "Metabolites",title = paste0("Results for patient ",pt), subtitle = paste0(stoftest,"\nZ-score > ",zscore_cutoff))+
      geom_vline(xintercept = 2, col = "grey", lwd = 0.5,lty=2)
    print(g)
  }
  
  if (pt=="all"){
    colors <- c("#b7355d", "#7435b7", "#3592b7","#35b740","#e4a010","#ffff00")
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

lapply(patient_list, function(pt) {
  # for each patient, go into the for loop for as many stofgroups there are.
  cat(paste0("\n","For ",pt,", done with plot nr. "))
  i_tot <- 48
  plot_height <- 0.5 * i_tot
  if (pt=="all"){
    pdf(paste0(output_dir,"/", pt, "_overview_",s,".pdf"),onefile = TRUE,
        width = 7, height = plot_height) # create the PDF device
  } else {
  pdf(paste0(output_dir,"/", pt, "_stoftesten_",s,".pdf"),onefile = TRUE,
      width = 5, height = plot_height) # create the PDF device
  }
  c = 0
  # Filter summed on the metabolites of interest (moi)
  for (metab.list in metab.list0){
    c = c + 1
    stoftest <- gsub(".txt","",stofgroup_files[c])
    metab.list$HMDB_name <- gsub(';','\n',metab.list$HMDB_name)
    make_plots(metab.list, stoftest, pt, zscore_cutoff, xaxis_cutoff)
    if (c%%3==0){
      cat(paste0(c," .."))
    }
  }
  s <- dev.off()
})
cat("### Step 5 # Make the violin plots is done.\n")

}


#############################
cat(paste0("\n All steps are executed, find output files here:\n",output_dir))
sink()
able_to_copy <- file.copy("log.txt",paste0(output_dir, "/log_",run_name,".txt"))
cat(paste0("\n log file successfully copied to ",output_dir," = ",able_to_copy))
beep(1)
