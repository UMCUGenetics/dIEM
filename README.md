# dIEM with visualization
For untargeted metabolomics, this tool calculates probability scores for metabolic disorders. In addition, it provides visual support with violin plots of the mass spectrometry (DI-HRMS) measurements for the lab specialists.


### Installation

Download zip-file or open the terminal and clone the repository:


```bash
cd ~ # go to folder
git clone https://github.com/UMCUGenetics/dIEM.git
```

### Folder Structure
```
|─── data (libraries for the dIEM algorithm)
|─── stofgroups (metabolite lists for the violin plots)
```

### Usage

The tool can be run in RStudio:
  -  Setting the working directory in the same directory as the `config.R` file. 
  -  Running the `packages_installation.R` script.
  -  Editing the `config.R` file with the appropriate paths and variables.
  -  Running the main script: `dIEM_violin_pipeline.R`

If the main script makes RStudio crash, there is probably not enough RAM on your device. 


### Input
  -  files listed under Folder Structure
  -  config.R
  -  excel output file from the inhouse DIMS pipeline with Z-scores

### Output
  -  violin plots overview pdf file (all patients with a Z-score higher than 5 are annotated)
  -  violin plots per patient, stofgroups pdf file
  -  output algorithm excel file (with probability scores per metabolic disease)
  -  input csv file for shiny app 
  -  log file

