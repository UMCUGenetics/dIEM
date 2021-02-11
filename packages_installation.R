list.of.packages <- c("beepr", "dplyr", "reshape2", "openxlsx", "ggplot2", "gghighlight", "data.table", "tidyr")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

