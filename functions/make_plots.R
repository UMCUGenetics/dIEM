make_plots <- function(metab.list, stoftest, pt, zscore_cutoff, xaxis_cutoff, ThisProbScore,
                       ratios_cutoff, ptcount) {
  stoftest <- as.character(stoftest)
  if (ncol(metab.list) > 2) {
    # metab list are the alarm values, so reduce the data frame to 2 columns and save metab.list4
    metab.list4 <- metab.list
    metab.list <- metab.list[c(1,2)]
  }
  # extract the top 20 highest and top 10 lowest scoring metabolites
  if (startsWith(stoftest, "top") & ptcount > 1) {
    if (endsWith(stoftest, "hoogst")){
      topX <- unique(summed[pt]) %>% slice_max(unique(summed[pt]),n = 20)
      topX <- inner_join(topX, summed[,c(1,2,(ptcount+1))], by = pt)
    } else {
      topX <- unique(summed[pt]) %>% slice_min(unique(summed[pt]),n = 10)
      topX <- inner_join(topX, summed[,c(1,2,(ptcount+1))], by = pt)
    }
    count = 0
    listrm <- c()
    b <- 100000.00000
    for (rows in c(1:nrow(topX))) {
      a <- topX[rows,1]
      if (a==b) { count = count + 1 }
      else { count = 0 }
      if (count >= 6) { listrm <- c(listrm, rows) }
      b <- a 
    }
    if (length(listrm) > 0) { 
      metab.list <- topX[-(listrm),-1]
      isos <- other_isobaric(metab.list,"infile",FALSE)
      metab.list <- rbind(metab.list, isos)
    } else {
      metab.list <- topX[-1]
    }
    # replace the metab.list dataframe with new values
    #metab.list <- topX[-1]
    
  }
  count <- 0
  for (metab in metab.list$HMDB_name){
    count <- count + 1 
    metab.list$HMDB_name[count] <- gsub("(.{45})", "\\1...;", metab, perl = TRUE)
  }
  # Filter summed on the metabolites of interest (moi)
  joined <- inner_join(metab.list, summed[-2], by = "HMDB_code")
  moi <- joined[,-2]
  j <- joined[,-1]
  jm <- reshape2::melt(j, id.vars = "HMDB_name")
  #jmo <- jm[ rev(order(match(jm$HMDB_name, joined$HMDB_name))), ]
  jma <- aggregate(HMDB_name ~ value+variable, jm, paste0, collapse = "_")
  jmas <- jma %>% separate(HMDB_name, into = c("HMDB_name", "isobar"), sep="_",extra = "merge", fill = "right")
  jmaso <- jmas[ rev(order(match(jmas$HMDB_name, joined$HMDB_name))), ]
  jmao <- jmaso %>% unite(HMDB_name,c("HMDB_name","isobar"),sep="_", na.rm = TRUE)
  #enters <- max(lengths(regmatches(jmao$HMDB_name, gregexpr(";|_", jma$HMDB_name))))
  jmao$HMDB_name <- gsub(';|_','\n',jmao$HMDB_name)
  jmao$HMDB_name <- factor(jmao$HMDB_name, levels=unique(jmao$HMDB_name))
  if (check.lists & (pt=="alle")) {
    write.xlsx(jmaso, paste0(output_dir,"/", stoftest, "_check2.xlsx"))
  }
  
  # Split the metabolite lists from isobaric compounds
  
  isos <- jmao %>% separate(HMDB_name, into = c("HMDB_name", "isobar"), sep="zzz",extra = "merge", fill = "right")
  b <- "yellow"
    count <- 0
    for (rownr in c(nrow(isos):1)) {
      a <- isos[rownr,4]
      if (!identical(a,b) & !is.na(a)) {
        count <- count + 1
      }
      isos[rownr,4] <- gsub("^", paste0(count,": "), isos[rownr,4], perl = TRUE)
      isos[rownr,3] <- gsub("\nqqq", paste0('^',count), isos[rownr,3], perl = TRUE) 
      b <- a
    }
    
    moi_m <- isos[,c(1:3)]
    moi_m$HMDB_name <- factor(moi_m$HMDB_name, levels=unique(moi_m$HMDB_name))
    footn <- isos[,4]
    footnu <- unique(footn[!is.na(footn)])
    footnus <- gsub("...\n", "", gsub("\nqqqzzz", "; ", footnu), perl = FALSE)
    for (i in c(1:length(footnus))){
      footnus[i] <- gsub("(.{120})", "\\1\n", gsub(":", ":\t", footnus[i]), perl = TRUE)
    }
    footnus <- footnus[length(footnus):1]
    if (pt=="alle") {
      isobarics_txt <<- c(isobarics_txt,"\n_____________________________________",stoftest,"\n",footnus)
    }
    
    namelist <- as.vector(unique(moi_m$HMDB_name))
    nrows <- length(namelist)
    if (ThisProbScore!=0) {
      stoftest.chunks <- list(moi_m)
    }
    if (split) {
      subSetSizes <- 15:30
      remainders <- nrows %% subSetSizes 
      minIndexes <- which(remainders == min(remainders))
      chunkSize <- max(subSetSizes[minIndexes])
      numberlists = nrows / chunkSize 
      moi_m.chunks <- list()
      b = 0
      c = 0
      stoftest.chunks <- c()
      for (index in c(1:numberlists)) { 
        b = c + 1
        stoftest.chunks <- c(stoftest.chunks,paste0(stoftest,"_",index))
        if (index==1){ a <- chunkSize+min(remainders)} else {a <- chunkSize}
        c = b + a - 1
        selectrow <- namelist[b:c]
        moi_m_cut <- moi_m[moi_m$HMDB_name %in% as.vector(namelist[b:c]),]
        moi_m.chunks[[index]] <- moi_m_cut
      }
      # reverse order of lists. 
      moi_m.chunks <- moi_m.chunks[length(moi_m.chunks):1]  
    }
    for (moi_m in moi_m.chunks) {
      i_tot <- length(unique(moi_m$HMDB_name))
      enters <- max(lengths(regmatches(unique(moi_m$HMDB_name), gregexpr("\n", unique(moi_m$HMDB_name)))))
      # adjust size for the font of y-axis labels and plotted dot sizes
      if (enters > 3) {
        fontsize1 <- -0.25*enters +1.75
      } else {
        fontsize1 <- 1
      }
      if (i_tot>15) {
        fontsize2 <- -0.02*i_tot + 1.30
        circlesize <- -0.013*i_tot + 1.26
      } else {
        fontsize2 <- 1
        circlesize <- 1
      }
      if (enters*i_tot>45){
        fontsize3 <- -0.006*(enters*i_tot) + 1.25
      } else {
        fontsize3 <- 1
      }
      #print(c(fontsize1, fontsize2, fontsize3, (enters*i_tot)))
      fontsize <- min(c(fontsize1, fontsize2, fontsize3))
      if (fontsize < 0.2) { fontsize <- 0.2 }
      if (circlesize < 0.3) { circlesize <- 0.3 }
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
      moi_m_max20$value <- as.numeric(lapply(moi_m_max20$value, function(x) ifelse(x < ratios_cutoff, as.numeric(ratios_cutoff), x)))
      
      # Get values that overlap with highZ and max20: i.e. 5 < z < 20
      group_highZ_max20 <- moi_m_max20 %>%
        group_by(value) %>% 
        filter(value > zscore_cutoff) %>%
        ungroup()
      
      if (stoftest=="0_alarmwaardes" & ThisProbScore==0 & ptcount > 1){
        # This will be front page with alarmvalues, no plots.
        # pt <- patient_list[3]
        # metab.list4 <- metab.list0[[1]]
        cat(pt)
        alarm <- inner_join(metab.list4, summed[,c("HMDB_code", pt)], by = "HMDB_code")
        colnames(alarm)[5] <- "zscore"
        #alarm$zscore <- as.numeric(round(alarm$zscore, 2))
        alarm$zscore[which(alarm$soort_grens=="onder")] <- as.numeric(lapply(alarm$zscore[which(alarm$soort_grens=="onder")], function(x) ifelse(x < alarm$alarmwaarde[which(alarm$zscore==x)], x, 0)))
        alarm$zscore[which(alarm$soort_grens=="boven")] <- as.numeric(lapply(alarm$zscore[which(alarm$soort_grens=="boven")], function(x) ifelse(x > alarm$alarmwaarde[which(alarm$zscore==x)], x, 0)))
      }
      
      if (!startsWith(pt,"all")){
        # patient one by one
        pt <- gsub("_IEM", "", pt)
        pt_colname <- pt
        # get values from patient
        pt_data <- moi_m[which(moi_m$variable==pt_colname),]
        pt_data_max20 <- moi_m_max20[which(moi_m_max20$variable==pt_colname),]
        pt_values <- pt_data$value # change this later, this gives warnings
        
        colors <- c("#22E4AC", "#00B0F0", "#504FFF","#A704FD","#F36265","#DA0641")
        #             green     blue      blue/purple purple    orange    red
        if (ThisProbScore==0 & !startsWith(stoftest,"top") & ptcount > 1){
          if (stoftest=="0_alarmwaardes"){
            plot.new()
            zeroes <- which(alarm$zscore!=0)
            if (length(zeroes)==0) {
              text(x=.5, y=.95, paste0("Dit zijn de alarmwaardes voor patient:\n\n",pt), font=1, cex=1, col="#F48024")
              text(x=.5, y=.85, paste0("Geen afwijkende waardes"), font=1, cex=1, col="#000000")
            } else {
              alarm <- alarm[zeroes,]
              find_cell <- function(table, row, col, name="core-fg"){
                l <- table$layout
                which(l$t==row & l$l==col & l$name==name)
              }
              alarm <- alarm[-1]
              alarm_table <- tableGrob(alarm)
              #color_this <- which(alarm$zscore!=0)
              #color_this <- ""
              high <- which(alarm$zscore!=0 & alarm$soort_grens=="boven")
              low <- which(alarm$zscore!=0 & alarm$soort_grens=="onder")
              # if (length(color_this)!=0){
              #   for (i in color_this){
              #     ind <- find_cell(alarm_table, i+1, 5, "core-bg")
              #     alarm_table$grobs[ind][[1]][["gp"]] <- gpar(fill="darkolivegreen1", col = "darkolivegreen4",fontface="bold", lwd=3)
              #   }
              #   
              # }
              if (length(high)!=0){
                for (i in high){
                  ind <- find_cell(alarm_table, i+1, 5, "core-bg")
                  alarm_table$grobs[ind][[1]][["gp"]] <- gpar(fill="#f2cdd7", col = colors[6],fontface="bold", lwd=3)
                }
              }
              if (length(low)!=0){
                for (i in low){
                  ind <- find_cell(alarm_table, i+1, 5, "core-bg")
                  alarm_table$grobs[ind][[1]][["gp"]] <- gpar(fill="#bcf6e6", col = colors[1],fontface="bold", lwd=3)
                }
              }
              text(x=.5, y=.95, paste0("Dit zijn de alarmwaardes voor patient:\n\n",pt), font=1, cex=1, col="#F48024")  # first 2 numbers are xy-coordinates within [0, 1]
              #text(x=.1, y=.8, paste0("Bovengrens:"), font=1, cex=0.5)  # first 2 numbers are xy-coordinates within [0, 1]
              grid.draw(alarm_table)
            }
            
            
          } else {
            # plot each stofgroup. NB: MPR changed shape from 21 to 22 (square)
            g <- ggplot(moi_m_max20, aes(x=value, y=HMDB_name))+
              theme(axis.text.y=element_text(size=rel(fontsize)), plot.caption = element_text(size=rel(fontsize)))+
              geom_violin(scale="width")+
              geom_point(data = pt_data_max20, aes(color=pt_data$value),size = 3.5*circlesize,shape=22, fill="white")+
              geom_jitter(data = group_highZ_max20, aes(color=group_highZ$value), size = 1.3*circlesize, position = position_dodge(1.5))+ #,colour = "#3592b7" 
              scale_fill_gradientn(colors = colors,values = NULL,space = "Lab",na.value = "grey50",guide = "colourbar",aesthetics = "colour")+
              labs(x = "Z-scores",y = "Metabolites",title = paste0("Results for patient ",pt), subtitle = stoftest, color = "z-score", caption = "Voor voetnoot zie 'isobarics.txt'")+
              geom_vline(xintercept = 2, col = "grey", lwd = 0.5,lty=2)+
              geom_vline(xintercept = -2, col = "grey", lwd = 0.5,lty=2)
            print(g)
          }
        }
        if (ThisProbScore==0 & startsWith(stoftest,"top") & ptcount > 1){
          # the top20 & top10 plots
          g <- ggplot(moi_m, aes(x=value, y=HMDB_name, color = value))+
            theme(axis.text.y=element_text(size=rel(fontsize)))+
            geom_violin(scale="width")+
            geom_point(data = pt_data, aes(color=pt_data$value),size = 3.5*circlesize,shape=21, fill="white")+
            geom_jitter(data = group_highZ, aes(color=group_highZ$value), size = 1.3*circlesize, position = position_dodge(1.5))+ #,colour = "#3592b7" 
            scale_fill_gradientn(colors = colors,values = NULL,space = "Lab",na.value = "grey50",guide = "colourbar",aesthetics = "colour")+
            labs(x = "Z-scores",y = "Metabolites",title = paste0("Results for patient ",pt), subtitle = stoftest, color = "z-score")+
            geom_vline(xintercept = 2, col = "grey", lwd = 0.5,lty=2)+
            geom_vline(xintercept = -2, col = "grey", lwd = 0.5,lty=2)
          #the plot without x-axis constraints
          print(g)
          
        }
        if (ThisProbScore > 0){
          # plot the metabolites of the top 5 IEMs
          g <- ggplot(moi_m_max20, aes(x=value, y=HMDB_name, color = value))+
            theme(axis.text.y=element_text(size=rel(fontsize)))+
            geom_violin(scale="width")+
            geom_point(data = pt_data_max20, aes(color=pt_data$value),size = 3.5*circlesize,shape=21, fill="white")+
            geom_jitter(data = group_highZ_max20, aes(color=group_highZ$value), size = 1.3*circlesize, position = position_dodge(1.5))+ #,colour = "#3592b7" 
            scale_fill_gradientn(colors = colors,values = NULL,space = "Lab",na.value = "grey50",guide = "colourbar",aesthetics = "colour")+
            labs(x = "Z-scores",y = "Metabolites",title = paste0("Algorithm results for patient ",pt), subtitle = paste0("Disease: ",stoftest,"\nProbability Score = ",format(round(ThisProbScore, 2), nsmall = 2)), color = "z-score")+
            geom_vline(xintercept = 2, col = "grey", lwd = 0.5,lty=2)+
            geom_vline(xintercept = -2, col = "grey", lwd = 0.5,lty=2)
          print(g)
        }
        
      }
      opencircle <- FALSE
      gradient <- FALSE
      # overview plots
      if (pt=="alle" & !startsWith(stoftest, "top") & !(stoftest=="0_alarmwaardes")){
        #overview plot with zscores of max 20 on the x-axis
        g <- ggplot(moi_m_max20, mapping = aes(x=value, y=HMDB_name))+
          theme(axis.text.y=element_text(size=rel(fontsize)))+
          geom_violin(scale="width")+
          {if(opencircle) geom_point(data = pt_data_max20, aes(color=pt_data$value),size = 3.5*circlesize,shape=21, fill="white")} +
          geom_jitter(data = group_highZ_max20,aes(color = group_highZ$variable), size = 2.5*circlesize, position = position_dodge(0.8))+ 
          {if (gradient) scale_fill_gradientn(colors = colors,values = NULL,space = "Lab",na.value = "grey50",guide = "colourbar",aesthetics = "colour")} +
          labs(x = "Z-scores",y = "Metabolites",title = "Overview plot", subtitle = stoftest, color = "patients")+
          geom_vline(xintercept = 2, col = "grey", lwd = 0.5,lty=2)+
          geom_vline(xintercept = -2, col = "grey", lwd = 0.5,lty=2)
        
        #delete_layers(g, "GeomText")
        print(g)
      }
      if (pt=="alle_volle_x-as" & !startsWith(stoftest, "top") & !(stoftest=="0_alarmwaardes")){
        #overview plot without x-axis constraints
        g <- ggplot(moi_m, mapping = aes(x=value, y=HMDB_name))+
          theme(axis.text.y=element_text(size=rel(fontsize)))+
          geom_violin(scale="width")+
          geom_jitter(data = group_highZ,aes(color = group_highZ$variable), size = 2.5*circlesize, position = position_dodge(0.8))+ 
          labs(x = "Z-scores",y = "Metabolites",title = "Overview plot", subtitle = stoftest, color = "patients")+
          geom_vline(xintercept = 2, col = "grey", lwd = 0.5,lty=2)+
          geom_vline(xintercept = -2, col = "grey", lwd = 0.5,lty=2)
        print(g)
      }
    }
    #return(isobarics_txt)
}