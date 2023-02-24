other_isobaric <- function(Combined, metab.list, infile, check.lists) {
  # not sure what this function does, find metabolites with same mass? Isomeric, not isobaric
  nrpat <- length(grep("P",names(Combined)))/2
  nrcontr <- length(grep("P",names(Combined)))/2
  same.mass <- Combined[c(1,2,(nrpat+nrcontr+5))] # HMDB columns plus first Z-score column
  if (ncol(metab.list) > 2) {
    cat("alarmvalues")
  } else { 
    names(same.mass) <- c("HMDB_code","HMDB_name","zscore")
    metab.list.int <- left_join(metab.list, same.mass, by = "HMDB_code")
    # if (check.lists) {
    #   write.xlsx(metab.list.int, paste0(output_dir,"/check_", infile, ".xlsx"))
    # }
    isobaric.metab.list <- unique(left_join(metab.list.int, same.mass, by = "zscore", keep = FALSE)[ , c(5,6)])
    names(isobaric.metab.list) <- c("HMDB_code", "HMDB_name")
    isobaric.rest <- anti_join(isobaric.metab.list, metab.list, by='HMDB_code')
    isobaric.rest <- isobaric.rest[order(isobaric.rest$HMDB_name, na.last = NA), ]
    isobaric.rest$HMDB_name <- lapply(isobaric.rest$HMDB_name, function(x) { gsub("^", "qqqzzz", x, perl = TRUE) })
    isobaric.rest[, 2] <- sapply(isobaric.rest[, 2], as.character)
    return(isobaric.rest)
  }
}
