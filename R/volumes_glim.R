mni.read.volume.index <- function(lobe=TRUE) {
  index <- NA
  if(lobe == TRUE) {
    index <- "/data/nih/nih1/quarantine_NIH_D/model_data/jacob/seg/jacob_atlas_brain_fine_lobes.guide"
  }
  else {
    index <- "/data/nih/nih1/quarantine_NIH_D/model_data/jacob/seg/jacob_atlas_brain_fine.guide"
  }
  
  index.table.tmp <- read.table(index, fill=TRUE)
  empty <- index.table.tmp$V3 == ""
  index.table.tmp$V3[empty] <- NA
  index.table <- data.frame(row.names=index.table.tmp$V1,
                            label=gsub(" ", ".", paste(index.table.tmp$V3, index.table.tmp$V2, sep=".")))
  
  return(index.table)
}

mni.read.volume.data <- function(glim.matrix, lobe=TRUE) {
  # assume that the filename is in GLIM format
  index.table <- mni.read.volume.index(lobe=lobe)

  r <- rownames(index.table)
 
  number.subjects <- length(glim.matrix[,1])
  number.entries <- nrow(index.table)

  volumes.table <- matrix(NA, nrow=number.entries, ncol=number.subjects)
  rownames(volumes.table) <- index.table$label
  for (i in 1:number.subjects) {
    tmp <- as.data.frame(read.table(as.character(glim.matrix[i, 1])))
    tmp2 <- matrix(NA, nrow=number.entries, ncol=1)
    rownames(tmp2) <- index.table$label
    for (j in 1:nrow(tmp)) {
      tmp2[as.character(index.table[as.character(tmp$V1[j]),]), ] <- tmp$V2[j]
    }
    volumes.table[,i] <- tmp2
  }
    return(volumes.table)
 }
