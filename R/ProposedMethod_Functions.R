#================================================================
# This function calculate the p-value of identify q genes from all k genes
# P(observed q-1 or more)
# (q: success-in-sample, m:success-in-bkgd, n:failure-in-bkgd, K:sample-size).
#================================================================
pValue <- function(q = NULL, m = NULL, n = NULL, k = NULL) {
  
  res = 1.0 - phyper(q-1, m, n, k, lower.tail = TRUE, log.p = FALSE)
  res
}

#================================================================
#' Adjust p-values
#================================================================
adjustpValues = function(results) {
  
  r <- results
  
  nR <- nrow(r)
  nC <- ncol(r)
  t <- as.vector(r)
  t <- p.adjust(t, method="fdr")
  r <- matrix(t, nrow = nR, ncol = nC)
  
  row.names(r) <- row.names(results)
  colnames(r) <- colnames(results)
  
  return(r)
}

#================================================================
#' Identify links among nodes from their expression data
#' @param data Expression data with rows being samples and columns being biological features
#' @param cause Range of cause
#' @param effect Range of effect
#' @param rootDir Root folder
#' @param f File name
#' @param usingpVal TRUE if using p-value
#' @param cutoff FDR cutoff
#' @return Edges from cause to effect
#================================================================
identifyEdges = function(data, cause, effect, rootDir, f = "temp", usingpVal = FALSE, cutoff = 0.05) {
  
  if(usingpVal == TRUE) {
    # Use Pearson to evaluate the relationship among cause and effect
    dataset <- paste(rootDir, "/Data/Output/dataset.csv", sep = "")
    write.csv(data[, c(cause, effect)], dataset, row.names = FALSE)
    results = Pearson_pValue(dataset, 1:length(cause), (length(cause)+1):(length(cause)+length(effect)))
    results <- t(results)
    
    # Get links which have p-value < cutoff
    results <- adjustpValues(results)
    results <- results < cutoff
    ind <- which(results, arr.ind=TRUE)
    edges <- ind
    edges[, 1] <- row.names(ind)
    edges[, 2] <- colnames(results)[ind[, 2]]
    
    # Remove dataset.csv file
    file.remove(paste(rootDir, "/Data/Output/dataset.csv", sep = ""))
  } else {
    # Use Pearson to evaluate the relationship among cause and effect
    dataset <- paste(rootDir, "/Data/Output/dataset.csv", sep = "")
    write.csv(data[, c(cause, effect)], dataset, row.names = FALSE)
    results = Pearson(dataset, 1:length(cause), (length(cause)+1):(length(cause)+length(effect)))
    results <- t(results)
    
    # Write file
    t <- paste(rootDir, "/Data/Output/PEDriver/Network/", f, ".csv", sep = "")
    write.csv(results, t, row.names = TRUE)
    
    # Get links which have the absolute coefficients more than
    # the average of all absolute coefficients
    a <- mean(abs(results))
    results <- abs(results) >= a
    ind <- which(results, arr.ind=TRUE)
    edges <- ind
    edges[, 1] <- row.names(ind)
    edges[, 2] <- colnames(results)[ind[, 2]]
    
    # Remove dataset.csv file
    file.remove(paste(rootDir, "/Data/Output/dataset.csv", sep = ""))
  }
  
  return(edges)
}

#================================================================
#' Build a network based on the expression data then filtering by existing datasets
#' @param interactions Interactions among mRNAs and TFs
#' @param nomR Number of mRNAs
#' @param noTF Number of transcription factors
#' @param data Expression data with rows being samples and 
#'    columns being biological features with the order mRNAs, TFs
#' @param rootDir Root folder
#' @return Edges of the network from cause to effect
#================================================================
buildNetworkForPersonalised = function(interactions, nomR, noTF, data, rootDir,usingpVal = FALSE, cutoff = 0.05){
  # Build network
  # TF => mRNA
  print("Build network...")
  print("TF => mRNA...")
  edges_cod_cod <- identifyEdges(data, (nomR+1):(nomR + noTF),
                                 1:nomR, rootDir,usingpVal, cutoff)
  # mRNA => mRNA
  print("mRNA => mRNA...")
  temp <- identifyEdges(data, 1:nomR, 1:nomR, rootDir,usingpVal, cutoff)
  temp[, 2] <- gsub("\\..*","",temp[, 2])
  edges_cod_cod <- rbind(edges_cod_cod, temp)
  # TF => TF
  print("TF => TF...")
  temp <- identifyEdges(data, (nomR+1):(nomR + noTF),
                        (nomR+1):(nomR + noTF), rootDir,usingpVal, cutoff)
  temp[, 2] <- gsub("\\..*","",temp[, 2])
  edges_cod_cod <- rbind(edges_cod_cod, temp)
  # Set colnames
  colnames(edges_cod_cod) <- c("cause", "effect")
  
  # Filter with existing datasets
  # TF => mRNA, mRNA => mRNA, TF => TF: Filter with interactions
  colnames(interactions) <- c("cause", "effect")
  interactions <- merge(interactions, edges_cod_cod)
  
  return(interactions)
}

#================================================================
#' Build a network based on the expression data then filtering by existing datasets
#' @param nomR Number of mRNAs
#' @param noTF Number of transcription factors
#' @param data Expression data with rows being samples and 
#'    columns being biological features with the order mRNAs, TFs
#' @param rootDir Root folder
#' @return Edges of the network from cause to effect
#================================================================
buildNetworkFromData = function(nomR, noTF, data, rootDir,usingpVal = FALSE, cutoff = 0.05){
  # Build network
  # TF => mRNA
  print("Build network...")
  print("TF => mRNA...")
  edges_cod_cod <- identifyEdges(data, (nomR+1):(nomR + noTF),
                                 1:nomR, rootDir,usingpVal, cutoff)
  # mRNA => mRNA
  print("mRNA => mRNA...")
  temp <- identifyEdges(data, 1:nomR, 1:nomR, rootDir,usingpVal, cutoff)
  temp[, 2] <- gsub("\\..*","",temp[, 2])
  edges_cod_cod <- rbind(edges_cod_cod, temp)
  # TF => TF
  print("TF => TF...")
  temp <- identifyEdges(data, (nomR+1):(nomR + noTF),
                        (nomR+1):(nomR + noTF), rootDir,usingpVal, cutoff)
  temp[, 2] <- gsub("\\..*","",temp[, 2])
  edges_cod_cod <- rbind(edges_cod_cod, temp)
  # Set colnames
  colnames(edges_cod_cod) <- c("cause", "effect")
  
  return(edges_cod_cod)
}

#================================================================
#' Build a network based on the expression data then filtering by existing datasets
#' @param nomR Number of mRNAs
#' @param noTF Number of transcription factors
#' @param data Expression data with rows being samples and 
#'    columns being biological features with the order mRNAs, TFs
#' @param rootDir Root folder
#' @return Edges of the network from cause to effect
#================================================================
getPearsonData = function(nomR, noTF, PE_data, rootDir){
  # Build network
  # TF => mRNA
  print("Calculate Pearson coeffience...")
  # res = matrix(nrow = nomR*noTF+nomR*nomR+noTF*noTF, ncol = 4)
  # nameOfGenes = colnames(PE_data)
  # k = 0
  
  print("TF => mRNA...")
  res = matrix(nrow = nomR*noTF, ncol = 4)
  nameOfGenes = colnames(PE_data)
  # Set colnames
  colnames(res) <- c("cause", "effect", "cor", "pvalue")
  k = 0
  for (i in (nomR+1):(nomR + noTF)) {
    for (j in 1:nomR) {
      k = k+1
      rt <- cor.test(PE_data[,i], PE_data[,j], method="pearson")
      res[k,] = c(nameOfGenes[i], nameOfGenes[j], rt$estimate, rt$p.value)
    }
  }
  write.csv(res, file = "TF_mRNA.csv",row.names = F)
  gc()
  
  # mRNA => mRNA
  print("mRNA => mRNA...")
  res = matrix(nrow = nomR*nomR-nomR, ncol = 2)
  #nameOfGenes = colnames(PE_data)
  # Set colnames
  colnames(res) <- c("cor", "pvalue")
  k = 0
  for (i in 1:(nomR-1)) {
    for (j in (i+1):nomR) {
      k = k+1
      rt <- cor.test(PE_data[,i], PE_data[,j], method="pearson")
      res[k,] = c(rt$estimate, rt$p.value)
      k = k+1
      res[k,] = c(rt$estimate, rt$p.value)
    }
  }
  write.csv(res, file = "mRNA_mRNA.csv",row.names = F)
  gc()

  # TF => TF
  print("TF => TF...")
  res = matrix(nrow = noTF*noTF-noTF, ncol = 2)
  #nameOfGenes = colnames(PE_data)
  # Set colnames
  colnames(res) <- c("cor", "pvalue")
  k = 0
  for (i in (nomR+1):(nomR + noTF-1)) {
    for (j in (i+1):(nomR + noTF)) {
      rt <- cor.test(PE_data[,i], PE_data[,j], method="pearson")
      k = k+1
      res[k,] = c(rt$estimate, rt$p.value)
      k = k+1
      res[k,] = c(rt$estimate, rt$p.value)
    }
  }
  write.csv(res, file = "TF_TF.csv",row.names = F)
  gc()
  
  return(1)
}

#================================================================
#' Analyse a network
#' @param nomR Number of mRNAs
#' @param noTF Number of transcription factors
#' @param network Edges of the network
#' @param data Expression data with rows being samples and 
#'    columns being biological features with the order mRNAs, TFs
#' @param fileName File to be written
#================================================================
analyseNetworkForPersonalised = function(nomR, noTF, network, data, fileName){
  # Prepare the data
  dat <- ""
  dat <- paste(dat, "Number of nodes:", nomR + noTF, "\n", sep = "\t")
  dat <- paste(dat, "Number of TFs:", noTF, "\n", sep = "\t")
  dat <- paste(dat, "Number of mRNAs:", nomR, "\n", sep = "\t")
  dat <- paste(dat, "\n", sep = "\t")
  mRs <- colnames(data)[1:nomR]
  TFs <- colnames(data)[(nomR+1):(nomR+noTF)]
  edge_type4 <- network[which(network$cause %in% TFs),]
  edge_type4 <- edge_type4[which(edge_type4$effect %in% TFs),]
  edge_type4 <- nrow(edge_type4)
  edge_type5 <- network[which(network$cause %in% TFs),]
  edge_type5 <- edge_type5[which(edge_type5$effect %in% mRs),]
  edge_type5 <- nrow(edge_type5)
  edge_type6 <- network[which(network$cause %in% mRs),]
  edge_type6 <- edge_type6[which(edge_type6$effect %in% mRs),]
  edge_type6 <- nrow(edge_type6)
  edge <- edge_type4 + edge_type5 + edge_type6
  dat <- paste(dat, "Number of edges:", edge, "\n", sep = "\t")
  dat <- paste(dat, "Number of TF-TF edges:", edge_type4, "\n", sep = "\t")
  dat <- paste(dat, "Number of TF-mRNA edges:", edge_type5, "\n", sep = "\t")
  dat <- paste(dat, "Number of mRNA-mRNA edges:", edge_type6, "\n", sep = "\t")
  
  # Write file
  writeLines(dat, fileName)
}

#================================================================
#' Analyse controllability of the network
#' @param result Result file of analysing controllability of the network
#' @param outFile File to be written
#================================================================
analyseControllability = function(result, outFile){
  # Read files
  result_output <- read.csv(result, header = FALSE)
  # Prepare the data
  items <- c("1: # of nodes having links", "2: # of edges", "3: average degree",
             "4: # of driver nodes", "5: fraction of driver nodes = Nd/N",
             "6: fraction of type-I critical nodes",
             "7: fraction of type-I redundant nodes", "8: fraction of type-I ordinary nodes",
             "9: fraction of type-II critical nodes",
             "10: fraction of type-II redundant nodes", "11: fraction of type-II ordinary nodes",
             "12: fraction of critical links",
             "13: fraction of redundant links", "14: fraction of ordinary links")
  dat <- ""
  for (i in 1:14) {
    dat <- paste(dat, items[i], result_output[1,i], "\n", sep = "\t")  
    if(i %in% c(3,5,8,11)) {
      dat <- paste(dat, "\n")
    }
  }
  
  # Write file
  writeLines(dat, outFile)
}

#================================================================
#' Convert a list to a string
#' @param aList A list
#' @return A string contains items of the list
#================================================================
getList = function(aList){
  len <- length(aList)
  l <- aList[1]
  if(len > 1) {
    for(i in 2:len) {
      l <- paste(l, aList[i], sep=",")
    }
  }
  
  return(l)
}


#================================================================
#' This function allows you to remove genes with low expression
#' @param d The data with rows being samples and columns being genes.
#' @param threshold The threshold to remove genes.
#' @return A matrix with rows being samples and columns being genes
#================================================================
removeLowExpGenes <- function(d, threshold) {
  
  # Rotate the matrix
  matrix <- t(d)
  
  # Remove all the rows that never have expression of threshold or higher
  r<- matrix[apply(matrix[,],1,function(X) length(X[X>=threshold])>0),]
  
  # Rotate the result
  r <- t(r)
  
  return(r)
}
