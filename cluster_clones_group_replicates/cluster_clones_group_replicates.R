###############################
##      0. Load Libraries    ##
###############################
# Load all required libraries for clustering, heatmap visualization, file I/O,
# and data manipulation.
library(cluster)
library(factoextra)
library(ggplot2)
library(xlsx)
library(rlist)
library(stringr)
library(gplots)
library(grDevices)
library(RColorBrewer)
library(pheatmap)
library(grid)
library(tidyverse)
library(dplyr)
library(icesTAF)

###############################
##   0.1 Define Dependencies ##
###############################
# Define the base directories and file paths used in the analysis.
baseDir            <- "~/Repositories/Frequency-dependentSelection_paper" # Base repository directory
dataDir            <- file.path(baseDir, "/cluster_clones_group_replicates/") # Directory containing data and code files
# TableS1Dir contains the single cell genome files as named_mat.csv.gz for each experiment. Named for Salehi, Shah (2021) paper from which the data originates
tableS1Dir         <- file.path(dataDir, "TableS1") # This needs to be added manually by the user as its too large to host on Git
metaDataPath       <- file.path(dataDir, "41586_2021_3648_MOESM3_ESM.xlsx") # Path to meta-data Excel file
replicatePathsOutputDir          <- file.path(dataDir, "paths")        # Directory for saving final output paths
shahPathsFile      <- file.path(dataDir, "shahExperimentPaths.xlsx") # Excel file containing experiment paths for reconstruction
lociFile           <- "hg38_cytoBand.txt"  # Cytoband file with genomic coordinates
origins            <- list('TNBC-SA609', 'TNBC-SA1035', 'TNBC-SA535', 'p53_KO1', 'p53_KO2') # List of experimental origins
shahCloneColorsDir <- "/color_arrays/"   # Directory with pre-defined color palettes for clone annotation (matches clone colors to other Figs)
clusterHeatmapsDir <- file.path(baseDir, "/shahClusterHeatmaps/") # Directory to store heatmap outputs

# Optional: provide a replicate file (with columns: origin, replicateID).
# If replicateFile is not desired or does not exist, set to NULL so the script computes replicate structure.
replicateFile <- file.path(dataDir, "replicateStructure.csv")

###############################
##      1. Copy-Number Caller  ##
###############################
# This section calls functions to convert segment-level copy number matrices
# into chromosome arm-level copy number matrices.
# The output 'arm_level_cn.rds' is later used by the clustering routine.

# ---- 1A. Define functions ----

# Function: get_loci_info()
# Description: Reads in a cytoband file for a human genome assembly and returns a data frame 
#              containing the start and end positions for each chromosome arm (p and q arms).
get_loci_info <- function(info_path){
  bands <- read.table(info_path, sep="\t")
  colnames(bands) <- c("chrom","start","end","band","idk")
  bands$arm <- substr(bands$band,1,1)         # Determine arm by taking first character of band name
  bands <- bands[bands$arm %in% c("p","q"),]  # Filter for p and q arms only
  
  bands <- split(bands, f=interaction(bands$chrom, bands$arm))
  
  df <- do.call(rbind, lapply(bands, function(bi){
    chrom = bi$chrom[1]
    chrom <- substr(chrom, 4, nchar(chrom))  # Remove any prefix (e.g., "chr") if present
    data.frame(chrom = chrom,
               arm = bi$arm[1],
               start = min(bi$start),
               end = max(bi$end))
  }))
  
  return(df)
}

# Function: Mode()
# Description: Returns the most frequent (modal) value from a numeric vector.
Mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}

# Set working directory to the data folder so that file operations occur in the correct location.
setwd(dataDir)

# Obtain genomic arm information from the cytoband file.
arms <- get_loci_info(lociFile)

# Get the list of file names from the TableS1 directory, excluding a README file if present.
fileNames <- list.files(tableS1Dir)
fileNames <- setdiff(fileNames, "README")

prevRownames <- NULL

# Loop through each experimental file directory in TableS1.
for(name in fileNames){
  
  # Change directory to the current experiment's folder.
  setwd(file.path(tableS1Dir, name))
  
  # Read in the segment-level copy number matrix (with row names as 'bin_name').
  x <- read.csv("named_mat.csv.gz", row.names = "bin_name")
  
  currentRownames <- rownames(x)
  
  # Check if the row names are numeric only (i.e., without underscores) and match the previous iteration.
  # If so, combine them with the previous row names to build a complete genomic coordinate.
  allNumeric <- all(grepl('^[0-9]+$', currentRownames))
  sameLength <- !is.null(prevRownames) && (length(prevRownames) == length(currentRownames))
  
  if(allNumeric && sameLength) {
    # Combine the previous row names with the current numeric ones (e.g., "1_100000001" + "100500000").
    combinedNames <- paste(prevRownames, currentRownames, sep = "_")
    rownames(x)   <- combinedNames
    
    # Update prevRownames with the combined names to use in subsequent iterations.
    prevRownames  <- combinedNames
    
  } else {
    # If row names are already in the correct format or lengths do not match,
    # update prevRownames with the current row names.
    prevRownames <- currentRownames
  }
  
  # Add chromosome arm information to each row.
  # Split the row names using "_" as a delimiter to separate chrom and position.
  df <- do.call(rbind, lapply(rownames(x), function(xi) unlist(strsplit(xi, split="_"))))
  df <- data.frame(df[,1:2])
  colnames(df) <- c("chrom","chrompos")
  df$arm <- sapply(1:nrow(df), function(i){
    chrom <- df$chrom[i]
    chrompos <- df$chrompos[i]
    # Determine if the position falls in the q arm by comparing to the start position from the arms data.
    is_q <- arms[arms$chrom == chrom & arms$arm == "q", "start"] < chrompos
    is_p <- !is_q
    c("q","p")[c(is_q, is_p)]
  })
  # Combine the chromosome and arm information with the copy number data.
  x <- cbind(df[, c("chrom", "arm")], x)
  
  # Split the matrix by unique chromosome-arm combinations.
  xcolnames <- colnames(x)
  x <- split(x, f=interaction(x$chrom, x$arm))
  lx <- sapply(x, nrow)
  x <- x[lx > 0]
  
  ## For each chromosome arm, assign the arm-level copy number as the mode of all its segment-level copy numbers.
  x <- data.frame(do.call(rbind, lapply(x,  function(xi){
    y <- c(xi[1, 1:2], as.numeric(apply(xi[, 3:ncol(xi)], 2, Mode)))
    y <- as.character(y)
    
  })))
  
  colnames(x) <- xcolnames
  # Save the arm-level copy number matrix for the current experiment.
  saveRDS(x, "arm_level_cn.rds")
  
}

##################################
##    2. Clustering and Heatmap ##
##################################
# This section clusters single-cell copy number data (arm-level) using k-means,
# then generates clone-frequency matrices and corresponding heatmaps.
# In the manuscript:
#   - 'times' or 'timepoint' corresponds to passage indices (t)
#   - 'km$cluster' corresponds to clone assignments (κᵢ)
#   - 'df_clust' summarizes per-time clone cell counts, which are then converted
#     into the frequency matrix Xᵣ used in the paper’s notation.

# Set working directory to dataDir and load additional libraries (repeating for clarity).
setwd(dataDir)
library(cluster)
library(factoextra)
library(ggplot2)
library(xlsx)
library(rlist)
library(stringr)
library(gplots)
library(grDevices)
library(RColorBrewer)
library(pheatmap)
library(grid)

# Read in meta-data from the Excel file.
meta <- read.xlsx("41586_2021_3648_MOESM3_ESM.xlsx", sheetIndex = 1)
meta$times <- sub("X", "", meta$timepoint)  # Remove prefix 'X' from timepoint values
meta$times <- as.numeric(meta$times)
meta$timepoint <- paste0(meta$timepoint, "_", meta$label)  # Combine timepoint and label for unique identification
meta$timepoint <- gsub("\\s+", " ", meta$timepoint)  # Clean extra whitespace if any
meta$origin <- sapply(str_split(meta$label, " ", n = 2), `[`, 1)  # Extract origin from label
meta$library_id <- trimws(as.character(meta$library_id))  # Ensure library IDs are clean strings

# Define a function to compute the average silhouette score for a given number of clusters.
silhouette_score <- function(k){
  km <- kmeans(x, centers = k, nstart = 25)
  ss <- silhouette(km$cluster, dist(x))
  mean(ss[, 3])
}
# Define a range of cluster numbers to test.
k <- 2:10

# Loop through each experimental origin to perform clustering and generate heatmaps.
for(org in unlist(origins)) {
  
  # Print progress for the current origin.
  print(paste("Clustering cells from:", org))
  
  # Identify unique dataset names (file folders) associated with the current origin.
  fileNames <- unique(meta[meta$origin == org, "datasetname"])
  losers <- c("SA004")  # Exclude any datasets that should not be processed.
  fileNames <- setdiff(fileNames, losers)
  fileNames <- na.omit(fileNames)
  
  arms <- list()   # List to store arm-level copy number matrices for each dataset.
  matches <- list()# List to store matching library IDs and timepoint info.
  valid_files <- logical(length(fileNames))  # Track which files have valid data for clustering
  
  # Process each file for the current origin.
  for (i in seq_along(fileNames)) {
    name <- fileNames[i]
    setwd(file.path(tableS1Dir, name))
    
    # Read in the pre-computed arm-level copy number matrix.
    x <- readRDS("arm_level_cn.rds")
    nx <- colnames(x)[-c(1, 2)]  # Extract cell/library names (ignoring chrom/arm columns)
    lookup <- sapply(strsplit(nx, ".", fixed = TRUE), "[[", 2)  # Parse to get matching IDs
    matched_idx <- match(lookup, meta$library_id)
    
    x <- x[, -c(1, 2)]  # Remove chromosome and arm columns; retain only copy number data
    x <- apply(x, 1, as.numeric)  # Convert data to numeric type
    arms[[name]] <- x
    matches[[name]] <- list(nx = nx, times = meta$timepoint[matched_idx])
    valid_files[i] <- TRUE
  }
  
  # If no valid data was found for the current origin, skip to the next origin.
  if (length(arms) == 0) {
    print(paste("No valid data for clustering in origin:", org))
    next
  }
  
  # Combine all arm-level matrices into a single matrix for clustering.
  arms <- list.rbind(arms)
  x <- arms
  
  # Use silhouette scores to determine the optimal number of clusters.
  avg_sil <- sapply(k, silhouette_score)
  best_k <- k[which.max(avg_sil)]
  
  # Perform k-means clustering using the optimal number of clusters.
  km <- kmeans(x, centers = best_k, nstart = 25)
  
  # Initialize lists to store timepoint information and heatmap data.
  tees <- list()
  enex <- list()
  heats <- list()
  
  # Loop again through each file to assign timepoint information to each cell.
  for(i in seq_along(fileNames)) {
    
    if(!valid_files[i]) next  # Skip any files that were marked as invalid
    
    name <- fileNames[i]
    setwd(file.path(tableS1Dir, name))
    
    # Read in the arm-level matrix and filter out sex chromosomes.
    x <- readRDS("arm_level_cn.rds")
    x <- x[!x$chrom %in% c("X", "Y"), ]
    x$chrom <- as.numeric(x$chrom)
    x <- x[order(x$chrom, x$arm), ]
    nx <- colnames(x)[-c(1:2)]
    lookup <- sapply(strsplit(nx, ".", fixed = TRUE), "[[", 2)
    
    # Debug: Print warning if any library IDs are missing in meta-data.
    missing_values <- setdiff(lookup, meta$library_id)
    if(length(missing_values) > 0) {
      print(paste("Mismatch detected. Missing lookup values in meta$library_id for file", name, ":", paste(missing_values, collapse = ", ")))
    }
    
    times <- meta$timepoint[match(lookup, meta$library_id)]
    
    # Prepare numeric copy number data for heatmap construction.
    x2 <- x[, -c(1, 2)]
    x2 <- apply(x2, 1, as.numeric)
    
    heats[[name]] <- x2
    tees[[name]] <- times
    enex[[name]] <- nx
  }
  
  # Combine timepoint and library names from all files.
  times <- unlist(tees)
  nx <- unlist(enex)
  
  # Create a data frame linking each cell (by name) with its cluster assignment and timepoint.
  df_clust <- data.frame(name = nx, cluster = km$cluster, times = times)
  
  # Aggregate clone counts per timepoint.
  df_clust$n <- 1
  df_clust <- aggregate(list(n = df_clust$n), by = list(cluster = df_clust$cluster, times = df_clust$times), sum)
  
  # Identify clusters that do not appear in every time point.
  allTimes <- unique(df_clust$times)
  unevenClusters <- rownames(plyr::count(df_clust$cluster)[plyr::count(df_clust$cluster)$freq != length(allTimes), ])
  
  # For clusters missing at a particular timepoint, add a dummy count (1 cell) so that every cluster is represented.
  for(clust in unevenClusters) {
    missedTimes <- setdiff(allTimes, df_clust[df_clust$cluster == clust, ]$times)
    for(t in missedTimes) {
      newRow <- data.frame(cluster = clust, times = t, n = 1)
      df_clust <- rbind(df_clust, newRow)
    }
  }
  
  # Convert the aggregated counts into a frequency matrix as expected by downstream MATLAB analysis.
  X_r <- as.data.frame(sapply(unique(df_clust$times), function(x) df_clust[df_clust$times == x, ]$n / sum(df_clust[df_clust$times == x, ]$n)))
  rownames(X_r) <- letters[1:nrow(X_r)]
  col <- rownames(X_r)
  
  # Save the clone frequency matrix to CSV.
  if(org == "p53-/-"){
    org <- name
  }
  setwd(dataDir)
  write.csv(X_r, file = paste0(org, "_cloneFreqs.csv"))
  
  # Create a combined heatmap object from all arm-level data.
  combined_heats <- do.call(rbind, heats)
  
  # If duplicate cell names exist, ensure they are made unique.
  if(anyDuplicated(nx)) {
    nx <- make.unique(nx)
  }
  
  # Set the row names of the heatmap data matrix.
  rownames(combined_heats) <- nx
  
  # Create a row annotation data frame for clusters.
  row_annotation <- data.frame(Cluster = factor(km$cluster))
  rownames(row_annotation) <- nx
  
  # Sort the heatmap data by cluster assignment.
  sorted_indices <- order(km$cluster)
  combined_heats <- combined_heats[sorted_indices, ]
  
  # Determine the number of clusters and create an appropriate color palette.
  num_clusters <- length(unique(km$cluster))
  colors <- brewer.pal(max(3, num_clusters), "Set3")
  ann_colors <- list(Cluster = setNames(colors[1:num_clusters], levels(row_annotation$Cluster)))
  
  # Read in a pre-defined clone colors CSV to match colors with other figures.
  csvFileName <- paste0(shahCloneColorsDir, 'clone_colors_', org, '.csv')
  cloneColorTable <- read.csv(csvFileName)
  
  # Extract RGB values and convert them to hexadecimal color codes.
  CloneID <- as.character(cloneColorTable$CloneID)
  R <- as.numeric(cloneColorTable$R)
  G <- as.numeric(cloneColorTable$G)
  B <- as.numeric(cloneColorTable$B)
  clone_colors_hex <- rgb(R, G, B)
  
  # Create a named vector of clone colors.
  clone_colors <- setNames(clone_colors_hex, CloneID)
  
  # If there are more clusters than available colors, generate additional colors.
  num_additional_colors <- best_k - length(clone_colors)
  if (num_additional_colors > 0) {
    additional_colors <- colorRampPalette(brewer.pal(8, "Dark2"))(num_additional_colors)
    clone_colors <- c(clone_colors, additional_colors)
  }
  ann_colors <- list(Cluster = setNames(clone_colors[1:num_clusters], unique(row_annotation$Cluster)))
  
  # Create a second annotation: assign time information.
  # The time string is split into two parts: the passage number and the sample ID.
  times_split <- str_split_fixed(times, "_", 2)
  times_before_underscore <- times_split[, 1]
  sampleID <- times_split[, 2]
  row_annotation$Time <- times_before_underscore[sorted_indices]
  num_times <- length(unique(times_before_underscore))
  gray_colors <- gray.colors(num_times, start = 0.9, end = 0.1)
  ann_colors$Time <- setNames(gray_colors, levels(factor(times_before_underscore)))
  
  # Add sample ID annotation and assign colors from a combined palette.
  row_annotation$Sample <- sampleID[sorted_indices]
  num_samples <- length(unique(sampleID))
  combined_palette <- c(brewer.pal(9, "Set1"), brewer.pal(8, "Set2"), brewer.pal(12, "Set3"))
  if (length(combined_palette) < num_samples) {
    combined_palette <- colorRampPalette(combined_palette)(num_samples)
  }
  sample_colors <- combined_palette[1:num_samples]
  ann_colors$Sample <- setNames(sample_colors, levels(factor(sampleID)))
  
  # Create a color palette for the heatmap itself based on the range of unique values.
  num_unique_values <- length(unique(combined_heats))
  spectral_palette <- colorRampPalette(brewer.pal(11, "Spectral"))(num_unique_values)
  
  # Generate a PDF file with the heatmap visualization.
  pdf(file = file.path(clusterHeatmapsDir, paste0("heatmap_", org, ".pdf")), width = 10, height = 8)
  
  pheatmap(combined_heats,
           annotation_row = row_annotation,
           annotation_colors = ann_colors,
           show_rownames = FALSE,
           show_colnames = TRUE,
           cluster_rows = FALSE,
           cluster_cols = TRUE,
           color = spectral_palette,
           main = paste("Copy Number Heatmap for", org))
  
  dev.off() # Close the PDF device
  
}

##############################################
##  3. Path Reconstruction
##############################################
# This section reconstructs evolutionary paths (cTP) from the clone frequency matrices.
# It reads in meta-data and clone frequency files, then for each origin builds a series
# of submatrices representing the evolution of clone frequencies. The output is saved
# in subdirectories under the replicatePathsOutputDir.

# Define the cTP function: clusters to paths.
# --------------------------------------------------------
# metaDataPath: Path to the meta-data Excel file from Shah et al. (2021)
# replicatePathsOutputDir: Directory under which replicate result subfolders are created
# replicateFile (optional): CSV file with columns 'origin' and 'replicateID'. If provided,
#                           it is used to define the replicate structure. If not, the script
#                           computes replicate information based on the data.
cTP <- function(metaDataPath, replicatePathsOutputDir, replicateFile = NULL){
  
  # Read in meta-data from the Shah et al. Excel file.
  meta <- read.xlsx(metaDataPath, sheetIndex = 1)
  meta$Passage <- sub("X", "", meta$timepoint)
  meta$Passage <- as.numeric(meta$Passage)
  meta$timepoint <- paste0(meta$timepoint, "_", meta$label)
  meta$origin <- sapply(str_split(meta$label, " ", n = 2), `[`, 1)
  meta$library_id <- trimws(as.character(meta$library_id))
  
  out <- list()
  
  # Loop over each experimental origin.
  for(org in unlist(origins)){
    message("Building cTP paths for origin: ", org)
    
    # Define the path to the clone frequency file generated in Section 2.
    freqFile <- file.path(dataDir, paste0(org, "_cloneFreqs.csv"))
    if(!file.exists(freqFile)) {
      message(" -- No cloneFreqs file found for ", org, ", skipping.")
      next
    }
    clusters <- read.csv(freqFile, row.names = 1)
    
    # Read the experiment paths from the Shah experiment paths Excel file.
    if(!file.exists(shahPathsFile)) {
      message(" -- No experiment path file found: ", shahPathsFile)
      next
    }
    
    paths <- tryCatch({
      read.xlsx(shahPathsFile, sheetIndex = org, header = FALSE)
    }, error = function(e) {
      message("Sheet not found for origin: ", org)
      return(NULL)
    })
    
    if(is.null(paths)) next
    
    # Determine the replicate structure based on the provided replicate file, if available.
    if(!is.null(replicateFile) && file.exists(replicateFile)){
      repStruct <- read.csv(replicateFile)
      repStruct <- repStruct[repStruct$origin == org, ]
      if(nrow(repStruct) != nrow(paths)){
        message("WARNING: Number of rows in replicate file (", nrow(repStruct),
                ") does not match number of rows in paths (", nrow(paths), ") for origin ", org)
      }
      replicateIDs <- repStruct$replicateID
      extraRep <- rep(NA, nrow(paths))  # No rxFrac data when using an external replicate file.
    } else {
      # If no external replicate file is provided, compute replicate info from the paths.
      repDecider <- c()
      fracker    <- c()
      
      # For each row in the paths sheet, calculate the fraction of "rx" (treated) cells.
      for(n in seq_len(nrow(paths))){
        entCount <- ncol(paths) - sum(is.na(paths[n,]))
        holidayCount <- length(unique(c(
          grep(".Holiday.Line", paths[n,]),
          grep(".Holiday", paths[n,])
        )))
        tmp <- unique(c(
          grep(".UnRx.Line", paths[n,]),
          grep(".Mixture.", paths[n,]),
          grep(".UnRx", paths[n,])
        ))
        unRxCount <- length(tmp)
        rxCount   <- entCount - unRxCount - holidayCount
        rxFrac    <- rxCount / entCount
        
        repDecider[n] <- rxFrac
        fracker[n]    <- rxFrac
      }
      
      # Map each unique rxFrac to a replicate ID.
      repDecider <- match(repDecider, unique(repDecider))
      replicateIDs <- repDecider
      extraRep <- fracker
    }
    
    newsies <- list()
    
    # For each path (row) in the experiment paths sheet, build a submatrix of clone frequencies.
    for(i in seq_len(nrow(paths))){
      newCluster <- matrix(NA, nrow(clusters), 1)
      rownames(newCluster) <- rownames(clusters)
      
      # Loop through the columns of the current path row to extract the corresponding data.
      for(k in seq_len(ncol(paths))){
        if(!is.na(paths[i, k])){
          colName <- as.character(paths[i, k])
          if(!colName %in% colnames(clusters)){
            message("DEBUG: In origin ", org, " row ", i, ", col ", k, ": colName '", colName, "' not found in clusters.")
            next
          }
          nc <- as.data.frame(clusters[, colName], row.names = rownames(clusters))
          names(nc) <- colName
          newCluster <- cbind(newCluster, nc)
        }
      }
      
      # Remove columns that are completely NA.
      newCluster <- newCluster[, colSums(is.na(newCluster)) < nrow(newCluster), drop = FALSE]
      
      # Filter out clones that never exceed 10% frequency (if desired).
      idx2 <- sapply(rownames(newCluster), function(x) any(newCluster[x,] > 0.10))
      newCluster <- newCluster[idx2, , drop = FALSE]
      
      newCluster <- as.data.frame(newCluster)
      newCluster$replicateID <- replicateIDs[i]
      newCluster$rxFrac      <- extraRep[i]  # Will be NA if using an external replicate file
      
      newsies[[i]] <- newCluster
      
      # Write the clone frequency submatrix to a CSV file.
      outDirPaths <- file.path(replicatePathsOutputDir, org, as.character(newCluster$replicateID[1]))
      mkdir(outDirPaths)  # Uses icesTAF's mkdir function (alternatively, use dir.create())
      
      outFile <- file.path(outDirPaths, paste0(org, "_", i, "_", newCluster$replicateID[1], "_cloneFreqs.csv"))
      write.csv(newCluster, file = outFile)
    }
    
    out[[org]] <- newsies
  }
  
  return(out)
}

# ---- 3B. Run the cTP function ----

# Ensure that the output directory for replicate paths exists.
if(!dir.exists(replicatePathsOutputDir)) dir.create(replicatePathsOutputDir, recursive = TRUE)

# Execute the cTP function, passing the optional replicate file if available.
results <- cTP(metaDataPath, replicatePathsOutputDir, replicateFile = replicateFile)

# Completion message.
message("All steps completed.")
