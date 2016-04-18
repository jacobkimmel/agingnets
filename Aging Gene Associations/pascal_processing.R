library(stringr)
library(ggplot2)

# Import a single Pascal pathway analysis output
# Plots p-values on a histogram and exports sig. enriched
# pathways as a list to enriched/
processPascal <- function(data_file){
  df <- read.table(data_file, stringsAsFactors = FALSE)
  names(df) <- df[1,]
  df <- df[-c(1),]
  df$chi2Pvalue <- as.numeric(df$chi2Pvalue)
  df$empPvalue <- as.numeric(df$empPvalue)
  # Sort based on each statistical mode
  df.chi2sort <- df[order(df[,2]),]
  df.empsort <- df[order(df[,3]),]
  
  # Plot histogram of chi-squared p-vals
  require(ggplot2)
  data_set_name <- str_sub(data_file, 6, -10)
  filename <- paste('plots/', data_set_name, '_', 'ChiSq_Freq_Dist.png', sep = '')
  png(filename, width = 7, height = 5, units = 'in', res = 300)
  print(
    ggplot(data=df, aes(df$chi2Pvalue)) + geom_histogram() + labs(x='Chi-squared p-value', y='Frequency', title='Pathway Analysis P-value distribution')
    
  )
  dev.off()
  
  # Determine which pathways are significantly enriched
  num_enriched <- length(which(df.chi2sort$chi2Pvalue<0.05))
  enriched <- df.chi2sort[1:num_enriched,1]
  filename <- paste('enriched/', data_set_name, '_', 'EnrichedPaths.txt', sep = '')
  write(enriched, file = filename, sep = ',')
}

# Runs processPascal in batch on all Pascal outputs in the data/ directory
batchprocess <- function(data_dir) {
  filenames <- list.files(data_dir, pattern = "*REACTOME--sum.txt", full.names = TRUE)
  res <- lapply(filenames, processPascal)
  return(enriched_pathways)
}

# Counts how many Pascal outputs identify a given pathway as enriched
# Exports results as a table in the enriched directory 
rich_paths <- function(enriched_dir){
  filenames <- list.files(enriched_dir, pattern = '*.txt', full.names = TRUE)
  for (file in filenames){
    
    # if the merged dataset doesn't exist, create it
    if (!exists("dataset")){
      dataset <- read.table(file, header=FALSE, sep="\t")
    }
    
    # if the merged dataset does exist, append to it
    if (exists("dataset")){
      temp_dataset <-read.table(file, header=FALSE, sep="\t")
      dataset<-rbind(dataset, temp_dataset)
      rm(temp_dataset)
    }
    
  }
  # Vectorize rownames and values, coerce into a data.frame, write to table 
  c1 <- rownames( sort( table( unlist(dataset)), decreasing = TRUE ) )
  c2 <- sort(table(unlist(enriched_paths)), decreasing = TRUE)
  return_table <- cbind( as.vector(c1), as.vector(c2) )
  write.table(return_table, file = 'enriched/EnrichedPathwaysGlobal.csv', sep = ',', col.names = FALSE, quote = FALSE)
}

