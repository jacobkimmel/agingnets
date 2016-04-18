# Read in p-vals file
require(stringr)

# Get list of GWAS's from magnum_output/
# Names in the format of Framingham Heart Study read-outs -- "pha####"
gwas <- list.files("magnum_output/")
gwas_no1788 <- gwas[-2]


# Create matrix of FHS-style names with their corresponding read-out
fhs_codes <- c("pha1754", "pha1788", "pha1970", "pha1982", "pha1986", "pha1990", "pha1994")
read_outs <- c("Bone Ultrasound attenuation", "Morbidity Free Survival at age 65", "Hand Grip", "Disability, walking speed", "Age at Natural Menopause", "Death Age", "Death Past Average Life Expectancy")
gwas_names <- cbind(fhs_codes, read_outs)

# Loop through all 32-high level network outputs 
# and produce bar plots for enrichment 
for (line in gwas){
  pvals = paste("magnum_output/", line, "/", "phs000007.pha00", str_sub(line, -4, -1), ".sum.genescores.pvals.txt", sep = '')
  df <- read.delim(file = pvals, sep = "\t", skip = 1)
  df.sort <- df[order(df$Pvalue, decreasing = TRUE),]
  df.sort$nlogPvalue <- -log10(df.sort$Pvalue)
  # Get read out description from gwas_names matrix
  read_out <- gwas_names[which(gwas_names[,1] == line), ][2]
  outname = paste("magnum_plots/", line, "_EnrichmentPlot.png", sep = '')
  png(outname, height = 8.5, width = 15, units = 'in', res = 300)
  # Set left margin to 4X usual value, las = 1 makes labeles parallel to bars 
  # Barplot -log10(Pvalues)
  par(mai=c(1,4,1,1))
  barplot(df.sort$nlogPvalue, main = read_out, horiz = TRUE, names.arg = df$Network, las = 1, xlab = "-log(P-value)", col = ifelse(df.sort$nlogPvalue < -log10(0.05), "gray", "blue"))
  # Adds a vertical red line for significance at p = 0.05 
  abline(v = -log10(0.05), col = 'red')
  dev.off()
}

for (line in gwas_no1788){
  specpvals <- paste("magnum_output/", line, '/', str_sub(line, -4, -1), '_specific_networks/', "phs000007.pha00", str_sub(line, -4, -1), ".sum.genescores.pvals.txt", sep = '')
  df <- read.delim(file = specpvals, sep = "\t", skip = 1)
  df.sort <- df[order(df$Pvalue, decreasing = TRUE),]
  df.sort$nlogPvalue <- -log10(df.sort$Pvalue)
  # Get read out description from gwas_names matrix
  read_out <- gwas_names[which(gwas_names[,1] == line), ][2]
  outname = paste("magnum_plots/", line, "_SpecificNetworks_EnrichmentPlot.png", sep = '')
  png(outname, height = 8.5, width = 15, units = 'in', res = 300)
  # Set left margin to 4X usual value, las = 1 makes labeles parallel to bars 
  # Barplot -log10(Pvalues)
  par(mai=c(1,3.9,1,1))
  barplot(df.sort$nlogPvalue, main = paste(read_out, "Specific Networks", sep = ' '), horiz = TRUE, names.arg = df$Network, las = 1, xlab = "-log(P-value)", col = ifelse(df.sort$nlogPvalue < -log10(0.05), "gray", "blue"))
  # Adds a vertical red line for significance at p = 0.05 
  abline(v = -log10(0.05), col = 'red')
  dev.off()
}

# Write enriched networks to seperate files 
for (line in gwas_no1788){
  pvals = paste("magnum_output/", line, "/", "phs000007.pha00", str_sub(line, -4, -1), ".sum.genescores.pvals.txt", sep = '')
  specpvals <- paste("magnum_output/", line, '/', str_sub(line, -4, -1), '_specific_networks/', "phs000007.pha00", str_sub(line, -4, -1), ".sum.genescores.pvals.txt", sep = '')
  
  df <- read.table(file = pvals, sep = "\t", skip = 1, header = TRUE)
  df.sort <- df[order(df$Pvalue, decreasing = FALSE),]
  df.spec <- read.table(file = specpvals, sep = "\t", skip = 1, header = TRUE)
  df.spec.sort <- df[order(df$Pvalue, decreasing = FALSE),]
  
  num_enriched <- length(which(df.sort$Pvalue < 0.05))
  enriched <- as.character(df.sort[1:num_enriched,1])
  num_enriched.spec <- length(which(df.spec.sort$Pvalue < 0.05))
  enriched.spec <- as.character(df.spec.sort[1:num_enriched.spec,1])
  enriched.all <- rbind(enriched, enriched.spec)
  filename <- paste('enriched_networks/', line, '_', 'EnrichedNetworks.txt', sep = '')
  write(enriched, file = filename, sep = ',')
  filename <- paste('enriched_networks/', line, '_', 'EnrichedSpecNetworks.txt', sep = '')
  write(enriched.all, file = filename, sep = ',')
  
}

# Find enriched networks detected by more than one read-out

rich_nets <- function(enriched_nets_dir){
  for (file in enriched_nets_dir){
    filenames <- list.files(enriched_nets_dir, pattern = '*.txt', full.names = TRUE)
    for (file in filenames){
      
      # if the merged dataset doesn't exist, create it
      if (!exists("dataset")){
        dataset <- read.table(file, header=FALSE, sep="\t")
      }
      
      # if the merged dataset does exist, append to it
      if (exists("dataset")){
        temp_dataset <-read.table(file, header=FALSE, sep="\t")
        dataset <- rbind(dataset, temp_dataset)
        rm(temp_dataset)
      }
      
    }
    # Vectorize rownames and values, coerce into a data.frame, write to table 
    #c1 <- rownames( sort( table( unlist(dataset)), decreasing = TRUE ) )
    #c2 <- sort(table(unlist(enriched_paths)), decreasing = TRUE)
    #return_table <- cbind( as.vector(c1), as.vector(c2) )
    #write.table(return_table, file = 'enriched_networks/EnrichedNetworksGlobal.csv', sep = ',', col.names = FALSE, quote = FALSE)
  }
  return(dataset)
}