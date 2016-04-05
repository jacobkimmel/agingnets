# New-style PASCAL output converter to old-style output

data_file = '/Applications/PASCAL/output/GIANT_HEIGHT_Wood_et_al_2014_pvals.sum.genescores.chr22.txt'
convert_output <- function(dir, dataset){
  
  name <- paste(dir, dataset, sep = '')
  df <- read.table(name, header = TRUE)
  df$Status <- "Simul"
  df$gene_id <- df$gene_symbol
  
  filename <- paste("converted_pascal/", dataset, sep = '')
  write.table(df, file = filename, quote = FALSE, sep = "\t", row.names = FALSE)
  
}