####################
# Load packages used
####################

# Check if present, if not: require, install, and load dplyr
if(!require(dplyr)){
  install.packages("dplyr")
  library(dplyr)
}

# Check if present, if not: require, install, and load openxlsx
if(!require(openxlsx)){
  install.packages("openxlsx")
  library(openxlsx)
}

############################################################
# Define directories (Change for your specific directories!)
############################################################

# Define working directory 
working_dir <- "C:/Users/aliad/OneDrive/Documents/R_PROJECTS/MCB-540_Final-Project_R"

# Define data directory 
data_dir <- "C:/Users/aliad/OneDrive/Documents/R_PROJECTS/MCB-540_Final-Project_R/data/"

# Define results directory
results_dir <- "C:/Users/aliad/OneDrive/Documents/R_PROJECTS/MCB-540_Final-Project_R/results/"

###################
# Read in DATA SETS
###################

# Read in the list of significant genes
sig_genes_mdamb231_data <- read.delim(paste(data_dir,"sig_genes_mdamb231.tsv", sep = "/"), header = TRUE, sep = "\t")

# Read in the data of the entire list of genes 
all_genes_data <- read.delim(paste(data_dir,"GSE212143.top.table.tsv", sep = "/"), header = TRUE, sep = "\t")

##############################
## Query All SIGNIFICANT GENES
##############################

# Create a look up criteria to query all significant genes
lookup_gene_bySymbol <- sig_genes_mdamb231_data$Symbol

# All significant genes data
mdamb231_all_sig_genes_data <- subset(all_genes_data, Symbol %in% lookup_gene_bySymbol)

########################
## Pre-process DATA SETS
########################

# Correct for mismatches between sig_genes_mdamb231_data and mdamb231_sig_genes_data data sets

# Rename log2FoldChange and n1Log10AdjPvalue column names in sig_genes_mdamb231_data for clarity
colnames(sig_genes_mdamb231_data)[4] <- "log2FoldChange"
colnames(sig_genes_mdamb231_data)[5] <- "n1Log10AdjPvalue"

# Create n1Log10AdjPvalue column in mdamb231_sig_genes_data
mdamb231_sig_genes_data <- mdamb231_all_sig_genes_data %>%
  mutate(n1Log10AdjPvalue = -log10(padj))

# Create a vector of column names, placing n1Log10AdjPvalue after padj column
new_order <- c(setdiff(names(mdamb231_sig_genes_data), "n1Log10AdjPvalue"),
               "n1Log10AdjPvalue")

# Find the position of padj column
target_index <- which(names(mdamb231_sig_genes_data) == "padj")

# Split the vector before and after the target column, and insert the new column
new_order <- c(new_order[1:target_index], "n1Log10AdjPvalue", new_order[(target_index + 1):length(new_order)])

# Reorder the columns in the data frame
mdamb231_sig_genes_data <- mdamb231_sig_genes_data[new_order]

# Exclude the duplicate last column
mdamb231_sig_genes_data <- mdamb231_sig_genes_data[, -ncol(mdamb231_sig_genes_data)]

############################################################
## Query Significant Genes of ALL SELECTED REGIONS (1,2,3,4)
## REGION 1: LogFC < -13
## REGION 2: -log10(Pvalue) > 22
## REGION 3: LogFC > 0 AND -log10(PValue) > 14
## REGION 4: LogFC > 20
############################################################

# Region 1: LogFC < -13
region1_data <- subset(mdamb231_sig_genes_data, log2FoldChange < -13)

# Region 2: n1log10PValue > 22
region2_data <- subset(mdamb231_sig_genes_data, n1Log10AdjPvalue > 22)

# Region 3: LogFC > 0 AND n1log10PValue > 14
region3_data <- subset(mdamb231_sig_genes_data, log2FoldChange > 0 & n1Log10AdjPvalue > 14)

# Region 4: LogFC > 20
region4_data <- subset(mdamb231_sig_genes_data, log2FoldChange > 20)

########################################################################################
# Write out all reformed data sets and findings into human-readable (Excel) spreadsheets 
########################################################################################

# Write out 'all genes' data (21,206 genes)
write.xlsx(all_genes_data, paste(results_dir,"all_genes.xlsx", sep = "/"))

# Write out 'all significant genes' data (1,262 genes)
write.xlsx(mdamb231_sig_genes_data, paste(results_dir,"sig_genes.xlsx", sep = "/"))

# Write out findings...

# Write out 'region 1' data 
write.xlsx(region1_data, paste(results_dir,"r1_genes.xlsx", sep = "/"))

# Write out 'region 2' data
write.xlsx(region2_data, paste(results_dir,"r2_genes.xlsx", sep = "/"))

# Write out 'region 3' data
write.xlsx(region3_data, paste(results_dir,"r3_genes.xlsx", sep = "/"))

# Write out 'region 4' data
write.xlsx(region4_data, paste(results_dir,"r4_genes.xlsx", sep = "/"))
