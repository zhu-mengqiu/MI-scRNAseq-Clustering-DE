# Download NBT_hiseq_linear_tpm_values.txt from the the following link
# https://www.pollenlab.org/datasets
# counts table of scRNA-seq data from Pollen et al., 2014

# set working directory to the folder where NBT_hiseq_linear_tpm_values.txt is 
# stored 

# load packages
library(data.table)

# load expression value data
data <- read.table('NBT_hiseq_linear_tpm_values.txt', header = T, sep  = '\t')
genes <- data[,1]
data <- data[,-1]

# extract cell type labels from column names
cells <- colnames(data)
ncells <- length(cells)
label <- vector(mode = 'character', length = ncells)
for (i in 1:ncells){
  label[i] <- strsplit(cells[i], '_')[[1]][2]
}

# remove genes that have zero expression values across all cells
idx <- which(apply(data, 1, function(r) mean(r==0)) == 1)
data <- data[-idx,]
genes <- genes[-idx]

# write expression values and cell type labels to text files
# expression values: Pollen_raw.txt
# cell type labels: Pollen_label.txt
data <- data.table(genes, data)
cells <- paste0('CELL_', 1:length(label))
colnames(data) <- c('GENE_ID', cells)
write.table(data, file = 'Pollen_raw.txt', sep = '\t', 
            row.names = F, col.names = T, quote = F)
write.table(label, 'Pollen_label.txt', quote = F, sep = '\n', col.names = F,
            row.names = F)