# Find and download NBT_hiseq_linear_tpm_values.txt from the the following link
# https://www.pollenlab.org/datasets
# Pollen et al., 2014

library(data.table)

# load expression value data
data <- read.table('NBT_hiseq_linear_tpm_values.txt', header = T, sep  = '\t')
genes <- data[,1]
data <- data[,-1]

# extract labels from column names
cells <- colnames(data)
ncells <- length(cells)
label <- vector(mode = 'character', length = ncells)
for (i in 1:ncells){
  label[i] <- strsplit(cells[i], '_')[[1]][2]
}

# remove genes that are zero across all cells
idx <- which(apply(data, 1, function(r) mean(r==0)) == 1)
data <- data[-idx,]
genes <- genes[-idx]

# write expression values and cell labels to text files
data <- data.table(genes, data)
cells <- paste0('CELL_', 1:length(label))
colnames(data) <- c('GENE_ID', cells)
write.table(data, file = 'Pollen_raw.txt', sep = '\t', 
            row.names = F, col.names = T, quote = F)
write.table(label, 'Pollen_label.txt', quote = F, sep = '\n', col.names = F,
            row.names = F)

# label distribution
label <- read.table('Pollen_label.txt')$V1
type <- unique(label)
sapply(type, function(x) sum(label == x))
