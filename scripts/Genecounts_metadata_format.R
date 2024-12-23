rm (list= ls())
setwd("~/Library/CloudStorage/Box-Box/UCSF Postdoc/Projects/MERLIN/Results/RNAseq analyses/DGE")

#### Read tpm counts and log transform data 
genecounts <- read.csv("gene counts_r1.csv", header = TRUE, stringsAsFactors = F)

####  set gene symbol as row names
genecounts$gene_symbol <- make.unique (as.character(genecounts$gene_symbol)) ###makes repeated gene names unique
row.names(genecounts) <- genecounts$gene_symbol ##set as row names
genecounts$gene_symbol <- NULL #remove row name

#upload metadata, use filtered list for PCA standard method
metadata <-  read.csv ("metadata_round1.csv", header = TRUE, stringsAsFactors = T)

### edit the column names in samples 
colnames (genecounts)

#remove the X from all the column names
colnames(genecounts) <- sub("^X","", colnames(genecounts))

##### check the head of the  metadata and samples to make sure the order of the samples are the same 
head(metadata$samplecode)
head(colnames(genecounts))

### only consider samples in the different groups
###First order samples 
metadata<- metadata[order(metadata$samplecode),]
genecounts <- genecounts[, order(colnames(genecounts))]

##### Some of the colnames on samples have zeros, remove the zeros
colnames(genecounts) <- gsub("^0+","", colnames(genecounts))
head(colnames(genecounts))

#### check if column names in samples match sample codes in metadata
all(colnames(genecounts) %in% metadata$samplecode) ## came out false so do some more checking 
#
#check for mismatches between the two datasets 
setdiff(metadata$samplecode, colnames(genecounts)) #which samplecodes in metadata are not in samples
setdiff(colnames(genecounts), metadata$samplecode) #which colnnames in samples are not in metadata

#remove all leading and trailing spaces 
metadata$samplecode <- trimws(metadata$samplecode)
colnames(genecounts) <- trimws(colnames(genecounts))

####remove all . in samplecode
colnames(genecounts) <- gsub("\\.","-", colnames(genecounts))

#remove the second underscores
colnames(genecounts) <- as.character(colnames(genecounts))

colnames(genecounts) <- gsub("^(.*?_.*?)(_)(.*)$", "\\1 \\3", colnames(genecounts))

#### work on the last column that has a mismatch
###identify column 
undesired_column <-"Nuclease_free water_" # Replace with specific column name 

#Remove underscores only for that column 
colnames(genecounts)[colnames(genecounts)==undesired_column] <- gsub("_"," ", undesired_column)

#remove all leading and trailing spaces 
metadata$samplecode <- trimws(metadata$samplecode)
colnames(genecounts) <- trimws(colnames(genecounts))

#check for mismatches between the two datasets 
setdiff(metadata$samplecode, colnames(genecounts)) #which samplecodes in metadata are not in samples
setdiff(colnames(genecounts), metadata$samplecode) #which colnnames in samples are not in metadata

###reorder samples 
metadata<- metadata[order(metadata$samplecode),]
genecounts <- genecounts[, order(colnames(genecounts))]

head(metadata$samplecode)
head(colnames(genecounts))

###
filteredList <- metadata$samplecode
genecounts <- genecounts[, colnames(genecounts) %in% filteredList]
 
##### order is the same
remaining_samples <-colnames(genecounts)
all(remaining_samples==metadata$samplecode)

### format the genecounts table and combine with metadata
genes <- t(genecounts)
genes_metadata<- cbind(genes, metadata)

### Add in a column for pre and post 
pre_HIV <- c("BL", "X-1")
post_HIV <- c("DX", "M12", "M24", "M24 ART", "M3", "M48", "M48 ART", "M6")#

#Add a new pre/post column
genes_metadata$pre_post <- ifelse(
  genes_metadata$Time.point %in% pre_HIV, 
  "Pre",
  ifelse(
    genes_metadata$Time.point %in% post_HIV,
    "Post",
    NA
  )
)

#### 
write.csv(genes_metadata,"~/Library/CloudStorage/Box-Box/UCSF Postdoc/Projects/MERLIN/Results/RNAseq analyses/DGE/genes_metadata.csv", row.names=T )
