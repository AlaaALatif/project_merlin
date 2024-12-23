# If not installed, uncomment and install as needed:
# install.packages("BiocManager")
# BiocManager::install("limma")
# BiocManager::install("edgeR")

library(limma)
library(edgeR)
library(tidyverse)


count_data <- read.csv("/Users/alaa/Documents/ucsf/data/suliman/merlin/counts/experiment1_raw_counts_2024-12-23.csv", header = TRUE, stringsAsFactors = FALSE)

# If the first column is "patient_id" or "sample_id", you might want to move it out:
# For example, if the first column is sample_id:
rownames(count_data) <- count_data$sample_id
count_data <- count_data[ , -1] 

# If genes are currently columns and samples are rows, transpose:
count_data <- t(count_data)
head(count_data)
# Load cohort metadata
metadata <- read.csv("/Users/alaa/Documents/ucsf/data/suliman/merlin/meta/experiment1_processed_metadata_2024-12-23.csv", 
                     header = TRUE, stringsAsFactors = FALSE)
# Verify alignment
all(colnames(count_data) == metadata$sample_id)
# Option 1: reorder metadata to match the count_data column names:
# metadata <- metadata[match(colnames(count_data), metadata$sample_id), ]
# store count data in a DGEList object
dge <- DGEList(counts = count_data)
# Example filter: keep genes with >= 10 counts in at least 3 samples
keep <- rowSums(cpm(dge) >= 1) >= 3
dge <- dge[keep, , keep.lib.sizes=FALSE]
# Compute normalization factors
dge <- calcNormFactors(dge)
# Setup the design matrix
metadata$hiv_diagnosis <- factor(metadata$hiv_diagnosis)
metadata$library_pool  <- factor(metadata$library_pool)
# Create the design matrix
design <- model.matrix(~ hiv_diagnosis + library_pool, data = metadata)
# The block argument is patient_id
patient_ids <- factor(metadata$patient_id)
# We first do a voom transformation *without* fitting the final model yet
v <- voom(dge, design, plot=FALSE)
# Then estimate the correlation
corfit <- duplicateCorrelation(v, design, block = patient_ids)
corfit$consensus
# Apply voom again, now specifying the correlation structure
v <- voom(dge, design, plot=FALSE, block=patient_ids, correlation=corfit$consensus)
# Fit the Linear Model
fit <- lmFit(v, design, block = patient_ids, correlation = corfit$consensus)
# Apply empirical Bayes moderation
fit <- eBayes(fit)
# Check design matrix columns to see which coefficient is for hiv_diagnosis
colnames(design)
# Typically, if we have (Intercept), hiv_diagnosisPositive, library_poolPool2...
# then hiv_diagnosisPositive might be coefficient 2.
deg_results <- topTable(fit, coef = 2, number = Inf, adjust.method = "BH")
deg_results
# Look at the first few
head(deg_results)
deg_results_output_filepath = "/Users/alaa/Documents/ucsf/data/suliman/merlin/dge/experiment1_dge_results_hivdx_library_pool_duplicationrate_patient_id_limma_voom_2024-12-23.csv"
write.csv(deg_results, deg_results_output_filepath)
