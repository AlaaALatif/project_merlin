if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
#biocLite("edgeR")
BiocManager::install("edgeR")

counts_file = "/Volumes/20200922 bulk RNAseq KKersten/kelly_experiments/tcell_exhaustion_2/bulk_rnaseq_output/kallisto/merged_counts.tsv"
#counts = read.csv(counts_file)
counts = read.delim(counts_file, row.names=1)
counts_reference = read.delim("https://raw.githubusercontent.com/ucdavis-bioinformatics-training/2018-June-RNA-Seq-Workshop/master/thursday/all_counts.txt")
d0 = DGEList(counts)
d0 = calcNormFactors(d0)
cutoff = 1
drop = which(apply(cpm(d0), 1, max) < cutoff)
d = d0[-drop,]
dim(d)
snames = colnames(counts)
snames