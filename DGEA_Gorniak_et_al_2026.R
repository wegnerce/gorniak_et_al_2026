#
# *****************************************************************************
#
# Gorniak et al., 2025 - From minerals to microbe: mobilization of lanthanides, 
# uptake, storage, and gene expression dynamics
# 
# doi:
# EBI/ENA accession:
# https://github.com/wegnerce/gorniak_et_al_2026
#
# *****************************************************************************
#

# 0. NEEDED LIBRARIES ----------------------------------------------------------
# Differential gene expression analysis was done using edgeR.
# To install this package and needed dependencies we use BiocManager
# and Bioconductor a repository for open source,  R-based, Bioinformatic tools

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install(c("edgeR", "limma", "mixOmics", "HTSFilter", 
                       "RColorBrewer", "svglite"))


# 0.a LOAD LIBRARIES -----------------------------------------------------------
library(edgeR)
library(limma)
library(mixOmics)
library(HTSFilter)
library(RColorBrewer)
library(svglite)

#library(bigPint)
library(dplyr)
library(ggplot2)
library(plotly)
library(data.table)

# 1. DATA INSPECTION -----------------------------------------------------------
# Set working directory
setwd("data_files")

# We load the count data, specifying that the id column contains the row names:
counts_AL1<-read.table("readcounts_featureCounts.txt", header=TRUE, row.names="Geneid", sep ="\t")

# edgeR methods work on  particular data structure called
# a 'DGEList'. We can create such an object from our data:

d_AL1<-DGEList(counts=dplyr::select(counts_AL1, Bas.1, Bas.2, Bas.3, Flint.2, 
                                    Gad.1, Gad.2, Gad.3, LaCl.1, LaCl.2, 
                                    LaCl.3, LaOx.1, LaOx.2, LaOx.3, LaP.1,
                                    LaP.2, LaP.3, Lop.1, Lop.2, Lop.3, 
                                    MinLaCl.1, MinLaCl.2, MinLaCl.3, Mon.1,
                                    Mon.2, Mon.3, Xen.1, Xen.2, Xen.3), 
           group=c("Bas", "Bas", "Bas", "Flint", "Gad", "Gad", "Gad", "LaCl",
                   "LaCl", "LaCl", "LaOx", "LaOx", "LaOx", "LaP", "LaP", "LaP",
                   "Lop", "Lop", "Lop", "MinLaCl", "MinLaCl",
                   "MinLaCl", "Mon", "Mon", "Mon", "Xen", "Xen",
                   "Xen"), 
           genes=data.frame(Length=dplyr::select(counts_AL1, Length)))

# The DGEList can be best imagined as container for our data, containing
# information such as gene IDs, counts, normalization fators (once they
# are computed...) and so on.

# We can have an overview of the DGEList objects:
d_AL1

# The next couple of lines of code are used to explore our datasets a bit
# Calculate pseudocounts, log2(K+1), and explore them with a histogram and a boxplot
# NOTE: These pseudocounts are generated based on the read counts and not applied
#       to the edgeR DGEList object, edgeR prefers "vanilla" data as it is
#       normalizing the data during differential analysis.
pseudoCounts <- log2(d_AL1$counts+1)
head(pseudoCounts)
hist(pseudoCounts[, "Bas.1"])
hist(pseudoCounts[, "Bas.2"])
hist(pseudoCounts[, "Bas.3"])
hist(pseudoCounts[, "Flint.2"])
hist(pseudoCounts[, "Gad.1"])
hist(pseudoCounts[, "Gad.2"])
hist(pseudoCounts[, "Gad.3"])
hist(pseudoCounts[, "LaCl.1"])
hist(pseudoCounts[, "LaCl.2"])
hist(pseudoCounts[, "LaCl.3"])
hist(pseudoCounts[, "LaOx.1"])
hist(pseudoCounts[, "LaOx.2"])
hist(pseudoCounts[, "LaOx.3"])
hist(pseudoCounts[, "LaP.1"])
hist(pseudoCounts[, "LaP.2"])
hist(pseudoCounts[, "LaP.3"])
hist(pseudoCounts[, "Lop.1"])
hist(pseudoCounts[, "Lop.2"])
hist(pseudoCounts[, "Lop.3"])
hist(pseudoCounts[, "MinLaCl.1"])
hist(pseudoCounts[, "MinLaCl.2"])
hist(pseudoCounts[, "MinLaCl.3"])
hist(pseudoCounts[, "Mon.1"])
hist(pseudoCounts[, "Mon.2"])
hist(pseudoCounts[, "Mon.3"])
hist(pseudoCounts[, "Xen.1"])
hist(pseudoCounts[, "Xen.2"])
hist(pseudoCounts[, "Xen.3"])

# Check inter-sample relationships
plotMDS(pseudoCounts)
sampleDists <- as.matrix(dist(t(pseudoCounts)))
sampleDists

# 2. Differential gene expression analysis -------------------------------------
# Let's do this differential gene expression analysis
# The edgeR User guide / vignette is really excellent with a lot
# of useful references, we strongly recommend reading it to get
# a basic understanding how data is processed and why.
# # https://bioconductor.org/packages/release/bioc/vignettes/edgeR/inst/doc/edgeRUsersGuide.pdf

d_AL1$counts <- d_AL1$counts[apply(d_AL1$counts, 1, sum) != 0, ]

d_AL1 <- calcNormFactors(d_AL1)
d_AL1 <- estimateCommonDisp(d_AL1)
d_AL1 <- estimateTagwiseDisp(d_AL1)

d_rpkm_AL1 <- rpkm(d_AL1, log = TRUE)
d_cpm_AL1 <- cpm(d_AL1, log = TRUE)

write.table(d_rpkm_AL1, file="rpkm_log_AL1.csv", sep = "\t", quote = FALSE)
write.table(d_cpm_AL1, file="cpm_log_AL1.csv", sep = "\t", quote = FALSE)

# A word about normalization factors (quoting the user guide):
# ============================================================
# The calcNormFactors function normalizes for RNA composition by finding a set of scaling factors
# for the library sizes that minimize the log-fold changes between the samples for most genes. The
# default method for computing these scale factors uses a trimmed mean of M-values (TMM) between
# each pair of samples. We call the product of the original library size and the scaling factor the
# effective library size. The effective library size replaces the original library size in all downstream
# analyses.
# TMM is the recommended for most RNA-Seq data where the majority (more than half) of
# the genes are believed not differentially expressed between any pair of the samples. 
# The normalization factors of all the libraries multiply to unity. A normalization factor below
# one indicates that a small number of high count genes are monopolizing the sequencing, causing
# the counts for other genes to be lower than would be usual given the library size. As a result, the
# library size will be scaled down, analogous to scaling the counts upwards in that library. Conversely,
# a factor above one scales up the library size, analogous to downscaling the counts.

# And a word about pseudo counts in edgeR (remember what I said above):
# ======================================================================
# The classic edgeR functions estimateCommonDisp and exactTest produce a matrix of pseudo-counts
# as part of the output object. The pseudo-counts are used internally to speed up computation of the
# conditional likelihood used for dispersion estimation and exact tests in the classic edgeR pipeline.
# The pseudo-counts represent the equivalent counts would have been observed had the library sizes
# all been equal, assuming the fitted model. The pseudo-counts are computed for a specific purpose,
# and their computation depends on the experimental design as well as the library sizes, so users
# are advised not to interpret the psuedo-counts as general-purpose normalized counts. They are
# intended mainly for internal use in the edgeR pipeline.
d_AL1

plotBCV(d_AL1)

# The biological coefficient of variation
# =======================================
# Biological CV (BCV) is the coefficient of variation with which the (unknown) true abundance
# of the gene varies between replicate RNA samples. It represents the CV that would remain between
# biological replicates if sequencing depth could be increased indefinitely. The technical CV decreases
# as the size of the counts increases. BCV on the other hand does not. BCV is therefore likely to be
# the dominant source of uncertainty for high-count genes, so reliable estimation of BCV is crucial
# for realistic assessment of differential expression in RNA-Seq experiments. If the abundance of each
# gene varies between replicate RNA samples in such a way that the genewise standard deviations are
# proportional to the genewise means, a commonly occurring property of measurements on physical
# quantities, then it is reasonable to suppose that BCV is approximately constant across genes.
# We allow however for the possibility that BCV might vary between genes and might also show a
# systematic trend with respect to gene expression or expected count.

# When you look at the plot, you see a lot of variation for genes with low logCPM, that makes perfect
# sense as we have a lot of "noise" here. Starting from a logCPM of 2-3 we have very low variation.
# And this is how it should be.

dtest_LaCl_LaOx <- exactTest(d_AL1, pair = c("LaCl", "LaOx"))
dtest_LaCl_LaP <- exactTest(d_AL1, pair = c("LaCl", "LaP"))

dtest_MinLa_Bas <- exactTest(d_AL1, pair = c("MinLaCl", "Bas"))
dtest_MinLa_Flint <- exactTest(d_AL1, pair = c("MinLaCl", "Flint"))
dtest_MinLa_Gad <- exactTest(d_AL1, pair = c("MinLaCl", "Gad"))
dtest_MinLa_Lop <- exactTest(d_AL1, pair = c("MinLaCl", "Lop"))
dtest_MinLa_Mon <- exactTest(d_AL1, pair = c("MinLaCl", "Mon"))
dtest_MinLa_Xen <- exactTest(d_AL1, pair = c("MinLaCl", "Xen"))

# The exactTest function is essentially a Fisher's exact test applied on RNA-Seq data. We have two
# nominal variables, namely our treatments/conditions (co-culture vs. alone) and our genes. We now
# test whether the treatment has an effect on our read counts per gene.
# For a better idea about Fisher's exact test I recommend:
# http://www.biostathandbook.com/fishers.html

d_AL1$genes = NULL

# ... Flint --> Gad
dfiltered <- HTSFilter(d_AL1, conds=c("Bas", "Bas", "Bas", "Gad", "Gad", 
                                      "Gad", "Gad", "LaCl", "LaCl", "LaCl",
                                      "LaOx", "LaOx", "LaOx", "LaP", "LaP",
                                      "LaP", "Lop", "Lop", "Lop", "MinLaCl",
                                      "MinLaCl", "MinLaCl", "Mon", "Mon", "Mon",
                                      "Xen", "Xen", "Xen"))$filteredData

dtest_filtered_LaCl_LaOx <- exactTest(dfiltered, pair = c("LaCl", "LaOx"))
dtest_filtered_LaCl_LaP <- exactTest(dfiltered, pair = c("LaCl", "LaP"))

dtest_filtered_MinLa_Bas <- exactTest(dfiltered, pair = c("MinLaCl", "Bas"))
dtest_filtered_MinLa_Flint <- exactTest(dfiltered, pair = c("MinLaCl", "Flint"))
dtest_filtered_MinLa_Gad <- exactTest(dfiltered, pair = c("MinLaCl", "Gad"))
dtest_filtered_MinLa_Lop <- exactTest(dfiltered, pair = c("MinLaCl", "Lop"))
dtest_filtered_MinLa_Mon <- exactTest(dfiltered, pair = c("MinLaCl", "Mon"))
dtest_filtered_MinLa_Xen <- exactTest(dfiltered, pair = c("MinLaCl", "Xen"))

# One issue with DGE is that hypothesis testing is done on a gene-by-gene basis. As a result determined
# p-values must be corrected for multiple testing. Common correction procedures to identify false positives
# lead to a loss of power to detect truly differentially expressed genes, due to high numbers of hypothesis
# tests performed.

# To prevent this we use HTSFilter to remove background noise (genes with a constant, low expression value).
# A quote from the HTSFilter manual:
# The HTSFilter package implements a novel data-based filtering procedure based
# on the calculation of a similarity index among biological replicates for read
# counts arising from replicated transcriptome sequencing (RNA-seq) data. This
# technique provides an intuitive data-driven way to filter high-throughput transcriptome
# sequencing data and to effectively remove genes with low, constant
# expression levels without incorrectly removing those that would otherwise have
# been identified as DE.

# Check the p-values
# Interesting is always the proportion of low p-value genes
hist(dtest_LaCl_LaOx$table[,"PValue"], breaks =50)
hist(dtest_LaCl_LaP$table[,"PValue"], breaks =50)
hist(dtest_MinLa_Bas$table[,"PValue"], breaks =50)
hist(dtest_MinLa_Flint$table[,"PValue"], breaks =50)
hist(dtest_MinLa_Gad$table[,"PValue"], breaks =50)
hist(dtest_MinLa_Lop$table[,"PValue"], breaks =50)
hist(dtest_MinLa_Mon$table[,"PValue"], breaks =50)
hist(dtest_MinLa_Xen$table[,"PValue"], breaks =50)

hist(dtest_filtered_LaCl_LaOx$table[,"PValue"], breaks =50)
hist(dtest_filtered_LaCl_LaP$table[,"PValue"], breaks =50)
hist(dtest_filtered_MinLa_Bas$table[,"PValue"], breaks =50)
hist(dtest_filtered_MinLa_Flint$table[,"PValue"], breaks =50)
hist(dtest_filtered_MinLa_Gad$table[,"PValue"], breaks =50)
hist(dtest_filtered_MinLa_Lop$table[,"PValue"], breaks =50)
hist(dtest_filtered_MinLa_Mon$table[,"PValue"], breaks =50)
hist(dtest_filtered_MinLa_Xen$table[,"PValue"], breaks =50)

# Extract a summary of DGE stats
resNoFilt_LaCl_LaOx <- topTags(dtest_LaCl_LaOx, n=nrow(dtest_LaCl_LaOx$table))
resNoFilt_LaCl_LaP <- topTags(dtest_LaCl_LaP, n=nrow(dtest_LaCl_LaP$table))

resNoFilt_MinLa_Bas <- topTags(dtest_MinLa_Bas, n=nrow(dtest_MinLa_Bas$table))
resNoFilt_MinLa_Flint <- topTags(dtest_MinLa_Flint, 
                                 n=nrow(dtest_MinLa_Flint$table))
resNoFilt_MinLa_Gad <- topTags(dtest_MinLa_Gad, n=nrow(dtest_MinLa_Gad$table))
resNoFilt_MinLa_Lop <- topTags(dtest_MinLa_Lop, n=nrow(dtest_MinLa_Lop$table))
resNoFilt_MinLa_Mon <- topTags(dtest_MinLa_Mon, n=nrow(dtest_MinLa_Mon$table))
resNoFilt_MinLa_Xen <- topTags(dtest_MinLa_Xen, n=nrow(dtest_MinLa_Xen$table))

resFilt_LaCl_LaOx <- topTags(dtest_filtered_LaCl_LaOx, 
                             n=nrow(dtest_filtered_LaCl_LaOx$table))
resFilt_LaCl_LaP <- topTags(dtest_filtered_LaCl_LaP, 
                            n=nrow(dtest_filtered_LaCl_LaP$table))

resFilt_MinLa_Bas <- topTags(dtest_filtered_MinLa_Bas, 
                             n=nrow(dtest_filtered_MinLa_Bas$table))
resFilt_MinLa_Flint <- topTags(dtest_filtered_MinLa_Flint, 
                                 n=nrow(dtest_filtered_MinLa_Flint$table))
resFilt_MinLa_Gad <- topTags(dtest_filtered_MinLa_Gad, 
                             n=nrow(dtest_filtered_MinLa_Gad$table))
resFilt_MinLa_Lop <- topTags(dtest_filtered_MinLa_Lop, 
                             n=nrow(dtest_filtered_MinLa_Lop$table))
resFilt_MinLa_Mon <- topTags(dtest_filtered_MinLa_Mon, 
                             n=nrow(dtest_filtered_MinLa_Mon$table))
resFilt_MinLa_Xen <- topTags(dtest_filtered_MinLa_Xen, 
                             n=nrow(dtest_filtered_MinLa_Xen$table))

# topTags extracts top DE tags ranked by p-value of LogFC, by default sorting is
# based on p-values.

# How many differentially expressed genes do we have?
sum(resNoFilt_LaCl_LaOx$table$FDR < 0.05)
sum(resNoFilt_LaCl_LaP$table$FDR < 0.05)

sum(resNoFilt_MinLa_Bas$table$FDR < 0.05)
sum(resNoFilt_MinLa_Flint$table$FDR < 0.05)
sum(resNoFilt_MinLa_Gad$table$FDR < 0.05)
sum(resNoFilt_MinLa_Lop$table$FDR < 0.05)
sum(resNoFilt_MinLa_Mon$table$FDR < 0.05)
sum(resNoFilt_MinLa_Xen$table$FDR < 0.05)

# How many differentially expressed genes do we have?
sum(resFilt_LaCl_LaOx$table$FDR < 0.05)
sum(resFilt_LaCl_LaP$table$FDR < 0.05)

sum(resFilt_MinLa_Bas$table$FDR < 0.05)
sum(resFilt_MinLa_Flint$table$FDR < 0.05)
sum(resFilt_MinLa_Gad$table$FDR < 0.05)
sum(resFilt_MinLa_Lop$table$FDR < 0.05)
sum(resFilt_MinLa_Mon$table$FDR < 0.05)
sum(resFilt_MinLa_Xen$table$FDR < 0.05)

# Extract up-/downregulated genes
# Here three parameters can be used for filtering:
# ** FDR (=p-value)
# ** logFC 
# ** logCPM
# Feel free to modify these values!
sigReg_LaCl_LaOx <- resFilt_LaCl_LaOx$table[(resFilt_LaCl_LaOx$table$FDR<0.01) &
                              (abs(resFilt_LaCl_LaOx$table$logFC)>1) & 
                              (resFilt_LaCl_LaOx$table$logCPM>4),]

sigReg_LaCl_LaP <- resFilt_LaCl_LaP$table[(resFilt_LaCl_LaP$table$FDR<0.01) &
                            (abs(resFilt_LaCl_LaP$table$logFC)>1) & 
                            (resFilt_LaCl_LaP$table$logCPM>4),]

sigReg_MinLa_Bas <- resFilt_MinLa_Bas$table[(resFilt_MinLa_Bas$table$FDR<0.01) &
                              (abs(resFilt_MinLa_Bas$table$logFC)>1) & 
                              (resFilt_MinLa_Bas$table$logCPM>4),]

sigReg_MinLa_Flint <- resFilt_MinLa_Flint$table[(resFilt_MinLa_Flint$table$FDR<0.01) &
                            (abs(resFilt_MinLa_Flint$table$logFC)>1) & 
                            (resFilt_MinLa_Flint$table$logCPM>4),]

sigReg_MinLa_Gad <- resFilt_MinLa_Gad$table[(resFilt_MinLa_Gad$table$FDR<0.01) &
                              (abs(resFilt_MinLa_Gad$table$logFC)>1) & 
                              (resFilt_MinLa_Gad$table$logCPM>4),]

sigReg_MinLa_Lop <- resFilt_MinLa_Lop$table[(resFilt_MinLa_Lop$table$FDR<0.01) &
                            (abs(resFilt_MinLa_Lop$table$logFC)>1) & 
                            (resFilt_MinLa_Lop$table$logCPM>4),]

sigReg_MinLa_Mon <- resFilt_MinLa_Mon$table[(resFilt_MinLa_Mon$table$FDR<0.01) &
                              (abs(resFilt_MinLa_Mon$table$logFC)>1) & 
                              (resFilt_MinLa_Mon$table$logCPM>4),]

sigReg_MinLa_Xen <- resFilt_MinLa_Xen$table[(resFilt_MinLa_Xen$table$FDR<0.01) &
                            (abs(resFilt_MinLa_Xen$table$logFC)>1) & 
                            (resFilt_MinLa_Xen$table$logCPM>4),]

head(sigReg_LaCl_LaOx)
head(sigReg_LaCl_LaP)

head(sigReg_MinLa_Bas)
head(sigReg_MinLa_Flint)
head(sigReg_MinLa_Gad)
head(sigReg_MinLa_Lop)
head(sigReg_MinLa_Mon)
head(sigReg_MinLa_Xen)

list_of_DEGs <- c(row.names(sigReg_LaCl_LaOx), row.names(sigReg_LaCl_LaP),
                  row.names(sigReg_MinLa_Bas), row.names(sigReg_MinLa_Flint),
                  row.names(sigReg_MinLa_Gad), row.names(sigReg_MinLa_Lop),
                  row.names(sigReg_MinLa_Mon), row.names(sigReg_MinLa_Xen))

# Lists of differentially expressed genes (up- and downregulated) are saved in our current working
# directory.
write.table(sigReg_LaCl_LaOx, file="sigReg_LaCl_vs_LaOx.tsv", sep ="\t", quote = FALSE)
write.table(sigReg_LaCl_LaP, file="sigReg_LaCl_vs_LaP.tsv", sep ="\t", quote = FALSE)

write.table(sigReg_MinLa_Bas, file="sigReg_MinLa_vs_Bas.tsv", sep ="\t", quote = FALSE)
write.table(sigReg_MinLa_Flint, file="sigReg_MinLa_vs_Flint.tsv", sep ="\t", quote = FALSE)
write.table(sigReg_MinLa_Gad, file="sigReg_MinLa_vs_Gad.tsv", sep ="\t", quote = FALSE)
write.table(sigReg_MinLa_Lop, file="sigReg_MinLa_vs_Lop.tsv", sep ="\t", quote = FALSE)
write.table(sigReg_MinLa_Mon, file="sigReg_MinLa_vs_Mon.tsv", sep ="\t", quote = FALSE)
write.table(sigReg_MinLa_Xen, file="sigReg_MinLa_vs_Xen.tsv", sep ="\t", quote = FALSE)
