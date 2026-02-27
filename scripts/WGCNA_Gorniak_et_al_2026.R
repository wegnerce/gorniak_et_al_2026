#
# *****************************************************************************
#
# Weighted Gene Correlation Network Analysis (WGCNA)
# based on the Horvath Lab tutorials
# ATTENTION: They are no longer online available as of 2025/07
#
# *****************************************************************************
#
# 0. NEEDED LIBRARIES ----------------------------------------------------------
# Libraries can be installed as needed through Bioconductur and the 
# BiocManager package
# --> https://www.bioconductor.org/install/
# install.packages(c("BiocManager", "tidyverse"))
# BiocManager::install("WGCNA")
# BiocManager::install("flashClust")
# devtools::install_github("kevinblighe/CorLevelPlot")
# Install clusterProfiler
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(c("clusterProfiler", "enrichplot", "cowplot", "edgeR", 
                       "WGCNA", "magrittr", "readr", "dplyr", "flashClust", 
                       "tidyverse", "CorLevelPlot", "ggsci",  "RColorBrewer"))

# 0.a LOAD LIBRARIES -----------------------------------------------------------
library(edgeR)
library(WGCNA)
library(magrittr)
library(readr)
library(dplyr)
library(flashClust)
library(tidyverse)
library(CorLevelPlot)
library(ggsci)
library(clusterProfiler)
library(enrichplot)
library(cowplot)
library(RColorBrewer)

options(stringsAsFactors = FALSE)

# 1. DATA PREPROCESSING --------------------------------------------------------
#
# A word about the do's and don'ts of data preprocessing
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Originally developed for microarray data analysis, WGCNA can be used as
# well for related data such as RNAseq-derived gene expression data.
#
# As many (if not most) statistical methods, WGCNA benefits from a larger
# sampling size. The authors recommend to use WGCNA only for datasets that
# include data from >= 30 samples. In general, the rule applies "more samples,
# more robustness".
#
# It is not recommended to pre-filter datasets based on differential gene
# expression.  Filtering genes by differential expression will lead to a set 
# of correlated genes that will essentially form a single (or a few highly 
# correlated) modules. It also completely invalidates the scale-free topology 
# assumption, so choosing soft thresholding power by scale-free topology fit 
# will fail.
#
# Instead filtering should focus on genes that introduce noise through their
# low read recruitment. A good starting point might be to remove genes with
# a mean row count < 20, or genes that recruit less than 10 reads in 90%
# of all samples.

# Set working directory
setwd("/home/wegnerce/Dropbox/tmp_work/AL1_minerals_RNAseq/submission/osf_io")

# We load the count data, specifying that the id column contains the row names:
counts_AL1<-read.table("readcounts_featureCounts_Gorniak-etal-2026.txt", header=TRUE, row.names="Geneid", sep ="\t")

# edgeR methods work on  particular data structure called
# a 'DGEList'. We can create such an object from our data:
d_AL1<-DGEList(counts=dplyr::select(counts_AL1, Bas.1, Bas.2, Bas.3, Ferro.3, 
                                    Gad.1, Gad.2, Gad.3, LaCl.1, LaCl.2, 
                                    LaCl.3, LaOx.1, LaOx.2, LaOx.3, LaP.1,
                                    LaP.2, LaP.3, Lop.1, Lop.2, Lop.3, 
                                    MinLaCl.1, MinLaCl.2, MinLaCl.3, Mon.1,
                                    Mon.2, Mon.3, Xen.1, Xen.2, Xen.3), 
               group=c("Bas", "Bas", "Bas", "Ferro", "Gad", "Gad", "Gad", "LaCl",
                       "LaCl", "LaCl", "LaOx", "LaOx", "LaOx", "LaP", "LaP", "LaP",
                       "Lop", "Lop", "Lop", "MinLaCl", "MinLaCl",
                       "MinLaCl", "Mon", "Mon", "Mon", "Xen", "Xen",
                       "Xen"), 
               genes=data.frame(Length=dplyr::select(counts_AL1, Length)))

# We remove lowly expressed genes, since they represent noise and do not
# contribute to DGE.
# Here we only keep reads that feature CPM values above 4 in all 28 datasets.
keep <- rowSums(cpm(d_AL1) > 4) >= 28
table(keep)#result: FALSE 385 and TRUE 3878
d_AL1 <- d_AL1[keep, , keep.lib.sizes=FALSE]

# We normalize the datasets, check the edgeR manual/vignette for details:
# --> https://www.bioconductor.org/packages/devel/bioc/vignettes/edgeR/inst/doc/edgeRUsersGuide.pdf
d_AL1 <- calcNormFactors(d_AL1)
d_AL1 <- estimateCommonDisp(d_AL1)
d_AL1 <- estimateTagwiseDisp(d_AL1)

# Next, the normalized data (=readcounts) are transformed in to CPM.
d_cpm_WGCNA_AL1 <- cpm(d_AL1, normalized.lib.sizes=TRUE, log=FALSE, 
                       prior.count=0.25)

# WGCNA requires us to transform the data: samples as rows, genes as columns.
expression.data <- as.data.frame(t(d_cpm_WGCNA_AL1))

# Optional: Enable multi-threading
allowWGCNAThreads()

# 2. WGCNA ---------------------------------------------------------------------
#
# --> https://bioinformaticsworkbook.org/tutorials/wgcna.html#gsc.tab=0
#
# 2.a Iterative filtering of samples and genes with too many missing entries----
gsg <- goodSamplesGenes(expression.data, verbose = 3)
if (!gsg$allOK) {
  expression.data <- expression.data[gsg$goodGenes, gsg$goodSamples]
}
# All good?
gsg$allOK

# 2.b Sample clustering---------------------------------------------------------
#Cluster samples based on distance 
sampleTree <- hclust(dist(expression.data), method = "average")
#Set the graphical parameters
par(cex = 0.6);
par(mar = c(0,4,2,0))
#Plot the cluster dendrogram
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,
     cex.axis = 1.5, cex.main = 2)

# 2.c Pick a soft threshold-----------------------------------------------------
# WGCNA, as the name implies, is a tool primarily intended for analyzing
# weighted networks. In a weighted network, you don't decide which nodes are 
# connected and which are not - all nodes are in principle connected, but the 
# strength varies (by convention) between 0 and 1. Soft thresholding really 
# means suppressing low correlations in a continuous ("soft") manner rather than
# the discontinuous ("hard") thresholding used in constructing unweighted 
# networks.
powers <- c(1:25)

sft = pickSoftThreshold(
  expression.data, # <= Input data
  #blockSize = 30,
  powerVector = powers,
  corFnc = bicor,
  verbose = 5
)

plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3]) * sft$fitIndices[,2],
     xlab="Soft Threshold (power)", ylab="Scale Free Topology Model Fit, signed R^2",
     type="n", main="Scale independence")
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3]) * sft$fitIndices[,2],
     labels=powers, cex=0.9, col="#691312")
abline(h=0.75, col="#003366")  # Aim for >0.75

picked_power = 15

# Check connectivity
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels= sft$fitIndices[,1],col="#691312")

# --> https://fuzzyatelin.github.io/bioanth-stats/module-F21-Group1/module-F21-Group1.html#Weighted_Gene_Correlation_Network_Analysis


# 2.d Creating the network------------------------------------------------------
temp_cor <- cor       
cor <- WGCNA::cor         # Force it to use WGCNA cor function (fix a namespace conflict issue)
netwk <- blockwiseModules(expression.data,                # <= input here
                          
                          # == Adjacency Function ==
                          power = picked_power,                # <= power here
                          networkType = "signed",
                          
                          # == Tree and Block Options ==
                          deepSplit = 2,
                          pamRespectsDendro = F,
                          # detectCutHeight = 0.75,
                          minModuleSize = 30,
                          maxBlockSize = 4000,
                          
                          # == Module Adjustments ==
                          reassignThreshold = 0,
                          mergeCutHeight = 0.30,
                          
                          # == TOM == Archive the run results in TOM file (saves time)
                          saveTOMs = T,
                          saveTOMFileBase = "ER",
                          
                          # == Output Options
                          numericLabels = T,
                          verbose = 3)

# extract the topological overlap matrix (TOM)----------------------------------
# We can either extract TOM from the previously generated network or calculate it manually
# Option 1; Extract TOM from previously generated network. Be aware that the command "blockwiseModules" sets diagonal of TOMsimilarity to zero
load("ER-block.1.RData")
class(TOM)

TOM_similarity <- as.matrix(TOM) 
rownames(TOM_similarity) <- colnames(expression.data)
colnames(TOM_similarity) <- colnames(expression.data)

# Option 2: Calculate TOM manually
#TOM_manual <-TOMsimilarityFromExpr(expression.data, power = picked_power, networkType = "signed")
#rownames(TOM_manual) <- colnames(expression.data)
#colnames(TOM_manual) <- colnames(expression.data)


# Plot the cluster dendrogram---------------------------------------------------
# Convert labels to colors
mergedColors = labels2colors(netwk$colors)#mergedColors is a character vector with entries like "turquoise", "blue", etc.
table(mergedColors)
# Plot the dendrogram and the module colors underneath
plotDendroAndColors(
  netwk$dendrograms[[1]],
  mergedColors[netwk$blockGenes[[1]]],
  "Module colors",
  dendroLabels = FALSE,
  hang = 0.03,
  addGuide = TRUE,
  guideHang = 0.05 )

module_df <- data.frame(
  gene_id = names(netwk$colors),
  colors = labels2colors(netwk$colors)
)

# Extract a table of genes and the assigned modules
#module_df is the mapping between each gene and its module assignment
module_df[1:5,]

write_delim(module_df,
            file = "gene_modules_filterByExpr.txt",
            delim = "\t")

# 3. ASSOCIATE MODULES WITH TRAITS ---------------------------------------------
#https://rpubs.com/natmurad/WGCNA

# 3.a Get Module Eigengenes per cluster ----------------------------------------
MEs0 <- moduleEigengenes(expression.data, mergedColors)$eigengenes

# Reorder modules so similar modules are next to each other
MEs0 <- orderMEs(MEs0)
module_order = names(MEs0) %>% gsub("ME","", .)

#clustering the eigengenes modules
#calculate the module dissimilarity eigengenes (without grey module)
MEs0selected <- MEs0[ , !(names(MEs0) %in% "MEgrey")]

selectedMEDiss = 1-cor(MEs0selected)

METree = hclust(as.dist(selectedMEDiss), method = "average")
plot(METree, main = "Clustering of module eigengenes",
     xlab = "", sub = "")

# 3.b Prepare metadata file-----------------------------------------------------
# Read metadata
metaData <- read.table("metadata_Gorniak-etal-2026.tsv", header = TRUE, stringsAsFactors = FALSE)
head(metaData)
samples <- rownames(expression.data)

# Binarize categorical variables in metadata
metaData_group.bin <- binarizeCategoricalColumns(metaData$group,
                                                 dropFirstLevelVsAll = FALSE,
                                                 includePairwise = FALSE,
                                                 includeLevelVsAll = TRUE,
                                                 minCount = 1) 

rownames(metaData_group.bin) <- samples
head(metaData_group.bin)

heatmap.data <- merge(MEs0, metaData_group.bin, by = "row.names")
head(heatmap.data)

heatmap.data <- heatmap.data %>%
  column_to_rownames(var = "Row.names")
head(heatmap.data)

metaData <- metaData %>%
  column_to_rownames(var = "sample")
head(metaData)

heatmap.data <- merge(heatmap.data, metaData[-c(1:3)], by = "row.names")
heatmap.data <- heatmap.data %>%
  column_to_rownames(var = "Row.names")
head(heatmap.data)

# Correlation analysis between categorical metadata (after binarization)
nGenes = ncol(expression.data)
nSamples = nrow(expression.data)
module.trait.correlation = cor(MEs0, metaData_group.bin, use = "p") #p for pearson 
#calculate the p-value associated with the correlation
module.trait.Pvalue = corPvalueStudent(module.trait.correlation, nSamples)

# Rename column names
colnames(heatmap.data)[colnames(heatmap.data) %in% c("data.Bas.vs.all", "data.Ferro.vs.all", "data.Gad.vs.all", "data.LaCl.vs.all", "data.LaOx.vs.all", "data.LaP.vs.all", "data.Lop.vs.all", "data.MinLaCl.vs.all", "data.Mon.vs.all", "data.Xen.vs.all")] <- c("Bastnaesite", "Ferrocerium", "Gadolinite", "LaCl3", "La2O3", "LaPO4", "Loparite", "Min_LaCl3", "Monazite", "Xenotime")
colnames(metaData_group.bin)[colnames(metaData_group.bin) %in% c("data.Bas.vs.all", "data.Ferro.vs.all", "data.Gad.vs.all", "data.LaCl.vs.all", "data.LaOx.vs.all", "data.LaP.vs.all", "data.Lop.vs.all", "data.MinLaCl.vs.all", "data.Mon.vs.all", "data.Xen.vs.all")] <- c("Bastnaesite", "Ferrocerium", "Gadolinite", "LaCl3", "La2O3", "LaPO4", "Loparite", "Min_LaCl3", "Monazite", "Xenotime")

CorLevelPlot(heatmap.data, 
             x = names(heatmap.data[11:23]),
             y = names(heatmap.data[1:10]),
             col = c("blue1", "skyblue", "white", "pink", "red"))

# 3.d Observe gene significances and module memberships of genes from selected modules-------------------------
#https://rpubs.com/natmurad/WGCNA

module.membership.measure <- cor(MEs0, expression.data, use = "p")
module.membership.measure.pvals <- corPvalueStudent(module.membership.measure, nSamples)

#names (colors) of the modules
modNames = substring(names(MEs0), 3)

names(module.membership.measure) = paste("MM", modNames, sep="")
names(module.membership.measure.pvals) = paste("p.MM", modNames, sep="")

write.table(as.data.frame(t(module.membership.measure)),
            file = "module_membership_measure_20250930.txt",
            sep = "\t")
write.table(as.data.frame(t(module.membership.measure.pvals)),
            file = "module_membership_pvals_20250930.txt",
            sep = "\t")

# Correlation between genes and selected traits
gene_significance <- data.frame(NA_col = rep(NA, 3878))
gene_significance_pval <- data.frame(NA_col = rep(NA, 3878))
rownames(gene_significance) <- colnames(expression.data)
rownames(gene_significance_pval) <- colnames(expression.data)
gene.signifList = list()
gene.signifList.pval = list()
#for(i in colnames(metaData_group.bin)) {
for(i in colnames(heatmap.data[-c(1:10)])) {
  #gene.signif <- cor(expression.data, metaData_group.bin[i], use = "p")
  gene.signif <- cor(expression.data, heatmap.data[-c(1:10)][i], use = "p")
  gene.signif.pvals <- corPvalueStudent(gene.signif, nSamples)
  colnames(gene.signif) <- i
  colnames(gene.signif.pvals) <- i
  gene.signifList[[i]] <- gene.signif
  gene.signifList.pval[[i]] <- gene.signif.pvals
}
gene_significance = as.data.frame(do.call(cbind, gene.signifList))
gene_significance_pval = as.data.frame(do.call(cbind, gene.signifList.pval))

names(gene_significance) = paste("GS.", names(gene_significance), sep ="")
names(gene_significance_pval) = paste("p.GS.", names(gene_significance_pval), 
                                      sep ="")

geneInfo_summary <- merge(gene_significance, gene_significance_pval, 
                          by = "row.names", sort = FALSE)
geneInfo_summary$module_color <- mergedColors

geneInfo_summary <- geneInfo_summary %>%
  column_to_rownames(var = "Row.names")
head(geneInfo_summary)

# combined data.frame containing gene significances and module memberships
geneInfo_summary <- merge(geneInfo_summary, t(module.membership.measure),
                  by = "row.names",
                  sort = FALSE)

write.csv(geneInfo_summary, file = "geneInfo_summary.csv")


# Plot gene significances and module memberships for selected modules
# Exemplary shown for module blue and trait Xenotime

module = "blue" #adjust
column = match(module, modNames)
moduleGenes = mergedColors==module


#add colors according to functions
gene_functions <- read.table("gene_function_list.txt", 
                             header = FALSE, sep = "\t", stringsAsFactors = FALSE)
colnames(gene_functions) <- c("gene", "annot")

# Build lookup vector: names = gene, values = annot
func_map <- setNames(gene_functions$annot, gene_functions$gene)

# Assign colors to each function
func_levels <- unique(gene_functions$annot)

func_colors <- c(
  "carbohydrate_metabolism" = "#34A0A4",
  "others"  = "#CED4DA",
  "lanthanome"  = "gold",
  "transport" = "#2C4268",
  "bacteriochlorophyll_biosynthesis" ="#7209B7",
  "C1_metabolism" = "#005F73",
  "motility_and_chemotaxis" = "#FF4800",
  "Fe_uptake_conversion_and_FeS_cluster_assembly" = "black",
  "methylolanthanin" = "gold",
  "N_metabolism" = "#A7C957",
  "PHA" = "#F72585",
  "PQQ" = "#ADE8F4",
  "Ribosome" ="#CD9777",
  "stress_associated" = "#6B705C",
  "S_metabolism" = "#780116",
  "Polysaccharide_synthesis_and_transport" = "#F28482",
  "hypotheticals" = "#CED4DA"
)

func_shapes <- c(
  "carbohydrate_metabolism" = "18",
  "others"  = "19",
  "lanthanome"  = "17",
  "transport" = "19",
  "bacteriochlorophyll_biosynthesis" ="15",
  "C1_metabolism" = "17",
  "motility_and_chemotaxis" = "19",
  "Fe_uptake_conversion_and_FeS_cluster_assembly" = "15",
  "methylolanthanin" = "18",
  "N_metabolism" = "19",
  "PHA" = "15",
  "PQQ" = "17",
  "Ribosome" ="15",
  "stress_associated" = "15",
  "S_metabolism" = "19",
  "Polysaccharide_synthesis_and_transport" = "18",
  "hypotheticals" = "22"
)

func_alpha <- c(
  "carbohydrate_metabolism" = "0.7",
  "others"  = "0.2",
  "lanthanome"  = "1",
  "transport" = "0.7",
  "bacteriochlorophyll_biosynthesis" ="0.7",
  "C1_metabolism" = "0.7",
  "motility_and_chemotaxis" = "0.7",
  "Fe_uptake_conversion_and_FeS_cluster_assembly" = "0.7",
  "methylolanthanin" = "0.8",
  "N_metabolism" = "1",
  "PHA" = "0.7",
  "PQQ" = "1",
  "Ribosome" = "0.7",
  "stress_associated" = "0.7",
  "S_metabolism" = "0.9",
  "Polysaccharide_synthesis_and_transport" = "0.9",
  "hypotheticals" = "0.5"
)

# Example: set custom plotting order (first = back, last = front)
priority_order <- c(
  "others", 
  "hypotheticals", 
  "carbohydrate_metabolism",
  "Ribosome",
  "transport",
  "N_metabolism",
  "stress_associated",
  "motility_and_chemotaxis",
  "S_metabolism",
  "Fe_uptake_conversion_and_FeS_cluster_assembly",
  "bacteriochlorophyll_biosynthesis",
  "PHA",
  "Polysaccharide_synthesis_and_transport",
  "C1_metabolism",
  "PQQ",
  "methylolanthanin",
  "lanthanome"   # <- these will end up on top
)

t_expression.data <- t(expression.data)

genes_in_module <- rownames(t_expression.data)[moduleGenes]
gene_funcs <- func_map[genes_in_module]  # lookup functional category for each gene

gene_cols  <- func_colors[gene_funcs]
gene_shapes <- as.numeric(func_shapes[gene_funcs])
gene_alpha <- as.numeric(func_alpha[gene_funcs])

plot_order <- order(match(gene_funcs, priority_order))
****
  
  
xvals <- abs(t(module.membership.measure)[moduleGenes, column])
yvals <- abs(gene_significance$GS.Xenotime[moduleGenes])# adjust

plot(xvals, yvals, type="n",
     xlab = paste("Module Membership in", module, "module"),
     ylab = "Gene significance for Gadolinite",
     main = "Module membership vs. gene significance",
     xlim = c(0, 1),     # set x-axis from 0 to 1
     ylim = c(0, max(yvals)))  # optional y-axis control

for (i in plot_order) {
  points(xvals[i], yvals[i],
         col = adjustcolor(gene_cols[i], alpha.f = gene_alpha[i]),
         pch = gene_shapes[i],
         cex = 2.5)  # dot size
}

# Compute correlation and p-value
corTest <- cor.test(xvals, yvals, use = "p")

r <- corTest$estimate
p <- corTest$p.value


# Add regression line (optional)
abline(lm(yvals ~ xvals), col = "black", lwd = 2)

# Add correlation text in top-left corner
text(
  x = 0.05, y = max(yvals) * 0.95,
  labels = paste0("r = ", round(r, 2), "\nP = ", signif(p, 2)),
  adj = 0
)

#add legend
plot.new()
plot.window(xlim = c(0, 1), ylim = c(0, 1))

# Add legend
legend(
  "center",
  legend = names(func_colors),
  col = func_colors,
  pch = as.numeric(func_shapes),
  pt.cex = 1.8,
  bty = "n",
  ncol = 1,
  title = "Functional category"
)




# 4. PLOT GENE EXPRESSION OF SELECTED MODULES ----------------------------------
#  Check the expression of modules of interest
d_log_cpm_WGCNA_AL1 <- cpm(d_AL1, normalized.lib.sizes=TRUE, log=TRUE, prior.count=0.25)

# pick out a few modules of interest here
modules_of_interest = c("turquoise" , "blue", "brown", "yellow", "green", "red",
                        "black", "pink", "magenta")

# Pull out list of genes in that module
submod = module_df %>%
  subset(colors %in% modules_of_interest)

row.names(module_df) = module_df$gene_id
subset_cpm = d_log_cpm_WGCNA_AL1[module_df$gene_id,]

subset_cpm = d_log_cpm_WGCNA_AL1[submod$gene_id,]

subset_cpm_df = data.frame(subset_cpm) %>%
  mutate(
    gene_id = row.names(.)
  ) %>%
  pivot_longer(-gene_id) %>%
  mutate(
    #module = module_df[gene_id,]$colors
    module = submod[gene_id,]$colors
  )

# Assume MEs is your module eigengene matrix (samples x modules)
MEs_scaled <- as.data.frame(scale(MEs0))  # scales each column (module) to mean=0, SD=1

MEs_scaled$sample <- rownames(MEs_scaled)

metaData <- as.data.frame(metaData)#LG
metaData$sample <- rownames(metaData)
metaData <- metaData[, c("sample", setdiff(names(metaData), "sample"))]

MEs_long <- MEs_scaled %>%
  pivot_longer(cols = starts_with("ME"), names_to = "Module", values_to = "Eigengene") %>%
  left_join(metaData, by = "sample")

# Order samples by mineral type and experiment
MEs_long$sample <- factor(MEs_long$sample,
                          levels = metaData %>%
                            arrange(experiment, mineral_type) %>%  # or `group`, `Timepoint`, etc.
                            pull(sample))

MEs_long_filtered <- MEs_long %>% filter(Module %in% c("MEturquoise" , "MEblue", 
                                                       "MEbrown", "MEyellow", 
                                                       "MEgreen", "MEred", 
                                                       "MEblack", "MEpink",
                                                       "MEmagenta"))

ggplot(MEs_long, aes(x = sample, y = Eigengene)) +
  geom_line(aes(group = Module), alpha = 0.2, color = "black") +  # optional: light connecting lines
  geom_point(aes(fill = mineral_type), color = "black", size = 3, stroke = 0.75, shape = 21) +
  facet_wrap(~ Module, scales = "free_y") +
  theme_classic() +
  scale_fill_manual(values = c(
    "carbonate" = "#7f0000ff",
    "chloride" = "black",
    "ferrocerium" = "#ffb717ff",
    "oxide" = "#145e82ff",
    "phosphate" ="#db5411ff",
    "silicate" = "#78d1e6ff"
  )) +
  labs(x = "Sample", y = "Module eigengene",
       title = "Module eigengene expression across samples") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

ggplot(MEs_long_filtered, aes(x = sample, y = Eigengene)) +
  geom_line(aes(group = Module), alpha = 0.2, color = "black") +  # optional: light connecting lines
  geom_point(aes(fill = mineral_type), color = "black", size = 3, stroke = 0.75, shape = 21) +
  facet_wrap(~ Module, scales = "free_y") +
  theme_classic() +
  scale_fill_manual(values = c(
    "carbonate" = "#7f0000ff",
    "chloride" = "black",
    "ferrocerium" = "#ffb717ff",
    "oxide" = "#145e82ff",
    "phosphate" ="#db5411ff",
    "silicate" = "#78d1e6ff"
  )) +
  labs(x = "Sample", y = "Module eigengene",
       title = "Module eigengene expression across samples") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

