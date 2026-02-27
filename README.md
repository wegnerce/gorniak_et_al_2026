# Gorniak et al., 2026 |  Source-independent enrichment of light lanthanides: microbial mobilization, selective uptake, and intracellular storage

## doi:

Linda Gorniak<sup>1,2</sup>   
Sophie M. Gutenthaler-Tietze<sup>3,4</sup>     
Alina Lobe<sup>4</sup>     
Lena J. Daumann<sup>3,4</sup>     
Robin Steudtner<sup>5</sup>     
Thorsten Schäfer<sup>2</sup>     
Frank Steiniger<sup>6</sup>     
Martin Westermann<sup>6</sup>     
Kirsten Küsel<sup>1,7,8</sup>     
Carl-Eric Wegner<sup>1,3</sup>     

<sup>1</sup>Institute of Biodiversity, Aquatic Geomicrobiology, Friedrich Schiller University, Jena, Germany  
<sup>2</sup>Applied Geology, Institute of Geosciences, Friedrich Schiller University, Burgweg 11, 07749 Jena, Germany    
<sup>3</sup>Bioinorganic Chemistry, Heinrich Heine University, Universitätsstr. 1, 40225 Düsseldorf, Germany    
<sup>4</sup>Department of Chemistry, Ludwig Maximilian University, 80539 München, Germany  
<sup>5</sup>Institute of Resource Ecology, Helmholtz-Zentrum Dresden-Rossendorf e.V, Bautzner Landstraße 400, 01328 Dresden, Germany    
<sup>6</sup>Electron Microscopy Center, Jena University Hospital, Ziegelmühlenweg 1, 07743 Jena, Germany  
<sup>7</sup>Cluster of Excellence Balance of the Microverse, Friedrich Schiller University Jena, Jena, Germany   
<sup>8</sup>German Center for Integrative Biodiversity Research (iDiv) Halle-Jena-Leipzig, Germany

Correspondence should be addressed to: [**Carl-Eric Wegner**][1]

------------------

This repository provides the following resources:

### Data files
* _[Read count data per gene][2]_: The carried out steps to obtain the read count data, starting from raw data are explained here: [`RNA-Seq data processing`][5] 
* _[Mapping statistics][3]_: The given statistics are based on mapping processed reads onto the  genome of _Beijerinckiaceae_ bacterium RH AL1 using _bbmap_ and generating mapping statistics with _featureCounts_. 
* _[log2CPM data for all coding genes][4]_: The provided log2CPM (counts per million) data are the basis for differential gene expression and related analyses.
* _[Oerview of differentially expressed genes][8]_: A spreadsheet with a combined overview of differentially expressed genes for comparisons between the individual experimental conditions (= different lanthanide sources).

### Walkthroughs and scripts for data processing and analysis

* _[RNASeq data processing][7]_: Here we deccribe how datasets have been pre-processed prior to differential gene expression analysis.
* _[Differential gene expression analysis][5]_: This R script allows you to reproduce the overviews of differentially expressed genes, which we used for subsequent in-depth analysis of transcriptome changes in response to the different provided lanthanide sources.
* _[Network analysis][6]_: The code in this script was used to identify modules of related genes based on weighted gene correlation network analysis.

[1]: http://mailto:%20carl-eric.wegner@hhu.de
[2]: https://github.com/wegnerce/gorniak_et_al_2026/blob/main/data_files/readcounts_featureCounts.txt
[3]: https://github.com/wegnerce/gorniak_et_al_2026/blob/main/data_files/readcounts_featureCounts.txt.summary
[4]: https://github.com/wegnerce/gorniak_et_al_2026/blob/main/data_files/log_cpm.csv
[5]: https://github.com/wegnerce/gorniak_et_al_2026/blob/main/scripts/DGEA_Gorniak_et_al_2026.R
[6]: https://github.com/wegnerce/gorniak_et_al_2026/blob/main/scripts/WGCNA_Gorniak_et_al_2026.R
[7]: https://github.com/wegnerce/gorniak_et_al_2026/blob/main/data_files/overview_differential_gene_expression_analysis.xlsx
