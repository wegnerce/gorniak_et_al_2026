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

RNASeq DATA PROCESSING
-------------------------------------------------

**1._Needed software_**

All listed tools can be easily installed and set up using [`conda`](https://docs.conda.io/en/latest/miniconda.html) and [`bioconda`](https://bioconda.github.io/).

- [_fastQC_](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) (v0.11.9)
- [_bbmap_ part of the BBtools suite](https://jgi.doe.gov/data-and-tools/software-tools/bbtools/) (v38.26)
- [_SortMeRNA_](https://bioinfo.lifl.fr/RNA/sortmerna/) (v2.1b)
- [_samtools_](http://www.htslib.org/) (v1.3.1)
- [_featureCounts_ part of the Subread package](http://subread.sourceforge.net/) (v1.6.3)

The described steps are carried out on the `commandline` using basic `BASH` scripting.

For the sake of reproducibility, we strongly recommend to use workflow management systems such as [_snakemake_](https://snakemake.github.io/) or [_nextflow_](https://www.nextflow.io/).

We provide an automated, _snakemake_-based workflow for RNAseq data processing on [_github_](https://github.com/wegnerce/smk_rnaseq).

We also provide a _snakemake_ profile to be used with HPC environments that rely on [_SLURM_](https://slurm.schedmd.com/documentation.html) for job scheduling on [_github_](https://github.com/wegnerce/smk_slurm_profile).

Please have a look at the repositories for details. Below you find an outline of the steps that are part of our data processing workflow, and which are automated by the mentioned _snakemake_ workflow. 

**2._Quality control,quality assessment and trimming_**

We assume that all _raw_ data files are located in a folder named `00_RAW`. Dataset names are listed in a `.txt` file called datasets.txt (one dataset per line). Forward and reverse read files end with ".1.fastq.gz" and ".2.fastq.gz", respectively. Quality control,quality assessment and trimming is done as follows using _bbduk_ and _fastQC_:

    #!/bin/bash
    for dataset in `cat datasets.txt`
        do
            fastqc -o ./00_RAW/ ./00_RAW/"$dataset".1.fastq.gz
            fastqc -o ./00_RAW/ ./00_RAW/"$dataset".2.fastq.gz
            bbduk.sh -Xmx8g in1=00_RAW/"$dataset".1.fastq.gz in2=00_RAW/"$dataset".2.fastq.gz out1=01_TRIMMED/"$dataset"_1_trimmed.fastq out2=01_TRIMMED/"$dataset"_2_trimmed.fastq stats="$dataset"_stats_QC.txt minlen=75 qtrim=rl trimq=20 ktrim=r k=25 mink=11 ref=adapters.fa hdist=1
            fastqc -o ./01_TRIMMED/ ./01_TRIMMED/"$dataset"_1_trimmed.fastq
            fastqc -o ./01_TRIMMED/ ./01_TRIMMED/"$dataset"_2_trimmed.fastq
        done

The command above does three things: 
1. it generates QC reports of the raw data, 
2. performs quality control + adapter trimming, 
3. and generates QC reports of the trimmed data. 

Trimmed read data and corresponding QC reports are stored in a folder called `01_TRIMMED`.

**3._Filtering of rRNA-derived sequences_**

For downstream analyses, we need to remove rRNA-derived sequences. We use  _SortMeRNA_ for doing this. Trimmed sequences are supposed to be in a folder called `01_TRIMMED`. We can re-use the aforementioned `datasets.txt` file.

    #!/bin/bash
    for dataset in `awk '{print $1}' datasets.txt`
        do merge-paired-reads.sh 01_TRIMMED/"$dataset"_R1_trimmed.fastq 01_TRIMMED/"$dataset"_R1_trimmed.fastq 02_rRNA_FILTERED/"$dataset"_merged.fastq
        sortmerna --ref rRNA_databases/silva-bac-16s-id90.fasta,index/silva-bac-16s-db:rRNA_databases/silva-bac-23s-id98.fasta,index/silva-bac-23s-db:rRNA_databases/silva-arc-16s-id95.fasta,index/silva-arc-16s-db:rRNA_databases/silva-arc-23s-id98.fasta,index/silva-arc-23s-db --reads 02_rRNA_FILTERED/"$dataset"_merged.fastq --other 02_rRNA_FILTERED/"$dataset"_non_rRNA --aligned 02_rRNA_FILTERED/"$dataset"_rRNA --fastx --best 5 --paired_in -a 4 -v --blast 1
        unmerge-paired-reads.sh 02_rRNA_FILTERED/"$dataset"_rRNA 02_rRNA_FILTERED/"$dataset"_rRNA_f.fastq 02_rRNA_FILTERED/"$dataset"_rRNA_r.fastq
        rm 02_rRNA_FILTERED/"$dataset"_merged.fastq
        rm 02_rRNA_FILTERED/"$dataset"_rRNA.fastq
        done

The database paths will vary dependent on how and where you installed _SortMeRNA_. After filtering, filtered sequences are kept in a folder called `02_rRNA_FILTERED`.

_SortMeRNA_ first merges read pairs, matches them against the specified databases, and then splits the data into rRNA- and non rRNA-derived sequences. The latter is kept for further processing.

**4._Read mapping_**

Trimmed and filtered RNA-Seq sequences are mapped onto the genome of _Beijerinckiaceae_ bacterium RH AL1 using _bbmap_. The genome of strain RH AL1 is available from this repository, as well as an annotation table containing annotation data from NCBI RefSeq and KEGG GENES.

[Genome .fna](https://github.com/wegnerce/gorniak_et_al_2026/blob/main/data_files/RHAL1_chromosome_plasmid.fa)

[Genome annotation](https://github.com/wegnerce/gorniak_et_al_2026/blob/main/data_files/RHAL1_annotation_table.csv)


    #!/bin/bash
    for dataset in `cat datasets.txt`
        do
            bbmap.sh slow k=11 in=./02_FILTERED/"$dataset"_non_rRNA_1.fastq in2=./02_FILTERED/"$dataset"_non_rRNA_2.fastq ref=RHAL1_chromosome_plasmid.fa out=./03_MAPPED/"$dataset"_mapped_AL1.bam statsfile="$dataset"_mapping_stats.txt
            samtools sort ./03_MAPPED/"$dataset"_mapped_AL1.bam -o ./03_MAPPED/"$dataset"_mapped_AL1_sorted.bam
            samtools index ./03_MAPPED/"$dataset"_mapped_AL1_sorted.bam
            rm ./03_MAPPED/"$dataset"_mapped_AL1.bam
        done

The resulting indexed, and sorted .bam files are kept in a folder called `03_MAPPED`.

**5._Generating read count data for downstream differential gene expression analysis_**

Read count data are generated with _featureCounts_. _featureCounts_ uses .gtf or .saf files as reference. A .saf file for the _Beijerinckiaceae_ bacterium RH AL1 genome can be found in this repository.

[Genome .saf](https://github.com/wegnerce/gorniak_et_al_2026/blob/main/data_files/RHAL1_chromosome_plasmid.saf)

    #!/bin/bash
    for dataset in `cat datasets.txt`
        do
            featureCounts -a "RHAL1_chromosome_plasmid.saf -o readcounts_featureCounts_AL1.txt -F SAF -t GeneID 03_MAPPED/*sorted.bam
        done

The resulting `readcounts_featureCounts_AL1.txt` file is the starting point for downstream differential gene expression analysis.

[Per gene read count data](https://github.com/wegnerce/gorniak_et_al_2026/blob/main/data_files/readcounts_featureCounts.txt)
