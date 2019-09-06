# Statement of tutorial objectives

The aim of this tutorial is to demonstrate a workflow for discovering and characterising genomic structural variations (SV). This is achieved by mapping Nanopore long sequence reads to a reference genome and evaluating discordant mapping characteristics using the **`Sniffles`** software. This tutorial is based on the Oxford Nanopore [pipeline-structural-variation project](https://github.com/nanoporetech/pipeline-structural-variation) and is available from the linked Github pages. The pipeline may be used to identify genomic insertion, duplication and deletion events.

The tutorial is packaged with example data from the [Genome in a Bottle](https://www.nist.gov/programs-projects/genome-bottle) project, and the workflow can be reproduced to address questions such as

* How many sequence reads map to the reference genome?
* What is the depth of sequence coverage across the genome?
* How many SVs can be identified, what is their size distribution, and how frequent are the different types of SV?
* How do the SVs called compare to a reference truthset when accuracy and precision of predictions are considered?

The tutorial workflow uses a  configuration file, **`config.yaml`**, that specifies the DNA sequences to use, the reference genome and analysis characteristics. This file can be modified to run the workflow using different DNA sequence collections, alternative reference genomes, and with parameters tuned to the needs of SV discovery. 

## What will this tutorial produce?

* A rich HTML format report containing summary statistics and figures highlighting the performance of the long sequence read mapping and SV detection analysis
* **`VCF`** format file of total SVs discovered and filtered files of just the high-confidence structural variants
* **`Excel`** format spreadsheet file identifying the SVs and genes that overlap within the genes' protein coding space
* **`Excel`** format spreadsheet file identifying the SVs and annotated repeat sequences that overlap
* **`IGV Coordinates`** and instructions for reviewing candidate genomic regions with **`IGV`**, the Integrated Genomics Viewer
* **`Instructions for Ribbon`**, an interactive web-application that provides functionality for visualising structural variants

## Methods utilised include: 

A number of different bioinformatics software packages will be used to analyse the provided DNA sequence data. The workflow used in the tutorial applies the [pipeline-structural-variation](https://github.com/nanoporetech/pipeline-structural-variation). The key software packages used include

* **`conda`** for management of bioinformatics software installations
* **`snakemake`** for managing the bioinformatics workflow
* **`git`** packages for downloading the tutorial from Github repository
* **`git-lfs`** is required to download the sequence and metadata files provided with the tutorial
* **`minimap2`** for mapping sequence reads to reference genome
* **`samtools`** is used to handle the `sam` and `bam` format mapping data
* **`sniffles`** for structural variation calling
* **`R`** is a statistical analysis software and is used for the analysis and reporting of the sequence summary data
* **`RSamtools`** and **`GenomicAlignments`**; R software for parsing BAM files
* **`IGV`** for visualising mapping characteristics at specific genomic regions
* **`Ribbon`** is a web-application for the presentation of structural variation data

## The computational requirements include: 

* Computer running Linux (Centos7, Ubuntu 18_10, Fedora 29)
* Multiple CPU cores are ideal; a 4 core CPU at minimum would be recommended 
* At least 24 Gb RAM - this is sufficient for mapping against the human genome and can report a full MinION flowcell worth of sequence data. The packaged dataset uses a region of human chromosome 20 and 8Gb RAM is sufficient
* At least 25 Gb spare disk space for analysis and data files; using a complete 30X human dataset would require approximately 500 Gb of available storage space.
* Runtime with provided example data - approximately 1 hour (assuming a recent 4 core processor with sufficient RAM)



# Software installation

1. Most software dependecies are managed using **`conda`**. Please install as described at  <br> [https://conda.io/projects/conda/en/latest/user-guide/install/index.html](https://conda.io/projects/conda/en/latest/user-guide/install/index.html). You will need to accept the license agreement during installation and we recommend that you allow the Conda installer to prepend its path to your `.bashrc` file when asked.
```
    wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
    bash Miniconda3-latest-Linux-x86_64.sh
    bash
```
2. Download Nanopore tutorial & example files into a folder named `ont_tutorial_sv`. This tutorial requires the **`git-lfs`** large file support capabilities, which should be installed through **`Conda`** first
```
    conda install -c conda-forge git-lfs
    git lfs install
    git clone --recursive https://github.com/nanoporetech/ont_tutorial_sv.git ont_tutorial_sv
```
3. Change your working directory into the new `ont_tutorial_sv` folder 
```
    cd ont_tutorial_sv
```
4. Install Conda software dependencies
```
    conda env create --name ont_tutorial_sv --file environment.yaml
```
5. Initialise the Conda environment
```
    conda activate ont_tutorial_sv
```

# Introduction

Structural variation (SV) is a type of genetic variation. SV has traditionally been challenging to analyse but Oxford Nanopore long sequence reads provide new opportunities for the discovery and characterisation of SV. Longer sequence reads map faithfully to the reference genome and broken mapping segments can be used to discover insertions, deletions and duplications (@sniffles2018).

This tutorial demonstrates a workflow for the analysis and exploration of SVs. The tutorial is provided with publicly available Nanopore sequence data from human chromosome 4 and a human chromosome 4 sequence reference suitable for the benchmarking of SV discovery.

There are five goals for this tutorial:

* To introduce a literate framework for analysing structural variation from Oxford Nanopore DNA sequence data
* To introduce and utilise best data-management practices to ensure reproducibility
* To map sequence reads to the reference genome and to discover structural varations using **`pipeline-structural-variation`**, an SV discovery workflow optimised for Oxford Nanopore DNA sequence data
* To annotate genomic insertion, deletion and duplications against the human reference genome to identify the variants that overlap genes and annotated genome repeat units
* To benchmark the performance of SV discovery through assessment of precision and recall using a structural variation truthset


# Getting started and best practices

This tutorial requires a computer workstation running a Linux operating system. The workflow described has been tested using **`Fedora 29`**, **`Centos 7`** and **`Ubuntu 18_04`**. This tutorial has been prepared in the **`Rmarkdown`** file format. This utilises *markdown* (an easy-to-write plain text format as used in many Wiki systems) - see @R-rmarkdown for more information about **`rmarkdown`**. The document template contains chunks of embedded **`R code`** that are dynamically executed during the report preparation. 

The described analytical workflow makes extensive use of the **`conda`** package management and the **`snakemake`** workflow software. These software packages and the functionality of **`Rmarkdown`** provide the source for a rich, reproducible and extensible tutorial document.

The workflow contained within this Tutorial performs an authentic bioinformatics analysis and uses human chromosome 4 as a reference sequence. There are some considerations in terms of memory and processor requirement. Indexing the whole human genome for sequence read mapping using **`minimap2`** will use at least **`18 Gb`** of memory. The minimal recommended hardware setup for this tutorial is a 4 threaded computer with at least 8 Gb of RAM and 25 Gb of storage space. 

There are few dependencies that need to be installed at the system level prior to running the tutorial. The **`conda`** package management software will coordinate the installation of the required bioinformatics software and their dependencies in user space - this is dependent on a robust internet connection.

As a best practice this tutorial will separate primary DNA sequence data (the base-called fastq files) from the **`Rmarkdown`** source and the genome reference data. The analysis results and figures will again be placed in a separate working directory. The required layout for the primary data is shown in the figure below.

![](Static/Images/SVFolderLayout.png) 

This minimal structure will be prepared over the next few sections of this tutorial. The fastq format DNA sequences must be placed within a folder called **`RawData`** and the reference genome must be placed in the folder named **`ReferenceData`**. All results will be placed in the folder named **`Analysis`** and different sub-folders will be created for different steps in the workflow. The **`Static`** folder contains accessory methods, texts, bibliography and graphics. The **`pipeline-structural-variation`** folder is a local copy of the [pipeline-structural-variation ](https://github.com/nanoporetech/pipeline-structural-variation) github project. The **`...`** means that there are other files in the folder not shown in the figure.


# Experimental setup

The first required step for running a structural variation discovery workflow is to define the experimental design. The design is configured within a YAML format file (**`config.yaml`**). The YAML file for the tutorial is shown in the figure below.

![](Static/Images/SVExperimentalDesign.png) 

There are several critical parameters that should be appropriately specified

* **`input_fastq`** specifies the fastq sequence file containing sequence reads that should be used for SV discovery. This file must be fastq format - the **`pipeline-structural-variation`** workflow can use fastq files that have been compressed using **`gzip`**.
* **`reference_fasta`** specifies the path to the **fasta** format reference genome to use. This sequence may be compressed using **`gzip`**
* **`sample_name`** This is a name for the analysis - multiple runs can be stored under the **`Analysis`** folder
* **`biocGenome`** is used to instruct the R/Bioconductor software which version of a genome is to be used for the analysis
* **`GeneIdMappings`** is the name of a Bioconductor package that provides the mappings of gene names for the reference genome used
* **`GenomeAnnotation`** is the name of a Bioconductor package containing the gene annotation for the reference genome used
* **`RepeatDB`** is a URL link to an appropriate **`RepeatMasker`** database of genomic repeats

Other parameters may be used to tune the workflow characteristics and reported content

* **`min_sv_length`** - The length of the shortest SVs to report
* **`max_sv_length`** - The length of the longest SVs to report
* **`min_read_length`** - Minimum read length - shorter reads will be discarded during the analysis
* **`min_read_mapping_quality`** - Minimum mapping quality. Reads with lower mapping qualities will be discarded
* **`min_read_support`** - Minimum read support required to call a SV (**auto** for auto-detect)
* **`tutorialText`** is used to define whether the full tutorial text is included in the report; set this to **False** for a more result focused report


\newpage

# Snakemake

This SV discovery tutorial uses **`snakemake`** (@snakemake2012). Snakemake is a workflow management system implemented in Python. The aim of the Snakemake tool is to enable reproducible and scalable data analyses. The workflow produced within this document should be portable between laptop resources, computer servers and other larger scale IT deployments. The Snakemake workflow additionally defines the sets of required software (and software versions where appropriate) and will automate the installation and deployment of this software through the **conda** package management system.

The **`snakemake`** workflow will call methods that include **`minimap2`** @minimap22018, **`samtools`** @samtools2009, **`Sniffles`** @sniffles2018 and **`pipeline-structural-variation`**. The defined workflow is summarised in the figure below. 

<!-- [^Comment]: # (DAG can be produced with the command = 
snakemake --snakefile ./pipeline-structural-variation/Snakefile --configfile ./config.yaml  --forceall --dag eval | grep -v target.bed | grep -v Working | dot -Tpng > structural_variation_workflow.png)
--> 

![](Static/Images/structural_variation_workflow.png) 

The commands within the **`Snakefile`** workflow illustrated above include

* index the reference genome sequence using **`minimap2`**
* map DNA sequence reads against the reference genome index using **`minimap2`** - prepare a sorted and indexed BAM file of sequence results
* use **`sniffles`** software to identify SVs
* define BED coordinates from the mapped genome and calculate depth-of-coverage across the mapped regions
* filter the called SVs using both mapping coordinate information and depth-of-coverage information to produce a high-confidence set of SVs.
* **`GIAB`** truthset data is downloaded; truthset variants that intersect with the experimental data will be identified
* **`Truvari`** script will be downloaded and used to assess performance of identified variants against the truthset


# Run the snakemake workflow file


<!-- [^Comment]: # (Include --stats snakemake.stats for logging runtime for the different steps) -->

\fontsize{8}{12}
```
# snakemake is used to run the workflow

snakemake --snakefile ./pipeline-structural-variation/Snakefile --configfile ./config.yaml --jobs 4 eval
```
\fontsize{10}{14}

The **`--jobs 4`** instructs the Snakemake workflow to use 4 CPU cores. Please adjust this value to best match the computer resources available. 

The **`eval`** task specifies that the analysis is performed and the results are compared to the GIAB truthset. If you are not using the dataset packaged with the tutorial (or other GM24385 dataset) then the workflow can be run with the **`call`** suffix instead. This will map reads and prepare a high-confidence set of SVs; the benchmarking will not be performed.

\pagebreak
 

# Prepare the analysis report

The tutorial document can also be prepared from the Linux command line with the following command. 

\fontsize{8}{12}
```
R --slave -e 'rmarkdown::render("ont_tutorial_sv.Rmd", "html_document")'
```
\fontsize{10}{14}

This will prepare an HTML format document called **`ont_tutorial_sv.html`**.
