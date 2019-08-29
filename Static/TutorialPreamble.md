# Statement of tutorial objectives

The aim of this tutorial is to demonstrate a workflow for discovering and characterising genomic structural variations (SV). This is achieved by mapping Nanopore long sequence reads to a reference genome and evaluating discordant mapping characteristics using the **`Sniffles`** software. This tutorial is based on the [Oxford Nanopore pipeline-structural-variation project](https://github.com/nanoporetech/pipeline-structural-variation) available from our Github pages and may be used to identify genomic insertion, duplication and deletion events.

The tutorial is packaged with example data from the [Genome in a Bottle](https://www.nist.gov/programs-projects/genome-bottle) project, and the workflow can be reproduced to address questions such as

* How many sequence reads map to the reference genome?
* What is the depth of sequence coverage across the genome, and across individual chromosomes?
* How many SVs can be identified, what is their size distribution, and  how frequent are the different types of SV?
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

* **`conda`** for management of bioinformatics software installations
* **`snakemake`** for managing the bioinformatics workflow
* **`minimap2`** for mapping sequence reads to reference genome
* **`sniffles`** for structural variation calling
* **`RSamtools`** and **`GenomicAlignments`**; R software for parsing BAM files
* **`IGV`** for visualising mapping characteristics at specific genomic regions
* **`Ribbon`** is a web-application for the presentation of structural variation data

## The computational requirements include: 

* Computer running Linux (Centos7, Ubuntu 18_10, Fedora 29, macOS)
* Multiple CPU cores are ideal; a 4 core CPU at minimum would be recommended 
* At least 24 Gb RAM - this is sufficient for mapping against the human genome and can report a full MinION flowcell worth of sequence data. The packaged dataset uses just human chromosome 4 and 8Gb RAM is sufficient
* At least 25 Gb spare disk space for analysis and data files; using a complete 30X human dataset would require approximately 500 Gb of available storage space.
* Runtime with provided example data - approximately 1 hour (assuming a recent 4 core processor and sufficient RAM)

\pagebreak

# Software installation

1. Most software dependencies are managed though **`conda`**. Install as described at  <br> [https://conda.io/docs/install/quick.html](https://conda.io/docs/install/quick.html). You will need to accept the license agreement during installation and we recommend that you allow the Conda installer to prepend its path to your `.bashrc` file when asked.
```
    wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
    bash Miniconda3-latest-Linux-x86_64.sh
    bash
```
2. Download Nanopore tutorial & example files into a folder named `ont_tutorial_sv`. This tutorial requires the **`git-lfs`** large file support capabilities, which should be installed through **`Conda`** first
```
    conda install -c conda-forge git-lfs
    git lfs install
    git clone https://github.com/nanoporetech/ont_tutorial_sv.git ont_tutorial_sv
```
3. Change your working directory into the new `ont_tutorial_sv` folder 
```
    cd ont_tutorial_sv
```
4. The Conda software dependencies are managed directly by the **`snakemake`** software that will be introduced in the following sections.


# Introduction

Structural variation is a type of genetic variation. Structural variation (SV) has traditionally been challenging to analyse but Oxford Nanopore long sequence reads provide new opportunities for the discovery and characterisation of SV. Longer sequence reads map faithfully to the reference genome and broken mapping segments can be used to discover insertions, deletions and duplications [reference to Sniffles manuscript].

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

As a best practice this tutorial will separate primary DNA sequence data (the base-called fastq files) from the **`Rmarkdown`** source and the genome reference data. The analysis results and figures will again be placed in a separate working directory. The required layout for the primary data is shown in the figure below. This minimal structure will be prepared over the next few sections of this tutorial. The DNA sequences and mapping BED file must be placed within a folder called **`RawData`** and the reference genome and annotation files must be placed in a folder named **`ReferenceData`**. All results will be placed in a folder named **`Analysis`** and different sub-folders will be created for different steps in the workflow. The **`Static`** folder contains accessory methods, texts, bibliography and graphics. The **`pipeline-structural-variation`** folder is a local copy of the [Oxford Nanopore pipeline-structural-variation ](https://github.com/nanoporetech/pipeline-structural-variation) github project.


![](Static/Images/FolderLayout.png) 

# Experimental setup

The first required step for running a structural variation discovery workflow is to define the experimental design. 

INSERT SOME TEXT ON WHAT CAN BE TUNED, REQUIRED FILE INPUTS AND RECOMMENDATIONS FOR DEPTH OF COVERAGE

The design for the tutorial is defined within a YAML format configuration file (**`config.yaml`**). The tutorial's file is shown in the figure below.


![](Static/Images/ExperimentalDesign.png) 




* **`pipeline`** identifies the workflow that this configuration file belongs to

SOME TEXT AND DESCRIPTION FOR WHAT THE INDIVIDUAL PARAMETERS ACTUALLY DO


\newpage

# Snakemake

This SV discovery tutorial uses **`snakemake`** (@snakemake2012). Snakemake is a workflow management system implemented in Python. The aim of the Snakemake tool is to enable reproducible and scalable data analyses. The workflow produced within this document should be portable between laptop resources, computer servers and other larger scale IT deployments. The Snakemake workflow additionally defines the sets of required software (and software versions where appropriate) and will automate the installation and deployment of this software through the **conda** package management system.

The **`snakemake`** workflow will call methods that include **`minimap2`** @minimap22018, **`samtools`** @samtools2009, **`Sniffles`** @sniffles2018 and **`pipeline-structural-variation`**. The defined workflow is summarised in the figure below. 

<!-- [^Comment]: # (DAG can be produced with the command = 
snakemake --snakefile ./pipeline-structural-variation/Snakefile --configfile ./config.yaml  --forceall --dag call | grep -v target.bed | grep -v Working | dot -Tpng > structural_variation_workflow.png)
--> 


![](Static/Images/structural_variation_workflow.png) 

The precise commands within the **`Snakefile`** include

* index the reference genome sequence using **``minimap2`**
* map DNA sequence reads against the reference genome index using **`minimap2`** - prepare a sorted and indexed BAM file of sequence results
* use the **`sniffles`** software to identify SVs
* define BED coordinates from the mapped genome and calculate depth-of-coverage across the mapped regions
* filter the called SVs using both mapping coordinate information and depth-of-coverage information to produce a high-confidence set of SVs.
* **`GIAB`** truthset data will be downloaded; truthset variants that intersect with the experimental data will be identified
* **`Truvari`** will be used to assess performance of identified variants against the truthset


# Run the snakemake workflow file


<!-- [^Comment]: # (Include --stats snakemake.stats for logging runtime for the different steps) -->

\fontsize{8}{12}
```
# snakemake is used to run the workflow

snakemake --snakefile ./pipeline-structural-variation/Snakefile --configfile ./config.yaml --jobs 8 eval
```
\fontsize{10}{14}

The **`eval`** task specifies that the analysis is performed and the results are compared to the GIAB truthset. If you are not using the dataset packaged with the tutorial (or other GM24385 dataset) then the workflow can be run with the **`call`** suffix instead. This will map reads and prepare a high-confidence set of SVs; the benchmarking will not be performed.

\pagebreak
 

# Prepare the analysis report

The tutorial document can also be prepared from the Linux command line with the following command. 

\fontsize{8}{12}
```
R --slave -e 'rmarkdown::render("ont_tutorial_sv.Rmd", "html_document")'
```
\fontsize{10}{14}

This will prepare an HTML format document called **`ont_tutorial_sv.html`** that, depending of modifications to the **`config.yaml`** file will closely resemble this tutorial vignette.
