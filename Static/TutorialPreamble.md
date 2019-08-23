# Statement of tutorial objectives

The aim of this tutorial is to demonstrate a workflow for discovering and characterising genomic structural variations. This is achieved by mapping Nanopore long sequence reads to a reference genome and evaluating discordant mapping coordinates using the **`Sniffles`** software. This tutorial is based on the [Oxford Nanopore pipeline-structural-variation project](https://github.com/nanoporetech/pipeline-structural-variation) available from out Github pages and may be used to identify large insertion and deletion events.

The tutorial is packaged with example data from the [Genome in a Bottle](https://www.nist.gov/programs-projects/genome-bottle), and the workflow can be reproduced to address questions such as

* How many sequence reads map to the reference genome?
* What is the depth of coverage across the genome
* How many structural variations can be identified, and what how frequent are the different types
* What is the precision and recall of structural variants from a reference truthset?

Editing the workflow's configuration file, **`config.yml`**, will allow the analyses to be run using different DNA sequence collections, reference genomes, and with different parameters to tune the structural variation discovery workflow.

## What will this tutorial produce?

* A rich HTML format report containing summary statistics and figures highlighting performance of the enrichment protocol
* **`Coordinates`** and instructions for reviewing candidate genomic regions with **`IGV`** the Integrated Genomics Viewer


## Methods utilised include: 

* **`conda`** for management of bioinformatics software installations
* **`snakemake`** for managing the bioinformatics workflow
* **`minimap2`** for mapping sequence reads to reference genome
* **`sniffles`** for structural variation calling
* **`RSamtools`** and **`GenomicAlignments`**; R software for parsing BAM files
* **`IGV`** for visualising mapping characteristics at specific genomic regions

## The computational requirements include: 

* Computer running Linux (Centos7, Ubuntu 18_10, Fedora 29, macOS)
* Multiple CPU cores are ideal; a 4 core CPU at minimum would be recommended 
* At least 16 Gb RAM - this is sufficient for mapping against the human genome and can report a full MinION flowcell worth of sequence data. The packaged dataset uses just human chromosome 4 and 8Gb RAM is sufficient
* At least 15 Gb spare disk space for analysis and indices
* Runtime with provided example data - approximately 45 minutes

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
4. Install Conda software dependencies
```
    conda env create --name ont_tutorial_sv --file environment.yaml
```
5. Initialise the Conda environment 
```
    source activate ont_tutorial_sv
    R CMD javareconf > /dev/null 2>&1
```


\pagebreak

# Introduction

Structural variation is a key variety of genetic variation. Structural variation has traditionally been challenging to analyse but Nanopore long sequence reads provide new opportunities for the discovery and characterisation of structural variation.


There are five goals for this tutorial:

* To introduce a literate framework for analysing structural variation from Oxford Nanopore DNA sequence data
* To utilise best data-management practices
* To map sequence reads to the reference genome and to discover structural varations using optimised SV discovery workflow
* To annotate insertion and deletion events for associated gene and repeat element content
* To benchmark the performance of SV discovery through assessment of accuracy and precision


# Getting started and best practices

This tutorial requires a computer workstation running a Linux operating system. The workflow described has been tested using **`Fedora 29`**, **`Centos 7`** and **`Ubuntu 18_04`**. This tutorial has been prepared in the **`Rmarkdown`** file format. This utilises *markdown* (an easy-to-write plain text format as used in many Wiki systems) - see @R-rmarkdown for more information about **`rmarkdown`**. The document template contains chunks of embedded **`R code`** that are dynamically executed during the report preparation. 

The described analytical workflow makes extensive use of the **`conda`** package management and the **`snakemake`** workflow software. These software packages and the functionality of **`Rmarkdown`** provide the source for a rich, reproducible and extensible tutorial document.

The workflow contained within this Tutorial performs an authentic bioinformatics analysis and using the human chromosome 4 as a reference sequence. There are some considerations in terms of memory and processor requirement. Indexing the whole human genome for sequence read mapping using **`minimap2`** will use at least **`18 Gb`** of memory. The minimal recommended hardware setup for this tutorial is a 4 threaded computer with at least 8 Gb of RAM and 10 Gb of storage space. 

There are few dependencies that need to be installed at the system level prior to running the tutorial. The **`conda`** package management software will coordinate the installation of the required bioinformatics software and their dependencies in user space - this is dependent on a robust internet connection.

As a best practice this tutorial will separate primary DNA sequence data (the base-called fastq files) from the **`Rmarkdown`** source and the genome reference data. The analysis results and figures will again be placed in a separate working directory. The required layout for the primary data is shown in the figure below. This minimal structure will be prepared over the next few sections of this tutorial. The DNA sequences and mapping BED file must be placed within a folder called **`RawData`** and the reference genome and annotation files must be placed in a folder named **`ReferenceData`**. All results will be placed in a folder named **`Analysis`** and different sub-folders will be created for different steps in the workflow. The **`Static`** folder contains accessory methods, texts, bibliography and graphics.


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

The **`snakemake`** workflow will call methods that include **`minimap2`** @minimap22018, **`samtools`** @samtools2009 and **`Sniffles`** @sniffles2018. The defined workflow is summarised in the figure below. 

<!-- [^Comment]: # (DAG can be produced with the command = snakemake --snakefile ./pipeline-structural-variation/Snakefile --configfile --> ./config.yaml  --forceall --dag call | grep -v target.bed | grep -v Working | dot -Tpng > structural_variation_workflow.png)

![](Static/Images/structural_variation_workflow.png) 

The precise commands within the **`Snakefile`** include

* download the specified reference genome
* use **`minimap2`** to index the reference genome
* map DNA sequence reads against the reference genome index using **`minimap2`**
* convert **`minimap2`** output (**`SAM`**) into a sorted **`BAM`** format using **`samtools`** and filter out the unmapped reads


WHAT OTHER STEPS ARE INCLUDED IN THE SNAKEFILE AND SHOULD BE INCLUDED HERE


# Run the snakemake workflow file


<!-- [^Comment]: # (Include --stats snakemake.stats for logging runtime for the different steps) -->

\fontsize{8}{12}
```
# snakemake is used to run the workflow
# --snakefile parameter is used to locate the Snakefile workflow
# --configfile is used to override the Snakefile default parameter file
# --jobs specifies the number of threads to use for the workflow (e.g. use 4 if you have 4 core PC)
# --stats collects information on the snakemake run that will be included later on in this report
# `call` is a snakemake target that specifies the sniffles workflow

snakemake --snakefile ./pipeline-structural-variation/Snakefile --configfile ./config.yaml --jobs 8 call
```
\fontsize{10}{14}


\pagebreak
 

# Prepare the analysis report

The **`Rmarkdown`** script can be run using the **`knit`** dialog in the **`Rstudio`** software. 

The document can also be rendered from the command line with the following command. This command is  run automatically during the Snakemake workflow.

\fontsize{8}{12}
```
R --slave -e 'rmarkdown::render("ont_tutorial_cas9.Rmd", "html_document")'
```
\fontsize{10}{14}

