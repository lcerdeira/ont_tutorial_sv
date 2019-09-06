![.](Static/Images/ONT_logo.png "Oxford Nanopore Technologies")

******************

# 1. Introduction:


### Overview:

The **Structural variation analysis using Nanopore long DNA sequence reads** tutorial demonstrates a workflow for discovering and characterising genomic structural variations (SV). This is achieved by mapping Nanopore long sequence reads to a reference genome and evaluating discordant mapping characteristics using the Sniffles software. This tutorial is based on the Oxford Nanopore [pipeline-structural-variation project](https://github.com/nanoporetech/pipeline-structural-variation) which is available from the linked Github page. The pipeline may be used to identify genomic insertion, duplication and deletion events.

### Features:

The tutorial is packaged with example data from the [Genome in a Bottle](https://www.nist.gov/programs-projects/genome-bottle) project (GM24385). [Truvari](https://github.com/spiralgenetics/truvari/) is used to benchmark SV calling performance against the GIAB truthset of sample-specific known SVs. The analysis considers a collection of sequence reads that map to a 50 Mb region from Chromosome 20. 

* How many sequence reads map to the reference genome?
* What is the depth of sequence coverage across the genome?
* How many SVs can be identified, what is their size distribution, and how frequent are the different types of SV?
* How do the SVs called compare to a reference truthset when accuracy and precision of predictions are considered?
* How can we visualise candidate SVs using the IGV software?


******************

# 2. Getting started:

### What are the bioinformatics objectives?

There are five goals for this tutorial:

* To introduce a literate framework for analysing structural variation from Oxford Nanopore DNA sequence data
* To introduce and utilise best data-management practices to ensure reproducibility
* To map sequence reads to the reference genome and to discover structural varations using **`pipeline-structural-variation`**, an SV discovery workflow optimised for Oxford Nanopore DNA sequence data
* To annotate genomic insertion, deletion and duplications against the human reference genome to identify the variants that overlap genes and annotated genome repeat units
* To benchmark the performance of SV discovery through assessment of precision and recall using a structural variation truthset


### Input and output: 

The tutorial workflow uses a  configuration file, **`config.yaml`**, that specifies the DNA sequences to use, the reference genome and analysis characteristics. This file can be modified to run the workflow using different DNA sequence collections, alternative reference genomes, and with parameters tuned to the needs of SV discovery. 

### What will this tutorial produce?

* A rich HTML format report containing summary statistics and figures highlighting the performance of the long sequence read mapping and SV detection analysis
* **`VCF`** format file of total SVs discovered and filtered files of just the high-confidence structural variants
* **`Excel`** format spreadsheet file identifying the SVs and genes that overlap within the genes' protein coding space
* **`Excel`** format spreadsheet file identifying the SVs and annotated repeat sequences that overlap
* **`IGV Coordinates`** and instructions for reviewing candidate genomic regions with **`IGV`**, the Integrated Genomics Viewer
* **`Instructions for Ribbon`**, an interactive web-application that provides functionality for visualising structural variants

### Dependencies:

* Computer running Linux (Centos7, Ubuntu 18_10, Fedora 29)
* Multiple CPU cores are ideal; a 4 core CPU at minimum would be recommended 
* At least 24 Gb RAM - this is sufficient for mapping against the human genome and can report a full MinION flowcell worth of sequence data. The packaged dataset uses a region of human chromosome 20 and 8Gb RAM is sufficient
* At least 25 Gb spare disk space for analysis and data files; using a complete 50X human dataset would require approximately 500 Gb of available storage space.
* Runtime with provided example data - approximately 1 hour (assuming a recent 4 core processor with sufficient RAM)


Software dependencies include

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


### Installation:

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

### Compilation From Source

This tutorial does not contain software that requires compilation. All packages are managed by the **`Conda`** package management software


### Usage: 

In your Conda environment, and in the tutorial working directory,

1. *optional* edit the provided **`config.yaml`** file to match your own sequence files, reference genome and annotation. 
2. Run the **`snakemake`** command to perform the bioinformatics analysis on the specified sequence files. Several analysis steps can benefit from multiple computer cores; use the `-j` parameter to parallise the analysis (this document is assuming that 8 cores are available).
```
    snakemake --snakefile ./pipeline-structural-variation/Snakefile --configfile ./config.yaml --jobs 4 eval
```
3. Render the tutorial report using the command below
```
   R --slave -e 'rmarkdown::render("ont_tutorial_sv.Rmd", "html_document")'
```



******************

# 3. Results

This tutorial workflow will produce a rich description of your sequence library characteristics and the results from the Structural variation analysis using Nanopore long DNA sequence reads. Please visit the tutorial page at [https://community.nanoporetech.com/knowledge/bioinformatics](https://community.nanoporetech.com/knowledge/bioinformatics) for further information

******************

# 4. Help:

### Licence and Copyright:

© 2019 Oxford Nanopore Technologies Ltd.

Bioinformatics-Tutorials is distributed by Oxford Nanopore Technologies under the terms of the MPL-2.0 license.

### FAQs:



### Abbreviations:


* __BAM__ is a compressed and binary version of the __SAM__ file format; please see __SAM__

* __BED__ is a tab-delimited file format and acronym for Browser Extensible Data. In this tutorial BED files are used to define genomic regions and the columns describe respectively [chromosome, start, end, name]

* __knit__ is the command to render an Rmarkdown file. The knitr package is used to embed code, the results of R analyses and their figures within the typeset text from the document. 

* __L50__  describes the number of sequences (or contigs) that are longer than, or equal to, the N50 length and therefore include half the bases of the assembly

* __N50__  describes the length (read length, contig length etc) where half the bases of the sequence collection are contained within reads/contigs of this length or longer

* __SAM__ is a file format and acronym for Sequence Alignment/Map file format. SAM files are a tab-delimited file and store information on the sequence reads that can be mapped to a genome and the confidence that the mapping is correct.

* __QV__  the quality value, -log10(p) that any given base is incorrect. QV may be either at the individual base level, or may be averaged across whole sequences

* __Rmarkdown__ is an extension to markdown. Functional R code can be embedded in a plain-text document and subsequently rendered to other formats including the PDF format of this report.

* __VCF__ (Variant Call Format) is a standard bioinformatics file format for storing gene sequence variations
