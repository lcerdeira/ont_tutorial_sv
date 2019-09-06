\pagebreak

# Customise the tutorial template

The **`R`** code to prepare your own report is included in the distributed **`Rmarkdown`** file. The **`Rmarkdown`** file can be loaded, viewed, and edited using the **`RStudio`** software. 

**Final thoughts.** Behind this **`Rmarkdown`** file is a modest amount of **`R code`** - please explore the **`Rmarkdown`** template; modify it and run with your own samples.

To extract the whole set of **`R code`** from the **`Rmarkdown`**, use the **`purl`** command - this will extract the R code into its own file.

```
knitr::purl("ont_tutorial_sv.Rmd", quiet=TRUE)
```


# Glossary of terms

* __BAM__ is a compressed and binary version of the __SAM__ file format; please see __SAM__

* __BED__ is a tab-delimited file format and acronym for Browser Extensible Data. In this tutorial BED files are used to define genomic regions and the columns describe respectively [chromosome, start, end, name]

* __knit__ is the command to render an Rmarkdown file. The knitr package is used to embed code, the results of R analyses and their figures within the typeset text from the document. 

* __L50__  describes the number of sequences (or contigs) that are longer than, or equal to, the N50 length and therefore include half the bases of the assembly

* __N50__  describes the length (read length, contig length etc) where half the bases of the sequence collection are contained within reads/contigs of this length or longer

* __SAM__ is a file format and acronym for Sequence Alignment/Map file format. SAM files are a tab-delimited file and store information on the sequence reads that can be mapped to a genome and the confidence that the mapping is correct.

* __QV__  the quality value, -log10(p) that any given base is incorrect. QV may be either at the individual base level, or may be averaged across whole sequences

* __Rmarkdown__ is an extension to markdown. Functional R code can be embedded in a plain-text document and subsequently rendered to other formats including the PDF format of this report.

* __VCF__ (Variant Call Format) is a standard bioinformatics file format for storing gene sequence variations

