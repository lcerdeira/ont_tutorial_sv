parsedBam1aryReads <- parseBamFile(bamFile)
parsedBam2aryReads <- parseBamFile(bamFile, level="SECONDARY")
parsedBamSpryReads <- parseBamFile(bamFile, level="SUPPLEMENTARY")
parsedBamUnmdReads <- harvestUnmappedReads(bamFile)


primaryMappings = parsedBam1aryReads
secondaryMappings=parsedBam2aryReads
supplementaryMappings=parsedBamSpryReads
failedMappings=parsedBamUnmdReads


primaryMappings = parsedBam1aryReads
secondaryMappings=NULL
supplementaryMappings=NULL
failedMappings=NULL








collateMappingCharacteristics <- function(primaryMappings, secondaryMappings=NULL, supplementaryMappings=NULL, failedMappings=NULL) {
  # handle empty secondaryMappings / supplementaryMappings / failedMappings
  # this does assume that parsedBamFile makes sense ...
  if (is.null(secondaryMappings)) {
    secondaryMappings <- primaryMappings[0,]
  }
  if (is.null(supplementaryMappings)) {
    supplementaryMappings <- primaryMappings[0,]
  }
  if (is.null(failedMappings)) {
    failedMappings <- data.frame(qname=character(),
                                 flag=character(), 
                                 qual=character(), 
                                 stringsAsFactors=FALSE)
  }
  
  # basic counts for #s of reads
  mappedSeqs <- sum(primaryMappings$readStarts) 
  unmappedSq <- nrow(failedMappings)
  totalReads <- mappedSeqs + unmappedSq
  
  # basic counts for #s of nucleotides
  mappedNts <- sum(primaryMappings$basesReadsStarted) 
  unmappedNts <- sum(width(failedMappings$qual))
  fastqNts <- mappedNts + unmappedNts
  mappedClippedNts <- sum(primaryMappings$cigarMapped)
  supplNts <- sum(sum(supplementaryMappings$basesReadsStarted, na.rm=TRUE),sum(secondaryMappings$basesReadsStarted, na.rm=TRUE), na.rm=TRUE)
  # basic counts for read lengths
  mappedLength <- mean(rep(primaryMappings$readLen, primaryMappings$readStarts))
  unmappedLength <- mean(width(failedMappings$qual))
  # basic counts for quality scores
  mappedQuality <- phredmean(rep(primaryMappings$readQ, primaryMappings$readStarts))
  unmappedQuality <- 0
  if (nrow(failedMappings) > 0)
    unmappedQuality <- phredmean(unlist(lapply(parsedBamUnmdReads$qual, qualToMeanQ)))
  
  #suppQuality <- mean(c(rep(secondaryMappings$readQ, secondaryMappings$readStarts), rep(supplementaryMappings$readQ, supplementaryMappings$readStarts)), na.rm=TRUE)
  # this information is not stored in this version...
  
  # mapping details ...
  mismatches <- sum(primaryMappings$mismatches, na.rm=TRUE)
  insertions <- sum(primaryMappings$cigarInsertionBases, na.rm=TRUE)
  deletions <- sum(primaryMappings$cigarDeletionBases, na.rm=TRUE)
  meanMapQ <- mean(primaryMappings$mapq, na.rm=TRUE)
  meanSSMapQ <- mean(c(secondaryMappings$mapq, supplementaryMappings$mapq), na.rm=TRUE)
  
  # reference genome characteristics
  refSize <- sum(primaryMappings$width)
  Ncount <- sum(primaryMappings$ncount)
  GCcount <- sum(primaryMappings$gccount)
  chromosomeCount <- length(unique(primaryMappings$chrId))
  meanCov <- sum(primaryMappings$meanCov * primaryMappings$width)/sum(primaryMappings$width)
  suppCov <- sum(secondaryMappings$meanCov * secondaryMappings$width)/sum(secondaryMappings$width) + 
    sum(supplementaryMappings$meanCov * supplementaryMappings$width)/sum(supplementaryMappings$width)
  
  
  mappingSummary <- data.frame(metric=character(), count=character(), percentage=character(), stringsAsFactors = FALSE)
  mappingSummary <- mappingSummary %>% add_row(metric="total sequence reads", 
                                               count=scales::comma_format()(totalReads), 
                                               percentage="")
  mappingSummary <- mappingSummary %>% add_row(metric="mapped reads (primary)", 
                                               count=scales::comma_format()(mappedSeqs), 
                                               percentage=paste0(round(mappedSeqs / totalReads * 100, digits=2),"%"))
  mappingSummary <- mappingSummary %>% add_row(metric=".....multi-mapping reads (secondary)", 
                                               count=scales::comma_format()(sum(secondaryMappings$readStarts)), 
                                               percentage="")
  mappingSummary <- mappingSummary %>% add_row(metric=".....multi-mapping reads (supplementary)", 
                                               count=scales::comma_format()(sum(supplementaryMappings$readStarts)), 
                                               percentage="")
  mappingSummary <- mappingSummary %>% add_row(metric="unmapped reads", 
                                               count=scales::comma_format()(unmappedSq), 
                                               percentage=paste0(round(unmappedSq / totalReads * 100, digits=2),"%"))
  mappingSummary <- mappingSummary %>% add_row(metric="total fastq nucleotides", 
                                               count=scales::comma_format()(fastqNts), 
                                               percentage="")
  mappingSummary <- mappingSummary %>% add_row(metric="sum of nucleotides from mapped reads", 
                                               count=scales::comma_format()(mappedNts), 
                                               percentage=paste0(round(mappedNts / fastqNts * 100, digits=2),"%"))
  mappingSummary <- mappingSummary %>% add_row(metric=".....nucleotides from clipped mapping", 
                                               count=scales::comma_format()(mappedClippedNts), 
                                               percentage=paste0(round(mappedClippedNts/fastqNts * 100, digits=2),"%"))  
  mappingSummary <- mappingSummary %>% add_row(metric="unmapped nucleotides", 
                                               count=scales::comma_format()(unmappedNts), 
                                               percentage=paste0(round(unmappedNts / fastqNts * 100, digits=2),"%"))    
  mappingSummary <- mappingSummary %>% add_row(metric="supp/secondary nucleotides", 
                                               count=scales::comma_format()(supplNts), 
                                               percentage="")  
  mappingSummary <- mappingSummary %>% add_row(metric="mean read length (mapped)", 
                                               count=scales::comma_format()(mappedLength), 
                                               percentage="")  
  mappingSummary <- mappingSummary %>% add_row(metric="mean read length (unmapped)", 
                                               count=scales::comma_format()(unmappedLength), 
                                               percentage="")  
  mappingSummary <- mappingSummary %>% add_row(metric="mean read length (supp/second)", 
                                               count=scales::comma_format()(mean(c(secondaryMappings$readLen, supplementaryMappings$readLen), na.rm=TRUE)), 
                                               percentage="")  
  mappingSummary <- mappingSummary %>% add_row(metric="mean read quality (mapped)", 
                                               count=round(mappedQuality, digits=2), 
                                               percentage="")  
  mappingSummary <- mappingSummary %>% add_row(metric="mean read quality (unmapped)", 
                                               count=round(unmappedQuality, digits=2), 
                                               percentage="")  
  mappingSummary <- mappingSummary %>% add_row(metric="mapping mismatches", 
                                               count=scales::comma_format()(mismatches), 
                                               percentage="")  
  mappingSummary <- mappingSummary %>% add_row(metric="base insertions (INS)", 
                                               count=scales::comma_format()(insertions), 
                                               percentage="")  
  mappingSummary <- mappingSummary %>% add_row(metric="base deletions (DEL)", 
                                               count=scales::comma_format()(deletions), 
                                               percentage="")  
  mappingSummary <- mappingSummary %>% add_row(metric="mean mapping quality (primary)", 
                                               count=round(meanMapQ, digits=2), 
                                               percentage="")  
  mappingSummary <- mappingSummary %>% add_row(metric="mean mapping quality (supp/second)", 
                                               count=round(meanSSMapQ, digits=2), 
                                               percentage="")  
  mappingSummary <- mappingSummary %>% add_row(metric="Chromosome Count", 
                                               count=chromosomeCount, 
                                               percentage="")  
  mappingSummary <- mappingSummary %>% add_row(metric="Reference genome size (nt)", 
                                               count=scales::comma_format()(refSize), 
                                               percentage="")  
  mappingSummary <- mappingSummary %>% add_row(metric="N count across reference", 
                                               count=scales::comma_format()(Ncount), 
                                               percentage=paste0(round(Ncount / refSize * 100, digits=2),"%"))  
  mappingSummary <- mappingSummary %>% add_row(metric="GC count", 
                                               count=scales::comma_format()(GCcount), 
                                               percentage=paste0(round(GCcount / (refSize - Ncount) * 100, digits=2),"%"))  
  mappingSummary <- mappingSummary %>% add_row(metric="Mean coverage (primary)", 
                                               count=round(meanCov, digits=2), 
                                               percentage="")  
  mappingSummary <- mappingSummary %>% add_row(metric="Mean coverage (supp/second)", 
                                               count=round(suppCov, digits=2), 
                                               percentage="")  
  
  rownames(mappingSummary) <- mappingSummary[,1]
  mappingSummary <- mappingSummary[,-1]       
  
  return(mappingSummary)
}


mappingCharacteristics <- collateMappingCharacteristics(parsedBam1aryReads, parsedBam2aryReads, parsedBamSpryReads, parsedBamUnmdReads) 


row.names(mappingCharacteristics)[6] <- paste0(row.names(mappingCharacteristics)[6], footnote_marker_symbol(1, "html"))
row.names(mappingCharacteristics)[8] <- paste0(row.names(mappingCharacteristics)[8], footnote_marker_symbol(2, "html"))
row.names(mappingCharacteristics)[16]<- paste0(row.names(mappingCharacteristics)[16], footnote_marker_symbol(3, "html"))
row.names(mappingCharacteristics)[21]<- paste0(row.names(mappingCharacteristics)[21], footnote_marker_symbol(4, "html"))
row.names(mappingCharacteristics)[24]<- paste0(row.names(mappingCharacteristics)[24], footnote_marker_symbol(5, "html"))
row.names(mappingCharacteristics)[25]<- paste0(row.names(mappingCharacteristics)[25], footnote_marker_symbol(6, "html"))

kable(mappingCharacteristics, format="html", col.names=rep(" ", ncol(mappingCharacteristics)), caption="Table summarising mapping characteristics", booktabs=TRUE, table.envir='table*', linesep="", escape = FALSE)  %>%
  #add_header_above(c(" ", "NGMLR"=2, "Minimap2"=2)) %>%
  kable_styling(c("striped", "condensed")) %>%
  pack_rows("Nucleotides mapped", 6, 10) %>%
  pack_rows("Sequence length characteristics (nt)", 11, 13) %>%
  pack_rows("Sequence quality characteristics (phred Qval)", 14, 15) %>%
  pack_rows("Mapping quality characteristics", 16, 20) %>%
  pack_rows("Reference genome characteristics", 21, 26) %>%
  footnote(symbol=c("fastq bases are calculated from the qwidth field of the mapped sequences and from the sequence length of unmapped sequences","clipped mapping calculated from CIGARquery coordinates", "mismatching bases reported by SAM NM edit distance tag", "standard workflow attempts to remove the mitochondrial genome", "The GC content is calculated as % of GC at positions where nucleotide is A/C/G/T", "depth of coverage based only on primary mapping reads"), symbol_title="please note: ", footnote_as_chunk = TRUE)


```