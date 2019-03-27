#' @export
#' @import GenomicRanges
#' @importFrom methods is
#' @importFrom utils installed.packages txtProgressBar setTxtProgressBar
#'
#' @title Assess Genes
#' @description Assess and categorize a set of genes for a genome using proteomics hits, evolutionarily conserved starts,
#' and evolutionarily conserved stops as evidence
#'
#' @usage
#' AssessGenes(geneLeftPos,
#'             geneRightPos = NA_integer_,
#'             geneStrand = NA_character_,
#'             inputMapObj,
#'             geneSource = "",
#'             minCovNum = 10,
#'             minCovPct = 5,
#'             minConCovRatio_Strong = 0.99,
#'             limConCovRatio_NotCon = 0.8,
#'             minConCovRatio_Stop = 0.5,
#'             noConStopsGeneFrac = 0.5,
#'             minMissORFLen = 0,
#'             allowNestedORFs = FALSE,
#'             useNTermProt = FALSE,
#'             verbose = TRUE)
#'
#' @param geneLeftPos An integer vector with the left positions of each gene, in terms of the forward strand.
#' Can also be a \code{GRanges} object from the \code{GenomicRanges} package that holds all of the positional
#' information (including strand) for the genes. In that case, the next two parameters should be left as NA.
#'
#' @param geneRightPos An integer vector with the right positions of each gene, in terms of the forward strand.
#' Should be left at the default value of \code{NA_integer_} if \code{geneLeftPos} is a GRanges object.
#'
#' @param geneStrand A character vector consisting of "+" and "-", specifying which strand each gene is on.
#' Should be left at the default value of \code{NA_character_} if \code{geneLeftPos} is a GRanges object.
#'
#' @param inputMapObj EITHER an object of class \code{Assessment} and subclass \code{DataMap} OR a character string
#' corresponding to the strain identifier for one of such objects from \code{AssessORFData}.
#' 
#' @param geneSource Optional character string that describes the source of the gene set, i.e. a database or gene prediction
#' program. Used when viewing and identifying the object returned by the function.
#' 
#' @param minCovNum Minimum number of related genomes required to have synteny to a position in the central genome.
#' Recommended to use the default value.
#' 
#' @param minCovPct Minimum percentage of related genomes required to have synteny to a position in the central genome.
#' Must be an integer ranging from 0 to 100. Recommended to use the default value.
#'
#' @param minConCovRatio_Strong Minimum value of the start codon conservation to coverage ratio needed to call a start strongly conserved.
#' Must range from 0 to 1. Lower values allow more conserved starts through. Recommended to use the default value.
#' 
#' @param limConCovRatio_NotCon Maximum, non-inclusive value of the conservation to coverage ratio needed to call a possible conserved
#' start not conserved. Used when making a decision on how to categorize the conserved start evidence. Must range from 0 to 1
#' Recommended to use the default value.
#' 
#' @param minConCovRatio_Stop Minimum value of the stop codon conservation to coverage ratio needed to say a position in the central
#' genome corresponds to a conserved stop across the related genomes. Must range from 0 to 1. Lower values allow more conserved stops
#' through. Recommended to use the default value.
#' 
#' @param noConStopsGeneFrac Value from 0 to 1 describing the fractional range of positions in a gene, starting from the start of the
#' gene and moving towards the stop of the gene, to use in searching for conserved stops. For example, a value of 0.25 means that the
#' first quarter of the gene is checked for conserved stops, a value of 0.5 correspond to the first half of the gene, etc. Recommended
#' to use the default value.
#' 
#' @param minMissORFLen Minimum ORF length required to include an ORF with protein hits
#' but no given/predicted gene start in the final output.
#' 
#' @param allowNestedORFs Logical indicating whether or not to include ORFs with protein hits but no given/predicted
#' gene starts that are completely nested within an ORF in another frame in the final output.
#'
#' @param verbose Logical indicating whether or not to display progress and status messages.
#' 
#' @param useNTermProt Logical indicating whether or not to treat proteomics evidence in the given mapping object as originating
#' from N-terminal proteomics experiments. The mapping object must be built with N-terminal proteomics data.
#' Default value is FALSE.
#'
#' @details
#' For each of the given genes, \code{AssessGene} assigns a category based on where conserved starts, conserved stops, and/or
#' proteomics hits are located in relation to the start of the gene. The category assignments for the genes are stored in the
#' \code{CategoryAssignments} vector in the \code{Results} object returned by the function. Please see
#' \code{\link{Assessment-class}} for a list of all possible categories and their descriptions.
#' 
#' If \code{geneLeftPos} is a \code{GRanges} object, then the left and right positions of each gene along with the strand of each
#' gene are extracted from the object. Any sequence names given for the genes within the \code{GRanges} object are ignored, and
#' the \code{CategoryAssignments} in the returned \code{Results} object follows the same order as to how the genes are listed
#' within the \code{GRanges} object.
#' 
#' If gene positional information is instead given as three vectors, then the three vectors, \code{geneLeftPos}, \code{geneRightPos},
#' and \code{geneStrand}, must all be of the same length. The same index within each vector must provide information on the same gene
#' (think of the vectors as columns of the same table). \code{geneLeftPos} and \code{geneRightPos} describe the upstream and downstream
#' positions (respectively) for each gene in terms of the forward strand. For genes on the forward strand, \code{geneLeftPos} corresponds
#' to the start positions and \code{geneRightPos} corresponds to stop positions. For genes on the reverse strand, \code{geneLeftPos}
#' corresponds to the stop positions and \code{geneRightPos} corresponds to the start positions. Gene positions on the reverse strand
#' must be relative to the 5' to 3' direction of the forward strand (as opposed to being relative to the 5' to 3' direction of the reverse
#' strand). This means that none of the elements of \code{geneLeftPos} can be greater than (or equal to) the corresponding element in
#' \code{geneRightPos}. The \code{CategoryAssignments} in the returned \code{Results} object has the same length as and aligns with the
#' indexing of the three given gene positional information vectors. 
#'
#' Please ensure that the same genome used in the mapping function is also used to derive the set of genes for this assessment function.
#' The function will only error if any gene positions are outside the bounds of the genome and does not make any other checks to make sure
#' the genes are valid for the genome.
#'
#' The maximum of either \code{minCovNum} (option 1) or \code{minCovPct} divided 100 then multiplied by the number of related genomes
#' (option 2) is used as the minimum coverage required in determining conserved starts and stops.
#' 
#' Additionally, open reading frames with proteomics evidence but no gene start are categorized based on whether or not there
#' is a conserved start upstream of the proteomic evidence. The positions and lengths of these open reading frames are included
#' in the \code{N_CS-_PE+_ORFs} and \code{N_CS+_PE+_ORFs} matrices within the final object that is returned.
#' 
#' If the proteomics evidence provided in the given mapping object comes from N-terminal proteomics experiments (i.e., if the value of
#' the \code{NTermProteomics} item within the mapping object is TRUE), the \code{useNTermProt} can be set to TRUE to impose stricter
#' requirements on the use of proteomics evidence in determining the correctness of the given genes. When \code{useNTermProt} is set to
#' TRUE, the start of first peptide mapping to an ORF where there is a given gene must directly align with the start of that gene or be
#' one codon off from the start (in cases where the protein product of the gene has undergone N-terminal methionine excision) in order for
#' the gene to be considered as having supporting protein evidence. If the first peptide hit does not align like that, the gene is
#' considered as having disproving protein evidence. Currently, N-terminal proteomics does not produce enough N-terminal peptides so
#' setting this flag as TRUE does not provide meaningful results. It is recommended to leave this flag as FALSE in all situations.
#'
#' @return An object of class \code{Assessment} and subclass \code{Results}
#' 
#' @seealso \code{\link{Assessment-class}}
#' 
#' @examples
#'
#' ## Example showing the minimum number of arguments that need to be specified:
#'
#' \dontrun{
#' myResObj <- AssessGenes(geneLeftPos = myGenesLeft,
#'                         geneRightPos = myGenesRight,
#'                         geneStrand = myGenesStrand,
#'                         inputMapObj = myMapObj)
#' }
#'
#'
#'
#' ## Example from vignette is shown below
#'
#' currMapObj <- readRDS(system.file("extdata",
#'                                   "MGAS5005_PreSaved_DataMapObj.rds",
#'                                   package = "AssessORF"))
#'
#' currProdigal <- readLines(system.file("extdata",
#'                                       "MGAS5005_Prodigal.sco",
#'                                       package = "AssessORF"))[-1:-2]
#'
#' prodigalLeft <- as.numeric(sapply(strsplit(currProdigal, "_", fixed=TRUE), `[`, 2L))
#' prodigalRight <- as.numeric(sapply(strsplit(currProdigal, "_", fixed=TRUE), `[`, 3L))
#' prodigalStrand <- sapply(strsplit(currProdigal, "_", fixed=TRUE), `[`, 4L)
#' 
#' currResObj <- AssessGenes(geneLeftPos = prodigalLeft,
#'                           geneRightPos = prodigalRight,
#'                           geneStrand = prodigalStrand,
#'                           inputMapObj = currMapObj,
#'                           geneSource = "Prodigal")
#'
#' print(currResObj)
#'
AssessGenes <- function(geneLeftPos,
                        geneRightPos = NA_integer_,
                        geneStrand = NA_character_,
                        inputMapObj,
                        geneSource = "",
                        minCovNum = 10L,
                        minCovPct = 5L,
                        minConCovRatio_Strong = 0.99,
                        limConCovRatio_NotCon = 0.8,
                        minConCovRatio_Stop = 0.5,
                        noConStopsGeneFrac = 0.5,
                        minMissORFLen = 0L,
                        allowNestedORFs = FALSE,
                        useNTermProt = FALSE,
                        verbose = TRUE) {
  
  ## Check inputs for error.
  
  if ((!is.logical(verbose)) ||(anyNA(verbose)) || (length(verbose) != 1)) {
    stop("'verbose' must be of type logical, be either TRUE or FALSE, and consist of only 1 element.")
  }
  
  ## --------------------------------------------------------------------------------------------------------------- ##
  
  if (anyNA(inputMapObj)) {
    stop("'inputMapObj' cannot have any NAs or have more than one element.")
  }
  
  if ((is(inputMapObj, "Assessment")) && (is(inputMapObj, "DataMap"))) {
    mapObj <- inputMapObj
  } else if (is.character(inputMapObj) && (length(inputMapObj) == 1)){
    if (!(any(installed.packages()[, 1] == "AssessORFData"))) {
      stop("The package 'AssessORFData' has not been installed.")
    }
    
    mapObj <- AssessORFData::GetDataMapObj(inputMapObj)
    
  } else {
    stop("'inputMapObj' must either be an object of class 'Assessment', subclass 'DataMap' ",
         "or be a single character string.")
  }
  
  ## Get the genome length to use for checking gene boundaries.
  genomeLength <- mapObj$GenomeLength
  
  ## --------------------------------------------------------------------------------------------------------------- ##
  
  ## First check to see if gene positional information is given as a GRanges object.
  if (is(geneLeftPos, "GRanges")) {
    ## If so, extract the gene positional information from the object.
    geneSet <- geneLeftPos
    geneLeftPos <-  start(geneSet)
    geneRightPos <- end(geneSet)
    geneStrand <- as.character(strand(geneSet))
  }
  
  ## The gene positional information is (now) given as three vectors.
  ## Check to see if all the information on the genes is accurate.
  ## The left and right position vectors must be of type integer.
  ## Positions cannot extend beyond the bounds of the genome.
  ## The strand vector must be of type character and contain only "+" or "-".
  
  ## The left position vector must be of type integer. Stop otherwise.
  if ((!is.numeric(geneLeftPos))  || (anyNA(geneLeftPos)) || (any(geneLeftPos %% 1 != 0))) {
    stop( "Left positions must be given as integer numbers")
  }
  
  ## Check left positions relative to the genome length. Stop if any positions are out of bounds.
  if ((any(geneLeftPos <= 0L)) || (any(geneLeftPos > genomeLength))) {
    stop("Left positions must be within the bounds of the genome.")
  }
  
  ## The right position vectors must be of type integer. Stop otherwise.
  if ((!is.numeric(geneRightPos))  || (anyNA(geneRightPos)) || (any(geneRightPos %% 1 != 0))) {
    stop("Right positions must be given as integer numbers.")
  }
  
  ## Check right positions relative to the genome length. Stop if any positions are out of bounds.
  if ((any(geneRightPos <= 0L)) || (any(geneRightPos > genomeLength))) {
    stop("Right positions must be within the bounds of the genome.")
  }
  
  ## The strand vector must be of type character. Stop otherwise.
  if ((!is.character(geneStrand)) || (anyNA(geneStrand))) {
    stop("Strand information must be given as a character vector.")
  }
  
  ## The strand vector must contain only "+" or "-". Stop otherwise.
  if (any((geneStrand != "+") & (geneStrand != "-"))) {
    stop("Strand information must consist only of + and -.")
  }
  
  ## All three vectors (left, right, and strand) must be of the same length. Stop otherwise.
  if ((length(geneLeftPos) != length(geneRightPos)) || (length(geneLeftPos) != length(geneStrand))) {
    stop("All three gene position vectors must be of the same length.")
  }
  
  ## Left positions must be less than the corresponding right positions. Stop otherwise.
  if (any(geneLeftPos >= geneRightPos)) {
    stop("Left positions fpr all genes must be strictly less than the corresponding right positions.")
  }
  
  
  ## --------------------------------------------------------------------------------------------------------------- ##
  
  if ((!is.character(geneSource)) || (anyNA(geneSource))) {
    stop("'geneSource' must be a valid character string.")
  }
  
  if (length(geneSource) != 1) {
    stop("'geneSource' must consist of exactly one character string.")
  }
  
  if ((!is.numeric(minCovNum))  || (anyNA(minCovNum))) {
    stop("The minimum number of related genomes required must be a valid real number.")
  }
  
  if (length(minCovNum) != 1) {
    stop("Exactly one number must be inputted as the minimum number of related genomes required.")
  }
  
  if ((minCovNum %% 1 != 0) || (minCovNum <= 0)) {
    stop("The minimum number of related genomes required must be a ",
         "non-negative integer that is greater than 0.")
  }
  
  if ((!is.numeric(minCovPct))  || (anyNA(minCovPct))) {
    stop("The minimum percentage of related genomes required must be a valid real number.")
  }
  
  if (length(minCovPct) != 1) {
    stop("Exactly one number must be inputted as the minimum percentage of related genomes required.")
  }
  
  if ((minCovPct %% 1 != 0) || (minCovPct < 0) || (minCovPct > 100)) {
    stop("The minimum percentage of related genomes required must be a ",
         "non-negative integer that is greater than 0 and less than 100 (both inclusive).")
  }
  
  if ((!is.numeric(minConCovRatio_Strong)) || (anyNA(minConCovRatio_Strong))) {
    stop("'minConCovRatio_Strong' must be a valid real number.")
  }
  
  if (length(minConCovRatio_Strong) != 1) {
    stop("'minConCovRatio_Strong' must consist of exactly one number.")
  }
  
  if ((minConCovRatio_Strong <= 0) || (minConCovRatio_Strong >= 1)) {
    stop("'minConCovRatio_Strong' must be greater than 0 and less than 1.")
  }
  
  if ((!is.numeric(limConCovRatio_NotCon)) || (anyNA(limConCovRatio_NotCon))) {
    stop("'limConCovRatio_NotCon' must be a valid real number.")
  }
  
  if (length(limConCovRatio_NotCon) != 1) {
    stop("'limConCovRatio_NotCon' must consist of exactly one number.")
  }
  
  if ((limConCovRatio_NotCon <= 0) || (limConCovRatio_NotCon >= 1)) {
    stop("'limConCovRatio_NotCon' must be greater than 0 and less than 1.")
  }
  
  if ((!is.numeric(minConCovRatio_Stop)) || (anyNA(minConCovRatio_Stop))) {
    stop("'minConCovRatio_Stop' must be a valid real number.")
  }
  
  if (length(minConCovRatio_Stop) != 1) {
    stop("'minConCovRatio_Stop' must consist of exactly one number.")
  }
  
  if ((minConCovRatio_Stop <= 0) || (minConCovRatio_Stop >= 1)) {
    stop("'minConCovRatio_Stop' must be greater than 0 and less than 1.")
  }
  
  if ((!is.numeric(noConStopsGeneFrac)) || (anyNA(noConStopsGeneFrac))) {
    stop("'noConStopsGeneFrac' must be a valid real number.")
  }
  
  if (length(noConStopsGeneFrac) != 1) {
    stop("'noConStopsGeneFrac' must consist of exactly one number.")
  }
  
  if ((noConStopsGeneFrac <= 0) || (noConStopsGeneFrac >= 1)) {
    stop("'noConStopsGeneFrac' must be greater than 0 and less than 1.")
  }
  
  if ((!is.numeric(minMissORFLen)) || (anyNA(minMissORFLen))) {
    stop("'minMissORFLen' must be a valid real number.")
  }
  
  if (length(minMissORFLen) != 1) {
    stop("'minMissORFLen' must consist of exactly one number.")
  }
  
  if ((minMissORFLen %% 1 != 0) || (minMissORFLen < 0)) {
    stop("'minMissORFLen' must be a non-negative integer greater than or equal to 0.")
  }
  
  if ((!is.logical(allowNestedORFs)) ||(anyNA(allowNestedORFs)) || (length(allowNestedORFs) != 1)) {
    stop("'allowNestedORFs' must be of type logical, be either TRUE or FALSE, and consist of only 1 element.")
  }
  
  if ((!is.logical(useNTermProt)) ||(anyNA(useNTermProt)) || (length(useNTermProt) != 1)) {
    stop("'useNTermProt' must be of type logical, be either TRUE or FALSE, and consist of only 1 element.")
  }
  
  ## --------------------------------------------------------------------------------------------------------------- ##
  
  ## Pull the necessary vectors and lists from the mapping object.
  
  stops <- mapObj$StopsByFrame
  
  fwdProt <- mapObj$FwdProtHits
  revProt <- mapObj$RevProtHits
  
  isNTerm <- mapObj$NTermProteomics && useNTermProt
  useProt <- mapObj$HasProteomics
  
  fwdCov <- mapObj$FwdCoverage
  fwdConStart <- mapObj$FwdConStarts
  fwdConStop <- mapObj$FwdConStops
  
  revCov <- mapObj$RevCoverage
  revConStart <- mapObj$RevConStarts
  revConStop <- mapObj$RevConStops
  
  useCons <- mapObj$HasConservation
  
  ## --------------------------------------------------------------------------------------------------------------- ##
  
  ## Check if there are enough related genomes to use conservation evidence when categorizing genes.
  if ((useCons) && (minCovNum > mapObj$NumRelatedGenomes)) {
    useCons <- FALSE
    
    warning("The value for the minimum number of related genomes required is greater than the number of related ",
            "genomes used to generate the mapping object. Conserved starts & stops will not be used in assessment.")
  }
  
  ## If conservation evidence can be used, determine the minimum coverage, using either the minimum
  ## number of related genomes required or the minimum percentage of related genomes required.
  if (useCons) {
    minCov <- max(minCovNum, minCovPct * mapObj$NumRelatedGenomes / 100)
  }
  
  ## --------------------------------------------------------------------------------------------------------------- ##
  
  ## If the predicted stop is downstream of the current in-frame region-ending
  ## stop, this function will check which of the next stop-to-stop regions need
  ## to be skipped, making sure no predicted genes will be skipped in the
  ## process. Returns the indices of the region-ending stops to skip.
  checkNextStops <- function(currPredStop) {
    skipIdxs <- which((frameStops > regStart) & (frameStops < currPredStop))
    skipIdxs <- c(skipIdxs, skipIdxs[length(skipIdxs)] + 1L)
    
    lastStop <- frameStops[skipIdxs[length(skipIdxs)]]
    
    skippedStarts <- strandFrameStarts[((strandFrameStarts > regEnd) &
                                          (strandFrameStarts < lastStop))]
    
    if (length(skippedStarts) > 0) {
      firstNoSkip <- (which(frameStops[skipIdxs] > skippedStarts[1]))[1]
      
      skipIdxs <- skipIdxs[seq(from = 1, to = firstNoSkip - 1, by = 1)]
    }
    
    return(skipIdxs)
  }
  
  ## This function checks if there is overlapping protein evidence in other
  ## frames. It is used when recording information for categories with protein
  ## evidence and no gene start.
  checkOtherFrames <- function() {
    currProtHitPos <- regRange[regionalProteinHits]
    
    sameRange <- (min(currProtHitPos)):(max(currProtHitPos))
    
    oppRange <- genomeLength - sameRange + 1L
    
    for (oFrame in (seq(1, 6))[-frameID]) {
      if (frameID <= 3) {
        if (oFrame <= 3) {
          oppStrandProt <- fwdProt[[oFrame]][[1]][sameRange]
        } else {
          oppStrandProt <- revProt[[oFrame - 3]][[1]][oppRange]
        }
        
      } else {
        if (oFrame <= 3) {
          oppStrandProt <- fwdProt[[oFrame]][[1]][sameRange]
        } else {
          oppStrandProt <- revProt[[oFrame - 3]][[1]][regRange]
        }
      }
      
      if (any(oppStrandProt > 0)) {
        return(oFrame)
      }
    }
    
    return(0L)
  }
  
  ## This function checks if the current ORF is completely within another ORF.
  ## It is used to filter out some of the ORFs with protein hits and no given start.
  isNestedORF <- function() {
    currProtHitPos <- regRange[regionalProteinHits]
    
    sameLeftPos <- min(currProtHitPos)
    sameRightPos <- max(currProtHitPos)
    
    oppLeftPos <- genomeLength - sameRightPos + 1L
    oppRightPos <- genomeLength - sameLeftPos + 1L
    
    for (oFrame in (seq(1, 6))[-frameID]) {
      
      if (((frameID <= 3) && (oFrame <= 3)) || ((frameID > 3) && (oFrame > 3))) {
        ## Use same
        nested <- any((stops[[oFrame]][-length(stops[[oFrame]])] <= sameLeftPos) &
                        (stops[[oFrame]][-1] >= sameRightPos))
      } else {
        ## Use opposite
        nested <- any((stops[[oFrame]][-length(stops[[oFrame]])] <= oppLeftPos) &
                        (stops[[oFrame]][-1] >= oppRightPos))
      }
      
      if (nested) {
        return(TRUE)
      }
    }
    
    return(FALSE)
  }
  
  ## This matrix stores information on ORFs with protein hits but no strong conserved start
  ## upstream of the protein hits and no predicted gene
  cat1A_ORFs <- matrix(0, nrow = sum(vapply(stops, length, integer(1))), ncol = 5,
                       dimnames = list(NULL,
                                       c("Start", "End", "Length", "Frame", "OtherProtFrame")))
  cat1A_Counter <- 0L
  
  ## This matrix stores information on ORFs with protein hits and at least one strong conserved
  ## start upstream of the protein hits but no predicted gene.
  cat1B_ORFs <- matrix(0, nrow = sum(vapply(stops, length, integer(1))), ncol = 5,
                       dimnames = list(NULL,
                                       c("Start", "End", "Length", "Frame", "OtherProtFrame")))
  cat1B_Counter <- 0L
  
  ## This vector stores the categorization for each predicted gene.
  catGeneAssignment <- character(length(geneLeftPos))
  
  if (verbose) {
    pBar <- txtProgressBar(style=ifelse(interactive(), 3, 1))
  }
  
  ## Iterate through each of the 6 frames of the genome, using
  ## the same set of vectors to work through each frame.
  for (frameID in seq(1, 6)) {
    if (frameID <= 3) {
      ## Get and set up the vectors to work through the current forward frame.
      strandID <- "+"
      cFrame <- frameID
      shift <- cFrame %% 3
      
      frameProteinHits <- fwdProt[[cFrame]][[1]]
      
      strandCoverage <- fwdCov
      strandConservedStarts <- fwdConStart
      strandConservedStops <- fwdConStop
      
      strandOnly <- which(geneStrand == strandID)
      
      strandStarts <- geneLeftPos[strandOnly]
      strandStops <- geneRightPos[strandOnly]
    } else {
      ## Get and set up the vectors to work through the current reverse frame.
      strandID <- "-"
      cFrame <- frameID - 3
      shift <- cFrame %% 3
      
      frameProteinHits <- revProt[[cFrame]][[1]]
      
      strandCoverage <- revCov
      strandConservedStarts <- revConStart
      strandConservedStops <- revConStop
      
      strandOnly <- which(geneStrand == strandID)
      
      strandStarts <- genomeLength - geneRightPos[strandOnly] + 1L
      strandStops <- genomeLength - geneLeftPos[strandOnly] + 1L
    }
    
    ## Add in names to the strand-specific gene position vectors so that the
    ## categorization of genes can be mapped back to the order the genes were
    ## inputted into the function.
    names(strandStarts) <- names(strandStops) <- as.character(strandOnly)
   
    ## Get and set up the predicted gene position vectors for the current frame. 
    strandFrameStarts <- strandStarts[(strandStarts %% 3 == shift)]
    strandFrameStops <- strandStops[(strandStarts %% 3 == shift)]
    
    ## In order to handle genes at the boundary of the genome, additional
    ## in-frame positions need to be added to the stops vector for the frame:
    ## the start of the first possible in-frame codon and the start(s) of the
    ## last one or two possible in-frame codon. The start of the last codon is
    ## always included. The start of the second to last possible codon is only
    ## included if the last codon is not completely within the genome.a
    finalCodonStart <- ifelse(genomeLength %% 3 == shift, genomeLength,
                              ifelse((genomeLength - 1L) %% 3 == shift,
                                     genomeLength - 1L, genomeLength - 2L))
    
    if ((finalCodonStart + 2L) == genomeLength) {
      frameStops <- sort(unique(c(cFrame, stops[[frameID]],
                                  finalCodonStart)))
    } else {
      frameStops <- sort(unique(c(cFrame, stops[[frameID]],
                                  finalCodonStart - 3L, finalCodonStart)))
    }
    
    ## If a predicted gene has a stop that is beyond the end of region, it will
    ## be necessary to skip subsequent genome stops until the predicted stop.
    skipStop <- logical(length(frameStops))
    
    ## For each frame, loop through each stop-to-stop region and categorize the
    ## gene(s) or proteomics hits within that region if there are any.
    for (idx in seq(1, length(frameStops) - 1)) {
      ## Skip the current region based on the region-ending stop. See below and
      ## the checkNextStops() function for how region skipping is determined.
      if (skipStop[idx + 1]) {
        next()
      }
      
      ## Set boundaries of the region to examine.
      ## The first variable corresponds to the region-starting stop
      ## The second variable corresponds to the region-ending stop.
      ## The third variable is the range between those two positions.
      regStart <- frameStops[idx]
      regEnd <- frameStops[idx + 1]
      regRange <- regStart:regEnd
      
      ## For each region, start by assuming there are no valid conserved
      ## starts, and no proteomics hits mapping to that region.
      noConStarts <- TRUE
      noProt <- TRUE
      
      ## Get the gene start(s) and gene stop(s) in that region.
      regionalPredictedStart <- strandFrameStarts[((strandFrameStarts >= regStart) &
                                                     (strandFrameStarts <= regEnd))]
      regionalPredictedStop <- strandFrameStops[((strandFrameStarts >= regStart) &
                                                   (strandFrameStarts <= regEnd))]
      
      ## Get the index for the gene(s).
      ## This is based on the order the genes were inputted into the function.
      geneIdx <- as.integer(names(regionalPredictedStart))
      
      ## If there is more than one gene in the region, warn the user and
      ## put the genes into the no evidence category.
      if (length(regionalPredictedStart) > 1) {
        ## Get the forward strand positions for both the region
        ## and the genes to use in the warning message.
        if (frameID <= 3) {
          fwdRegStart <- regStart
          fwdRegEnd <- regEnd
          fwdPredStart <- regionalPredictedStart
          fwdPredStop <- regionalPredictedStop
        } else {
          fwdRegStart <- genomeLength - regStart + 1L
          fwdRegEnd <- genomeLength - regEnd + 1L
          fwdPredStart <- genomeLength - regionalPredictedStart + 1L
          fwdPredStop <- genomeLength - regionalPredictedStop + 1L
        }
        
        ## Warn the user about multiple genes within a particular region.
        warning("More than one in-frame gene start between stops ",
                fwdRegStart, " and ", fwdRegEnd, " in frame ", frameID, ".",
                paste0(" A gene start is at ", fwdPredStart, " with stop at ",
                       fwdPredStop, " (gene ", geneIdx, ")."),
                " These genes will be placed into the no evidence category.")
        
        lastPredStop <- regionalPredictedStop[length(regionalPredictedStart)]
        
        ## If the most downstream predicted stop in the set of genes is
        ## downstream of the region-ending stop, skip subsequent genome
        ## stops as long as no predicted genes will be skipped.
        if ((lastPredStop - 2L) > regEnd) {
          skipStop[checkNextStops(lastPredStop)] <- TRUE
        }
        
        ## Assign those genes to the no evidence category.
        catGeneAssignment[geneIdx] <- "Y CS- PE-"
        next()
      }
      
      ## If there is one gene in the region, make sure the stop makes sense.
      ## Note that the region-ending stop position corresponds to the first
      ## position of the stop codon while the predicted stop position should
      ## correspond to the last position of the stop codon.
      if (length(regionalPredictedStart) == 1) {
        ## Check if the predicted stop is beyond the end of the region.
        beyondEnd <- ((regionalPredictedStop - 2L) > regEnd)
        
        ## Check if the predicted stop is out of frame.
        outOfFrame <- ((regionalPredictedStop - 2L) %% 3 != shift)
        
        ## Check if the predicted stop is the before the end of the region.
        beforeEnd <- ((regionalPredictedStop - 2L) < regEnd)
        
        ## If the predicted stop is not correct for the gene, warn the user
        ## appropriately and place the gene in the no evidence category.
        if (beyondEnd || outOfFrame || beforeEnd) {
          ## Get the forward strand positions for both the correct stop
          ## and the predicted stop to use in the warning message(s).
          if (frameID <= 3) {
            fwdRightStop <- regEnd + 2L
            fwdPredStop <- regionalPredictedStop
          } else {
            fwdRightStop <- genomeLength - regEnd + 1L
            fwdPredStop <- genomeLength - regionalPredictedStop + 1L
          }
          
          ## Print the warning message.
          ## Skip subsequent genome stops if necessary.
          if (beyondEnd) {
            warning("The predicted stop for gene ", geneIdx, ", position ", fwdPredStop,
                    ", is downstream of the correct stop. The stop should be at ",
                    fwdRightStop, ". The gene is in frame ", frameID, ".",
                    " The gene will be placed into the no evidence category.")
            
            ## If the predicted stop is downstream of the end of the region, skip
            ## subsequent genome stops between the correct stop and the predicted
            ## one as long as no predicted genes will be skipped.
            skipStop[checkNextStops(regionalPredictedStop)] <- TRUE
            
          } else if (outOfFrame) {
            warning("The predicted stop for gene ", geneIdx, ", position ",
                    fwdPredStop, ", is out of frame. The stop should be at ",
                    fwdRightStop, ". The gene is in frame ", frameID, ".",
                    " The gene will be placed into the no evidence category.")
            
          } else if (beforeEnd) {
            warning("The predicted stop for gene ", geneIdx, ", position ", fwdPredStop,
                    ", is upstream of the correct stop. The stop should be at ",
                    fwdRightStop, ". The gene is in frame ", frameID, ".",
                    " The gene will be placed into the no evidence category.")
          }
          
          ## Assign the gene to the no evidence category.
          catGeneAssignment[geneIdx] <- "Y CS- PE-"
          next()
        }
      }
      
      if (useProt) {
        ## Get the protein hits in that region.
        regionalProteinHits <- frameProteinHits[regRange] > 0
        regionalProteinHitStarts <- which(diff(c(0, regionalProteinHits)) == 1) # start
        
        ## If there are protein hits, determine the start of the first protein hit.
        if (length(regionalProteinHitStarts) >= 1) {
          firstProteinHitStart <- min(regRange[regionalProteinHitStarts])
          
          if ((isNTerm) && (length(regionalPredictedStart) > 0)) {
            protAligned <- (firstProteinHitStart == regionalPredictedStart) ||
              (firstProteinHitStart == (regionalPredictedStart + 3L))
          }
          
          noProt <- FALSE
        }
      }
      
      ## Skip the current region if there is no proteomics hit
      ## or gene start mapping to it.
      if ((noProt) && (length(regionalPredictedStart) <= 0)) {
        next()
      }
      
      ## If evolutionary conservation is being used in assessement,
      ## find the conserved starts in that region.
      if (useCons) {
        ## As long as the predicted start maps to a designated start codon, ...
        if ((length(regionalPredictedStart) <= 0) ||
            !(is.na(strandConservedStarts[regionalPredictedStart]))) {
          ## .... find which start positions have some conservation and coverage in the region.
          allConStartIdxs <- which((strandConservedStarts[regRange] / strandCoverage[regRange] > 0) &
                                     (strandCoverage[regRange] >= minCov) &
                                     (((regRange - cFrame) %% 3) == 0))
          
          if (length(allConStartIdxs) >= 1) {
            ## If there are any such positions, get all possible conserved starts.
            allConservedStarts <- regRange[allConStartIdxs]
            
            ## Only the conserved starts upstream of the protein hits can be used as evidence.
            if (!noProt) {
              conStarts <- allConservedStarts[(allConservedStarts <= firstProteinHitStart)]
            } else {
              conStarts <- allConservedStarts
            }
            
            ## If there are conserved starts (upstream of any protein hits), ...
            if (length(conStarts) >= 1) {
              ## ... get their scores.
              conStartScores <- strandConservedStarts[conStarts] / strandCoverage[conStarts]
              
              ## Check if there are any above the strong conserved start threshold.
              strongConStartIdxs <- which(conStartScores >= minConCovRatio_Strong)
              
              ## If at least one of the scores is above the threshold,
              ## there is at least one valid conserved start
              if (length(strongConStartIdxs) >= 1) {
                noConStarts <- FALSE
              }
            }
          }
        }
      }
      
      ## If evolutionary conservation is being used in assessement, find the conserved stops
      ## in that region. Requires a predicted gene in that region.
      if ((length(regionalPredictedStart) > 0) && (useCons)) {
        ## Find which positions have some coverage and stop codon conservation in the region.
        allConStopIdxs <- which((strandConservedStops[regRange] / strandCoverage[regRange] > minConCovRatio_Stop) &
                                  (strandCoverage[regRange] >= minCov) &
                                  (((regRange - cFrame) %% 3) == 0))
        
        if (length(allConStopIdxs) >= 1) {
          allConservedStops <- regRange[allConStopIdxs]
          
          ## Only the conserved stops upstream of the protein hits can be used as evidence.
          if (!noProt) {
            currConStops <- allConservedStops[(allConservedStops <= firstProteinHitStart)]
          } else {
            currConStops <- allConservedStops
          }
          
          geneLen <- (regEnd + 2L) - (regionalPredictedStart) + 1L
          
          ## Only the conserved stops early on in the gene can be used as evidence.
          condition1 <- (currConStops >= regionalPredictedStart)
          condition2 <- (currConStops < (regionalPredictedStart + (geneLen * noConStopsGeneFrac)))
          
          if (any(condition1 & condition2)) {
            ## If there is no protein evidence and there are valid conserved stop(s),
            ## put the gene in the potential false positive category.
            if (noProt) {
              ## At least one predicted start in that region
              ## Conserved stops disprove the gene --> "Y CS! PE-"
              catGeneAssignment[geneIdx] <- "Y CS! PE-"
              next()
            }
            
            ## If there is protein evidence and valid conserved stop(s) upstream,
            ## put the gene in the wrong start by conserved stop category.
            if (!noProt) {
              ## At least one predicted start in that region
              ## Protein evidence proves gene's existence
              ## Conserved stops disprove the gene start --> "Y CS! PE+"
              catGeneAssignment[geneIdx] <- "Y CS! PE+"
              next()
            }
          }
        }
      }
      
      ## If there is only a gene start mapping to the current region...
      if ((noProt) && (noConStarts) && (length(regionalPredictedStart) > 0)) {
        ## At least one predicted start in that region
        ## No protein hits and no conserved starts --> Category 7, "Y CS- PE-"
        catGeneAssignment[geneIdx] <- "Y CS- PE-"
        next()
      }
      
      ## Check if there are any protein hits within that region (section for categories 1-3).
      if (!noProt) {
        ## There are protein hits in that region.
        ## Check if there is any predicted start within that region.
        if (length(regionalPredictedStart) <= 0) {
          ## No predicted starts within that region.
          
          if ((regEnd - regStart) >= minMissORFLen) {
            ## Check if region is nested, if necessary.
            if ((allowNestedORFs) || !(isNestedORF())) {
              if (frameID <= 3) {
                cat1Pos <- c(regStart + 3L, regEnd + 2L)
              } else {
                cat1Pos <- genomeLength - c(regStart + 3L, regEnd + 2L) + 1L
              }
              
              ## Check if there is any conserved start (upstream of the first protein hit) in that region.
              if (noConStarts) {
                ## Protein hit with no conserved start upstream and no predicted start --> Category 1A
                ## "N CS- PE+"
                cat1A_Counter <- cat1A_Counter + 1L
                
                cat1A_ORFs[cat1A_Counter, ] <- c(cat1Pos, regEnd - regStart, frameID, checkOtherFrames())
              } else {
                ## Protein hit with conserved start upstream but no predicted start --> Category 1B
                ## "N CS< PE+"
                cat1B_Counter <- cat1B_Counter + 1L
                
                cat1B_ORFs[cat1B_Counter, ] <- c(cat1Pos, regEnd - regStart, frameID, checkOtherFrames())
              }
            }
          }
          
          next()
        } else {
          ## At least one predicted start in that region
          
          ## Check if there is any conserved start upstream of the first protein hit in that region.
          if (noConStarts) {
            ## No conserved starts upstream of the first protein hit in that region
            
            ## Check if the gene start is upstream or downstream of the first protein hit.
            if((isNTerm && protAligned) || ((!isNTerm) && (regionalPredictedStart <= firstProteinHitStart))) {
              ## Upstream
              ## Protein evidence with upstream gene start and no conserved start --> Category 2
              ## "Y CS- PE+"
              catGeneAssignment[geneIdx] <- "Y CS- PE+"
            } else {
              ## Downstream
              ## Protein evidence with downstream gene start and no conserved start --> Category 3A
              ## "Y CS- PE!"
              catGeneAssignment[geneIdx] <- "Y CS- PE!"
            }
            
            next()
          } else {
            ## At least one conserved start upstream of the first protein hit in that region
            
            ## Check if the gene start is downstream of the first protein hit in that region.
            if ((isNTerm && (!protAligned)) || ((!isNTerm) && (regionalPredictedStart > firstProteinHitStart))) {
              ## Protein evidence with downstream gene start
              ## and at least one conserved start upstream --> Category 3B
              ## "Y CS< PE!"
              catGeneAssignment[geneIdx] <- "Y CS< PE!"
              next()
            }
          }
        }
      }
      
      ## (Section for categories 4-6)
      ## There should be strong conserved start(s) and a given gene start by this point.
      
      ## Determine which of the all possible conserved starts in the region have bad scores.
      allConStartScores <- strandConservedStarts[allConservedStarts] / strandCoverage[allConservedStarts]
      badAllConStartIdxs <- which(allConStartScores < limConCovRatio_NotCon)
      
      ## Determine which of the conserved starts that are upstream of the protein evidence have bad scores.
      badConStartIdxs <- which(conStartScores < limConCovRatio_NotCon)
      
      ## Determine which of the conserved starts is the best one (or are the best ones).
      bestConStartIdx <- which(conStartScores == max(conStartScores))
      
      ## If there is only one conserved start with the best score,
      ## and that strongest conserved start is aligned with the given start,
      ## and there is at least one bad conserved start somewhere in the region, ...
      if ((length(bestConStartIdx) == 1) &&
          (conStarts[bestConStartIdx] == regionalPredictedStart) &&
          (length(badAllConStartIdxs) >= 1)) {
        ## ..., there is good reason to pick the aligned conserved start as evidence. --> CS+
        
        if (!noProt) {
          ## Protein evidence downstream --> Category 4A, "Y CS+ PE+"
          catGeneAssignment[geneIdx] <- "Y CS+ PE+"
        } else {
          ## No protein evidence --> Category 4B, "Y CS+ PE-"
          catGeneAssignment[geneIdx] <- "Y CS+ PE-"
        }
        
        next()
      }
      
      ## If there is a bad conserved start aligned with the given start, ... 
      if ((length(badConStartIdxs) >= 1) && (any(conStarts[badConStartIdxs] == regionalPredictedStart))) {
        ## The strongest conserved starts are somewhere in the region but not aligned with given start.
        ## Check where the most upstream of the strongest conserved starts is relative to the given start.
        pickUpstream <- (conStarts[bestConStartIdx[1]] < regionalPredictedStart)
        
        if (!pickUpstream) {
          ## All of the strongest conserved starts are downstream of the gene start.
          if (!noProt) {
            ## Protein evidence downstream --> Category 5A, "Y CS> PE+"
            catGeneAssignment[geneIdx] <- "Y CS> PE+"
          } else {
            ## No protein evidence --> Category 5B, "Y CS> PE-"
            catGeneAssignment[geneIdx] <- "Y CS> PE-"
          }
        } else {
          ## At least one of the strongest conserved starts is upstream of the gene start.
          if (!noProt) {
            ## Protein evidence downstream --> Category 6A, "Y CS< PE+"
            catGeneAssignment[geneIdx] <- "Y CS< PE+"
          } else {
            ## No protein evidence --> Category 6B, "Y CS< PE-"
            catGeneAssignment[geneIdx] <- "Y CS< PE-"
          }
        }
        
        next()
      }
      
      ## A good decision cannot be made on conserved starts. --> CS-
      if (!noProt) {
        ## Protein evidence with upstream gene start and no valid conserved start
        ## --> Category 2, "Y CS- PE+"
        catGeneAssignment[geneIdx] <- "Y CS- PE+"
      } else {
        ## No protein hits and no valid conserved starts --> Category 7, "Y CS- PE-"
        catGeneAssignment[geneIdx] <- "Y CS- PE-"
      }
    }
    
    
    if (verbose) {
      setTxtProgressBar(pBar, frameID / 6)
      
      if (frameID == 6) {
        cat("\n")
      }
    }
  }
  
  
  ## --------------------------------------------------------------------------------------------------------------- ##
  
  if (cat1A_Counter <= 0) {
    cat1A_ORFs <- matrix(0, nrow = 0, ncol = 5,
                         dimnames = list(NULL,
                                         c("Start", "End", "Length", "Frame", "OtherProtFrame")))
  } else {
    cat1A_ORFs <- cat1A_ORFs[seq_len(cat1A_Counter), , drop = FALSE]
  }
  
  if (cat1B_Counter <= 0) {
    cat1B_ORFs <- matrix(0, nrow = 0, ncol = 5,
                         dimnames = list(NULL,
                                         c("Start", "End", "Length", "Frame", "OtherProtFrame")))
  } else {
    cat1B_ORFs <- cat1B_ORFs[seq_len(cat1B_Counter), , drop = FALSE]
  }
  
  return(structure(list("StrainID" = mapObj$StrainID,
                        "Species" = mapObj$Species,
                        "GenomeLength" = genomeLength,
                        "GeneLeftPos" = geneLeftPos,
                        "GeneRightPos" = geneRightPos,
                        "GeneStrand" = geneStrand,
                        "GeneSource" = geneSource,
                        "NumGenes" = length(geneLeftPos),
                        "N_CS-_PE+_ORFs" = cat1A_ORFs,
                        "N_CS<_PE+_ORFs" = cat1B_ORFs,
                        "CategoryAssignments" = catGeneAssignment),
                   class = c("Assessment", "Results")))
}