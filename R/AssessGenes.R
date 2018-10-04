#' @export
#'
#' @title Assess Genes
#' @description Assess and categorize a set of genes for a genome using proteomics hits, evolutionarily conserved starts,
#' and evolutionarily conserved stops as evidence
#'
#' @usage
#' AssessGenes(geneLeftPos,
#'             geneRightPos,
#'             geneStrand,
#'             inputMapObj,
#'             geneSource = "",
#'             minCovNum = 10,
#'             minCovPct = 5,
#'             minConCovRatio_Best = 0.99,
#'             limConCovRatio_NotCon = 0.8,
#'             minConCovRatio_Stop = 0.5,
#'             noConStopsGeneFrac = 0.5,
#'             minNumStops = 2,
#'             minMissORFLen = 0,
#'             allowNestedORFs = FALSE,
#'             useNTermProt = FALSE,
#'             verbose = TRUE)
#'
#' @param geneLeftPos An integer vector with the left positions of each gene, in terms of the forward strand.
#'
#' @param geneRightPos An integer vector with the right positions of each gene, in terms of the forward strand.
#'
#' @param geneStrand A character vector consisting of "+" and "-", specifying which strand each gene is on.
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
#' @param minConCovRatio_Best Minimum value of the start codon conservation to coverage ratio needed to call a start conserved.
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
#' gene and moving towards the stop of gene, to use in searching for conserved stops. For example, a value of 0.25 means that the first
#' quarter of the gene is checked for conserved stops, a value of 0.5 correspond to the first half of the gene, etc. Recommended to use
#' the default value.
#' 
#' @param minNumStops Minimum number of conserved stop positions required to be within a gene (with no mapped proteomics hits)
#' in order to categorize that gene as an overprediction. Recommended to use the default value.
#' 
#' @param minMissORFLen Minimum ORF length required to include an ORF with protein hits but no gene start in the final results.
#' 
#' @param allowNestedORFs Logical indicating whether or not to include ORFs with protein hits but no gene starts that are
#' completely nested within an ORF in another frame in the final results.
#'
#' @param verbose Logical indicating whether or not to display progress and status messages.
#' 
#' @param useNTermProt Logical indicating whether or not to require the given gene starts to align with (or be one codon off from)
#' the start of the first proteomics hit in the ORF. The mapping object must be built with N-terminal proteomics data.
#' Default value is FALSE.
#'
#' @details
#' For each of the given genes, \code{AssessGene} assigns a category based on where the conserved starts, conserved stops, and
#' proteomics hits are located in relation to the start of the gene. The category assignments for the genes are stored in the
#' \code{CategoryAssignments} vector in the \code{Results} object returned by the function. This vector has the same length as
#' and aligns with the indexing of the three given gene positional information vectors. Please see
#' \code{\link{Assessment-class}} for a list of all possible categories and their descriptions.
#' 
#' The vectors \code{geneLeftPos}, \code{geneRightPos}, and \code{geneStrand} must all be of the same length. Using the same
#' index with each vector must provide information on the same gene (think of the vectors as columns of the same table).
#'
#' \code{geneLeftPos} and \code{geneRightPos} describe the upstream and downstream positions (respectively) for each gene in
#' terms of the forward strand. For genes on the forward strand, \code{geneLeftPos} corresponds to the start positions and
#' \code{geneRightPos} corresponds to stop positions. For genes on the reverse strand, \code{geneLeftPos} corresponds to the
#' stop positions and \code{geneRightPos} corresponds to the start positions. Gene positions on the reverse strand must be
#' relative to the 5' to 3' direction of the forward strand (as opposed to being relative to the 5' to 3' direction of the
#' reverse strand). This means that none of the elements of \code{geneLeftPos} can be greater than (or equal to) the
#' corresponding element in \code{geneRightPos}.
#'
#' Please ensure that the same genome that is used in the mapping function is also used to derive the set of genes for this
#' assessment function. The function will only error if any gene positions are outside the bounds of the genome.
#'
#' The maximum of either \code{minCovNum} or (\code{minCovPct} / 100) multiplied by the number of related genomes is used as
#' the minimum coverage required in determining conserved starts and stops.
#' 
#' Additionaly, open reading frames with proteomics evidence but no gene start are categorized based on whether or not there
#' is a conserved start upstream of the proteomic evidence. The positions and lengths of these open reading frames is included
#' in the \code{N_CS-_PE+_ORFs} and \code{N_CS+_PE+_ORFs} matrices within the final object that is returned.
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
AssessGenes <- function(geneLeftPos,
                        geneRightPos,
                        geneStrand,
                        inputMapObj,
                        geneSource = "",
                        minCovNum = 10L,
                        minCovPct = 5L,
                        minConCovRatio_Best = 0.99,
                        limConCovRatio_NotCon = 0.8,
                        minConCovRatio_Stop = 0.5,
                        noConStopsGeneFrac = 0.5,
                        minNumStops = 2L,
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
    if (!("package:AssessORFData" %in% search())) {
      if (!(any(installed.packages()[, 1] == "AssessORFData"))) {
        stop("The package 'AssessORFData' has not been installed.")
      }
      stop("The package 'AssessORFData' has not been loaded.")
    }
    
    mapObj <- GetDataMapObj(inputMapObj)
    
  } else {
    stop("'inputMapObj' must either be an object of class 'Assessment', subclass 'DataMap' ",
         "or be a single character string.")
  }
  
  ## Get the genome length to use for checking gene boundaries.
  genomeLength <- mapObj$GenomeLength
  
  ## --------------------------------------------------------------------------------------------------------------- ##
  
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
  
  if ((!is.numeric(minConCovRatio_Best)) || (anyNA(minConCovRatio_Best))) {
    stop("'minConCovRatio_Best' must be a valid real number.")
  }
  
  if (length(minConCovRatio_Best) != 1) {
    stop("'minConCovRatio_Best' must consist of exactly one number.")
  }
  
  if ((minConCovRatio_Best <= 0) || (minConCovRatio_Best >= 1)) {
    stop("'minConCovRatio_Best' must be greater than 0 and less than 1.")
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
  
  if ((!is.numeric(minNumStops))  || (anyNA(minNumStops))) {
    stop("'minNumStops must be a valid real number.")
  }
  
  if (length(minNumStops) != 1) {
    stop("Exactly one number must be inputted for 'minNumStops'.")
  }
  
  if ((minNumStops %% 1 != 0) || (minNumStops <= 0)) {
    stop("'minNumStops' must be a non-negative integer that is greater than 0.")
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
  
  if ((useCons) && (minCovNum > mapObj$NumRelatedGenomes)) {
    useCons <- FALSE
    
    warning("The value for the minimum number of related genomes required is greater than the number of related",
            " genomes used to generate the map object. Conserved starts & stops will not be used in assessment.")
  }
  
  if (useCons) {
    minCov <- max(minCovNum, minCovPct * mapObj$NumRelatedGenomes / 100)
  }
  
  ## --------------------------------------------------------------------------------------------------------------- ##
  
  checkRange <- function(x) {
    if (all(x < 0))
      stop("X-axis out of bounds.")
    if (all(x > genomeLength))
      stop("X-axis out of bounds.")
    if (x[1] < 0)
      x[1] <- 0
    if (x[2] > genomeLength)
      x[2] <- genomeLength
    x <- floor(x)
    return(x[1]:x[2])
  }
  
  ## This function checks if there is overlapping protein evidence in other frames. It is used
  ## when recording information for categories with protein evidence and no gene start.
  checkOtherFrames <- function() {
    currProtHitPos <- orfRange[regionalProteinHits]
    
    sameRange <- (min(currProtHitPos)):(max(currProtHitPos))
    
    oppRange <- genomeLength - sameRange + 1L
    
    for (oFrame in (1:6)[-frameID]) {
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
          oppStrandProt <- revProt[[oFrame - 3]][[1]][orfRange]
        }
      }
      
      if (any(oppStrandProt > 0)) {
        return(oFrame)
      }
    }
    
    return(0L)
  }
  
  ## This function checks if the current ORF is completely within another ORF. It is used
  ## to filter out some of the ORFs with protein hits and no given start.
  isNestedORF <- function() {
    currProtHitPos <- orfRange[regionalProteinHits]
    
    sameLeftPos <- min(currProtHitPos)
    sameRightPos <- max(currProtHitPos)
    
    oppLeftPos <- genomeLength - sameRightPos + 1L
    oppRightPos <- genomeLength - sameLeftPos + 1L
    
    for (oFrame in (1:6)[-frameID]) {
      
      if (((frameID <= 3) && (oFrame <= 3)) || ((frameID > 3) && (oFrame > 3))) {
        ## Use same
        nested <- any((stops[[oFrame]][-length(stops[[oFrame]])] <= sameLeftPos) &
                        (stops[[oFrame]][-1] >= sameRightPos))
      } else {
        ## USe opposite
        nested <- any((stops[[oFrame]][-length(stops[[oFrame]])] <= oppLeftPos) &
                        (stops[[oFrame]][-1] >= oppRightPos))
      }
      
      if (nested) {
        return(TRUE)
      }
    }
    
    return(FALSE)
  }
  
  cat1A_Counter <- 0L
  
  cat1A_ORFs <- matrix(0, nrow = sum(sapply(stops, length)), ncol = 5,
                       dimnames = list(NULL,
                                       c("Start", "End", "Length", "Frame", "OtherProtFrame")))
  
  cat1B_Counter <- 0L
  
  cat1B_ORFs <- matrix(0, nrow = sum(sapply(stops, length)), ncol = 5,
                       dimnames = list(NULL,
                                       c("Start", "End", "Length", "Frame", "OtherProtFrame")))
  
  catGeneAssignment <- character(length(geneLeftPos))
  
  if (verbose) {
    pBar <- txtProgressBar(style=ifelse(interactive(), 3, 1))
  }
  
  for (frameID in 1:6) {
    if (frameID <= 3) {
      strandID <- "+"
      cFrame <- frameID
      shift <- cFrame %% 3
      
      frameProteinHits <- fwdProt[[cFrame]][[1]]
      
      strandCoverage <- fwdCov
      strandConservedStarts <- fwdConStart
      strandConservedStops <- fwdConStop
      
      strandOnly <- which(geneStrand == strandID)
      
      strandStarts <- geneLeftPos[strandOnly]
      names(strandStarts) <- as.character(strandOnly)
      
      strandFrameStarts <- strandStarts[(strandStarts %% 3 == shift)]
      
    } else {
      strandID <- "-"
      cFrame <- frameID - 3
      shift <- cFrame %% 3
      
      frameProteinHits <- revProt[[cFrame]][[1]]
      
      strandCoverage <- revCov
      strandConservedStarts <- revConStart
      strandConservedStops <- revConStop
      
      strandOnly <- which(geneStrand == strandID)
      
      strandStarts <- genomeLength - geneRightPos[strandOnly] + 1
      names(strandStarts) <- as.character(strandOnly)
      
      strandFrameStarts <- strandStarts[(strandStarts %% 3 == shift)]
    }
    
    ## In order to handle genes at the boundary of the genome, two additional
    ## in-frame positions need to be added to the stops vector for the frame:
    ## the start of the first possible in-frame codon and the start of the
    ## last possible in-frame codon.
    
    finalCodonStart <- ifelse(genomeLength %% 3 == shift, genomeLength,
                              ifelse((genomeLength - 1L) %% 3 == shift,
                                     genomeLength - 1L, genomeLength - 2L))
    
    frameStops <- c(cFrame, stops[[frameID]], finalCodonStart)
    
    for (idx in seq(1, length(frameStops) - 1)) {
      ## Set boundaries of the region to examine.
      orfStart <- frameStops[idx]
      orfEnd <- frameStops[idx + 1]
      orfRange <- checkRange(c(orfStart, orfEnd))
      
      ## For each ORF, start by assuming there are no valid conserved starts,
      ## and no proteomics hits mapping to that ORF.
      noConStarts <- TRUE
      noProt <- TRUE
      
      ## Get the gene start(s) in that region.
      regionalPredictedStart <-strandFrameStarts[((strandFrameStarts >= orfStart) &
                                                    (strandFrameStarts <= orfEnd))]
      geneIdx <- as.integer(names(regionalPredictedStart))
      
      if (length(regionalPredictedStart) > 1) {
        warning("More than one in-frame gene start between stops ",
                orfStart, " and ", orfEnd, " in frame ", frameID, ".",
                paste(" A gene start is at ", regionalPredictedStart,
                      " (gene ", geneIdx, ").", sep = ""))
        
        catGeneAssignment[geneIdx] <- "Y CS- PE-"
        next()
      }
      
      if (useProt) {
        ## Get the protein hits in that region.
        regionalProteinHits <- frameProteinHits[orfRange] > 0
        regionalProteinHitStarts <- which(diff(c(0, regionalProteinHits)) == 1) # start
        
        ## If there are protein hits, determine the start of the first protein hit.
        if (length(regionalProteinHitStarts) >= 1) {
          firstProteinHitStart <- min(orfRange[regionalProteinHitStarts])
          
          if ((isNTerm) && (length(regionalPredictedStart) > 0)) {
            protAligned <- (firstProteinHitStart == regionalPredictedStart) ||
              (firstProteinHitStart == (regionalPredictedStart + 3L))
          }
          
          noProt <- FALSE
        }
      }
      
      ## If evolutionary conservation is being used in assessement,
      ## find the conserved starts in that region.
      if (useCons) {
        ## As long as the predicted start maps to a designated start codon, ...
        if ((length(regionalPredictedStart) <= 0) ||
            !(is.na(strandConservedStarts[regionalPredictedStart]))) {
          ## .... find which start positions have some conservation and coverage in the region.
          allConStartIdxs <- which((strandConservedStarts[orfRange] / strandCoverage[orfRange] > 0) &
                                     (strandCoverage[orfRange] >= minCov) &
                                     (((orfRange - cFrame) %% 3) == 0))
          
          if (length(allConStartIdxs) >= 1) {
            ## If there are any such positions, get all possible conserved starts.
            allConservedStarts <- orfRange[allConStartIdxs]
            
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
              
              ## Check if there are any above the conserved start threshold.
              bestConStartIdxs <- which(conStartScores >= minConCovRatio_Best)
              
              if (length(bestConStartIdxs) >= 1) {
                ## If so, there is at least one valid conserved start.
                noConStarts <- FALSE
              }
            }
          }
        }
      }
      
      ## Skip the current ORF if there is no proteomics hit or gene start mapping to it.
      if ((noProt) && (length(regionalPredictedStart) <= 0)) {
        next()
      }
      
      ## If evolutionary conservation is being used in assessement, find the conserved stops
      ## in that region. Requires a predicted gene in that region.
      if ((length(regionalPredictedStart) > 0) && (useCons)) {
        ## Find which positions have some coverage and stop codon conservation in the region.
        allConStopIdxs <- which((strandConservedStops[orfRange] / strandCoverage[orfRange] > minConCovRatio_Stop) &
                                  (strandCoverage[orfRange] >= minCov) &
                                  (((orfRange - cFrame) %% 3) == 0))
        
        if (length(allConStopIdxs) >= 1) {
          allConservedStops <- orfRange[allConStopIdxs]
          
          ## Only the conserved stops upstream of the protein hits can be used as evidence.
          if (!noProt) {
            currConStops <- allConservedStops[(allConservedStops <= firstProteinHitStart)]
          } else {
            currConStops <- allConservedStops
          }
          
          geneLen <- (orfEnd + 2L) - (regionalPredictedStart) + 1L
          
          ## Only the conserved stops early on in the gene can be used as evidence.
          condition1 <- (currConStops >= regionalPredictedStart)
          condition2 <- (currConStops < (regionalPredictedStart + (geneLen * noConStopsGeneFrac)))
          
          if (any(condition1 & condition2)) {
            validConStops <- currConStops[(condition1 & condition2)]
            
            ## Get the most downstream conserved stop to check if there is a conserved start downstream.
            mostDownstreamStop <- max(validConStops)
            
            multiStops <- (length(validConStops) >= minNumStops)
            
            conStopConStart <- (any(mostDownstreamStop < conStarts[bestConStartIdxs]))
            
            ## If there is no protein evidence and either there are (multiple) valid conserved stops or there is
            ## a strong conserved start downstream of the stops, put the gene in the conserved stop category.
            if ((noProt) && (multiStops || conStopConStart)) {
              ## At least one Prodigal start in that region
              ## Conserved stops disprove the gene --> "Y CS! PE-"
              catGeneAssignment[geneIdx] <- "Y CS! PE-"
              next()
            }
            
            ## If there is protein evidence and there is a strong conserved start downstream of the stops,
            ## put the gene in the conserved stop category.
            if ((!noProt) && conStopConStart) {
              ## At least one Prodigal start in that region
              ## Conserved stops disprove the gene --> "Y CS! PE+"
              catGeneAssignment[geneIdx] <- "Y CS! PE+"
              next()
            }
          }
        }
      }
      
      ## If there is only a gene start mapping to the current ORF...
      if ((noProt) && (noConStarts) && (length(regionalPredictedStart) > 0)) {
        ## At least one Prodigal start in that region
        ## No protein hits and no conserved starts --> Category 7, "Y CS- PE-"
        catGeneAssignment[geneIdx] <- "Y CS- PE-"
        next()
      }
      
      ## Check if there are any protein hits within that region (section for categories 1-3).
      if (!noProt) {
        ## There are protein hits in that region.
        ## Check if there is any Prodigal start within that region.
        if (length(regionalPredictedStart) <= 0) {
          ## No prodigal starts within that region.
          
          if ((orfEnd - orfStart) >= minMissORFLen) {
            ## Check if ORF is nested, if necessary.
            if ((allowNestedORFs) || !(isNestedORF())) {
              if (frameID <= 3) {
                cat1Pos <- c(orfStart + 3L, orfEnd + 2L)
              } else {
                cat1Pos <- genomeLength - c(orfStart + 3L, orfEnd + 2L) + 1L
              }
              
              ## Check if there is any conserved start (upstream of the first protein hit) in that region.
              if (noConStarts) {
                ## Protein hit with no conserved start upstream and no Prodigal start --> Category 1A
                ## "N CS- PE+"
                cat1A_Counter <- cat1A_Counter + 1L
                
                cat1A_ORFs[cat1A_Counter, ] <- c(cat1Pos, orfEnd - orfStart, frameID, checkOtherFrames())
              } else {
                ## Protein hit with conserved start upstream but no Prodigal start --> Category 1B
                ## "N CS< PE+"
                cat1B_Counter <- cat1B_Counter + 1L
                
                cat1B_ORFs[cat1B_Counter, ] <- c(cat1Pos, orfEnd - orfStart, frameID, checkOtherFrames())
              }
            }
          }
          
          next()
        } else {
          ## At least one Prodigal start in that region
          
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
      ## There should definitely be good conserved start(s) and a given gene start by this point.
      
      ## Determine which of the all possible conserved starts in the region have bad scores.
      allConStartScores <- strandConservedStarts[allConservedStarts] / strandCoverage[allConservedStarts]
      badAllConStartIdxs <- which(allConStartScores < limConCovRatio_NotCon)
      
      ## Determine which of the conserved starts that are upstream of the protein evidence have bad scores.
      badConStartIdxs <- which(conStartScores < limConCovRatio_NotCon)
      
      
      ## If there is a good conserved start aligned with the given start,
      ##  and there is at least one bad conserved start, ...
      if ((any(conStarts[bestConStartIdxs] == regionalPredictedStart)) && (length(badAllConStartIdxs) >= 1)) {
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
      
      ## There is at least one good conserved start somewhere in the ORF but it is not aligned with given start.
      ## If there is a bad conserved start aligned with the given start, ...
      if ((length(badConStartIdxs) >= 1) && (any(conStarts[badConStartIdxs] == regionalPredictedStart))) {
        ## ... check where the best conserved starts are relative to the given gene start.
        pickDownstream <- (all(conStarts[bestConStartIdxs] > regionalPredictedStart))
        
        if (pickDownstream) {
          ## All the best nearby conserved starts are downstream of the gene start.
          if (!noProt) {
            ## Protein evidence downstream --> Category 5A, "Y CS> PE+"
            catGeneAssignment[geneIdx] <- "Y CS> PE+"
          } else {
            ## No protein evidence --> Category 5B, "Y CS> PE-"
            catGeneAssignment[geneIdx] <- "Y CS> PE-"
          }
        } else {
          ## At least one of the best nearby conserved starts is upstream of the gene start.
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
    cat1A_ORFs <- cat1A_ORFs[1:cat1A_Counter, , drop = FALSE]
  }
  
  if (cat1B_Counter <= 0) {
    cat1B_ORFs <- matrix(0, nrow = 0, ncol = 5,
                         dimnames = list(NULL,
                                         c("Start", "End", "Length", "Frame", "OtherProtFrame")))
  } else {
    cat1B_ORFs <- cat1B_ORFs[1:cat1B_Counter, , drop = FALSE]
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