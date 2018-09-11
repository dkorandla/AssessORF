#' @export
#' @import DECIPHER
#'
#' @title Map Evidence to a Genome
#' @description Maps proteomics hits and evolutionarily conserved starts to a central genome
#'
#' @usage
#' MapAssessmentData(
#' genomes_DBFile,
#' tblName = "Seqs",
#' central_ID,
#' related_IDs,
#' protHits_Seqs,
#' protHits_Scores = rep.int(1, length(protHits_Seqs)),
#' strainID = "",
#' speciesName = "",
#' protHits_Threshold = 0,
#' protHits_IsNTerm = FALSE,
#' related_KMerLen = 8,
#' related_MinDist = 0.01,
#' related_MaxDistantN = 1000,
#' startCodons = c("ATG", "GTG", "TTG"),
#' useProt = TRUE,
#' useCons = TRUE,
#' verbose = TRUE)
#'
#' @param genomes_DBFile A SQLite connection object or a character string specifying the path to the database file.
#'
#' @param tblName Character string specifying the table where the genome sequences are located.
#'
#' @param central_ID Character string specifying which identifier corresponds to the central genome, the genome
#' to which the proteomics data and evolutionary conservation data will be mapped.
#'
#' @param related_IDs Character vector of strings specifying identifiers that correspond to related genomes, the genomes
#' that will be used to determine which start codons (ATG, GTG, and TTG) are evolutionarily conserved.
#'
#' @param protHits_Seqs Character vector of amino acid strings that correspond to the sequences for the proteomics hits.
#'
#' @param protHits_Scores Numeric vector of (confidence) scores for the proteomics hits. Scores cannot be negative.
#' The default option assigns a score of one to each proteomics hit.
#'
#' @param strainID Optional character string that specifies the strain identifier that the central genome corresponds to.
#'
#' @param speciesName Optional character string that specifies the name of the species that the central genome corresponds to.
#'
#' @param protHits_Threshold Optional number that specifies what percent of the lowest scoring proteomics hits should be dropped.
#' Must be a non-negative integer less than 100.
#' 
#' @param protHits_IsNTerm Logical describing whether or not the proteomics hits come from N-terminal proteomics.
#' Default value is false.
#'
#' @param related_KMerLen The k-mer length to be used when measuring distances between the central genome and related genomes.
#' Default value is 8.
#' 
#' @param related_MinDist The minimum fractional distance required for a related genome to be used in finding
#' evolutionary conservation. Default value is 0.01.
#'
#' @param related_MaxDistantN The maximum number of related genomes to use in finding evolutionary conservation after the
#' related genomes have been sorted from most distantly related to most closely related in relation to the central genome.
#' Default value is 1000.
#' 
#' @param startCodons A charcter vector consisting of three-letter DNA strings to use as the start codons when finding
#' evolutionarily conserved starts.
#' 
#' @param useProt Logical indicating whether or not proteomics evidence should be mapped to the genome.
#' Default value is true. Cannot be false if \code{useCons} is false.
#' 
#' @param useCons Logical indicating whether or not evolutionary conservation evidence should be mapped to the genome.
#' Default value is true. Cannot be false if \code{useProt} is false.
#'
#' @param verbose Logical indicating whether or not to display progress and status messages.
#'
#' @details
#' All genomes used inside this function, including the central genome, must be inside the specified table of the specified
#' database. If the central genome is not found, the function returns an error. Please see the Using AssessORF vignette
#' for details on how to populate a database with genomic sequences.
#'
#' Information on the proteomics hits is primarily given by \code{protHits_Seqs} and \code{protHits_Scores}. The sequences
#' (\code{protHits_Seqs}) are mapped to the six-frame translations of the central genome, and the scores (\code{protHits_Scores})
#' are used in thresholding and plotting the proteomics hits.
#'
#' \code{protHits_Scores} can be a single number. In that case, that number is used the as the score for all proteomics hits.
#' Otherwise, the \code{protHits_Scores} must be of the same length as \code{protHits_Seqs}.
#'
#' Only proteomics hits with a score greater than the value of the percentile that corresponds to the value of \code{protHits_Threshold}
#' will be kept and the rest of the hits will be dropped. If all the proteomics hits have the same score or if \code{protHits_Threshold}
#' is zero, no thresholding will occur and no hits will be dropped.
#' 
#' Please note that \code{protHits_IsNTerm} has no affect on how the proteomics evidence is mapped to the central genome but it can be
#' used to affect how genes are assessed and categorized in \code{\link{AssessGenes}}.
#'
#' Evolutionarily conserved starts and conserved stop are found by first measuring how far the related genomes are from the central
#' genome using k-mer frequencies. Next, the most distant related genomes are aligned to the central genome. This provides information
#' on how often each position in the central genome is covered by syntenic matches to related genomes (coverage), how often those
#' positions correspond to the start codons (start codon conservation) in both genomes, and how often those positions correspond to stop
#' codons in related genomes (stop codon conservation). A ratio of conservation to coverage is used in downstream functions to measure
#' the strength of both conserved starts and conserved stops.
#'
#' Related genomes should be genomes from species that are closely related to the given strain. \code{related_IDs} specifies the
#' identifiers for the sequences of the related genomes inside the database. Any related genome identifier (each element of
#' \code{related_IDs}) is considered invalid and not used when finding evolutionary conservation if it is not found in the
#' databse. Please note that the function will only error when none of the related genomes are found. \code{related_MinDist} is used
#' to prevent the inclusion of related genomes that are too similar to the central genome. It is recommended to keep
#' \code{related_MinDist} at the default value. If there are less valid related genomes in the sequence database than value of
#' \code{related_MaxDistantN}, all related genomes will be used in finding evolutionary conservation.
#' 
#' The logical flag \code{useProt} is used to indicate whether or not proteomics evidence has been provided and should be mapped to
#' the genome. Error checking will not occur for any arguments that involve proteomics if it is false.
#' 
#' The logical flag \code{useCons} is used to indicate whether or not evolutionary conservation evidence has been provided and should be
#' mapped to the genome. Error checking will not occur for any arguments that involve evolutionary conservation if it is false. 
#'
#' @return An object of class \code{Assessment} and subclass \code{DataMap}
#' 
#' @seealso \code{\link{Assessment-class}}
#'
#'@examples
#'
#' ## Example showing the minimum number of arguments that need to be specified
#' ## to map both proteomics and evolutionary conservation data:
#'
#' \dontrun{
#' myMapObj <- MapAssessmentData(myDBFile, central_ID = "1",
#'                               related_IDs = as.character(2:1001),
#'                               protHits_Seqs = myProtSeqs)
#' }
#' 
#' 
#' ## Runnable example that uses evolutionary conservation data only:
#' ## Human adenovirus 1 is the strain of interest, and the set of Adenoviridae
#' ## genomes will serve as the set of genome. The cenral genome, also known as
#' ## the genome of human adenovirus 1, is at identifier 1. The related genomes
#' ## are at identifiers 2 - 105.
#' 
#' myMapObj <- MapAssessmentData(system.file("extdata",
#'                                           "Adenoviridae.sqlite",
#'                                           package = "AssessORF"),
#'                               central_ID = "1",
#'                               related_IDs = as.character(2:53),
#'                               speciesName = "Human adenovirus 1",
#'                               useProt = FALSE)
#' 
#'
MapAssessmentData <- function(genomes_DBFile,
                              tblName = "Seqs",
                              central_ID,
                              related_IDs,
                              protHits_Seqs,
                              protHits_Scores = rep.int(1L, length(protHits_Seqs)),
                              strainID = "",
                              speciesName = "",
                              protHits_Threshold = 0L,
                              protHits_IsNTerm = FALSE,
                              related_KMerLen = 8L,
                              related_MinDist = 0.01,
                              related_MaxDistantN = 1000L,
                              startCodons = c("ATG", "GTG", "TTG"),
                              useProt = TRUE,
                              useCons = TRUE,
                              verbose = TRUE) {
  
  ## Check inputs for error.
  
  if ((!is.logical(verbose)) ||(anyNA(verbose)) || (length(verbose) != 1)) {
    stop("'verbose' must be of type logical, be either TRUE or FALSE, and consist of only 1 element.")
  }
  
  if ((!is.logical(useProt)) ||(anyNA(useProt)) || (length(useProt) != 1)) {
    stop("'useProt' must be of type logical, be either TRUE or FALSE, and consist of only 1 element.")
  }
  
  if ((!is.logical(useCons)) ||(anyNA(useCons)) || (length(useCons) != 1)) {
    stop("'useCons' must be of type logical, be either TRUE or FALSE, and consist of only 1 element.")
  }
  
  if ((!useProt) && (!useCons)) {
    stop("Either 'useProt' or 'useCons' must be true. Must specify at least one type of evidence.")
  }
  
  if ((!is.character(central_ID)) || (anyNA(central_ID))) {
    stop("The ID for the central genome must be a valid character string.")
  }
  
  if (length(central_ID) != 1) {
    stop("Exactly one character string must be inputted as the ID for the central genome.")
  }
  
  if (central_ID == "") {
    stop("The ID for the central genome cannot be an empty character string.")
  }
  
  if (useProt) {
    if (!is.character(protHits_Seqs)) {
      stop("Sequences for the proteomics hits must be inputted as a character vector.")
    }
    
    if (length(protHits_Seqs) <= 0) {
      stop("'protHits_Seqs' must have a length greater than 0.")
    }
    
    if ((anyNA(protHits_Seqs)) || (any(protHits_Seqs == ""))) {
      stop("Sequences for the proteomics hits must be inputted as valid, non-empty character strings.")
    }
    
    ## --------------------------------------------------------------------------------------------------------------- ##
    
    if ((!is.numeric(protHits_Scores)) || (anyNA(protHits_Scores))) {
      stop("Scores for the proteomics hits must be inputted as valid numbers.")
    }
    
    thresholdProt <- FALSE
    
    if (length(protHits_Scores) <= 0) {
      
      warning("'protHits_Scores' does not have a length greater than 0. ",
              "Assigning a score of 1 to all proteomics hits. ",
              "No threshold will be applied.")
      
      protHits_Scores <- rep.int(1L, length(protHits_Seqs))
      
    } else if ((length(protHits_Seqs) != length(protHits_Scores)) && ((length(protHits_Score) != 1))) {
      
      stop("The number of proteomic hit sequences does not match up with the number of proteomic hit scores. ",
           "The number of scores must either be one or the same as the number of sequences.")
      
    } else if (all(protHits_Scores == 0)) {
      warning("'protHits_Scores' consists entirely of zeroes. ",
              "Assigning a score of 1 to all proteomics hits. ",
              "No threshold will be applied.\n")
      
      protHits_Scores <- rep.int(1L, length(protHits_Seqs))
      
    } else if (any(protHits_Scores <= 0)) {
      stop("Scores for the proteomics hits must all be positive.")
    } else if (length(protHits_Seqs) == length(protHits_Scores)) {
      
      if (any(protHits_Scores != protHits_Scores[1])) {
        thresholdProt <- TRUE
        
        if (verbose) {
          cat("'protHits_Scores' has distinct scores.",
              "A threshold will be applied (if possible).\n")
        }
      }
      
    } else {
      ## 'protHits_Score' has a length of 1.
      
      if (verbose) {
        cat("'protHits_Scores' has length of 1.",
            "Assigning the same score to all proteomics hits.",
            "No threshold will be applied.\n")
      }
      
      protHits_Scores <- rep(protHits_Scores, length(protHits_Seqs))
    }
    
    ## --------------------------------------------------------------------------------------------------------------- ##
    
    if (thresholdProt) {
      if ((!is.numeric(protHits_Threshold)) || (anyNA(protHits_Threshold))) {
        stop("The threshold for proteomic hits must be a valid real number.")
      }
      
      if (length(protHits_Threshold) != 1) {
        stop("Exactly one number must be inputted as the threshold for proteomic hits.")
      }
      
      if (protHits_Threshold == 0) {
        
        if (verbose) {
          cat("'protHits_Threshold' is zero so it is not possible to apply a threshold.\n")
        }
        
        thresholdProt <- FALSE
        
      } else if ((protHits_Threshold %% 1 != 0) || (protHits_Threshold < 0) || (protHits_Threshold >= 100)) {
        stop("The threshold for proteomic hits must be a non-negative integer ",
             "that is greater than (or equal to) 0 and less than 100.")
      }
    }
    
    ## --------------------------------------------------------------------------------------------------------------- ##
    
    if ((!is.logical(protHits_IsNTerm)) ||(anyNA(protHits_IsNTerm)) || (length(protHits_IsNTerm) != 1)) {
      stop("'protHits_IsNTerm' must be of type logical, be either TRUE or FALSE, and consist of only 1 element.")
    }
  }
  
  ## --------------------------------------------------------------------------------------------------------------- ##
  
  if (useCons) {
    if ((!is.character(related_IDs)) || (anyNA(related_IDs))) {
      stop("The IDs for the related genomes must be valid character strings.")
    }
    
    if (length(related_IDs) <= 0) {
      stop("One or more character strings must be inputted as the IDs for the related genomes.")
    }
    
    if (any(related_IDs == "")) {
      stop("None of the IDs for the related genomes can be empty character strings.")
    }
    
    if ((!is.numeric(related_KMerLen))  || (anyNA(related_KMerLen))) {
      stop("The length of k-mers to use in distance measuring must be a valid real number.")
    }
    
    if (length(related_KMerLen) != 1) {
      stop("Exactly one number must be inputted as the length of k-mers to use in distance measuring.")
    }
    
    if ((related_KMerLen %% 1 != 0) || (related_KMerLen <= 0) ||(related_KMerLen >= 11)) {
      stop("The length of k-mers to use in distance measuring must be a ",
           "non-negative integer that is greater than 0 and less than 11.")
    }
    
    if ((!is.numeric(related_MinDist))  || (anyNA(related_MinDist))) {
      stop("The minimum fractional distance required for a related genome to use in ",
           "finding evolutionary conservation must be a valid real number.")
    }
    
    if (length(related_MinDist) != 1) {
      stop("Exactly one number must be inputted as the minimum fractional distance ",
           "required for a related genome to use in finding evolutionary conservation.")
    }
    
    if ((related_MinDist <= 0) || (related_MinDist >= 1)) {
      stop("The minimum fractional distance required for a related genome to use in ",
           "finding evolutionary conservation must be greater than 0 and less than 1.")
    }
    
    if ((!is.numeric(related_MaxDistantN))  || (anyNA(related_MaxDistantN))) {
      stop("The maximum number of most distantly related genomes to use in ",
           "finding evolutionary conservation must be a valid real number.")
    }
    
    if (length(related_MaxDistantN) != 1) {
      stop("Exactly one number must be inputted as the maximum number ",
           "of related genomes to use in finding evolutionary conservation.")
    }
    
    if ((related_MaxDistantN %% 1 != 0) || (related_MaxDistantN <= 0)) {
      stop("The maximum number of related genomes to use in finding evolutionary ",
           "conservation must be a non-negative integer that is greater than 0.")
    }
    
    if ((!is.character(startCodons)) || (anyNA(startCodons))) {
      stop("'startCodons' must consist of valid character strings.")
    }
    
    if (length(startCodons) <= 0) {
      stop("'startCodons' must consist of one or more character strings.")
    }
    
    ## Make all the letters in start codons upper case so it's easier to check.
    startCodons <- toupper(startCodons)
    
    if (any(nchar(startCodons) != 3L) ||
        !(all(unlist(strsplit(startCodons, split = "")) %in% c("A", "C", "G", "T")))) {
      stop("'startCodons' must consist only of three-letter DNA strings.")
    }
  }
  
  ## --------------------------------------------------------------------------------------------------------------- ##
  
  if ((!is.character(strainID)) || (anyNA(strainID))) {
    stop("The strain ID must be a valid character string.")
  }
  
  if (length(strainID) != 1) {
    stop("Exactly one character string must be inputted as the strain ID.")
  }
  
  if ((!is.character(speciesName)) || (anyNA(speciesName))) {
    stop("The species name must be a valid character string.")
  }
  
  if (length(speciesName) != 1) {
    stop("Exactly one character string must be inputted as the species name.")
  }
  
  ## --------------------------------------------------------------------------------------------------------------- ##
  
  fwdGenome <- SearchDB(dbFile = genomes_DBFile, identifier = central_ID,
                        type = "DNAStringSet", verbose = FALSE)
  
  if (length(fwdGenome) <= 0) {
    stop("Central genome not found in database. ",
         "Please check the following inputs: 'genomes_DBFile' and 'central_ID'.")
  }
  
  fwdGenome <- DNAStringSet(unlist(fwdGenome))
  
  ## --------------------------------------------------------------------------------------------------------------- ##
  
  revGenome <- reverseComplement(fwdGenome)
  
  fwdFrame1AA <- suppressWarnings(translate(fwdGenome, if.fuzzy.codon="solve"))
  fwdFrame2AA <- suppressWarnings(translate(subseq(fwdGenome, 2), if.fuzzy.codon="solve"))
  fwdFrame3AA <- suppressWarnings(translate(subseq(fwdGenome, 3), if.fuzzy.codon="solve"))
  
  fwdFrameAA <- list(fwdFrame1AA, fwdFrame2AA, fwdFrame3AA)
  
  fwdProtMap <- lapply(width(fwdGenome), integer)
  fwdProtMap <- list(fwdProtMap, fwdProtMap, fwdProtMap)
  
  revFrame1AA <- suppressWarnings(translate(revGenome, if.fuzzy.codon="solve"))
  revFrame2AA <- suppressWarnings(translate(subseq(revGenome, 2), if.fuzzy.codon="solve"))
  revFrame3AA <- suppressWarnings(translate(subseq(revGenome, 3), if.fuzzy.codon="solve"))
  
  revFrameAA <- list(revFrame1AA, revFrame2AA, revFrame3AA)
  
  revProtMap <- lapply(width(revGenome), integer)
  revProtMap <- list(revProtMap, revProtMap, revProtMap)
  
  ## --------------------------------------------------------------------------------------------------------------- ##
  
  genomeLen <- width(fwdGenome)
  
  fwdConStops <- revConStops <- fwdCov <- revCov <- integer(genomeLen)
  
  fwdConStarts <- revConStarts <- rep.int(NA_integer_, genomeLen)
  
  ## --------------------------------------------------------------------------------------------------------------- ##
  
  # find all stops in each frame
  stopsByFrame <- vector("list", 6)
  
  for (frame in seq_along(fwdFrameAA)) {
    s <- vmatchPattern("*", fwdFrameAA[[frame]])
    s <- unlist(startIndex(s))
    s <- (s - 1)*3 + frame
    stopsByFrame[[frame]] <- s
  }
  
  for (frame in seq_along(revFrameAA)) {
    s <- vmatchPattern("*", revFrameAA[[frame]])
    s <- unlist(startIndex(s))
    s <- (s - 1)*3 + frame
    stopsByFrame[[frame + 3]] <- s
  }
  
  ## --------------------------------------------------------------------------------------------------------------- ##
  
  if (useProt) {
    ## Threshold proteomics hits if necessary.
    
    if (thresholdProt) {
      quantProbs <- seq(0, 1, 0.01)
      scoreQuants <- quantile(protHits_Scores, quantProbs)
      names(scoreQuants) <- as.character(quantProbs)
      
      threshold <- quantProbs[which(quantProbs <= (protHits_Threshold / 100))]
      threshold <- threshold[length(threshold)]
      
      protInds <- which(protHits_Scores > scoreQuants[as.character(threshold)])
      
      protHits_Scores <- protHits_Scores[protInds]
      protHits_Seqs <- protHits_Seqs[protInds]
      
      if (verbose) {
        cat(protHits_Threshold, "percentile threshold applied.\n")
      }
    }
    
    ## --------------------------------------------------------------------------------------------------------------- ##
    
    ## Map proteomics hits.
    
    if (verbose) {
      pBar_Prot <- txtProgressBar(style=ifelse(interactive(), 3, 1))
    }
    
    for (rIdx in seq_along(protHits_Seqs)) {
      peptide <- protHits_Seqs[rIdx]
      
      completeAASeq <- peptide
      startShift <- 0L
      
      ## Methionine (M) is the typical starting amino acid for newly translated proteins,
      ## regardless of what amino acid the first codon of the gene codes for. Usually,
      ## this initial M is cleaved off in post-translation processing but can remain for
      ## a variety of reasons including the enrichment for N-terminal (the starting end
      ## of the protein) peptides in N-terminal proteomcis. Therefore, the first amino
      ## acid of a peptide can potentially be ignored when mapping proteomics hits back
      ## to the central genome.
      
      ## Get the first amino acid of the peptide.
      firstAA <- unlist(strsplit(peptide, ""))[1]
      
      ## Check if it is an M. 
      isStartAA <- (firstAA == "M")
      
      if (isStartAA) {
        peptide <- substring(peptide, 2)
        startShift <- -1L
      }
      
      # Match to the 6-frame translations.
      matchFrameAA <- vector("list", 6)
      
      ## Check the forward frames (frames 1-3).
      for (frame in seq_along(fwdFrameAA)) {
        matchFrameAA[[frame]] <- vmatchPattern(peptide, fwdFrameAA[[frame]])
      }
      
      ## Check the reverse frames (frames 4-6).
      for (frame in seq_along(revFrameAA)) {
        matchFrameAA[[frame + length(fwdFrameAA)]] <- vmatchPattern(peptide, revFrameAA[[frame]])
      }
      
      matchLens <- lapply(matchFrameAA, lengths)
      
      if (sum(unlist(matchLens)) != 1) {
        if (isStartAA) {
          # Match to the 6-frame translations.
          matchFrameAA <- vector("list", 6)
          
          ## Check the forward frames (frames 1-3).
          for (frame in seq_along(fwdFrameAA)) {
            matchFrameAA[[frame]] <- vmatchPattern(completeAASeq, fwdFrameAA[[frame]])
          }
          
          ## Check the reverse frames (frames 4-6).
          for (frame in seq_along(revFrameAA)) {
            matchFrameAA[[frame + length(fwdFrameAA)]] <- vmatchPattern(completeAASeq, revFrameAA[[frame]])
          }
          
          matchLens <- lapply(matchFrameAA, lengths)
          
          if (sum(unlist(matchLens)) != 1) {
            next
          }
          
          startShift <- 0
          
        } else {
          next
        }
      }
      
      cFrame <- which(sapply(matchLens, sum)==1)
      cSeq <- which(matchLens[[cFrame]]==1)
      
      begIdx <- unlist(startIndex(matchFrameAA[[cFrame]][cSeq])) + startShift
      endIdx <- unlist(endIndex(matchFrameAA[[cFrame]][cSeq]))
      
      if (cFrame <= length(fwdFrameAA)) { # forward hit
        begIdx <- (begIdx - 1)*3 + cFrame
        endIdx <- endIdx*3  + cFrame - 1
        
        fwdProtMap[[cFrame]][[cSeq]][begIdx:endIdx] <-
          fwdProtMap[[cFrame]][[cSeq]][begIdx:endIdx] + as.numeric(protHits_Scores[rIdx])
        
      } else { # reverse hit
        cFrame <- cFrame - 3
        begIdx <- (begIdx - 1)*3 + cFrame
        endIdx <- endIdx*3  + cFrame - 1
        
        revProtMap[[cFrame]][[cSeq]][begIdx:endIdx] <-
          revProtMap[[cFrame]][[cSeq]][begIdx:endIdx] + as.numeric(protHits_Scores[rIdx])
      }
      
      if (verbose) {
        setTxtProgressBar(pBar_Prot, rIdx / length(protHits_Seqs))
      }
    }
  } else {
    protHits_IsNTerm <- FALSE
  }
  
  ## --------------------------------------------------------------------------------------------------------------- ##
  
  if (useCons) {
    ## Determine how distant the related genomes are and which ones to use.
    
    central_KMers <- oligonucleotideFrequency(fwdGenome,
                                              width = related_KMerLen,
                                              as.prob = TRUE)
    
    numRelated <- length(related_IDs)
    distToCentral <- rep(NA_real_, numRelated)
    
    if (verbose) {
      pBar_Freq <- txtProgressBar(style=ifelse(interactive(), 3, 1))
    }
    
    for (rIdx in seq_along(related_IDs)) {
      currRGenome <- SearchDB(genomes_DBFile,
                              identifier = related_IDs[rIdx],
                              type = "DNAStringSet",
                              verbose = FALSE)
      
      currRGenome <- DNAStringSet(unlist(currRGenome))
      
      currFreqs <- oligonucleotideFrequency(currRGenome,
                                            width = related_KMerLen,
                                            as.prob = TRUE)
      
      distToCentral[rIdx] <- sum(abs(central_KMers - currFreqs)) / min(sum(central_KMers), sum(currFreqs))
      
      if (verbose) {
        setTxtProgressBar(pBar_Freq, rIdx / numRelated)
      }
    }
    
    ## --------------------------------------------------------------------------------------------------------------- ##
    
    ## Determine which related genomes to use.
    
    related_NoNA_Order <- order(distToCentral,
                                decreasing = TRUE,
                                na.last = NA)
    
    if (length(related_NoNA_Order) <= 0) {
      stop("No related genomes found in database. ",
           "Please check the following inputs: 'genomes_DBFile' and 'related_IDs'.")
    }
    
    validRelated_IDs <- related_IDs[related_NoNA_Order]
    validDistToCentral <- distToCentral[related_NoNA_Order]
    
    dVal <- 1 - ((1 - related_MinDist) ^ related_KMerLen)
    
    distantRelatedIdxs <- which(validDistToCentral >= dVal)
    
    if (length(distantRelatedIdxs) <= 0) {
      stop("Related genomes are too similar. ",
           "Please check the following inputs: 'genomes_DBFile' and 'related_IDs'.")
    }
    
    validRelated_IDs <- validRelated_IDs[distantRelatedIdxs]
    
    if (length(validRelated_IDs) > related_MaxDistantN) {
      validRelated_IDs <- validRelated_IDs[seq_len(related_MaxDistantN)]
    }
    
    numTopR <- length(validRelated_IDs)
    
    ## --------------------------------------------------------------------------------------------------------------- ##
    
    ## Set up conservation and coverage vectors.
    
    for (gPos in seq_len(genomeLen - 2)) {
      fwdCodon <- as.character(subseq(fwdGenome, gPos, gPos + 2))
      
      if (fwdCodon %in% startCodons) {
        fwdConStarts[gPos] <- 0L
      }
      
      revCodon <- as.character(subseq(revGenome, gPos, gPos + 2))
      
      if (revCodon %in% startCodons) {
        revConStarts[gPos] <- 0L
      }
    }
    
    stopCodons <- c("TAG", "TGA", "TAA")
    
    revCompStartCodons <- as.character(reverseComplement(DNAStringSet(startCodons)))
    revCompStopCodons <- as.character(reverseComplement(DNAStringSet(stopCodons)))
    
    ## --------------------------------------------------------------------------------------------------------------- ##
    
    ## Find evolutionarily conserved starts.
    
    if (verbose) {
      pBar_Synteny <- txtProgressBar(style=ifelse(interactive(), 3, 1))
    }
    
    for (rIdx in seq_along(validRelated_IDs)) {
      
      synteny <- FindSynteny(genomes_DBFile,
                             identifier = c(central_ID, validRelated_IDs[rIdx]),
                             verbose = FALSE)
      
      dna <- AlignSynteny(synteny,
                          genomes_DBFile,
                          verbose = FALSE)
      
      dna <- unlist(dna[[1]])
      
      dna <- as.character(dna)
      
      s <- strsplit(dna, "", fixed=TRUE)
      
      for (iIdx in seq_len(length(s)/2)) {
        w <- synteny[2,1][[1]][iIdx, "start1"]:synteny[2,1][[1]][iIdx, "end1"]
        
        fwdCov[w] <- fwdCov[w] + 1L
        revCov[genomeLen - w + 1] <- revCov[genomeLen - w + 1] + 1L
        
        pos <- which(s[[iIdx*2 - 1]] != "-")
        pos1 <- pos[-c(length(pos) - 1, length(pos))]
        pos2 <- pos[-c(1, length(pos))]
        pos3 <- pos[-c(1, 2)]
        
        # which nucleotides are identical
        #identical <- which(s[[j*2 - 1]][pos] == s[[j*2]][pos])
        # which codons are both ATG, GTG, or TTG
        
        centralCodons <- paste(s[[(iIdx * 2) - 1]][pos1],
                               s[[(iIdx * 2) - 1]][pos2],
                               s[[(iIdx * 2) - 1]][pos3],
                               sep = "")
        
        relatedCodons <- paste(s[[iIdx * 2]][pos1],
                               s[[iIdx * 2]][pos2],
                               s[[iIdx * 2]][pos3],
                               sep = "")
        
        ## Forward conserved starts
        identical <- which((centralCodons %in% startCodons) & (relatedCodons %in% startCodons))
        
        if (length(identical) > 0) {
          fwdConStarts[w[identical]] <- fwdConStarts[w[identical]] + 1L
        }
        
        ## Reverse conserved starts
        identical <- which((centralCodons %in% revCompStartCodons) & (relatedCodons %in% revCompStartCodons))
        
        if (length(identical) > 0) {
          revConStarts[genomeLen - w[identical] - 1] <- revConStarts[genomeLen - w[identical] - 1] + 1L
        }
        
        ## Forward conserved stops
        hasStop <- which(relatedCodons %in% stopCodons)
        
        if (length(hasStop) > 0){
          fwdConStops[w[hasStop]] <- fwdConStops[w[hasStop]] + 1L
        }
        
        ## Reverse conserved stops
        hasStop <- which(relatedCodons %in% revCompStopCodons)
        
        if (length(hasStop) > 0){
          revConStops[genomeLen - w[hasStop] - 1] <- revConStops[genomeLen - w[hasStop] - 1] + 1L
        }
      }
      
      if (verbose) {
        setTxtProgressBar(pBar_Synteny, rIdx / numTopR)
      }
    }
  } else {
    validRelated_IDs <- character()
  }
  
  if (verbose) {
    cat("\n")
  }
  
  ## --------------------------------------------------------------------------------------------------------------- ##
  
  ## Save data to an assessment object.
  
  return(structure(list("StrainID" = strainID,
                        "Species" = speciesName,
                        "GenomeLength" = genomeLen,
                        "StopsByFrame" = stopsByFrame,
                        "NTermProteomics" = protHits_IsNTerm,
                        "FwdProtHits" = fwdProtMap,
                        "RevProtHits" = revProtMap,
                        "FwdCoverage" = fwdCov,
                        "FwdConStarts" = fwdConStarts,
                        "FwdConStops" = fwdConStops,
                        "RevCoverage" = revCov,
                        "RevConStarts" = revConStarts,
                        "RevConStops" = revConStops,
                        "NumRelatedGenomes" = length(validRelated_IDs),
                        "HasProteomics" = useProt,
                        "HasConservation" = useCons),
                   class = c("Assessment", "DataMap")))
}