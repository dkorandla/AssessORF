#' @export
#' @importFrom methods is
#' @importFrom utils installed.packages
#' 
#' @title Compare Assessment Results
#' @description Compare two objects of class Assessment. subclass Results to determine how their gene sets and the
#' corresponding category assignments vary
#'
#' @usage
#' CompareAssessmentResults(obj1,
#'                          obj2,
#'                          printSummary = TRUE,
#'                          returnDetails = FALSE)
#'
#' @param obj1,obj2 Objects of class \code{Assessment} and subclass \code{Results} to compare against each other.
#' Alternatively, either \code{obj1} or \code{obj2} (or both) can be a two-element character vector that specifies one
#' of such objects from \code{AssessORFData}. The first element in the vector should be the strain identifier, and the
#' second element should be the gene source identifier. Both objects should have been generated from the same mapping
#' object.
#' 
#' @param printSummary Logical indicating whether or not to print out a summary of the comparison analysis.
#' 
#' @param returnDetails Logical indicating whether or not to return a list of details from the comparison analysis.
#' See the next section for what items are in the list.
#'
#' @details
#' Since the same mapping object (an object class \code{Assessment} and subclass \code{DataMap}) can be used to assess
#' multiple sets of genes for genome, it is meaningful to compare how those gene sets and their category assignments
#' from \code{AssessGenes} vary from one another. To make describing this function easier, let us assume that one set
#' of genes consists of the complete set of predictions made by a gene-finding program on a particular strain's genome
#' and that the other set of genes consists of the complete set of predictions from a second gene-finding program.
#' 
#' When gene-finding programs predict genes for a genome, they make a decision on which regions of the genome code for
#' proteins. There is (usually) only one option for the stop codon that ends a particular coding region, but there are
#' typically multiple options available for the start codon that will mark the beginning of a region. It is therefore
#' useful to find out which (general) coding regions the two programs agree on by determining which stops are found in
#' both sets of predicted genes. From there, the starts each program picked for those shared coding regions can be
#' compared to see whether they agree or not. If the same start is chosen by both programs for a particular shared
#' stop / coding region, then that is an example of a gene predicted by both programs. If the starts chosen by the two
#' programs for a particular shared stop / coding region are different, then that it is an example of both programs
#' agreeing that that particular region of the genome codes for protein but disagreeing on where in the genome that
#' region starts. It would be interesting to see what category was assigned to each start, especially if one start has
#' evolutionary conservation and the other does not.
#' 
#' This function compares the set of genes and corresponding category assignments in the object specified by
#' \code{obj1} (object 1) to the set of genes and corresponding category assignments in the object specified by
#' \code{obj2} (object 2). It then reports the results of the comparison analysis in the format specified by the
#' logical parameters \code{printSummary} and \code{returnDetails}.
#' 
#' If \code{printSummary} is true, the function prints out the following information: the number of shared coding
#' regions (i.e., the number of stops in both gene sets), the number of shared genes (i.e., the number of times both a
#' start and its corresponding stop are found in both sets), and the number of instances where a stop is found in both
#' gene sets but the corresponding starts in each set disagree. For the shared stop - different start set, the function
#' also prints the number of instances where the start from one object has conservation evidence while the
#' corresponding start in the other object does not.
#' 
#' If \code{returnDetails} is true, the function returns a 11-item list. Each item of the list is described below. The
#' contents of the object 1 and object 2 gene vectors correspond to the ordering of the genes inside object 1 or object
#' 2, respectively. For the category assignment matrix for shared stop - different start set, it is possible for the
#' gene in object 1 to be assigned to the same category as the corresponding gene from object 2, and the table reflects
#' that.
#' \itemize{
#'  \item{"StrainID"}{Same as the strain identifier inside \code{obj1} and \code{obj2}}
#'  \item{Species}{Same as the strain identifier inside objects 1 and 2}
#'  \item{Obj1_GeneSource}{Same as the gene source identifier inside \code{obj1}}
#'  \item{Obj2_GeneSource}{Same as the gene source identifier inside \code{obj2}}
#'  \item{Obj1_Genes_SharedCodingRegions}{The genes from object 1 that share a stop with a gene from the object 2}
#'  \item{Obj2_Genes_SharedCodingRegions}{The genes from object 2 that share a stop with a gene from the object 1}
#'  \item{Obj1_Genes_SharedGenes}{The genes from object 1 that share a start and stop with a gene from the object 2}
#'  \item{Obj2_Genes_SharedGeness}{The genes from object 2 that share a start and stop with a gene from the object 1}
#'  \item{Obj1_Genes_SharedStopDiffStart}{The genes from object 1 that share a stop with a gene from the object 2 but
#'  have a different start from the corresponding object 2 gene}
#'  \item{Obj2_Genes_SharedStopDiffStart}{The genes from object 2 that share a stop with a gene from the object 1 but
#'  have a different start from the corresponding object 1 gene}
#'  \item{CategoryTable_SharedStopDiffStart}{A 12-by-12 matrix describing the number of times the gene from object 1
#'  was assigned one category and the gene from object 2 was assigned some other category for the
#'  shared stop - different start set}
#' }
#' 
#' Please ensure that \code{obj1} and \code{obj2} come from the same strain / mapping object. The function will do its
#' best to make sure the identifying information for \code{obj1} and \code{obj2} match.
#' 
#' \code{printSummary} and \code{returnDetails} cannot both be FALSE.
#'
#' @return If \code{returnDetails} is true, the function returns a 11-item list. Otherwise, the function invisibly
#' returns object 1.
#' 
#' @seealso \code{\link{Assessment-class}}, \code{\link{AssessGenes}}
#'
#' @examples
#'
#' ## Example showing how to use the function with the AssessORFData package:
#'
#' \dontrun{
#' compare1 <- CompareAssessmentResults(obj1 = c("MGAS5005", "Prodigal"),
#'                                      obj2 = c("MGAS5005", "GeneMarkS2"),
#'                                      printSummary = TRUE,
#'                                      returnDetails = TRUE)
#' }
#' 
#' resObj1 <- readRDS(system.file("extdata",
#'                                "MGAS5005_PreSaved_ResultsObj_Prodigal.rds",
#'                                package = "AssessORF"))
#'                                
#' resObj2 <- readRDS(system.file("extdata",
#'                                "MGAS5005_PreSaved_ResultsObj_GeneMarkS2.rds",
#'                                package = "AssessORF"))
#'                                
#' compare2 <- CompareAssessmentResults(obj1 = resObj1,
#'                                      obj2 = resObj2,
#'                                      printSummary = TRUE,
#'                                      returnDetails = TRUE)
#' 
CompareAssessmentResults <- function(obj1,
                                     obj2,
                                     printSummary = TRUE,
                                     returnDetails = FALSE) {
  
  ## Check the logical inputs for errors.
  
  if ((!is.logical(printSummary)) || (anyNA(printSummary)) || (length(printSummary) != 1)) {
    stop("'printSummary' must be of type logical, be either TRUE or FALSE, and consist of only 1 element.")
  }
  
  if ((!is.logical(returnDetails)) || (anyNA(returnDetails)) || (length(returnDetails) != 1)) {
    stop("'returnDetails' must be of type logical, be either TRUE or FALSE, and consist of only 1 element.")
  }
  
  if ((!printSummary) && (!returnDetails)) {
    stop("Either 'printSummary' or 'returnDetails' must be TRUE.")
  }
  
  ## --------------------------------------------------------------------------------------------------------------- ##
  
  ## Check the two objects that are to be compared for errors.
  ## If either one is a two-element character vector, get the corresponding results object from AssessORFData.
  
  if (anyNA(obj1)) {
    stop("The first object cannot have any NAs.")
  }
  
  if ((is(obj1, "Assessment")) && (is(obj1, "Results"))) {
    ## Do nothing
  } else if ((is.character(obj1)) && (length(obj1) == 2)) {
    if (!(any(installed.packages()[, 1] == "AssessORFData"))) {
      stop("The package 'AssessORFData' has not been installed.")
    }
    
    obj1_Char <- obj1
    
    obj1 <-AssessORFData::GetResultsObj(obj1[1], obj1[2])
    
  } else {
    stop("The first object must either be an object of class 'Assessment', subclass 'Results' ",
         "or be a character vector of length 2.")
  }
  
  if (anyNA(obj2)) {
    stop("The second object cannot have any NAs.")
  }
  
  if ((is(obj2, "Assessment")) && (is(obj2, "Results"))) {
    ## Do nothing
  } else if ((is.character(obj2)) && (length(obj2) == 2)) {
    if (!(any(installed.packages()[, 1] == "AssessORFData"))) {
      stop("The package 'AssessORFData' has not been installed.")
    }
    
    obj2_Char <- obj2
    
    obj2 <-AssessORFData::GetResultsObj(obj2[1], obj2[2])
    
  } else {
    stop("The second object must either be an object of class 'Assessment', subclass 'Results' ",
         "or be a character vector of length 2.")
  }
  
  ## --------------------------------------------------------------------------------------------------------------- ##
  
  ## Make sure that the two objects that are to be compared come from the same mapping object.
  
  if ((obj1$StrainID != obj2$StrainID) || (obj1$Species != obj2$Species) || (obj1$GenomeLength != obj1$GenomeLength)) {
    stop("The two results objects to be compared against each other do not seem to come from the same mapping object.",
         " Please check your inputs and try again.")
  }
  
  ## --------------------------------------------------------------------------------------------------------------- ##
  
  ## For each results object, group the genes and their category assignments by strand (forward followed by reverse).
  ## Next, condense the positional information into one character vector for gene starts and one charcter vector for
  ## gene stops so that is easier to compare the two sets of genes against each other (format is "position_strand").
  
  obj1_Starts <- c(paste((obj1$GeneLeftPos)[(obj1$GeneStrand == "+")], "+", sep = "_"),
                   paste((obj1$GeneRightPos)[(obj1$GeneStrand == "-")], "-", sep = "_"))
  
  obj1_Stops <-  c(paste((obj1$GeneRightPos)[(obj1$GeneStrand == "+")], "+", sep = "_"),
                   paste((obj1$GeneLeftPos)[(obj1$GeneStrand == "-")], "-", sep = "_"))
  
  obj1_Cats <- c((obj1$CategoryAssignments)[(obj1$GeneStrand == "+")],
                 (obj1$CategoryAssignments)[(obj1$GeneStrand == "-")])
  
  obj2_Starts <- c(paste((obj2$GeneLeftPos)[(obj2$GeneStrand == "+")], "+", sep = "_"),
                   paste((obj2$GeneRightPos)[(obj2$GeneStrand == "-")], "-", sep = "_"))
  
  obj2_Stops <-  c(paste((obj2$GeneRightPos)[(obj2$GeneStrand == "+")], "+", sep = "_"),
                   paste((obj2$GeneLeftPos)[(obj2$GeneStrand == "-")], "-", sep = "_"))
  
  obj2_Cats <- c((obj2$CategoryAssignments)[(obj2$GeneStrand == "+")],
                 (obj2$CategoryAssignments)[(obj2$GeneStrand == "-")])
  
  ## Store the grouping order for use in generating the output.
  obj1_Indices <- c(which(obj1$GeneStrand == "+"), which(obj1$GeneStrand == "-"))
  obj2_Indices <- c(which(obj2$GeneStrand == "+"), which(obj2$GeneStrand == "-"))
  
  ## Generate an identifying string for each object.
  obj1_String <- ifelse(obj1$GeneSource != "", paste("Object 1", obj1$GeneSource, sep = " - "), "Object 1")
  obj2_String <- ifelse(obj2$GeneSource != "", paste("Object 2", obj2$GeneSource, sep = " - "), "Object 2")
  
  ## --------------------------------------------------------------------------------------------------------------- ##
  
  ## For the set of genes where both objects have the same stop but different starts, set up a matrix to store the
  ## number of times the gene from object 1 was assigned one particular category and the gene from object 2 was
  ## assigned another category (possibly the same one).
  
  catNames <- c("Y CS+ PE+", "Y CS+ PE-", "Y CS- PE+", "Y CS- PE-",
                "Y CS< PE!", "Y CS- PE!", "Y CS! PE+", "Y CS! PE-",
                "Y CS> PE+", "Y CS> PE-", "Y CS< PE+", "Y CS< PE-")
  
  catByCatTable <- matrix(0L, nrow = 12L, ncol = 12L,
                          dimnames = list(paste(obj1_String, catNames),
                                          paste(obj2_String, catNames)))
  
  ## Determine the indices of the genes from the object 1 that have the same stop but different start as genes
  ## from object 2. Note that this set of indices is specific to how the genes are grouped inside the function
  ## and not how they are stored inside either object.
  idxs_ShareStopDiffStart <- which((obj1_Stops %in% obj2_Stops) & !(obj1_Starts %in% obj2_Starts))
  
  ## Loop through the same stop, different start genes and increment the appropriate cell in the matrix based on
  ## what category was assigned to the gene in object 1 and what category was assigned to the gene in object 2.
  for (gIdx_1 in idxs_ShareStopDiffStart) {
    currStop_1 <- obj1_Stops[gIdx_1]
    
    gIdx_2 <- which(obj2_Stops == currStop_1)
    
    currCatIdx_1 <- which(catNames == obj1_Cats[gIdx_1])
    currCatIdx_2 <- which(catNames == obj2_Cats[gIdx_2])
    
    catByCatTable[currCatIdx_1, currCatIdx_2] <- catByCatTable[currCatIdx_1, currCatIdx_2] + 1L
  }
  
  ## --------------------------------------------------------------------------------------------------------------- ##
  
  ## This vector stores the indices of the genes in object 1 that have the same stop as a gene in object 2.
  obj1_ShareStop <- obj1_Indices[(obj1_Stops %in% obj2_Stops)]
  
  ## This vector stores the indices of the genes in object 1 that share the same start and stop as a gene in object 2.
  obj1_ShareStartStop <- obj1_Indices[(obj1_Stops %in% obj2_Stops) & (obj1_Starts %in% obj2_Starts)]
  
  ## This vector stores the indices of the genes in object 1 that have the same stop but different start as a gene
  ## in object 2. Note that this set of indices is specific to how the genes are grouped inside object 1.
  obj1_ShareStopDiffStart <- obj1_Indices[idxs_ShareStopDiffStart]
  
  ## --------------------------------------------------------------------------------------------------------------- ##
  
  if (printSummary) {
    ## Put together and print the summary if requested. 
    speciesName <- obj1$Species
    strainID <- obj1$StrainID
    
    if ((speciesName != "") && (strainID != "")) {
      genomeID <- paste(speciesName, strainID)
    } else if ((speciesName != "")) {
      genomeID <- speciesName
    } else if ((strainID != "")) {
      genomeID <- strainID
    } else {
      genomeID <- ""
    }
    
    printOut <- "Comparison of Two Assessment Results Objects"
    
    if (genomeID == "") {
      printOut <- paste0(printOut, "\n")
    } else {
      printOut <- paste0(printOut, " from ", genomeID, "\n")
    }
    
    printOut <- paste0(printOut, obj1_String, " vs ", obj2_String, "\n\n")
    
    printOut <- paste0(printOut,
                       "Number of Shared Coding Regions (number of stops found in both sets): ",
                       length(obj1_ShareStop), "\n\n")
    
    printOut <- paste0(printOut,
                       "Number of Shared Genes (number of start-stop pairs found in both sets): ",
                       length(obj1_ShareStartStop), "\n\n")
    
    printOut <- paste0(printOut, "Number of Shared Stop - Different Start Coding Regions: ",
                       length(obj1_ShareStopDiffStart),
                       "\n(This is where the stop for a coding region is found in both sets,\n",
                       "but the corresponding starts in each set differ.)\n\n")
    
    printOut <- paste0(printOut, "For the shared stop - different start set, there were ",
                       sum(rowSums(catByCatTable[c(1, 2), c(3:6, 9:12)])), " instances\nwhere the start from ",
                       "object 1 has conservation evidence ('CS+')\n",
                       "and the start from object 2 does not ('CS-', 'CS>', 'CS<').\n")
    
    printOut <- paste0(printOut, "For the shared stop - different start set, there were ",
                       sum(colSums(catByCatTable[c(3:6, 9:12), c(1, 2)])), " instances\nwhere the start from ",
                       "object 2 has conservation evidence ('CS+')\n",
                       "and the start from object 1 does not ('CS-', 'CS>', 'CS<').\n")
    
    cat(printOut)
    
  }
  
  if (returnDetails) {
    ## Put together and return the list if requested.
    
    ## These object 2 vectors are similar to the object 1 vectors described in the previous section.
    obj2_ShareStop <- obj2_Indices[(obj2_Stops %in% obj1_Stops)]
    obj2_ShareStartStop <- obj2_Indices[(obj2_Stops %in% obj1_Stops) & (obj2_Starts %in% obj1_Starts)]
    obj2_ShareStopDiffStart <- obj2_Indices[(obj2_Stops %in% obj1_Stops) & !(obj2_Starts %in% obj1_Starts)]
    
    return(list("StrainID" = obj1$StrainID,
                "Species" = obj1$Species,
                "Obj1_GeneSource" = obj1$GeneSource,
                "Obj2_GeneSource" = obj2$GeneSource,
                "Obj1_Genes_SharedCodingRegions" = obj1_ShareStop,
                "Obj2_Genes_SharedCodingRegions" = obj2_ShareStop,
                "Obj1_Genes_SharedGenes" = obj1_ShareStartStop,
                "Obj2_Genes_SharedGeness" = obj2_ShareStartStop,
                "Obj1_Genes_SharedStopDiffStart" = obj1_ShareStopDiffStart,
                "Obj2_Genes_SharedStopDiffStart" = obj2_ShareStopDiffStart,
                "CategoryTable_SharedStopDiffStart" = catByCatTable))
  } else {
    invisible(obj1)
  }
}