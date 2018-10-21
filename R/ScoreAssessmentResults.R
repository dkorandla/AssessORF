#' @export
#' @importFrom methods is
#'
#' @title Score Gene Assessment Results
#' @description Scores the results from the assessment of a set of genes using one of three modes
#'
#' @param x An object of class \code{Assessment} and subclass \code{Results}.
#'
#' @param mode Must either be "a" (use all evidence), "p" (use proteomics evidence only),
#' "c" (use evolutionary conservation evidence only), or "w" (use all evidence but with weights).
#' 
#' @details
#' \code{ScoreAssessmentResults} calculates an accuracy-like score for the categorization of genes within the given
#' \code{Results} object using the given mode of calculation. The accuracy score for a mode is equal to the number of genes
#' that were categorized to be correct for that mode divided by the total number of genes that could have been categorized
#' as correct for that mode (i.e. a count of the number of genes that had available and useable evidence for that particular
#' mode).
#' 
#' Open reading frames with proteomics evidence but no predicted start are included in the total gene count when calculating
#' the accuracy score for the proteomics mode and for both all evidence modes.
#' 
#' In the weighted, all evidence mode, weights for each category are determined by the number of types of evidence that are
#' supporting or against genes in the category. Counts are multiplied by their corresponding weight, and the maximum value for
#' a weight is 2.
#' 
#' @return A numeric vector of length one containing the calculated accuracy-like score. 
#' 
#' @seealso \code{\link{Assessment-class}}
#' 
#' @examples
#'
#' currResObj <- readRDS(system.file("extdata",
#'                                   "MGAS5005_PreSaved_ResultsObj_Prodigal.rds",
#'                                   package = "AssessORF"))
#'
#' ScoreAssessmentResults(currResObj, "a")
#' 
#' ScoreAssessmentResults(currResObj, "c")
#' 
#' ScoreAssessmentResults(currResObj, "p")
#' 
ScoreAssessmentResults <- function(x, mode = "a") {
  if (!(is(x, "Assessment") && is(x, "Results"))) {
    stop("'x' must be an object of class 'Assessment' and subclass 'Results'.")
  }
  
  if ((!is.character(mode)) || (anyNA(mode)) || (length(mode) != 1) || nchar(mode) != 1) {
    stop("'mode' must be a valid, single character.")
  }
  
  catSumTable <- table(x$CategoryAssignments)
  
  allCatSums <- integer(14L)
  
  names(allCatSums) <- c("Y CS+ PE+", "Y CS+ PE-", "Y CS- PE+", "Y CS- PE-",
                         "Y CS< PE!", "Y CS- PE!", "Y CS! PE+", "Y CS! PE-",
                         "Y CS> PE+", "Y CS> PE-", "Y CS< PE+", "Y CS< PE-",
                         "N CS< PE+", "N CS- PE+")
  
  allCatSums[names(catSumTable)] <- catSumTable
  allCatSums["N CS< PE+"] <- sum(nrow(x$'N_CS<_PE+_ORFs'), na.rm = TRUE)
  allCatSums["N CS- PE+"] <- sum(nrow(x$'N_CS-_PE+_ORFs'), na.rm = TRUE)
  
  if (tolower(mode) == "p") {
    
    accScore <- sum(allCatSums[c("Y CS+ PE+", "Y CS- PE+", "Y CS! PE+", "Y CS> PE+", "Y CS< PE+")]) /
      sum(allCatSums[c("Y CS+ PE+", "Y CS- PE+", "Y CS> PE+", "Y CS< PE+", "Y CS! PE+",
                       "Y CS< PE!", "Y CS- PE!",
                       "N CS< PE+", "N CS- PE+")])
    
  } else if (tolower(mode) == "c") {
    
    accScore <- sum(allCatSums[c("Y CS+ PE+", "Y CS+ PE-")]) /
      sum(allCatSums[c("Y CS+ PE+", "Y CS+ PE-",
                       "Y CS< PE!", "Y CS! PE+", "Y CS! PE-",
                       "Y CS> PE+", "Y CS> PE-", "Y CS< PE+", "Y CS< PE-",
                       "N CS< PE+")])
    
  } else if (tolower(mode) == "a") {
    
    accScore <- sum(allCatSums[c("Y CS+ PE+", "Y CS+ PE-", "Y CS- PE+")]) /
      sum(allCatSums[c("Y CS+ PE+", "Y CS+ PE-", "Y CS- PE+",
                       "Y CS< PE!", "Y CS- PE!", "Y CS! PE+", "Y CS! PE-",
                       "Y CS> PE+", "Y CS> PE-", "Y CS< PE+", "Y CS< PE-",
                       "N CS< PE+", "N CS- PE+")])
    
  } else if (tolower(mode) == "w") {
    
    allCatSums[c("Y CS+ PE+", "Y CS< PE!", "N CS< PE+")] <- allCatSums[c("Y CS+ PE+", "Y CS< PE!", "N CS< PE+")] * 2L
    
    accScore <- sum(allCatSums[c("Y CS+ PE+", "Y CS+ PE-", "Y CS- PE+")]) /
      sum(allCatSums[c("Y CS+ PE+", "Y CS+ PE-", "Y CS- PE+",
                       "Y CS< PE!", "Y CS- PE!", "Y CS! PE+", "Y CS! PE-",
                       "Y CS> PE+", "Y CS> PE-", "Y CS< PE+", "Y CS< PE-",
                       "N CS< PE+", "N CS- PE+")])
    
  } else {
    stop("Invalid mode")
  }
  
  return(accScore)
}