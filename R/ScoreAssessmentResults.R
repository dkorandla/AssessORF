#' @export
#'
#' @title Score Gene Assessment Results
#' @description Scores the results from the assessment of a set of genes using one of three modes
#'
#' @param x An object of class \code{Assessment} and subclass \code{Results}
#'
#' @param mode Must either be "a" (use all evidence), "p" (use proteomics evidence only),
#' or "c" (use evolutionary conservation evidence only)
#' 
#' @seealso \code{\link{Assessment-class}}
#' 
#' @examples
#'
#' ScoreAssessmentResults(readRDS(system.file("extdata", "ATCC17978_PreSaved_ResultsObj_Prodigal.rds", package = "AssessORF")), "a")
#' 
#' ScoreAssessmentResults(readRDS(system.file("extdata", "ATCC17978_PreSaved_ResultsObj_Prodigal.rds", package = "AssessORF")), "c")
#' 
#' ScoreAssessmentResults(readRDS(system.file("extdata", "ATCC17978_PreSaved_ResultsObj_Prodigal.rds", package = "AssessORF")), "p")
#' 
ScoreAssesmentResults <- function(x, mode = "a") {
  if (class(x)[1] != "Assessment") {
    stop("x must be an object of class 'Assessment'.")
  }
  
  if (class(x)[2] != "Results") {
    stop("x must be of subclass 'Results'.")
  }
  
  if ((!is.character(mode)) || (anyNA(mode)) || (length(mode) != 1) || nchar(mode) != 1) {
    stop("mode must be a valid, single character.")
  }
  
  catSumTable <- table(x$CategoryAssignments)
  
  allCatSums <- integer(13L)
  
  names(allCatSums) <- c("Y CS+ PE+", "Y CS+ PE-", "Y CS- PE+", "Y CS- PE-",
                         "Y CS! PE-", "Y CS< PE!", "Y CS- PE!",
                         "Y CS> PE+", "Y CS> PE-", "Y CS< PE+", "Y CS< PE-",
                         "N CS< PE+", "N CS- PE+")
  
  allCatSums[names(catSumTable)] <- catSumTable
  allCatSums["N CS< PE+"] <- sum(nrow(x$'N_CS<_PE+_ORFs'), na.rm = TRUE)
  allCatSums["N CS- PE+"] <- sum(nrow(x$'N_CS-_PE+_ORFs'), na.rm = TRUE)
  
  if (tolower(mode) == "p") {
    
    accScore <- sum(allCatSums[c("Y CS+ PE+", "Y CS- PE+", "Y CS> PE+", "Y CS< PE+")]) /
      sum(allCatSums[c("Y CS+ PE+", "Y CS- PE+", "Y CS> PE+", "Y CS< PE+",
                       "Y CS< PE!", "Y CS- PE!", "N CS< PE+", "N CS- PE+")])
    
  } else if (tolower(mode) == "c") {
    
    accScore <- sum(allCatSums[c("Y CS+ PE+", "Y CS+ PE-")]) /
      sum(allCatSums[c("Y CS+ PE+", "Y CS+ PE-",
                       "Y CS! PE-", "Y CS< PE!",
                       "Y CS> PE+", "Y CS> PE-", "Y CS< PE+", "Y CS< PE-",
                       "N CS< PE+")])
    
  } else if (tolower(mode) == "a") {
    
    accScore <- sum(allCatSums[c("Y CS+ PE+", "Y CS+ PE-", "Y CS- PE+")]) /
      sum(allCatSums[c("Y CS+ PE+", "Y CS+ PE-", "Y CS- PE+",
                       "Y CS! PE-", "Y CS< PE!", "Y CS- PE!",
                       "Y CS> PE+", "Y CS> PE-", "Y CS< PE+", "Y CS< PE-",
                       "N CS< PE+", "N CS- PE+")])
    
  } else {
    stop("Invalid mode")
  }
  
  return(accScore)
}
