#' @export
#' 
#' @method as.matrix Assessment
#'
#' @title Tabulate the Category Assignments for Assessment Results Objects
#' @description The \code{as.matrix} method for \code{Assessment} and subclass \code{Results} objects
#'
#' @param x An object of class \code{Assessment} and subclass \code{Results}.
#'
#' @param ... Additional arguments.
#'
#' @details
#' \code{as.matrix.Assessment} tabulates and returns the number of times each category appears in the \code{CategoryAssignments}
#' vector within the given \code{Results} object. If the number of genes for any the 12 main gene / ORF categories is zero, a
#' count (of zero) will still be included for that category.
#'
#' @return A one-row matrix with the counts for the number of genes/ORFs that fall into each category. The corresponding
#' category codes serve as the column names, and the name of the row is the strain ID.
#'
#' @seealso \code{\link{Assessment-class}}
#' 
#' @examples
#'
#' as.matrix(readRDS(system.file("extdata",
#'                               "MGAS5005_PreSaved_ResultsObj_Prodigal.rds",
#'                               package = "AssessORF")))
#'
as.matrix.Assessment <- function(x, ...) {
  if (class(x)[1] != "Assessment") {
    stop("'x' must be an object of class 'Assessment'.")
  }
  
  if (class(x)[2] == "Results") {
    speciesName <- x$Species
    strainID <- x$StrainID
    
    if ((speciesName != "") && (strainID != "")) {
      genomeID <- paste(speciesName, strainID, sep = "_")
    } else if ((speciesName != "")) {
      genomeID <- speciesName
    } else if ((strainID != "")) {
      genomeID <- strainID
    } else {
      genomeID <- ""
    }
    
    geneSource <- x$GeneSource
    
    if ((genomeID != "") && (geneSource != "")) {
      rowTitle <- paste(genomeID, geneSource)
    } else if (genomeID != "") {
      rowTitle <- paste(genomeID, "UnknownSource")
    } else if (geneSource != "") {
      rowTitle <- paste("UnknownGenome", geneSource)
    } else {
      rowTitle <- "UnknownGenome UnknownSource"
    }
    
    ## --------------------------------------------------------------------------------------------------------------- ##
    
    ## Get the number of genes in each category.
    catSumTable <- table(x$CategoryAssignments)
    
    allCatSums <- integer(14L)
    
    names(allCatSums) <- c("Y CS+ PE+", "Y CS+ PE-", "Y CS- PE+", "Y CS- PE-",
                           "Y CS< PE!", "Y CS- PE!", "Y CS! PE+", "Y CS! PE-",
                           "Y CS> PE+", "Y CS> PE-", "Y CS< PE+", "Y CS< PE-",
                           "N CS< PE+", "N CS- PE+")
    
    allCatSums[names(catSumTable)] <- catSumTable
    
    allCatSums["N CS< PE+"] <- sum(nrow(x$'N_CS<_PE+_ORFs'), na.rm = TRUE)
    allCatSums["N CS- PE+"] <- sum(nrow(x$'N_CS-_PE+_ORFs'), na.rm = TRUE)
    
    return(matrix(allCatSums, nrow = 1, ncol = length(allCatSums),
                  dimnames = list(rowTitle, names(allCatSums))))
  } else {
    stop("x must be of subclass 'Results'.")
  }
}

#' @export
#'
#' @title Print Assessment Objects
#' @description The \code{print} method for \code{Assessment} objects
#'
#' @param x An object of class \code{Assessment} and of either subclass \code{DataMap} or subclass \code{Results}.
#'
#' @param ... Further printing parameters.
#'
#' @details
#' If \code{x} is of subclass \code{DataMap}, the length of the genome is printed along with any supplied identifying
#' information for the genome.
#'
#' If \code{x} is of subclass \code{Results}, the number of genes in each category and the accuracy scores are printed out
#' along with any supplied identifying information.
#' 
#' @return Invisibly returns the input object \code{x}
#'
#' @seealso \code{\link{Assessment-class}}
#' 
#' @examples
#'
#' print(readRDS(system.file("extdata",
#'                           "MGAS5005_PreSaved_DataMapObj.rds",
#'                           package = "AssessORF")))
#' 
#' print(readRDS(system.file("extdata",
#'                           "MGAS5005_PreSaved_ResultsObj_Prodigal.rds",
#'                           package = "AssessORF")))
#'
print.Assessment <- function(x, ...) {
  if (class(x)[1] != "Assessment") {
    stop("'x' must be an object of class 'Assessment'.")
  }
  
  if (class(x)[2] == "DataMap") {
    speciesName <- x$Species
    strainID <- x$StrainID
    
    if ((speciesName != "") && (strainID != "")) {
      genomeID <- paste(speciesName, strainID)
    } else if ((speciesName != "")) {
      genomeID <- speciesName
    } else if ((strainID != "")) {
      genomeID <- strainID
    } else {
      genomeID <- ""
    }
    
    printOut <- "An Assessment object that maps proteomics data and evolutionary data\n"
    
    if (genomeID != "") {
      printOut <- paste(printOut, genomeID, "\n", sep = "")
    }
    
    cat(paste(printOut, "Genome Length: ", x$GenomeLength, "\n", sep = ""))
    
  } else if (class(x)[2] == "Results") {
    
    speciesName <- x$Species
    strainID <- x$StrainID
    
    if ((speciesName != "") && (strainID != "")) {
      genomeID <- paste(speciesName, strainID)
    } else if ((speciesName != "")) {
      genomeID <- speciesName
    } else if ((strainID != "")) {
      genomeID <- strainID
    } else {
      genomeID <- ""
    }
    
    printOut <- "An Assessment object that categorizes predicted genes\n"
    
    if (genomeID != "") {
      printOut <- paste(printOut, "Strain: ", genomeID, "\n", sep = "")
    }
    
    if (x$GeneSource != "") {
      printOut <- paste(printOut, "Number of Genes Provided from ", x$GeneSource, ": ", x$NumGenes, "\n\n", sep = "")
    } else {
      printOut <- paste(printOut, "Number of Genes Provided: ", x$NumGenes, "\n\n", sep = "")
    }
    
    ## --------------------------------------------------------------------------------------------------------------- ##
    
    ## Get the number of genes in each category.
    catSumTable <- table(x$CategoryAssignments)
    
    allCatSums <- integer(14L)
    
    names(allCatSums) <- c("Y CS+ PE+", "Y CS+ PE-", "Y CS- PE+", "Y CS- PE-",
                           "Y CS< PE!", "Y CS- PE!", "Y CS! PE+", "Y CS! PE-",
                           "Y CS> PE+", "Y CS> PE-", "Y CS< PE+", "Y CS< PE-",
                           "N CS< PE+", "N CS- PE+")
    
    allCatSums[names(catSumTable)] <- catSumTable
    allCatSums["N CS< PE+"] <- sum(nrow(x$'N_CS<_PE+_ORFs'), na.rm = TRUE)
    allCatSums["N CS- PE+"] <- sum(nrow(x$'N_CS-_PE+_ORFs'), na.rm = TRUE)
    
    ## --------------------------------------------------------------------------------------------------------------- ##
    
    printOut <- paste(printOut, "Score Using All Evidence: ",
                      round(ScoreAssessmentResults(x), 4), "\n", sep = "")
    
    printOut <- paste(printOut, "Score Using Only Proteomic Evidence: ",
                      round(ScoreAssessmentResults(x, "p"), 4), "\n", sep = "")
    
    printOut <- paste(printOut, "Score Using Only Conserved Start Evidence: ",
                      round(ScoreAssessmentResults(x, "c"), 4), "\n", sep = "")
    
    printOut <- paste(printOut, "Score Using All Evidence With Weights: ",
                      round(ScoreAssessmentResults(x, "w"), 4), "\n\n", sep = "")
    
    ## --------------------------------------------------------------------------------------------------------------- ##
    
    catIndex <- names(allCatSums) == "Y CS+ PE+"
    printOut <- paste0(printOut, names(allCatSums[catIndex]),
                       " (correct with strong evidence): ", allCatSums[catIndex], "\n")
    
    catIndex <- (names(allCatSums) == "Y CS+ PE-") | (names(allCatSums) == "Y CS- PE+")
    printOut <- paste0(printOut, paste0(names(allCatSums[catIndex]),
                                        " (correct with some evidence): ",
                                        allCatSums[catIndex], "\n" , collapse = ""))
    
    catIndex <- names(allCatSums) == "Y CS- PE-"
    printOut <- paste0(printOut, names(allCatSums[catIndex]),
                       " (no evidence): ", allCatSums[catIndex], "\n" )
    
    catIndex <- names(allCatSums) == "Y CS< PE!"
    printOut <- paste0(printOut, names(allCatSums[catIndex]),
                       " (definitely incorrect): ", allCatSums[catIndex], "\n")
    
    catIndex <- (names(allCatSums) == "Y CS- PE!") | (names(allCatSums) == "Y CS! PE+") |
      (names(allCatSums) == "Y CS! PE-")
    printOut <- paste0(printOut, paste0(names(allCatSums[catIndex]),
                                        " (likely incorrect): ",
                                        allCatSums[catIndex], "\n" , collapse = ""))
    
    catIndex <- (names(allCatSums) == "Y CS> PE+") | (names(allCatSums) == "Y CS> PE-") |
      (names(allCatSums) == "Y CS< PE+") | (names(allCatSums) == "Y CS< PE-")
    printOut <- paste0(printOut, paste0(names(allCatSums[catIndex]),
                                        " (potentially incorrect): ",
                                        allCatSums[catIndex], "\n" , collapse = ""))
    
    catIndex <- names(allCatSums) == "N CS< PE+"
    printOut <- paste0(printOut, names(allCatSums[catIndex]),
                       " (likely missing genes): ", allCatSums[catIndex], "\n")
    
    catIndex <- names(allCatSums) == "N CS- PE+"
    printOut <- paste0(printOut, names(allCatSums[catIndex]),
                       " (potentially missing genes): ", allCatSums[catIndex], "\n")
    
    ## --------------------------------------------------------------------------------------------------------------- ##
    
    printOut <- paste(printOut, "\nNumber of Genes with Supporting Evidence: ",
                      sum(allCatSums[c("Y CS+ PE+", "Y CS+ PE-", "Y CS- PE+")]), sep = "")
    
    printOut <- paste(printOut, "\nNumber of Genes with Contradictory Evidence: ",
                      sum(allCatSums[c("Y CS> PE+", "Y CS> PE-", "Y CS< PE+", "Y CS< PE-",
                                       "Y CS! PE+", "Y CS! PE-", "Y CS< PE!", "Y CS- PE!")]),
                      sep = "")
    
    printOut <- paste(printOut, "\nNumber of ORFs with Protein Evidence and No Given Start: ",
                      sum(allCatSums[c("N CS< PE+", "N CS- PE+")]), "\n", sep = "")
    
    cat(printOut)
    
  } else {
    stop("'x' is an object of unrecognized subclass '", class(x)[2], "'.")
  }
  
  invisible(x)
}

#' @export
#' @import graphics
#' @import grDevices
#'
#' @title Plot Assessment Objects
#' @description The \code{plot} method for \code{Assessment} objects
#'
#' @param x An object of class \code{Assessment} and of either subclass \code{DataMap} or subclass \code{Results}.
#'
#' @param y An optional object of class \code{Assessment} and of either subclass \code{DataMap} or subclass \code{Results}. Its
#' subclass must be different than the subclass of \code{x}
#'
#' @param related_MinConStart Minimum value of the conservation to coverage ratio needed to call a start conserved. Must range
#' from 0 to 1. Lower values allow more conserved starts through. Recommended to use default value.
#'
#' @param ... Further plotting parameters.
#'
#' @details
#' If only \code{x} is specified and \code{x} is of subclass \code{DataMap}, an interactive genome viewer showing how the
#' proteomics data and evolutionary conservation data maps to the central genome is plotted.
#'
#' If only \code{x} is specified and \code{x} is of subclass \code{Results}, a bar chart describing the number of genes in each
#' category is plotted.
#'
#' If both \code{x} and \code{y} are specified, an interactive genome viewer showing how the proteomics data, evolutionary
#' evolutionary conservation data, and gene set map to the central genome is plotted.
#' 
#' @return Invisibly returns the input object \code{x}
#'
#' @seealso \code{\link{Assessment-class}}, \code{\link{locator}}
#' 
#' @examples
#'
#' currMapObj <- readRDS(system.file("extdata",
#'                                   "MGAS5005_PreSaved_DataMapObj.rds",
#'                                   package = "AssessORF"))
#' 
#' currResObj <- readRDS(system.file("extdata",
#'                                   "MGAS5005_PreSaved_ResultsObj_Prodigal.rds",
#'                                   package = "AssessORF"))
#'
#' plot(currMapObj)
#' 
#' plot(currResObj)
#' 
#' plot(currMapObj, currResObj)
#' 
#' plot(currResObj, currMapObj)
#'
plot.Assessment <- function(x, y = NULL,
                            related_MinConStart = 0.8, ...) {
  
  if (class(x)[1] != "Assessment") {
    stop("'x' must be an object of class 'Assessment'.")
  }
  
  if (!is.null(y)) {
    if (class(y)[1] != "Assessment") {
      stop("'y' must be an object of class 'Assessment'.")
    }
    
    if (!is.numeric(related_MinConStart)) {
      stop("The value for the conserved start threshold must be a real number.")
    }
    
    if (class(x)[2] == "DataMap") {
      mapObj <- x
      
      if (class(y)[2] != "Results") {
        stop("'y' must be an object of subclass 'Results' when x is an object of subclass 'DataMap'.")
      }
      
      predObj <- y
    } else if (class(x)[2] == "Results") {
      predObj <- x
      
      if (class(y)[2] != "DataMap") {
        stop("'y' must be an object of subclass 'DataMap' when x is an object of subclass 'Results'.")
      }
      
      mapObj <- y
    } else {
      stop("'x' is an unrecognized object of class '", class(x)[2], "'.")
    }
    
    PlotAssessmentMapping(mapObj, predObj, related_MinConStart)
    
  } else if (class(x)[2] == "DataMap") {
    if (!is.numeric(related_MinConStart)) {
      stop("The value for the conserved start threshold must be a real number.")
    }
    
    PlotAssessmentMapping(x, NULL, related_MinConStart)
    
  } else if (class(x)[2] == "Results") {
    geneSource <- x$GeneSource
    
    part2 <- "Gene Category Assignments"
    
    if (geneSource != "") {
      part2 <- paste(geneSource, part2)
    }
    
    ## --------------------------------------------------------------------------------------------------------------- ##
    
    speciesName <- x$Species
    strainID <- x$StrainID
    
    if ((speciesName != "") && (strainID != "")) {
      plotTitle <- bquote(italic(.(speciesName))~.(strainID)~.(part2))
    } else if ((speciesName != "")) {
      plotTitle <- bquote(italic(.(speciesName))~.(part2))
    } else if ((strainID != "")) {
      plotTitle <- bquote(~.(strainID)~.(part2))
    } else {
      plotTitle <- part2
    }
    
    ## --------------------------------------------------------------------------------------------------------------- ##
    
    ## Get the number of genes in each category
    catSumTable <- table(x$CategoryAssignments)
    
    nonCollapseIdxs <- which((names(catSumTable) != "Y CS! PE+") &
                               (names(catSumTable) != "Y CS! PE-") &
                               (names(catSumTable) != "Y CS> PE+") &
                               (names(catSumTable) != "Y CS> PE-") &
                               (names(catSumTable) != "Y CS< PE+") &
                               (names(catSumTable) != "Y CS< PE-"))
    
    allCatSums <- integer(11L)
    
    names(allCatSums) <-  c("Y CS+ PE+", "Y CS+ PE-", "Y CS- PE+", "Y CS- PE-",
                            "Y CS< PE!", "Y CS- PE!", "Y CS! PE\u00B1",
                            "Y CS> PE\u00B1", "Y CS< PE\u00B1",
                            "N CS< PE+", "N CS- PE+")
    
    allCatSums[names(catSumTable)[nonCollapseIdxs]] <- catSumTable[nonCollapseIdxs]
    
    allCatSums["Y CS! PE\u00B1"] <- sum(catSumTable[which((names(catSumTable) == "Y CS! PE+") |
                                                            (names(catSumTable) == "Y CS! PE-"))], na.rm = TRUE)
    
    allCatSums["Y CS> PE\u00B1"] <- sum(catSumTable[which((names(catSumTable) == "Y CS> PE+") |
                                                            (names(catSumTable) == "Y CS> PE-"))], na.rm = TRUE)
    
    allCatSums["Y CS< PE\u00B1"] <- sum(catSumTable[which((names(catSumTable) == "Y CS< PE+") |
                                                            (names(catSumTable) == "Y CS< PE-"))], na.rm = TRUE)
    
    allCatSums["N CS< PE+"] <- sum(nrow(x$'N_CS<_PE+_ORFs'), na.rm = TRUE)
    allCatSums["N CS- PE+"] <- sum(nrow(x$'N_CS-_PE+_ORFs'), na.rm = TRUE)
    
    ## --------------------------------------------------------------------------------------------------------------- ##
    
    tenFactor <- 10 ^ floor(log10(max(allCatSums)))
    digitMult <- max(allCatSums) / tenFactor
    lastDigit <- floor(digitMult)
    
    if ((digitMult - lastDigit) <= 0.5) {
      currYMax <- (lastDigit + 0.5) * tenFactor
    } else {
      currYMax <- (lastDigit + 1) * tenFactor
    }
    
    ## --------------------------------------------------------------------------------------------------------------- ##
    
    catColors <- c("darkgreen", "lightgreen", "lightgreen", "white",
                   "darkred", "indianred1", "indianred1",
                   "gray", "gray",
                   "darkblue", "lightblue", "darkred")
    
    par(mfrow=c(1,1), mar = c(6, 4, 6, 2))
    
    barplot(allCatSums, las = 2, ylim = c(0, currYMax),
            col = catColors,
            family = "mono")
    
    title(main = plotTitle, line = 5)
    
    legend(x = ceiling(par("usr")[1]), y = 1.2 * currYMax, bty = "n", ncol = 3, xpd = TRUE,
           fill = c("darkgreen", "lightgreen", "white", "darkred", "indianred1", "gray", "darkblue", "lightblue"),
           legend = c("Definitely Correct", "Likely Correct", "No Evidence",
                      "Definitely Incorrect", "Likely Incorrect", "Potentially Incorrect",
                      "Likely Missing", "Potentially Missing"))
    
  } else {
    stop("'x' is an object of unrecognized subclass '", class(x)[2], "'.")
  }
  
  invisible(x)
}

#' @export
#' @importFrom graphics mosaicplot
#' @importFrom stats quantile
#'
#' @title Plot Genes by Category and Length
#' @description The \code{mosaicplot} method for \code{Assessment} object
#'
#' @param x An object of class \code{Assessment} and subclass \code{Results}.
#'
#' @param ... Further \code{mosaicplot} parameters.
#'
#' @details
#' \code{mosaicplot.Assessment} plots all the genes in the given \code{Results} object by category and length. This set of genes
#' includes both the supplied predicted genes as well as open reading frames with proteomics evidence but no predicted start.
#' 
#' The set of genes are separated into ten quantile bins based on the length of the gene/open reading frame. The genes are then
#' plotted by length bin and category in a mosaic format, with each column representing a length bin and each row/block
#' representing a category.
#' 
#' @return Invisibly returns the input object \code{x}
#'
#' @seealso \code{\link{Assessment-class}}
#' 
#' @examples
#'
#' mosaicplot(readRDS(system.file("extdata",
#'                                "MGAS5005_PreSaved_ResultsObj_Prodigal.rds",
#'                                package = "AssessORF")))
#'
mosaicplot.Assessment <- function(x, ...) {
  if (class(x)[1] != "Assessment") {
    stop("'x' must be an object of class 'Assessment'.")
  }
  
  if (class(x)[2] == "Results") {
    geneSource <- x$GeneSource
    
    part2 <- "Gene Category Assignments by Nucleotide Length"
    
    if (geneSource != "") {
      part2 <- paste(geneSource, part2)
    }
    
    ## --------------------------------------------------------------------------------------------------------------- ##
    
    speciesName <- x$Species
    strainID <- x$StrainID
    
    if ((speciesName != "") && (strainID != "")) {
      plotTitle <- bquote(italic(.(speciesName))~.(strainID)~.(part2))
    } else if ((speciesName != "")) {
      plotTitle <- bquote(italic(.(speciesName))~.(part2))
    } else if ((strainID != "")) {
      plotTitle <- bquote(~.(strainID)~.(part2))
    } else {
      plotTitle <- part2
    }
    
    ## --------------------------------------------------------------------------------------------------------------- ##
    
    geneLengths <- x$GeneRightPos - x$GeneLeftPos + 1
    names(geneLengths) <- x$CategoryAssignments
    
    if (!is.null(nrow(x$'N_CS<_PE+_ORFs'))){
      cat1BLens <- x$'N_CS<_PE+_ORFs'[, 'Length']
      names(cat1BLens) <- rep("N CS< PE+", nrow(x$'N_CS<_PE+_ORFs'))
      
      geneLengths <- c(geneLengths, cat1BLens)
    }
    
    if (!is.null(nrow(x$'N_CS-_PE+_ORFs'))){
      cat1ALens <- x$'N_CS-_PE+_ORFs'[, 'Length']
      names(cat1ALens) <- rep("N CS- PE+", nrow(x$'N_CS-_PE+_ORFs'))
      
      geneLengths <- c(geneLengths, cat1ALens)
    }
    
    quantBins <- quantile(geneLengths, seq(0, 1, 0.1))
    
    catByLenTable <- table(cut(geneLengths, quantBins), names(geneLengths))
    
    catCodeNames <- c("Y CS+ PE+", "Y CS+ PE-", "Y CS- PE+", "Y CS- PE-",
                      "Y CS< PE!", "Y CS- PE!", "Y CS! PE+", "Y CS! PE-",
                      "Y CS> PE+", "Y CS> PE-", "Y CS< PE+", "Y CS< PE-",
                      "N CS< PE+", "N CS- PE+", "Y CS> PE!")
    
    catColors <- c("darkgreen", "lightgreen", "lightgreen", "white",
                   "darkred", "indianred1", "indianred1", "indianred1",
                   "gray", "gray", "gray", "gray",
                   "darkblue", "lightblue", "darkred")
    
    catCodeNames <- catCodeNames[catCodeNames %in% colnames(catByLenTable)]
    
    usedCats <- which(colSums(catByLenTable[, catCodeNames], na.rm = TRUE) != 0)
    
    catCodeNames <- catCodeNames[usedCats]
    
    catColors <- catColors[usedCats]
    
    plot(catByLenTable[, catCodeNames], main = plotTitle, col = catColors, las = 1)
    
  } else {
    stop("'x' must be of subclass 'Results'.")
  }
  
  invisible(x)
}
