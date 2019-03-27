\name{Assessment}
\alias{Assessment-class}

\title{Assessment objects}

\description{
In order to assess the quality of a set of (predicted) genes for a genome, evidence must first be mapped to that genome.
Next, each gene must be categorized based on how strong the evidence is for that gene or against that gene. Class \code{Assessment}
furnishes objects that can store the necessary information for assessing a set of genes for a genome and also provides
functions for viewing and visualizing assessment information. Specifically, class \code{Assessment} objects utilize proteomic hits
and evolutionarily conserved start & stop codons as evidence to determine the correctness for each gene in a given set.
}

\section{\code{DataMap} Objects}{
  Objects of class \code{Assessment} and subclass \code{DataMap} are used to store the mapping of proteomics and evolutionary
  conservation to the genome of interest (central genome). They are generated through the function \code{\link{MapAssessmentData}},
  and they have a list structure containing the following elements:
  
  \describe{
    \item{\code{StrainID}}{Equal to \code{strainID} if it was specified; otherwise ""}
    \item{\code{Species}}{Equal to \code{speciesName} if it was specified; otherwise ""}
    \item{\code{GenomeLength}}{Length of the central genome}
    \item{\code{StopsByFrame}}{Where the stops are in each frame, used to bound open reading frames in downstream functions}
    \item{\code{N-TermProteomics}}{Logical describing whether or not the proteomics hits are from N-terminal proteomics}
    \item{\code{FwdProtHits}}{Proteomic hit information that maps to the three forward frames of the central genome}
    \item{\code{RevProtHits}}{Proteomic hit information that maps to the three reverse frames of the central genome}
    \item{\code{FwdCoverage}}{Coverage of the forward strand of the central genome}
    \item{\code{FwdConStarts}}{Start codon conservation of the forward strand of the central genome}
    \item{\code{FwdConStops}}{Stop codon conservation of the forward strand of the central genome}
    \item{\code{RevCoverage}}{Coverage of the reverse strand of the central genome}
    \item{\code{RevConStarts}}{Start codon conservation of the reverse strand of the central genome}
    \item{\code{RevConStops}}{Stop codon conservation of the reverse strand of the central genome}
    \item{\code{NumRelatedGenomes}}{Final number of related genomes that were mapped to the central genome}
    \item{\code{HasProteomics}}{Logical describing whether or not proteomics evidence has been mapped to the central genome}
    \item{\code{HasConservation}}{Logical describing whether or not evolutionary conservation evidence has been mapped to the central genome}
  }
}

\section{\code{Results} Objects}{
  Objects of class \code{Assessment} and subclass \code{Results} are used to store how correct a set of genes for a given genome.
  The function \code{\link{AssessGenes}} generates \code{Results} using a \code{DataMap} object and information on a set of genes
  for the genome corresponding to the \code{DataMap} object. \code{Results} objects have a list structure containing the following
  elements:
  
  \describe{
    \item{\code{StrainID}}{Equal to the \code{strainID} of the corresponding \code{DataMap} object}
    \item{\code{Species}}{Equal to \code{speciesName} of the corresponding \code{DataMap} object}
    \item{\code{GenomeLength}}{Length of the genome}
    \item{\code{GeneLeftPos}}{Left positions of the given set of genes (in forward strand terms)}
    \item{\code{GeneRightPos}}{Right positions of the given set of genes (in forward strand terms)}
    \item{\code{GeneStrand}}{Strand information of the given set of genes ("+" or "-")}
    \item{\code{GeneSource}}{The source of the given set of genes}
    \item{\code{NumGenes}}{Number of genes given}
    \item{\code{N_CS-_PE+_ORFs}}{Data for open reading frames with no gene start but with proteomics evidence}
    \item{\code{N_CS<_PE+_ORFs}}{Data for open reading frames with no gene start but with proteomics evidence and at least one
    valid evolutionarily conserved start}
    \item{\code{CategoryAssignments}}{A character vector that stores the category assignment for each of the given genes in the
    same order as the gene information (please see below for a list of all possible categories, their descriptions, and their
    character string codes)}
  }
}

\section{Gene Categories}{
  The \code{CategoryAssignments} vector in \code{Results} objects describes how the proteomics evidence and evolutionarily
  conserved start/stop codon evidence support or disprove the corresponding set of genes. In the vector, each gene is assigned
  a character string code that has the following format: "Y CS[_] PE[_]". The first part, the "Y", signifies that for this ORF
  contains a predicted gene. The second part, the "CS[_]", describes how the conserved start(s) lines up with the given gene
  start. The third part, the "PE[_]", describes how the proteomics hits line up with the given gene start.
  
  \describe{
    \item{\code{Y CS+ PE+}}{There is a good conserved start aligned with the gene start with protein evidence downstream.}
    \item{\code{Y CS+ PE-}}{There is a good conserved start aligned with the gene start without protein evidence downstream.}
    \item{\code{Y CS- PE+}}{There is no good conserved start aligned with the predicted start, and there is protein evidence downstream of the gene start.}
    \item{\code{Y CS- PE-}}{There is no good conserved start aligned with the predicted start, and there is no protein evidence downstream of the gene start.}
    \item{\code{Y CS! PE-}}{There are either multiple good conserved stops in the middle of the gene, or the most downstream, good conserved stop is followed by a good conserved start. There is no protein evidence downstream of the gene start}
    \item{\code{Y CS! PE+}}{The most downstream, good conserved stop is followed by a good conserved start, and there is protein evidence downstream of the gene start.}
    \item{\code{Y CS< PE!}}{The protein evidence disagrees with/is upstream of the gene start, and there is a good conserved start upstream of the protein evidence.}
    \item{\code{Y CS- PE!}}{The protein evidence disagrees with/is upstream of the gene start, and there is no good conserved start upstream of the protein evidence.}
    \item{\code{Y CS> PE+}}{The best conserved starts are downstream of the predicted start, and there is protein evidence downstream of the gene start.}
    \item{\code{Y CS> PE-}}{The best conserved starts are downstream of the predicted start, and there is no protein evidence downstream of the gene start.}
    \item{\code{Y CS< PE+}}{At least one of the best conserved starts is upstream of the predicted start, and there is protein evidence downstream of the gene start.}
    \item{\code{Y CS< PE-}}{At least one of the best conserved starts is upstream of the predicted start, and there is no protein evidence downstream of the gene start.}
  }
}

\section{S3 Methods}{
  \code{\link{as.matrix.Assessment}} (only works with objects of class \code{Results})

  \code{\link{print.Assessment}}

  \code{\link{plot.Assessment}}

  \code{\link{mosaicplot.Assessment}} (only works with objects of class \code{Results})
}