<<<<<<< HEAD
PlotAssessmentMapping <- function(mapObj,
                                 predObj = NULL,
                                 minConStart) {
  genomeLength <- mapObj$GenomeLength
  stops <- mapObj$StopsByFrame

  fwdProt <- mapObj$FwdProtHits
  revProt <- mapObj$RevProtHits

  fwdCov <- mapObj$FwdCoverage
  fwdConStart <- mapObj$FwdConStarts

  revCov <- mapObj$RevCoverage
  revConStart <- mapObj$RevConStarts

  addPred <- FALSE

  if (!is.null(predObj)) {
    addPred <- TRUE

    predLeft <- predObj$GeneLeftPos
    predRight <- predObj$GeneRightPos
    predStrand <- predObj$GeneStrand
  }

  speciesName <- mapObj$Species
  strainID <- mapObj$StrainID

  if ((speciesName != "") && (strainID != "")) {
    plotTitle <- bquote(italic(.(speciesName))~.(strainID)~"Genome Viewer")
  } else if ((speciesName != "")) {
    plotTitle <- bquote(italic(.(speciesName))~"Genome Viewer")
  } else if ((strainID != "")) {
    plotTitle <- bquote(~.(strainID)~"Genome Viewer")
  } else {
    plotTitle <- "Genome Viewer"
  }

  ## This function determines the saturation for proteomic hit blocks based on their scores.
  chooseColor <- function(color,
                          protScores,
                          start,
                          end,
                          max=20) {

    stopifnot(length(start)==length(end),
              color > 1 && color < 5,
              color==floor(color),
              all(start <= length(protScores)),
              all(end <= length(protScores)),
              all(start > 0),
              all(end > 0),
              all(end >= start))

    maxes <- integer(length(start))

    for (i in seq_along(start)) {
      maxes[i] <- max(protScores[start[i]:end[i]])
      if (maxes[i] > max)
        maxes[i] <- max
      if (maxes[i] < 0)
        maxes[i] <- 0
    }

    if (color==2) {
      colors <- rgb(1 - maxes/max, 1 - maxes/max, 1)
    } else if (color==3) {
      colors <- rgb(1, 1 - maxes/max, 1 - maxes/max)
    } else if (color==4) {
      colors <- rgb(1 - maxes/max, 1, 1 - maxes/max)
    }

    return(colors)
  }

  plotRange <- function(fwdRange) {
    revRange <- genomeLength - fwdRange + 1L

    layout(matrix(1:2))

    for (frame in c(4:6, 1:3)) {

      if (frame==1L) {

        par(mar=c(4.1, 4.1, 0.1, 0.5))

        plot(NA,
             xlim=range(fwdRange),
             ylim= c(0, 3),
             xlab= "Position",
             ylab= "Reading frame (forward)",
             xaxs= "i",
             yaxs= "i",
             yaxt= "n")
        axis(2, 0:2 + 0.5, 1:3)

        if (addPred) {
          # add predicted genes
          predFwd <- which(predStrand == "+")
          abline(v = predLeft[predFwd], col = "magenta")
          abline(v = predRight[predFwd], col = "cyan")
        }

        abline(h = 0:3)
      } else if (frame == 4L) {

        par(mar=c(2.1, 4.1, 2.1, 0.5))

        plot(NA,
             xlim= rev(range(revRange)),
             ylim= c(0, 3),
             xlab= "",
             ylab= "Reading frame (reverse)",
             main= plotTitle,
             xaxs= "i",
             yaxs= "i",
             yaxt= "n")
        axis(2, 0:2 + 0.5, 1:3)

        if (addPred) {
          # add predicted genes
          predRev <- which(predStrand == "-")
          abline(v = genomeLength - predLeft[predRev] + 1, col="cyan")
          abline(v = genomeLength - predRight[predRev] + 1, col="magenta")
        }

        abline(h = 0:3)
      }
      if (frame <= 3) {
        range <- fwdRange

        y <- fwdProt[[frame]][[seq]][range] > 0

        w1 <- which(diff(c(0, y))==1) # start
        w2 <- which(diff(c(y, 0))==-1) # end

        if (length(w1) > 0) {
          rect(range[w1] - 0.5, frame - 1, range[w2] + 0.5, frame,
               col=chooseColor(frame + 1, fwdProt[[frame]][[seq]][range], w1, w2),
               border=NA)
        }

        if (length(stops[[frame]]) > 0) {
          rect(stops[[frame]], frame - 1, stops[[frame]] + 2, frame,
               col="yellow", border=NA)
        }

        w <- which(fwdConStart[range]/fwdCov[range] > minConStart &
                     ((range - frame) %% 3)==0)

        if (length(w) > 0) {
          segments(range[w], frame - 1, range[w], frame,
                   col=gray(1 - fwdConStart[range[w]]/fwdCov[range[w]], alpha=0.5))
        }
      } else {
        range <- revRange

        y <- revProt[[frame - 3]][[seq]][range] > 0

        w1 <- which(diff(c(0, y))==1) # start
        w2 <- which(diff(c(y, 0))==-1) # end

        if (length(w1) > 0){
          rect(range[w1] - 0.5, frame - 4, range[w2] + 0.5, frame - 3,
               col=chooseColor(frame - 2, revProt[[frame - 3]][[seq]][range], w1, w2),
               border=NA)
        }

        if (length(stops[[frame]]) > 0){
          rect(stops[[frame]], frame - 4, stops[[frame]] + 2, frame - 3,
               col="yellow", border=NA)
        }

        w <- which(revConStart[range]/revCov[range] > minConStart &
                     ((range - frame - 3) %% 3)==0)

        if (length(w) > 0) {
          segments(range[w], frame - 4, range[w], frame - 3,
                   col = gray(1 - revConStart[range[w]]/revCov[range[w]], alpha=0.5))
        }
      }
    }
  }

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

  seq <- 1L
  l <- list(x=integer(4))

  cat("How to Interact With the Genome Viewer Via the Locator:\n",
      "Click the graphics window one or more times and then terminate the locator.\n",
      "1 click : Scroll the viewer to the left or the right.\n",
      "2 clicks: Zoom into the horizontal range between the two click points.\n",
      "3 clicks: Zoom out 10-fold.\n",
      "4 clicks: Zoom out completely and view the entire genome.\n",
      "To stop interaction, click zero times then terminate the locator.\n",
      "Depending on the graphical device, terminating the locator can either be ",
      "done by pressing the 'Stop' button, hitting the 'Esc' key, or right-clicking.\n")

  repeat {
    if (length(l$x)==0) { # no clicks
      break
    } else if (length(l$x)==1) { # scroll
      temp <- list(x=par("usr")[1:2])
      mid <- (temp$x[1] + temp$x[2])/2
      if (l$x - mid > 0) { # move right on forward
        temp$x <- temp$x + (temp$x[2] - temp$x[1])/2
      } else { # move left on forward
        temp$x <- temp$x - (temp$x[2] - temp$x[1])/2
      }
      plotRange(checkRange(c(temp$x[1], temp$x[2])))
    } else if (length(l$x)==2) { # zoom to range
      plotRange(checkRange(c(l$x[1], l$x[2])))
    } else if (length(l$x)==3) { # zoom out 10-fold
      temp <- par("usr")[1:2]
      dtemp <- diff(temp)
      temp <- c(temp[1] - dtemp*4.5,
                temp[2] + dtemp*4.5)
      plotRange(checkRange(temp))
    } else { # plot all
      plotRange(seq_len(genomeLength))
    }
    l <- locator()
  }
  
  par(mfrow=c(1,1))

  invisible(mapObj)
}
=======
PlotAssessmentMapping <- function(mapObj,
                                 predObj = NULL,
                                 minConStart) {
  genomeLength <- mapObj$GenomeLength
  stops <- mapObj$StopsByFrame

  fwdProt <- mapObj$FwdProtHits
  revProt <- mapObj$RevProtHits

  fwdCov <- mapObj$FwdCoverage
  fwdConStart <- mapObj$FwdConStarts

  revCov <- mapObj$RevCoverage
  revConStart <- mapObj$RevConStarts

  addPred <- FALSE

  if (!is.null(predObj)) {
    addPred <- TRUE

    predLeft <- predObj$GeneLeftPos
    predRight <- predObj$GeneRightPos
    predStrand <- predObj$GeneStrand
  }

  speciesName <- mapObj$Species
  strainID <- mapObj$StrainID

  if ((speciesName != "") && (strainID != "")) {
    plotTitle <- bquote(italic(.(speciesName))~.(strainID)~"Genome Viewer")
  } else if ((speciesName != "")) {
    plotTitle <- bquote(italic(.(speciesName))~"Genome Viewer")
  } else if ((strainID != "")) {
    plotTitle <- bquote(~.(strainID)~"Genome Viewer")
  } else {
    plotTitle <- "Genome Viewer"
  }

  ## This function determines the saturation for proteomic hit blocks based on their scores.
  chooseColor <- function(color,
                          protScores,
                          start,
                          end,
                          max=20) {

    stopifnot(length(start)==length(end),
              color > 1 && color < 5,
              color==floor(color),
              all(start <= length(protScores)),
              all(end <= length(protScores)),
              all(start > 0),
              all(end > 0),
              all(end >= start))

    maxes <- integer(length(start))

    for (i in seq_along(start)) {
      maxes[i] <- max(protScores[start[i]:end[i]])
      if (maxes[i] > max)
        maxes[i] <- max
      if (maxes[i] < 0)
        maxes[i] <- 0
    }

    if (color==2) {
      colors <- rgb(1 - maxes/max, 1 - maxes/max, 1)
    } else if (color==3) {
      colors <- rgb(1, 1 - maxes/max, 1 - maxes/max)
    } else if (color==4) {
      colors <- rgb(1 - maxes/max, 1, 1 - maxes/max)
    }

    return(colors)
  }

  plotRange <- function(fwdRange) {
    revRange <- genomeLength - fwdRange + 1L

    layout(matrix(1:2))

    for (frame in c(4:6, 1:3)) {

      if (frame==1L) {

        par(mar=c(4.1, 4.1, 0.1, 0.5))

        plot(NA,
             xlim=range(fwdRange),
             ylim= c(0, 3),
             xlab= "Position",
             ylab= "Reading frame (forward)",
             xaxs= "i",
             yaxs= "i",
             yaxt= "n")
        axis(2, 0:2 + 0.5, 1:3)

        if (addPred) {
          # add predicted genes
          predFwd <- which(predStrand == "+")
          abline(v = predLeft[predFwd], col = "magenta")
          abline(v = predRight[predFwd], col = "cyan")
        }

        abline(h = 0:3)
      } else if (frame == 4L) {

        par(mar=c(2.1, 4.1, 2.1, 0.5))

        plot(NA,
             xlim= rev(range(revRange)),
             ylim= c(0, 3),
             xlab= "",
             ylab= "Reading frame (reverse)",
             main= plotTitle,
             xaxs= "i",
             yaxs= "i",
             yaxt= "n")
        axis(2, 0:2 + 0.5, 1:3)

        if (addPred) {
          # add predicted genes
          predRev <- which(predStrand == "-")
          abline(v = genomeLength - predLeft[predRev] + 1, col="cyan")
          abline(v = genomeLength - predRight[predRev] + 1, col="magenta")
        }

        abline(h = 0:3)
      }
      if (frame <= 3) {
        range <- fwdRange

        y <- fwdProt[[frame]][[seq]][range] > 0

        w1 <- which(diff(c(0, y))==1) # start
        w2 <- which(diff(c(y, 0))==-1) # end

        if (length(w1) > 0) {
          rect(range[w1] - 0.5, frame - 1, range[w2] + 0.5, frame,
               col=chooseColor(frame + 1, fwdProt[[frame]][[seq]][range], w1, w2),
               border=NA)
        }

        if (length(stops[[frame]]) > 0) {
          rect(stops[[frame]], frame - 1, stops[[frame]] + 2, frame,
               col="yellow", border=NA)
        }

        w <- which(fwdConStart[range]/fwdCov[range] > minConStart &
                     ((range - frame) %% 3)==0)

        if (length(w) > 0) {
          segments(range[w], frame - 1, range[w], frame,
                   col=gray(1 - fwdConStart[range[w]]/fwdCov[range[w]], alpha=0.5))
        }
      } else {
        range <- revRange

        y <- revProt[[frame - 3]][[seq]][range] > 0

        w1 <- which(diff(c(0, y))==1) # start
        w2 <- which(diff(c(y, 0))==-1) # end

        if (length(w1) > 0){
          rect(range[w1] - 0.5, frame - 4, range[w2] + 0.5, frame - 3,
               col=chooseColor(frame - 2, revProt[[frame - 3]][[seq]][range], w1, w2),
               border=NA)
        }

        if (length(stops[[frame]]) > 0){
          rect(stops[[frame]], frame - 4, stops[[frame]] + 2, frame - 3,
               col="yellow", border=NA)
        }

        w <- which(revConStart[range]/revCov[range] > minConStart &
                     ((range - frame - 3) %% 3)==0)

        if (length(w) > 0) {
          segments(range[w], frame - 4, range[w], frame - 3,
                   col = gray(1 - revConStart[range[w]]/revCov[range[w]], alpha=0.5))
        }
      }
    }
  }

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

  seq <- 1L
  l <- list(x=integer(4))

  cat("How to Interact With the Genome Viewer Via the Locator:\n",
      "Click the graphics window one or more times and then terminate the locator.\n",
      "1 click : Scroll the viewer to the left or the right.\n",
      "2 clicks: Zoom into the horizontal range between the two click points.\n",
      "3 clicks: Zoom out 10-fold.\n",
      "4 clicks: Zoom out completely and view the entire genome.\n",
      "To stop interaction, click zero times then terminate the locator.\n",
      "Depending on the graphical device, terminating the locator can either be ",
      "done by pressing the 'Stop' button, hitting the 'Esc' key, or right-clicking.\n")

  repeat {
    if (length(l$x)==0) { # no clicks
      break
    } else if (length(l$x)==1) { # scroll
      temp <- list(x=par("usr")[1:2])
      mid <- (temp$x[1] + temp$x[2])/2
      if (l$x - mid > 0) { # move right on forward
        temp$x <- temp$x + (temp$x[2] - temp$x[1])/2
      } else { # move left on forward
        temp$x <- temp$x - (temp$x[2] - temp$x[1])/2
      }
      plotRange(checkRange(c(temp$x[1], temp$x[2])))
    } else if (length(l$x)==2) { # zoom to range
      plotRange(checkRange(c(l$x[1], l$x[2])))
    } else if (length(l$x)==3) { # zoom out 10-fold
      temp <- par("usr")[1:2]
      dtemp <- diff(temp)
      temp <- c(temp[1] - dtemp*4.5,
                temp[2] + dtemp*4.5)
      plotRange(checkRange(temp))
    } else { # plot all
      plotRange(seq_len(genomeLength))
    }
    l <- locator()
  }
  
  par(mfrow=c(1,1))

  invisible(mapObj)
}
>>>>>>> master-holder
