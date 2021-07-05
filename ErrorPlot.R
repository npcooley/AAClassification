# A script for plotting results similar to figure 2 in the manuscript:
# Accurate annotation of protein coding sequences with IDTAXA
# author: Nicholas Cooley
# contact: npc19@pitt.edu

# Requirements:
# R version > 4.0
# DECIPHER > 2.19
# gridExtra v ?

# Argument form:
# RScript <RDatafile1> <RDatafile2> <RDatafile3> <SEP> <RDatafile4> <RDatafileetc...>
# This script saves two outputs:
# errorplot.pdf
# errortable.pdf

# All results follow the same form and have the same performance summaries
# supports 1 panel, 2 panels, 3 panels, or 4 panels
# 1 panel format: 7v7 pdf
# 2 panel format: 2 stacked plots 7x3.5 pdf
# 3 panel format: 3 stacked plots 10.5x3.5 pdf
# 4 panel format: standard 4 panel figure 7x7 pdf
# if SEP appears once, two panels will be created with the data files on each size
# being binned together into panels
# if SEP appears twice, three panels will be generated in a similar manner and so on
# if any panel appears to be empty, script will complain and stop

# IDTaxa results are expected to contain an 'IPerf' object
# BLAST results are expected to contain an 'EPerf' object
# a 'PPerf' object, and a Res object
# HMMER results are expected to contain 'Res', and 'HPerf'

###### -- libraries -----------------------------------------------------------

suppressMessages(library(DECIPHER))
suppressMessages(library(gridExtra))

###### -- arguments -----------------------------------------------------------

Args <- commandArgs(trailingOnly = TRUE)

###### -- code body -----------------------------------------------------------
# argument form:
# RScript <scriptname> <CVdata.RData> <CVdata.RData> <SEP> <CVdata.RData> <etc>
# script will divide arguments on `SEP` to build panels of data to plot together
# data naming conventions `METHOD_Level_Other`
# both BLAST and HMMER are not expected to have an `other` category,
# BLASTP and BLASTN should be specified in METHOD
# all HMMER results are in AA space

# Hard code color vector, lifted from:
# Sasha's List of 20 simple distinct colors
# https://sashamaps.net/docs/resources/20-colors/
ColVector <- c('#e6194b',
               '#3cb44b',
               '#ffe119',
               '#4363d8',
               '#f58231',
               '#911eb4',
               '#46f0f0',
               '#f032e6',
               '#bcf60c',
               '#fabebe',
               '#008080',
               '#e6beff',
               '#9a6324',
               '#fffac8',
               '#800000',
               '#aaffc3',
               '#808000',
               '#ffd8b1',
               '#000075',
               '#808080',
               '#ffffff',
               '#000000')

PanelCounts <- sum(grepl(pattern = "SEP",
                         x = Args,
                         fixed = TRUE)) + 1L
PanelFiles <- sapply(seq(PanelCounts),
                     function(x) rep(NA,
                                     length(Args)),
                     simplify = FALSE,
                     USE.NAMES = FALSE)
L1 <- 1L
L2 <- 1L
for (m1 in seq_along(Args)) {
  if (Args[m1] == "SEP") {
    # argument is seperator,
    # move to next panel,
    # reset iterator
    L1 <- L1 + 1L
    L2 <- 1L
  } else if (Args[m1] != "SEP") {
    # argument isn't seperator, 
    # add file to current iterator
    PanelFiles[[L1]][L2] <- Args[m1]
    L2 <- L2 + 1L
  }
}

PanelFiles <- sapply(PanelFiles,
                     function(x) x[!is.na(x)],
                     USE.NAMES = FALSE,
                     simplify = FALSE)

PanelFiles <- PanelFiles[lengths(PanelFiles) > 0L]

if (length(PanelFiles) == 0L |
    length(PanelFiles) >= 5L) {
  stop ("Please check inputs, panel arrangement is not supported.")
}

# panel arrangement
# if 1:
# 1 par, both axes
# if 2:
# 2 pars, first with both axes, second with only x
# if 3:
# 2 pars, first with both axes, second with only x
# third with both
# if 4:
# 4 pars, first with y axis only, second with no axes
# third with both axes, fourth with x axis only

# 4 possible arrangements
MARLIST <- AXESLIST <- vector(mode = "list",
                              length = 4L)

# uniform across all plots
# par(mgp = c(2.2, 1.0, 0.0))

for (m1 in seq_along(MARLIST)) {
  if (m1 == 1L) {
    MARLIST[[m1]] <- list(c(3.5, 3.5, 0.5, 1.0)) # both x and y axes
    AXESLIST[[m1]] <- list(c("Classification Rate",
                             "Error Rate"))
  } else if (m1 == 2L) {
    MARLIST[[m1]] <- list(c(3.5, 3.5, 0.5, 1.0), # both x and y axes
                          c(3.5, 2.5, 0.5, 2.0)) # x axis only
    AXESLIST[[m1]] <- list(c("Classification Rate",
                             "Error Rate"),
                           c("",
                             "Error Rate"))
  } else if (m1 == 3L) {
    MARLIST[[m1]] <- list(c(3.5, 3.5, 0.5, 1.0), # both x and y axes
                          c(3.5, 2.5, 0.5, 2.0), # x axis only
                          c(3.5, 3.5, 0.5, 1.0)) # both x and y axes
    AXESLIST[[m1]] <- list(c("Classification Rate",
                             "Error Rate"),
                           c("",
                             "Error Rate"),
                           c("Classification Rate",
                             "Error Rate"))
  } else if (m1 == 4L) {
    MARLIST[[m1]] <- list(c(2.0, 3.5, 2.0, 1.0), # y axis only
                          c(2.0, 2.5, 2.0, 2.0), # no axes at all
                          c(3.5, 3.5, 0.5, 1.0), # x axis only
                          c(3.5, 2.5, 0.5, 2.0)) # both x and y axes
    AXESLIST[[m1]] <- list(c("",
                             "Error Rate"),
                           c("",
                             ""),
                           c("Classification Rate",
                             "Error Rate"),
                           c("Classification Rate",
                             ""))
  }
}

# for the table
# script will only create a two column table
# if a user wishes to create this format with different levels of a different heirarchy
# these must be edited to reflect the user's new schemes

###### -- plot stuff ----------------------------------------------------------

# create vector for table row names
# remove rows where nchar == 0 at end
TableRowNames <- vector(mode = "character",
                        length = length(Args) *2L)
KOOverConf <- KOMissConf <- KOOverClass <- KOMissClass <- rep(NA_real_,
                                                              length(Args) *2L)
TAXOverConf <- TAXMissConf <- TAXOverClass <- TAXMissClass <- rep(NA_real_,
                                                                  length(Args) *2L)
KOConfConv <- TAXConfConv <- rep(NA_real_,
                                 length(Args) *2L)

PanelParams <- length(PanelFiles)
# open graphics device
# set panel arrangement
if (PanelParams == 1L) {
  pdf("ErrorPlot.pdf",
      height = 7,
      width = 7)
} else if (PanelParams == 2L) {
  pdf("ErrorPlot.pdf",
      height = 7,
      width = 3.5)
  layout(mat = matrix(data = 1:2,
                      nrow = 2,
                      ncol = 1))
} else if (PanelParams == 3L) {
  pdf("ErrorPlot.pdf",
      height = 7,
      width = 7)
  layout(mat = matrix(data = c(1, 0, 2, 3),
                      nrow = 2,
                      ncol = 2))
} else if (PanelParams == 4L) {
  pdf("ErrorPlot.pdf",
      height = 7,
      width = 7)
  layout(mat = matrix(data = c(1,3,2,4),
                      nrow = 2,
                      ncol = 2))
}

# expect files to be named `Method_Level_Mode.RData`
# in HMMERs case this is Method_Level_.RData
# legend name will come from first and last _ separated characters
# table row will be designated first and last
# table column will be designated by middle

LegendVector <- vector(mode = "list",
                       length = length(PanelFiles))
for (m1 in seq_along(PanelFiles)) {
  
  # uniform across all plots
  par(mgp = c(2.2, 1.0, 0.0))
  par(mar = MARLIST[[PanelParams]][[m1]])
  plot(x = 0,
       y = 0,
       type = "n",
       xlim = c(0, 1),
       ylim = c(0, 1),
       xlab = AXESLIST[[PanelParams]][[m1]][1L],
       ylab = AXESLIST[[PanelParams]][[m1]][2L])
  LegendPos <- 1L
  LegendVector[[m1]] <- rep(NA_character_,
                            length(PanelFiles[[m1]]) * 2L)
  for (m2 in seq_along(PanelFiles[[m1]])) {
    # within each panel generate a legend vector
    # legend vector is generated from the names of the files
    # expecting a format of `ID_otherID_otherID_etc...`
    # name will be grabbed from first and last `_` separated characters
    FilesLoaded <- load(file = PanelFiles[[m1]][m2],
                        verbose = FALSE)
    if ("IPerf" %in% FilesLoaded) {
      # IDTaxa result
      # plot IDTAXA confidence, extract 60% confidence for table
      # create legend tag and assign color
      D <- get("IPerf")
      lines(x = D[[1]],
            y = D[[2]], # OC
            lty = 2,
            col = ColVector[LegendPos])
      lines(x = D[[1]],
            y = D[[3]], # MC
            lty = 1,
            col = ColVector[LegendPos])
      # points(x = D[[4]][c(1, 1), c(4, 4)],
      #        y = D[[4]][2:3, c(4, 4)],
      #        col = ColVector[LegendPos])
      ph <- strsplit(x = PanelFiles[[m1]][m2],
                     split = "/",
                     fixed = TRUE)[[1]]
      ph <- ph[length(ph)]
      ph <- gsub(pattern = ".RData",
                 replacement = "",
                 x = ph,
                 fixed = TRUE)
      ph1 <- strsplit(x = ph,
                      split = "_",
                      fixed = TRUE)[[1]][2]
      ph2 <- paste0(strsplit(x = ph,
                             split = "_",
                             fixed = TRUE)[[1]][c(1, 3)],
                    collapse = " ")
      LegendVector[[m1]][LegendPos] <- ph2
      if (ph2 %in% TableRowNames) {
        # row name has already been specified
        # ask which column(s) to fill in and add the err(s) at 60 % confidence
        # and 60 % classified
        w <- which(TableRowNames == ph2)
        if (ph1 == "KO") {
          w1 <- which.min(abs(D[[1]] - .6))
          KOOverClass[w] <- D[[2]][w1] # over class where 60 % of classifiables classified
          KOMissClass[w] <- D[[3]][w1] # miss class where 60 % of classifiables classified
          KOOverConf[w] <- D[[4]][2, 4] # over class where Conf is 60%
          KOMissConf[w] <- D[[4]][3, 4] # miss class where Conf is 60%
          KOConfConv[w] <- w1 * 0.1
        } else {
          w1 <- which.min(abs(D[[1]] - .6))
          TAXOverClass[w] <- D[[2]][w1] # over class where 60 % of classifiables classified
          TAXMissClass[w] <- D[[3]][w1] # miss class where 60 % of classifiables classified
          TAXOverConf[w] <- D[[4]][2, 4] # over class where Conf is 60%
          TAXMissConf[w] <- D[[4]][3, 4] # miss class where Conf is 60%
          TAXConfConv[w] <- w1 * 0.1
        }
      } else {
        # row name has not been specified, add it
        # ask which column(s) to fill in and add the err(s) at 60 % confidence
        # and 60 % classified
        # assign to first row that is yet to be occupied
        w <- which(nchar(TableRowNames) == 0)[1L]
        if (ph1 == "KO") {
          w1 <- which.min(abs(D[[1]] - .6))
          KOOverClass[w] <- D[[2]][w1] # over class where 60 % of classifiables classified
          KOMissClass[w] <- D[[3]][w1] # miss class where 60 % of classifiables classified
          KOOverConf[w] <- D[[4]][2, 4] # over class where Conf is 60%
          KOMissConf[w] <- D[[4]][3, 4] # miss class where Conf is 60%
          KOConfConv[w] <- w1 * 0.1
        } else {
          w1 <- which.min(abs(D[[1]] - .6))
          TAXOverClass[w] <- D[[2]][w1] # over class where 60 % of classifiables classified
          TAXMissClass[w] <- D[[3]][w1] # miss class where 60 % of classifiables classified
          TAXOverConf[w] <- D[[4]][2, 4] # over class where Conf is 60%
          TAXMissConf[w] <- D[[4]][3, 4] # miss class where Conf is 60%
          TAXConfConv[w] <- w1 * 0.1
        }
        TableRowNames[nchar(TableRowNames) == 0L][1L] <- ph2
      }
      LegendPos <- LegendPos + 1L
    } else if ("HPerf" %in% FilesLoaded) {
      # HMMER result
      # plot HMMER e-value as confidence
      # convert 60 % confidence as evalue back to evalue for table
      # create legend tag and assign color
      D <- get("HPerf")
      lines(x = D[[1]],
            y = D[[2]], # OC
            lty = 2,
            col = ColVector[LegendPos])
      lines(x = D[[1]],
            y = D[[3]], # MC
            lty = 1,
            col = ColVector[LegendPos])
      # points(x = D[[4]][c(1, 1), c(4, 4)],
      #        y = D[[4]][2:3, c(4, 4)],
      #        col = ColVector[LegendPos])
      ph <- strsplit(x = PanelFiles[[m1]][m2],
                     split = "/",
                     fixed = TRUE)[[1]]
      ph <- ph[length(ph)]
      ph <- gsub(pattern = ".RData",
                 replacement = "",
                 x = ph,
                 fixed = TRUE)
      ph1 <- strsplit(x = ph,
                      split = "_",
                      fixed = TRUE)[[1]][2]
      ph2 <- strsplit(x = ph,
                      split = "_",
                      fixed = TRUE)[[1]][1]
      LegendVector[[m1]][LegendPos] <- ph2
      if (ph2 %in% TableRowNames) {
        # row name has already been specified
        # ask which column(s) to fill in and add the err(s) at 60 % confidence
        # and 60 % classified
        w <- which(TableRowNames == ph2)
        x <- as.numeric(Res[, 5])
        x <- min(x[x > 0])
        if (ph1 == "KO") {
          w1 <- which.min(abs(D[[1]] - .6))
          KOOverClass[w] <- D[[2]][w1] # over class where 60 % of classifiables classified
          KOMissClass[w] <- D[[3]][w1] # miss class where 60 % of classifiables classified
          KOOverConf[w] <- D[[4]][2, 4] # over class where Conf is 60%
          KOMissConf[w] <- D[[4]][3, 4] # miss class where Conf is 60%
          KOConfConv[w] <- w1 * 0.001
          KOConfConv[w] <- 10^(log10(x * KOConfConv[w]))
        } else {
          w1 <- which.min(abs(D[[1]] - .6))
          TAXOverClass[w] <- D[[2]][w1] # over class where 60 % of classifiables classified
          TAXMissClass[w] <- D[[3]][w1] # miss class where 60 % of classifiables classified
          TAXOverConf[w] <- D[[4]][2, 4] # over class where Conf is 60%
          TAXMissConf[w] <- D[[4]][3, 4] # miss class where Conf is 60%
          TAXConfConv[w] <- w1 * 0.001
          TAXConfConv[w] <- 10^(log10(x * TAXConfConv[w]))
        }
      } else {
        # row name has not been specified, add it
        # ask which column(s) to fill in and add the err(s) at 60 % confidence
        # and 60 % classified
        # assign to first row that is yet to be occupied
        w <- which(nchar(TableRowNames) == 0)[1L]
        x <- as.numeric(Res[, 5])
        x <- min(x[x > 0])
        if (ph1 == "KO") {
          w1 <- which.min(abs(D[[1]] - .6))
          KOOverClass[w] <- D[[2]][w1] # over class where 60 % of classifiables classified
          KOMissClass[w] <- D[[3]][w1] # miss class where 60 % of classifiables classified
          KOOverConf[w] <- D[[4]][2, 4] # over class where Conf is 60%
          KOMissConf[w] <- D[[4]][3, 4] # miss class where Conf is 60%
          KOConfConv[w] <- w1 * 0.001
          KOConfConv[w] <- 10^(log10(x * KOConfConv[w]))
        } else {
          w1 <- which.min(abs(D[[1]] - .6))
          TAXOverClass[w] <- D[[2]][w1] # over class where 60 % of classifiables classified
          TAXMissClass[w] <- D[[3]][w1] # miss class where 60 % of classifiables classified
          TAXOverConf[w] <- D[[4]][2, 4] # over class where Conf is 60%
          TAXMissConf[w] <- D[[4]][3, 4] # miss class where Conf is 60%
          TAXConfConv[w] <- w1 * 0.001
          TAXConfConv[w] <- 10^(log10(x * TAXConfConv[w]))
        }
        TableRowNames[nchar(TableRowNames) == 0L][1L] <- ph2
      }
      LegendPos <- LegendPos + 1L
    } else if ("PPerf" %in% FilesLoaded) {
      # BLAST result
      # plot BLAST PID as confidence AND
      # BLAST e-value as confidence 
      # convert 60 % confidence as evalue back to evalue for table
      # create legend tag and assign color for both PID and EVL
      # EVL plotting and stuff
      D <- get("EPerf")
      lines(x = D[[1]],
            y = D[[2]], # OC
            lty = 2,
            col = ColVector[LegendPos])
      lines(x = D[[1]],
            y = D[[3]], # MC
            lty = 1,
            col = ColVector[LegendPos])
      # points(x = D[[4]][c(1, 1), c(4, 4)],
      #        y = D[[4]][2:3, c(4, 4)],
      #        col = ColVector[LegendPos])
      ph <- strsplit(x = PanelFiles[[m1]][m2],
                     split = "/",
                     fixed = TRUE)[[1]]
      ph <- ph[length(ph)]
      ph <- gsub(pattern = ".RData",
                 replacement = "",
                 x = ph,
                 fixed = TRUE)
      ph1 <- strsplit(x = ph,
                      split = "_",
                      fixed = TRUE)[[1]][2]
      ph2 <- strsplit(x = ph,
                      split = "_",
                      fixed = TRUE)[[1]][1]
      ph3 <- paste(ph2, "EVL", sep = " ")
      LegendVector[[m1]][LegendPos] <- ph3
      if (ph3 %in% TableRowNames) {
        # row name has already been specified
        # ask which column(s) to fill in and add the err(s) at 60 % confidence
        # and 60 % classified
        w <- which(TableRowNames == ph3)
        if (ph1 == "KO") {
          w1 <- which.min(abs(D[[1]] - .6))
          KOOverClass[w] <- D[[2]][w1] # over class where 60 % of classifiables classified
          KOMissClass[w] <- D[[3]][w1] # miss class where 60 % of classifiables classified
          KOOverConf[w] <- D[[4]][2, 4] # over class where Conf is 60%
          KOMissConf[w] <- D[[4]][3, 4] # miss class where Conf is 60%
          KOConfConv[w] <- w1 * 0.001
          KOConfConv[w] <- 10^(log10(min(Res$evalue[Res$evalue > 0]) * KOConfConv[w]))
        } else {
          w1 <- which.min(abs(D[[1]] - .6))
          TAXOverClass[w] <- D[[2]][w1] # over class where 60 % of classifiables classified
          TAXMissClass[w] <- D[[3]][w1] # miss class where 60 % of classifiables classified
          TAXOverConf[w] <- D[[4]][2, 4] # over class where Conf is 60%
          TAXMissConf[w] <- D[[4]][3, 4] # miss class where Conf is 60%
          TAXConfConv[w] <- w1 * 0.001
          KOConfConv[w] <- 10^(log10(min(Res$evalue[Res$evalue > 0]) * TAXConfConv[w]))
        }
      } else {
        # row name has not been specified, add it
        # ask which column(s) to fill in and add the err(s) at 60 % confidence
        # and 60 % classified
        # assign to first row that is yet to be occupied
        w <- which(nchar(TableRowNames) == 0)[1L]
        if (ph1 == "KO") {
          w1 <- which.min(abs(D[[1]] - .6))
          KOOverClass[w] <- D[[2]][w1] # over class where 60 % of classifiables classified
          KOMissClass[w] <- D[[3]][w1] # miss class where 60 % of classifiables classified
          KOOverConf[w] <- D[[4]][2, 4] # over class where Conf is 60%
          KOMissConf[w] <- D[[4]][3, 4] # miss class where Conf is 60%
          KOConfConv[w] <- w1 * 0.001
          KOConfConv[w] <- 10^(log10(min(Res$evalue[Res$evalue > 0]) * KOConfConv[w]))
        } else {
          w1 <- which.min(abs(D[[1]] - .6))
          TAXOverClass[w] <- D[[2]][w1] # over class where 60 % of classifiables classified
          TAXMissClass[w] <- D[[3]][w1] # miss class where 60 % of classifiables classified
          TAXOverConf[w] <- D[[4]][2, 4] # over class where Conf is 60%
          TAXMissConf[w] <- D[[4]][3, 4] # miss class where Conf is 60%
          TAXConfConv[w] <- w1 * 0.001
          TAXConfConv[w] <- 10^(log10(min(Res$evalue[Res$evalue > 0]) * TAXConfConv[w]))
        }
        TableRowNames[nchar(TableRowNames) == 0L][1L] <- ph3
      }
      
      # PID plotting and stuff
      D <- get("PPerf")
      lines(x = D[[1]],
            y = D[[2]], # OC
            lty = 2,
            col = ColVector[LegendPos + 1L])
      lines(x = D[[1]],
            y = D[[3]], # MC
            lty = 1,
            col = ColVector[LegendPos + 1L])
      # points(x = D[[4]][c(1, 1), c(4, 4)],
      #        y = D[[4]][2:3, c(4, 4)],
      #        col = ColVector[LegendPos + 1L])
      ph3 <- paste(ph2, "PID", sep = " ")
      LegendVector[[m1]][LegendPos + 1L] <- ph3
      if (ph3 %in% TableRowNames) {
        # row name has already been specified
        # ask which column(s) to fill in and add the err(s) at 60 % confidence
        # and 60 % classified
        w <- which(TableRowNames == ph3)
        # x <- as.numeric(Res[, 5])
        # x <- min(x[x > 0])
        if (ph1 == "KO") {
          w1 <- which.min(abs(D[[1]] - .6))
          KOOverClass[w] <- D[[2]][w1] # over class where 60 % of classifiables classified
          KOMissClass[w] <- D[[3]][w1] # miss class where 60 % of classifiables classified
          KOOverConf[w] <- D[[4]][2, 4] # over class where Conf is 60%
          KOMissConf[w] <- D[[4]][3, 4] # miss class where Conf is 60%
          KOConfConv[w] <- w1 * 0.001
        } else {
          w1 <- which.min(abs(D[[1]] - .6))
          TAXOverClass[w] <- D[[2]][w1] # over class where 60 % of classifiables classified
          TAXMissClass[w] <- D[[3]][w1] # miss class where 60 % of classifiables classified
          TAXOverConf[w] <- D[[4]][2, 4] # over class where Conf is 60%
          TAXMissConf[w] <- D[[4]][3, 4] # miss class where Conf is 60%
          TAXConfConv[w] <- w1 * 0.001
        }
      } else {
        # row name has not been specified, add it
        # ask which column(s) to fill in and add the err(s) at 60 % confidence
        # and 60 % classified
        # assign to first row that is yet to be occupied
        w <- which(nchar(TableRowNames) == 0)[1L]
        if (ph1 == "KO") {
          w1 <- which.min(abs(D[[1]] - .6))
          KOOverClass[w] <- D[[2]][w1] # over class where 60 % of classifiables classified
          KOMissClass[w] <- D[[3]][w1] # miss class where 60 % of classifiables classified
          KOOverConf[w] <- D[[4]][2, 4] # over class where Conf is 60%
          KOMissConf[w] <- D[[4]][3, 4] # miss class where Conf is 60%
          KOConfConv[w] <- w1 * 0.001
        } else {
          w1 <- which.min(abs(D[[1]] - .6))
          TAXOverClass[w] <- D[[2]][w1] # over class where 60 % of classifiables classified
          TAXMissClass[w] <- D[[3]][w1] # miss class where 60 % of classifiables classified
          TAXOverConf[w] <- D[[4]][2, 4] # over class where Conf is 60%
          TAXMissConf[w] <- D[[4]][3, 4] # miss class where Conf is 60%
          TAXConfConv[w] <- w1 * 0.001
        }
        TableRowNames[nchar(TableRowNames) == 0L][1L] <- ph3
      }
      
      LegendPos <- LegendPos + 2L
    } # data type logicals
    
    rm(list = FilesLoaded)
  } # loop through intra panel data
  LegendVector[[m1]] <- LegendVector[[m1]][!is.na(LegendVector[[m1]])]
  L <- legend(x = 0.005,
              y = 0.9,
              legend = rep(NA, length(LegendVector[[m1]]) * 2),
              col = rep(ColVector[seq(length(LegendVector[[m1]]))],
                        2),
              lty = c(rep(1, length(LegendVector[[m1]])),
                      rep(2, length(LegendVector[[m1]]))),
              pch = rep(NA_integer_, length(length(LegendVector[[m1]])) * 2L),
              ncol = 2,
              bty = 'n',
              x.intersp = 0.5,
              inset = 0.02,
              cex = 0.75)
  legend(x = L$rect$left,
         y = L$rect$top,
         legend = LegendVector[[m1]],
         col = rep(NA,2),
         lty = c(1,1),
         ncol = 1,
         x.intersp = 3,
         bg = NA,
         bty = "n",
         cex = 0.75)
  mtext("MC",
        cex = 0.55,
        at = 0.055,
        line = -2.65)
  mtext("OC",
        cex = 0.55,
        at = 0.155,
        line = -2.65)
  # explanation of bquote behavior
  # https://stackoverflow.com/questions/15074127/use-expression-with-a-variable-r
  mtext(bquote(paste(bold(.(LETTERS[m1])))),
        cex = 1.25,
        at = 0.055,
        line = -1.65)
  cat("\nPanel Completed!\n")
} # loops through panels

# close graphics device
dev.off()

RowsOccupied <- nchar(TableRowNames) > 0
TableRowNames <- TableRowNames[RowsOccupied]
KOOverClass <- KOOverClass[RowsOccupied]
KOMissClass <- KOMissClass[RowsOccupied]
TAXMissClass <- TAXMissClass[RowsOccupied]
TAXOverClass <- TAXOverClass[RowsOccupied]

KOOverConf <- KOOverConf[RowsOccupied]
KOMissConf <- KOMissConf[RowsOccupied]
TAXMissConf <- TAXMissConf[RowsOccupied]
TAXOverConf <- TAXOverConf[RowsOccupied]

ClassTable <- cbind(KOOverClass,
                    KOMissClass,
                    TAXOverClass,
                    TAXMissClass)

ConfTable <- cbind(KOOverConf,
                   KOMissConf,
                   TAXOverConf,
                   TAXMissConf)

rownames(ClassTable) <- rownames(ConfTable) <- TableRowNames
colnames(ClassTable) <- colnames(ConfTable) <- c("KO-OC",
                                                 "KO-MC",
                                                 "TAX-OC",
                                                 "TAX-MC")

pdf(file = "ErrByConf60.pdf",
    height = 4,
    width = 7)
grid.table(ConfTable)
dev.off()

pdf(file = "ErrByClass60.pdf",
    height = 4,
    width = 7)
grid.table(ClassTable)
dev.off()



