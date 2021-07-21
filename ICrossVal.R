# A script for generating the results that are shown in Figure 2 in the manuscript:
# Accurate annotation of protein coding sequences with IDTAXA
# author: Nicholas Cooley
# contact: npc19@pitt.edu

# Requirements:
# R version > 4.0
# DECIPHER > 2.19
# BLAST - please see BLAST installation instructions
# HMMER  - please see HMMER installation instructions
# R must have access to these executables
# which means that they must be in the PATH

# This script takes in two arguments
# 1) an RData file that contains 3 associated arguments
# 1a) an XStringSet names Seqs
# 1b) a list named Holdouts
# 1c) a list named TaxVectors
# 2) an integer
# 3) an output file name

# the usage of this script should look like:
# $Rscript <thisscript.R> <my_rdatafile.RData> <some_integer> <myresult.RData>

# Seqs is an xstringset that has unique KEGG IDs for sequence names

# Holdouts is a list of target holdout sequences for cross validation.
# In the supplied files, Holdouts has 5 initial entries, one for every tax level
# that can be tested, each level has 10 cvs which are represented by integers that
# relate to the positions in TaxVectors

# TaxVectors is a list that contains the taxonomies for all sequences present in Seqs.
# In the provided files, TaxVector is a list of length 5 list entry is populated
# by taxonomies for Seqs with increasing levels of assignment, i.e. TaxVectors[[1]]
# is simply the assigned KO number, where TaxVectors[[2]] takes the form of:
# "Root; KO_Number; Highest_Level_KEGG_Tax_Assignment", and so on

###### -- Libraries -----------------------------------------------------------

suppressMessages(library(DECIPHER))

###### -- Arguments -----------------------------------------------------------

Args <- commandArgs(trailingOnly = TRUE)
InFile <- Args[1L]
TargetTax <- as.integer(Args[2L])
OutFile <- Args[3L]

if (length(Args) < 4) {
  UseStandardAlphabet <- FALSE
} else {
  UseStandardAlphabet <- as.logical(Args[4L])
}


###### -- Functions -----------------------------------------------------------

OC_MC_FromIDTaxa <- function(ResList,
                             ResPositions = NULL,
                             Taxonomy,
                             EnforceSingleton = NULL,
                             Level = NULL,
                             ShowPlot = FALSE,
                             IncludeConf = FALSE,
                             ConfPoints = seq(10, 100, 10),
                             Verbose = FALSE) {
  if (Verbose) {
    FunStart <- Sys.time()
    pBar <- txtProgressBar(style = 1L)
  }
  
  PosMatrix <- do.call(rbind,
                       strsplit(x = Taxonomy,
                                split = "; ",
                                fixed = TRUE))
  IsSinglton <- matrix(data = TRUE,
                       ncol = ncol(PosMatrix),
                       nrow = nrow(PosMatrix))
  
  # return(PosMatrix)
  
  if (is.null(EnforceSingleton)) {
    
    # pull singletons naturally from the data
    
    for (m1 in seq_len(ncol(IsSinglton))) {
      
      UniqueRows <- unique(PosMatrix[duplicated(PosMatrix[, 1:m1, drop = FALSE]), 1:m1, drop = FALSE])
      
      # for lists...
      # PresRows <- do.call(paste0,
      #                     PosMatrix[, 1:m1, drop = FALSE]) %in% do.call(paste0,
      #                                                                   UniqueRows)
      PresRows <- apply(X = PosMatrix[, 1:m1, drop = FALSE],
                        MARGIN = 1L,
                        FUN = function(x) paste0(x,
                                                 collapse = "")) %in% apply(X = UniqueRows,
                                                                            MARGIN = 1L,
                                                                            FUN = function(x) paste0(x,
                                                                                                     collapse = ""))
      IsSinglton[, m1] <- !PresRows
      if (Verbose) {
        setTxtProgressBar(pb = pBar,
                          value = m1 / ncol(IsSinglton))
      }
    }
    
    # return(IsSinglton)
    
    if (Verbose) {
      cat("\n")
      print("Singleton Matrix complete.")
      cat("\n")
    }
  } else {
    
    # enforce singleton status for modified holdout sets
    
    for (m1 in seq_len(ncol(IsSinglton))) {
      IsSinglton[, m1] <- EnforceSingleton
    }
    
    if (Verbose) {
      cat("\n!!Singleton status set by user!!\n")
    }
    
  }
  
  
  
  if (is.null(Level)) {
    Level <- ncol(IsSinglton)
  }
  
  # collect confidences
  # collect correct predictions
  
  if (is(ResList, "list")) {
    # print("Here")
    AllConfs <- sapply(ResList,
                       function(x) x[[1]][[2]][Level],
                       simplify = TRUE,
                       USE.NAMES = FALSE)
    AllAssign <- sapply(ResList,
                        function(x) paste(x[[1]][[1]][1:Level],
                                          collapse = "; "),
                        USE.NAMES = FALSE,
                        simplify = TRUE)
  } else if (is(ResList, "Taxa")) {
    # print("There")
    AllConfs <- sapply(ResList,
                       function(x) x[[2]][Level],
                       simplify = TRUE,
                       USE.NAMES = FALSE)
    AllAssign <- sapply(ResList,
                        function(x) paste(x[[1]][1:Level],
                                          collapse = "; "),
                        USE.NAMES = FALSE,
                        simplify = TRUE)
  } else {
    stop("Class Not Compatible.")
  }
  
  
  CheckVal <- apply(X = PosMatrix[, 1:Level, drop = FALSE],
                    MARGIN = 1L,
                    FUN = function(x) paste0(x,
                                             collapse = "; "))
  
  if (is.null(ResPositions)) {
    AssignmentsCorrect <- AllAssign == CheckVal
  } else {
    AssignmentsCorrect <- AllAssign == CheckVal[ResPositions]
  }
  
  
  
  # return(list(AllConfs,
  #             AllAssign,
  #             AssignmentsCorrect,
  #             CheckVal))
  
  # x[c] = (number of non-singletons >= c) / (total number of non-singletons)
  # OCR[c] = (number of singletons >= c) / (total number of singletons)
  # MCR[c] = (number of non-singletons incorrectly classified >= c) / (total number of non-singletons)
  
  ConfScale <- seq(0, 100, 0.1)
  
  ClassRate <- vector(mode = "integer",
                      length = length(ConfScale))
  OverRate <- vector(mode = "integer",
                     length = length(ConfScale))
  MissRate <- vector(mode = "integer",
                     length = length(ConfScale))
  
  # return(list(IsSinglton,
  #             AllConfs,
  #             ConfScale,
  #             AssignmentsCorrect,
  #             AllAssign,
  #             CheckVal))
  
  
  if (is.null(ResPositions)) {
    for (m1 in seq_along(ConfScale)) {
      ClassRate[m1] <- length(which(!IsSinglton[, Level] &
                                      AllConfs >= ConfScale[m1])) / sum(!IsSinglton[, Level])
      OverRate[m1] <- length(which(IsSinglton[, Level] &
                                     AllConfs >= ConfScale[m1])) / sum(IsSinglton[, Level])
      MissRate[m1] <- length(which(!IsSinglton[, Level] &
                                     AllConfs >= ConfScale[m1] &
                                     !AssignmentsCorrect)) / sum(!IsSinglton[, Level])
      
      if (Verbose) {
        setTxtProgressBar(pb = pBar,
                          value = m1 / length(ConfScale))
      }
    }
  } else {
    for (m1 in seq_along(ConfScale)) {
      ClassRate[m1] <- length(which(!IsSinglton[ResPositions, Level] &
                                      AllConfs >= ConfScale[m1])) / sum(!IsSinglton[ResPositions, Level])
      OverRate[m1] <- length(which(IsSinglton[ResPositions, Level] &
                                     AllConfs >= ConfScale[m1])) / sum(IsSinglton[ResPositions, Level])
      MissRate[m1] <- length(which(!IsSinglton[ResPositions, Level] &
                                     AllConfs >= ConfScale[m1] &
                                     !AssignmentsCorrect)) / sum(!IsSinglton[ResPositions, Level])
      
      if (Verbose) {
        setTxtProgressBar(pb = pBar,
                          value = m1 / length(ConfScale))
      }
    }
  }
  
  
  if (Verbose) {
    cat("\n")
    FunEnd <- Sys.time()
    print(FunEnd - FunStart)
    cat("\n")
  }
  
  if (IncludeConf) {
    ConfMat <- matrix(data = NA_real_,
                      nrow = 3L,
                      ncol = length(ConfPoints))
    ConfMat[1L, ] <- ClassRate[ConfScale %in% ConfPoints]
    ConfMat[2L, ] <- OverRate[ConfScale %in% ConfPoints]
    ConfMat[3L, ] <- MissRate[ConfScale %in% ConfPoints]
  }
  
  if (ShowPlot) {
    plot(x = 0,
         y = 0,
         type = "n",
         xlim = c(0, 1),
         ylim = c(0, 1),
         main = "OC/MC",
         xlab = "Classification Rate",
         ylab = "Error Rate")
    lines(x = ClassRate,
          y = OverRate,
          col = "red",
          lty = 2,
          lwd = 1.5)
    lines(x = ClassRate,
          y = MissRate,
          col = "blue",
          lty = 2,
          lwd = 1.5)
    if (IncludeConf) {
      points(x = ConfMat[1, ],
             y = ConfMat[2, ],
             pch = 20,
             col = "red")
      points(x = ConfMat[1, ],
             y = ConfMat[3, ],
             pch = 20,
             col = "blue")
    }
  }
  
  if (IncludeConf) {
    return(list(ClassRate,
                OverRate,
                MissRate,
                ConfMat))
  } else {
    return(list(ClassRate,
                OverRate,
                MissRate))
  }
  
}

###### -- Code Body -----------------------------------------------------------

TimeStart <- Sys.time()

cat("\nScript starting at:\n")
print(TimeStart)

load(file = InFile,
     verbose = TRUE)

L <- length(Holdouts[[TargetTax]])
pBar <- txtProgressBar(style = 1L)
Res <- vector(mode = "list",
              length = L)

cat("\nBeginning CVs.\n")
for (m1 in seq_len(L)) {
  
  # create a placeholder vector for holdout ids
  ph <- Holdouts[[TargetTax]][[m1]]
  
  # learn the taxonomy minus the holdouts
  if (UseStandardAlphabet) {
    train <- LearnTaxa(train = Seqs[-ph],
                       taxonomy = TaxVectors[[TargetTax]][-ph],
                       maxChildren = 1, # for the provided data
                       alphabet = AA_STANDARD,
                       verbose = FALSE)
    
    Res[[m1]] <- IdTaxa(test = Seqs[ph],
                        trainingSet = train,
                        fullLength = 0.95,
                        threshold = 0,
                        verbose = FALSE)
  } else {
    train <- LearnTaxa(train = Seqs[-ph],
                       taxonomy = TaxVectors[[TargetTax]][-ph],
                       maxChildren = 1, # for the provided data
                       verbose = FALSE)
    
    Res[[m1]] <- IdTaxa(test = Seqs[ph],
                        trainingSet = train,
                        fullLength = 0.95,
                        threshold = 0,
                        verbose = FALSE)
  }
  
  
  setTxtProgressBar(pb = pBar,
                    value = m1 / L)
}

TimeInt <- Sys.time()

cat("\n")

cat("\nCVs completed in:\n")
print(TimeInt - TimeStart)

IRes <- do.call(c,
                Res)

IPerf <- OC_MC_FromIDTaxa(ResList = IRes,
                          Taxonomy = TaxVectors[[TargetTax]],
                          ResPositions = unlist(Holdouts[[TargetTax]]),
                          IncludeConf = TRUE,
                          Verbose = TRUE)

TimeStop <- Sys.time()

###### -- Save Results --------------------------------------------------------

save(IRes, # IDTaxa result
     IPerf, # IDTaxa performance
     file = OutFile,
     compress = "xz")

cat("\nScript completed in:\n")
print(TimeStop - TimeStart)








