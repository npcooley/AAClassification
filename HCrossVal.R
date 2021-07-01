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

Args <- c("~/Box Sync/Cooley/20200923_AATraining.RData",
          "5",
          "~/Box Sync/AAClassification/Figure2/HMMResTaxLevel.RData",
          "10",
          "TRUE")

Args <- commandArgs(trailingOnly = TRUE)
InFile <- Args[1L]
TargetTax <- as.integer(Args[2L])
OutFile <- Args[3L]
if (length(Args) < 4L) {
  CoresToUse <- 0L
} else {
  CoresToUse <- as.integer(Args[4L])
}

if (length(Args) < 5L) {
  SaveModels <- FALSE
} else {
  SaveModels <- as.logical(Args[5L])
}

t1 <- try(system(command = "hmmscan -h", intern = TRUE))
if (is(object = t1,
       class2 = "try-error")) {
  stop ("R appears to not have access to HMMER.")
}

# give HMMER access to cores
if (CoresToUse == 0) {
  AllowCores <- detectCores()
  if (is.na(AllowCores)) {
    AllowCores <- 1L
  }
} else {
  AllowCores <- CoresToUse
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

# create two vectors that are the length of all labels present in the taxonomy
# for each label create a model for all sequences
# and create a model for any holdout given
# this relies on each label only having 1 holdout
# which should be a constraint that is clearly stated in the manuscript
cat("\nGenerating Models.\n")

CurrentTax <- TaxVectors[[TargetTax]]
PlaceHolderSize <- length(gregexpr(pattern = ";",
                                   text = TaxVectors[[TargetTax]][1L],
                                   fixed = TRUE)[[1]])
UniqueTax <- unique(CurrentTax)
HMMTax <- gsub(pattern = " ",
               replacement = "_",
               x = UniqueTax)
AllHoldouts <- unlist(Holdouts[[TargetTax]])
L <- length(UniqueTax)
NonHoldoutHMMs <- HoldoutHMMs <- vector(mode = "list",
                                        length = L)
pBar <- txtProgressBar(style = 1)

for (m1 in seq_len(L)) {
  
  TEMP01 <- tempfile()
  TEMP02 <- tempfile()
  TEMP03 <- tempfile()
  
  w <- which(CurrentTax == UniqueTax[m1])
  CurrentSeqs <- Seqs[w]
  if (length(CurrentSeqs) == 1L) {
    # singleton, build model, leave holdout null
    
    writeXStringSet(x = CurrentSeqs,
                    filepath = TEMP01)
    BuildHMM <- paste("hmmbuild ",
                      "-n ",
                      '\"',
                      HMMTax[m1],
                      '\"',
                      " --cpu ",
                      AllowCores,
                      " ",
                      TEMP02,
                      " ",
                      TEMP01,
                      sep = "")
    system(command = BuildHMM,
           intern = TRUE)
    NonHoldoutHMMs[[m1]] <- readLines(TEMP02)
    
  } else {
    # non-singleton build model
    # if a holdout is designated by any of the holdout vectors
    # build that model as well
    
    # first ask if a holdout is present in the provided labels
    w1 <- any(w %in% AllHoldouts)
    if (w1 & length(w) == 2L) {
      # holdout present, when withheld only one remains, no alignment necessary
      # for non-holdout
      Ali <- AlignSeqs(CurrentSeqs,
                       processors = AllowCores,
                       verbose = FALSE)
      writeXStringSet(x = Ali,
                      filepath = TEMP01)
      BuildHMM <- paste("hmmbuild ",
                        "-n ",
                        '\"',
                        HMMTax[m1],
                        '\"',
                        " --cpu ",
                        AllowCores,
                        " ",
                        TEMP02,
                        " ",
                        TEMP01,
                        sep = "")
      system(command = BuildHMM,
             intern = TRUE)
      NonHoldoutHMMs[[m1]] <- readLines(TEMP02)
      
      # for holdout
      w2 <- w[!(w %in% AllHoldouts)]
      HoldoutSeqs <- Seqs[w2]
      writeXStringSet(x = HoldoutSeqs,
                      filepath = TEMP03)
      BuildHMM <- paste("hmmbuild ",
                        "-n ",
                        '\"',
                        HMMTax[m1],
                        '\"',
                        " --cpu ",
                        AllowCores,
                        " ",
                        TEMP03,
                        " ",
                        TEMP01,
                        sep = "")
      system(command = BuildHMM,
             intern = TRUE)
      HoldoutHMMs[[m1]] <- readLines(TEMP03)
      
    } else if (w1 & length(w) > 2L) {
      # holdout present, when withheld, more than one seq remains,
      # alignment necessary
      Ali <- AlignSeqs(CurrentSeqs,
                       processors = AllowCores,
                       verbose = FALSE)
      writeXStringSet(x = Ali,
                      filepath = TEMP01)
      BuildHMM <- paste("hmmbuild ",
                        "-n ",
                        '\"',
                        HMMTax[m1],
                        '\"',
                        " --cpu ",
                        AllowCores,
                        " ",
                        TEMP02,
                        " ",
                        TEMP01,
                        sep = "")
      system(command = BuildHMM,
             intern = TRUE)
      NonHoldoutHMMs[[m1]] <- readLines(TEMP02)
      
      # for holdout
      w2 <- w[!(w %in% AllHoldouts)]
      HoldoutSeqs <- Seqs[w2]
      Ali <- AlignSeqs(HoldoutSeqs,
                       processors = AllowCores,
                       verbose = FALSE)
      writeXStringSet(x = Ali,
                      filepath = TEMP03)
      BuildHMM <- paste("hmmbuild ",
                        "-n ",
                        '\"',
                        HMMTax[m1],
                        '\"',
                        " --cpu ",
                        AllowCores,
                        " ",
                        TEMP03,
                        " ",
                        TEMP01,
                        sep = "")
      system(command = BuildHMM,
             intern = TRUE)
      HoldoutHMMs[[m1]] <- readLines(TEMP03)
      
    } else {
      # no holdout present, and singletons already accounted for,
      # align seqs designated by w
      Ali <- AlignSeqs(CurrentSeqs,
                       processors = AllowCores,
                       verbose = FALSE)
      writeXStringSet(x = Ali,
                      filepath = TEMP01)
      BuildHMM <- paste("hmmbuild ",
                        "-n ",
                        '\"',
                        HMMTax[m1],
                        '\"',
                        " --cpu ",
                        AllowCores,
                        " ",
                        TEMP02,
                        " ",
                        TEMP01,
                        sep = "")
      system(command = BuildHMM,
             intern = TRUE)
      NonHoldoutHMMs[[m1]] <- readLines(TEMP02)
    }
  }
  
  unlink(TEMP01)
  unlink(TEMP02)
  unlink(TEMP03)
  
  setTxtProgressBar(pb = pBar,
                    value = m1 / L)
}
cat("\nModel construction completed in:\n")
TimeInt <- Sys.time()
print(TimeInt - TimeStart)

save(HoldoutHMMs,
     NonHoldoutHMMs,
     file = OutFile,
     compress = "xz")

L <- length(Holdouts[[TargetTax]])
pBar <- txtProgressBar(style = 1L)
Res <- HRes <- vector(mode = "list",
                      length = L)

cat("\nBeginning CVs.\n")
for (m1 in seq_len(L)) {
  
  TEMP01 <- tempfile()
  TEMP02 <- tempfile()
  TEMP03 <- tempfile()
  
  # two large lists of character vectors exist that the HMM libraries for
  # cross validation can be constructed from
  # for each CV a library will be constructed and written out to a temp file
  # HMMER will press and then scan that time file
  # HMMER will need the testing sequences as well, those are written out to another
  # temp file
  
  # writing test sequences out to a tempfile:
  writeXStringSet(x = Seqs[Holdouts[[TargetTax]][[m1]]],
                  filepath = TEMP01,
                  append = FALSE)
  CurrentSeqIDs <- names(Seqs[Holdouts[[TargetTax]][[m1]]])
  # determine which unique taxons are being pulled from the holdout
  # and which are being pulled from the non-holdout
  w <- which(UniqueTax %in% TaxVectors[[TargetTax]][Holdouts[[TargetTax]][[m1]]])
  
  DummyList <- vector(mode = "list",
                      length = length(UniqueTax))
  for (m2 in seq_along(UniqueTax)) {
    if (m2 %in% w) {
      DummyList[[m2]] <- HoldoutHMMs[[m2]]
    } else {
      DummyList[[m2]] <- NonHoldoutHMMs[[m2]]
    }
  }
  
  writeLines(text = unlist(DummyList),
             con = TEMP02)
  
  system(command = paste("hmmpress -f ",
                         TEMP02,
                         sep = ""),
         intern = FALSE)
  
  system(command = paste("hmmscan ",
                         "--noali ",
                         "--tblout ",
                         TEMP03,
                         " --cpu ",
                         AllowCores,
                         " ",
                         TEMP02,
                         " ",
                         TEMP01,
                         sep = ""),
         intern = FALSE)
  
  Res[[m1]] <- readLines(TEMP03)
  
  x <- strsplit(x = Res[[m1]],
                split = "[ ]+",
                fixed = FALSE)
  x <- x[lengths(x) == 19L]
  x <- x[-1L]
  x <- do.call(rbind,
               x)
  colnames(x) <- c("Target_name",
                   "accession",
                   "query_name",
                   "accession",
                   "full_seq_Eval",
                   "full_seq_score",
                   "full_seq_bias",
                   "best_1_Eval",
                   "best_1_score",
                   "best_1_bias",
                   "exp",
                   "reg",
                   "clu",
                   "ov",
                   "env",
                   "dom",
                   "rep",
                   "inc",
                   "description")
  Res[[m1]] <- x
  x <- data.frame("Target_name" = x[, 1L],
                  "Target_ACC" = x[, 2L],
                  "Query_name" = x[, 3L],
                  "Query_ACC" = x[, 4L],
                  "Full_Seq_EVal" = as.numeric(x[, 5L]),
                  "Full_Seq_Score" = as.numeric(x[, 6L]),
                  "Full_Seq_Bias" = as.numeric(x[, 7L]))
  
  # put an IDTaxa like object into HRes[[m1]]
  # capture minimum HMMER evalue
  # set any non-found seqs to an evalue of 1
  # set any evalues greater than 1 to 1
  minEvalue <- min(x$Full_Seq_EVal[x$Full_Seq_EVal > 0])
  HEval <- vector(mode = "numeric",
                  length = length(CurrentSeqIDs))
  HTarget <- HID <- vector(mode = "character",
                           length = length(CurrentSeqIDs))
  for (m2 in seq_along(CurrentSeqIDs)) {
    ph <- x[x$Query_name == CurrentSeqIDs[m2], , drop = FALSE]
    if (nrow(ph) == 0L) {
      # HMMER found no significant hits, assign placeholder at zero confidence
      HEval[m2] <- 1
      HTarget[m2] <- paste("Root",
                           paste(rep("Placeholder",
                                     PlaceHolderSize),
                                 collapse = ";_"),
                           sep = ";_")
      HID[m2] <- CurrentSeqIDs[m2]
    } else {
      # HMMER found hits, capture most significant and it's 'Target_name'
      ph1 <- which.min(ph$Full_Seq_EVal)
      HEval[m2] <- ph$Full_Seq_EVal[ph1]
      HTarget[m2] <- ph$Target_name[ph1]
      HID[m2] <- ph$Query_name[ph1]
    }
  }
  
  if (any(HEval > 1)) {
    HEval[HEval > 1] <- 1
  }
  
  if (any(HEval == 0)) {
    HEval[HEval == 0] <- minEvalue
  }
  
  HConf <- (-log10(HEval) / -log10(minEvalue)) * 100
  HAssign <- gsub(pattern = "_",
                  replacement = " ",
                  x = HTarget)
  
  HRes[[m1]] <- vector(mode = "list",
                       length = length(CurrentSeqIDs))
  for (m2 in seq_along(HRes[[m1]])) {
    HRes[[m1]][[m2]] <- list("taxon" = strsplit(HAssign[m2],
                                                split = "; ",
                                                fixed = TRUE)[[1]],
                             "confidence" = rep(HConf[m2],
                                                (PlaceHolderSize + 1L)))
  }
  
  names(HRes[[m1]]) <- HID
  class(HRes[[m1]]) <- c("Taxa", "Test")
  
  # b <- OC_MC_FromIDTaxa(ResList = HRes[[m1]],
  #                       ResPositions = Holdouts[[TargetTax]][[m1]],
  #                       Taxonomy = CurrentTax,
  #                       ShowPlot = TRUE,
  #                       Verbose = TRUE)
  
  
  unlink(TEMP01)
  unlink(TEMP02)
  unlink(TEMP03)
  
  # HMMER prints out statements during this loop
  # setTxtProgressBar(pb = pBar,
  #                   value = m1 / L)
}

TimeInt <- Sys.time()

cat("\n")

cat("\nCVs completed in:\n")
print(TimeInt - TimeStart)

Res <- do.call(rbind,
               Res)

HRes <- do.call(c,
                HRes)

HPerf <- OC_MC_FromIDTaxa(ResList = HRes,
                          Taxonomy = TaxVectors[[TargetTax]],
                          ResPositions = unlist(Holdouts[[TargetTax]]),
                          IncludeConf = TRUE,
                          Verbose = TRUE)

TimeStop <- Sys.time()

###### -- Save Results --------------------------------------------------------

if (SaveModels) {
  save(HoldoutHMMs,
       NonHoldoutHMMs,
       Res,
       HRes,
       HPerf,
       file = OutFile,
       compress = "xz")
} else {
  save(Res,
       HRes,
       HPerf,
       file = OutFile,
       compress = "xz")
}


cat("\nScript completed in:\n")
print(TimeStop - TimeStart)








