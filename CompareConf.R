###### -- ask how BLAST and IDTAXA differ on annotations ----------------------
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

###### -- Libraries -----------------------------------------------------------

suppressMessages(library(DECIPHER))

###### -- Arguments -----------------------------------------------------------

Args <- commandArgs(trailingOnly = TRUE)
# Args <- c("~/Box Sync/Cooley/20200923_AATraining.RData",
#           "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/011/074/885/GCF_011074885.1_ASM1107488v2/GCF_011074885.1_ASM1107488v2_protein.faa.gz",
#           "bbrx")
TrainingData <- Args[1L]
ProtFASTA <- Args[2L]
TaxCode <- Args[3L]

###### -- Functions -----------------------------------------------------------

GenericBLAST <- function(Seqs,
                         DB,
                         UseCores = FALSE) {
  
  SubjectSeqs <- DB
  QuerySeqs <- Seqs
  
  if (UseCores) {
    UseCores <- detectCores()
    if (is.na(UseCores)) {
      UseCores <- 1L
    }
  } else {
    UseCores <- 1L
  }
  
  dir.create(paste("./TempBLAST"))
  writeXStringSet(x = SubjectSeqs,
                  filepath = "./TempBLAST/SubjectSeqs.fa",
                  append = FALSE)
  writeXStringSet(x = QuerySeqs,
                  filepath = "./TempBLAST/QuerySeqs.fa",
                  append = FALSE)
  
  if (is(object = Seqs,
         class2 = "AAStringSet")) {
    system(command = paste("makeblastdb -in ",
                           "./TempBLAST/SubjectSeqs.fa ",
                           "-dbtype prot ",
                           "-out ./TempBLAST/DB",
                           sep = ""))
    system(command = paste("blastp -query ",
                           "./TempBLAST/QuerySeqs.fa",
                           " -outfmt 6 -evalue ",
                           1e-3,
                           " -num_threads ",
                           UseCores,
                           " -use_sw_tback -db ",
                           "./TempBLAST/DB",
                           " -out ./TempBLAST/OutTable.txt",
                           sep = ""))
  } else if (is(object = Seqs,
                class2 = "DNAStringSet")) {
    system(command = paste("makeblastdb -in ",
                           "./TempBLAST/SubjectSeqs.fa ",
                           "-dbtype nucl ",
                           "-out ./TempBLAST/DB",
                           sep = ""))
    system(command = paste("blastn -query ",
                           "./TempBLAST/QuerySeqs.fa",
                           " -outfmt 6 -evalue ",
                           1e-3,
                           " -num_threads ",
                           UseCores,
                           " -db ",
                           "./TempBLAST/DB",
                           " -out ./TempBLAST/OutTable.txt",
                           sep = ""))
  }
  
  
  if (length(readLines("./TempBLAST/OutTable.txt")) > 0L) {
    
    BLASTResults <- read.table("./TempBLAST/OutTable.txt",
                               sep = "\t",
                               header = FALSE,
                               stringsAsFactors = FALSE)
    
    colnames(BLASTResults) <- c("qseqid",
                                "sseqid",
                                "pident",
                                "length",
                                "mismatch",
                                "gapopen",
                                "qstart",
                                "qend",
                                "sstart",
                                "send",
                                "evalue",
                                "bitscore")
    
    unlink("./TempBLAST/OutTable.txt")
    system("rm -rf ./TempBLAST")
    return(BLASTResults)
    
  } else {
    
    system("rm -rf ./TempBLAST")
    return(NULL)
  }
  
}

###### -- Code Body -----------------------------------------------------------

load(file = TrainingData,
     verbose = TRUE)

TargetSeqs <- readAAStringSet(ProtFASTA)
FullSeqNames <- names(TargetSeqs)
SeqIDs <- unname(sapply(FullSeqNames,
                        function(x) strsplit(x = x,
                                             split = " ",
                                             fixed = TRUE)[[1]][1],
                        USE.NAMES = FALSE,
                        simplify = TRUE))
names(TargetSeqs) <- SeqIDs

SeqTax <- unlist(regmatches(x = names(Seqs),
                            m = gregexpr(pattern = "(?<=\\.)([^:]+)(?=:)",
                                         text = names(Seqs),
                                         perl = TRUE)))

Seqs <- Seqs[!(SeqTax %in% TaxCode)]
Tax <- TaxVectors[[1L]][!(SeqTax %in% TaxCode)]

train <- LearnTaxa(train = Seqs,
                   taxonomy = Tax,
                   maxChildren = 1L,
                   verbose = TRUE)

IRes <- IdTaxa(test = TargetSeqs,
               trainingSet = train,
               fullLength = 0.95,
               threshold = 0,
               verbose = TRUE)

Res <- GenericBLAST(Seqs = TargetSeqs,
                    DB = Seqs,
                    UseCores = TRUE)

EVL <- PID <- vector(mode = "numeric",
                     length = length(TargetSeqs))
EVLClass <- PIDClass <- EVLID <- PIDID <- vector(mode = "character",
                                                 length = length(TargetSeqs))

pBar <- txtProgressBar(style = 1L)
L <- length(TargetSeqs)
PlaceHolderSize <- length(gregexpr(pattern = ";",
                                   text = Tax[1L],
                                   fixed = TRUE)[[1]])
ph1 <- names(TargetSeqs)
ph3 <- min(Res$evalue[Res$evalue > 0])
for (m2 in seq_along(TargetSeqs)) {
  if (ph1[m2] %in% Res$qseqid) {
    # id hit in table
    ph2 <- Res[Res$qseqid == ph1[m2], , drop = FALSE]
    w1 <- which.max(ph2$pident)
    PID[m2] <- ph2$pident[w1]
    PIDClass[m2] <- Tax[names(Seqs) == ph2$sseqid[w1]]
    PIDID[m2] <- ph2$sseqid[w1]
    w1 <- which.min(ph2$evalue)
    if (ph2$evalue[w1] == 0L) {
      EVL[m2] <- ph3
    } else {
      EVL[m2] <- ph2$evalue[w1]
    }
    EVLClass[m2] <- Tax[names(Seqs) == ph2$sseqid[w1]]
    EVLID[m2] <- ph2$sseqid[w1]
  } else {
    # id hit not in table
    EVL[m2] <- 1
    PID[m2] <- 0
    EVLClass[m2] <- PIDClass[m2] <- paste("Root",
                                          paste(rep("placeholder_classification",
                                                    PlaceHolderSize),
                                                collapse = "; "),
                                          sep = "; ")
    EVLID[m2] <- PIDID[m2] <- paste("placeholder_",
                                    m2,
                                    sep = "")
  }
  
  setTxtProgressBar(pb = pBar,
                    value = m2 / length(ph1))
}

EVL <- (-log10(EVL) / max(-log10(EVL))) * 100
PRes <- ERes <- vector(mode = "list",
                       length = length(ph1))
for (m2 in seq_along(PRes)) {
  PRes[[m2]] <- list("taxon" = strsplit(PIDClass[m2],
                                              split = "; ",
                                              fixed = TRUE)[[1]],
                           "confidence" = rep(PID[m2],
                                              (PlaceHolderSize + 1L)))
  ERes[[m2]] <- list("taxon" = strsplit(EVLClass[m2],
                                              split = "; ",
                                              fixed = TRUE)[[1]],
                           "confidence" = rep(EVL[m2],
                                              (PlaceHolderSize + 1L)))
}

names(PRes) <- PIDID
names(ERes) <- EVLID

class(PRes) <- c("Taxa", "Test")
class(ERes) <- c("Taxa", "Test")

save(PRes,
     ERes,
     IRes,
     file = paste0(TaxCode,
                   "_Result.RData"),
     compress = "xz")

x <- sapply(IRes,
            function(x) x$confidence[2],
            USE.NAMES = FALSE,
            simplify = TRUE)
y <- sapply(PRes,
            function(x) x$confidence[2],
            USE.NAMES = FALSE,
            simplify = TRUE)
z <- sapply(ERes,
            function(x) x$confidence[2],
            USE.NAMES = FALSE,
            simplify = TRUE)
xt <- sapply(IRes,
             function(x) x$taxon[2],
             USE.NAMES = FALSE,
             simplify = TRUE)
yt <- sapply(PRes,
             function(x) x$taxon[2],
             USE.NAMES = FALSE,
             simplify = TRUE)
zt <- sapply(ERes,
             function(x) x$taxon[2],
             USE.NAMES = FALSE,
             simplify = TRUE)

ColVector <- c(rgb(red = 0.95,
                   blue = 0.05,
                   green = 0.05,
                   alpha = 0.25),
               rgb(red = 0.05,
                   blue = 0.95,
                   green = 0.05,
                   alpha = 0.25))

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
               '#000000')[c(4,2)]

AVal <- "40"
ColVector <- paste0(ColVector,
                    AVal)


pdf(paste0(TaxCode,
           "_PIDComparison.pdf"),
    height = 7,
    width = 7)
par(mgp = c(2.2, 1.0, 0.0))
par(mar=c(10, 10, 1, 1))
plot(x = x,
     y = y,
     pch = 16,
     cex = 0.75,
     col = ifelse(test = xt == yt,
                  yes = ColVector[2L],
                  no = ColVector[1L]),
     xlab = "IDTAXA Confidence (%)",
     ylab = "BLAST PID (%)")
legend("bottomright",
       legend = c("Classifications disagree",
                  "Classifications agree"),
       border = NA,
       pch = 16,
       col = gsub(x = ColVector,
                  pattern = "40$",
                  replacement = ""),
       cex = 0.75)
dev.off()

pdf(paste0(TaxCode,
           "_EVLComparison.pdf"),
    height = 7,
    width = 7)
par(mgp = c(2.2, 1.0, 0.0))
par(mar=c(10, 10, 1, 1))

AxisPosition <- c(100,
                  1e-18,
                  1e-38,
                  1e-58,
                  1e-78,
                  1e-99,
                  1e-119,
                  1e-139,
                  1e-159,
                  1e-180)
AxisLabels <- c("Unassigned",
                as.character(c(1e-18,
                               1e-38,
                               1e-58,
                               1e-78,
                               1e-99,
                               1e-119,
                               1e-139,
                               1e-159,
                               1e-180)))
MinEVL <- 1e-180
yy <- 10^(log10(MinEVL) * (z * 0.01))
yy[yy == 1] <- 100
plot(x = x,
     y = yy,
     pch = 16,
     cex = 0.5,
     log = "y",
     ylim = c(1, 1e-180),
     col = ifelse(test = xt == zt,
                  yes = ColVector[2L],
                  no = ColVector[1L]),
     xlab = "IDTAXA Confidence (%)",
     ylab = "BLAST EVL",
     axes = FALSE,
     frame.plot = TRUE)
axis(side = 2,
     at = AxisPosition,
     labels = AxisLabels)
axis.break(axis = 2,
           breakpos = 1e-3)
axis(side = 1,
     at = seq(from = 0,
              to = 100,
              by = 20),
     labels = seq(from = 0,
                  to = 100,
                  by = 20))
legend("bottomright",
       legend = c("Classifications disagree",
                  "Classifications agree"),
       border = NA,
       pch = 16,
       col = gsub(x = ColVector,
                  pattern = "40$",
                  replacement = ""),
       cex = 0.75)
dev.off()



