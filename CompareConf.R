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

pdf(paste0(TaxCode,
           "_Comparison.pdf"),
    height = 3.5,
    width = 7)

layout(mat = matrix(data = 1:2,
                    ncol = 2))
par(mgp = c(2.2, 1.0, 0.0))
par(mar = c(3.5, 3.5, 0.5, 1.0))
plot(x = x,
     y = y,
     pch = 20,
     cex = 0.5,
     col = ifelse(test = xt == yt,
                  yes = ColVector[2L],
                  no = ColVector[1L]),
     xlab = "IDTAXA Confidence",
     ylab = "BLAST PID")
legend("bottomright",
       legend = c("Classifications agree",
                  "Classifications disagree"),
       border = NA,
       pch = 20,
       col = c("blue",
               "red"),
       cex = 0.75)

par(mgp = c(2.2, 1.0, 0.0))
par(mar = c(3.5, 3.5, 0.5, 1.0))
plot(x = x,
     y = z,
     pch = 20,
     cex = 0.5,
     col = ifelse(test = xt == zt,
                  yes = ColVector[2L],
                  no = ColVector[1L]),
     xlab = "IDTAXA Confidence",
     ylab = "BLAST EVL as Confidence")
legend("bottomright",
       legend = c("Classifications agree",
                  "Classifications disagree"),
       border = NA,
       pch = 20,
       col = c("blue",
               "red"),
       cex = 0.75)
dev.off()



