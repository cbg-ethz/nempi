
library(snowfall)
library(TCGAbiolinks)

if (is.null(types)) {
    types <- TCGAbiolinks:::getGDCprojects()$project_id
    ## special "FM-AD"
    dont <- paste0(unique(gsub("-.*", "", types)), "-")
    dont <- dont[-grep("TCGA", dont)]
    samplenr <- matrix(NA, length(types), 2)
    rownames(samplenr) <- types
    donr <- TRUE
    snrcount <- 0
} else {
    dont <- "AGCT"
    donr <- FALSE
}

sizemat <- matrix(0, 1, 2)
colnames(sizemat) <- c("Tumor", "Normal")
rownames(sizemat) <- ""

path <- "mutclust/"

for (type in types) {
    if (donr) {
        snrcount <- snrcount + 1
    }
    if (length(grep(paste(dont, collapse = "|"), type)) > 0) { next() }
    print(type)
    if (file.exists(paste0(path, type, "_final.rda")) & !newmut & !newllr & !newsave) {
        load(paste0(path, type, "_final.rda"))
    } else {

        summ <- TCGAbiolinks:::getProjectSummary(type)
        library(SummarizedExperiment)

        ## get methylation:
        if (file.exists(paste0(path, type, "_meth.rda"))) {
            load(paste0(path, type, "_meth.rda"))
        } else {
            data <- GDCquery(project = paste(type, sep = ""),
                             sample.type = "Primary Tumor",
                             data.category = "DNA methylation",
                             data.type = "Methylation beta value",
                             legacy = TRUE
                             )
            GDCdownload(data)
            data <- GDCprepare(data, summarizedExperiment = 0)
            ## map methy sites to genes:
            library(methyAnalysis)
            if (is.data.frame(data)) {
                meth2 <- as.matrix(data[, -(1:3)])
                rownames(meth2) <- gsub(";.*", "", data[, 1])
            } else {
                meth <- data@rowRanges
                library(TxDb.Hsapiens.UCSC.hg19.knownGene)
                probe2gene <- annotateDMRInfo(meth, 'TxDb.Hsapiens.UCSC.hg19.knownGene')
                meth2 <- assay(data)
                rownames(meth2) <- probe2gene$sigDMRInfo@elementMetadata@listData$GeneSymbol
            }
            meth <- meth2
            meth <- meth[which(apply(meth, 1, function(x) return(any(is.na(x)))) == FALSE), ]
            methm <- meth[which(duplicated(rownames(meth)) == FALSE), ]
            count <- 0
            for (i in which(duplicated(rownames(meth)) == FALSE)) {
                count <- count + 1
                methm[count, ] <- apply(meth[which(rownames(meth) %in% rownames(methm)[i]), , drop = FALSE], 2, median)
            }
            meth <- methm
            meth <- meth[order(rownames(meth)), order(colnames(meth))]
            colnames(meth) <- gsub("A$", "", lapply(strsplit(colnames(meth), "-"), function(x) { y <- x[1:4]; y <- paste(y, collapse = "-"); return(y) }))
            methm <- meth[, which(duplicated(colnames(meth)) == FALSE)]
            for (i in which(duplicated(colnames(meth)) == TRUE)) {
                j <- which(colnames(methm) == colnames(meth)[i])
                methm[, j] <- apply(meth[, which(colnames(meth) %in% colnames(methm)[i]), drop = FALSE], 2, median)
            }
            meth <- methm
            meth <- meth[order(rownames(meth)), order(colnames(meth))]
            save(meth, file = paste0(path, type, "_meth.rda"))
        }

        print("meth done")

        ## get copy number variation:
        if (file.exists(paste0(path, type, "_cnv.rda"))) {
            load(paste0(path, type, "_cnv.rda"))
        } else {
            data <- getGistic(gsub("TCGA-", "", type), type = "thresholded")
            ## data <- GDCquery(project = paste(type, sep = ""),
            ##                  sample.type = "Primary solid Tumor",
            ##                  data.category = "Copy Number Variation",
            ##                  data.type = "Copy Number Segment",
            ##                  )
            ## GDCdownload(data)
            ## data <- GDCprepare(data)
            cnv <- data[, -(1:3)]
            cnv <- apply(cnv, c(1,2), as.numeric)
            rownames(cnv) <- data[, 1]
            colnames(cnv) <- gsub("A$", "", lapply(strsplit(colnames(cnv), "-"), function(x) { y <- x[1:4]; y <- paste(y, collapse = "-"); return(y) }))
            cnv <- cnv[order(rownames(cnv)), order(colnames(cnv))]
            save(cnv, file = paste0(path, type, "_cnv.rda"))
        }

        print("cnv done")

        ## get expression
        if (file.exists(paste0(path, type, "_query.rda"))) {
            load(paste0(path, type, "_query.rda"))
        } else {
            if (length(grep("-AML$|-LAML$", type)) > 0) {
                sampletype <- c("Primary Blood Derived Cancer - Peripheral Blood", "Solid Tissue Normal")
            } else {
                sampletype <- c("Primary solid Tumor", "Solid Tissue Normal")
            }
            data <- GDCquery(project = paste(type, sep = ""),
                             sample.type = sampletype,
                             data.category = "Transcriptome Profiling",
                             data.type = "Gene Expression Quantification",
                             workflow.type = "HTSeq - Counts"
                             )
            GDCdownload(data)
            if (is.null(data)) {
                print("is null")
                next()
            }
            data <- GDCprepare(data)
            save(data, file = paste0(path, type, "_query.rda"))
        }
        ## process expression data:
        D <- assay(data)
        class <- data@colData@listData$definition

        print("gene expression done")

        ## get mutation
        if (file.exists(paste0(path, type, "_mut.rda"))) {
            load(paste0(path, type, "_mut.rda"))
        } else {
            type2 <- gsub(paste(paste0(unique(gsub("-.*", "", types)), "-"), collapse = "|"), "", type)
            library(data.table)
            GDCquery_Maf(tumor = type2, save.csv = TRUE, pipeline = "varscan2") # mutation file
            GDCquery_Maf(tumor = type2, save.csv = TRUE, pipeline = "muse") # mutation file
            GDCquery_Maf(tumor = type2, save.csv = TRUE, pipeline = "somaticsniper") # mutation file
            GDCquery_Maf(tumor = type2, save.csv = TRUE, pipeline = "mutect2") # mutation file
            mut <- list()
            count <- 1
            type3 <- gsub("-", "\\.", type)
            for (i in list.files("GDCdata")) {
                if (length(grep(type3, i)) > 0 & length(grep("csv", i)) > 0) {
                    mut[[count]] <- fread(paste0("GDCdata/", i))
                    count <- count + 1
                }
            }
            save(mut, file = paste0(path, type, "_mut.rda"))
        }

        ## clinical
        if (file.exists(paste0(path, type, "_clin.rda"))) {
            load(paste0(path, type, "_clin.rda"))
        } else {
            clinical <- GDCquery_clinic(project = type, type = "clinical")
            save(clinical, file = paste0(path, type, "_clin.rda"))
        }

        ## process mutation data: (https://cancer.sanger.ac.uk/cosmic/help/mutation/overview)
        if (file.exists(paste0(path, type, "_mut0.rda")) & !newmut) {
            load(paste0(path, type, "_mut0.rda"))
        } else {
            n <- nmut # try something to get all patients with at least one mutation
            library(snowfall)
            countSamples <- function(i, mut.org, genes, mut.mat, coln) {
                i <- which(rownames(mut.mat) %in% genes[i])
                tmp <- mut.org[which(mut.org$Hugo_Symbol %in% rownames(mut.mat)[i]), coln]
                tmp2 <- mut.mat[i, ]
                tmp2[which(colnames(mut.mat) %in% tmp)] <- 1
                return(tmp2)
            }
            typeSamples <- function(i, mut.org, genes, type.mat, coln, coln2) {
                i <- which(rownames(type.mat) %in% genes[i])
                tmp <- mut.org[which(mut.org$Hugo_Symbol %in% rownames(type.mat)[i]), coln]
                tmp3 <- mut.org[which(mut.org$Hugo_Symbol %in% rownames(type.mat)[i]), coln2]
                tmp2 <- type.mat[i, ]
                tmp2[which(colnames(type.mat) %in% tmp)] <- tmp3
                return(tmp2)
            }
            biggl <- list()
            for (i in length(mut)) {
                mutation <- mut[[i]]
                biggl[[i]] <- mutation$Hugo_Symbol
            }
            freq <- sort(table(unlist(biggl)), decreasing = TRUE)
            if (n == 0) {
                allsub <- names(freq)
            } else {
                allsub <- names(freq)[1:n]
            }
            M <- Mtype <- list()
            for (i in 1:length(mut)) {
                mutation <- mut[[i]]
                mut.mat <- matrix(0, length(allsub), length(unique(mutation$Tumor_Sample_Barcode)))
                type.mat <- matrix("", length(allsub), length(unique(mutation$Tumor_Sample_Barcode)))
                colnames(type.mat) <- colnames(mut.mat) <- sort(unique(mutation$Tumor_Sample_Barcode))
                rownames(type.mat) <- rownames(mut.mat) <- allsub
                coln <- which(colnames(mutation) %in% "Tumor_Sample_Barcode")
                coln2 <- which(colnames(mutation) %in% "Variant_Classification")
                mut.org <- mutation[which(mutation$Hugo_Symbol %in% allsub), ]
                sfInit(parallel = T, cpus = 4)
                sfExport("mut.mat", "coln")
                tmp <- sfLapply(as.list(1:length(allsub)), countSamples, mut.org, allsub, mut.mat, coln)
                sfStop()
                tmp <- do.call("rbind", tmp)
                rownames(tmp) <- allsub
                colnames(tmp) <- colnames(mut.mat)
                M[[i]] <- tmp
                sfInit(parallel = T, cpus = 4)
                sfExport("type.mat", "coln2")
                tmp <- sfLapply(as.list(1:length(allsub)), typeSamples, mut.org, allsub, type.mat, coln, coln2)
                sfStop()
                tmp <- do.call("rbind", tmp)
                rownames(tmp) <- allsub
                colnames(tmp) <- colnames(mut.mat)
                Mtype[[i]] <- tmp
            }
            samples <- intersect(intersect(intersect(colnames(M[[1]]), colnames(M[[2]])), colnames(M[[3]])), colnames(M[[4]]))
            M0 <- M[[1]][, which(colnames(M[[1]]) %in% samples)]
            Mtype0 <- Mtype[[1]][, which(colnames(Mtype[[1]]) %in% samples)]
            for (i in 2:length(M)) {
                M0 <- M0 + M[[i]][, which(colnames(M[[i]]) %in% samples)]
                Mtype0 <- matrix(paste(Mtype0, Mtype[[i]][, which(colnames(Mtype[[i]]) %in% samples)]), nrow(Mtype0))
            }
            rownames(Mtype0) <- rownames(M0)
            colnames(Mtype0) <- colnames(M0)
            save(M0, Mtype0, file = paste0(path, type, "_mut0.rda"))
        }

        ## process expression data:
        D <- assay(data)
        class <- data@colData@listData$definition
        M <- M0
        Mtype <- Mtype0
        colnames(M) <- lapply(colnames(M), function(x) {
            y <- unlist(strsplit(x, "-"))
            y <- paste(y[1:4], collapse = "-")
            y <- unlist(strsplit(y, ""))
            y <- paste(y[1:(length(y)-1)], collapse = "")
            return(y)
        })
        colnames(D) <- lapply(colnames(D), function(x) {
            y <- unlist(strsplit(x, "-"))
            y <- paste(y[1:4], collapse = "-")
            y <- unlist(strsplit(y, ""))
            y <- paste(y[1:(length(y)-1)], collapse = "")
            return(y)
        })
        colnames(M) <- gsub("A$", "", lapply(strsplit(colnames(M), "-"), function(x) { y <- x[1:4]; y <- paste(y, collapse = "-"); return(y) }))
        M <- M[order(rownames(M)), order(colnames(M))]
        Mtype <- Mtype[order(rownames(Mtype)), order(colnames(Mtype))]

        print("mutation done")

        ## log odds:
        if (file.exists(paste0(path, type, "_llr.rda")) & !newllr) {
            load(paste0(path, type, "_llr.rda"))
        } else {
            library(edgeR)
            if (sum(class %in% "Solid Tissue Normal") < 10) {
                distrParNC <- function(i, data) {
                    data[i, ] <- data[i, ] - median(data[, i])
                    llrcol <- numeric(ncol(data))
                    div0 <- quantile(data[i, ], 0.25)
                    div1 <- quantile(data[i, ], 0.75)
                    sigma <- sd(data[i, ])
                    for (j in 1:ncol(data)) {
                        if (data[i, j] <= 0) {
                            llrcol[j] <- log2(div0/data[i, j])
                        } else {
                            llrcol[j] <- log2(data[i, j]/div1)
                        }
                    }
                    return(llrcol)
                }
                highcounts <- which(apply(D, 1, median) >= 10)
                DC <- D[highcounts, ]
                genenames <- rownames(D)[highcounts, ]
                nf <- calcNormFactors(DC)
                DC <- t(t(DC)/nf) # DC <- DC2
                ## this is very adventurous:
                ## DC <- t(scale(t(DC)))
                ## DC <- abs(DC)
                ## pc <- 1-min(DC)
                ## tmp <- log2(DC+pc)
                ## hist(tmp)
                sfInit(parallel = T, cpus = 4)
                sfExport("DC")
                tmp <- do.call("rbind", sfLapply(1:nrow(DC),
                                                 distrParNC, DC))
                colnames(tmp) <- colnames(DC)
                rownames(tmp) <- genenames
                sfStop()
                DF <- DC
            } else {
                library(ks)
                mykcde <- function(x, k) {
                    a <- which.min(abs(k$eval.points - x))
                    b <- k$estimate[a]
                    b <- min(b, 1-b)
                    return(b)
                }
                distrParKs <- function(i, data, C) {
                    llrcol <- numeric(ncol(data))
                    ddistr <- list()
                    dogenes <- unique(colnames(data))[which(!(unique(colnames(data)) %in% ""))]
                    for (j in dogenes) {
                        D <- which(colnames(data) %in% j)
                        ddistr[[j]] <- kcde(data[i, D])
                    }
                    cdistr <- kcde(data[i, C])
                    for (j in which(!(colnames(data) %in% ""))) {
                        gene <- colnames(data)[j]
                        llrcol[j] <- log2(mykcde(data[i, j], ddistr[[gene]])/mykcde(data[i,j], cdistr))
                    }
                    llrcol <- llrcol[-C]
                    return(llrcol)
                }
                DN <- D[, which(class %in% "Solid Tissue Normal")]
                DT <- D[, which(class %in% "Primary solid Tumor")]
                DF <- cbind(DT, DN)
                C <- (ncol(DT)+1):ncol(DF)
                highcounts <- which(apply(DF, 1, median) >= 10)
                DF <- DF[highcounts, ]
                genenames <- rownames(DF)
                colnames(DF)[1:ncol(DT)] <- "P" # not knock-down specific!
                colnames(DF)[grep(paste(unique(gsub("-.*", "", types)), collapse = "|"), colnames(DF))] <- ""
                nf <- calcNormFactors(DF)
                DF <- t(t(DF)/nf)
                sfInit(parallel = T, cpus = 4)
                sfExport("DF", "C", "mykcde")
                sfLibrary(ks)
                tmp <- do.call("rbind", sfLapply(1:nrow(DF),
                                                 distrParKs, DF, C))
                sfStop()
                colnames(tmp) <- colnames(DT)
            }
            save(tmp, DF, file = paste0(path, type, "_llr.rda"))
        }

        rownames(tmp) <- rownames(DF)

        par(mfrow=c(1,3))
        hist(tmp)

        tmp[which(is.na(tmp) | is.infinite(tmp))] <- 0

        hist(tmp)

        D <- tmp

        print("expression done")

        ## sd.glob <- sd(tmp)

        ## tmp <- tmp[-which(apply(tmp, 1, sd) < sd.glob), ]

        ## hist(tmp)

        ## prep clinical data:
        clinical[which(clinical$vital_status%in% "dead"), which(colnames(clinical) %in% "vital_status")] <- 1
        clinical[which(clinical$vital_status%in% "alive"), which(colnames(clinical) %in% "vital_status")] <- 0
        count <- 0
        for (stage in sort(unique(clinical$tumor_stage))) {
            clinical[which(clinical$tumor_stage%in% stage), which(colnames(clinical) %in% "tumor_stage")] <- count
            count <- count + 1
        }
        clinical$tumor_stage <- as.numeric(clinical$tumor_stage)
        clinical[which(is.na(clinical$days_to_death)), which(colnames(clinical) %in% "days_to_death")] <- clinical[which(is.na(clinical$days_to_death)), which(colnames(clinical) %in% "days_to_last_follow_up")]
        clinical$vital_status<- as.numeric(clinical$vital_status)

        print("clinical done")

        save(clinical, D, M, Mtype, DF, class, meth, cnv, file = paste0(path, type, "_final.rda"))
    }
    print(table(class))
    sizemat <- rbind(sizemat, table(class))
    rownames(sizemat)[nrow(sizemat)] <- type
    if (donr) {
        samplenr[snrcount, 1] <- sum(class %in% "Primary solid Tumor")
        samplenr[snrcount, 2] <- sum(class %in% "Solid Tissue Normal")
    }
}

stop("done")

samplenr2 <- samplenr[-which(is.na(samplenr[, 1]) == TRUE), ]
barplot(t(samplenr2[order(apply(samplenr2, 1, sum)), ]), horiz = 1, space = 1, las = 2)

newllr <- 1; newmut <- 1; nmut <- 1; newsave <- 1; types <- c("TCGA-SKCM","TCGA-UVM"); source("~/Documents/testing/general/TCGA.r")

nonesolid <- c("TCGA-LAML")
solidnonormal <- c()

## analysis:

type <- "TCGA-BRCA"

path <- "mutclust/"

load(paste0(path, type, "_final.rda"))

M[which(M < 3)] <- 0
M[which(M > 0)] <- 1

M[grep("Silent", Mtype)] <- 0
M <- M[order(rownames(M)), order(colnames(M))]

cnv[which(cnv == 0)] <- 0
cnv[which(cnv != 0)] <- 1
cnv <- cnv[order(rownames(cnv)), order(colnames(cnv))]

meth[which(meth > 0.5)] <- 1
meth[which(meth <= 0.5)] <- 0
meth <- meth[order(rownames(meth)), order(colnames(meth))]
meth[is.na(meth)] <- 0

P <- matrix(0, length(unique(c(rownames(M), rownames(cnv), rownames(meth)))), length(unique(c(colnames(M), colnames(cnv), colnames(meth)))))
rownames(P) <- sort(unique(c(rownames(M), rownames(cnv), rownames(meth))))
colnames(P) <- sort(unique(c(colnames(M), colnames(cnv), colnames(meth))))
colnames(P) <- gsub("A$", "", lapply(strsplit(colnames(P), "-"), function(x) { y <- x[1:4]; y <- paste(y, collapse = "-"); return(y) }))
P <- P[, which(duplicated(colnames(P)) == FALSE)]

P[which(rownames(P) %in% rownames(M)), which(colnames(P) %in% colnames(M))] <- M[which(rownames(M) %in% rownames(P)), which(colnames(M) %in% colnames(P))]

Pmut <- P

P <- Pmut*0

P[which(rownames(P) %in% rownames(meth)), which(colnames(P) %in% colnames(meth))] <- P[which(rownames(P) %in% rownames(meth)), which(colnames(P) %in% colnames(meth))] + meth[which(rownames(meth) %in% rownames(P)), which(colnames(meth) %in% colnames(P))]

Pmeth <- P

P <- Pmeth*0

P[which(rownames(P) %in% rownames(cnv)), which(colnames(P) %in% colnames(cnv))] <- P[which(rownames(P) %in% rownames(cnv)), which(colnames(P) %in% colnames(cnv))] + cnv[which(rownames(cnv) %in% rownames(P)), which(colnames(cnv) %in% colnames(P))]

Pcnv <- P

P <- Pmut+Pmeth+Pcnv

P2 <- P # full abberations including cnv and meth

P <- Pmut

## Bailey et al Cell 2018

goi <- c("MAP2K4", "GATA3", "GPS2", "TBX3", "PTPRD", "NCOR1", "CBFB", "CDKN1B") # BRCA

P <- P[which(rownames(P) %in% goi), ]

P[which(P > 1)] <- 1
P <- apply(P, 2, function(x) return(x/sum(x)))
P[is.na(P)] <- 0

## data imputation:

library(naturalsort)
library(nem)
library(cluster)
library(Rcpp)
library(Rgraphviz)
library(mnem)

source("~/Documents/mnem/R/mnems.r")
source("~/Documents/mnem/R/mnems_low.r")
sourceCpp("~/Documents/mnem/src/mm.cpp")

source("~/Documents/nempi/R/nempi_main.r")
source("~/Documents/nempi/R/nempi_low.r")

Rho <- cbind(P, matrix(0, nrow(P), sum(!(colnames(D) %in% colnames(P)))))

colnames(Rho) <- c(colnames(P), colnames(D)[which(!(colnames(D) %in% colnames(P)))])

Rho <- Rho[, colnames(D)]

if (sum(apply(Rho, 1, sum) == 0) > 0) {
    Rho <- Rho[-which(apply(Rho, 1, sum) == 0), ]
}

Rho[is.na(Rho)] <- 0

sum(apply(Rho, 2, sum) == 0)/ncol(Rho) # unlabelled

pdf("temp.pdf", width = 12, height = 6)
tmp <- Rho
colnames(tmp) <- NULL
epiNEM::HeatmapOP(tmp, col = "RdBu", Rowv = 0, bordercol = "transparent")
dev.off()

inc <- sort(apply(Rho, 1, sum))

D2 <- D[which(apply(D, 1, median) != 0), ]

D3 <- D2[, which(duplicated(colnames(D2)) == FALSE)]
Rho <- Rho[, which(duplicated(colnames(D2)) == FALSE)]
for (i in which(duplicated(colnames(D2)) == TRUE)) {
    j <- which(colnames(D3) %in% colnames(D2)[i])
    D3[, j] <- apply(D2[, which(colnames(D2) %in% colnames(D2)[i]), drop = FALSE], 1, median)
}

D2 <- D3

colnames(D2) <- c(rownames(Rho), sample(rownames(Rho), ncol(D2)-nrow(Rho), replace = TRUE))

converged <- 10

start <- Sys.time()
nempires <- nempi(D2, Gamma = Rho, full = TRUE, converged = converged)
end <- Sys.time()

print(end - start)

ures <- nempires

sum(ures$lls[2:length(ures$lls)] - ures$lls[1:(length(ures$lls)-1)] < 0)

pdf("temp.pdf", width = 12, height = 6)
epiNEM::HeatmapOP(ures$Gamma, bordercol = rgb(0,0,0,0), col = "RdBu")
#plot(ures, edgewidth = 30)
dev.off()

pdf("temp.pdf", width = 9, height = 6)
par(mfrow=c(2,3))
plotConvergence(ures, type = "b", col = "blue")
dev.off()

source("~/Documents/nempi/R/nempi_main.r")
D4 <- D2
colnames(D4) <- apply(Rho, 2, function(x) {
    Sgenes <- paste(sort(rownames(Rho)[which(x > 0)]), collapse = "_")
    return(Sgenes)
})

## run on hpc cluster:

path <- ""
type <- "TCGA-BRCA"

if (do == 1) {

    library(e1071)

    load(paste0(path, type, "_nempi.rda"))

    D4 <- D2
    colnames(D4) <- apply(Rho, 2, function(x) {
        Sgenes <- paste(sort(rownames(Rho)[which(x > 0)]), collapse = "_")
        return(Sgenes)
    })

    svmres <- classpi(D4, full = TRUE, method = "svm")

    save(svmres, file = paste0("temp_", do, "_", as.numeric(Sys.time()), ".rda"))

}

if (do == 2) {

    library(nnet)

    load(paste0(path, type, "_nempi.rda"))

    D4 <- D2
    colnames(D4) <- apply(Rho, 2, function(x) {
        Sgenes <- paste(sort(rownames(Rho)[which(x > 0)]), collapse = "_")
        return(Sgenes)
    })

    nnres <- classpi(D4, full = TRUE, method = "nnet", MaxNWts = 50000, size = 5) # takes forever

    save(nnres, file = paste0("temp_", do, "_", as.numeric(Sys.time()), ".rda"))

}

if (do == 3) {

    library(CALIBERrfimpute)

    load(paste0(path, type, "_nempi.rda"))

    D4 <- D2
    colnames(D4) <- apply(Rho, 2, function(x) {
        Sgenes <- paste(sort(rownames(Rho)[which(x > 0)]), collapse = "_")
        return(Sgenes)
    })

    mfdata <- cbind(as.data.frame(t(D4)), colnames(D4))
    mfdata[which(mfdata == "", arr.ind = TRUE)] <- NA

    micedata <- mfdata
    colnames(micedata) <- paste0(LETTERS[1:ncol(micedata)], 1:ncol(micedata))
    miceres <- mice(micedata, method = c(rep('rfcont', ncol(micedata)-1), 'rfcat'), m = 2, maxit = 2)

    save(miceres, file = paste0("temp_", do, "_", as.numeric(Sys.time()), ".rda"))

}

if (do == 4) {

    library(e1071)

    load(paste0(path, type, "_nempi.rda"))

    D4 <- D2
    colnames(D4) <- apply(Rho, 2, function(x) {
        Sgenes <- paste(sort(rownames(Rho)[which(x > 0)]), collapse = "_")
        return(Sgenes)
    })

    rfres <- classpi(D4, full = TRUE, method = "randomForest")

    save(rfres, file = paste0("temp_", do, "_", as.numeric(Sys.time()), ".rda"))

}

if (do == 5) {

    library(e1071)

    load(paste0(path, type, "_nempi.rda"))

    D4 <- D2
    colnames(D4) <- apply(Rho, 2, function(x) {
        Sgenes <- paste(sort(rownames(Rho)[which(x > 0)]), collapse = "_")
        return(Sgenes)
    })

    mfdata <- cbind(as.data.frame(t(D4)), colnames(D4))
    mfdata[which(mfdata == "", arr.ind = TRUE)] <- NA
    library(missForest)
    mfimp <- missForest(mfdata)
    D4 <- D2
    colnames(D4) <- mfimp$ximp[, ncol(mfimp$ximp)]
    tmp <- mynem(D4, multi = TRUE)
    Gamma <- getGamma(D4)
    ures <- list()
    ures$Gamma <- apply(Gamma, 2, function(x) return(x/sum(x)))
    ures$res <- list()
    ures$res$adj <- tmp$adj
    ures$null <- TRUE
    ures$combi <- 1
    mfres <- ures

    save(mfres, file = paste0("temp_", do, "_", as.numeric(Sys.time()), ".rda"))

}

## knn

library(class)
train <- t(D4[, which(colnames(D4) != "")])
test <- t(D4[, which(colnames(D4) == "")])
knn0 <- 0
if (knn0) {
    train <- rbind(train, NULL = rep(0, ncol(train)))
    cl <- c(colnames(D4)[which(colnames(D4) != "")], "NULL")
} else {
    cl <- colnames(D4)[which(colnames(D4) != "")]
}
knnres <- knn(train, test, cl, prob = TRUE)
D3 <- D4
colnames(D3)[which(colnames(D3) %in% "")] <- as.character(knnres)
tmp <- mynem(D3, multi = TRUE)
Gamma <- getGamma(D3)
ures <- list()
ures$Gamma <- Gamma # apply(Gamma, 2, function(x) return(x/sum(x)))
ures$res <- list()
ures$res$adj <- tmp$adj
ures$null <- TRUE
ures$combi <- 1

knnres <- ures

## save(nempires, knnres, rfres, mfres, svmres, nnres, Rho, D2, Pmut, Pmeth, Pcnv, file = paste0(path, type, "_nempi.rda"))

path <- "mutclust/"; type <- "TCGA-BRCA"

load(paste0(path, type, "_nempi.rda"))

## ## load("~/Mount/Euler/temp_1_1573814974.40659.rda") # old

## load("~/Mount/Euler/temp_1_1574076694.65703.rda")

## ## load("~/Mount/Euler/temp_4_1573819581.9749.rda") # old

## load("~/Mount/Euler/temp_4_1574080528.67352.rda")

## ## load("~/Mount/Euler/temp_2_1573821112.2412.rda") # old

## load("~/Mount/Euler/temp_2_1574084503.12431.rda")

## load("~/Mount/Euler/temp_5_1574081415.91547.rda")

## ures <- rfres

## ures <- mfres

## ures <- svmres

## ures <- nnres

## ures <- knnres

ures <- nempires

## check against methylation and cnvs:

pdf("nempi_gamma.pdf", width = 12, height = 6)
tmp <- ures$Gamma
colnames(tmp) <- NULL
epiNEM::HeatmapOP(tmp, bordercol = rgb(0,0,0,0), col = "RdBu", colorkey = NULL)
dev.off()

pdf("nempi_phi.pdf", width = 6, height = 6)
Pgenes <- sort(unique(colnames(D2)))
adjtmp <- ures$res$adj
colnames(adjtmp) <- rownames(adjtmp) <- Pgenes
plotDnf(adjtmp, edgelwd = 2)
dev.off()

## cnv/meth enrichment:

methods <- list("NEM$\\pi$" = nempires, knn = knnres, mf = mfres, nn = nnres, rf = rfres, svm = svmres)

mutinc <- 1

Lall <- Lcnv <- Lmeth <- Lmut <- list()

for (i in 1:length(methods)) {
    if (i != 8) {
        print(names(methods)[i])
        ures <- methods[[i]]
        newGamma <- ures$Gamma
    } else {
        print("random")
        newGamma <- newGamma*0
        newGamma[sample(1:length(newGamma), floor(0.45*length(newGamma)))] <- 1 # well that is included into the test...
    }
    hist(newGamma)
    if (i == 1) {
        rntmp <- rownames(newGamma); newGamma <- t(mytc(ures$res$adj))%*%newGamma; rownames(newGamma) <- rntmp
    }
    P <- Pmut+Pmeth+Pcnv
    if (!mutinc) { # include mutations or not (0)
        P[which(Pmut == 1)] <- 0
    }
    P <- P[which(rownames(P) %in% rownames(Rho)), which(colnames(P) %in% colnames(Rho))]
    P <- P[order(rownames(P)), order(colnames(P))]
    P[which(P > 1)] <- 1
    Ptmp <- cbind(P, matrix(0, nrow(P), sum(!(colnames(Rho) %in% colnames(P)))))
    colnames(Ptmp) <- c(colnames(P), colnames(Rho)[which(!(colnames(Rho) %in% colnames(P)))])
    P <- Ptmp[, colnames(Rho)]
    ## fisher:
    if (i %in% c(4,5)) {
        cut <- 0.07
    } else if (i == 6) {
        cut <- 0.1
    } else {
        cut <- 1/8
    }
    newGamma[which(newGamma > cut)] <- 1
    newGamma[which(newGamma <= cut)] <- 0
    pmeth <- newGamma
    pmeth[which(newGamma == 1 & P == 1)] <- 2
    pmeth[which(newGamma == 0 & P == 1)] <- -2
    if (!mutinc) { # include mutations or not (0)
        pmeth[which(Rho > 0)] <- 0
    }
    colnames(pmeth) <- NULL
    ##pdf(paste0("FigS_", names(methods)[i], ".pdf"), height = 6, width = 12)
    setEPS()
    postscript(paste0("FigS_", names(methods)[i], ".eps"), height = 6, width = 12)
    print(epiNEM::HeatmapOP(pmeth, bordercol = rgb(0,0,0,0), col = "RdBu", colorkey = NULL))
    dev.off()
    print("cnv")
    P <- Pcnv
    P <- P[which(rownames(P) %in% rownames(Rho)), which(colnames(P) %in% colnames(Rho))]
    P <- P[order(rownames(P)), order(colnames(P))]
    P[which(P > 1)] <- 1
    Ptmp <- cbind(P, matrix(0, nrow(P), sum(!(colnames(Rho) %in% colnames(P)))))
    colnames(Ptmp) <- c(colnames(P), colnames(Rho)[which(!(colnames(Rho) %in% colnames(P)))])
    P <- Ptmp[, colnames(Rho)]
    F <- matrix(c(sum(pmeth >= 1 & P == 1), sum(pmeth >= 1 & P == 0), sum(pmeth == -2 & P == 1), sum(pmeth == 0 & P == 0)), 2)
    print(1 - phyper(F[1,1]-1, sum(F[, 1]), sum(F[, 2]), sum(F[1, ])))
    Lcnv[[i]] <- F
    print("meth")
    P <- Pmeth
    P <- P[which(rownames(P) %in% rownames(Rho)), which(colnames(P) %in% colnames(Rho))]
    P <- P[order(rownames(P)), order(colnames(P))]
    P[which(P > 1)] <- 1
    Ptmp <- cbind(P, matrix(0, nrow(P), sum(!(colnames(Rho) %in% colnames(P)))))
    colnames(Ptmp) <- c(colnames(P), colnames(Rho)[which(!(colnames(Rho) %in% colnames(P)))])
    P <- Ptmp[, colnames(Rho)]
    F <- matrix(c(sum(pmeth >= 1 & P == 1), sum(pmeth >= 1 & P == 0), sum(pmeth == -2 & P == 1), sum(pmeth == 0 & P == 0)), 2)
    print(1 - phyper(F[1,1]-1, sum(F[, 1]), sum(F[, 2]), sum(F[1, ])))
    Lmeth[[i]] <- F
    print("mut")
    P <- Pmut
    P <- P[which(rownames(P) %in% rownames(Rho)), which(colnames(P) %in% colnames(Rho))]
    P <- P[order(rownames(P)), order(colnames(P))]
    P[which(P > 1)] <- 1
    Ptmp <- cbind(P, matrix(0, nrow(P), sum(!(colnames(Rho) %in% colnames(P)))))
    colnames(Ptmp) <- c(colnames(P), colnames(Rho)[which(!(colnames(Rho) %in% colnames(P)))])
    P <- Ptmp[, colnames(Rho)]
    F <- matrix(c(sum(pmeth >= 1 & P == 1), sum(pmeth >= 1 & P == 0), sum(pmeth == -2 & P == 1), sum(pmeth == 0 & P == 0)), 2)
    print(1 - phyper(F[1,1]-1, sum(F[, 1]), sum(F[, 2]), sum(F[1, ])))
    Lmut[[i]] <- F
    Fmat <- matrix(c(sum(pmeth == 2), sum(pmeth == 1), sum(pmeth == -2), sum(pmeth == 0)), 2)
    Lall[[i]] <- Fmat
    ## print(fisher.test(Fmat, alternative = "greater"))
    print("p-value")
    print(1 - phyper(Fmat[1,1]-1, sum(Fmat[, 1]), sum(Fmat[, 2]), sum(Fmat[1, ])))
}

## create tables:

for (i in 1:length(methods)) {
    cat(paste0(names(methods)[i], " & ", Lcnv[[i]][1,1], " & ", Lcnv[[i]][2,1], " & ", Lcnv[[i]][2,2], " & ", Lcnv[[i]][1,2], "\\\\\n"))
}

for (i in 1:length(methods)) {
    cat(paste0(names(methods)[i], " & ", Lmeth[[i]][1,1], " & ", Lmeth[[i]][2,1], " & ", Lmeth[[i]][2,2], " & ", Lmeth[[i]][1,2], "\\\\\n"))
}

for (i in 1:length(methods)) {
    cat(paste0(names(methods)[i], " & ", Lmut[[i]][1,1], " & ", Lmut[[i]][2,1], " & ", Lmut[[i]][2,2], " & ", Lmut[[i]][1,2], "\\\\\n"))
}

for (i in 1:length(methods)) {
    if (names(methods)[i] == "nn") { next() }
    Fmat <- Lall[[i]]
    ptmp <- 1 - phyper(Fmat[1,1]-1, sum(Fmat[, 1]), sum(Fmat[, 2]), sum(Fmat[1, ]))
    if (ptmp == 0) {
        ptmp <- "$< 2.2\\times10^{-16}$"
    } else {
        ptmp <- paste0("$", signif(ptmp), "$")
    }
    cat(paste0(names(methods)[i], " & ", Lall[[i]][1,1], " & ", Lall[[i]][2,1], " & ", Lall[[i]][2,2], " & ", Lall[[i]][1,2], " & ", ptmp, "\\\\\n"))
}


## check correlation

## P4 <- apply(P3, 2, function(x) return(x/sum(x)))
## P4[is.na(P4)] <- 0

## cor(as.vector(newGamma), as.vector(P3))

##

cormat <- matrix(0, nrow(pmeth), 2)
fishres <- numeric(nrow(pmeth))
names(fishres) <- rownames(pmeth)
for (i in 1:nrow(pmeth)) {
    Fmat <- matrix(c(sum(pmeth[i, ] == 2), sum(pmeth[i, ] == 1), sum(pmeth[i, ] == -2), sum(pmeth[i, ] == 0)), 2)
    fishres[i] <- fisher.test(Fmat, alternative = "g")$p.value
    cormat[i, ] <- c(sum(Fmat[1, ]), sum(Fmat[, 1]))
}

## GATA3 & PTPRD:

Fmat <- matrix(c(sum(pmeth[3, ] %in% c(2,-2) & pmeth[7, ] %in% c(2,-2)),
                 sum(pmeth[3, ] %in% c(1,0) & pmeth[7, ] %in% c(-2,2)),
                 sum(pmeth[3, ] %in% c(-2,2) & pmeth[7, ] %in% c(1,0)),
                 sum(pmeth[3, ] %in% c(1,0) & pmeth[7, ] %in% c(1,0))), 2)

1 - phyper(Fmat[1,1]-1, sum(Fmat[, 1]), sum(Fmat[, 2]), sum(Fmat[1, ]))


## pca:

pca <- princomp(D2)

col <- apply(newGamma, 2, sum)

col <- col/max(col)

K <- kmeans(t(newGamma), nrow(newGamma))

plot(pca$loadings[, 1:2], col = K$cluster)#rgb(col,0,0,1))

## tsne:

sne <- tsne(t(newGamma))

plot(sne, col = K$cluster)

## R profiling:

Rprof("temp.txt", line.profiling=TRUE)
ures <- nempi(D2[1:20, ], Gamma = Gamma, complete = 1, full = TRUE, converged = converged, combi = combi)
Rprof(NULL)
summaryRprof("temp.txt", lines = "show")$sampling.time
head(summaryRprof("temp.txt", lines = "show")$by.self, 10)

##

type <- "TCGA-BRCA"
load(paste0(path, type, "_nempi.rda"))

pdf("Fig5.pdf", width = 10, height = 10)
plot(uresn, edgelwd = 2)
dev.off()

pdf("Fig6.pdf", width = 10, height = 5)
tmp <- uresn$Gamma
colnames(tmp) <- NULL
epiNEM::HeatmapOP(tmp, bordercol = rgb(0,0,0,0), col = "RdBu", colorkey = NULL)
dev.off()

pdf("Fig7.pdf", width = 10, height = 10)
phitmp <- mytc(uresn$res$adj)
tmp <- t(phitmp)%*%uresn$Gamma
colnames(tmp) <- NULL
rownames(tmp) <- rownames(uresn$Gamma)
tmp2 <- tmp
colnames(tmp) <- NULL
tmp <- Gamma
colnames(tmp) <- NULL
tmp3 <- tmp
tmp4 <- tmp2
p1 <- epiNEM::HeatmapOP(tmp2, bordercol = rgb(0,0,0,0), col = "RdBu", clusterx = tmp2, colorkey = NULL)
p2 <- epiNEM::HeatmapOP(tmp, bordercol = rgb(0,0,0,0), col = "RdBu", clusterx = tmp2, colorkey = NULL)
print(p1, position=c(0, .5, 1, 1), more=TRUE)
print(p2, position=c(0, 0, 1, .5))
sum(tmp == 1 & tmp2 == 1)/sum(tmp == 1)
dev.off()

pdf("Fig8.pdf", width = 10, height = 10)
plot(uresf, edgelwd = 2)
dev.off()

pdf("Fig9.pdf", width = 10, height = 5)
tmp <- uresf$Gamma
colnames(tmp) <- NULL
epiNEM::HeatmapOP(tmp, bordercol = rgb(0,0,0,0), col = "RdBu", colorkey = NULL)
dev.off()

pdf("Fig10.pdf", width = 10, height = 10)
phitmp <- mytc(uresf$res$adj)
tmp <- t(phitmp)%*%uresf$Gamma
colnames(tmp) <- NULL
rownames(tmp) <- rownames(uresf$Gamma)
tmp2 <- tmp
colnames(tmp) <- NULL
tmp <- Gamma
colnames(tmp) <- NULL
p1 <- epiNEM::HeatmapOP(tmp2, bordercol = rgb(0,0,0,0), col = "RdBu", clusterx = tmp2, colorkey = NULL)
p2 <- epiNEM::HeatmapOP(tmp, bordercol = rgb(0,0,0,0), col = "RdBu", clusterx = tmp2, colorkey = NULL)
print(p1, position=c(0, .5, 1, 1), more=TRUE)
print(p2, position=c(0, 0, 1, .5))
sum(tmp == 1 & tmp2 == 1)/sum(tmp == 1)
dev.off()

## new figure:

P <- t(mytc(uresf$res$adj))%*%uresf$Gamma
rownames(P) <- rownames(uresf$Gamma)
PM <- P
PM[which(PM > 1/6 & Gamma == 1)] <- P[which(PM > 1/6 & Gamma == 1)] + 1
PM[which(PM <= 1/6 & Gamma == 1)] <- P[which(PM <= 1/6 & Gamma == 1)] - 1
epiNEM::HeatmapOP(PM, bordercol = rgb(0,0,0,0), col = "RdYlBu", breaks = seq(-1,2,length.out=5), clusterx = tmp2)

## other plots:

pdf("temp.pdf", width = 10, height = 10)
plotDnf(c("M1=M4", "M2=M4", "M3=M5", "M1=M5"), edgelwd = 2)
dev.off()

M <- matrix(0, 5, 10)
rownames(M) <- paste0("M", 1:5)
colnames(M) <- paste0("S", 1:10)

M[1, 1:4] <- 1
M[2, c(2,7:9)] <- 1
M[3, 5:6] <- 1
M[3, 10] <- 1

phi <- matrix(0, 5, 5)
diag(phi) <- 1
phi[1, 1:5] <- 1
phi[2, 3] <- phi[4, 5] <- 1

## M <- t(phi)%*%M; M[M > 1] <- 1

rownames(M) <- paste0("M", 1:5)

pdf("temp.pdf", width = 8, height = 4)
epiNEM::HeatmapOP(M, Colv = 0, Rowv = 0, col = "RdBu", colorkey = NULL)
dev.off()

## check for mutation type (das führt zu nichts)

colnames(Mtype) <- unlist(lapply(strsplit(colnames(Mtype), "-"), function(x) {
    y <- x[1:3]
    y <- paste(c(y, "01"), collapse = "-")
    return(y)
    }))

checkgene <- "CBFB"

van <- intersect(colnames(Gamma)[which(PM[checkgene, ] < 0)], colnames(Mtype))

A <- sum(unlist(lapply(strsplit(Mtype[checkgene, van], " "), function(x) if ("Silent" %in% x) { return(1) } else { return(0) })))
B <- length(van) - A

van <- intersect(colnames(Gamma)[which(PM[checkgene, ] > 0)], colnames(Mtype))

C <- sum(unlist(lapply(strsplit(Mtype[checkgene, van], " "), function(x) if ("Silent" %in% x) { return(1) } else { return(0) })))
D <- length(van) - C

table(unlist(lapply(strsplit(Mtype[checkgene, van], " "), function(x) return(names(table(x))[which.max(table(x))]))))

## table(unlist(strsplit(as.vector(Mtype), " ")))

## pdf("Fig10.pdf", width = 10, height = 10)
## p1 <- epiNEM::HeatmapOP(tmp4, bordercol = rgb(0,0,0,0), col = "RdBu", clusterx = tmp2, colorkey = NULL)
## p2 <- epiNEM::HeatmapOP(tmp3, bordercol = rgb(0,0,0,0), col = "RdBu", clusterx = tmp2, colorkey = NULL)
## p3 <- epiNEM::HeatmapOP(tmp2, bordercol = rgb(0,0,0,0), col = "RdBu", clusterx = tmp2, colorkey = NULL)
## p4 <- epiNEM::HeatmapOP(tmp, bordercol = rgb(0,0,0,0), col = "RdBu", clusterx = tmp2, colorkey = NULL)
## print(p1, position=c(0, .5, .5, 1), more=TRUE)
## print(p2, position=c(0, 0, .5, .5), more=TRUE)
## print(p3, position=c(.5, .5, 1, 1), more=TRUE)
## print(p4, position=c(.5, 0, 1, .5))
## sum(tmp == 1 & tmp2 == 1)/sum(tmp == 1)
## dev.off()

source("~/Documents/testing/R/nempi.r")
source("~/Documents/testing/R/nempi_low.r")

bsres <- unembs(D2, Gamma = Gamma, complete = 1, full = TRUE, converged = converged, combi = combi, bsruns = 10, bssize = 0.5)

pdf("temp.pdf", width = 10, height = 5)
epiNEM::HeatmapOP(bsres$Gamma, bordercol = rgb(0,0,0,0), col = "RdBu", colorkey = NULL)
dev.off()

pdf("temp.pdf", width = 20, height = 10)
par(mfrow=c(1,2))
tmp <- bsres$phi
tmp[which(tmp > 0)] <- 1
diag(tmp) <- 0
freqtop <- as.vector(t(bsres$phi)[which(lower.tri(bsres$phi) == TRUE)])/10
freqtop <- freqtop[which(freqtop != 0)]
tmptop <- tmp
tmptop[lower.tri(tmptop)] <- 0
tmptop <- adj2dnf(tmptop)
tmptop <- tmptop[grep("=", tmptop)]
plotDnf(tmptop, edgelwd = 2, edgelab = freqtop, freq = freqtop)
freqbot <- as.vector(t(bsres$phi)[which(upper.tri(bsres$phi) == TRUE)])/10
freqbot <- freqbot[which(freqbot != 0)]
tmpbot <- tmp
tmpbot[upper.tri(tmpbot)] <- 0
tmpbot <- adj2dnf(tmpbot)
tmpbot <- tmpbot[grep("=", tmpbot)]
plotDnf(tmpbot, edgelwd = 2, edgelab = freqbot, freq = freqbot)
dev.off()

##

source("testing/vignettes/TCGA_cluster.r")

pcares <- prcomp(tmp)

library(naturalsort)
library(nem)
library(cluster)
library(Rcpp)

source("~/Documents/mnem/R/mnems.r")
source("~/Documents/mnem/R/mnems_low.r")
sourceCpp("~/Documents/mnem/src/mm.cpp")

res <- mnem(tmp, starts = 10, search = "greedy", type = "cluster3", complete = 1, multi = 1)

tmp2 <- tmp
Rprof("temp.txt", line.profiling=TRUE)
res <- mnem(tmp2, starts = 10, search = "greedy", type = "cluster3", complete = 1, multi = 1, k = 2)
Rprof(NULL)
summaryRprof("temp.txt", lines = "show")$sampling.time
head(summaryRprof("temp.txt", lines = "show")$by.self)

## resk <- mnemk(tmp2, starts = 10, search = "estimate", type = "cluster3", complete = 1, multi = 1)

cluster <- apply(getAffinity(res$probs, mw = res$mw, complete = TRUE), 2, which.max)

par(mfrow=c(1,2))
plot(pcares$rotation[, 1:2], col = cluster)

names(cluster) <- gsub("-01$|-03$", "", colnames(M))

cluster <- cluster[which(names(cluster) %in% clinical$submitter_id)]

cluster <- cluster[order(names(cluster))]

clinical <- rbind(clinical, clinical[which(clinical[, 1] %in% names(cluster)[which(duplicated(names(cluster)))]), ])

clinical <- clinical[order(clinical[, 1]), ]

print(all(clinical[, 1] == names(cluster)))

fit <- survfit(Surv(days_to_death, vital_status) ~ cluster, clinical)
plot(fit, col = 1:length(table(cluster)), lty = 1:length(table(cluster)))
legend(max(clinical$days_to_death, na.rm = TRUE), 1, 1:length(table(cluster)), lty = 1:length(table(cluster)), col = 1:length(table(cluster)), xjust = 1, yjust = 1)

fit <- coxph(Surv(days_to_death, vital_status) ~ cluster + age_at_diagnosis, clinical)
print(fit)

fit <- coxph(Surv(days_to_death, vital_status) ~ cluster + age_at_diagnosis + tumor_stage, clinical)
print(fit)

fit <- coxph(Surv(days_to_death, vital_status) ~ cluster + tumor_stage, clinical)
print(fit)

fit <- coxph(Surv(days_to_death, vital_status) ~ cluster, clinical)
print(fit)

kres <- clustNEM(tmp, nem = 0, k = length(res$comp), nstart = 10)

## kres <- clustNEM(tmp, nem = 0, nstart = 10)

plot(pcares$rotation[, 1:2], col = kres$cluster)

kcluster <- kres$cluster

names(kcluster) <- gsub("-01$|-03$", "", colnames(M))

kcluster <- kcluster[which(names(kcluster) %in% clinical$submitter_id)]

kcluster <- kcluster[order(names(kcluster))]

fit <- survfit(Surv(days_to_death, vital_status) ~ kcluster, clinical)
plot(fit, col = 1:length(table(cluster)), lty = 1:length(table(cluster)))
legend(max(clinical$days_to_death, na.rm = TRUE), 1, 1:length(table(cluster)), lty = 1:length(table(cluster)), col = 1:length(table(cluster)), xjust = 1, yjust = 1)

fit <- coxph(Surv(days_to_death, vital_status) ~ kcluster + age_at_diagnosis, clinical)
print(fit)

fit <- coxph(Surv(days_to_death, vital_status) ~ kcluster, clinical)
print(fit)

















