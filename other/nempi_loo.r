
pipe <- as.character(commandArgs(TRUE)[1])

type <- as.character(commandArgs(TRUE)[2])

bsruns <- as.numeric(commandArgs(TRUE)[3]) # bsruns <- 10

startrun <- as.numeric(commandArgs(TRUE)[4]) # startrun <- 10

if (type %in% c("rf", "svm", "nn")) {
    options(expressions = 500000)
}

safeSave <- function(..., list = character(), file = stop("'file' must be specified"),
                     ascii = FALSE, version = NULL, envir = parent.frame(), compress = isTRUE(!ascii),
                     compression_level, eval.promises = TRUE, precheck = TRUE) {

    i <- 1
    while(file.exists(paste0(file, "_", i))) {
        i <- i + 1
    }

    save(..., list = list, file = paste0(file, "_", i),
         ascii = ascii, version = version, envir = envir, compress = compress,
         compression_level = compression_level, eval.promises = eval.promises, precheck = precheck)
}

aucs <- function(a,b) {
    auc <- roc <- 0
    ppv <- rec <- spec <- NULL
    for (cut in c(2,seq(1,0, length.out = 100),-1)) {
        tmp <- a*0
        tmp[which(a > cut)] <- 1
        tp <- sum(tmp == 1 & b == 1)
        fp <- sum(tmp == 1 & b == 0)
        tn <- sum(tmp == 0 & b == 0)
        fn <- sum(tmp == 0 & b == 1)
        ppvtmp <-  tp/(tp+fp)
        if (is.na(ppvtmp)) { ppvtmp <- 0.5 }
        rectmp <- tp/(tp+fn)
        spectmp <- 1-tn/(tn+fp)
        if (length(ppv) > 0) {
            auc <- auc + (rectmp-rec[length(rec)])*(ppvtmp+ppv[length(ppv)])/2
            roc <- roc + (spectmp-spec[length(spec)])*(rectmp+rec[length(rec)])/2
        }
        ppv <- c(ppv, ppvtmp)
        rec <- c(rec, rectmp)
        spec <- c(spec, spectmp)
    }
    return(list(auc=auc,roc=roc,ppv=ppv,rec=rec,spec=spec))
}

path <- "/cluster/work/bewi/members/mpirkl/"

if (pipe %in% "loo") { 
    library(naturalsort)
    library(nem)
    library(cluster)
    library(Rcpp)
    library(Rgraphviz)
    library(randomForest)
    library(e1071)
    source("mnems.r")
    source("mnems_low.r")
    sourceCpp(code=readChar("mm.cpp", file.info("mm.cpp")$size))
    source("nempi_main.r")
    source("nempi_low.r")
    load("TCGA-BRCA_nempi.rda")
    ## sourceCpp("~/Documents/mnem/src/mm.cpp"); source("~/Documents/mnem/R/mnems.r"); source("~/Documents/mnem/R/mnems_low.r"); source("~/Documents/nempi/R/nempi_main.r"); source("~/Documents/nempi/R/nempi_low.r")
    ## load("nempi_backup/backup3/TCGA-BRCA_nempi.rda")
    P <- Pmut
    P <- P[which(rownames(P) %in% rownames(Rho)), which(colnames(P) %in% colnames(Rho))]
    P <- P[order(rownames(P)), order(colnames(P))]
    P[which(P > 1)] <- 1
    Ptmp <- cbind(P, matrix(0, nrow(P), sum(!(colnames(Rho) %in% colnames(P)))))
    colnames(Ptmp) <- c(colnames(P), colnames(Rho)[which(!(colnames(Rho) %in% colnames(P)))])
    P <- Ptmp[, colnames(Rho)]
    D4 <- D2
    colnames(D4) <- apply(Rho, 2, function(x) {
        Sgenes <- paste(sort(rownames(Rho)[which(x > 0)]), collapse = "_")
        return(Sgenes)
    })
    D5 <- D4[, which(apply(P, 2, sum) != 0)]
    genes.med <- apply(D5, 1, median)
    D5 <- D5[which(genes.med != 0), ]
    PM <- P[, which(apply(P, 2, sum) != 0)]
    tp <- fp <- fn <- tn <- prauc <- numeric(ncol(D5))
    if (type %in% "nempi") {
        converged <- 1
        for (i in startrun:(startrun+bsruns-1)) {
            if (i > ncol(D5)) { break() }
            D6 <- D5[, -i]
            Rho2 <- PM[, -i]
            Rho3 <- PM
            Rho3[, i] <- 0
            nempitmp <- nempi(D6, Gamma = Rho2, full = TRUE, converged = converged)
            Rho4 <- Rho3
            Rho4[, -i] <- nempitmp$Gamma
            nempitmp <- nempi(D5, Gamma = Rho4, full = FALSE, converged = converged)
            Rho5 <- t(mytc(nempitmp$res$adj))%*%nempitmp$Gamma
            tp[i] <- sum(PM[, i] == 1 & Rho5[, i] >= 1/8)
            fp[i] <- sum(PM[, i] == 0 & Rho5[, i] >= 1/8)
            fn[i] <- sum(PM[, i] == 1 & Rho5[, i] < 1/8)
            tn[i] <- sum(PM[, i] == 0 & Rho5[, i] < 1/8)
            prauc[i] <- aucs(Rho5[, i],PM[, i])$auc
        }
        safeSave(prauc, tp, fp, fn, tn, file = paste0(path, "nempi_bs/nempi_loo.rda"))
    } else if (type %in% "svm") {
        for (i in startrun:(startrun+bsruns-1)) {
            if (i > ncol(D5)) { break() }
            D6 <- D5[, -i]
            Rho2 <- PM[, -i]
            Rho3 <- PM
            Rho3[, i] <- 0
            svmres <- classpi(D6, full = TRUE, method = "svm")
            Rho4 <- Rho3
            Rho4[, -i] <- svmres$Gamma
            D7 <- D5
            colnames(D7)[i] <- ""
            svmres <- classpi(D7, full = FALSE, method = "svm")
            Rho5 <- svmres$Gamma
            tp[i] <- sum(PM[, i] == 1 & Rho5[, i] >= 1/8)
            fp[i] <- sum(PM[, i] == 0 & Rho5[, i] >= 1/8)
            fn[i] <- sum(PM[, i] == 1 & Rho5[, i] < 1/8)
            tn[i] <- sum(PM[, i] == 0 & Rho5[, i] < 1/8)
            prauc[i] <- aucs(Rho5[, i],PM[, i])$auc
        }
        safeSave(prauc, tp, fp, fn, tn, file = paste0(path, "nempi_bs/svm_loo.rda"))
    } else if (type %in% "rf") {
        for (i in startrun:(startrun+bsruns-1)) {
            if (i > ncol(D5)) { break() }
            D6 <- D5[, -i]
            Rho2 <- PM[, -i]
            Rho3 <- PM
            Rho3[, i] <- 0
            svmres <- classpi(D6, full = TRUE, method = "randomForest")
            Rho4 <- Rho3
            Rho4[, -i] <- svmres$Gamma
            D7 <- D5
            colnames(D7)[i] <- ""
            svmres <- classpi(D7, full = FALSE, method = "randomForest")
            Rho5 <- svmres$Gamma
            tp[i] <- sum(PM[, i] == 1 & Rho5[, i] >= 1/8)
            fp[i] <- sum(PM[, i] == 0 & Rho5[, i] >= 1/8)
            fn[i] <- sum(PM[, i] == 1 & Rho5[, i] < 1/8)
            tn[i] <- sum(PM[, i] == 0 & Rho5[, i] < 1/8)
            prauc[i] <- aucs(Rho5[, i],PM[, i])$auc
        }
        safeSave(prauc, tp, fp, fn, tn, file = paste0(path, "nempi_bs/rf_loo.rda"))
    } else if (type %in% "nn") {
        for (i in startrun:(startrun+bsruns-1)) {
            if (i > ncol(D5)) { break() }
            D6 <- D5[, -i]
            Rho2 <- PM[, -i]
            Rho3 <- PM
            Rho3[, i] <- 0
            library(nnet)
            svmres <- classpi(D6, full = TRUE, method = "nnet", MaxNWts = 50000, size = 1)
            Rho4 <- Rho3
            Rho4[, -i] <- svmres$Gamma
            D7 <- D5
            colnames(D7)[i] <- ""
            svmres <- classpi(D7, full = FALSE, method = "nnet", MaxNWts = 50000, size = 1)
            Rho5 <- svmres$Gamma
            tp[i] <- sum(PM[, i] == 1 & Rho5[, i] >= 1/8)
            fp[i] <- sum(PM[, i] == 0 & Rho5[, i] >= 1/8)
            fn[i] <- sum(PM[, i] == 1 & Rho5[, i] < 1/8)
            tn[i] <- sum(PM[, i] == 0 & Rho5[, i] < 1/8)
            prauc[i] <- aucs(Rho5[, i],PM[, i])$auc
        }
        safeSave(prauc, tp, fp, fn, tn, file = paste0(path, "nempi_bs/nn_loo.rda"))
    } else if (type %in% "mf") {
        for (i in startrun:(startrun+bsruns-1)) {
            if (i > ncol(D5)) { break() }
            D6 <- D5
            colnames(D6)[i] <- ""
            Rho2 <- PM
            Rho2[, i] <- 0
            Rho3 <- PM
            Rho3[, i] <- 0
            mfdata <- cbind(as.data.frame(t(D6)), factor(colnames(D6)))
            mfdata[which(mfdata == "", arr.ind = TRUE)] <- NA
            library(missForest)
            mfimp <- missForest(mfdata)
            colnames(D6) <- mfimp$ximp[, ncol(mfimp$ximp)]
            tmp <- nem(D6)
            ures <- list()
            ures$Gamma <- apply(getGamma(D6), 2, function(x) return(x/sum(x)))
            ures$res <- list()
            ures$res$adj <- tmp$adj
            ures$null <- TRUE
            ures$combi <- 1
            Rho5 <- ures$Gamma
            tp[i] <- sum(PM[, i] == 1 & Rho5[, i] >= 1/8)
            fp[i] <- sum(PM[, i] == 0 & Rho5[, i] >= 1/8)
            fn[i] <- sum(PM[, i] == 1 & Rho5[, i] < 1/8)
            tn[i] <- sum(PM[, i] == 0 & Rho5[, i] < 1/8)
            prauc[i] <- aucs(Rho5[, i],PM[, i])$auc
        }
        safeSave(prauc, tp, fp, fn, tn, file = paste0(path, "nempi_bs/mf_loo.rda"))
    } else if (type %in% "knn") {
        for (i in startrun:(startrun+bsruns-1)) {
            if (i > ncol(D5)) { break() }
            D6 <- D5
            colnames(D6)[i] <- ""
            Rho2 <- PM
            Rho2[, i] <- 0
            Rho3 <- PM
            Rho3[, i] <- 0
            library(class)
            train <- t(D6[, which(colnames(D6) != "")])
            test <- t(D6[, which(colnames(D6) == "")])
            knn0 <- 0
            cl <- colnames(D6)[which(colnames(D6) != "")]
            knnres <- knn(train, test, cl, prob = TRUE)
            D3 <- D6
            colnames(D3)[which(colnames(D3) %in% "")] <- as.character(knnres)
            tmp <- nem(D3)
            ures <- list()
            ures$Gamma <- apply(getGamma(D3), 2, function(x) return(x/sum(x)))
            ures$res <- list()
            ures$res$adj <- tmp$adj
            ures$null <- TRUE
            ures$combi <- 1
            Rho5 <- ures$Gamma
            tp[i] <- sum(PM[, i] == 1 & Rho5[, i] >= 1/8)
            fp[i] <- sum(PM[, i] == 0 & Rho5[, i] >= 1/8)
            fn[i] <- sum(PM[, i] == 1 & Rho5[, i] < 1/8)
            tn[i] <- sum(PM[, i] == 0 & Rho5[, i] < 1/8)
            prauc[i] <- aucs(Rho5[, i],PM[, i])$auc
        }
        safeSave(prauc, tp, fp, fn, tn, file = paste0(path, "nempi_bs/knn_loo.rda"))
    } else if (type %in% "rand") {
        for (i in startrun:(startrun+bsruns-1)) {
            if (i > ncol(D5)) { break() }
            Rho5 <- PM
            Rho5[,i] <- runif(nrow(Rho5))
            tp[i] <- sum(PM[, i] == 1 & Rho5[, i] >= 1/8)
            fp[i] <- sum(PM[, i] == 0 & Rho5[, i] >= 1/8)
            fn[i] <- sum(PM[, i] == 1 & Rho5[, i] < 1/8)
            tn[i] <- sum(PM[, i] == 0 & Rho5[, i] < 1/8)
            prauc[i] <- aucs(Rho5[, i],PM[, i])$auc
        }
        safeSave(prauc, tp, fp, fn, tn, file = paste0(path, "nempi_bs/rand_loo.rda"))
    }
    stop("successfull loo")
} else {

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
                                 sample.type = "Primary solid Tumor",
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
                    sampletype <- c("Primary Tumor", "Solid Tissue Normal")
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

    stop("tcga done")

}

## commands for hpc:

system("scp ~/Documents/mnem/src/mm.cpp euler:")
system("scp ~/Documents/mnem/R/mnems.r euler:")
system("scp ~/Documents/mnem/R/mnems_low.r euler:")
system("scp ~/Documents/nempi/R/nempi_main.r euler:")
system("scp ~/Documents/nempi/R/nempi_low.r euler:")

system("scp testing/nempi/other/nempi_loo.r euler:")

do="loo"

rm error.txt
rm output.txt

ram=8000 # nempi 8000, rf 64000, svm 20000, mf 8000, knn 8000, nn 20000
queue=4 # nempi 2min, rf 90min, svm 10min, mf 2min, knn 1min, nn 40min
bsruns=10
startrun=1
R="R/bin/R"

type="nn" # nempi, svm, rf, nn, knn, mf

if [ ${type} == 'rf' ] || [ ${type} == 'svm' ] || [ ${type} == 'nn' ]
then
maxpp=500000
else
maxpp=50000
fi
if [ ${type} == 'rf' ]
then
ram=64000
bsruns=1
fi
if [ ${type} == 'svm' ]
then
ram=20000
fi
if [ ${type} == 'nn' ]
then
ram=20000
bsruns=1
fi

bsub -M ${ram} -q normal.${queue}h -n 1 -e error.txt -o output.txt -R "rusage[mem=${ram}]" "${R} --vanilla --max-ppsize ${maxpp} --silent --no-save --args '${do}' '${type}' '${bsruns}' '${startrun}' < nempi_loo.r"

end=$(expr 91 + $bsruns )
end=$(expr $end / $bsruns )
## end=85

for (( i=2; i<=$end; i++ )); do
a=$(expr ${i} - 1 )
b=$(expr ${a} \* ${bsruns} )
startrun=$(expr ${b} + 1 )
echo $startrun
bsub -M ${ram} -q normal.${queue}h -n 1 -e error.txt -o output.txt -R "rusage[mem=${ram}]" "${R} --vanilla --max-ppsize ${maxpp} --silent --no-save --args '${do}' '${type}' '${bsruns}' '${startrun}' < nempi_loo.r"
done

for i in {2..10}; do
bsub -M ${ram} -q normal.${queue}h -n 1 -e error.txt -o output.txt -R "rusage[mem=${ram}]" "${R} --vanilla --max-ppsize ${maxpp} --silent --no-save --args '${do}' '${type}' '${bsruns}' '${startrun}' < nempi_loo.r"
done

## read in loo results:

methods <- c("nempi", "svm", "nn", "rf", "mf", "knn", "rand")
methlist <- list()
comparePre <- compareRec <- comparePrauc <- NULL
for (method in methods) {
    print(method)
    methlist[[method]] <- list()
    tp2 <- tn2 <- fp2 <- fn2 <- prauc2 <- numeric(91)
    for (file in list.files("~/Mount/Eulershare/nempi_bs")) {
        if (length(grep(paste0("^", method, "_"), file)) > 0 &
            length(grep(paste0("_loo\\."), file)) > 0) {
            load(paste0("~/Mount/Eulershare/nempi_bs/", file))
            tp2 <- tp2 + tp
            tn2 <- tn2 + tn
            fp2 <- fp2 + fp
            fn2 <- fn2 + fn
            prauc2 <- prauc2 + prauc
        }
    }
    methlist[[method]][["tp"]] <- tp2
    methlist[[method]][["tn"]] <- tn2
    methlist[[method]][["fp"]] <- fp2
    methlist[[method]][["fn"]] <- fn2
    methlist[[method]][["pre"]] <- tp2/(tp2+fp2)
    methlist[[method]][["rec"]] <- tp2/(tp2+fn2)
    methlist[[method]][["rec"]][which(is.na(methlist[[method]][["rec"]]))] <- 0
    methlist[[method]][["prauc"]] <- prauc2
    comparePre <- cbind(comparePre, methlist[[method]][["pre"]])
    compareRec <- cbind(compareRec, methlist[[method]][["rec"]])
    comparePrauc <- cbind(comparePrauc, methlist[[method]][["prauc"]])
}

setEPS()
postscript("temp.eps", height = 7, width = 10)
layout.mat <- matrix(c(rep(c(1,4,7),each=2),10,rep(c(2,5,8),each=2),11,rep(c(3,6,9),each=2),12),7)
par(mfrow=c(1,2))
col <- c("red", "blue", "darkgreen", "brown", "orange", "turquoise","grey")
mnem:::myboxplot(comparePre, dens = 0, col = col, xaxt = "n", main = "precision of mutation calls", ylab = "precision",box=1,border = col,medcol="black")
axis(1, 1:7, round(apply(comparePre,2,mean), 2))
mnem:::myboxplot(compareRec, dens = 0, col = col, xaxt = "n", main = "recall of mutation calls", ylab = "recall",box=1,border = col,medcol="black")
axis(1, 1:7, round(apply(compareRec,2,mean), 2))
dev.off()


setEPS()
postscript("temp.eps", height = 5, width = 5)
layout.mat <- matrix(c(rep(rep(1,each=3),3),2:4),4,byrow=1)
layout(layout.mat)
wtest <- numeric(ncol(comparePrauc))
for (i in 1:length(wtest)) {
    wtest[i] <- wilcox.test(comparePrauc[,1],comparePrauc[,i],alternative="greater")$p.value
}
col <- c("red", "blue", "darkgreen", "brown", "orange", "turquoise","grey")
mnem:::myboxplot(comparePrauc, dens = 0, col = col, xaxt = "n", main = "area under the precision-recall curve of mutation calls", ylab = "area under the precision-recall curve",box=1,border = col,medcol="black")
axis(1, 1:7, round(apply(comparePrauc,2,mean), 2))
par(mar=rep(0, 4))
plot.new()
legend("topleft",legend=c(expression(NEM~pi), "svm", "neural net"),col=c("red", "blue", "darkgreen"),fill=c("red", "blue", "darkgreen"),box.col = "transparent",cex = cex.leg,ncol = 1,border="transparent")
plot.new()
legend("topleft",legend=c("random forest", "missForest", "knn"),col=c("brown", "orange", "turquoise"),fill=c("brown", "orange", "turquoise"),box.col = "transparent",cex = cex.leg,ncol = 1,border="transparent")
plot.new()
legend("topleft",legend=c("random"),col=c("grey"),fill=c("grey"),box.col = "transparent",cex = cex.leg,ncol = 1,border="transparent")
dev.off()
