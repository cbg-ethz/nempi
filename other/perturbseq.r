data <- as.character(commandArgs(TRUE)[1])
set <- as.numeric(commandArgs(TRUE)[2])
path <- as.character(commandArgs(TRUE)[3])
cores <- as.numeric(commandArgs(TRUE)[4])
subset <- as.numeric(commandArgs(TRUE)[5])

## path <- "~/Mount/Eulershare/perturbseq2/"

print(data)
print(set)
print(path)
print(cores)

library(Rcpp)
library(naturalsort)

source("mnems.r")
source("mnems_low.r")
sourceCpp(code=readChar("mm.cpp", file.info("mm.cpp")$size))
devtools::load_all("CausalPathways/dce/")

## source("~/Documents/mnem/R/mnems.r"); source("~/Documents/mnem/R/mnems_low.r"); sourceCpp("~/Documents/mnem/src/mm.cpp"); library(Rgraphviz);

if (data == "csv") {
    M <- readRDS(paste0(path, "Counts.rds"))
    for (i in unique(colnames(M))) {
        write.csv(M[,which(colnames(M)==i)],file=paste0(path, "Counts_", i, ".csv"))
    }
} else if (data == "mnem") {
    D <- mnem::createApp(sets = set, types = "data", n = 100)
    saveRDS(D, file = paste0("Perturbseq_data_", set, ".rds"))
} else if (data == "perturbseq2") {
    M <- list()
    exp.id <- list(e1 = c("75","01"), e2 = c("77","05"), e3 = c("81","10"))
    for (i in 1:3) {
        id1 <- exp.id[[i]][1]
        id2 <- exp.id[[i]][2]
        M[[i]] <- Matrix::readMM(paste0(path, "GSM24066", id1, "_10X0", id2, "_matrix.mtx.txt"))
        gc()
        M[[i]] <- as(M[[i]], "matrix")
        gc()
        genes <- read.table(paste0(path, "GSM24066", id1, "_10X0", id2, "_genes.tsv"))
        cells <- read.csv(paste0(path, "GSM24066", id1, "_10X0", id2, "_cell_identities.csv"))
        barcodes <- read.table(paste0(path, "GSM24066", id1, "_10X0", id2, "_barcodes.tsv"))
        rownames(M[[i]]) <- genes[,2]
        colnames(M[[i]]) <- barcodes[,1]
        cells <- cells[which(cells$good.coverage == TRUE | cells$good.coverage == "True"),]
        M[[i]] <- M[[i]][,which(colnames(M[[i]]) %in% cells$cell.BC)]
        colnames(M[[i]]) <- cells$guide.identity[match(colnames(M[[i]]), cells$cell.BC)]
        if (i == 1) {
            colnames(M[[i]]) <- gsub("_pD.*|_pM.*|_pB.*|_only.*", "", colnames(M[[i]]))
            colnames(M[[i]]) <- gsub("\\*.*", "Unknwn", colnames(M[[i]]))
            colnames(M[[i]])[grep("mod",colnames(M[[i]]))] <- "Ctrl_Pilot"
        } else if (i == 2) {
            treatments <- gsub(".*-","",barcodes[,1])
            treatments <- treatments[which(barcodes[,1] %in% cells$cell.BC)]
            colnames(M[[i]]) <- gsub("_pD.*|_pM.*|_pB.*|_only.*", "", colnames(M[[i]]))
            colnames(M[[i]]) <- gsub("\\*.*", "Unknwn", colnames(M[[i]]))
            colnames(M[[i]]) <- gsub("3x_neg_ctrl","Ctrl_Epi",colnames(M[[i]]))
            colnames(M[[i]]) <- paste0(colnames(M[[i]]), "_", treatments)
        } else if (i == 3) {
            treatments <- gsub(".*-","",barcodes[,1])
            treatments <- treatments[which(barcodes[,1] %in% cells$cell.BC)]
            colnames(M[[i]]) <- gsub("_pD.*|_pM.*|_pB.*|_only.*", "", colnames(M[[i]]))
            colnames(M[[i]]) <- gsub("\\*.*", "Unknwn", colnames(M[[i]]))
            colnames(M[[i]])[grep("mod",colnames(M[[i]]))] <- "Ctrl_Main"
            colnames(M[[i]]) <- paste0(colnames(M[[i]]), "_", treatments)
        }
        gc()
    }
    saveRDS(M[[1]], file = paste0(path, "Counts_Pilot.rds"))
    saveRDS(M[[2]], file = paste0(path, "Counts_Epistasis.rds"))
    saveRDS(M[[3]], file = paste0(path, "Counts_Main.rds"))
    genes <- rownames(M[[1]])
    for (i in 2:3) {
        genes <- intersect(genes, rownames(M[[i]]))
    }
    for (i in 1:3) {
        M[[i]] <- M[[i]][genes,]
    }
    M <- do.call("cbind", M)
    gc()
    saveRDS(M, file = paste0(path, "Counts.rds"))
} else if (data == "edger") {
    print("start")
    for (i in subset) {
        print(i)
        if (i == 4) {
            M <- readRDS(paste0(path, "Counts_Main.rds"))
            genes.med <- apply(M, 1, median)
            M <- M[which(genes.med > 0),]
            treatments <- gsub(".*_","",colnames(M))
            colnames(M) <- gsub("_.*","",colnames(M))
            colnames(M) <- gsub("Ctrl","0_Ctrl",colnames(M))
            if (set != 2) {
                design <- model.matrix(~colnames(M)+factor(treatments))
                colnames(design) <- gsub("colnames\\(M\\)", "", colnames(design))
                genes <- unique(colnames(M))
                genes <- genes[-grep("Ctrl",genes)]
                library(edgeR)
                for (j in genes) {
                    print(j)
                    cells <- grep(paste0("Ctrl|",j),colnames(M))
                    Mtmp <- M[,cells]
                    designTmp <- design[cells,]
                    designTmp <- designTmp[,which(apply(designTmp,2,sum)>0)]
                    y <- DGEList(counts=Mtmp)
                    y <- calcNormFactors(y)
                    y <- estimateDisp(y,design=designTmp)
                    fit <- glmQLFit(y,design=designTmp)
                    P <- matrix(NA, nrow(Mtmp), 1)
                    colnames(P) <- j
                    rownames(P) <- rownames(M)
                    F <- P
                    qlf <- glmQLFTest(fit,coef=2)
                    P[,1] <- qlf$table$PValue
                    F[,1] <- qlf$table$logFC
                    saveRDS(P, file = paste0(path, "P_edgeR_", i, "_", j, ".rds"))
                    saveRDS(F, file = paste0(path, "F_edgeR_", i, "_", j, ".rds"))
                }
            }
        } else if (i %in% 1:3) {
            M <- readRDS(paste0(path, "Counts_Epistasis.rds"))
            M <- M[,grep(paste0("_",i),colnames(M))]
            genes.med <- apply(M, 1, median)
            M <- M[which(genes.med > 0),]
            colnames(M) <- gsub("Ctrl","0_Ctrl",colnames(M))
            design <- model.matrix(~colnames(M))
            colnames(design) <- gsub("colnames\\(M\\)", "", colnames(design))
        } else if (i == 5) {
            M <- readRDS(paste0(path, "Counts_Pilot.rds"))
            genes.med <- apply(M, 1, median)
            M <- M[which(genes.med > 0),]
            colnames(M) <- gsub("Ctrl","0_Ctrl",colnames(M))
            design <- model.matrix(~colnames(M))
            colnames(design) <- gsub("colnames\\(M\\)", "", colnames(design))
        }
        gc()
        if (set == 1 & i %in% c(1:3,5)) {
            print("bulk")
            library(edgeR)
            y <- DGEList(counts=M)
            y <- calcNormFactors(y)
            y <- estimateDisp(y,design=design)
            fit <- glmQLFit(y,design=design)
            P <- matrix(NA, nrow(M), ncol(design)-1)
            colnames(P) <- colnames(design)[-1]
            rownames(P) <- rownames(M)
            F <- P
            for (j in 1:ncol(P)) {
                print(colnames(P)[j])
                qlf <- glmQLFTest(fit,coef=j+1)
                P[,j] <- qlf$table$PValue
                F[,j] <- qlf$table$logFC
            }
            saveRDS(P, file = paste0(path, "P_edgeR_", i, ".rds"))
            saveRDS(F, file = paste0(path, "F_edgeR_", i, ".rds"))
        } else if (set == 2 & i %in% c(1:5)) {
            print("single cell")
            colnames(M) <- gsub("_Epi.*|_Pilot.*|_Main.*","",colnames(M))
            colnames(M)[which(!(colnames(M) == "0_Ctrl"))] <- paste0(colnames(M)[which(!(colnames(M) == "0_Ctrl"))], "_", 1:sum(!(colnames(M) == "0_Ctrl")))
            design <- model.matrix(~colnames(M))
            colnames(design) <- gsub("colnames\\(M\\)", "", colnames(design))
            library(edgeR)
            doCell <- function(x) {
                designTmp <- design[,c(1,x+1)]
                idx <- c(which(apply(design,1,sum)==1),which(designTmp[,-1]==1))
                Mtmp <- M[,idx]
                designTmp <- designTmp[idx,]
                y <- DGEList(counts=Mtmp)
                y <- calcNormFactors(y)
                y <- estimateDisp(y,design=designTmp)
                fit <- glmFit(y, design=designTmp)
                qlf <- glmLRT(fit,coef=2)
                P <- qlf$table$PValue
                F <- qlf$table$logFC
                return(cbind(P,F))
            }
            if (cores > 1) {
                library(snowfall)
                sfInit(parallel = TRUE, cpus = cores)
                sfExport("design","M")
                sfLibrary(edgeR)
                tmp <- sfLapply(1:(ncol(design)-1),doCell)
                sfStop()
            } else {
                tmp <- lapply(1:(ncol(design)-1),doCell)
            }
            PF <- do.call("cbind",tmp)
            F <- PF[,grep("F",colnames(PF))]
            P <- PF[,grep("P",colnames(PF))]
            colnames(P) <- colnames(F) <- colnames(design)[-1]
            saveRDS(P, file = paste0(path, "P_edgeR_", i, "_sc.rds"))
            saveRDS(F, file = paste0(path, "F_edgeR_", i, "_sc.rds"))
        }
    }
} else if (data == "pca") {
    if (subset %in% c(1,2,3)) {
        M <- readRDS(paste0(path, "Counts_Epistasis.rds"))
        M <- M[,grep(paste0("_",i),colnames(M))]
    } else if (subset == 4) {
        M <- readRDS(paste0(path, "Counts_Main.rds"))
    } else if (subset == 5) {
        M <- readRDS(paste0(path, "Counts_Pilot.rds"))
    }
    genes.med <- apply(M, 1, median)
    M <- M[which(genes.med > 0),]
    pca <- prcomp(t(log(M+1)))
    col <- unique(colnames(M))
    col <- match(colnames(M),col)
    varpct <- pca$sdev^2
    varpct <- round(varpct/sum(varpct),2)
    comps <- c(1,2)
    pdf(paste0(path, "pca_", subset, ".pdf"), width = 20, height = 20)
    par(mfrow=c(1,2))
    plot(pca$x[, comps], col = col, xlab = varpct[comps[1]], ylab = varpct[comps[2]], main = "log(x+1) transformed raw data")
    plot(sort(unique(col)),col=sort(unique(col)))
    legend("bottomright",legend=unique(colnames(M)),fill=sort(unique(col)),col=sort(unique(col)))
    dev.off()
} else if (data == "lods") {
    print("start")
    for (i in subset) {
        print(i)
        if (i == 4) {
            M <- readRDS(paste0(path, "Counts_Main.rds"))
            genes.med <- apply(M, 1, median)
            M <- M[which(genes.med > 0),]
            colnames(M) <- gsub("_.*","",colnames(M))
        } else if (i %in% 1:3) {
            M <- readRDS(paste0(path, "Counts_Epistasis.rds"))
            M <- M[,grep(paste0("_",i),colnames(M))]
            genes.med <- apply(M, 1, median)
            M <- M[which(genes.med > 0),]
        } else if (i == 5) {
            M <- readRDS(paste0(path, "Counts_Pilot.rds"))
            genes.med <- apply(M, 1, median)
            M <- M[which(genes.med > 0),]
        }
        M <- Linnorm::Linnorm(M,DataImputation=FALSE)
        print("log odds")
        colnames(M) <- gsub("_Epi.*|_Pilot.*|_Main.*","",colnames(M))
        ctrls <- which(colnames(M)=="Ctrl")
        doGene <- function(x,ctrls=ctrls) {
            ctrlcde <- ks::kcde(M[x,ctrls])
            kdidx <- which(!(colnames(M) == "Ctrl"))
            Z0 <- length(ctrlcde$eval.points)
            L <- numeric(length(kdidx))
            for (i in 1:length(kdidx)) {
                kd <- colnames(M)[kdidx[i]]
                kds <- which(colnames(M)==kd)
                if (all(M[x,kds]==0)) {
                    L[i] <- 0
                    next()
                }
                kdcde <- ks::kcde(M[x,kds])
                Z1 <- length(kdcde$eval.points)
                pkd <- min(sum(kdcde$eval.points >= M[x,i])/Z1,sum(kdcde$eval.points <= M[x,i])/Z1)
                pctrl <- min(sum(ctrlcde$eval.points >= M[x,i])/Z0,sum(ctrlcde$eval.points <= M[x,i])/Z0)
                if (pkd == 0) {
                    pkd <- min(0.01,pctrl)
                }
                if (pctrl == 0) {
                    pctrl <- min(0.01,pkd)
                }
                if (pctrl == 0 & pkd == 0) {
                    L[i] <- 0
                } else {
                    L[i] <- log(pkd/pctrl)
                }
            }
            print(x)
            return(L)
        }
        if (cores > 1) {
            library(snowfall)
            sfInit(parallel = TRUE, cpus = cores)
            sfExport("ctrls","M")
            sfLibrary(ks)
            tmp <- sfLapply(1:nrow(M),doGene,ctrls)
            sfStop()
        } else {
            tmp <- lapply(1:nrow(M),doGene,ctrls)
        }
        L <- do.call("rbind",tmp)
        colnames(L) <- colnames(M)[which(!(colnames(M) == "Ctrl"))]
        rownames(L) <- rownames(M)
        saveRDS(L, file = paste0(path, "L_linnorm_", i, "_sc.rds"))
    }
}
stop("success")

## hpc commands:

system("scp ~/Documents/mnem/R/* euler.ethz.ch:")
system("scp ~/Documents/testing/nest/nest_app.r euler.ethz.ch:")

rm error.txt
rm output.txt

ram=16000 # data: 64gb, wilcox sc: 8gb for 32 nodes, edger: 64gb, voom: 64gb
queue=24
cores=8

data="lods"
set=2
subset=5
path=/cluster/work/bewi/members/mpirkl/perturbseq2/

bsub -M ${ram} -q normal.${queue}h -n ${cores} -e error.txt -o output.txt -R "rusage[mem=${ram}]" "R/bin/R --vanilla --silent --no-save --args '${data}' '${set}' '${path}' '${cores}' '${subset}' < nest_app.r"

## temp:

rm error.txt
rm output.txt

queue=24
ram=64000
cores=4

bsub -M ${ram} -q normal.${queue}h -n ${cores} -e error.txt -o output.txt -R "rusage[mem=${ram}]" "R/bin/Rscript -e 'library(DawnRank); load(\"DawnRank_temp.rda\"); dawnRankScore <- DawnRank(adjMatrix = pathway2, mutationMatrix = M, expressionMatrix = normalizedDawn, mu = 3, parallel = ${cores}); saveRDS(dawnRankScore,file=\"nempi_dawnrank2.rds\")'"

## escape quotes like \"

rm error.txt
rm output.txt

queue=24
ram=4000
cores=1

Rcmd=/cluster/work/bewi/members/mpirkl/miniconda3/bin/Rscript
## Rcmd=~/R/bin/Rscript

variable=latent
values=5
reps=1

bsub -M ${ram} -q normal.${queue}h -n ${cores} -e error.txt -o output.txt -R "rusage[mem=${ram}]" "${Rcmd} main.R --variable ${variable} --values ${values} --append FALSE --replicates ${reps} --link identity --output test.csv"

for i in {1..100}
do
bsub -M ${ram} -q normal.${queue}h -n ${cores} -e error.txt -o output.txt -R "rusage[mem=${ram}]" "${Rcmd} main.R --variable ${variable} --values ${values} --append TRUE --replicates ${reps} --link identity --output test.csv"
done


bsub -M ${ram} -q normal.${queue}h -n ${cores} -e error.txt -o output.txt -R "rusage[mem=${ram}]" "${Rcmd} main.R --variable node.num --values 100 --methods \"dce,dce.log,dce.tpm,dce.tpmlog\""

for i in {1..99}
do
bsub -M ${ram} -q normal.${queue}h -n ${cores} -e error.txt -o output.txt -R "rusage[mem=${ram}]" "${Rcmd} main.R --variable node.num --values 100 --methods \"dce,dce.log,dce.tpm,dce.tpmlog\""
done

## results:

## load functions:

library(Rcpp)
library(naturalsort)
source("~/Documents/mnem/R/mnems.r"); source("~/Documents/mnem/R/mnems_low.r"); sourceCpp("~/Documents/mnem/src/mm.cpp"); library(Rgraphviz);

##

for (i in 1:5) {
    L <- readRDS(paste0("~/Mount/Eulershare/perturbseq2/L_linnorm_", i, "_sc.rds"))
    print(dim(L))
}

## edger results:

P <- readRDS("~/Mount/Eulershare/perturbseq2/P_edgeR.rds")
F <- readRDS("~/Mount/Eulershare/perturbseq2/F_edgeR.rds")

##

P <- V$p.value[,-1]
F <- V$coefficients[,-1]

##

discFDR <- function(P) {
    Q <- P*0
    Q[1:length(Q)] <- p.adjust(P)
    E <- P*0
    E[which(Q < 0.1)] <- 1
    return(E)
}

llrFc <- function(F) {
    E <- abs(F)-mean(abs(F))
    return(E)
}

##

D <- E
colnames(D) <- gsub("_[1-3]","",colnames(D))
D <- D[,-grep("_",colnames(D))]
D <- D[,which(!(colnames(D)=="Unknwn"))]
colnames(D)[which(colnames(D) == "Gal4-4(mod)")] <- "GAL4"
Sgenes <- getSgenes(D)

res <- mnem(D,search="estimate",k=1,complete=TRUE,starts=10)
res$ll
pdf("temp.pdf", width = 8, height = 7)
plot(res, cells = 0, bestCell = 0, egenes = 0, showdata = TRUE)
dev.off()

##

phi <- res$comp[[1]]$phi*0

for (i in 1:1000) {
    Ds <- D[sample(1:nrow(D), nrow(D), replace = TRUE),]
    ress <- nemEst(Ds,method=method)
    phi <- phi + mnem::transClose(ress$phi)
    ## cat(paste0(i,"."))
}

phi2 <- phi/1010
phi2[which(phi2 < 0.99)] <- 0

pdf("temp.pdf", width = 20, height = 40)
mnem::plotDnf(transitive.reduction(phi2))
dev.off()

## bum model:

hist(P,freq=FALSE)
lines(density(punif(seq(0,1,length.out=1000))),col=2)
lines(density(pbeta(seq(0,1,length.out=1000),10,1)),col=2)

lods <- as.vector(P)
names(lods) <- rownames(P)
bum <- BioNet::bumOptim(lods)
dbeta <- dbeta(lods, bum$a, 1)
lods <- log(dbeta)
lods <- matrix(lods, nrow = nrow(P))
rownames(lods) <- rownames(P)
colnames(lods) <- colnames(P)

lods <- dbeta(P, 0.1, 1)
E <- lods
method <- "llr"

## gene length:

counts.csv <- read.csv("~/Mount/Eulershare/perturbseq2/Counts_COPB1.csv")

library("biomaRt")
mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
gene.info <- getBM(filters = "hgnc_symbol",
                     attributes = c("entrezgene_id", "hgnc_symbol","transcript_length"),
                     values = counts.csv[,1], mart = mart)

D <- as(tmp[,-1],"matrix")
rownames(D) <- tmp[,1]
rownames(D) <- gsub("-.*","",gsub("\\..*","",rownames(D)))
D <- rowsum(D, rownames(D))

gene.info <- getBM(filters = "hgnc_symbol",
                     attributes = c("entrezgene_id", "hgnc_symbol","transcript_length"),
                     values = rownames(D), mart = mart)

## logical model:

t <- 2
P <- readRDS(paste0("~/Mount/Eulershare/perturbseq2/P_edgeR_",t,"_.rds"))
F <- readRDS(paste0("~/Mount/Eulershare/perturbseq2/F_edgeR_",t,"_.rds"))
rownames(P) <- rownames(F) <- 1:nrow(F)

source("~/Documents/testing/booltest/booltest.r")

E <- discFDR(P)
method <- "disc"

D <- E
colnames(D) <- gsub("_[1-3]","",colnames(D))
D <- D[,-grep("_",colnames(D))]
D <- D[,which(!(colnames(D)=="Unknwn"))]
colnames(D)[which(colnames(D) == "Gal4-4(mod)")] <- "GAL4"
Sgenes <- getSgenes(D)

res <- mnem(D,search="exhaustive",k=1,complete=TRUE,starts=10)
res$ll
pdf("temp.pdf", width = 8, height = 7)
plot(res, cells = 0, bestCell = 0, egenes = 0, showdata = TRUE)
dev.off()

D <- E
colnames(D) <- gsub("_[1-3]","",colnames(D))
D <- D[,which(!(colnames(D)=="Unknwn"))]
if (method == "disc") {
    D <- D[which(apply(D,1,sum) > 0),]
} else {
    D <- D[which(apply(D,1,max) > 0),]
}
Sgenes <- getSgenes(D)

pdf("temp.pdf",width=5,height=40)
epiNEM::HeatmapOP(D,bordercol="transparent",aspect="iso")
dev.off()

lres <- ltest(D,method="exhaustive",search=search)
lres$pvals

for (i in 1:4) {
    print(lres$nemC[[i]]$adj)
    print(lres$nems[[i]]$adj)
}

for (i in 1:4) {
    lres$nemC[[i]]$adj <- transClose(lres$nemC[[i]]$adj)
    lres$nems[[i]]$adj <- transClose(lres$nems[[i]]$adj)
}

for (bs in 1:10) {
    lresbs <- ltest(D[sample(1:nrow(D), nrow(D), replace = TRUE), ],method = method,search=search)
    for (i in 1:4) {
        lres$nemC[[i]]$adj <- lres$nemC[[i]]$adj + transClose(lresbs$nemC[[i]]$adj)
        lres$nems[[i]]$adj <- lres$nems[[i]]$adj + transClose(lresbs$nems[[i]]$adj)
        lres$pvals <- lres$pvals + lresbs$pvals
    }
}

for (i in 1:4) {
    print(lres$nemC[[i]]$adj)
    print(lres$nems[[i]]$adj)
}

