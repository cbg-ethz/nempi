
method <- as.character(commandArgs(TRUE)[1])

type <- as.character(commandArgs(TRUE)[2])

goi <- as.numeric(commandArgs(TRUE)[3])

goitype <- as.numeric(commandArgs(TRUE)[4])

path <- "/cluster/work/bewi/members/mpirkl/"

## method <- "nempi"; type <- "TCGA-BRCA"; goi <- 0; goitype <- 1; path <- "~/Mount/Eulershare/"

library(naturalsort)
library(cluster)
library(Rcpp)
library(Rgraphviz)
library(mnem)
library(randomForest)
library(e1071)
library(missForest)
library(class)
library(nnet)

source("mnems.r")
source("mnems_low.r")
sourceCpp(code=readChar("mm.cpp", file.info("mm.cpp")$size))
source("nempi_main.r")
source("nempi_low.r")

## sourceCpp("~/Documents/mnem/src/mm.cpp"); source("~/Documents/mnem/R/mnems.r"); source("~/Documents/mnem/R/mnems_low.r"); source("~/Documents/nempi/R/nempi_main.r"); source("~/Documents/nempi/R/nempi_low.r")
## load("mutclust/TCGA-BRCA_final.rda")

type <- "TCGA-BRCA"
load(paste0(type, "_final.rda"))

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
Ps <- list()
Ps$all <- P2
Ps$meth <- Pmeth
Ps$cnv <- Pcnv
Ps$mut <- Pmut

if (method=="auc") {
    ## goitype <- 4
    methods <- c("nempi", "svm", "nn", "rf", "mf", "knn")
    mutinc <- 1
    files <- list.files(paste0(path,"nempi_tcga/"))
    if (goitype==1) {
        files <- files[grep("_pan",files)]
    } else if (goitype==2) {
        files <- files[grep("_all",files)]
    } else if (goitype==3) {
        files <- files[grep("_mut",files)]
    } else if (goitype==4) {
        files <- files[grep("_org",files)]
    }
    pr <- roc <- list()
    aucrows <- max(unlist(lapply(methods,function(x) {
        y <- length(grep(paste0("^",x),files))
        return(y)
    })))
    pr$all <- pr$mut <- pr$cnv <- pr$meth <-
        pr$allp <- pr$mutp <- pr$cnvp <- pr$methp <-
            pr$time <- matrix(NA,aucrows,length(methods)+2,
                              dimnames=list(rownames=NULL,colnames=c(methods,c("random bin","random cont"))))
    roc$all <- roc$mut <- roc$cnv <- roc$meth <-
        roc$allp <- roc$mutp <- roc$cnvp <- roc$methp <-
            roc$time <- matrix(NA,aucrows,length(methods)+2,
                              dimnames=list(rownames=NULL,colnames=c(methods,c("random bin","random cont"))))
    pr$missing <- roc$missing <- numeric(aucrows)
    for (method in methods) {
        ## method <- "nempi"
        print(method)
        files.method <- files[grep(paste0("^",method),files)]
        idx <- which(methods==method)
        idx2 <- 0
        for (file in files.method) {
            if (length(grep(paste0("^",method),file))==0) { next() }
            idx2 <- idx2 + 1
            cat(paste0(idx2,"."))
            ## file <- "nempi_1.rds"
            ures <- readRDS(paste0(path,"nempi_tcga/",file))
            pr$time[idx2,idx] <- ures$time
            goi <- rownames(ures$Gamma)
            P <- Pmut
            P <- P[which(rownames(P) %in% goi), ]
            P[which(P > 1)] <- 1
            ## P <- apply(P, 2, function(x) return(x/sum(x)))
            ## P[is.na(P)] <- 0
            Rho <- cbind(P, matrix(0, nrow(P), sum(!(colnames(D) %in% colnames(P)))))
            colnames(Rho) <- c(colnames(P), colnames(D)[which(!(colnames(D) %in% colnames(P)))])
            Rho <- Rho[, colnames(D)]
            if (sum(apply(Rho, 1, sum) == 0) > 0) {
                Rho <- Rho[-which(apply(Rho, 1, sum) == 0), ]
            }
            Rho[is.na(Rho)] <- 0
            if (method=="nempi") {
                pr$missing[idx2] <- roc$missing[idx2] <- sum(apply(Rho, 2, sum) == 0)/ncol(Rho) # unlabelled
            }
            inc <- sort(apply(Rho, 1, sum))
            D2 <- D[which(apply(D, 1, median) != 0), ]
            D3 <- D2[, which(duplicated(colnames(D2)) == FALSE)]
            D4 <- D3
            Rho <- Rho[, which(duplicated(colnames(D2)) == FALSE)]
            colnames(D4) <- apply(Rho, 2, function(x) {
                Sgenes <- paste(sort(rownames(Rho)[which(x > 0)]), collapse = "_")
                return(Sgenes)
            })
            for (mtype in c("all")) { # ,"mut","cnv","meth")) {
                ## mtype <- "all"
                P <- Ps[[mtype]] # P2
                if (!mutinc) { # include mutations or not (0)
                    P[which(Pmut == 1)] <- 0
                }
                P <- P[which(rownames(P) %in% rownames(Rho)), which(colnames(P) %in% colnames(Rho))]
                P <- P[order(rownames(P)), order(colnames(P))]
                Ptmp <- cbind(P, matrix(0, nrow(P), sum(!(colnames(Rho) %in% colnames(P)))))
                colnames(Ptmp) <- c(colnames(P), colnames(Rho)[which(!(colnames(Rho) %in% colnames(P)))])
                P <- Ptmp[, colnames(Rho)]
                P[which(P > 1)] <- 1
                meas <- list()
                meas$data <- P
                colnames(meas$data) <- apply(P, 2, function(x) {
                    y <- paste(sort(rownames(P)[which(x==1)]), collapse = "_")
                    return(y)
                })
                measnet <- matrix(0, nrow(P), nrow(P))
                diag(measnet) <- 1
                meas$Nem[[1]] <- measnet
                meas$Theta[[1]] <- ures$res$subtopo
                if (idx == 1) {
                    pitmp <- pifit(ures, meas, D4)
                    randres <- ures
                    randres$Gamma[1:length(randres$Gamma)] <- sample(c(0,1),length(randres$Gamma),replace=TRUE)
                    rand <- pifit(randres, meas, D4, propagate=FALSE)
                    pr[[mtype]][idx2,length(methods)+1] <- rand$auc
                    roc[[mtype]][idx2,length(methods)+1] <- rand$roc
                    randres2 <- randres
                    randres2$Gamma[1:length(randres2$Gamma)] <- runif(length(randres2$Gamma),0,1)
                    randres2$Gamma <- apply(randres2$Gamma, 2, function(x) {
                        y <- runif(1,0,1)
                        return(x/(sum(x)+y))
                    })
                    rand <- pifit(randres2, meas, D4, propagate=FALSE)
                    pr[[mtype]][idx2,length(methods)+2] <- rand$auc
                    roc[[mtype]][idx2,length(methods)+2] <- rand$roc
                } else {
                    pitmp <- pifit(ures, meas, D4, propagate=FALSE)
                }
                pr[[mtype]][idx2,idx] <- pitmp$auc
                roc[[mtype]][idx2,idx] <- pitmp$roc
                meas$Nem[[1]] <- ures$res$adj
                if (idx == 1) {
                    pitmp <- pifit(ures, meas, D4)
                    rand <- pifit(randres, meas, D4, propagate=FALSE)
                    pr[[paste0(mtype,"p")]][idx2,length(methods)+1] <- rand$auc
                    roc[[paste0(mtype,"p")]][idx2,length(methods)+1] <- rand$roc
                    rand <- pifit(randres2, meas, D4, propagate=FALSE)
                    pr[[paste0(mtype,"p")]][idx2,length(methods)+2] <- rand$auc
                    roc[[paste0(mtype,"p")]][idx2,length(methods)+2] <- rand$roc
                } else {
                    pitmp <- pifit(ures, meas, D4, propagate=FALSE)
                }
                pr[[paste0(mtype,"p")]][idx2,idx] <- pitmp$auc
                roc[[paste0(mtype,"p")]][idx2,idx] <- pitmp$roc
            }
        }
    }
    if (goitype==1) {
        saveRDS(pr,file=paste0(path,"nempi_tcga/pr_results_pan.rds"))
        saveRDS(roc,file=paste0(path,"nempi_tcga/roc_results_pan.rds"))
    } else if (goitype==2) {
        saveRDS(pr,file=paste0(path,"nempi_tcga/pr_results_all.rds"))
        saveRDS(roc,file=paste0(path,"nempi_tcga/roc_results_all.rds"))
    } else if (goitype==3) {
        saveRDS(pr,file=paste0(path,"nempi_tcga/pr_results_mut.rds"))
        saveRDS(roc,file=paste0(path,"nempi_tcga/roc_results_mut.rds"))
    } else if (goitype==4) {
        saveRDS(pr,file=paste0(path,"nempi_tcga/pr_results_org.rds"))
        saveRDS(roc,file=paste0(path,"nempi_tcga/roc_results_org.rds"))
    }
    stop("success")
}

## genes of interest

if (goi!=0) {
    ## goi <- 1
    goi0 <- goi
    gois <- rownames(Pmut)[which(apply(Pmut,1,sum)>0)]
    prob <- rep(1/length(gois),length(gois))
    if (goitype==1) {
        gois <- intersect(c("JAK1","FGFR1","MSH2","HIST1H1C","MGMT","SMARCA1","ZC3H12A","SETBP1","PIK3CB","KRT222","HGF","IL7R","DAZAP1","LEMD2","MACF1","MYC","THRAP3","RPL5","H3F3C","AXIN2","POLQ","AR","UNCX","PIK3CG","ABL1","PMS1","IRF6","ALK","ESR1","PGR","CHD8","H3F3A","PTPRC","CACNA1A","TBL1XR1","EGR3","JAK3","GRIN2D","POLE","MSH3","RARA","ATXN3","PMS2","ACVR1B","TLR4","PLXNB2","RNF111","CDKN2C","MLH1","BRCA2","ARAF","IRF2","CHEK2","USP9X","TRAF3","EPHA3","NOTCH2","CBWD3","PAX5","RHEB","KMT2A","BRAF1","BCL2","JAK2"),gois) # pan-cancer drivers: too mutated/cnved/methed ?? no, don't think so.
    } else if (goitype==3) {
        prob0 <- apply(Pmut[which(apply(Pmut,1,sum)>0),],1,function(x) return(sum(x!=0)))
        prob <- prob0
        prob[which(prob<20)] <- 0
        prob <- prob/sum(prob)
    }
    ## goi <- 1
    set.seed(goi*9247)
    goi <- sample(gois,10,prob=prob)
    print(goi)
    epiNEM::HeatmapOP(Pmut[goi,])
    apply(Pmut[goi,],1,sum)
} else {
    goi0 <- goi
    goi <- c("MAP2K4", "GATA3", "GPS2", "TBX3", "PTPRD", "NCOR1", "CBFB", "CDKN1B") # BRCA
}

P <- Pmut
P <- P[which(rownames(P) %in% goi), ]
P[which(P > 1)] <- 1
P <- apply(P, 2, function(x) return(x/sum(x)))
P[is.na(P)] <- 0

## data imputation:

Rho <- cbind(P, matrix(0, nrow(P), sum(!(colnames(D) %in% colnames(P)))))
colnames(Rho) <- c(colnames(P), colnames(D)[which(!(colnames(D) %in% colnames(P)))])
Rho <- Rho[, colnames(D)]
if (sum(apply(Rho, 1, sum) == 0) > 0) {
    Rho <- Rho[-which(apply(Rho, 1, sum) == 0), ]
}
Rho[is.na(Rho)] <- 0
print("unlabelled:")
print(sum(apply(Rho, 2, sum) == 0)/ncol(Rho))
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
D4 <- D2
colnames(D4) <- apply(Rho, 2, function(x) {
    Sgenes <- paste(sort(rownames(Rho)[which(x > 0)]), collapse = "_")
    return(Sgenes)
})

## pdf("temp.pdf", width = 12, height = 6);tmp <- Rho;colnames(tmp) <- NULL;epiNEM::HeatmapOP(tmp, col = "RdBu", Rowv = 0, bordercol = "transparent");dev.off()

start <- as.numeric(Sys.time())
if (method=="nempi") {
    converged <- 10
    nempires <- nempi(D2, Gamma = Rho, full = TRUE, converged = converged)
    ures <- nempires
} else if (method=="svm") {
    svmres <- classpi(D4, full = TRUE, method = "svm")
    ures <- svmres
} else if (method=="nn") {
    nnres <- classpi(D4, full = TRUE, method = "nnet", MaxNWts = 50000, size = 2)
    ures <- nnres
} else if (method=="mf") {
    mfdata <- cbind(as.data.frame(t(D4)), factor(colnames(D4)))
    mfdata[which(mfdata == "", arr.ind = TRUE)] <- NA
    mfimp <- missForest(mfdata)
    D4 <- D2
    colnames(D4) <- mfimp$ximp[, ncol(mfimp$ximp)]
    tmp <- nem(D4, multi = TRUE)
    Gamma <- getGamma(D4)
    ures <- list()
    ures$Gamma <- apply(Gamma, 2, function(x) return(x/sum(x)))
    ures$res <- list()
    ures$res$adj <- tmp$adj
    ures$null <- TRUE
    ures$combi <- 1
    mfres <- ures
} else if (method=="rf") {
    rfres <- classpi(D4, full = TRUE, method = "randomForest")
    ures <- rfres
} else if (method=="knn") {
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
    tmp <- nem(D3, multi = TRUE)
    Gamma <- getGamma(D3)
    ures <- list()
    ures$Gamma <- Gamma # apply(Gamma, 2, function(x) return(x/sum(x)))
    ures$res <- list()
    ures$res$adj <- tmp$adj
    ures$null <- TRUE
    ures$combi <- 1
    knnres <- ures
} else if (method=="mice") { # takes forever
    mfdata <- cbind(as.data.frame(t(D4)), factor(colnames(D4)))
    mfdata[which(mfdata == "", arr.ind = TRUE)] <- NA

    micedata <- mfdata
    colnames(micedata) <- paste0(LETTERS[1:ncol(micedata)], 1:ncol(micedata))
    miceres <- mice(micedata, method = c(rep('rfcont', ncol(micedata)-1), 'rfcat'), m = 2, maxit = 2)
    ures <- miceres
}
end <- as.numeric(Sys.time())
print(end - start)
ures$time <- end-start

if (goi0!=0) {
    if (goitype==1) {
        saveRDS(ures, file=paste0(path, "nempi_tcga/", method, "_", goi0, "_pan.rds"))
    } else if (goitype==2) {
        saveRDS(ures, file=paste0(path, "nempi_tcga/", method, "_", goi0, "_all.rds"))
    } else if (goitype==3) {
        saveRDS(ures, file=paste0(path, "nempi_tcga/", method, "_", goi0, "_mut.rds"))
    }
} else {
    saveRDS(ures, file=paste0(path, "nempi_tcga/", method, "_org.rds"))
}
stop("tcga done")

## commands for hpc:

system("scp ~/Documents/mnem/R/mnems.r euler.ethz.ch:")
system("scp ~/Documents/mnem/R/mnems_low.r euler.ethz.ch:")
system("scp ~/Documents/nempi/R/nempi_main.r euler.ethz.ch:")
system("scp ~/Documents/nempi/R/nempi_low.r euler.ethz.ch:")

system("scp testing/nempi/other/TCGA_nempi.r euler.ethz.ch:")

rm error.txt
rm output.txt
rm .RData

ram=8000 # 4h: nempi 8, mf 8, knn 8, svm 32, nn 32, rf 32, auc 16
queue=4
method="auc"
type="TCGA-BRCA"
goi=1
goitype=3

bsub -M ${ram} -q normal.${queue}h -e error.txt -o output.txt -R "rusage[mem=${ram}]" "R/bin/R --max-ppsize=500000 --vanilla --silent --no-save --args '${method}' '${type}' '${goi}' '${goitype}' < TCGA_nempi.r"

for goi in {2..100}; do
bsub -M ${ram} -q normal.${queue}h -e error.txt -o output.txt -R "rusage[mem=${ram}]" "R/bin/R --max-ppsize=500000 --vanilla --silent --no-save --args '${method}' '${type}' '${goi}' '${goitype}' < TCGA_nempi.r"
done

for i in {141716598..141716731}; do
bkill ${i}
done

## analysis:

method <- "nempi"
ures <- readRDS(paste0("~/Mount/Eulershare/nempi_tcga/",method,"_org.rds"))

sum(ures$lls[2:length(ures$lls)] - ures$lls[1:(length(ures$lls)-1)] < 0)

pdf("temp.pdf", width = 12, height = 6)
epiNEM::HeatmapOP(ures$Gamma, bordercol = rgb(0,0,0,0), col = "RdBu")
#plot(ures, edgewidth = 30)
dev.off()

pdf("temp.pdf", width = 9, height = 6)
par(mfrow=c(2,3))
plotConvergence(ures, type = "b", col = "blue")
dev.off()

pdf("temp.pdf", width = 6, height = 6)
Pgenes <- sort(unique(rownames(ures$Gamma)))
adjtmp <- ures$res$adj
colnames(adjtmp) <- rownames(adjtmp) <- Pgenes
plotDnf(adj2dnf(adjtmp), edgelwd = 2)
dev.off()

pdf("temp.pdf", width = 12, height = 6)
epiNEM::HeatmapOP(t(mytc(ures$res$adj))%*%ures$Gamma, bordercol = rgb(0,0,0,0), col = "RdBu")
#plot(ures, edgewidth = 30)
dev.off()

source("~/Documents/nempi/R/nempi_main.r")
D4 <- D2
colnames(D4) <- apply(Rho, 2, function(x) {
    Sgenes <- paste(sort(rownames(Rho)[which(x > 0)]), collapse = "_")
    return(Sgenes)
})

## new analysis:

goitype <- "pan"
stat <- readRDS(paste0("~/Mount/Eulershare/nempi_tcga/pr_results_",goitype,".rds"))

cex.main <- 1.5
cex.lab <- 1.5
idx <- c(1:6,8)
setEPS()
postscript("temp.eps", height = 5, width = 10)
par(mar=c(3.85,4,2.75,1),oma=c(0,0,0,0))
layout.mat <- matrix(c(rep(rep(1:3,each=1),3),4:6),4,byrow=1)
layout(layout.mat)
methods <- c("nempi", "svm", "nn", "rf", "mf", "knn")
wc.ps <- numeric((length(idx)-1))
wc.ps[1:length(wc.ps)] <- NA
for (i in 1:(length(idx)-1)) {
    wc.ps[i] <- wilcox.test(stat$all[,1],stat$all[,idx[1+i]],alternative="greater")$p.value
}
wc.psp <- numeric((length(idx)-1))
wc.psp[1:length(wc.psp)] <- NA
for (i in 1:(length(idx)-1)) {
    wc.psp[i] <- wilcox.test(stat$allp[,1],stat$allp[,idx[1+i]],alternative="greater")$p.value
}
col <- c("red", "blue", "darkgreen", "brown", "orange", "turquoise","grey","darkgrey")
ylim <- c(min(stat$all,na.rm=TRUE),1)
mnem:::myboxplot(stat$all[,idx],col=col,dens=0,xaxt = "n",main = "no perturbation propagation",ylab ="area under the precision-recall curve",box=1,border=col,medcol="black",ylim=ylim,cex.main=cex.main,cex.lab=cex.lab)
max.val <- apply(stat$all[,idx[-1]],2,max,na.rm=TRUE) + 0.05
if (goitype=="pan") {
    max.val <- max.val + rep(0.05,length(max.val))
}
text(2:(length(idx[-1])+1),max.val,labels=round(wc.ps,7),srt=90)
mnem:::myboxplot(stat$allp[,idx],col=col,dens=0,xaxt = "n",main = "perturbation propagation",ylab ="area under the precision-recall curve",box=1,border=col,medcol="black",ylim=ylim,cex.main=cex.main,cex.lab=cex.lab)
max.valp <- apply(stat$allp[,idx[-1]],2,min,na.rm=TRUE) - 0.05
if (goitype=="pan") {
    max.valp <- max.valp + c(0,0,0,-0.01,-0.02,-0.05)
}
##text(2:(length(idx[-1])+1),max.valp,labels=round(wc.psp,7))
text((length(idx[-1])+1),max.valp[length(max.valp)],labels=round(wc.psp[length(wc.psp)],7),srt=90)
mnem:::myboxplot(stat$time[,1:6],col=col,dens=0,xaxt = "n",main = "running time",ylab ="seconds",box=1,border=col,medcol="black",log="y",cex.main=cex.main,cex.lab=cex.lab)
##mnem:::myboxplot(stat$missing,col=col,dens=0,xaxt = "n",main = "fraction of unlabelled samples",ylab ="fraction",box=1,border=col,medcol="black",log="y")
cex.leg <- 1.5
par(mar=rep(0, 4))
plot.new()
legend("topleft",legend=c(expression(NEM~pi), "svm", "neural net"),col=c("red", "blue", "darkgreen"),fill=c("red", "blue", "darkgreen"),box.col = "transparent",cex = cex.leg,ncol = 1,border="transparent")
plot.new()
legend("topleft",legend=c("random forest", "missForest", "knn"),col=c("brown", "orange", "turquoise"),fill=c("brown", "orange", "turquoise"),box.col = "transparent",cex = cex.leg,ncol = 1,border="transparent")
plot.new()
legend("topleft",legend=c("random"),col=c("grey"),fill=c("grey"),box.col = "transparent",cex = cex.leg,ncol = 1,border="transparent")
dev.off()

## colnames
## rownames     nempi       svm        nn        rf        mf       knn random bin
##     [1,] 0.8830668 0.8682249 0.7571648 0.7754192 0.7362392 0.7826912  0.6653634
##         colnames
## rownames random cont
## [1,]   0.7385025

## colnames
## rownames     nempi       svm      nn        rf        mf       knn random bin
##     [1,] 0.6007286 0.6215111 0.57883 0.6321455 0.6049589 0.5877279  0.5210848
##         colnames
## rownames random cont
##     [1,]   0.5418864

## old analysis:

## save(nempires, knnres, rfres, mfres, svmres, nnres, Rho, D2, Pmut, Pmeth, Pcnv, file = paste0(path, type, "_nempi.rda"))

path <- "mutclust/"; type <- "TCGA-BRCA"

load(paste0(path, type, "_nempi.rda"))

source("~/Documents/nempi/R/nempi_main.r")
D4 <- D2
colnames(D4) <- apply(Rho, 2, function(x) {
    Sgenes <- paste(sort(rownames(Rho)[which(x > 0)]), collapse = "_")
    return(Sgenes)
})

ures <- nempires

## check against methylation and cnvs:

pdf("temp.pdf", width = 12, height = 6)
tmp <- ures$Gamma
colnames(tmp) <- NULL
epiNEM::HeatmapOP(tmp, bordercol = rgb(0,0,0,0), col = "RdBu", colorkey = NULL)
dev.off()

pdf("temp.pdf", width = 6, height = 6)
Pgenes <- sort(unique(colnames(D2)))
adjtmp <- ures$res$adj
colnames(adjtmp) <- rownames(adjtmp) <- Pgenes
plotDnf(adjtmp, edgelwd = 2)
dev.off()

## cnv/meth enrichment:

# load("~/Documents/nempi_backup/backup3/TCGA-BRCA_nempi.rda")

source("~/Documents/nempi/R/nempi_main.r")
source("~/Documents/nempi/R/nempi_low.r")

methods <- list("NEM$\\pi$" = nempires)#, knn = knnres, mf = mfres, nn = nnres, rf = rfres, svm = svmres)

mutinc <- 1

Lall <- Lcnv <- Lmeth <- Lmut <- Laucor <- list()

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
    Ptmp <- cbind(P, matrix(0, nrow(P), sum(!(colnames(Rho) %in% colnames(P)))))
    colnames(Ptmp) <- c(colnames(P), colnames(Rho)[which(!(colnames(Rho) %in% colnames(P)))])
    P <- Ptmp[, colnames(Rho)]
    P <- t(mytc(methods[[1]]$res$adj))%*%P
    P[which(P > 1)] <- 1
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
    F <- matrix(c(sum(pmeth >= 1 & P == 1), sum(pmeth >= 1 & P == 0), sum(pmeth <= 0 & P == 1), sum(pmeth <= 0 & P == 0)), 2)
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
    F <- matrix(c(sum(pmeth >= 1 & P == 1), sum(pmeth >= 1 & P == 0), sum(pmeth <= 0 & P == 1), sum(pmeth <= 0 & P == 0)), 2)
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
    F <- matrix(c(sum(pmeth >= 1 & P == 1), sum(pmeth >= 1 & P == 0), sum(pmeth <= 0 & P == 1), sum(pmeth <= 0 & P == 0)), 2)
    print(1 - phyper(F[1,1]-1, sum(F[, 1]), sum(F[, 2]), sum(F[1, ])))
    Lmut[[i]] <- F
    Fmat <- matrix(c(sum(pmeth == 2), sum(pmeth == 1), sum(pmeth == -2), sum(pmeth == 0)), 2)
    Lall[[i]] <- Fmat
    ## print(fisher.test(Fmat, alternative = "greater"))
    print("p-value")
    print(1 - phyper(Fmat[1,1]-1, sum(Fmat[, 1]), sum(Fmat[, 2]), sum(Fmat[1, ])))
    meas <- list()
    meas$data <- matrix(0, nrow(pmeth), ncol(pmeth))
    colnames(meas$data) <- apply(pmeth, 2, function(x) {
        y <- paste(sort(rownames(pmeth)[which(x %in% c(2,-2))]), collapse = "_")
        return(y)
    })
    measnet <- matrix(0, nrow(pmeth), nrow(pmeth))
    diag(measnet) <- 1
    meas$Nem[[1]] <- measnet
    meas$Theta[[1]] <- ures$res$subtopo
    print("AUC")
    if (i == 1) {
        pitmp <- pifit(ures, meas, D4)
    } else {
        pitmp <- pifit(ures, meas, D4, propagate=FALSE)
    }
    print(pitmp$auc)
    print("Cor")
    print(pitmp$cor)
    Laucor[[i]] <- c(pitmp$auc, pitmp$cor)
}

source("~/Documents/nempi/R/nempi_main.r")
source("~/Documents/nempi/R/nempi_low.r")

P <- Pmut+Pmeth+Pcnv
if (!mutinc) { # include mutations or not (0)
    P[which(Pmut == 1)] <- 0
}
P <- P[which(rownames(P) %in% rownames(Rho)), which(colnames(P) %in% colnames(Rho))]
P <- P[order(rownames(P)), order(colnames(P))]
Ptmp <- cbind(P, matrix(0, nrow(P), sum(!(colnames(Rho) %in% colnames(P)))))
colnames(Ptmp) <- c(colnames(P), colnames(Rho)[which(!(colnames(Rho) %in% colnames(P)))])
P <- Ptmp[, colnames(Rho)]
# P <- t(mytc(methods[[1]]$res$adj))%*%P
P[which(P > 1)] <- 1

## plot aucs:
setEPS()
postscript(paste0("Fig_auc.eps"), height = 5, width = 10)
cols <- c("red","turquoise","orange","darkgreen","brown","blue")
par(mfrow=c(1,2))
for (prop in 1:2) {
    for (i in 1:length(methods)) {
        print(names(methods)[i])
        ures <- methods[[i]]
        meas <- list()
        meas$data <- matrix(0, nrow(P), ncol(P))
        colnames(meas$data) <- apply(P, 2, function(x) {
            y <- paste(sort(rownames(P)[which(x == 1)]), collapse = "_")
            return(y)
        })
        meas$Theta[[1]] <- methods[[1]]$res$subtopo
        print("AUC")
        if (prop == 1) {
            meas$Nem[[1]] <- methods[[1]]$res$adj
        } else {
            meas$Nem[[1]] <- diag(nrow(pmeth))
        }
        if (i == 1) {
            pitmp <- pifit(ures, meas, D4)
        } else {
            pitmp <- pifit(ures, meas, D4, propagate = FALSE)
        }
        if (i == 1) {
            if (prop == 1) {
                main <- "precision-recall curves (propagated)"
            } else {
                main <- "precision-recall curves (not propagated)"
            }
            plot(pitmp$rec, pitmp$ppv, type="l", xlim = c(0,1), ylim = c(0,1), main = main, xlab = "recall", ylab = "precision", col = cols[i])
            abline(h=0.5,col=rgb(0,0,0,1),lty=3)
        } else {
            lines(pitmp$rec, pitmp$ppv, type="l", col = cols[i])
        }
        print(pitmp$auc)
    }
    if (prop == 1) {
        cols <- c("red", "blue", "darkgreen", "brown", "orange", "turquoise")
        legend <- c("svm","neural net","random forest","missForest","knn")
        legend("bottomright",legend=c(expression(NEM~pi),legend),fill=cols,col=cols)
        cols <- c("red","turquoise","orange","darkgreen","brown","blue")
    }
}
dev.off()

source("~/Documents/nempi/R/nempi_main.r")
source("~/Documents/nempi/R/nempi_low.r")

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
    cat(paste0(names(methods)[i], " & ", Lall[[i]][1,1], " & ", Lall[[i]][2,1], " & ", Lall[[i]][2,2], " & ", Lall[[i]][1,2], " & ", round(Lall[[i]][1,1]/(Lall[[i]][1,1]+Lall[[i]][2,1])*100), "\\% & ", round(Lall[[i]][1,1]/(Lall[[i]][1,1]+Lall[[i]][1,2])*100), "\\% & ", ptmp, " & ", round(Laucor[[i]][1], 2), "\\\\\n"))
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

















