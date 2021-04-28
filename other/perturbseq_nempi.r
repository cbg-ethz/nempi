set <- as.numeric(commandArgs(TRUE)[1])
domice <- as.numeric(commandArgs(TRUE)[2])
Pgenes <- as.numeric(commandArgs(TRUE)[3])
unkpct <- as.numeric(commandArgs(TRUE)[4])

library(naturalsort)
library(cluster)
library(Rcpp)
library(e1071)
library(nnet)
library(randomForest)
library(missForest)
library(class)
library(CALIBERrfimpute)

## source("~/Documents/mnem/R/mnems.r"); source("~/Documents/mnem/R/mnems_low.r"); sourceCpp("~/Documents/mnem/src/mm.cpp"); source("~/Documents/nempi/R/nempi_main.r"); source("~/Documents/nempi/R/nempi_low.r"); path <- "~/Mount/Eulershare/perturbseq2/"; set <- 3; domice <- 0; Pgenes <- 10

## uncomment for leo/euler:
source("mnems.r")
source("mnems_low.r")
## sourceCpp("mm.cpp")
sourceCpp(code=readChar("mm.cpp", file.info("mm.cpp")$size))
source("nempi_main.r")
source("nempi_low.r")

path <- "/cluster/work/bewi/members/mpirkl/perturbseq2/"

acc <- function(a,b) {
    gamsave <- a
    gamtmp2 <- b
    auc <- roc <- 0
    ppv <- rec <- spec <- NULL
    for (cut in c(2,seq(1,0, length.out = 100),-1)) {
        gamtmp <- apply(gamsave, 2, function(x) {
            y <- x*0
            y[which(x > cut)] <- 1
            return(y)
        })
        tp <- sum(gamtmp == 1 & gamtmp2 == 1)
        fp <- sum(gamtmp == 1 & gamtmp2 == 0)
        tn <- sum(gamtmp == 0 & gamtmp2 == 0)
        fn <- sum(gamtmp == 0 & gamtmp2 == 1)
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

lods <- readRDS(paste0(path, "L_linnorm_", set, "_sc.rds")) # lods <- lods.bckp
colnames(lods) <- gsub("_[0-9].*","",colnames(lods))
if (set %in% c(4,5)) {
    lods <- lods[,which(!(colnames(lods) == "Unknwn"))]
    search <- "greedy"
    if (set == 4) {
        Sgenes <- sort(sample(getSgenes(lods),Pgenes))
        Gamma0 <- getGamma(lods)
        Gamma0 <- Gamma0[which(rownames(Gamma0) %in% Sgenes),]
        colnames(lods) <- apply(Gamma0,2,function(x) {
            y <- paste(sort(rownames(Gamma0)[which(x==1)]),collapse="_")
            return(y)
        })
        lods <- lods[,which(!(colnames(lods) == ""))]
    }
} else {
    search <- "exhaustive"
}
fres <- nem(lods,search=search)
phi <- fres$adj
n <- ncol(phi)
Gamma <- t(mytc(phi))%*%getGamma(lods)
Gamma[which(Gamma > 1)] <- 1
prevalence <- sum(Gamma==1)/length(Gamma)
phi_gtn_dens <- (sum(mytc(phi)==1)-n)/((n*(n-1))/2)

if (set == 4) {
    if (domice) {
        file <- paste0("nempi_epistasis_", set, "_", domice, "_", Pgenes, "_", unkpct, ".csv")
    } else {
        file <- paste0("nempi_epistasis_", set, "_", Pgenes, "_", unkpct, ".csv")
    }
} else {
    if (domice) {
        file <- paste0("nempi_epistasis_", set, "_", domice, "_", unkpct, ".csv")
    } else {
        file <- paste0("nempi_epistasis_", set, "_", unkpct, ".csv")
    }
}
if (!file.exists(file)) {
    if (domice) {
        write.table(matrix(c("nempi","svm","randomForest","neuralNet","knn","missForest","mice","random","phi_acc","phi_gtn_dens","nempi","svm","randomForest","neuralNet","knn","missForest","mice","random2","prevelance"),1),file=file,sep=",",col.names=FALSE,row.names=FALSE)
    } else {
        write.table(matrix(c("nempi","svm","randomForest","neuralNet","knn","missForest","random","phi_acc","phi_gtn_dens","nempi","svm","randomForest","neuralNet","knn","missForest","random2","prevelance"),1),file=file,sep=",",col.names=FALSE,row.names=FALSE)
    }
}

lods.sub <- lods
unknowns <- sample(1:ncol(lods),floor(unkpct*ncol(lods)))
colnames(lods.sub)[unknowns] <- ""
## nempi
nres <- nempi(lods.sub,search=search)
n <- nrow(phi)
phiacc <- 1-sum(abs(mytc(phi)-mytc(nres$res$adj)))/(n*(n-1))
Gamma2 <- t(mytc(nres$res$adj))%*%nres$Gamma
nacc <- acc(Gamma2,Gamma)
nacc2 <- acc(nres$Gamma,getGamma(lods))
## svm
svmres <- classpi(lods.sub)
svmacc <- acc(svmres$Gamma,Gamma)
svmacc2 <- acc(svmres$Gamma,getGamma(lods))
## random forest
rfres <- classpi(lods.sub,method="randomForest")
rfacc <- acc(rfres$Gamma,Gamma)
rfacc2 <- acc(rfres$Gamma,getGamma(lods))
## neural net
nnres <- classpi(lods.sub,method="nnet",size=2)
nnacc <- acc(nnres$Gamma,Gamma)
nnacc2 <- acc(nnres$Gamma,getGamma(lods))
## knn
train <- t(lods.sub[, which(colnames(lods.sub) != "")])
test <- t(lods.sub[, which(colnames(lods.sub) == "")])
cl <- colnames(lods.sub)[which(colnames(lods.sub) != "")]
tmp <- knn(train, test, cl, prob=TRUE)
lods.sub2 <- lods.sub
colnames(lods.sub2)[which(colnames(lods.sub2) %in% "")] <- as.character(tmp)
knnres <- list()
knnres$Gamma <- getGamma(lods.sub2)
knnres$Gamma <- apply(knnres$Gamma, 2, function(x) return(x/sum(x)))
knnacc <- acc(knnres$Gamma,Gamma)
knnacc2 <- acc(knnres$Gamma,getGamma(lods))
## missForest
mfdata <- cbind(as.data.frame(t(lods.sub)), factor(colnames(lods.sub)))
mfdata[which(mfdata == "", arr.ind = TRUE)] <- NA
sink("NUL")
mfimp <- missForest(mfdata)
sink()
lods.sub2 <- lods.sub
colnames(lods.sub2) <- mfimp$ximp[, ncol(mfimp$ximp)]
mfres <- list()
mfres$Gamma <- getGamma(lods.sub2)
mfres$Gamma <- apply(mfres$Gamma, 2, function(x) return(x/sum(x)))
mfacc <- acc(mfres$Gamma,Gamma)
mfacc2 <- acc(mfres$Gamma,getGamma(lods))
## mice:
if (domice) {
    micedata <- mfdata
    colnames(micedata) <- paste0(LETTERS[1:ncol(micedata)], 1:ncol(micedata))
    sink("NUL")
    miceres <- mice(micedata, method = c(rep('', ncol(micedata)-1), 'rfcat'), m = 1, maxit = 1)
    sink()
    lods.sub2 <- lods.sub
    colnames(lods.sub2)[which(colnames(lods.sub2) %in% "")] <- as.character(miceres$imp[[length(miceres$imp)]][, 1])
    miceres <- list()
    miceres$Gamma <- getGamma(lods.sub2)
    miceres$Gamma <- apply(miceres$Gamma, 2, function(x) return(x/sum(x)))
    miceacc <- acc(miceres$Gamma,Gamma)
    miceacc2 <- acc(miceres$Gamma,getGamma(lods))
}
## random:
rand <- getGamma(lods.sub)
rand[,unknowns] <- runif(nrow(rand)*length(unknowns),0,1) # sample(c(0,1),nrow(rand)*length(unknowns),replace=TRUE)
random <- acc(rand,Gamma)
random2 <- acc(rand,getGamma(lods))

if (domice) {
    write.table(matrix(c(nacc$auc,svmacc$auc,rfacc$auc,nnacc$auc,knnacc$auc,mfacc$auc,miceres$auc,random$auc,phiacc,phi_gtn_dens,nacc2$auc,svmacc2$auc,rfacc2$auc,nnacc2$auc,knnacc2$auc,mfacc2$auc,miceres2$auc,random2$auc,prevalence),1),file=file,append=TRUE,sep=",",col.names=FALSE,row.names=FALSE)
} else {
    write.table(matrix(c(nacc$auc,svmacc$auc,rfacc$auc,nnacc$auc,knnacc$auc,mfacc$auc,random$auc,phiacc,phi_gtn_dens,nacc2$auc,svmacc2$auc,rfacc2$auc,nnacc2$auc,knnacc2$auc,mfacc2$auc,random2$auc,prevalence),1),file=file,append=TRUE,sep=",",col.names=FALSE,row.names=FALSE)
}

stop("success")

## cluster:

system("scp nempi/R/nempi_main.r euler.ethz.ch:")
system("scp nempi/R/nempi_low.r euler.ethz.ch:")
system("scp mnem/R/mnems_low.r euler.ethz.ch:")
system("scp mnem/R/mnems.r euler.ethz.ch:")
system("scp testing/nempi/other/perturbseq.r euler.ethz.ch:")

rm error.txt
rm output.txt
rm .RData

queue=4
ram=8000

set=4
domice=0
Pgenes=10
unkpct=0.1

if [ ${set} == '4' ] && [ ${Pgenes} == '15' ]
then
queue=24
fi

if [ ${set} == '4' ] && [ ${Pgenes} == '10' ] && [ ${unkpct} == '0.1' ]
then
queue=24
fi

bsub -M ${ram} -q normal.${queue}h -n 1 -e error.txt -o output.txt -R "rusage[mem=${ram}]" "R/bin/R --vanilla --silent --no-save --args '${set}' '${domice}' '${Pgenes}' '${unkpct}' < perturbseq.r"

for i in {2..71}; do
bsub -M ${ram} -q normal.${queue}h -n 1 -e error.txt -o output.txt -R "rusage[mem=${ram}]" "R/bin/R --vanilla --silent --no-save --args '${set}' '${domice}' '${Pgenes}' '${unkpct}' < perturbseq.r"
done

for i in {145596019..145596204}; do
bkill ${i}
done

## figures:

supp <- 1
if (supp) {
    layout.mat <- matrix(c(rep(rep(1:3,each=1),2),rep(rep(4:6,each=1),2),
                           rep(rep(7:9,each=1),2),rep(rep(10:12,each=1),2),13:15),9,byrow=1)
    height <- 12
    width <- 8
} else {
    layout.mat <- matrix(c(rep(rep(1:3,each=1),2),rep(rep(4:6,each=1),2),7:9),5,byrow=1)
    height <- 7
    width <- 8
}
setEPS()
postscript("temp.eps", height = height, width = width)
par(mar=c(3.85,4,4,1),oma=c(0,0,0,0))
layout(layout.mat)
path <- "~/Mount/Euler/"
cols <- c("red", "blue", "darkgreen", "brown", "orange", "turquoise", "grey");
data.names <- c("Epistasis 1","Epistasis 2","Epistasis 3","Main (10 P-genes)","Pilot","Main (15 P-genes)")
phiacc <- array(NA,c(100,6,3),list(rows=1:100,data=data.names,unkpct=c(0.1,0.5,0.9)))
ylab <- "area under the precision-recall curve"
for (i in c(5,4,6,1,2,3)) {
    count <- 0
    for (unkpct in c(0.1,0.5,0.9)) {
        count <- count + 1
        if (i == 4) {
            ylim.min <- 0.2
            res <- read.csv(paste0(path, "nempi_epistasis_4_10_", unkpct, ".csv"))
        } else if (i == 6) {
            res <- read.csv(paste0(path, "nempi_epistasis_4_15_", unkpct, ".csv"))
        } else {
            if (i == 5) {
                ylim.min <- 0.4
            } else {
                ylim.min <- 0.75
            }
            res <- read.csv(paste0(path, "nempi_epistasis_", i, "_", unkpct, ".csv"))
        }
        if (length(res$phi_acc)>0) {
            phiacc[1:100,i,count] <- res$phi_acc[1:100]
        }
        print(dim(res))
        res2 <- res[,c(1:7)]
        if ((supp & unkpct != 0.5) | (!supp & unkpct == 0.5)) {
            mnem:::myboxplot(res2,col = cols,ylim=c(ylim.min,1),main=paste0(data.names[i], "\nunobserved: ", unkpct),ylab=ylab,box=1,scatter=1,dens=0,xaxt = "n",border = cols,medcol="black",cex.main=1.5,cex.lab=1.25)
        }
    }
}
if (supp) { cex.leg <- 2 } else { cex.leg <- 1.5 }
par(mar=rep(0, 4))
plot.new()
legend("topleft",legend=c(expression(NEM~pi), "svm", "neural net"),col=c("red", "blue", "darkgreen"),fill=c("red", "blue", "darkgreen"),box.col = "transparent",cex = cex.leg,ncol = 1,border="transparent")
plot.new()
legend("topleft",legend=c("random forest", "missForest","knn"),col=c("brown", "orange", "pink"),fill=c("brown", "orange", "turquoise"),box.col = "transparent",cex = cex.leg,ncol = 1,border="transparent")
plot.new()
legend("topleft",legend=c("random"),col=c("grey"),fill=c("grey"),box.col = "transparent",cex = cex.leg,ncol = 1,border="transparent")
dev.off()

layout.mat <- matrix(c(rep(rep(1:3,each=1),2),rep(rep(4:6,each=1),2),
                       rep(rep(7:9,each=1),2),rep(rep(10:12,each=1),2),13:15),9,byrow=1)
height <- 8
width <- 8
setEPS()
postscript("temp.eps", height = height, width = width)
par(mar=c(3.85,4,4,1),oma=c(0,0,0,0))
layout(layout.mat)
path <- "~/Mount/Euler/"
cols <- c("red", "blue", "darkgreen", "brown", "orange", "turquoise", "grey");
data.names <- c("Epistasis 1","Epistasis 2","Epistasis 3","Main (10 P-genes)","Pilot","Main (15 P-genes)")
ylab <- "area under the precision-recall curve"
for (i in c(4,6)) {
    count <- 0
    for (unkpct in c(0.1,0.5,0.9)) {
        count <- count + 1
        if (i == 4) {
            res <- read.csv(paste0(path, "nempi_epistasis_4_10_", unkpct, ".csv"))
        } else if (i == 6) {
            res <- read.csv(paste0(path, "nempi_epistasis_4_15_", unkpct, ".csv"))
        }
        print(dim(res))
        res2 <- res[,c(1:7)]
        res3 <- res2[which(res[,9]>=median(res[,9])),]
        mnem:::myboxplot(res3,col = cols,ylim=c(0.2,1),main=paste0(data.names[i], "\nunobserved: ", unkpct),ylab=ylab,box=1,scatter=1,dens=0,xaxt = "n",border = cols,medcol="black",xlab=expression(dense ~ phi),cex.main=1.5,cex.lab=1.25)
        res3 <- res2[which(res[,9]<median(res[,9])),]
        mnem:::myboxplot(res3,col = cols,ylim=c(0.2,1),main=paste0(data.names[i], "\nunobserved: ", unkpct),ylab=ylab,box=1,scatter=1,dens=0,xaxt = "n",border = cols,medcol="black",xlab=expression(sparse ~ phi),cex.main=1.5,cex.lab=1.25)
    }
}
cex.leg <- 2
par(mar=rep(0, 4))
plot.new()
legend("topleft",legend=c(expression(NEM~pi), "svm", "neural net"),col=c("red", "blue", "darkgreen"),fill=c("red", "blue", "darkgreen"),box.col = "transparent",cex = cex.leg,ncol = 1,border="transparent")
plot.new()
legend("topleft",legend=c("random forest", "missForest","knn"),col=c("brown", "orange", "pink"),fill=c("brown", "orange", "turquoise"),box.col = "transparent",cex = cex.leg,ncol = 1,border="transparent")
plot.new()
legend("topleft",legend=c("random"),col=c("grey"),fill=c("grey"),box.col = "transparent",cex = cex.leg,ncol = 1,border="transparent")
dev.off()

pdf("FigS_crispr_phi.pdf",width=9,height=3)
unkpcts <- c(0.1,0.5,0.9)
par(mfrow=c(1,3))
for (i in 1:3) {
    #mnem:::myboxplot(phiacc[,c(1:3,5,4,6),i],col = cols[1],ylim=c(0,1),main=expression(accuracy ~ of ~ phi),ylab="normalised hamming distance",box=1,scatter=1,dens=0,xaxt = "n",border = cols[1],medcol="black")
    mnem:::myboxplot(phiacc[,c(1:3,5,4,6),i],col = cols[1],ylim=c(0,1),main=bquote(accuracy ~ of ~ phi ~ (unknowns: ~ .(unkpcts[i]))),ylab="normalised hamming distance",box=1,scatter=1,dens=0,xaxt = "n",border = cols[1],medcol="black",cex.main=1.5,cex.lab=1.25)
    axis(1,1:6,c(expression(Epistasis ~ 1),"Epistasis 2\n","Epistasis 3","Pilot\n","Main (10)","Main (15)\n"),padj=1,las=1)
}
dev.off()

## old:

pdf("temp.pdf",width=10,height=18)
par(mfrow=c(6,5))
phiacc <- list()
count <- 0
for (unkpct in c(0.1,0.5,0.9)) {
    count <- count + 1
    path <- "~/Mount/Euler/"
    cols <- c("red", "blue", "darkgreen", "brown", "orange", "turquoise", "grey");
    data.names <- c("Epistasis 1","Epistasis 2","Epistasis 3","Main (10 P-genes)","Pilot","Main (15 P-genes)")
    ylab <- "area under the precision-recall curve"
    for (i in c(1,2,4,3,5,6)) {
        if (i == 4) {
            res <- read.csv(paste0(path, "nempi_epistasis_4_10_", unkpct, ".csv"))
        } else if (i == 6) {
            res <- read.csv(paste0(path, "nempi_epistasis_4_15_", unkpct, ".csv"))
        } else {
            res <- read.csv(paste0(path, "nempi_epistasis_", i, "_", unkpct, ".csv"))
        }
        if (i == 1) {
            phiacc[[count]] <- res$phi_acc
        } else {
            phiacc[[count]] <- cbind(phiacc[[count]],res$phi_acc)
        }
        print(dim(res))
        res2 <- res[,c(1:7)]
        mnem:::myboxplot(res2,col = cols,ylim=c(0,1),main=data.names[i],ylab=ylab,box=1,scatter=1,dens=0,xaxt = "n",border = cols,medcol="black")
        if (i %in% c(4,6)) {
            res3 <- res2[which(res[,9]>=median(res[,9])),]
            mnem:::myboxplot(res3,col = cols,ylim=c(0,1),main=data.names[i],ylab=ylab,box=1,scatter=1,dens=0,xaxt = "n",border = cols,medcol="black",xlab=expression(dense ~ phi))
            res3 <- res2[which(res[,9]<median(res[,9])),]
            mnem:::myboxplot(res3,col = cols,ylim=c(0,1),main=data.names[i],ylab=ylab,box=1,scatter=1,dens=0,xaxt = "n",border = cols,medcol="black",xlab=expression(sparse ~ phi))
        }
    }
}
dev.off()

## mnem:

library(Rgraphviz)

path <- "~/Mount/Eulershare/perturbseq2/"

lods <- readRDS(paste0(path,"L_linnorm_1_sc.rds"))
colnames(lods) <- gsub("_[0-9].*","",colnames(lods))

mres <- mnem(lods,search="exhaustive",k=3,complete=TRUE)

par(mfrow=c(2,3))
mnem::plotConvergence(mres)

plot(mres)

mkres <- mnemk(lods,search="exhaustive",complete=TRUE,start=10)

par(mfrow=c(1,2))
plot(mkres$lls,type="b")
plot(mkres$ics,type="b")

par(mfrow=c(2,2))
mnem::plotConvergence(mkres$best)

plot(mkres$best,legend=TRUE,showdata=TRUE)

pca <- prcomp(lods)

cols <- numeric(ncol(lods))
for (i in 1:length(unique(colnames(lods)))) {
    cols[which(colnames(lods) %in% unique(colnames(lods))[i])] <- i
}

plot(pca$rotation[,1:2], col = cols)

