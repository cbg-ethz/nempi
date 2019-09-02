## Sgenes <- 100; wrong <- 0; highnoise <- 1; doclass <- 1; knowns <- 8; runs2 <- NA

library(naturalsort)
library(nem)
library(cluster)
library(Rcpp)
library(e1071)
library(nnet)
library(randomForest)
library(missForest)
library(class)
library(CALIBERrfimpute)

## source("~/Documents/mnem/R/mnems.r"); source("~/Documents/mnem/R/mnems_low.r"); sourceCpp("~/Documents/mnem/src/mm.cpp"); source("~/Documents/testing/nempi/R/nempi_main.r"); source("~/Documents/testing/nempi/R/nempi_low.r")

## uncomment for leo/euler:
source("mnems.r")
source("mnems_low.r")
## sourceCpp("mm.cpp")
sourceCpp(code=readChar("mm.cpp", file.info("mm.cpp")$size))
source("nempi_main.r")
source("nempi_low.r")
Sgenes <- as.numeric(commandArgs(TRUE)[1])
wrong <- as.numeric(commandArgs(TRUE)[2])
highnoise <- as.numeric(commandArgs(TRUE)[3])
doclass <- as.numeric(commandArgs(TRUE)[4])
knowns <- as.numeric(commandArgs(TRUE)[5])
runs2 <- as.numeric(commandArgs(TRUE)[6])
##

if (is.na(knowns)) { knowns <- Sgenes }

Egenes <- 10
if (knowns < Sgenes) {
    nCells <- 1000
} else {
    nCells <- Sgenes*10*2
}

runs <- 100
noises <- c(1,3,5)
lost <- c(0.1, 0.5, 0.9)
multi <- c(0.2, 0.1)

complete <- 1
combi <- 1
combisvm <- NULL
null <- 1
badCells <- 0.1

if (wrong) { lost <- lost[1:2] }

if (knowns < Sgenes) { lost <- lost[c(1,3)] }

result <- array(0, c(runs, length(noises), length(lost), 6, 13), list(rep("runs", runs), noises, lost, c("net", "known", "unknown", "cor", "time", "nullsamples"), c("unem_s_n", "unem_s_f_n", "svm", "svm_f", "nnet", "nnet_f", "rf", "rf_f", "empty", "random", "missForest", "mice", "knn")))

getfalse <- 0

for (i in 1:runs) {
    if (!is.na(runs2)) {
        if (i != runs2) {
            next()
        }
    }
    cat(i)
    for (j in 1:length(noises)) {

        ## if (knowns < Sgenes & j == 2) { next() }

        for (k in 1:length(lost)) {

            ## i <- j <- 1; k <- 2

            if (knowns < Sgenes) {
                Sgenes2 <- sort(sample(Sgenes, knowns))
                Sgenes3 <- seq_len(Sgenes)[-Sgenes2]
            } else {
                Sgenes2 <- Sgenes
            }

            noise1 <- noises[j]
            if (!wrong) {
                noise2 <- lost[k]
            } else {
                noise2 <- 0.5
            }

            sim <- simData(Sgenes = Sgenes, Egenes = Egenes, mw = 1, nCells = nCells, Nem = 1, multi = multi, badCells = floor(nCells*badCells))
            sdata <- sim$data
            sdata[which(sim$data == 1)] <- rnorm(sum(sim$data == 1), 1, noise1)
            sdata[which(sim$data == 0)] <- rnorm(sum(sim$data == 0), -1, noise1)

            lost2 <- sort(sample(1:ncol(sdata), floor(noise2*ncol(sdata))))

            stop <- FALSE
            while(!stop) {
                lost2 <- sort(sample(1:ncol(sdata), floor(noise2*ncol(sdata))))
                if (highnoise) {
                    if (all(Sgenes2 %in% unique(unlist(strsplit(colnames(sdata)[-lost2], "_"))))) { stop <- TRUE }
                } else {
                    if (all(Sgenes2 %in% unique(colnames(sdata)[-lost2]))) { stop <- TRUE }
                }
            }
            colnames(sdata)[lost2] <- ""
            if (wrong) {
                multi2 <- NULL
                if (multi[1] != FALSE) { multi2 <- multi }
                sdata2 <- sdata
                stop <- FALSE
                while(!stop) {
                    sdata <- sdata2
                    colnames(sdata)[sample(which(!(colnames(sdata) %in% "")), floor(lost[k]*sum(!(colnames(sdata) %in% ""))))] <- paste(naturalsort(sample(1:Sgenes, sample(1:(length(multi2)+1), 1))), collapse = "_")
                    lost2 <- sort(sample(1:ncol(sdata), floor(noise2*ncol(sdata))))
                    if (highnoise) {
                        if (all(Sgenes2 %in% unique(unlist(strsplit(colnames(sdata)[-lost2], "_"))))) { stop <- TRUE }
                    } else {
                        if (all(Sgenes2 %in% unique(colnames(sdata)[-lost2]))) { stop <- TRUE }
                    }
                }
            }
            D <- sdata

            if (multi[1] != FALSE) {
                multi2 <- TRUE
            }

            paras <- expand.grid(c(0,1), c(0,1), c(0,1))

            if (knowns < Sgenes) {
                colnames(D) <- gsub(paste(paste0("_", Sgenes3, "_"), collapse = "|"), "_", colnames(D))
                colnames(D) <- gsub(paste(c(paste0("^", Sgenes3, "_"),
                                            paste0("_", Sgenes3, "$"),
                                            paste0("^", Sgenes3, "$")), collapse = "|"), "", colnames(D))
                colnames(D) <- gsub(paste(c(paste0("^", Sgenes3, "_"),
                                            paste0("_", Sgenes3, "$"),
                                            paste0("^", Sgenes3, "$")), collapse = "|"), "", colnames(D))
                full <- mytc(sim$Nem[[1]])[Sgenes3, Sgenes2]
                Sgenes4 <- rownames(full)[which(apply(full, 1, sum) == 0)]
            }

            ## source("~/Documents/testing/nempi/R/nempi_main.r"); source("~/Documents/testing/nempi/R/nempi_low.r");

            for (s in 2) {
                type <- "null"
                start <- as.numeric(format(Sys.time(), "%s"))
                ures <- nempi(D, full = paras[s, 1], type = type, soft = paras[s+2, 2], combi = combi, multi = multi2, complete = complete, null = null)
                result[i, j, k, 5, s] <- as.numeric(format(Sys.time(), "%s")) - start
                if (knowns < Sgenes) {
                    tmp <- pifit(ures, sim, D, knowns = Sgenes2)
                } else {
                    tmp <- pifit(ures, sim, D)
                }
                par(mfrow=c(2,3))
                result[i, j, k, 1, s] <- tmp$net
                result[i, j, k, 2, s] <- tmp$known[1]
                result[i, j, k, 3, s] <- tmp$known[2]
                result[i, j, k, 4, s] <- tmp$cor
                result[i, j, k, 6, s] <- sum(apply(ures$Gamma, 2, sum) < 0.5 & colnames(sim$data) %in% Sgenes4)/sum(colnames(sim$data) %in% Sgenes4)
                if (!all(ures$lls[2:length(ures$lls)] - ures$lls[1:(length(ures$lls)-1)] >= 0)) { getfalse <- getfalse + 1 }
            }
            print("nempi")

            ## result[i, j, k, , s]

            ## if (getfalse > 0) { stop() }

            ## source("~/Documents/testing/nempi/R/nempi_main.r"); source("~/Documents/testing/nempi/R/nempi_low.r");

            if (doclass) {
                for (s in c(2,4,6)) {
                    start <- as.numeric(format(Sys.time(), "%s"))
                    if (s %in% 1:2) {
                        ures <- classpi(D, full = paras[s, 1])
                    }
                    if (s %in% 3:4) {
                        ures <- classpi(D, full = paras[s, 1], method = "nnet", size = 2)
                    }
                    if (s %in% 5:6) {
                        ures <- classpi(D, full = paras[s, 1], method = "randomForest")
                    }
                    result[i, j, k, 5, s+2] <- as.numeric(format(Sys.time(), "%s")) - start
                    if (knowns < Sgenes) {
                        tmp <- pifit(ures, sim, D, propagate = 0, knowns = Sgenes2)
                    } else {
                        tmp <- pifit(ures, sim, D, propagate = 0)
                    }
                    result[i, j, k, 1, s+2] <- tmp$net
                    result[i, j, k, 2, s+2] <- tmp$known[1]
                    result[i, j, k, 3, s+2] <- tmp$known[2]
                    result[i, j, k, 4, s+2] <- tmp$cor
                    result[i, j, k, 6, s+2] <- sum(apply(ures$Gamma, 2, sum) < 0.5 & colnames(sim$data) %in% Sgenes4)/sum(colnames(sim$data) %in% Sgenes4)

                    sum(apply(ures$Gamma, 2, max) < 1 - apply(ures$Gamma, 2, sum) & colnames(sim$data) %in% Sgenes4)/sum(colnames(sim$data) %in% Sgenes4)


                }
                cat(s)
            }
            print("class")

            ## empty:

            D2 <- D[, which(colnames(D) != "")]
            Rho <- getRho(D2)
            tmp <- mynem(D2, Rho = Rho, multi = TRUE)
            n <- Sgenes
            A <- mytc(tmp$adj)
            B <- mytc(sim$Nem[[1]])
            if (knowns < Sgenes) {
                B <- B[which(rownames(B) %in% Sgenes2), which(colnames(B) %in% Sgenes2)]
            }
            result[i, j, k, 1, 9] <- (n*(n-1) - sum(abs(A - B)))/(n*(n-1))

            ## random:

            start <- as.numeric(format(Sys.time(), "%s"))
            D2 <- D
            if (knowns < Sgenes) {
                colnames(D2)[which(colnames(D) %in% "")] <- sample(Sgenes2, sum(colnames(D) %in% ""), replace = TRUE)
            } else {
                colnames(D2)[which(colnames(D) %in% "")] <- sample(1:Sgenes, sum(colnames(D) %in% ""), replace = TRUE)
            }
            tmp <- mynem(D2, multi = TRUE)
            Gamma <- getGamma(D2)
            ures <- list()
            ures$Gamma <- apply(Gamma, 2, function(x) return(x/sum(x)))
            ures$res <- list()
            ures$res$adj <- tmp$adj
            ures$null <- TRUE
            ures$combi <- 1
            result[i, j, k, 5, 10] <- as.numeric(format(Sys.time(), "%s")) - start
            if (knowns < Sgenes) {
                tmp <- pifit(ures, sim, D, propagate = 0, knowns = Sgenes2)
            } else {
                tmp <- pifit(ures, sim, D, propagate = 0)
            }
            result[i, j, k, 1, 10] <- tmp$net
            result[i, j, k, 2, 10] <- tmp$known[1]
            result[i, j, k, 3, 10] <- tmp$known[2]
            result[i, j, k, 4, 10] <- tmp$cor
            result[i, j, k, 6, 10] <- sum(apply(ures$Gamma, 2, sum) < 0.5 & colnames(sim$data) %in% Sgenes4)/sum(colnames(sim$data) %in% Sgenes4)

            ## missForest

            start <- as.numeric(format(Sys.time(), "%s"))
            mfdata <- cbind(as.data.frame(t(D)), colnames(D))
            mfdata[which(mfdata == "", arr.ind = TRUE)] <- NA
            sink("NUL")
            mfimp <- missForest(mfdata)
            sink()
            D2 <- D
            colnames(D2) <- mfimp$ximp[, ncol(mfimp$ximp)]
            tmp <- mynem(D2, multi = TRUE)
            Gamma <- getGamma(D2)
            ures <- list()
            ures$Gamma <- apply(Gamma, 2, function(x) return(x/sum(x)))
            ures$res <- list()
            ures$res$adj <- tmp$adj
            ures$null <- TRUE
            ures$combi <- 1
            result[i, j, k, 5, 11] <- as.numeric(format(Sys.time(), "%s")) - start
            if (knowns < Sgenes) {
                tmp <- pifit(ures, sim, D, propagate = 0, knowns = Sgenes2)
            } else {
                tmp <- pifit(ures, sim, D, propagate = 0)
            }
            result[i, j, k, 1, 11] <- tmp$net
            result[i, j, k, 2, 11] <- tmp$known[1]
            result[i, j, k, 3, 11] <- tmp$known[2]
            result[i, j, k, 4, 11] <- tmp$cor
            result[i, j, k, 6, 11] <- sum(apply(ures$Gamma, 2, sum) < 0.5 & colnames(sim$data) %in% Sgenes4)/sum(colnames(sim$data) %in% Sgenes4)
            print("missForest")

            ## mice:

            if (Sgenes <= 100) {
                start <- as.numeric(format(Sys.time(), "%s"))
                micedata <- mfdata
                colnames(micedata) <- paste0(LETTERS[1:ncol(micedata)], 1:ncol(micedata))
                sink("NUL")
                miceres <- mice(micedata, method = c(rep('', ncol(micedata)-1), 'rfcat'), m = 1, maxit = 1)
                sink()
                D2 <- D
                colnames(D2)[which(colnames(D2) %in% "")] <- as.character(miceres$imp[[length(miceres$imp)]][, 1])
                tmp <- mynem(D2, multi = TRUE)
                Gamma <- getGamma(D2)
                ures <- list()
                ures$Gamma <- apply(Gamma, 2, function(x) return(x/sum(x)))
                ures$res <- list()
                ures$res$adj <- tmp$adj
                ures$null <- TRUE
                ures$combi <- 1
                result[i, j, k, 5, 12] <- as.numeric(format(Sys.time(), "%s")) - start
                if (knowns < Sgenes) {
                    tmp <- pifit(ures, sim, D, propagate = 0, knowns = Sgenes2)
                } else {
                    tmp <- pifit(ures, sim, D, propagate = 0)
                }
                result[i, j, k, 1, 12] <- tmp$net
                result[i, j, k, 2, 12] <- tmp$known[1]
                result[i, j, k, 3, 12] <- tmp$known[2]
                result[i, j, k, 4, 12] <- tmp$cor
                result[i, j, k, 6, 12] <- sum(apply(ures$Gamma, 2, sum) < 0.5 & colnames(sim$data) %in% Sgenes4)/sum(colnames(sim$data) %in% Sgenes4)
                print("mice")
            }

            ## knn:

            start <- as.numeric(format(Sys.time(), "%s"))
            train <- t(D[, which(colnames(D) != "")])
            test <- t(D[, which(colnames(D) == "")])
            cl <- colnames(D)[which(colnames(D) != "")]
            knnres <- knn(train, test, cl, k = 3, prob=TRUE)
            D2 <- D
            colnames(D2)[which(colnames(D2) %in% "")] <- as.character(knnres)
            tmp <- mynem(D2, multi = TRUE)
            Gamma <- getGamma(D2)
            ures <- list()
            ures$Gamma <- apply(Gamma, 2, function(x) return(x/sum(x)))
            ures$res <- list()
            ures$res$adj <- tmp$adj
            ures$null <- TRUE
            ures$combi <- 1
            result[i, j, k, 5, 13] <- as.numeric(format(Sys.time(), "%s")) - start
            if (knowns < Sgenes) {
                tmp <- pifit(ures, sim, D, propagate = 0, knowns = Sgenes2)
            } else {
                tmp <- pifit(ures, sim, D, propagate = 0)
            }
            result[i, j, k, 1, 13] <- tmp$net
            result[i, j, k, 2, 13] <- tmp$known[1]
            result[i, j, k, 3, 13] <- tmp$known[2]
            result[i, j, k, 4, 13] <- tmp$cor
            result[i, j, k, 6, 13] <- sum(apply(ures$Gamma, 2, sum) < 0.5 & colnames(sim$data) %in% Sgenes4)/sum(colnames(sim$data) %in% Sgenes4)
            print("knn")

            ## result[i, j, k, , c((1:4)*2, 11,12,13,10)]

        }
    }
}

if (is.na(runs2)) {
    if (wrong) {
        save(result, getfalse, lost, file = paste("unem_misslabeled", highnoise, complete, noise2, Sgenes, Egenes, nCells, paste(c(multi, ".rda"), collapse = ""), sep = "_"))
    } else {
        if (knowns < Sgenes) {
            save(result, getfalse, lost, file = paste("unem_missing_unknowns", knowns, highnoise, complete, Sgenes, Egenes, nCells, paste(c(multi, ".rda"), collapse = ""), sep = "_"))
        } else {
            save(result, getfalse, lost, file = paste("unem_missing", highnoise, complete, Sgenes, Egenes, nCells, paste(c(multi, ".rda"), collapse = ""), sep = "_"))
        }
    }
} else {
    if (wrong) {
        save(result, getfalse, lost, file = paste("nempi/unem_misslabeled", runs2, highnoise, complete, noise2, Sgenes, Egenes, nCells, paste(c(multi, ".rda"), collapse = ""), sep = "_"))
    } else {
        if (knowns < Sgenes) {
            save(result, getfalse, lost, file = paste("nempi/unem_missing_unknowns", runs2, knowns, highnoise, complete, Sgenes, Egenes, nCells, paste(c(multi, ".rda"), collapse = ""), sep = "_"))
        } else {
            save(result, getfalse, lost, file = paste("nempi/unem_missing", runs2, highnoise, complete, Sgenes, Egenes, nCells, paste(c(multi, ".rda"), collapse = ""), sep = "_"))
        }
    }
}

stop("done")

## analyse ll decrease:

source("~/Documents/testing/R/nempi/nempi.r")
source("~/Documents/testing/R/nempi/nempi_low.r")
ures <- nempi(D, full = paras[s, 1], type = type, soft = paras[s, 3], combi = combi, multi = multi2, complete = complete, null = 0)
pifit(ures, sim, D)

sum(ures$lls[2:length(ures$lls)] - ures$lls[1:(length(ures$lls)-1)] < 0)

par(mfrow=c(2,3))
plotConvergence.nempi(ures, type = "b", col = "red")

## start pipeline:

Sgenes <- 5; wrong <- 1; highnoise <- 1; source("~/Documents/testing/nempi/vignettes/nempi_app.r");

Sgenes <- 5; wrong <- 0; highnoise <- 1; source("~/Documents/testing/nempi/vignettes/nempi_app.r");

Sgenes <- 10; wrong <- 1; highnoise <- 1; source("~/Documents/testing/nempi/vignettes/nempi_app.r");

Sgenes <- 10; wrong <- 0; highnoise <- 1; source("~/Documents/testing/nempi/vignettes/nempi_app.r");

Sgenes <- 20; wrong <- 1; highnoise <- 1; source("~/Documents/testing/nempi/vignettes/nempi_app.r");

Sgenes <- 20; wrong <- 0; highnoise <- 1; source("~/Documents/testing/nempi/vignettes/nempi_app.r");

## ## on leo/euler:

## leo:

module load r/3.5.1

module load curl/7.53.1

module load gmp/6.1.2

## euler:

module load bioconductor/3.6

module load curl/7.49.1

module load gmp/5.1.3

##

## BiocManager::install("mnem")

ram=10000

rm error.txt
rm output.txt
rm .RData

bsub -M ${ram} -q normal.24h -n 1 -e error.txt -o output.txt -R "rusage[mem=${ram}]" "R/bin/R --silent --no-save --args '15' '0' '1' '1' < nempi_app.r"

bsub -M ${ram} -q normal.24h -n 1 -e error.txt -o output.txt -R "rusage[mem=${ram}]" "R/bin/R --silent --no-save --args '15' '1' '1' '1' < nempi_app.r"

bsub -M ${ram} -q normal.24h -n 1 -e error.txt -o output.txt -R "rusage[mem=${ram}]" "R/bin/R --silent --no-save --args '10' '0' '1' '1' < nempi_app.r"

bsub -M ${ram} -q normal.24h -n 1 -e error.txt -o output.txt -R "rusage[mem=${ram}]" "R/bin/R --silent --no-save --args '10' '1' '1' '1' < nempi_app.r"

bsub -M ${ram} -q normal.24h -n 1 -e error.txt -o output.txt -R "rusage[mem=${ram}]" "R/bin/R --silent --no-save --args '5' '0' '1' '1' < nempi_app.r"

bsub -M ${ram} -q normal.24h -n 1 -e error.txt -o output.txt -R "rusage[mem=${ram}]" "R/bin/R --silent --no-save --args '5' '1' '1' '1' < nempi_app.r"

## new with unknowns

ram=10000

rm error.txt
rm output.txt
rm .RData

bsub -M ${ram} -q normal.24h -n 1 -e error.txt -o output.txt -R "rusage[mem=${ram}]" "R/bin/R --silent --no-save --args '50' '0' '1' '1' '8' < nempi_app.r"

ram=10000

rm error.txt
rm output.txt
rm .RData

bsub -M ${ram} -q normal.120h -n 1 -e error.txt -o output.txt -R "rusage[mem=${ram}]" "R/bin/R --silent --no-save --args '100' '0' '1' '1' '8' < nempi_app.r"

## parallel:

ram=30000
queue=4h
Sgenes=200

rm error.txt
rm output.txt
rm .RData

bsub -M ${ram} -q normal.${queue} -n 1 -e error.txt -o output.txt -R "rusage[mem=${ram}]" "R/bin/R --silent --no-save --args '${Sgenes}' '0' '1' '1' '8' '1' < nempi_app.r"

for i in {2..100}; do
        bsub -M ${ram} -q normal.${queue} -n 1 -e error.txt -o output.txt -R "rusage[mem=${ram}]" "R/bin/R --silent --no-save --args '${Sgenes}' '0' '1' '1' '8' '${i}' < nempi_app.r"
done

## load results:

## combine parallel:

path <- "~/Mount/Euler/nempi/"

Sgenes <- 200
highnoise <- 1
complete <- 1
Egenes <- 10
nCells <- 1000 # Sgenes*10*2
multi <- c(0.2, 0.1)
noise2 <- 0.5
runs <- 100

knowns <- 8

for (runs2 in 1:runs) {
    if (!file.exists(paste(paste0(path, "unem_missing_unknowns"), runs2, knowns, highnoise, complete, Sgenes, Egenes, nCells, paste(c(multi, ".rda"), collapse = ""), sep = "_"))) { cat(runs2); next() }
    load(paste(paste0(path, "unem_missing_unknowns"), runs2, knowns, highnoise, complete, Sgenes, Egenes, nCells, paste(c(multi, ".rda"), collapse = ""), sep = "_"))
    if (runs2 == 1) {
        result2 <- result
    }
    result2[runs2,,,,] <- result[runs2,,,,]
}
result <- result2

result[which(is.na(result) == TRUE)] <- 0

## save(result, lost, file = paste("unem_missing_unknowns", knowns, highnoise, complete, Sgenes, Egenes, nCells, paste(c(multi, ".rda"), collapse = ""), sep = "_"))

## normal pipeline:

Sgenes <- 5
highnoise <- 1
complete <- 1
Egenes <- 10
nCells <- Sgenes*10*2
multi <- c(0.2, 0.1)

## Documents/nempi/old2/

noise2 <- 0.5
load(paste("~/Documents/unem_misslabeled", highnoise, complete, noise2, Sgenes, Egenes, nCells, paste(c(multi, ".rda"), collapse = ""), sep = "_"))

## load(paste("~/Documents/unem_missing", highnoise, complete, Sgenes, Egenes, nCells, paste(c(multi, ".rda"), collapse = ""), sep = "_"))

result[which(is.na(result) == TRUE)] <- 0

wrong <- 0
paras <- expand.grid(c(0,1), c(0,1), c(0,1))
box <- 1
scatter <- ""#"random"
dens <- 0
show <- c((1:4)*2, 11, 12, 13, 10)
show2 <- c(3,7) # 3
if (wrong) {
    width0 <- length(lost)*(10/3)
    height0 <- 3*length(show2)
    main <- "miss-labeled"
} else {
    width0 <- length(lost)*(10/3)
    height0 <- 3*length(show2)
    main <- "missing"
}
ymin <- 0.5
ymin2 <- 0
pdf("temp.pdf", width = width0, height = height0)
par(mfcol=c(length(show2), length(lost)))
for (k in 1:length(lost)) {
    cols <- rgb(c(1,0,0,0.5),c(0,1,0,0.5),c(0,0,1,0.5), 0.75)
    cols <- rep(cols[1:length(show)], 3)
    if (1 %in% show2) {
        myboxplot(cbind(result[,1,k,2,show], result[,2,k,2,show], result[,3,k,2,show]), col = cols, ylim = c(ymin,1), main = paste0(main, ": ", lost[k]), xlab = expression(sigma), ylab = "known prediction", box = box, scatter = scatter, dens = dens, xaxt = "n")
        axis(1, length(show)/2 + 0.5 + c(0, length(show), 2*length(show)), c(1,3,5))
        abline(v=c(length(show)*(1:2)+0.5), col = 1)
        addgrid(x=c(0,1,0.2), y=c(0,length(show)*3+1,1))
    }
    if (6 %in% show2) {
        myboxplot((cbind(result[,1,k,2,show], result[,2,k,2,show], result[,3,k,2,show])+cbind(result[,1,k,3,show], result[,2,k,3,show], result[,3,k,3,show]))/2,
                  col = cols, ylim = c(ymin,1), main = paste0(main, ": ", lost[k]), xlab = expression(sigma), ylab = "balanced accuracy", box = box, scatter = scatter, dens = dens, xaxt = "n")
        axis(1, length(show)/2 + 0.5 + c(0, length(show), 2*length(show)), c(1,3,5))
        abline(v=c(length(show)*(1:2)+0.5), col = 1)
        addgrid(x=c(0,1,0.2), y=c(0,length(show)*3+1,1))
    }
    if (2 %in% show2) {
        myboxplot(cbind(result[,1,k,3,show], result[,2,k,3,show], result[,3,k,3,show]), col = cols, ylim = c(ymin,1), main = paste0(main, ": ", lost[k]), xlab = expression(sigma), ylab = "unknown prediction", box = box, scatter = scatter, dens = dens, xaxt = "n")
        axis(1, length(show)/2 + 0.5 + c(0, length(show), 2*length(show)), c(1,3,5))
        abline(v=c(length(show)*(1:2)+0.5), col = 1)
        addgrid(x=c(0,1,0.2), y=c(0,length(show)*3+1,1))
    }
    if (3 %in% show2) {
        myboxplot(cbind(result[,1,k,4,show], result[,2,k,4,show], result[,3,k,4,show]), col = cols, ylim = c(ymin2,1), main = paste0(main, ": ", lost[k]), xlab = expression(sigma), ylab = "perturbation profile correlation", box = box, scatter = scatter, dens = dens, xaxt = "n")
        axis(1, length(show)/2 + 0.5 + c(0, length(show), 2*length(show)), c(1,3,5))
        abline(v=c(length(show)*(1:2)+0.5), col = 1)
        addgrid(x=c(-1,1,0.5), y=c(0,length(show)*3+1,1))
    }
    if (4 %in% show2) {
        myboxplot(cbind(result[,1,k,1,show], result[,2,k,1,show], result[,3,k,1,show]), col = cols, ylim = c(ymin,1), main = paste0(main, ": ", lost[k]), xlab = expression(sigma), ylab = "network accuracy", box = box, scatter = scatter, dens = dens, xaxt = "n")
        axis(1, length(show)/2 + 0.5 + c(0, length(show), 2*length(show)), c(1,3,5))
        abline(v=c(length(show)*(1:2)+0.5), col = 1)
        addgrid(x=c(0,1,0.2), y=c(0,length(show)*3+1,1))
    }
    if (5 %in% show2) {
        ymax <- max(cbind(result[,1,k,5,show], result[,2,k,5,show], result[,3,k,5,show]))
        myboxplot(cbind(result[,1,k,5,show], result[,2,k,5,show], result[,3,k,5,show]), ylim = c(0, ymax), col = cols, main = paste0(main, ": ", lost[k]), xlab = expression(sigma), ylab = "seconds", box = box, scatter = scatter, dens = dens, xaxt = "n")
        axis(1, length(show)/2 + 0.5 + c(0, length(show), 2*length(show)), c(1,3,5))
        abline(v=c(length(show)*(1:2)+0.5), col = 1)
        addgrid(x=c(0,1,0.2), y=c(0,length(show)*3+1,1))
    }
    if (7 %in% show2) {
        ymax <- max(cbind(result[,1,k,6,show], result[,2,k,6,show], result[,3,k,6,show]))
        myboxplot(cbind(result[,1,k,6,show], result[,2,k,6,show], result[,3,k,6,show]), ylim = c(0, ymax), col = cols, main = paste0(main, ": ", lost[k]), xlab = expression(sigma), ylab = "ratio of null samples to uninformative samples", box = box, scatter = scatter, dens = dens, xaxt = "n")
        axis(1, length(show)/2 + 0.5 + c(0, length(show), 2*length(show)), c(1,3,5))
        abline(v=c(length(show)*(1:2)+0.5), col = 1)
        addgrid(x=c(0,1,0.2), y=c(0,length(show)*3+1,1))
    }
}
dev.off()

## several Sgenes in one plot:

wrong <- 0 # 0,1
knowns <- 8 # 800, 8
show2 <- c(4) # 4,6
shownoise <- c(1,2,3)

highnoise <- 1
complete <- 1
Egenes <- 10
multi <- c(0.2, 0.1)
noise2 <- 0.5
if (wrong) {
    lost <- c(0.1, 0.5)
} else {
    lost <- c(0.1, 0.5, 0.9)
}
if (knowns == 0) {
    doSgenes <- c(5,10,15)#,20)
} else {
    doSgenes <- c(50,100)#,200)
    lost <- c(0.1, 0.9)
}

show <- c((1:4)*2, 11, 12, 13, 10)
cols <- c("red", "blue", "darkgreen", "brown", "orange", "pink", "turquoise", "grey")
cols <- rep(cols[1:length(show)], 3)

if (wrong | knowns > 0) {
    parcols <- 2
    height0 <- 5.5
    lost2 <- lost[1:2]
    leg <- 0
} else {
    parcols <- 3
    height0 <- 8
    lost2 <- lost
    leg <- 2
}
paras <- expand.grid(c(0,1), c(0,1), c(0,1))
box <- 1
scatter <- ""#"random"
dens <- 0
ymin <- 0.5
ymin2 <- 0

if (6 %in% show2 | knowns == 0) {
    pdf("temp.pdf", height = ((8/3)*length(doSgenes)+leg)*length(show2), width = height0)
} else {
    a <- (5/4)
    pdf("temp.pdf", height = (((8/3)*length(doSgenes)+leg)*length(show2))*a, width = height0*a)
}
if (wrong == 0 & knowns == 0) {
    par(mfrow=c((length(doSgenes)+1)*length(show2), parcols))
} else {
    par(mfrow=c(length(doSgenes)*length(show2), parcols))
}
for (Sgenes in doSgenes) {
    for (k in 1:length(lost2)) {
        for (s in show2) {
            nCells <- 1000 # Sgenes*10*2
            if (wrong) {
                noise2 <- 0.5
                load(paste("~/Documents/unem_misslabeled", highnoise, complete, noise2, Sgenes, Egenes, nCells, paste(c(multi, ".rda"), collapse = ""), sep = "_"))
            } else {
                main <- paste0(Sgenes, " P-genes, unobserved")
                if (knowns < Sgenes) {
                    load(paste("~/Documents/unem_missing_unknowns", knowns, highnoise, complete, Sgenes, Egenes, nCells, paste(c(multi, ".rda"), collapse = ""), sep = "_"))
                } else {
                    load(paste("~/Documents/unem_missing", highnoise, complete, Sgenes, Egenes, nCells, paste(c(multi, ".rda"), collapse = ""), sep = "_"))
                }
            }
            ylab <- "perturbation profile correlation"
            ylim <- c(ymin2,1)
            if (dimnames(result)[[4]][s] == "time") {
                ylim <- NULL
            }
            if (dimnames(result)[[4]][s] != "cor") {
                ylab <- dimnames(result)[[4]][s]
            }
            if (s == 6) {
                ylab <- "fraction of identified null samples"
            }
            main <- paste0(main, ": ", lost[k])
            if (knowns < Sgenes) { main <- paste0(main, "\n unknowns: ", Sgenes-knowns) }
            databox <- NULL
            for (i in shownoise) {
                databox <- cbind(databox, result[,i,k,s,show])
            }
            ## databox <- cbind(result[,1,k,show2,show], result[,2,k,show2,show], result[,3,k,4,show])
            myboxplot(databox, col = cols, ylim = ylim, main = main, xlab = expression(sigma), ylab = ylab, box = box, scatter = scatter, dens = dens, xaxt = "n", border = cols, #notch = 1,
                      medcol = "black")
            axis(1, length(show)/2 + 0.5 + c(0, length(show), 2*length(show))[1:length(shownoise)], c(1,3,5)[shownoise])
            abline(v=c(length(show)*(1:(length(shownoise)-1))+0.5), col = 1)
            addgrid(x=c(-1,1,0.5), y=c(0,length(show)*3+1,1))
        }
    }
}
if (wrong == 0 & knowns >= Sgenes) {
    hist(rnorm(1000), border = "white", freq = 0, ylim = c(0,1), xlim = c(0,3), main = "", yaxt = "n", xlab = "", ylab = "", xaxt = "n")
    ## legend(0,1, c(expression(NEM~pi), "svm", "neural net", "random forest", "missForest", "random"), c("red", "blue", "darkgreen", "brown", "grey", "white"), box.col = "transparent", cex = 1, ncol = 3)
    legend(0,1, c(expression(NEM~pi), "svm", "neural net"), c("red", "blue", "darkgreen"), box.col = "transparent", cex = 1.5, ncol = 1)
    hist(rnorm(1000), border = "white", freq = 0, ylim = c(0,1), xlim = c(0,3), main = "", yaxt = "n", xlab = "", ylab = "", xaxt = "n")
    legend(0,1, c("random forest", "missForest", "mice"), c("brown", "orange", "pink"), box.col = "transparent", cex = 1.5, ncol = 1)
    hist(rnorm(1000), border = "white", freq = 0, ylim = c(0,1), xlim = c(0,3), main = "", yaxt = "n", xlab = "", ylab = "", xaxt = "n")
    legend(0,1, c("knn", "random"), c("turquoise", "grey"), box.col = "transparent", cex = 1.5, ncol = 1)
}
dev.off()

## nempi scheme:

D <- matrix(rnorm(10*20), 10)

rownames(D) <- paste0("E", 1:10)

Rho <- matrix(runif(5*20), 5)

Rho[which(Rho[, 5] != max(Rho[, 5]))[1], 5] <- max(Rho[, 5])

Rho <- apply(Rho, 2, function(x) return(x/sum(x)))

Rho[which(Rho[, 5] != max(Rho[, 5]))[1], 5] <- max(Rho[, 5])

phi <- matrix(c(1,0,0,0,0,
                1,0,0,0,0,
                0,1,0,0,0,
                1,0,0,0,0,
                0,0,0,1,0), 5)

rownames(Rho) <- colnames(phi) <- rownames(phi) <- paste0("P", 1:5)

colnames(D) <- colnames(Rho) <- paste0("S", 1:20)

Rhod <- apply(Rho, 2, function(x) { y <- x; y[which(y != max(x))] <- 0; y[which(y > 0)] <- 1; return(y) })

Rhod[, 11:20] <- 0

pdf("temp.pdf", width = 10, height = 4)
epiNEM::HeatmapOP(Rhod, col = "RdBu", colorkey = NULL, aspect = "iso", Rowv = 0, Colv = 0, cexRow = 2, cexCol = 2)
dev.off()

pdf("temp.pdf", width = 10, height = 8)
epiNEM::HeatmapOP(D, col = "RdBu", colorkey = NULL, aspect = "iso", Rowv = 0, Colv = 0, cexRow = 2, cexCol = 2)
dev.off()

pdf("temp.pdf", width = 5, height = 5)
plotDnf(phi, edgelwd = 3, fontsize = 10)
dev.off()

Rhos <- apply(Rho, 2, function(x) { y <- x; y[which(y != max(x))] <- 0; y[which(y > 0)] <- 1; return(y) })

Rhos <- Rhos + runif(length(Rhos), 0, 0.1)

Rhos <- apply(Rhos, 2, function(x) return(x/sum(x)))

pdf("temp.pdf", width = 10, height = 4)
epiNEM::HeatmapOP(Rhos, col = "RdBu", colorkey = NULL, aspect = "iso", Rowv = 0, Colv = 0, cexRow = 2, cexCol = 2)
dev.off()

## combinatorials:

library(Rgraphviz)

phi <- matrix(c(1,0,0,0,0,
                1,0,0,0,0,
                0,1,0,0,0,
                1,0,0,0,0,
                0,0,0,1,0), 5)

colnames(phi) <- rownames(phi) <- paste0("P", 1:5)

Rho <- matrix(0, 5, 10)

rownames(Rho) <- rownames(phi)

colnames(Rho) <- paste("S", 1:10, sep = "")

Rho[1, 1:3] <- 1
Rho[2, c(2,4:7)] <- 1
Rho[3, c(8:10, 10)] <- 1

pdf("temp.pdf", width = 5, height = 5)
p <- plotDnf(phi, edgelwd = 3, nodecol = list(P2 = "red", P3 = "red"), fontsize = 10)
epiNEM::HeatmapOP(Rho, col = "RdBu", Rowv = FALSE, aspect = "iso", colorkey = NULL, cexCol = 2, cexRow = 2)
Rho2 <- t(mnem:::mytc(phi))%*%Rho
Rho2[which(Rho2 > 1)] <- 1
epiNEM::HeatmapOP(Rho2, clusterx = Rho, col = "RdBu", Rowv = FALSE, aspect = "iso", colorkey = NULL, cexCol = 2, cexRow = 2)
dev.off()




















## warshall and other figures:

library(nanotime)
runs <- 1000
time <- matrix(0, runs, 3)
for (i in 1:runs) {
    cat(i)
    A <- simData(Sgenes = 10, Egenes = Egenes, mw = 1, nCells = nCells, Nem = 1, multi = multi)$Nem[[1]]
    A <- mytc(A)
    idx = which(A == 0)[1]
    nn <- dim(A)
    uv <- arrayInd(idx, nn)
    Phinew = A
    Phinew[idx[1]] = 1
    start <- nanotime(Sys.time())
    B <- mytc(Phinew, uv[1], uv[2])
    time[i, 1] <- as.numeric(nanotime(Sys.time()) - start)
    start <- nanotime(Sys.time())
    C <- transitive.closure(Phinew, mat = TRUE)
    time[i, 2] <- as.numeric(nanotime(Sys.time()) - start)
    start <- nanotime(Sys.time())
    D <- mytc(Phinew)
    time[i, 3] <- as.numeric(nanotime(Sys.time()) - start)
    if (any(B != C | C != D)) { stop() }
}

source("~/Documents/mnem/R/mnems_low.r")
pdf("temp.pdf", width = 6, height = 8)
myboxplot(time[1:i, ]/10^6, ylim = c(0.2,20), log = "y", xaxt = "n", main = "Wharshall vs Exponential", dens = 1, scatter = "", ylab = "seconds", lwd = 1, col = rgb(c(1,0,1), c(0,0,0.5), c(0,1,0), 0.5))
axis(1, 1:3, c("Warshall I.", "Exponential", "Warshall"))
dev.off()

microbenchmark(mytc(A), transitive.closure(A, mat = TRUE), time = 10)

Rho <- getRho(D)
pdf("temp.pdf", width = 20, height = 3)
epiNEM::HeatmapOP(Rho, col = "RdBu", cexRow = 0.75, cexCol = 0.75, rot = 0, aspect = "iso", colorkey = NULL)

epiNEM::HeatmapOP(ures2$Rho, col = "RdBu", cexRow = 0.75, cexCol = 0.75, rot = 0, aspect = "iso", clusterx = Rho, colorkey = NULL)

epiNEM::HeatmapOP(getRho(sim$data), col = "RdBu", cexRow = 0.75, cexCol = 0.75, rot = 0, aspect = "iso", clusterx = Rho, colorkey = NULL)

epiNEM::HeatmapOP(ures$Rho, col = "RdBu", cexRow = 0.75, cexCol = 0.75, rot = 0, aspect = "iso", clusterx = Rho, colorkey = NULL)

epiNEM::HeatmapOP(Rhosoft, col = "RdBu", cexRow = 0.75, cexCol = 0.75, rot = 0, aspect = "iso", clusterx = Rho, colorkey = NULL)
dev.off()

pdf("temp.pdf", width = 5, height = 5)
plotDnf(sim$Nem[[1]], edgelwd = 4, lwd = 4, fontsize = 14)
plotDnf(sim$Nem[[1]], edgelwd = 4, lwd = 4, fontsize = 14, nodecol = list("2" = "red", "9" = "red", "6" = "red"))
plotDnf(sim$Nem[[1]], edgelwd = 4, lwd = 4, fontsize = 14, nodecol = list("7" = "red", "4" = "red"))
dev.off()

pdf("temp.pdf", width = 9, height = 3)
par(mfrow=c(1,3))
pca <- rbind(cbind(rnorm(20, 0, 0.2), rnorm(20, 0, 0.2)),
             cbind(rnorm(30, 1, 0.2), rnorm(30, 1, 0.2)))
col <- c(sample(2:6, 15, replace = TRUE, prob = c(rep(0.1, 4), 0.6)), rep(1, 5),
                  sample(2:6, 20, replace = TRUE, prob = c(rep(0.1, 2), 0.6, rep(0.1, 2))), rep(1, 10))
plot(pca, col = col, xlab = "", ylab = "", main = "clusters")
col2 <- col
col2[which(col == 1)[1:5]] <- names(which.max(table(col[1:20])[-1]))
col2[which(col == 1)[6:15]] <- names(which.max(table(col[21:50])[-1]))
plot(pca, col = col2, xlab = "", ylab = "", main = "re-labeled")
col2[which(col != 6 & 1:50 <= 20)] <- names(which.max(table(col[1:20])[-1]))
col2[which(col != 4 & 1:50 > 20)] <- names(which.max(table(col[21:50])[-1]))
plot(pca, col = col2, xlab = "", ylab = "", main = "fully re-labeled")
dev.off()

pdf("temp.pdf", width = 8, height = 4)
par(mfrow=c(1,2))
plot(pca, col = col, xlab = "", ylab = "", main = "clusters")
d <- as.matrix(dist(pca))
col2 <- col
col2[which(col == 1)] <- col[which(col != 1)[apply(d[which(col == 1), which(col != 1)], 1, which.min)]]
plot(pca, col = col2, xlab = "", ylab = "", main = "re-labeled")
dev.off()

pdf("temp.pdf", width = 8, height = 4)
m <- matrix(c(c(1,1,0,0,1), c(0,0,1,1,0), c(1,1,1,1,1), rep(NA, 5*5), c(1,1,1,1,1)), 5)
rownames(m) <- 1:5
colnames(m) <- c(1,2,"1_2",rep("", 5),"sample")
epiNEM::HeatmapOP(m, Colv = 0, Rowv = 0, xrot = 0, colNA = "white", aspect = "iso", colorkey = NULL, col = "RdBu")
dev.off()

pdf("temp.pdf", width = 8, height = 4)
m <- matrix(c(c(1,1,0,0,1), c(0,0,1,1,0), c(0,1,1,0,0), rep(NA, 5*5), c(0,1,1,0,0)), 5)
rownames(m) <- 1:5
colnames(m) <- c(1,2,"1_2",rep("", 5),"sample")
epiNEM::HeatmapOP(m, Colv = 0, Rowv = 0, xrot = 0, colNA = "white", aspect = "iso", colorkey = NULL, col = "RdBu")
dev.off()

pdf("temp.pdf", width = 8, height = 4)
m <- matrix(c(c(1,1,0,0,1), c(1,1,1,1,1), c(1,1,1,1,1), rep(NA, 5*5), c(1,1,1,1,1)), 5)
rownames(m) <- 1:5
colnames(m) <- c(1,2,"1_2",rep("", 5),"sample")
epiNEM::HeatmapOP(m, Colv = 0, Rowv = 0, xrot = 0, colNA = "white", aspect = "iso", colorkey = NULL, col = "RdBu")
dev.off()







## try again the expected with respect to data creation: (does not seem to work.. how is it giving better results in unem?)

library(naturalsort)
library(nem)
library(cluster)
library(Rcpp)

source("~/Documents/mnem/R/mnems.r")
source("~/Documents/mnem/R/mnems_low.r")
sourceCpp("~/Documents/mnem/src/mm.cpp")

Sgenes <- 5
Egenes <- 10
samples <- 100

sim <- simData(Sgenes = Sgenes, Egenes = Egenes, mw = 1, nCells = Sgenes)

Rho <- matrix(0, Sgenes, samples)

Rho <- apply(Rho, 2, function(x) {
    y <- sample(seq_len(Sgenes), 1)
    x[y] <- 1
    return(x)
})

colnames(Rho) <- colnames(sim$data)[apply(Rho, 2, which.max)]

epiNEM::HeatmapOP(Rho)

ones <- which(Rho == 1)
zeros <- which(Rho == 0)

Rho[ones] <- runif(length(ones), 0.25, 1)

Rho[zeros] <- runif(length(zeros), 0, 0.75)

Rho <- apply(Rho, 2, function(x) return(x/sum(x)))

sdata <- sim$data[, naturalorder(colnames(sim$data))]%*%Rho

sdata <- (sdata - 0.5)/0.5

sdata <- sdata + rnorm(length(sdata), 0, 1)

epiNEM::HeatmapOP(sdata, col = "RdBu", cexRow = 0.75, cexCol = 0.75, bordercol = "transparent", rot = 0, dendrogram = "both")

D <- sdata

colnames(D) <- sample(seq_len(Sgenes), samples, replace = TRUE)

source("~/Documents/mnem/R/mnems_low.r")
res <- mynem(D, Rho = Rho, domean = FALSE)#, start = sim$Nem[[1]])

res$comp <- list()
res$comp[[1]] <- list()
res$comp[[1]]$phi <- res$adj

fitacc(res, sim)

Rho2 <- apply(Rho, 2, function(x) { y <- x*0; y[which.max(x)] <- 1; return(y) })

res <- mynem(D, Rho = Rho2, domean = FALSE)#, start = sim$Nem[[1]])

res$comp <- list()
res$comp[[1]] <- list()
res$comp[[1]]$phi <- res$adj

fitacc(res, sim)





res <- mnem(sdata, k = 2, search = "greedy", type = "random")

res <- mnem(sdata, k = 2, search = "greedy", type = "random")

fitacc(res, sim, strict = TRUE)

egenes <- numeric()

for (i in 1:5) {
    egenes <- c(egenes, sample(which(rownames(sdata) %in% i), 10))
}

sdataR <- sdata[egenes, ]

res <- mnem(sdataR, k = 2, search = "greedy", type = "random")

fitacc(res, sim, strict = TRUE)

## bugfixing:

#  D <- D_bkup; Gamma <- NULL


