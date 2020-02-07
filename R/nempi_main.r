#' Main function for NEM based perturbation imputation.
#'
#' Infers perturbations profiles based on a sparse perturbation
#' matrix and differential gene expression as log odds
#' @param D either a binary effects matrix or log odds matrix as
#' for Nested Effects Models (see package 'nem')
#' @param unknown colname of samples without mutation data, E.g. ""
#' @param Gamma matrix with expectations of perturbations, e.g. if you
#' have a binary mutation matrix, just normalize the columns to have sum 1
#' @param type "null": does not use the unknown samples for inference
#' at the start, "random" uses them in a random fashion (not recommended)
#' @param full if FALSE, does not change the known profiles
#' @param verbose if TRUE gives more output during inference
#' @param logtype log type for the log odds
#' @param null if FALSE does not use a NULL node for uninformative samples
#' @param soft if FALSE discretizes Gamma during the inference
#' @param combi if combi > 1, uses a more complex algorithm to infer
#' combinatorial perturbations (experimental)
#' @param converged the absolute difference of log likelihood till convergence
#' @param complete if TRUE uses the complete-data
#' logliklihood (recommended for many E-genes)
#' @param mw if NULL infers mixture weights, otherwise keeps them fixed
#' @param max_iter maximum iterations of the EM algorithm
#' @param keepphi if TRUE, uses the previous phi for the next inference,
#' if FALSE always starts with start network (and empty and full)
#' @param start starting network as adjacency matrix
#' @param phi if not NULL uses only this phi and does not infer a new one
#' @param ... additional parameters for the nem
#' function (see package mnem, function mnem or mnem:::mynem)
#' @return nempi object
#' @author Martin Pirkl
#' @export
#' @import mnem naturalsort nem matrixStats
#' @examples
#' D <- matrix(rnorm(1000*100), 1000, 100)
#' colnames(D) <- sample(seq_len(5), 100, replace = TRUE)
#' Gamma <- matrix(sample(c(0,1), 5*100, replace = TRUE, p = c(0.9, 0.1)), 5,
#' 100)
#' Gamma <- apply(Gamma, 2, function(x) return(x/sum(x)))
#' Gamma[is.na(Gamma)] <- 0
#' rownames(Gamma) <- seq_len(5)
#' result <- nempi(D, Gamma = Gamma)
nempi <- function(D, unknown = "", Gamma = NULL, type = "null", full = TRUE,
                 verbose = FALSE, logtype = 2, null = TRUE, soft = TRUE,
                 combi = 1, converged = 0.1, complete = TRUE,
                 mw = NULL, max_iter = 100, keepphi = TRUE, start = NULL,
                 phi = NULL, ...) {
    samplenames <- colnames(D)
    if (is.null(Gamma)) {
        realnames <- getSgenes(D)
    } else {
        realnames <- rownames(Gamma)
    }
    realnames <- naturalsort(realnames)
    D <- modData(D)
    colnames(D) <- gsub("\\..*", "", colnames(D))
    Sgenes <- getSgenes(D)
    n <- length(Sgenes)
    if (is.null(Gamma)) {
        Gamma <- getGamma(D)
        if (type %in% "random") {
            Gamma[cbind(sample(seq_len(n), sum(colnames(Gamma) %in% unknown),
                             replace = TRUE),
                      which(colnames(Gamma) %in% unknown))] <- 1
        } else if (type %in% "null") {
            Gamma[, which(colnames(Gamma) %in% unknown)] <- 0
        }
    }
    if (!full & sum(colnames(D) %in% unknown) == 0) {
        unk <- which(apply(Gamma, 2, sum) == 0)
        if (length(unk) > 0) {
            colnames(D)[unk] <- unknown
            colnames(Gamma) <- colnames(D)
        } else {
            stop(paste(c("No unlabeled samples available. ",
                         "Either use full inference, change ",
                         "column names or Gamma."), collapse = ""))
        }
    }
    if (type %in% "random") {
        Gamma[cbind(sample(seq_len(n), sum(apply(Gamma, 2, sum) == 0),
                         replace = TRUE),
                  which(apply(Gamma, 2, sum) == 0))] <- 1
    }
    if (soft) {
        Gamma <- apply(Gamma, 2, function(x) return(x/sum(x)))
        Gamma[which(is.na(Gamma))] <- 0
    }
    if (combi >= n) {
        combi <- n - 1
    }
    llold <- ll <- -Inf
    time0 <- TRUE
    stop <- FALSE
    Gammaold <- Gamma
    phiold <- thetaold <- 0
    lls <- numeric()
    if (is.null(mw)) {
        if (null) { nn <- n+1 } else { nn <- n }
        pi <- piold <- rep(1/nn, nn)
    } else {
        pi <- piold <- mw
    }
    olls <- numeric()
    evoGamma <- evophi <- evotheta <- evopi <- numeric()
    count <- 0
    if (is.null(start)) {
        start <- matrix(1, n, n)
        rownames(start) <- colnames(start) <- seq_len(n)
    }
    while(!stop | time0) {
        count <- count + 1
        miss <- which(apply(Gamma, 1, sum) == 0)
        if (length(miss) >= length(Sgenes)-1) {
            res <- list()
            res$adj <- matrix(0, length(Sgenes), length(Sgenes))
            colnames(res$adj) <- rownames(res$adj) <- Sgenes
            break()
        } else {
            if (is.null(phi)) {
                res0 <- mynem(D, modified = TRUE, Rho = Gamma,
                                     verbose = verbose,
                              logtype = logtype, start = NULL, ...)
                res <- mynem(D, modified = TRUE, Rho = Gamma,
                                    verbose = verbose,
                             logtype = logtype, start = start, ...)
                ones <- start*0 + 1
                res1 <- mynem(D, modified = TRUE, Rho = Gamma,
                                     verbose = verbose,
                              logtype = logtype, start = ones, ...)
                if (res0$score > res$score) { res <- res0 }
                if (res1$score > res$score) { res <- res1 }
                if (keepphi) {
                    start <- res$adj
                }
            } else {
                res <- scoreAdj(D, phi, Rho = Gamma, dotopo = TRUE)
                res$adj <- phi
            }
        }
        evophi <- c(evophi, sum(abs(mytc(res$adj) - phiold)))
        phiold <- mytc(res$adj)
        evotheta <- c(evotheta, sum(thetaold != res$subtopo))
        thetaold <- res$subtopo
        theta <- theta2theta(res$subtopo, res$adj)
        F <- mytc(res$adj)%*%theta
        S <- getulods(F, D, combi)
        single <- 1
        if (single) {
            G <- S$G
        } else {
            G <- S$I%*%D
            H <- S$H
        }
        if (null) {
            G <- rbind(0, G)
            sub <- 1
        } else {
            sub <- 0
        }
        expG <- function(G, complete = FALSE) {
            maxlog <- 1023
            if (complete & any(G > maxlog)) {
                maxnum <- 2^maxlog
                shrinkfac <- log(maxnum)/log(logtype)
                G <- apply(G, 2, function(x) {
                    xmax <- max(x)
                    if (xmax > maxlog) {
                        x <- x - (xmax - shrinkfac/length(x))
                    }
                    return(x)
                })
            }
            Z <- logtype^G
            return(Z)
        }
        if (any(pi == 0) & complete) {
            pi[which(pi == 0)] <- 1
        }
        if (!full) {
            G <- G[, which(colnames(Gamma) %in% unknown), drop = FALSE]
        }
        if (complete) {
            Z <- expG(G, complete = complete)
            Z <- Z/colSums(Z)[col(Z)]
            ll <- sum(apply(Z*(G + log(pi)/log(logtype)), 2, sum))
            probs <- Z
        } else {
            probs <- expG(G, complete = complete)
            probs <- probs*pi
            ll <- sum(log(apply(probs, 2, sum))/log(logtype))
            probs <- apply(probs, 2, function(x) return(x/sum(x)))
        }
        if (is.null(mw)) {
            pi <- apply(probs, 1, sum)
            pi <- pi/sum(pi)
        }
        ## ll <- res$score
        evopi <- c(evopi, sum(abs(sort(pi) - sort(piold))))
        piold <- pi
        if (soft) { G2 <- probs } else { G2 <- G }
        if (!full) {
            Gamma[, which(colnames(Gamma) %in% unknown)] <- 0
            if (soft) {
                if (null) {
                    G2 <- G2[-1, ]
                }
                Gamma[, which(colnames(Gamma) %in% unknown)] <- G2
            } else {
                idx <- as.list(apply(G2, 2, function(x)
                    return(which(x == max(x, na.rm = TRUE)) - sub)))
                aidx <- cbind(unlist(idx),
                          unlist(lapply(seq_len(length(idx)),
                                        function(x, idx) {
                                            y <- rep(x, length(idx[[x]]))
                                            return(y) }, idx)))
                aidx <- aidx[which(aidx[, 2] %in%
                                   which(colnames(Gamma) %in% unknown)), ]
                Gamma[aidx] <- 1
            }
        } else {
            Gamma <- Gamma*0
            if (soft) {
                if (null) {
                    G2 <- G2[-1, ]
                }
                Gamma <- G2
            } else {
                idx <- as.list(apply(G2, 2, function(x)
                    return(which(x == max(x, na.rm = TRUE)) - sub)))
                aidx <- cbind(unlist(idx),
                          unlist(lapply(seq_len(length(idx)),
                                        function(x, idx) {
                                            y <- rep(x, length(idx[[x]]))
                                            return(y) }, idx)))
                Gamma[aidx] <- 1
            }
        }
        if (!single) {
            Gamma <- lapply(seq_len(n), function(i) {
                x <- apply(Gamma[which(H[i, ] == 1), , drop = FALSE], 2, max,
                           na.rm = TRUE)
                return(x)
            })
            Gamma <- do.call("rbind", Gamma)
            rownames(Gamma) <- rownames(F)
            Gamma <- apply(Gamma, 2, function(x) return(x/sum(x)))
        }
        evoGamma <- c(evoGamma, sum(abs(Gamma - Gammaold)))
        if ((evoGamma[length(evoGamma)] == 0 & evophi[length(evophi)] == 0 &
             evotheta[length(evotheta)] == 0 &evopi[length(evopi)] == 0 &
             !time0)
            | abs(ll - llold) <= converged
            | count == max_iter | !(ll > llold)) {
            stop <- TRUE
        }
        llold <- ll
        Gammaold <- Gamma
        lls <- c(lls, ll)
        if (verbose) { print(ll) }
        time0 <- FALSE
    }
    if (null) { pi <- pi[-1] }
    names(pi) <- rownames(Gamma)
    rownames(Gamma) <- realnames
    colnames(Gamma) <- samplenames
    ures <- list(res = res, Gamma = Gamma, lls = lls, probs = probs,
                 full = full,
                 pi = pi, evoGamma = evoGamma, evophi = evophi,
                 evotheta = evotheta,
                 evopi = evopi, combi = combi, null = null)
    class(ures) <- "nempi"
    return(ures)
}
#' Bootstrapping function
#'
#' Bootstrap algorithm to get a more stable result.
#' @param D either a binary effects matrix or log odds matrix as
#' @param bsruns number of bootstraps
#' @param bssize number of E-genes for each boostrap
#' @param replace if TRUE, actual bootstrap, if False sub-sampling
#' @param ... additional parameters for the function nempi
#' @return list with aggregate Gamma and aggregate causal network phi
#' @author Martin Pirkl
#' @export
#' @examples
#' D <- matrix(rnorm(1000*100), 1000, 100)
#' colnames(D) <- sample(seq_len(5), 100, replace = TRUE)
#' Gamma <- matrix(sample(c(0,1), 5*100, replace = TRUE, p = c(0.9, 0.1)), 5,
#' 100)
#' Gamma <- apply(Gamma, 2, function(x) return(x/sum(x)))
#' Gamma[is.na(Gamma)] <- 0
#' rownames(Gamma) <- seq_len(5)
#' result <- nempibs(D, bsruns = 3, Gamma = Gamma)
nempibs <- function(D, bsruns = 100, bssize = 0.5, replace = TRUE, ...) {
    n <- getSgeneN(D)
    Sgenes <- getSgenes(D)
    aggGamma <- matrix(0, n, ncol(D))
    aggPhi <- matrix(0, n, n)
    rownames(aggGamma) <- rownames(aggPhi) <- colnames(aggPhi) <- Sgenes
    colnames(aggGamma) <- colnames(D)
    for (i in seq_len(bsruns)) {
        Dr <- D[sample(seq_len(nrow(D)), ceiling(nrow(D)*bssize),
                       replace = replace), ]
        tmp <- nempi(Dr, ...)
        aggGamma <- aggGamma + tmp$Gamma
        aggPhi <- aggPhi + mytc(tmp$res$adj)
    }
    return(list(Gamma = aggGamma, phi = aggPhi))
}
#' Plot convergence
#'
#' Produces different convergence plots based on a nempi object
#' @param x nempi object
#' @param ... additional parameters for plot
#' @return plot
#' @author Martin Pirkl
#' @export
#' @import graphics mnem
#' @examples
#' D <- matrix(rnorm(1000*100), 1000, 100)
#' colnames(D) <- sample(seq_len(5), 100, replace = TRUE)
#' Gamma <- matrix(sample(c(0,1), 5*100, replace = TRUE, p = c(0.9, 0.1)), 5,
#' 100)
#' Gamma <- apply(Gamma, 2, function(x) return(x/sum(x)))
#' Gamma[is.na(Gamma)] <- 0
#' rownames(Gamma) <- seq_len(5)
#' result <- nempi(D, Gamma = Gamma)
#' par(mfrow=c(2,3))
#' plotConvergence(result)
plotConvergence <- function(x, ...) {
    plot(x$lls, main = "log odds evolution", ylab = "log odds",
         xlab = "iterations", ...)
    hist(x$probs, main = "sample probabilities", xlab = "probabilities")
    plot(x$evoGamma, main = expression(evolution ~ of ~ Gamma),
         ylab = "sum of absolute distance", xlab = "iterations", ...)
    plot(x$evopi, main= expression(evolution ~ of ~ pi),
         ylab = "sum of absolute distance", xlab = "iterations", ...)
    plot(x$evotheta, main= expression(evolution ~ of ~ theta),
         ylab = "hamming distance", xlab = "iterations", ...)
    plot(x$evophi, main = expression(evolution ~ of ~ phi),
         ylab = "hamming distance", xlab = "iterations", ...)
}
#' Classification
#'
#' Builds and uses different classifiers to infer perturbation profiles
#' @param D either a binary effects matrix or log odds matrix as
#' for Nested Effects Models (see package 'nem')
#' @param unknown colname of samples without mutation data, E.g. ""
#' @param full if FALSE, does not change the known profiles
#' @param method either one of svm, nn, rf
#' @param size parameter for neural network (see package 'nnet')
#' @param MaxNWts parameters for neural network (see package 'nnet')
#' @param ... additional parameters for mnem:::mynem, see mnem::mnem
#' @return plot
#' @author Martin Pirkl
#' @export
#' @import e1071 nnet randomForest mnem nem
#' @importFrom stats predict
#' @examples
#' D <- matrix(rnorm(1000*100), 1000, 100)
#' colnames(D) <- sample(seq_len(5), 100, replace = TRUE)
#' Gamma <- matrix(sample(c(0,1), 5*100, replace = TRUE, p = c(0.9, 0.1)), 5,
#' 100)
#' Gamma <- apply(Gamma, 2, function(x) return(x/sum(x)))
#' Gamma[is.na(Gamma)] <- 0
#' rownames(Gamma) <- seq_len(5)
#' result <- classpi(D)
classpi <- function(D, unknown = "", full = TRUE,
                    method = "svm", size = NULL, MaxNWts = 10000, ...) {
    samplenames <- colnames(D)
    realnames <- getSgenes(D)
    realnames <- naturalsort(realnames)
    if (is.null(size)) { size <- length(realnames) }
    D_bkup <- D
    D <- modData(D)
    DK <- D[, which(colnames(D) != unknown)]
    if (!full) {
        Gammafull <- matrix(0, length(realnames)*2,
                            sum(colnames(D) %in% unknown))
    } else {
        Gammafull <- matrix(0, length(realnames)*2, ncol(D))
    }
    n <- length(realnames)
    count <- 0
    for (Sgene in as.character(seq_len(length(realnames)))) {
        labels <- rep(0, ncol(DK))
        labels[grep(paste0("^", Sgene, "$|_", Sgene, "$|^",
                           Sgene, "_|_", Sgene, "_"), colnames(DK))] <- Sgene
        train <- as.data.frame(t(DK))
        colnames(train) <- paste0("Var", seq_len(ncol(train)))
        train <- cbind(train, label = labels)
        if (method %in% "svm") {
            sres <- svm(label ~ ., data = train, probability = TRUE)
        }
        if (method %in% "nnet") {
            sres <- nnet(label ~ ., data = train, size = size, trace = FALSE,
                         MaxNWts = MaxNWts)
        }
        if (method %in% "randomForest") {
            sres <- randomForest(label ~ ., data = train)
        }
        DU <- D[, which(colnames(D) == unknown), drop = FALSE]
        test <- as.data.frame(t(DU))
        colnames(test) <- paste0("Var", seq_len(ncol(test)))
        traintest <- as.data.frame(t(D))
        colnames(traintest) <- paste0("Var", seq_len(ncol(traintest)))
        ll <- -Inf
        llold <- -Inf
        lls <- ll
        start <- TRUE
        unident <- FALSE
        test2 <- test
        while (ll - llold > 0 | start) {
            start <- FALSE
            llold <- ll
            if (!full) {
                if (method %in% "svm") {
                    p <- predict(sres, test2, probability = TRUE)
                    p <- attr(p, "probabilities")
                    p <- p[, naturalsort(colnames(p)), drop = FALSE]
                }
                if (method %in% "nnet") {
                    p <- predict(sres, test2)
                }
                if (method %in% "randomForest") {
                    p <- predict(sres, test2,
                                 type = "prob")
                    p <- p[, naturalsort(colnames(p)), drop = FALSE]
                }
                Gamma <- t(p)
                colnames(Gamma) <- rep("", ncol(Gamma))
            } else {
                if (method %in% "svm") {
                    p <- predict(sres, traintest,
                                 probability = TRUE)
                    p <- attr(p, "probabilities")
                    p <- p[, naturalsort(colnames(p)), drop = FALSE]
                }
                if (method %in% "nnet") {
                    p <-  predict(sres, traintest)
                }
                if (method %in% "randomForest") {
                    p <- predict(sres, traintest,
                                 type = "prob")
                    p <- p[, naturalsort(colnames(p)), drop = FALSE]
                }
                Gamma <- t(p)
                colnames(Gamma) <- colnames(D)
            }
            if (method %in% "nnet") {
                Gamma <- rbind("0" = (1 - Gamma), Sgene = Gamma)
            }
            Gamma[which(Gamma == 0)] <- 10^-323
            ll <- sum(apply(log(Gamma), 2, max))
            llold <- ll
            lls <- c(lls, ll)
        }
        count <- count + 1
        Gammafull[count, ] <- Gamma[2, ]
        Gammafull[count+length(realnames), ] <- Gamma[1, ]
    }
    Gamma <- Gammafull
    Gamma <- apply(Gamma, 2, function(x) return(x/sum(x)))
    Gamma <- Gamma[seq_len(length(realnames)), ]
    if (!full) {
        GammaD <- getGamma(D)
        GammaD[, which(colnames(D) %in% unknown)] <- Gamma
        Gamma <- GammaD
    }
    rownames(Gamma) <- realnames
    colnames(Gamma) <- samplenames
    if (unident) {
        res <- list()
        res$adj <- diag(1, n)
        colnames(res$adj) <- rownames(res$adj) <- seq_len(n)
    } else {
        res <- mynem(D_bkup, Rho = Gamma, ...)
    }
    ures <- list(res = res, Gamma = Gamma, probs = Gamma, full = full,
                 null = FALSE, ll = ll, lls = lls)
    return(ures)
}
