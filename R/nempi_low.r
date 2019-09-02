#' @noRd
#' @importFrom stats cor
pifit <- function(x, y, D, unknown = "", balanced = FALSE, propagate = TRUE,
                  knowns = NULL) {
    Gamma <- getGamma(y$data)
    Gamma[which(Gamma > 1)] <- 1
    x$Gamma[which(is.na(x$Gamma) == TRUE)] <- 0
    Sgenes <- rownames(x$Gamma)
    kept <- Sgenes[which(Sgenes %in%
                                  unlist(strsplit(colnames(D), "_")))]
    miss <- which(!(Sgenes %in% kept))
    combi <- min(x$combi, nrow(Gamma) - 1)
    null <- x$null
    if (length(miss) > 0) {
        x$Gamma <- rbind(x$Gamma, matrix(0, length(miss), ncol(x$Gamma)))
        rownames(x$Gamma)[-seq_len(length(kept))] <- miss
        x$Gamma <- x$Gamma[naturalorder(rownames(x$Gamma)), ]
        A <- x$res$adj
        rownames(A) <- colnames(A) <- kept
        A <- cbind(A, matrix(0, nrow(A), length(miss)))
        A <- rbind(A, matrix(0, length(miss), ncol(A)))
        rownames(A)[-seq_len(length(kept))] <-
            colnames(A)[-seq_len(length(kept))] <- miss
        A <- A[naturalorder(rownames(A)), naturalorder(colnames(A))]
    }
    A <- mytc(x$res$adj)
    B <- mytc(y$Nem[[1]])
    ## Gammasoft <- rhosoft(y, D, combi = combi, null = null)
    Gammasoft <- apply(Gamma, 2, function(x) return(x/sum(x)))
    Gammasoft[which(is.nan(Gammasoft) == TRUE)] <- 0
    if (propagate) {
        x$Gamma <- t(A)%*%x$Gamma
    }
    Gammasoft <- t(B)%*%Gammasoft
    if (!is.null(knowns)) {
        Gammasoft <- Gammasoft[which(rownames(Gammasoft) %in% knowns), ]
    }
    glob <<- Gammasoft
    corres <- cor(as.vector(x$Gamma), as.vector(Gammasoft))
    x$Gamma <- apply(x$Gamma, 2, function(x) {
        y <- x*0
        y[which(x > 1/nrow(Gamma))] <- 1
        return(y)
    })
    Gamma <- t(B)%*%Gamma
    Gamma[which(Gamma > 1)] <- 1
    colnames(Gamma) <- colnames(x$Gamma) <- colnames(D)
    n <- nrow(A)
    if (!is.null(knowns)) {
        B <- B[which(rownames(B) %in% knowns), which(colnames(B) %in% knowns)]
        Gamma <- Gamma[which(rownames(Gamma) %in% knowns), ]
    }
    net <- (n*(n-1) - sum(abs(A - B)))/(n*(n-1))
    tp <- sum(x$Gamma[, which(colnames(x$Gamma) != unknown)] == 1 &
              Gamma[, which(colnames(x$Gamma) != unknown)] == 1)
    fp <- sum(x$Gamma[, which(colnames(x$Gamma) != unknown)] == 1 &
              Gamma[, which(colnames(x$Gamma) != unknown)] == 0)
    tn <- sum(x$Gamma[, which(colnames(x$Gamma) != unknown)] == 0 &
              Gamma[, which(colnames(x$Gamma) != unknown)] == 0)
    fn <- sum(x$Gamma[, which(colnames(x$Gamma) != unknown)] == 0 &
              Gamma[, which(colnames(x$Gamma) != unknown)] == 1)
    if (balanced) {
        known <- (((tp)/(tp+fn))+((tn)/(tn+fp)))/2
    } else {
        known <- (tp+tn)/(tp+tn+fp+fn)
    }
    tp <- sum(x$Gamma[, which(colnames(x$Gamma) == unknown)] == 1 &
              Gamma[, which(colnames(x$Gamma) == unknown)] == 1)
    fp <- sum(x$Gamma[, which(colnames(x$Gamma) == unknown)] == 1 &
              Gamma[, which(colnames(x$Gamma) == unknown)] == 0)
    tn <- sum(x$Gamma[, which(colnames(x$Gamma) == unknown)] == 0 &
              Gamma[, which(colnames(x$Gamma) == unknown)] == 0)
    fn <- sum(x$Gamma[, which(colnames(x$Gamma) == unknown)] == 0 &
              Gamma[, which(colnames(x$Gamma) == unknown)] == 1)
    if (balanced) {
        known2 <- (((tp)/(tp+fn))+((tn)/(tn+fp)))/2
    } else {
        known2 <- (tp+tn)/(tp+tn+fp+fn)
    }
    known <- c(known, known2)
    names(known) <- c("known", "unknown")
    return(list(net = net, known = known, cor = corres))
}
#' @noRd
#' @importFrom utils combn
getulods <- function(F, D, combi) {
    n <- nrow(F)
    if (combi == 1) {
        G <- F%*%D
        J <- G
        IF <- F
        HF <- diag(1, nrow(F))
    } else {
        if (combi == n) {
            H <- expand.grid(rep(list(c(0,1)), n))
        } else {
            Hfull <- list()
            for (i in seq_len(combi)) {
                H <- matrix(0, n, choose(n,i))
                combis <- combn(seq_len(n), i)
                for (j in seq_len(ncol(combis))) {
                    H[combis[, j], j] <- 1
                }
                Hfull[[i]] <- H
            }
            H <- do.call("cbind", Hfull)
        }
        I <- t(H)%*%F
        I[which(I > 1)] <- 1
        didx <- which(duplicated(apply(I, 1, paste, collapse = "")))
        didx <- didx[which(didx > n)]
        IF <- I
        HF <- H
        if (length(didx) > 0) {
            I <- I[-didx, ]
            H <- H[, -didx]
        }
        J <- I%*%D
        G <- lapply(seq_len(n), function(i) {
            x <- apply(J[which(H[i, ] == 1), , drop = FALSE], 2, max,
                       na.rm = TRUE)
            return(x)
        })
        G <- do.call("rbind", G)
        rownames(G) <- rownames(F)
    }
    return(list(G = G, J = J, I = IF, H = HF))
}
#' @noRd
rhosoft <- function(x, D, logtype = 2, combi = 1, null = TRUE) {
    phi <- mytc(x$Nem[[1]])
    theta <- theta2theta(x$theta[[1]], phi)
    F <- phi%*%theta
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
        Gamma <- logtype^G
        Gamma <- apply(Gamma, 2, function(x) return(x/sum(x)))
        Gamma <- Gamma[-1, ]
    } else {
        Gamma <- logtype^G
        Gamma <- apply(Gamma, 2, function(x) return(x/sum(x)))
    }
    return(Gamma)
}
#' @import mnem
#' @noRd
getGamma <- function(data) {

    Sgenes <- getSgenes(data)

    Rho <- matrix(0, length(Sgenes), ncol(data))

    for (i in seq_len(length(Sgenes))) {
        Rho[i, grep(paste0("^", Sgenes[i], "_|_", Sgenes[i],
                           "$|_", Sgenes[i], "_|^", Sgenes[i], "$"),
                    colnames(data))] <- 1
    }

    rownames(Rho) <- Sgenes
    colnames(Rho) <- colnames(data)
    Rho <- Rho[naturalsort(rownames(Rho)), ]
    return(Rho)
}

