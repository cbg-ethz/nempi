#' @noRd
auc <- function(a,b) {
    n <- length(a)
    c <- sum((a-c(0,a[-n]))*(b+c(0,b[-n]))/c(1,rep(2,n-1)))
    return(c)
}
#' @noRd
#' @import mnem
getSgenes <- mnem:::getSgenes
#' @noRd
#' @import mnem
modData <- mnem:::modData
#' @noRd
#' @import mnem
mynem <- mnem:::mynem
#' @noRd
#' @import mnem
scoreAdj <- mnem:::scoreAdj
#' @noRd
#' @import mnem
mytc <- mnem:::mytc
#' @noRd
#' @import mnem
theta2theta <- mnem:::theta2theta
#' @noRd
#' @import mnem
getSgeneN <- mnem:::getSgeneN
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

