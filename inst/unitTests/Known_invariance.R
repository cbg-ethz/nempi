sim <- mnem::simData(Sgenes=5,Nems=1,nCells=100)
D <- (sim$data-0.5)/0.5
D <- D + rnorm(length(D))
known <- sample(seq_len(ncol(D)),floor(0.5*ncol(D)))
colnames(D)[!seq_len(ncol(D)) %in% known] <- ""
result <- nempi(D,full=FALSE)
checkTrue(all(colnames(D)[known]==colnames(sim$data)[known]))
