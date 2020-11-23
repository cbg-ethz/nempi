sim <- mnem::simData(Sgenes=5,Nems=1)
D <- (sim$data-0.5)/0.5
D <- D + rnorm(length(D))
result <- nempi(D,phi=sim$Nem[[1]])
checkTrue(all(result$res$adj==sim$Nem[[1]]))
