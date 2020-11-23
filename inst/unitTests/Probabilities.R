D <- matrix(rnorm(1000*100), 1000, 100)
colnames(D) <- sample(seq_len(5), 100, replace = TRUE)
result <- nempi(D)
checkEquals(result$Gamma, result$probs[-1,], checkNames = TRUE)
