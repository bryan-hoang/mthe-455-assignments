xneig4 <- function(x, a, b, val) { # compute the local energy at (a,b)
  n <- dim(x)[1]
  m <- dim(x)[2]
  nei <- rep(0, 4)
  if (a != 1) {
    nei[1] <- (x[a - 1, b] == val)
  } # Left
  if (a != n) {
    nei[2] <- (x[a + 1, b] == val)
  } # Right
  if (b != 1) {
    nei[3] <- (x[a, b - 1] == val)
  } # North
  if (b != m) {
    nei[4] <- (x[a, b + 1] == val)
  } # South

  return(sum(nei)) # Sum
}

isinghm <- function(niter, n, m = n, beta) {
  # generate X_1
  int_x <- sample(c(0, 1), n * m, prob = c(0.5, 0.5), rep = TRUE)
  x <- matrix(int_x, n, m)


  # generate X_t
  for (i in 1:(niter * n * m)) {
    rind <- sample(1:n, 1)
    cind <- sample(1:m, 1)

    # compute the engergy of the current x
    n0 <- xneig4(x, rind, cind, x[rind, cind])

    # compute the engergy if we flip the (sampl1[k],sampl2[l]) pixel
    n1 <- xneig4(x, rind, cind, 1 - x[rind, cind])

    if (runif(1) < exp(2 * beta * (n1 - n0))) { # if accept
      x[rind, cind] <- 1 - x[rind, cind]
    }
    if (i %% (40 * n * m) == 0) {
      print(paste0("Iter: ", i / (n * m)))
    }
  }
  return(x)
}

set.seed(7)
n <- 100 # number of rows, as well columns
niter <- 5000 # number of iterations in MH

# try different beta
for (beta in c(0.2, 0.3, 0.4, 0.5, 0.6)) {
  gen_img <- isinghm(niter, n, n, beta)

  # save image
  png(filename = paste0("genImage_", beta, ".png"))
  image(1:n, 1:n, gen_img, col = gray.colors(2))
  dev.off()
}
