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


likelihood <- function(x_val, obs_val, flip_p) {
  if (x_val == obs_val) {
    return(1 - flip_p)
  } else {
    return(flip_p)
  }
}



denoising_mh <- function(obs_im, niter, beta, flip_p) {
  x <- obs_im
  n <- nrow(x)
  m <- ncol(x)
  for (i in 1:(niter * n * m)) {
    rind <- sample(1:n, 1)
    cind <- sample(1:m, 1)
    x_val <- x[rind, cind]
    obs_val <- obs_im[rind, cind]

    # compute the engergy of the current x
    n0 <- xneig4(x, rind, cind, x[rind, cind])
    l0 <- likelihood(x_val, obs_val, flip_p)

    # compute the engergy if we flip the (sampl1[k],sampl2[l]) pixel
    n1 <- xneig4(x, rind, cind, 1 - x[rind, cind])
    l1 <- likelihood(1 - x_val, obs_val, flip_p)

    if (runif(1) < (exp(2 * beta * (n1 - n0)) * l1 / l0)) {
      x[rind, cind] <- 1 - x[rind, cind]
    }
    if (i %% (40 * n * m) == 0) {
      print(paste0("Iter: ", i / (n * m)))
    }
  }
  return(x)
}

set.seed(3)

niter <- 5000 # number of iterations in MH
beta <- 0.5



library(imager)
# read the noisy image and turn it into binary matrix
obs_im <- grayscale(load.image("A_noisy.png"))
obs_im <- as.matrix(round(obs_im / max(obs_im), 0))

# Denoising using M-H algorithm

flip_p <- 0.1
one_sample01 <- denoising_mh(obs_im, niter, beta, flip_p)
save.image(as.cimg(one_sample01), "denoised_A_01.png")

flip_p <- 0.05
one_sample005 <- denoising_mh(obs_im, niter, beta, flip_p)
save.image(as.cimg(one_sample005), "denoised_A_005.png")

flip_p <- 0.2
one_sample02 <- denoising_mh(obs_im, niter, beta, flip_p)
save.image(as.cimg(one_sample02), "denoised_A_02.png")
