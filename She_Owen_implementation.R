
#################################################
### function definition
#################################################
library(quantreg)
set.seed(2023)

# soft thresholding
soft.threshold <- function (x, lam = 1) ifelse(abs(x) < lam, 0, x - sign(x)*lam)

# hard thresholding
hard.threshold <- function (x, lam = 1) ifelse(abs(x) < lam, 0, x)

# thresholding wrapper
threshold <- function (x, lam = 1, type = "soft") {
  if (type == "soft") res <- soft.threshold(x, lam)
  if (type == "hard") res <- hard.threshold(x, lam)
  return( res )
}

# algorithm 1 (She and Owen)
ipod.1 <- function (x, y, b = NULL, type = "soft", lam = 1) {
  
  X <- cbind(1, x)
  decomp <- qr(X)
  Q <- qr.Q(decomp)
  R <- qr.R(decomp)
  
  ## initialization
  if (is.null(b)) b <- rq(y ~ x, tau = 0.5, data = data.frame(x = x, y = y))$coefficients
  
  ## loop
  j <- 0; dif <- Inf
  gam <- y - drop(X%*%b)
  repeat {
    y.adj <- y - gam
    b <- drop(solve(R) %*% t(Q) %*% y.adj)
    r <- y - drop(X%*%b)
    gam.new <- threshold(r, lam = lam, type = type)
    if (max(abs(gam.new - gam)) < 1e-5 | j == 100) break
    j <- j + 1
    gam <- gam.new
  }
  
  return(list(b = b, gam = gam))
  
}


# ipod algorithm wrapper
ipod <- function (x, y, algorithm = 1, b = NULL, type = "soft", lam = 1) {
  if (algorithm == 1) res <- ipod.1(x, y, b = b, type = type, lam = lam)
  if (algorithm == 2) res <- ipod.2(x, y, b = b, type = type, lam = lam)
  return(res)
}


# BIC selection
BIC <- function (x, y, lambda = NULL, algorithm = 1, type = "soft") {
  X <- cbind(1, x)
  H <- X %*% solve(t(X)%*%X) %*% t(X)
  eig <- eigen(H)
  c.idx <- which(abs(eig$values) < 1e-8)
  U.c <- eig$vectors[,c.idx]
  y.tilde <- drop(t(U.c)%*%y)
  m <- nrow(X) - ncol(X)
  
  init <- rq(y ~ x, tau = 0.5, data = data.frame(x = x, y = y))
  resid <- init$residuals
  sig2 <- mad(resid)^2
  
  if (is.null(lambda)) {
    lam.max <- max(abs(drop((diag(nrow(X))-H)%*%y) / sqrt(1-diag(H))))
    #lambda <- seq(lam.max, 0.5, by = -0.1) ## She and Owen's suggestion
    lambda <- 2^seq(log2(5), -4, length = 50)
  }
  
  score <- beta.mat <- gam.mat <- NULL
  for (a in seq_along(lambda)) {
    res <- ipod(x, y, algorithm = algorithm, type = type, lam = lambda[a])
    nzr <- sum(res$gam != 0)
    k <- nzr + 1
    #bic1 <- m * log(sum((y.tilde - drop(t(U.c)%*%res$gam))^2)/m) + k * (log(m) + 1)
    bic2 <- n * log(sig2) + sum((y.tilde - drop(t(U.c)%*%res$gam))^2) / sig2 + (nzr + 2) * log(n)
    score <- rbind(score, c(lam = lambda[a], nzr = nzr, bic1 = bic1, bic2 = bic2))
    beta.mat <- rbind(beta.mat, c(lam = lambda[a], res$b))
    gam.mat <- rbind(gam.mat, c(lam = lambda[a], res$gam))
  }
  colnames(beta.mat) <- c("lam", paste0("b", 1:ncol(X) - 1))
  colnames(gam.mat) <- c("lam", paste0("obs", 1:nrow(X)))
  
  return(list(score = score, beta.mat = beta.mat, gam.mat = gam.mat))
}




#################################################
# test with an example
#################################################

library(tidyverse)

n <- 100
p <- 1
beta <- c(-1, 1)

x <- rnorm(n)
#x <- sort(rnorm(n))
y <- beta[1] + beta[2]*x + rnorm(n)*0.5
out.idx <- 1:10
out <- rep(1, n); out[out.idx] <- 2
y[out.idx] <- y[out.idx] + rnorm(length(out.idx), 2, 0.3)
#y[out.idx] <- rnorm(length(out.idx), -2, 1)
dat <- data.frame(x = x, y = y, out = factor(out))
lse.est <- lm(y ~ x, data = dat)$coefficients
lad.est <- rq(y ~ x, tau = 0.5, data = dat)$coefficients
res1.soft <- ipod(x, y, algorithm = 1, type = "soft", lam = 1)
(res1.soft.b <- res1.soft$b)
res1.soft.gam <- res1.soft$gam
res1.hard <- ipod(x, y, algorithm = 1, type = "hard", lam = 1)
(res1.hard.b <- res1.hard$b)
res1.hard.gam <- res1.hard$gam

#method <- c("true", "lse", "lad", "soft-ipot", "hard-ipot")
#df <- data.frame(method = factor(1:length(method), labels = method), rbind(beta, lse.est, lad.est, res1.soft.b, res1.hard.b))
method <- c("true", "least sqaures", "soft-ipot", "hard-ipot")
df <- data.frame(method = factor(1:length(method), labels = method), rbind(beta, lse.est, res1.soft.b, res1.hard.b))
names(df) <- c("method", "b0", "b1")

ggplot(dat, aes(x = x, y = y, shape = out)) +
  geom_abline(data = df, aes(intercept = b0, slope = b1, color = method, linetype = method), size = 1) +
  geom_point() + scale_shape_manual(values = c(21, 19))
  scale_color_brewer(palette = "Set1")



## BIC selection

res1.soft <- BIC(x, y, lambda = NULL, algorithm = 1, type = "soft")
res1.hard <- BIC(x, y, lambda = NULL, algorithm = 1, type = "hard")
#res2.soft <- BIC(x, y, lambda = NULL, algorithm = 2, type = "soft")
#res2.hard <- BIC(x, y, lambda = NULL, algorithm = 2, type = "hard")


## bic curve - She and Owen's suggestion
df <- data.frame(lam = res1.soft$score[,"lam"], soft = res1.soft$score[,"bic1"], hard = res1.hard$score[,"bic1"]) %>%
  pivot_longer(col = -lam, names_to = "penalty", values_to = "bic") %>%
  arrange(lam)
ggplot(df, aes(x = lam, y = bic, color = penalty)) + geom_line() + scale_x_log10()


## bic curve - my suggestion
lam.opt.soft <- res1.soft$score[which.min(res1.soft$score[,"bic2"]), "lam"]
lam.opt.hard <- res1.hard$score[which.min(res1.hard$score[,"bic2"]), "lam"]
df <- data.frame(lam = res1.soft$score[,"lam"], soft = res1.soft$score[,"bic2"], hard = res1.hard$score[,"bic2"]) %>%
  pivot_longer(col = -lam, names_to = "penalty", values_to = "bic") %>%
  arrange(lam)
ggplot(df, aes(x = lam, y = bic, color = penalty)) + geom_line() + scale_x_log10() +
  geom_vline(xintercept = c(lam.opt.soft, lam.opt.hard), linetype = "dashed")



## gamma solution path - soft
df <- as.data.frame(res1.soft$gam.mat) %>% 
  pivot_longer(col = -lam, names_to = "gamma", values_to = "estimate") %>%
  arrange(lam)
ggplot(df, aes(x = lam, y = estimate, color = gamma)) + geom_line() + scale_x_log10() + guides(color = FALSE) +
  geom_vline(xintercept = lam.opt.soft)


## gamma solution path - hard
df <- as.data.frame(res1.hard$gam.mat) %>% 
  pivot_longer(col = -lam, names_to = "gamma", values_to = "estimate") %>%
  arrange(lam)
ggplot(df, aes(x = lam, y = estimate, color = gamma)) + geom_line() + scale_x_log10() + guides(color = FALSE) +
  geom_vline(xintercept = lam.opt.hard)


## beta solution path - soft
df <- as.data.frame(res1.soft$beta.mat) %>% 
  pivot_longer(col = -lam, names_to = "beta", values_to = "estimate") %>%
  arrange(lam)
ggplot(df, aes(x = lam, y = estimate, color = beta)) + geom_line() + scale_x_log10() +
  geom_vline(xintercept = lam.opt.soft)


## beta solution path - hard
df <- as.data.frame(res1.hard$beta.mat) %>% 
  pivot_longer(col = -lam, names_to = "beta", values_to = "estimate") %>%
  arrange(lam)
ggplot(df, aes(x = lam, y = estimate, color = beta)) + geom_line() + scale_x_log10() +
  geom_vline(xintercept = lam.opt.hard)
 



### refit with an optimal lambda

res1.soft <- ipod(x, y, algorithm = 1, type = "soft", lam = lam.opt.soft)
res1.soft$b
res1.soft$gam

res1.hard <- ipod(x, y, algorithm = 1, type = "hard", lam = lam.opt.hard)
res1.hard$b
res1.hard$gam


# result : soft
df <- subset(dat, res1.soft$gam != 0)
ggplot(dat, aes(x = x, y = y, shape = out)) + geom_point() + geom_abline(intercept = beta[1], slope = beta[2], size = 1) +
  geom_abline(intercept = res1.soft$b[1], slope = res1.soft$b[2], size = 1, linetype = "dashed", color = "red") +
  geom_point(data = df, aes(x = x, y = y), shape = 21, size = 4, color = "red") + scale_shape_manual(values = c(21, 19))
  
# result : hard
df <- subset(dat, res1.hard$gam != 0)
ggplot(dat, aes(x = x, y = y, shape = out)) + geom_point() + geom_abline(intercept = beta[1], slope = beta[2], size = 1) +
  geom_abline(intercept = res1.hard$b[1], slope = res1.hard$b[2], size = 1, linetype = "dashed", color = "blue") +
  geom_point(data = df, aes(x = x, y = y), shape = 21, size = 4, color = "blue") + scale_shape_manual(values = c(21, 19))

