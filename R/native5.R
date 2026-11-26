# Natively implemented mySORT (Chen et al. 2018, BMC Bioinformatics; SVR() of
# svr_without_rmse_corr.r in the lsbnb/mysort_2022 image), replacing e1071 and
# preprocessCore. Same conventions as native.R: X is genes x samples, results
# are samples x features, each function self-contained (see bd_isolate()).

# Quantile normalisation across the columns of a matrix, as
# preprocessCore::normalize.quantiles(): every column gets the distribution given
# by the row means of the sorted columns; tied values share the mean of the
# quantiles of their ranks.
bd_native$quantile_normalize <- function(X) {
  X <- data.matrix(X)
  target <- rowMeans(apply(X, 2L, sort))
  out <- X
  for (j in seq_len(ncol(X))) {
    r <- rank(X[, j], ties.method = "average")
    out[, j] <- (target[floor(r)] + target[ceiling(r)]) / 2
  }
  out
}

# Linear nu-support vector regression, solved in the primal. libsvm (behind
# e1071::svm) solves the dual with SMO to a tolerance of 0.001; here the same convex
# problem
#   min 1/2 |w|^2 + C n nu eps + C sum(xi + xi*)
#   s.t. |y_i - w.x_i - b| <= eps + xi_i (one-sided slacks xi, xi* >= 0), eps >= 0
# is solved with a log-barrier method. The slacks have a closed form given
# (w, b, eps), so Newton's method runs on the r + 2 variables (w, b, eps) only and
# one solve takes milliseconds. Returns w (the intercept is not needed by mySORT).
bd_native$nusvr_linear_w <- function(S, y, nu, C = 1, mu_end = 1e-13, shrink = 0.2) {
  n <- nrow(S); r <- ncol(S)
  cc <- C * n * nu
  # Slack at the barrier optimum for a margin distance p: the positive root of
  # C xi^2 + (C p - 2 mu) xi - mu p = 0, i.e. xi = (2 mu - C p + d) / (2 C) with
  # d = sqrt(C^2 p^2 + 4 mu^2). Both xi and p + xi are written without
  # cancellation. Returns list(xi, pxi = p + xi).
  slack <- function(p, mu) {
    d <- sqrt((C * p)^2 + 4 * mu^2)
    pos <- p > 0
    xi <- ifelse(pos, (2 * mu + 4 * mu^2 / (C * p + d)) / (2 * C), (2 * mu - C * p + d) / (2 * C))
    pxi <- ifelse(pos, (2 * mu + C * p + d) / (2 * C), (2 * mu + 4 * mu^2 / (d - C * p)) / (2 * C))
    list(xi = xi, pxi = pxi)
  }
  parts <- function(v, mu) {
    w <- v[seq_len(r)]; b <- v[r + 1L]; eps <- v[r + 2L]
    f <- drop(S %*% w) + b
    p <- eps - y + f
    q <- eps + y - f
    sp <- slack(p, mu); sq <- slack(q, mu)
    hv <- function(s) C * s$xi - mu * log(s$pxi) - mu * log(s$xi)
    list(sp = sp, sq = sq,
         phi = 0.5 * sum(w^2) + cc * eps - mu * log(eps) + sum(hv(sp)) + sum(hv(sq)))
  }
  v <- c(rep(0, r), 0, 1)
  mu <- 1
  repeat {
    for (it in seq_len(60L)) {
      pr <- parts(v, mu)
      hp <- -mu / pr$sp$pxi; hq <- -mu / pr$sq$pxi
      Hp <- mu / (pr$sp$pxi^2 + pr$sp$xi^2); Hq <- mu / (pr$sq$pxi^2 + pr$sq$xi^2)
      w <- v[seq_len(r)]; eps <- v[r + 2L]
      g <- c(w + crossprod(S, hp - hq), sum(hp - hq), cc - mu / eps + sum(hp + hq))
      # Hessian = sum_i Hp_i a_i a_i' + Hq_i c_i c_i' + diag(1_r, 0, mu / eps^2)
      A <- cbind(S, 1, 1); Cm <- cbind(-S, -1, 1)
      H <- crossprod(A * sqrt(Hp), A * sqrt(Hp)) + crossprod(Cm * sqrt(Hq), Cm * sqrt(Hq))
      H[cbind(seq_len(r), seq_len(r))] <- H[cbind(seq_len(r), seq_len(r))] + 1
      H[r + 2L, r + 2L] <- H[r + 2L, r + 2L] + mu / eps^2
      step <- tryCatch(-solve(H, g), error = function(e) -solve(H + diag(1e-10, r + 2L), g))
      dec <- -sum(g * step)
      if (dec < 1e-12 * (1 + abs(pr$phi))) break
      t <- 1
      if (v[r + 2L] + t * step[r + 2L] <= 0) t <- 0.99 * v[r + 2L] / -step[r + 2L]
      repeat {
        new <- v + t * step
        ph <- parts(new, mu)$phi
        if (is.finite(ph) && ph <= pr$phi - 1e-4 * t * dec) break
        t <- t / 2
        if (t < 1e-14) break
      }
      v <- new
    }
    if (mu <= mu_end) break
    mu <- mu * shrink
  }
  v[seq_len(r)]
}

# mySORT: quantile-normalised mixture, standardised signature, linear nu-SVR with
# nu chosen among 0.25, 0.5, 0.75 by RMSE, negative weights set to 0, weights
# scaled to 100 (percent).
bd_native$mysort_fit <- function(X, S) {
  Y <- quantile_normalize(as.matrix(X))
  dimnames(Y) <- dimnames(X)
  genes <- intersect(rownames(S), rownames(Y))
  S <- S[genes, , drop = FALSE]; Y <- Y[genes, , drop = FALSE]
  S <- (S - mean(S)) / stats::sd(as.vector(S))
  fit_w <- function(y, nu) {
    w <- nusvr_linear_w(S, y, nu)
    w[w < 0] <- 0
    w / sum(w)
  }
  out <- t(vapply(seq_len(ncol(Y)), function(i) {
    y <- Y[, i]; y <- (y - mean(y)) / stats::sd(y)
    nus <- c(0.25, 0.5, 0.75)
    ws <- lapply(nus, function(nu) fit_w(y, nu))
    rmse <- vapply(ws, function(w) sqrt(mean((drop(S %*% w) - y)^2)), 0)
    ws[[which.min(rmse)]] * 100
  }, numeric(ncol(S))))
  dimnames(out) <- list(colnames(X), colnames(S))
  out
}
