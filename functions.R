# Functions for clustered connectedness (time and frequency domain)

# lagged regressors
build_var_data <- function(Y, p) {
  Y <- as.matrix(Y)
  TT <- nrow(Y); N <- ncol(Y)
  if (TT <= p) stop("Sample too short for p = ", p)
  y <- Y[(p + 1):TT, , drop = FALSE]
  X <- do.call(cbind, lapply(seq_len(p), function(l) Y[(p + 1 - l):(TT - l), , drop = FALSE]))
  vn <- colnames(Y); if (is.null(vn)) vn <- paste0("x", seq_len(N))
  colnames(X) <- as.vector(outer(vn, paste0(".L", seq_len(p)), paste0))
  colnames(y) <- vn
  list(y = y, X = X)
}

# AIC/BIC/HQ for VAR(p), OLS
select_lag_order <- function(Y, pmax = 12, intercept = TRUE) {
  Y <- as.matrix(Y); N <- ncol(Y)
  Teff <- nrow(Y) - pmax
  res <- data.frame(p = 1:pmax, AIC = NA_real_, BIC = NA_real_, HQ = NA_real_)
  for (p in 1:pmax) {
    d <- build_var_data(Y[(pmax - p + 1):nrow(Y), , drop = FALSE], p)
    X <- if (intercept) cbind(1, d$X) else d$X
    if (ncol(X) >= nrow(X)) next
    B <- qr.solve(X, d$y)
    U <- d$y - X %*% B
    S <- crossprod(U) / Teff
    k <- ncol(X) * N
    ld <- as.numeric(determinant(S, logarithm = TRUE)$modulus)
    res$AIC[p] <- ld + 2 * k / Teff
    res$BIC[p] <- ld + log(Teff) * k / Teff
    res$HQ[p]  <- ld + 2 * log(log(Teff)) * k / Teff
  }
  attr(res, "selected") <- c(AIC = which.min(res$AIC), BIC = which.min(res$BIC),
                              HQ = which.min(res$HQ))
  res
}

# one adaptive elastic net equation
aenet_fit <- function(X, y, w, l1 = 0.5, l2 = 0.5, intercept = TRUE,
                      foldid = NULL, nfolds = 10, s = "lambda.min", lambda = NULL) {
  sy <- if (intercept) sqrt(mean((y - mean(y))^2)) else sqrt(mean(y^2))
  if (sy == 0) sy <- 1
  ys <- y / sy
  alpha <- l1 / (l1 + 2 * l2 * sy)
  if (is.null(lambda)) {
    if (is.null(foldid)) foldid <- sample(rep(seq_len(nfolds), length.out = length(y)))
    f <- glmnet::cv.glmnet(X, ys, alpha = alpha, penalty.factor = w, foldid = foldid,
                           intercept = intercept, standardize = FALSE)
    cf <- as.numeric(coef(f, s = s)); lam <- f[[s]]
  } else {
    f <- glmnet::glmnet(X, ys, alpha = alpha, lambda = lambda, penalty.factor = w,
                        intercept = intercept, standardize = FALSE, thresh = 1e-14)
    cf <- as.numeric(coef(f)); lam <- lambda
  }
  list(a = cf[1] * sy, b = cf[-1] * sy, lambda_glmnet = lam, alpha = alpha, sy = sy)
}

# VAR by adaptive elastic net, equation by equation
estimate_var_aenet <- function(Y, p, l1 = 0.5, l2 = 0.5, nfolds = 10, intercept = TRUE,
                               standardize = FALSE, s = "lambda.min",
                               seed = 1234, max_weight = 1e6) {
  if (!requireNamespace("glmnet", quietly = TRUE)) stop("Install package 'glmnet'.")
  Y <- as.matrix(Y)
  N <- ncol(Y)
  d <- build_var_data(Y, p)
  y <- d$y; X <- d$X
  Teff <- nrow(X); K <- ncol(X)

  mx <- colMeans(X); sx <- apply(X, 2, sd)
  sx[sx == 0] <- 1
  Xs <- if (standardize) sweep(sweep(X, 2, mx), 2, sx, "/") else X

  set.seed(seed)
  foldid <- sample(rep(seq_len(nfolds), length.out = Teff))

  coef_mat <- matrix(0, K, N, dimnames = list(colnames(X), colnames(y)))
  const <- setNames(numeric(N), colnames(y))
  lambda_sel <- setNames(numeric(N), colnames(y))

  for (i in seq_len(N)) {
    yi <- y[, i]
    if (Teff > K + intercept) {
      Xo <- if (intercept) cbind(1, Xs) else Xs
      b0 <- qr.solve(Xo, yi)
      if (intercept) b0 <- b0[-1]
    } else {
      cvr <- glmnet::cv.glmnet(Xs, yi, alpha = 0, foldid = foldid,
                               intercept = intercept, standardize = FALSE)
      b0 <- as.numeric(coef(cvr, s = "lambda.min"))[-1]
    }
    w <- pmin(1 / abs(b0), max_weight)

    f <- aenet_fit(Xs, yi, w, l1, l2, intercept, foldid = foldid, s = s)
    b <- f$b; a <- f$a
    if (standardize) { b <- b / sx; a <- a - sum(b * mx) }
    coef_mat[, i] <- b
    const[i] <- a
    lambda_sel[i] <- f$lambda_glmnet
  }

  U <- y - sweep(X %*% coef_mat, 2, const, "+")
  Sigma <- crossprod(U) / Teff

  Phi <- lapply(seq_len(p), function(l) t(coef_mat[((l - 1) * N + 1):(l * N), , drop = FALSE]))

  list(Phi = Phi, const = const, Sigma = Sigma, resid = U, coef = coef_mat,
       lambda = lambda_sel, p = p, N = N, Teff = Teff, varnames = colnames(y))
}

# MA coefficients A_0..A_{H-1}
ma_coefficients <- function(Phi, H) {
  N <- nrow(Phi[[1]]); P <- length(Phi)
  A <- vector("list", H)
  A[[1]] <- diag(N)
  if (H > 1) for (h in 2:H) {
    Ah <- matrix(0, N, N)
    for (l in seq_len(P)) if (h - l >= 1) Ah <- Ah + Phi[[l]] %*% A[[h - l]]
    A[[h]] <- Ah
  }
  A
}

# block-recursive Q, Omega block diagonal
cluster_Q <- function(Sigma, cluster_id, order = unique(cluster_id)) {
  N <- nrow(Sigma)
  stopifnot(length(cluster_id) == N, setequal(order, unique(cluster_id)))
  Qinv <- diag(N)
  prev <- integer(0)
  for (cl in order) {
    idx <- which(cluster_id == cl)
    if (length(prev) > 0) {
      Qinv[idx, prev] <- -Sigma[idx, prev, drop = FALSE] %*%
        solve(Sigma[prev, prev, drop = FALSE])
    }
    prev <- c(prev, idx)
  }
  Omega <- Qinv %*% Sigma %*% t(Qinv)
  Omega <- (Omega + t(Omega)) / 2
  same <- outer(cluster_id, cluster_id, "==")
  max_offblock <- if (any(!same)) max(abs(Omega[!same])) else 0
  Omega[!same] <- 0
  list(Qinv = Qinv, Q = solve(Qinv), Omega = Omega, max_offblock = max_offblock)
}

# clustered variance decomposition (row-normalized)
vd_clustered <- function(A, Sigma, cluster_id, order = unique(cluster_id)) {
  N <- nrow(Sigma)
  qo <- cluster_Q(Sigma, cluster_id, order)
  QO <- qo$Q %*% qo$Omega
  num <- matrix(0, N, N); den <- numeric(N)
  for (h in seq_along(A)) {
    num <- num + (A[[h]] %*% QO)^2
    den <- den + diag(A[[h]] %*% Sigma %*% t(A[[h]]))
  }
  theta_raw <- sweep(num, 2, diag(qo$Omega), "/") / den
  theta <- theta_raw / rowSums(theta_raw)
  dimnames(theta_raw) <- dimnames(theta) <- dimnames(Sigma)
  list(theta = theta, theta_raw = theta_raw, Omega = qo$Omega, Qinv = qo$Qinv)
}

# cluster net: (to - from) / N_k
cluster_net <- function(theta, cluster_id) {
  cl <- unique(cluster_id)
  sapply(cl, function(c) {
    ic <- cluster_id == c
    to <- sum(theta[!ic, ic]); from <- sum(theta[ic, !ic])
    100 * (to - from) / sum(ic)
  })
}

# X-13 seasonal adjustment
seasonal_adjust <- function(Y, dates, ...) {
  if (!requireNamespace("seasonal", quietly = TRUE)) stop("Install package 'seasonal'.")
  Y <- as.matrix(Y)
  st <- c(as.numeric(format(dates[1], "%Y")), as.numeric(format(dates[1], "%m")))
  Ysa <- Y
  diag <- data.frame(series = colnames(Y), adjusted = FALSE, qs_pval_original = NA_real_,
                     qs_pval_adjusted = NA_real_, note = "", stringsAsFactors = FALSE)
  for (j in seq_len(ncol(Y))) {
    m <- tryCatch(seasonal::seas(stats::ts(Y[, j], frequency = 12, start = st), ...),
                  error = function(e) e)
    if (inherits(m, "error")) { diag$note[j] <- conditionMessage(m); next }
    Ysa[, j] <- as.numeric(seasonal::final(m))
    diag$adjusted[j] <- TRUE
    q <- tryCatch(seasonal::qs(m), error = function(e) NULL)
    if (!is.null(q)) {
      if ("qsori" %in% rownames(q)) diag$qs_pval_original[j] <- q["qsori", "p-val"]
      if ("qssadj" %in% rownames(q)) diag$qs_pval_adjusted[j] <- q["qssadj", "p-val"]
    }
  }
  list(Y = Ysa, diag = diag)
}

# DFT indices in [down, up)
freq_indices <- function(n, up, down) {
  space <- (0:floor(n / 2)) / floor(n / 2) * pi
  out <- (space >= down & space < up) * 1
  out <- if (n %% 2 == 0) c(out, rev(out[2:(length(out) - 1)])) else c(out, rev(out[2:length(out)]))
  which(out == 1)
}

freq_partition <- function(partition, n, names = NULL) {
  stopifnot(all(diff(partition) < 0))
  part <- lapply(seq_len(length(partition) - 1),
                 function(i) freq_indices(n, partition[i], partition[i + 1]))
  if (any(lengths(part) == 0)) stop("Empty band; increase n.")
  if (length(unique(unlist(part))) != n) warning("Partition does not cover all frequencies.")
  if (is.null(names)) {
    per <- round(2 * pi / partition, 1)
    names <- sapply(seq_along(part), function(i) sprintf("%s-%s", per[i], per[i + 1]))
  }
  setNames(part, names)
}

# clustered VD by frequency band
vd_clustered_freq <- function(A, Sigma, cluster_id, order = unique(cluster_id),
                              partition = c(pi + 1e-5, pi / 3, 0),
                              band_names = c("short", "long")) {
  N <- nrow(Sigma); n <- length(A)
  qo <- cluster_Q(Sigma, cluster_id, order)
  QO <- qo$Q %*% qo$Omega; om <- diag(qo$Omega)
  bands <- freq_partition(partition, n, band_names)
  rng <- sort(unique(unlist(bands)))
  Amat <- t(sapply(A, as.vector))
  Fh <- stats::mvfft(Amat)
  Psi <- lapply(seq_len(n), function(s) matrix(Fh[s, ], N, N))
  den <- Reduce(`+`, lapply(Psi[rng], function(P) Re(diag(P %*% Sigma %*% t(Conj(P)))))) / n
  th_s <- lapply(Psi, function(P) sweep(Mod(P %*% QO)^2 / n, 2, om, "/") / den)
  tot <- rowSums(Reduce(`+`, th_s[rng]))
  theta_bands <- lapply(bands, function(b) {
    m <- Reduce(`+`, th_s[b]) / tot; dimnames(m) <- dimnames(Sigma); m })
  theta_total <- Reduce(`+`, theta_bands)
  share <- sapply(theta_bands, rowSums)
  list(theta_total = theta_total, theta_bands = theta_bands, band_share = share,
       bands = bands, partition = partition)
}

# summary measures (%), rows of theta = receivers
conn_measures <- function(th, id) {
  th <- 100 * th; dg <- diag(th)
  iC <- id == "core"; iN <- id == "noncore"; iX <- id == "controls"; n <- sum(iC | iN)
  part <- function(i, k) c(own = mean(dg[i]), comove = (sum(th[i, i]) - sum(dg[i])) / sum(i),
                           contag = sum(th[i, k]) / sum(i))
  w <- (sum(th[iC, iC]) - sum(dg[iC]) + sum(th[iN, iN]) - sum(dg[iN])) / n
  b <- (sum(th[iC, iN]) + sum(th[iN, iC])) / n
  x <- if (any(iX)) sum(th[!iX, iX]) / n else 0
  c(setNames(part(iC, iN), c("c_own", "c_comove", "c_contag")),
    setNames(part(iN, iC), c("n_own", "n_comove", "n_contag")),
    total = w + b + x, within = w, between = b, from_ctrl = x,
    net_n = (sum(th[iC, iN]) - sum(th[iN, iC])) / sum(iN),
    net_c = (sum(th[iN, iC]) - sum(th[iC, iN])) / sum(iC))
}

# component net, inflation block only
node_net <- function(th) { m <- 100 * th; diag(m) <- 0; colSums(m) - rowSums(m) }
