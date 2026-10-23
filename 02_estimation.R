# Rolling estimation. All models share the same estimator.
suppressMessages(library(glmnet))
source("R/functions.R")
d <- readRDS("output/data.rds")
Y <- d$Y; id <- d$id; TT <- nrow(Y); N <- ncol(Y)

va <- list(l1 = 0.5, l2 = 0.5, nfolds = 10, intercept = TRUE, standardize = FALSE,
           s = "lambda.min", seed = 1234)
win <- 120; H <- 12; nf <- 101
p2 <- c(pi + 1e-5, pi / 3, 0); b2 <- c("short", "long")
p3 <- c(pi + 1e-5, pi / 3, pi / 12, 0); b3 <- c("short", "medium", "long")
ob <- c("noncore", "core"); oa <- rev(ob)

ls <- select_lag_order(Y, pmax = 12)
write.csv(ls, "output/lag_selection.csv", row.names = FALSE)
p <- unname(attr(ls, "selected")["BIC"])

YX <- cbind(d$X, Y); idX <- c(setNames(rep("controls", ncol(d$X)), colnames(d$X)), id)
YI <- cbind(d$X, ipmi = d$ipmi, Y); idI <- c(idX[1:5], ipmi = "controls", id)
g1 <- rep("all", N)

fit <- function(Z, p) do.call(estimate_var_aenet, c(list(Y = Z, p = p), va))
vd <- function(f, cl, o, h = H) vd_clustered(ma_coefficients(f$Phi, h), f$Sigma, cl, o)$theta
vdf <- function(f, cl, o, pt = p2, bn = b2)
  vd_clustered_freq(ma_coefficients(f$Phi, nf), f$Sigma, cl, o, pt, bn)$theta_bands

one <- function(w) {
  r <- (w - win + 1):w
  f <- fit(Y[r, ], p); f3 <- fit(Y[r, ], 3); fx <- fit(YX[r, ], p); fi <- fit(YI[r, ], p)
  list(cc = vd(f, id, ob), cc_alt = vd(f, id, oa), gc = vd(f, g1, "all"),
       fq_cc = vdf(f, id, ob), fq_gc = vdf(f, g1, "all"), fq3 = vdf(f, id, ob, p3, b3),
       p3 = vd(f3, id, ob),
       x = vd(fx, idX, c("controls", ob)), x_alt = vd(fx, idX, c("controls", oa)),
       fq_x = vdf(fx, idX, c("controls", ob)),
       xi = vd(fi, idI, c("controls", ob)),
       w150 = if (w >= 150) vd(fit(Y[(w - 149):w, ], p), id, ob))
}

nc <- if (.Platform$OS.type == "windows") 1 else max(1, parallel::detectCores())
ends <- win:TT
t0 <- Sys.time()
res <- parallel::mclapply(ends, one, mc.cores = nc)
print(Sys.time() - t0)

stk <- function(k) simplify2array(lapply(res, `[[`, k))
stkb <- function(k, bn) setNames(lapply(bn, function(b) simplify2array(lapply(res, function(z) z[[k]][[b]]))), bn)
out <- list(dates = d$dates[ends], p = p, id = id, idX = idX, idI = idI,
            cc = stk("cc"), cc_alt = stk("cc_alt"), gc = stk("gc"), p3 = stk("p3"),
            x = stk("x"), x_alt = stk("x_alt"), xi = stk("xi"),
            fq_cc = stkb("fq_cc", b2), fq_gc = stkb("fq_gc", b2), fq_x = stkb("fq_x", b2),
            fq3 = stkb("fq3", b3),
            w150 = simplify2array(lapply(res[ends >= 150], `[[`, "w150")),
            dates150 = d$dates[ends[ends >= 150]])
saveRDS(out, "output/rolling.rds")
