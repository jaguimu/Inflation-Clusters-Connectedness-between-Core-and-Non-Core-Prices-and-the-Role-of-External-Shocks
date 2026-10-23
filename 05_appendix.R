# Appendix figures: data (B), ordering (C), robustness (D)
source("R/functions.R"); source("R/plots.R")
dd <- readRDS("output/data.rds"); r <- readRDS("output/rolling.rds"); d <- r$dates
od <- "output/figures"; dir.create(od, recursive = TRUE, showWarnings = FALSE)
M <- function(a, id) t(apply(a, 3, conn_measures, id = id))
sk <- c("total", "within", "between"); c3 <- COL[c(1, 3, 2)]

# B.1-B.2
lab <- c(ALIM_NSUB = "Food", COMB_NSUB = "Fuels", TRAN_NSUB = "Transport", PUB_NSUB = "Public utilities",
         AYB_SUB = "Food and beverages", TEX_CALZ_SUB = "Textiles and footwear", ELEC_SUB = "Household appliances",
         INDUSTR_SUB = "Other industrial goods", CFH_SUB = "Food away from home", EDUC_SUB = "Education",
         SALUD_SUB = "Health", ALQ_SUB = "Rents", SERV_OTROS_SUB = "Other services",
         energy = "Energy prices (World Bank)", food = "Food prices (World Bank)",
         fert = "Fertilizer prices (World Bank)", tc = "Exchange rate (PEN/USD)", gscpi = "GSCPI (standard deviations)")
grid <- function(Z, nr, nc) function() {
  par(mfrow = c(nr, nc), mar = c(2.2, 3, 1.8, 0.8), mgp = c(2, 0.6, 0))
  yrs <- seq(as.Date("2000-01-01"), as.Date("2025-01-01"), by = "5 years")
  for (v in colnames(Z)) {
    plot(dd$dates, Z[, v], type = "n", xlab = "", ylab = "", main = lab[v], las = 1, bty = "l",
         cex.main = 0.95, font.main = 1, cex.axis = 0.8, xaxt = "n")
    axis(1, at = yrs, labels = format(yrs, "%Y"), cex.axis = 0.8)
    abline(h = 0, col = "grey70", lwd = 0.8); lines(dd$dates, Z[, v], col = COL[1], lwd = 0.9)
  }
}
ns <- names(dd$id)[dd$id == "noncore"]; cs <- names(dd$id)[dd$id == "core"]
save_fig(file.path(od, "figB1_noncore"), grid(dd$Y[, ns], 2, 2), 7, 4.4)
save_fig(file.path(od, "figB1_core"), grid(dd$Y[, cs], 3, 3), 7, 6.2)
save_fig(file.path(od, "figB2_controls"), grid(dd$X, 2, 3), 7, 4.4)

# C.1
m_av <- M((r$cc + r$cc_alt) / 2, r$id); m_gc <- M(r$gc, r$id)
fc1 <- function() {
  par(mfrow = c(2, 1), mar = c(2.3, 3, 2, 0.8), mgp = c(2, 0.6, 0))
  yl <- range(0, m_av[, sk], m_gc[, sk])
  for (k in list(list(m_av, "Clustered identification, average of the two cluster orderings"),
                 list(m_gc, "Generalized identification"))) {
    ts_frame(d, yl, k[[2]])
    for (j in 1:3) lines(d, k[[1]][, sk[j]], col = c3[j], lwd = c(1.8, 1.5, 1.5)[j])
    legend("topleft", c("Total", "Within clusters", "Between clusters"), col = c3, lwd = c(1.8, 1.5, 1.5),
           bty = "n", cex = 0.8, horiz = TRUE)
  }
}
save_fig(file.path(od, "figC1_avg_orderings"), fc1, 8, 6)

# C.3
D0 <- M(r$cc, r$id)[, "net_n"] - M(r$cc_alt, r$id)[, "net_n"]
D1 <- M(r$x, r$idX)[, "net_n"] - M(r$x_alt, r$idX)[, "net_n"]
fc3 <- function() {
  par(mar = c(2.3, 3, 1, 0.8), mgp = c(2, 0.6, 0))
  ts_frame(d, range(0, D0, D1), zero = TRUE)
  lines(d, D0, col = COL[3], lwd = 1.8); lines(d, D1, col = COL[2], lwd = 1.8)
  legend("topleft", c("Without controls", "With controls"), col = COL[3:2], lwd = 1.8, bty = "n", cex = 0.8, horiz = TRUE)
}
save_fig(file.path(od, "figC3_ordering_controls"), fc3, 8, 3.4)

# D.1-D.2
sp <- list(base = list(d = d, m = M(r$cc, r$id)), p3 = list(d = d, m = M(r$p3, r$id)),
           w150 = list(d = r$dates150, m = M(r$w150, r$id)))
sl <- c(base = "Baseline (p = 1, 120 months)", p3 = "p = 3", w150 = "Window of 150 months")
sc <- c(base = COL[1], p3 = COL[2], w150 = COL[3])
multi <- function(keys, titles, zero = FALSE) function() {
  par(mfrow = c(length(keys), 1), mar = c(2.3, 3, 2, 0.8), mgp = c(2, 0.6, 0))
  for (k in seq_along(keys)) {
    ts_frame(d, range(0, unlist(lapply(sp, function(s) s$m[, keys[k]]))), titles[k], zero)
    for (s in names(sp)) lines(sp[[s]]$d, sp[[s]]$m[, keys[k]], col = sc[s], lwd = if (s == "base") 1.8 else 1.4)
    if (k == 1) legend("topleft", sl, col = sc, lwd = c(1.8, 1.4, 1.4), bty = "n", cex = 0.8, horiz = TRUE)
  }
}
save_fig(file.path(od, "figD1_rob_system"), multi(sk, c("Total connectedness", "Within clusters", "Between clusters")), 8, 7.5, 150)
save_fig(file.path(od, "figD2_rob_net"), multi(c("net_n", "net_c"), c("Non-core cluster", "Core cluster"), TRUE), 8, 5.5, 150)

# D.3
b3 <- c(short = "Short run (2-6 months)", medium = "Medium run (6-24 months)", long = "Long run (more than 24 months)")
F3 <- lapply(r$fq3, M, id = r$id)
fd3 <- function() {
  par(mfrow = c(1, 3), mar = c(2.3, 3, 2, 0.8), mgp = c(2, 0.6, 0))
  yl <- range(0, sapply(F3, function(m) m[, c("net_n", "net_c")]))
  for (b in names(b3)) {
    ts_frame(d, yl, b3[b], zero = TRUE)
    lines(d, F3[[b]][, "net_n"], col = COL[2], lwd = 1.6); lines(d, F3[[b]][, "net_c"], col = COL[1], lwd = 1.6)
    if (b == "short") legend("topleft", c("Non-core", "Core"), col = COL[2:1], lwd = 1.6, bty = "n", cex = 0.85)
  }
}
save_fig(file.path(od, "figD3_rob_freq3"), fd3, 10, 3.4, 150)

# D.4
n0 <- sp$base$m[, "net_n"]; n5 <- M(r$x, r$idX)[, "net_n"]; n6 <- M(r$xi, r$idI)[, "net_n"]
fd4 <- function() {
  par(mar = c(2.3, 3, 1, 0.8), mgp = c(2, 0.6, 0))
  ts_frame(d, range(0, n0, n5, n6), zero = TRUE)
  lines(d, n0, col = COL[3], lwd = 1.8); lines(d, n5, col = COL[2], lwd = 1.4, lty = 2); lines(d, n6, col = COL[2], lwd = 1.8)
  legend("topleft", c("Without controls", "Baseline controls", "Baseline controls + imported wholesale prices"),
         col = COL[c(3, 2, 2)], lwd = c(1.8, 1.4, 1.8), lty = c(1, 2, 1), bty = "n", cex = 0.8)
}
save_fig(file.path(od, "figD4_rob_ipmi"), fd4, 8, 3.6, 150)
