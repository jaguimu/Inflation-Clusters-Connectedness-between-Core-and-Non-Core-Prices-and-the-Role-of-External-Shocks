# Figures 1-8
suppressMessages(library(igraph))
source("R/functions.R"); source("R/plots.R")
r <- readRDS("output/rolling.rds"); d <- r$dates
od <- "output/figures"; dir.create(od, recursive = TRUE, showWarnings = FALSE)
M <- function(a, id) t(apply(a, 3, conn_measures, id = id))
isN <- r$id == "noncore"; ia <- d < as.Date("2022-01-01"); ib <- !ia
avgw <- function(a, i) { m <- 100 * apply(a[, , i, drop = FALSE], c(1, 2), mean); dimnames(m) <- dimnames(a)[1:2]; m }
sk <- c("total", "within", "between"); c3 <- COL[c(1, 3, 2)]
m_cc <- M(r$cc, r$id); m_gc <- M(r$gc, r$id)

# Fig 1
f1 <- function() {
  par(mfrow = c(2, 1), mar = c(2.3, 3, 2, 0.8), mgp = c(2, 0.6, 0))
  yl <- range(0, m_cc[, sk], m_gc[, sk])
  for (k in list(list(m_cc, "Clustered identification"), list(m_gc, "Generalized identification"))) {
    ts_frame(d, yl, k[[2]])
    for (j in 1:3) lines(d, k[[1]][, sk[j]], col = c3[j], lwd = c(1.8, 1.5, 1.5)[j])
    legend("topleft", c("Total", "Within clusters", "Between clusters"), col = c3, lwd = c(1.8, 1.5, 1.5),
           bty = "n", cex = 0.8, horiz = TRUE)
  }
}
# Fig 2
f2 <- function() {
  par(mfrow = c(2, 1), mar = c(2.3, 3, 2, 0.8), mgp = c(2, 0.6, 0))
  yl <- range(m_cc[, c("net_n", "net_c")], m_gc[, c("net_n", "net_c")])
  for (k in list(list(m_cc, "Clustered identification"), list(m_gc, "Generalized identification"))) {
    ts_frame(d, yl, k[[2]], zero = TRUE)
    lines(d, k[[1]][, "net_n"], col = COL[2], lwd = 1.8); lines(d, k[[1]][, "net_c"], col = COL[1], lwd = 1.8)
    legend("topleft", c("Non-core", "Core"), col = COL[2:1], lwd = 1.8, bty = "n", cex = 0.8, horiz = TRUE)
  }
}
save_fig(file.path(od, "fig1_system"), f1, 8, 6)
save_fig(file.path(od, "fig2_net"), f2, 8, 6)

# Fig 3 (and C.2)
L <- net_layout(avgw(r$cc, rep(TRUE, length(d))))
TH <- list(a = avgw(r$cc, ia), b = avgw(r$cc, ib), ga = avgw(r$gc, ia), gb = avgw(r$gc, ib))
sc <- net_scale(TH)
fnet <- function(k1, k2) function() {
  layout(matrix(c(1, 2, 3, 3), 2, byrow = TRUE), heights = c(1, 0.08)); par(mar = c(0.5, 0.5, 2, 0.5))
  net_panel(TH[[k1]], isN, L, sc, "Windows ending in 2009-2021")
  net_panel(TH[[k2]], isN, L, sc, "Windows ending in 2022-2025"); net_legend()
}
save_fig(file.path(od, "fig3_network"), fnet("a", "b"), 12, 6)
save_fig(file.path(od, "figC2_network_gc"), fnet("ga", "gb"), 12, 6)
write.csv(round(sapply(TH, function(m) colSums(offd(m)) - rowSums(offd(m))), 1), file.path(od, "net_components.csv"))

# Figs 4-5
bands <- c(short = "Short run (2-6 months)", long = "Long run (more than 6 months)")
FQ <- list(cc = lapply(r$fq_cc, M, id = r$id), gc = lapply(r$fq_gc, M, id = r$id))
idl <- c(cc = "Clustered", gc = "Generalized")
f4 <- function() {
  par(mfrow = c(2, 2), mar = c(2.3, 3, 2, 0.8), mgp = c(2, 0.6, 0))
  yl <- range(0, unlist(lapply(FQ, function(z) lapply(z, function(m) m[, sk]))))
  for (i in names(FQ)) for (b in names(bands)) {
    m <- FQ[[i]][[b]]; ts_frame(d, yl, paste0(idl[i], ": ", bands[b]), by = "4 years", cex.main = 0.9)
    for (j in 1:3) lines(d, m[, sk[j]], col = c3[j], lwd = c(1.8, 1.5, 1.5)[j])
    if (i == "cc" && b == "short") legend("topleft", c("Total", "Within clusters", "Between clusters"),
                                          col = c3, lwd = 1.6, bty = "n", cex = 0.75)
  }
}
f5 <- function() {
  par(mfrow = c(2, 2), mar = c(2.3, 3, 2, 0.8), mgp = c(2, 0.6, 0))
  yl <- range(0, unlist(lapply(FQ, function(z) lapply(z, function(m) m[, c("net_n", "net_c")]))))
  for (i in names(FQ)) for (b in names(bands)) {
    m <- FQ[[i]][[b]]; ts_frame(d, yl, paste0(idl[i], ": ", bands[b]), zero = TRUE, by = "4 years", cex.main = 0.9)
    lines(d, m[, "net_n"], col = COL[2], lwd = 1.8); lines(d, m[, "net_c"], col = COL[1], lwd = 1.8)
    if (i == "cc" && b == "short") legend("topleft", c("Non-core", "Core"), col = COL[2:1], lwd = 1.8, bty = "n", cex = 0.75)
  }
}
save_fig(file.path(od, "fig4_freq_system"), f4, 9, 6)
save_fig(file.path(od, "fig5_freq_net"), f5, 9, 6)

# Fig 6
TB <- list(sa = avgw(r$fq_cc$short, ia), sb = avgw(r$fq_cc$short, ib), la = avgw(r$fq_cc$long, ia), lb = avgw(r$fq_cc$long, ib))
sb <- net_scale(TB)
f6 <- function() {
  layout(matrix(c(1, 2, 3, 4, 5, 5), 3, byrow = TRUE), heights = c(1, 1, 0.07)); par(mar = c(0.3, 0.3, 1.6, 0.3))
  tt <- c(sa = "Short run, windows ending in 2009-2021", sb = "Short run, windows ending in 2022-2025",
          la = "Long run, windows ending in 2009-2021", lb = "Long run, windows ending in 2022-2025")
  for (k in names(TB)) net_panel(TB[[k]], isN, L, sb, tt[k], lcex = 0.58, asz = 0.3, tline = 0.3, tcex = 0.95)
  net_legend()
}
save_fig(file.path(od, "fig6_freq_network"), f6, 12, 11.5, res = 170)

# Figs 7-8
m_x <- M(r$x, r$idX); FX <- lapply(r$fq_x, M, id = r$idX)
c0 <- COL[3]; c1 <- COL[2]
f7 <- function() {
  par(mar = c(2.3, 3, 1, 0.8), mgp = c(2, 0.6, 0))
  ts_frame(d, range(0, m_cc[, "net_n"], m_x[, "net_n"]), zero = TRUE)
  lines(d, m_cc[, "net_n"], col = c0, lwd = 1.8); lines(d, m_x[, "net_n"], col = c1, lwd = 1.8)
  legend("topleft", c("Without controls", "With controls"), col = c(c0, c1), lwd = 1.8, bty = "n", cex = 0.8, horiz = TRUE)
}
f8 <- function() {
  par(mfrow = c(1, 2), mar = c(2.3, 3, 2, 0.8), mgp = c(2, 0.6, 0))
  yl <- range(0, sapply(c("short", "long"), function(b) range(FQ$cc[[b]][, "net_n"], FX[[b]][, "net_n"])))
  for (b in c("short", "long")) {
    ts_frame(d, yl, c(short = "(a) Short run (2-6 months)", long = "(b) Long run (more than 6 months)")[b], zero = TRUE)
    lines(d, FQ$cc[[b]][, "net_n"], col = c0, lwd = 1.8); lines(d, FX[[b]][, "net_n"], col = c1, lwd = 1.8)
    if (b == "short") legend("topleft", c("Without controls", "With controls"), col = c(c0, c1), lwd = 1.8, bty = "n", cex = 0.8)
  }
}
save_fig(file.path(od, "fig7_controls_net"), f7, 8, 3.4)
save_fig(file.path(od, "fig8_controls_freq_net"), f8, 10, 3.6)
