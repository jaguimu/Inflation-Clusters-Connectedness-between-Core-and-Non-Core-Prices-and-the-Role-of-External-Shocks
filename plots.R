# Plot helpers (base R + igraph)
COL <- c("#1b3a5c", "#c8553d", "#2a9d8f", "#d4a017", "#6d4c8d")
COL_NC <- "#e07b00"; COL_CN <- "#2a9d8f"; COL_X <- "#6a4c93"
PAL <- colorRampPalette(c("#2c5d8a", "#9fbbd6", "#f4f2ee", "#e8a58f", "#b8432f"))(101)

LAB <- c(AYB_SUB = "Food &\nbeverages", TEX_CALZ_SUB = "Textiles &\nfootwear", ELEC_SUB = "Household\nappliances",
         INDUSTR_SUB = "Other ind.\ngoods", CFH_SUB = "Food away\nfrom home", EDUC_SUB = "Education",
         SALUD_SUB = "Health", ALQ_SUB = "Rents", SERV_OTROS_SUB = "Other\nservices",
         ALIM_NSUB = "Food", COMB_NSUB = "Fuels", TRAN_NSUB = "Transport", PUB_NSUB = "Public\nutilities")

save_fig <- function(file, fun, w, h, res = 200) {
  pdf(paste0(file, ".pdf"), width = w, height = h); fun(); dev.off()
  png(paste0(file, ".png"), width = w, height = h, units = "in", res = res); fun(); dev.off()
}

ts_frame <- function(d, ylim, main = "", zero = FALSE, by = "2 years", xlim = range(d), cex.main = 0.95) {
  plot(d, rep(NA_real_, length(d)), type = "n", ylim = ylim, xlim = xlim, xlab = "", ylab = "", xaxt = "n",
       las = 1, bty = "l", main = main, cex.main = cex.main, font.main = 1, cex.axis = 0.8)
  yrs <- seq(as.Date("2010-01-01"), as.Date("2026-01-01"), by = by)
  axis(1, at = yrs, labels = format(yrs, "%Y"), cex.axis = 0.8)
  abline(h = pretty(ylim), col = "grey92", lwd = 0.8)
  if (zero) abline(h = 0, col = "grey50")
}

offd <- function(m) { diag(m) <- 0; m }
acol <- function(col, a) { v <- col2rgb(col) / 255; rgb(v[1], v[2], v[3], alpha = a) }

# node positions from average VD
net_layout <- function(m) {
  S <- offd(m); S <- (S + t(S)) / 2
  g <- igraph::graph_from_adjacency_matrix(S, mode = "undirected", weighted = TRUE, diag = FALSE)
  set.seed(7); igraph::layout_with_fr(g, weights = igraph::E(g)$weight, niter = 5000)
}
net_scale <- function(TH) list(tor = range(sapply(TH, function(m) colSums(offd(m)))),
                               wmax = max(sapply(TH, function(m) max(offd(m)))))

# m in %, m[i, j] = j -> i
net_panel <- function(m, isN, L, sc, main, lcex = 0.62, asz = 0.35, tline = 0.5, tcex = 1) {
  dd <- offd(m); vn <- rownames(dd)
  g <- igraph::graph_from_adjacency_matrix(t(dd), mode = "directed", weighted = TRUE, diag = FALSE)
  el <- igraph::as_edgelist(g, names = FALSE); s <- igraph::E(g)$weight / sc$wmax
  fn <- isN[el[, 1]]; tn <- isN[el[, 2]]
  type <- ifelse(fn & !tn, "nc", ifelse(!fn & tn, "cn", "within"))
  a <- 0.05 + 0.95 * pmin(1, 1.8 * s)
  ecol <- ifelse(type == "nc", acol(COL_NC, a), ifelse(type == "cn", acol(COL_CN, a),
                 gray(0.8 - 0.65 * s, alpha = 0.08 + 0.7 * s)))
  ewid <- ifelse(type == "within", 0.2 + 4.7 * s, 0.3 + 6.5 * s)
  nv <- colSums(dd) - rowSums(dd); to <- colSums(dd)
  ci <- round(50 * pmax(-1, pmin(1, nv / 25))) + 51
  vs <- 12 + 18 * (to - sc$tor[1]) / diff(sc$tor)
  o <- order(type != "within")
  g2 <- igraph::make_graph(as.vector(t(el[o, , drop = FALSE])), n = length(vn), directed = TRUE)
  plot(g2, layout = L, vertex.color = PAL[ci], vertex.frame.color = "grey30",
       vertex.shape = ifelse(isN, "square", "circle"), vertex.size = vs,
       vertex.label = paste0(LAB[vn], "\n", sprintf("%+.1f", nv)), vertex.label.cex = lcex,
       vertex.label.color = "black", vertex.label.degree = pi / 2, vertex.label.dist = 0.9 + vs / 9,
       vertex.label.family = "sans", edge.width = ewid[o], edge.color = ecol[o], edge.arrow.size = asz,
       edge.curved = 0.15, margin = c(0.1, 0, 0, 0))
  title(main, cex.main = tcex, font.main = 1, line = tline)
}
net_legend <- function() {
  par(mar = c(0, 0, 0, 0)); plot.new()
  legend("center", horiz = TRUE, bty = "n", cex = 0.72, pt.cex = 1.5, seg.len = 1.6,
         legend = c("Non-core", "Core", "Net transmitter", "Net receiver", "Non-core to core", "Core to non-core", "Within cluster"),
         pch = c(22, 21, 21, 21, NA, NA, NA), lty = c(NA, NA, NA, NA, 1, 1, 1), lwd = c(NA, NA, NA, NA, 3, 3, 3),
         pt.bg = c("white", "white", "#b8432f", "#2c5d8a", NA, NA, NA), col = c(rep("grey30", 4), COL_NC, COL_CN, "grey55"))
}
