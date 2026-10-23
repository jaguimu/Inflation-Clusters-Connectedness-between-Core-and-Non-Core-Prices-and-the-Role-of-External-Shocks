# Tables 1-4
source("R/functions.R")
r <- readRDS("output/rolling.rds"); d <- r$dates
dir.create("output/tables", recursive = TRUE, showWarnings = FALSE)
g <- list(all = rep(TRUE, length(d)), y0921 = d < as.Date("2022-01-01"), y2225 = d >= as.Date("2022-01-01"))
M <- function(a, id) t(apply(a, 3, conn_measures, id = id))
avg <- function(m) sapply(g, function(i) colMeans(m[i, , drop = FALSE]))

k1 <- c("c_own", "c_comove", "c_contag", "n_own", "n_comove", "n_contag", "total", "within", "between", "net_n")
k3 <- c("c_contag", "n_contag", "between", "net_n")
pair <- function(a, b) { A <- avg(a)[k3, ]; B <- avg(b)[k3, ]
  z <- do.call(cbind, lapply(names(g), function(p) cbind(A[, p], B[, p])))
  colnames(z) <- paste(rep(names(g), each = 2), c("without", "with"), sep = "_"); z }

m_cc <- M(r$cc, r$id); m_x <- M(r$x, r$idX)
fq_cc <- lapply(r$fq_cc, M, id = r$id); fq_x <- lapply(r$fq_x, M, id = r$idX)

t1 <- avg(m_cc)[k1, ]
t2 <- cbind(short = avg(fq_cc$short)[k1, ], long = avg(fq_cc$long)[k1, ])
t3 <- pair(m_cc, m_x)
t4 <- rbind(short = pair(fq_cc$short, fq_x$short), long = pair(fq_cc$long, fq_x$long))
rownames(t4) <- paste(rep(c("short", "long"), each = 4), k3, sep = "_")

for (k in c("t1", "t2", "t3", "t4")) {
  x <- round(get(k), 2); write.csv(x, sprintf("output/tables/table%s.csv", substr(k, 2, 2)))
  cat("\nTable", substr(k, 2, 2), "\n"); print(x)
}
