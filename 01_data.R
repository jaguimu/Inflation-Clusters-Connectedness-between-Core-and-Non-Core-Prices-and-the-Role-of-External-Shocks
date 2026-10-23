# Data: CPI components (SA) and controls
library(readxl)
source("R/functions.R")
dir.create("output", showWarnings = FALSE)

rd <- function(s) {
  d <- as.data.frame(read_excel("data/data.xlsx", sheet = s, .name_repair = "minimal"))
  names(d)[1] <- "date"; d <- d[!is.na(d$date), ]; d$date <- as.Date(d$date)
  d[order(d$date), ]
}
core <- rd("cluster1"); ncore <- rd("cluster2"); ctl <- rd("controls")

pct <- function(x) c(NA, 100 * (x[-1] / x[-length(x)] - 1))
for (v in c("energy", "food", "fert", "tc", "IPMI")) ctl[[v]] <- pct(ctl[[v]])

dat <- merge(merge(core, ncore, by = "date"), ctl, by = "date")
dat <- dat[dat$date >= as.Date("2000-01-01") & dat$date <= as.Date("2025-12-01"), ]
dates <- dat$date

vi <- c(names(core)[-1], names(ncore)[-1])
vx <- c("energy", "food", "fert", "tc", "gscpi")
stopifnot(!anyNA(dat[, c(vi, vx, "IPMI")]))

Y_raw <- as.matrix(dat[, vi])
sa <- seasonal_adjust(Y_raw, dates)
Y <- sa$Y
id <- setNames(c(rep("core", ncol(core) - 1), rep("noncore", ncol(ncore) - 1)), vi)

X <- as.matrix(dat[, vx])
ipmi <- dat$IPMI

write.csv(sa$diag, "output/seasonal_diagnostics.csv", row.names = FALSE)
saveRDS(list(Y = Y, Y_raw = Y_raw, X = X, ipmi = ipmi, id = id, dates = dates), "output/data.rds")
