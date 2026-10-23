# Runs the full replication
for (f in c("01_data.R", "02_estimation.R", "03_tables.R", "04_figures.R", "05_appendix.R")) {
  message(f); source(f, echo = FALSE, local = new.env())
}
