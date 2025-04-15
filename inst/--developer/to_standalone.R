
# path to standalone
outpath <- "../household-epidemics/inst/standalone/"

developer_export_to_standalone <- function(hhmodel, outpath) {
  write.csv(hhmodel$dataset,
            file.path(outpath, "data.tmp"),
            quote = FALSE, row.names = FALSE)
  write.csv(hhmodel$augm_defs,
            file.path(outpath, "augm.tmp"),
            quote = FALSE, row.names = FALSE)
  write.csv(hhmodel$param_defs,
            file.path(outpath, "param.tmp"),
            quote = FALSE, row.names = FALSE)
  write.table(t(as.data.frame(hhmodel$setup)),
              file.path(outpath, "setup.tmp"),
              row.names = TRUE, col.names = FALSE, sep = ",", quote = FALSE)
  write.table(t(as.data.frame(hhmodel$infec_profile_defs)),
              file.path(outpath, "infec.tmp"),
              row.names = TRUE, col.names = FALSE, sep = ",", quote = FALSE)
}

#developer_export_to_standalone(hhmodel, outpath)

