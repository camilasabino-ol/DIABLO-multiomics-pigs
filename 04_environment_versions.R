# Read-only environment report. Captures CURRENT server environment only.
# Does not prove historical versions used for archived analysis.
pkgs <- c("mixOmics","BiocParallel","igraph","cowplot","ggplot2","viridis",
          "dplyr","readxl","readr","writexl")
cat("Session information on date of this check:\n")
print(sessionInfo())
cat("\nPackage versions available in this environment:\n")
for (pkg in pkgs) {
  installed <- requireNamespace(pkg, quietly=TRUE)
  cat(sprintf("%-15s %s\n", pkg,
              if (installed) as.character(packageVersion(pkg))
              else "NOT INSTALLED"))
}
if (file.exists("final_diablo_model.rds")) {
  cat("\nArchived RDS serialization metadata:\n")
  print(infoRDS("final_diablo_model.rds"))
}
