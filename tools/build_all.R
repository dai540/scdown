.libPaths(c(Sys.getenv("R_LIBS_USER"), .libPaths()))

options(repos = c(CRAN = "https://cloud.r-project.org"))

cran_pkgs <- c("pkgload", "testthat", "pkgdown", "rmarkdown", "knitr", "Matrix", "SeuratObject", "uwot", "roxygen2")
for (pkg in cran_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, lib = Sys.getenv("R_LIBS_USER"))
  }
}

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager", lib = Sys.getenv("R_LIBS_USER"))
}
bioc_pkgs <- c("SingleCellExperiment", "SummarizedExperiment", "S4Vectors")
missing_bioc <- bioc_pkgs[!vapply(bioc_pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_bioc)) {
  BiocManager::install(missing_bioc, ask = FALSE, update = FALSE, lib = Sys.getenv("R_LIBS_USER"))
}

pkg_dir <- getwd()
Sys.setenv(RSTUDIO_PANDOC = "C:/Program Files/RStudio/resources/app/bin/quarto/bin/tools")

roxygen2::roxygenise(pkg_dir)
pkgload::load_all(pkg_dir, quiet = TRUE)
testthat::test_dir(file.path(pkg_dir, "tests", "testthat"), reporter = "summary")
pkgload::unload("scdown")
pkgdown::build_site(pkg_dir, new_process = FALSE, install = TRUE)

cat("scdown rebuild complete\n")
