#!/usr/bin/env Rscript
# install_R_packages.R
# Installs required R packages for the GenePrioritiser project.
# Usage: Rscript install_R_packages.R

required_cran <- c(
  'tidyverse',
  'magrittr',
  'data.table',
  'enrichR',
  'stringr',
  'ggplot2',
  'reshape2'
)

required_bioc <- c(
  'org.Hs.eg.db'
)

install_if_missing <- function(pkgs, repos = getOption('repos')) {
  to_install <- pkgs[!pkgs %in% rownames(installed.packages())]
  if (length(to_install) == 0) {
    message('All requested CRAN packages already installed.')
    return(invisible(TRUE))
  }
  message('Installing CRAN packages: ', paste(to_install, collapse = ', '))
  install.packages(to_install, repos = repos, dependencies = TRUE)
}

install_bioc_if_missing <- function(pkgs) {
  if (!requireNamespace('BiocManager', quietly = TRUE)) {
    install.packages('BiocManager', dependencies = TRUE)
  }
  to_install <- pkgs[!pkgs %in% rownames(installed.packages())]
  if (length(to_install) == 0) {
    message('All requested Bioconductor packages already installed.')
    return(invisible(TRUE))
  }
  message('Installing Bioconductor packages: ', paste(to_install, collapse = ', '))
  BiocManager::install(to_install, ask = FALSE)
}

main <- function() {
  # ensure a CRAN repo is set (default to RStudio CRAN mirror if unset)
  if (is.null(getOption('repos')) || getOption('repos')['CRAN'] == '@CRAN@') {
    options(repos = c(CRAN = 'https://cloud.r-project.org'))
  }

  install_if_missing(required_cran)
  install_bioc_if_missing(required_bioc)

  message('\nInstallation complete. You can now run the R scripts in src/data_preprocessing/.')
}

if (!interactive()) {
  main()
}
