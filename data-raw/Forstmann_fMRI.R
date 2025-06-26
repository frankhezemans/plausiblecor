library(dplyr, warn.conflicts = FALSE)
library(purrr)
library(tibble)
library(readr)
library(usethis)

# Download data from Forstmann et al. (2008), released with Ly et al. (2017) on
# OSF:
# Forstmann, B.U., Dutilh, G., Brown, S., Neumann, J., Von Cramon, D.Y.,
#   Ridderinkhof, K.R., & Wagenmakers, E.J. (2008). Striatum and pre-SMA
#   facilitate decision-making under time pressure. Proceedings of the National
#   Academy of Sciences, 105(45), 17538-17542.
#   https://doi.org/10.1073/pnas.0805903105
#
# Ly, A., Boehm, U., Heathcote, A., Turner, B.M., Forstmann, B., Marsman, M. &
#   Matzke, D. (2017). A Flexible and Efficient Hierarchical Bayesian Approach
#   to the Exploration of Individual Differences in Cognitive-model-based
#   Neuroscience. In Computational Models of Brain and Behavior, A.A. Moustafa
#   (Ed.). https://doi.org/10.1002/9781119159193.ch34

load_forstmann_osf <- function(osf_link) {
  temp_file <- tempfile(fileext = ".RData")
  utils::download.file(
    url = osf_link,
    destfile = temp_file,
    mode = "wb"
  )
  temp_env <- new.env(parent = emptyenv())
  load(file = temp_file, envir = temp_env)
  result <- as.list(temp_env)
  return(result)
}

Forstmann_fMRI <- load_forstmann_osf("https://osf.io/download/f7553/") %>%
  # select only the covariates data (other list element is behavioural data)
  purrr::pluck("covariates") %>%
  tibble::rownames_to_column(var = "subjects") %>%
  dplyr::mutate(subjects = factor(subjects)) %>%
  dplyr::rename(
    # use "striatum" as more accurate column name, that corresponds to
    # terminology used by Forstmann et al. (2008)
    striatum = basalganglia,
    pre_sma = presma
  )

readr::write_csv(Forstmann_fMRI, "data-raw/Forstmann_fMRI.csv")
usethis::use_data(Forstmann_fMRI, overwrite = TRUE)
