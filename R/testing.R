#' Test Rwrapper
#'
#' Test whether the R module being called matches one of the implemented modules.
#'
#' @export
test_rwrapper <- function() {
  # Test for valid options of rwrapper
  valid_wrappers <- c(
    "pipeline_alignment",
    "pipeline_coelution_detection",
    "pipeline_aggregate_split_peaks",
    "pipeline_coelution_labeling",
    "pipeline_standard_search",
    "pipeline_qc",
    "pipeline_eda"
  )

  # test that rwrapper is in an object in global environment
  if (!("rwrapper" %in% objects(envir = globalenv()))) {
    stop("\"rwrapper\" must exist in the global environment")
  }

  # verify that rwrapper's value is one of the reserved values in valid_wrappers
  stopifnot(length(rwrapper) == 1, class(rwrapper) == "character")
  if (!(rwrapper %in% valid_wrappers)) {
    stop("Unrecognized value for rwrapper; possible arguments are: ", paste(valid_wrappers, collapse = ", "))
  }
}

test_output_folder <- function(dots) {
  # test that output_folder is present in dots (list)
  if (!("output_folder" %in% names(dots))) {
    stop("output_folder is a required argument")
  }

  # test that output_folder exists on the file system
  if (!file.exists(dots$output_folder)) {
    stop("output_folder not found at: ", dots$output_folder)
  }
}

test_mzroll_db_file <- function(dots) {
  # test that mzroll_db_file is present in dots (list)
  if (!("mzroll_db_file" %in% names(dots))) {
    stop("mzroll_db_file is a required argument")
  }

  # test output_folder since this mzroll_db_path is generated from output_folder / mzroll_db_file
  test_output_folder(dots)

  # test that mzroll_db_file is present
  mzroll_db_path <- file.path(dots$output_folder, dots$mzroll_db_file)

  if (!file.exists(mzroll_db_path)) {
    stop("mzroll file not found at: ", mzroll_db_path, "\nThis mzroll_db_path is generated from:\n    output_folder: ", dots$output_folder, "\n    mzroll_db_file: ", dots$mzroll_db_file)
  }
}

#' Test Quahog Packages
#'
#' Determine whether the R packages required to run the quahog pipeline are installed.
#'
#' @export
test_quahog_packages <- function() {

  # check packages which are required to run the core quahog pipeline

  required_packages <- c(
    "clamshell",
    "calicomics",
    "mzkitcpp",
    "clamr",
    "clamdb",
    "clamqc",
    # cran
    "fastcluster",
    "furrr",
    "future",
    "fuzzyjoin",
    "lubridate",
    "Rcpp",
    "RMySQL",
    # bioconductor
    "IRanges",
    "mzR",
    "debugr"
  )

  required_packages_installed <- as.logical(lapply(required_packages, function(x) {
    requireNamespace(x, quietly = TRUE)
  }))

  if (any(!required_packages_installed)) {
    stop("required R packages are not installed: ", paste(required_packages[!required_packages_installed], collapse = ", "), "\n",
      "please install ", paste(required_packages[!required_packages_installed], collapse = " & "), " and continue",
      call. = FALSE
    )
  }
}
