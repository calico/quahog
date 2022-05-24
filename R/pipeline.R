#' Run Quahog
#'
#' Run the quahog pipeline directly from R with a configurable config file.
#'
#' This is only meant as a convenience function for running mzkit with
#' an appropriate config during development. The primary methods
#' for running quahog is via Docker.
#'
#' @param data_folder directory containing sample mzXML/mzML files.
#' @param output_folder directory where output files should be saved.
#' @param methodId mass spectrometry method being run (e.g., M001A)
#' @param config_args list of config arguments to change from the method-specific defaults
#' @param path_to_mass_spec path to mass_spec github repo
#' @param conda_env conda environment to use when running quahog
#' @param conda_path path to conda directory
#' @param background run process in the background?
#'
#' @return a proc object for the pipeline running in the background
#'
#' @export
run_quahog <- function(data_folder, output_folder = "/tmp/metabolomics_test", methodId,
                       config_args = list(), path_to_mass_spec = "~/mass_spec",
                       conda_env, conda_path, background = FALSE) {
  stopifnot("character" %in% class(data_folder), length(data_folder) == 1)
  stopifnot("character" %in% class(output_folder), length(output_folder) == 1)
  stopifnot("character" %in% class(methodId), length(methodId) == 1)
  stopifnot("list" %in% class(config_args))

  stopifnot("character" %in% class(path_to_mass_spec), length(path_to_mass_spec) == 1)
  if (!all(c("analysis", "lib", "docker", "info") %in% list.files(path_to_mass_spec))) {
    stop(path_to_mass_spec, " is not a path to the mass_spec github repo")
  }

  stopifnot("character" %in% class(conda_env), length(conda_env) == 1)
  stopifnot("character" %in% class(conda_path), length(conda_path) == 1)

  # setup paths

  mzkit_py_path <- normalizePath(file.path(path_to_mass_spec, "lib", "python", "mzkit_v2.py"))
  default_config <- file.path(path_to_mass_spec, "lib", "assets", "pipeline_configs", "default_config.json")

  # setup python environment to be used in mzkit_v2.py call

  # create config
  tmp_config_path <- file.path("/tmp", "mzkit_config.json")

  # build a relevant config

  config_args_available <- names(formals(mzkit_config_file_from_default)) %>%
    setdiff(c("default_config_path", "export_file_path", "methodId"))

  config_args_provided <- names(config_args)
  config_args_provided_valid <- intersect(config_args_provided, config_args_available)
  config_args_provided_invalid <- setdiff(config_args_provided, config_args_available)

  if (length(config_args_provided_invalid) != 0) {
    warning(length(config_args_provided_invalid), " args provided in \"config_args\" can't be used: ", paste(config_args_provided_invalid, collapse = ", "))
  }

  config_args <- append(
    config_args,
    list(
      default_config_path = default_config,
      methodId = methodId,
      export_file_path = tmp_config_path
    )
  )

  # generate config at tmp_config_path
  do.call(mzkit_config_file_from_default, config_args)

  # run the pipeline in foreground or background

  sh_script_path <- normalizePath("~/mass_spec/lib/sh/job_submission/call_mzkit_from_R.sh")

  conda_root_path <- dirname(dirname(conda_path))
  sh_args <- c(conda_root_path, conda_env, mzkit_py_path, data_folder, output_folder, tmp_config_path)

  if (background) {
    proc <- process$new(sh_script_path, sh_args, stdout = "|", "stderr" = "|")
    return(proc)
  } else {
    processx::run(sh_script_path, args = sh_args, spinner = TRUE)
    return(invsible(0))
  }
}
