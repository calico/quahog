#' Create mzkit config file from defaults
#'
#' @param default_config_path path to a default config file
#' @inheritParams run_quahog
#' @param export_file_path path to new config file
#' @param rt_alignment_file_override path to an alignment file or an empty string to ignore
#' @param spline_ridge_penalty spline_ridge_penalty Ridge penalty on sample-level splines which model deviations between observed retention times and the consensus rt of a compound.
#' @param is_alignment should alignment be performed?
#' @param is_splitting should splitting be corrected?
#' @param is_coelution should coelutions be found?
#' @param is_search should search be performed?
#' @param is_qc should QC be performed?
#' @param is_eda should EDA be performed?
#' @param reference_sample_override NULL to skip or the name of a reference sample
#' @param is_exclude_isotopes should isotopologues be excluded?
#'
#' @return writes a json config to \code{export_file_path}
#'
#' @export
mzkit_config_file_from_default <- function(default_config_path,
                                           methodId,
                                           export_file_path,
                                           rt_alignment_file_override = "",
                                           spline_ridge_penalty = 20,
                                           is_alignment = TRUE,
                                           is_splitting = FALSE,
                                           is_coelution = TRUE,
                                           is_search = TRUE,
                                           is_qc = TRUE,
                                           is_eda = TRUE,
                                           reference_sample_override = NULL,
                                           is_exclude_isotopes = FALSE) {
  params <- jsonlite::fromJSON(default_config_path)

  params$pipeline$alignment$use <- is_alignment
  params$pipeline$splitting$use <- is_splitting
  params$pipeline$coelution$use <- is_coelution
  params$pipeline$search$use <- is_search
  params$pipeline$qc$use <- is_qc
  params$pipeline$eda$use <- is_eda

  # Use uploaded RT alignment file in peakdetector step
  if (rt_alignment_file_override != "") {
    params$modules$peakdetector$parameters$alignSamples <- 1
    params$modules$peakdetector$parameters$alignmentFile <- rt_alignment_file_override
  }

  if (is_exclude_isotopes) {
    params$modules$pipeline_standard_search$parameters$exclude_isotopes <- "T"
  } else {
    params$modules$pipeline_standard_search$parameters$exclude_isotopes <- "F"
  }

  # update ridge penalty
  stopifnot(class(spline_ridge_penalty) %in% c("numeric", "integer"), length(spline_ridge_penalty) == 1, spline_ridge_penalty >= 0)
  params$modules$pipeline_alignment$parameters$spline_ridge_penalty <- spline_ridge_penalty

  # udpate methods
  params$globals$methodId <- methodId

  if (methodId == "M001A") {
    params$globals$chemical_class <- "polar"
    params$globals$mode <- "negative"
    params$globals$chromatographic_method <- "'MS-Chrom-001-A C18-TBA'"
    params$globals$collision_energies <- "20,50,100"
    params$modules$pipeline_standard_search$parameters$matching_model <- "polar_forest.Rds"
    params$globals$dbname <- "mass_spec_standards"
  } else if (methodId == "M002A") {
    params$globals$chemical_class <- "polar"
    params$globals$mode <- "positive"
    params$globals$chromatographic_method <- "'MS-Chrom-004-A ZIC-pHILIC'"
    params$globals$collision_energies <- "20,40,80"
    params$modules$pipeline_standard_search$parameters$matching_model <- "polar_forest.Rds"
    params$globals$dbname <- "mass_spec_standards"
  } else if (methodId == "M003A") {
    params$globals$chemical_class <- "polar"
    params$globals$mode <- "negative"
    params$globals$chromatographic_method <- "'MS-Chrom-004-A ZIC-pHILIC'"
    params$modules$pipeline_standard_search$parameters$matching_model <- "polar_forest.Rds"
    params$globals$collision_energies <- "20,50,100"
    params$globals$dbname <- "mass_spec_standards"
  } else if (methodId == "M004A") {
    params$globals$chemical_class <- "lipid"
    params$globals$mode <- "negative"
    params$globals$chromatographic_method <- "'MS-Chrom-002-A C30-25cm 20mM NH4HCO2'"
    params$globals$collision_energies <- "20,30,40"
    params$modules$pipeline_standard_search$parameters$matching_model <- "lipids_fast_match"
    params$globals$dbname <- "mass_spec_standards_extended"
  } else if (methodId == "M005A") {
    params$globals$chemical_class <- "lipid"
    params$globals$mode <- "positive"
    params$globals$chromatographic_method <- "'MS-Chrom-002-A C30-25cm 20mM NH4HCO2'"
    params$globals$collision_energies <- "20,30,40"
    params$modules$pipeline_standard_search$parameters$matching_model <- "lipids_fast_match"
    params$globals$dbname <- "mass_spec_standards_extended"
  } else if (methodId == "M006A") {
    params$globals$chemical_class <- "fatty acid"
    params$globals$mode <- "negative"
    params$globals$chromatographic_method <- "'MS-Chrom-003-A C8-TBA'"
    params$globals$collision_energies <- "50,60,70"
    params$modules$pipeline_standard_search$parameters$matching_model <- "baseline"
    params$globals$dbname <- "mass_spec_standards"
  }

  if (!is.null(reference_sample_override)) {
    params$modules$pipeline_qc$parameters$referenceSample <- reference_sample_override
  }

  params_json_formatted <- capture.output(cat(jsonlite::toJSON(params, pretty = TRUE)))
  params_json_formatted_cleaned <- character(length(params_json_formatted))


  # clean json file returned by jsonlite::toJSON to match mzkit_v2.py specifications

  next_line <- ""
  for (x in 1:length(params_json_formatted)) {
    line <- params_json_formatted[x]

    if (x < length(params_json_formatted)) {
      next_line <- params_json_formatted[x + 1]
    }

    if (grepl("\\[*\\]", line) && !grepl("\"modules\":", line)) {
      start_coord <- regexpr("\\[", line)
      stop_coord <- regexpr("\\]", line)

      bracketed_substring <- substring(line, start_coord + 1, stop_coord - 1)

      if (grepl(",", bracketed_substring) && !grepl("collision_energies", line)) { # collision_energies is a string containing commas
        params_json_formatted_cleaned[x] <- line # pass line along without modification
      } else {
        pre_bracketed <- substring(line, 1, start_coord - 1)
        post_bracketed <- substring(line, stop_coord + 1, length(line))

        if (grepl("\\}", next_line) && !grepl("\\{", next_line)) {
          cleaned_line <- paste0(pre_bracketed, bracketed_substring, post_bracketed, sep = "")
        } else {
          cleaned_line <- paste0(pre_bracketed, bracketed_substring, post_bracketed, ",", sep = "")
        }

        params_json_formatted_cleaned[x] <- cleaned_line
      }
    } else {
      params_json_formatted_cleaned[x] <- line # pass line along without modification
    }
  }

  # export cleaned json structure to export file path
  fileConn <- file(export_file_path)
  writeLines(params_json_formatted_cleaned, fileConn)
  close(fileConn)

  invisible(0)
}


#' Sideload json config
#'
#' @description
#' Load an quahog configuration file
#'
#' @param rwrapper string corresponding to R module to run
#' @inheritParams run_quahog
#' @param config_path path to a config file (or NULL if config is output_folder/config.json)
#' @param r_scripts_path path to mass_spec/lib/R
#'
#' @return a named character vector of rwrapper arguments
#'
#' @export
sideload_json_config <- function(rwrapper, output_folder, config_path = NULL, r_scripts_path) {
  stopifnot(file.exists(output_folder))
  if (class(config_path) == "NULL") {
    json_config <- jsonlite::fromJSON(file.path(output_folder, "config.json"))
  } else {
    stopifnot("character" %in% class(config_path), length(config_path) == 1, file.exists(config_path))
    json_config <- jsonlite::fromJSON(config_path)
  }

  formatted_flags <- list()

  formatted_flags[["rwrapper"]] <- rwrapper
  formatted_flags[["output_folder"]] <- output_folder
  formatted_flags[["r_scripts_path"]] <- r_scripts_path

  # add globals

  formatted_flags <- append(formatted_flags, json_config$globals)

  # add module-specific parameters

  formatted_flags <- append(formatted_flags, json_config$modules[[rwrapper]]$parameters)

  # Need to remove single quotes from string to avoid error in clamdb::match_method(dots, mass_spec_standards_con)
  formatted_flags$chromatographic_method <- gsub("'", "", formatted_flags$chromatographic_method)

  formatted_flags
}
