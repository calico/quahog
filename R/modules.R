#' Pipeline Alignment
#'
#' Align retention time across samples using M2-based alignment
#'
#' @param ... all arguments are passed in from mzkit.R using \code{do.call(rwrapper, as.list(formatted_flags))}
#' @return 0 invisibly, function called for side effects.
#'
#' @export
pipeline_alignment <- function(...) {
  debugr::dwatch(msg = "Starting pipeline_alignment. [mzkit<pipeline_wrappers.R>::pipeline_alignment]")

  # load an mzroll db and align samples (across retention time & full scan m/z)

  dots <- list(...)

  # output_folder, mzroll_db_file, MS1tol and MS2tol are required arguments which will be passed as ...
  # along with any optional arguments ms2_driven_aligner

  test_output_folder(dots)
  output_folder <- normalizePath(dots$output_folder)

  test_mzroll_db_file(dots)
  mzroll_db_path <- normalizePath(file.path(dots$output_folder, dots$mzroll_db_file))

  debugr::dwatch(msg = paste0("mzroll_db_path= ", mzroll_db_path, " [mzkit<pipeline_wrappers.R>::pipeline_alignment]", sep = " "))

  # test mzrollDB's schema and open connection

  mzroll_db_con <- clamr::mzroll_db_sqlite(mzroll_db_path)
  clamr_config <- clamr::build_clamr_config(dots, quietly = FALSE)

  debugr::dwatch(msg = paste0("mzroll_db_connection= ", slot(mzroll_db_con, "dbname"), " [mzkit<pipeline_wrappers.R>::pipeline_alignment]"))

  # Combine specified parameters with defaults to create function's call arguments

  if (!("maximum_mz_set_size" %in% names(dots))) {
    n_samples <- dplyr::tbl(mzroll_db_con, dbplyr::sql("SELECT COUNT (DISTINCT sampleId) FROM samples")) %>%
      dplyr::collect() %>%
      unname() %>%
      unlist()
    dots$maximum_mz_set_size <- 4 * n_samples
  }

  aligner_args <- reclass_cmd_args(
    chr_list = dots, call_fxn = clamr::ms2_driven_aligner,
    default_args = list(
      peak_quality_cutoff = 0.5,
      cosine_cutoff = 0.95,
      sd_rt_resid_cutoff = 0.1,
      spline_ridge_penalty = 200,
      return_plot_data = TRUE
    )
  )
  if (class(aligner_args$spline_ridge_penalty) == "integer") {
    aligner_args$spline_ridge_penalty <- as.numeric(aligner_args$spline_ridge_penalty)
  }

  aligner_args$mzroll_db_con <- mzroll_db_con
  aligner_args$clamr_config <- clamr_config

  # Perform the sample-to-sample alignment guided by peaks that match in precursor M/Z and MS2 fragmentation profiles

  debugr::dwatch(msg = "Prepared data for ms2_driven_aligner step. [mzkit<pipeline_wrappers.R>::pipeline_alignment]")

  sample_aligner <- do.call(clamr::ms2_driven_aligner, aligner_args)

  debugr::dwatch(msg = "finished ms2_driven_aligner. [mzkit<pipeline_wrappers.R>::pipeline_alignment]")

  # save retention time and ms1 drift deviations

  overwritten_tables <- intersect(DBI::dbListTables(mzroll_db_con), c("rt_update_key", "mz_update_key"))
  if (length(overwritten_tables) != 0) {
    warning(length(overwritten_tables), " tables (", paste(overwritten_tables, collapse = " & "), ") are already present in the mzroll_db_file and will be overwritten")
  }

  debugr::dwatch(msg = "created overwritten_tables. [mzkit<pipeline_wrappers.R>::pipeline_alignment]")

  clamshell::sqlite_write(mzroll_db_con, "rt_update_key", sample_aligner$rt_update_key, overwrite = TRUE)
  clamshell::sqlite_write(mzroll_db_con, "mz_update_key", sample_aligner$MS1_drift, overwrite = TRUE)
  DBI::dbDisconnect(mzroll_db_con)

  debugr::dwatch(msg = "disconnected from db. [mzkit<pipeline_wrappers.R>::pipeline_alignment]")

  # save plots (if present)

  if (aligner_args$return_plot_data) {

    # Create the QC directory if it does not already exist.
    qc_subfolder <- file.path(output_folder, "QC")
    dir.create(qc_subfolder)

    # Set permissions to facilitate easy deletion by mzkitchen
    rt_plot_data_file <- file.path(qc_subfolder, "rt_plot_data.Rds")

    saveRDS(sample_aligner$plot_data, rt_plot_data_file)

    system(glue::glue("chmod 777 {rt_plot_data_output_file}", rt_plot_data_output_file = rt_plot_data_file))
  }

  cat("MS2-based alignment DONE - rt_update_key & mz_update_key added to mzrollDB\n")
  invisible(0)
}

#' Pipeline Coelution Detection
#'
#' Summarizes an mzDeltas output to identify coelution relationships between peaks.
#' Output tables are added to the mzrollDB database.
#'
#' @param ... all arguments are passed in from mzkit.R using \code{do.call(rwrapper, as.list(formatted_flags))}
#' @return 0 invisibly, function called for side effects.
#'
#' @export
pipeline_coelution_detection <- function(...) {
  dots <- list(...)

  # output_folder, mzroll_db_file, adduct_file and MS1tol are required arguments which will be passed as ...

  test_output_folder(dots)
  output_folder <- normalizePath(dots$output_folder)

  test_mzroll_db_file(dots)
  mzroll_db_path <- normalizePath(file.path(dots$output_folder, dots$mzroll_db_file))
  mzroll_db_con <- clamr::mzroll_db_sqlite(mzroll_db_path)

  if (!("adduct_file" %in% names(dots))) {
    stop("adduct_file is a required flag when calling coelutionPeaks. e.g., adduct_file=mzdeltas.out")
  }
  adduct_path <- file.path(output_folder, dots$adduct_file)
  if (!file.exists(adduct_path)) {
    stop("adduct_file not found at ", adduct_path)
  }

  # test and format mass accuracy tolerances

  missing_instrument_tolerances <- setdiff("MS1tol", names(dots))
  if (length(missing_instrument_tolerances) != 0) {
    stop("Instrument tolerances must be specified with ", paste(missing_instrument_tolerances, collapse = ", "))
  }

  # test mode & analyte_type validity (for subsetting the adduct table).
  if (!"mode" %in% names(dots) || !dots$mode %in% c("negative", "positive")) {
    stop("mode is a required flag when calling coelutionPeaks; valid options are negative and positive")
  }
  if (!"chemical_class" %in% names(dots) || !dots$chemical_class %in% c("polar", "lipid", "fatty acid")) {
    stop("chemical_class is a required flag when calling coelutionPeaks; valid options are polar, lipid, fatty acid")
  }

  clamr_config <- clamr::build_clamr_config(dots)

  # Combine specified parameters with defaults to create function arguments

  find_common_coelutions_args <- reclass_cmd_args(dots, clamr::find_common_coelutions,
    default_args = list(coelution_n_cutoff = 1000L)
  )
  find_common_coelutions_args$adduct_path <- adduct_path

  # Find m/z differences between peaks which are common and the associated [mz, scans] of these peaks

  common_coelution_events <- do.call(clamr::find_common_coelutions, find_common_coelutions_args)

  # define a set of coelutions (isotopologues, adducts, etc.) that are expected in this mode / chemical class
  coelution_library <- clamr::standard_coelution_library %>%
    dplyr::filter(is.na(mode_specific) | mode_specific == dots$mode) %>%
    dplyr::filter(is.na(analyte_specific) | analyte_specific == dots$chemical_class)

  # match coelutions to a library of high confidence coelutions (e.g., +13C, +Na, ...)

  coelution_labels <- clamr::identify_coelutions(common_coelution_events,
    type = "library",
    coelution_n_cutoff = find_common_coelutions_args$coelution_n_cutoff,
    quietly = TRUE,
    coelution_library = coelution_library,
    amu_tolerance = 0.003
  ) %>%
    dplyr::group_by(formula) %>%
    dplyr::arrange(match_distance) %>%
    # assign all library formulas to the mz difference which best matches
    dplyr::slice(1) %>%
    dplyr::group_by(mz_set, mz_subset) %>%
    dplyr::summarize(
      formula = paste0(formula, collapse = " | "),
      type = type[1],
      is_loss = is_loss[1]
    ) %>%
    dplyr::ungroup() %>%
    dplyr::rename(library_formula = formula, library_is_loss = is_loss, library_type = type) %>%
    dplyr::filter(library_type %in% c("isotopologue", "adduct")) %>%
    dplyr::select(mz_set, mz_subset, library_formula, library_is_loss, library_type)

  reduced_coelution_events <- clamr::apply_coelution_labels(common_coelution_events, coelution_labels)

  # save a list of coeluting masses and types of masses as a binary
  overwritten_tables <- intersect(DBI::dbListTables(mzroll_db_con), c("coelutions", "coelution_mz_diffs", "coelution_peaks"))
  if (length(overwritten_tables) != 0) {
    warning(length(overwritten_tables), " tables (", paste(overwritten_tables, collapse = " & "), ") are already present in the mzroll_db_file and will be overwritten")
  }

  clamshell::sqlite_write(mzroll_db_con, "coelutions", reduced_coelution_events$coelutions, overwrite = TRUE)
  clamshell::sqlite_write(mzroll_db_con, "coelution_mz_diffs", reduced_coelution_events$coelution_mz_diffs, overwrite = TRUE)

  if (nrow(reduced_coelution_events$coelutions) != 0) {
    coelution_peaks <- reduced_coelution_events$coelutions %>%
      tidyr::nest(coelutions = -sample) %>%
      # for each sample, determine coelution peaks [mz, scans] and then match these to clamr db peaks
      dplyr::mutate(matches = purrr::map(coelutions,
        clamr::generate_coelution_peaks,
        MS1tol = clamr_config$MS1tol
      ) %>%
        purrr::transpose() %>%
        {
          .[["unique_coeluting_peaks"]]
        }) %>%
      dplyr::select(-coelutions) %>%
      tidyr::unnest(matches)
  } else {
    # return an empty tibble with correct classes
    coelution_peaks <- tibble::tibble(
      sample = "", mz_scan_group = "", minscan = 1L, maxscan = 1L, peakMz = 0.1,
      mzmin = 0.1, mzmax = 0.1, n = 1L
    ) %>%
      dplyr::slice(-1)
  }

  clamshell::sqlite_write(mzroll_db_con, "coelution_peaks", coelution_peaks, overwrite = TRUE)
  DBI::dbDisconnect(mzroll_db_con)

  cat("Coelution detection DONE - coelutions, coelution_mz_diffs & coelution_peaks added to mzrollDB\n")
  invisible(0)
}

#' Pipeline Aggregate Split Peaks
#'
#' Identifies and merges split peaks caused by the same analyte being picked up in
#' different peakgroups across samples.
#'
#' @param ... all arguments are passed in from mzkit.R using \code{do.call(rwrapper, as.list(formatted_flags))}
#' @return 0 invisibly, function called for side effects.
#'
#' @export
pipeline_aggregate_split_peaks <- function(...) {
  debugr::dwatch(msg = "started pipeline_aggregate_split_peaks [mass_spec<pipeline_wrappers.R>::pipeline_aggregate_split_peaks]")
  dots <- list(...)

  # output_folder, mzroll_db_file, MS1tol and MS2tol are required arguments which will be passed as ...
  # along with any optional arguments for peakgrp_coelution_annotation

  test_output_folder(dots)
  output_folder <- normalizePath(dots$output_folder)

  test_mzroll_db_file(dots)
  mzroll_db_path <- normalizePath(file.path(dots$output_folder, dots$mzroll_db_file))

  # test mzrollDB's schema and open connection

  mzroll_db_con <- clamr::mzroll_db_sqlite(mzroll_db_path)
  debugr::dwatch(msg = "connected to mzroll_db. [mass_spec<pipeline_wrappers.R>::pipeline_aggregate_split_peaks]")

  # remove redundant peakgroups
  clamr_config <- clamr::build_clamr_config(dots)

  identify_split_peaks_args <- reclass_cmd_args(dots, clamr::identify_split_peaks)
  identify_split_peaks_args$mzroll_db_con <- mzroll_db_con
  identify_split_peaks_args$clamr_config <- clamr_config

  # Perform the sample-to-sample alignment guided by peaks that match in precursor M/Z and MS2 fragmentation profiles

  debugr::dwatch(msg = "prepared data for split peak identification step. [mzkit<pipeline_wrappers.R>::pipeline_aggregate_split_peaks]")

  split_peaks <- do.call(clamr::identify_split_peaks, identify_split_peaks_args)

  debugr::dwatch(msg = "finished identify_split_peaks_args. [mzkit<pipeline_wrappers.R>::pipeline_aggregate_split_peaks]")

  clamr::merge_peakgroups(mzroll_db_con, split_peaks)

  cat("Split peak aggregation DONE - peakgroups and peaks updated in mzrollDB\n")
  invisible(0)
}

#' Pipeline Coelution Labeling
#'
#' Groups coeluting peakgroups which are likely to be adducts or isotopologues of one another.
#'
#' @param ... all arguments are passed in from mzkit.R using \code{do.call(rwrapper, as.list(formatted_flags))}
#' @return 0 invisibly, function called for side effects.
#'
#' @export
pipeline_coelution_labeling <- function(...) {
  debugr::dwatch(msg = "started pipeline_coelution_labeling. [mass_spec<pipeline_wrappers.R>::pipeline_coelution_labeling]")
  dots <- list(...)

  # output_folder, mzroll_db_file, MS1tol and MS2tol are required arguments which will be passed as ...
  # along with any optional arguments for peakgrp_coelution_annotation

  test_output_folder(dots)
  output_folder <- normalizePath(dots$output_folder)

  test_mzroll_db_file(dots)
  mzroll_db_path <- normalizePath(file.path(dots$output_folder, dots$mzroll_db_file))

  # test mzrollDB's schema and open connection

  missing_instrument_tolerances <- setdiff(c("MS1tol"), names(dots))
  if (length(missing_instrument_tolerances) != 0) {
    stop("Instrument tolerances must be specified with ", paste(missing_instrument_tolerances, collapse = ", "))
  }

  if (!"mode" %in% names(dots) || !dots$mode %in% c("negative", "positive")) {
    stop("mode is a required flag when calling coelutionPeaks; valid options are negative and positive")
  }

  mzroll_db_con <- clamr::mzroll_db_sqlite(mzroll_db_path)
  debugr::dwatch(msg = "connected to mzroll_db. [mass_spec<pipeline_wrappers.R>::pipeline_coelution_labeling]")

  # remove redundant peakgroups
  clamr_config <- clamr::build_clamr_config(dots)
  clamr::remove_redundant_peakgroups(mzroll_db_con, clamr_config, rt_grouping_tol = 0.5, group_corr_cutoff = 0.98)

  debugr::dwatch(msg = "removed redundant peak groups. [mass_spec<pipeline_wrappers.R>::pipeline_coelution_labeling]\n")

  required_mzroll_tables <- c("coelution_mz_diffs", "coelutions", "peakgroups", "peaks", "scans", "samples")
  missing_req_tables <- setdiff(required_mzroll_tables, DBI::dbListTables(mzroll_db_con))
  if (length(missing_req_tables) != 0) {
    stop(paste(missing_req_tables, collapse = " & "), " tables are missing from the .mzrollDB")
  }

  # TO DO - functionalize this coelution extraction
  coelutions <- list(
    coelutions = dplyr::tbl(mzroll_db_con, "coelutions") %>% dplyr::collect(),
    coelution_mz_diffs = dplyr::tbl(mzroll_db_con, "coelution_mz_diffs") %>% dplyr::collect()
  )
  class(coelutions) <- "coelution"

  debugr::dwatch(msg = "determined coelutions. [mass_spec<pipeline_wrappers.R>::pipeline_coelution_labeling]")

  # define search options

  peakgrp_coelution_annotation_args <- reclass_cmd_args(dots, clamr::peakgrp_coelution_annotation,
    default_args = list(
      min_edge_inconsistency_for_root = 1,
      root_weight_constant = 1.25
    )
  )
  peakgrp_coelution_annotation_args$mzroll_db_con <- mzroll_db_con
  peakgrp_coelution_annotation_args$clamr_config <- clamr_config
  peakgrp_coelution_annotation_args$coelutions <- coelutions

  debugr::dwatch(msg = "defined search options in preparation for grouping clusters.\n[mass_spec<pipeline_wrappers.R>::pipeline_coelution_labeling]")

  # group clusters of peakgroups linked via adduct / isotopologue relationships

  baseline_label <- ifelse(dots$mode == "positive", "[M+H]+", "[M-H]-")

  # Guard possible empty output
  peakgroup_coelution_annotation_results <- do.call(clamr::peakgrp_coelution_annotation, peakgrp_coelution_annotation_args)

  if (is.null(peakgroup_coelution_annotation_results)) {
    cat("Coelution labeling DONE - no coelutions detected.\n")
    invisible(0)
  }

  adduct_labels <- peakgroup_coelution_annotation_results %>%
    dplyr::ungroup() %>%
    # add baseline adduct type (M +/- H) to parent peaks
    dplyr::mutate(adductName = ifelse(adductName != "", adductName, baseline_label)) %>%
    dplyr::rename(parentGroupId_new = parentGroupId, adductName_new = adductName, metaGroupId_new = metaGroupId)

  debugr::dwatch(msg = "created adduct_labels. [mass_spec<pipeline_wrappers.R>::pipeline_coelution_labeling]")

  # Write updated adduct labels to peakgroups
  updated_peakgroups <- dplyr::tbl(mzroll_db_con, "peakgroups") %>%
    dplyr::collect() %>%
    dplyr::left_join(adduct_labels, by = "groupId") %>%
    dplyr::mutate(parentGroupId = parentGroupId_new, adductName = adductName_new, metaGroupId = metaGroupId_new) %>%
    dplyr::select(-parentGroupId_new, -adductName_new, -metaGroupId_new) %>%
    # If a group has no parent, its metaGroupId should be 0
    dplyr::mutate(metaGroupId = ifelse(parentGroupId == 0L, 0L, metaGroupId))

  debugr::dwatch(msg = "created updated_peakgroups. [mass_spec<pipeline_wrappers.R>::pipeline_coelution_labeling]")

  clamshell::sqlite_write(mzroll_db_con, "peakgroups", updated_peakgroups, overwrite = TRUE)
  DBI::dbDisconnect(mzroll_db_con)

  cat("Coelution labeling DONE - parentGroupId, adductName and metaGroupId variables of the peakgroups table have been overwritten\n")
  invisible(0)
}

#' Pipeline Direct Infusion
#'
#' Identify compounds in direct infusion data
#'
#' @param ... all arguments are passed in from mzkit.R using \code{do.call(rwrapper, as.list(formatted_flags))}
#' @return 0 invisibly, function called for side effects.
#'
#' @export
pipeline_direct_infusion <- function(...) {
  debugr::dwatch(msg = "Entered pipeline_direct_infusion function. [mass_spec<pipeline_wrappers.R>::pipeline_direct_infusion]")

  # based on peak group features (mz, rt, fragmentation) and adduct type
  # search for standard ions which match the peak well

  dots <- list(...)

  # output_folder, mzroll_db_file, MS1tol and MS2tol are required arguments which will be passed as ...
  # along with any optional arguments for identify_peakgroups

  test_output_folder(dots)
  output_folder <- normalizePath(dots$output_folder)

  debugr::dwatch(msg = "Resolved output_folder. [mass_spec<pipeline_wrappers.R>::pipeline_direct_infusion]\n")
  debugr::dwatch(msg = paste("output_folder=", output_folder))

  test_mzroll_db_file(dots)
  mzroll_db_path <- normalizePath(file.path(dots$output_folder, dots$mzroll_db_file))

  debugr::dwatch(msg = "Just before test_searchStandards_dots. [mass_spec<pipeline_wrappers.R>::pipeline_direct_infusion]")

  clamdb::test_searchStandards_dots(dots)

  debugr::dwatch(msg = "Just before matching model. [mass_spec<pipeline_wrappers.R>::pipeline_direct_infusion]\n")

  # determine which matching method to use and (if needed) find its scoring parameters
  matching_model <- setup_matching_model(dots)

  # load standard data for the mode/chemical class of interest

  mass_spec_standards_con <- DBI::dbConnect(
    RMySQL::MySQL(),
    user = dots$standard_db_user,
    password = Sys.getenv(dots$standard_db_passwd_key),
    dbname = dots$dbname,
    host = Sys.getenv(dots$standard_db_host_key)
  )

  DBI::dbSendQuery(mass_spec_standards_con, "SET SESSION WAIT_TIMEOUT=9999999")
  DBI::dbSendQuery(mass_spec_standards_con, "SET SESSION MAX_EXECUTION_TIME=9999999")

  debugr::dwatch(msg = "Succesfully defined connection to mass_spec standards [mass_spec<pipeline_wrappers.R>::pipeline_direct_infusion]")

  # test chromatography method validity based on standards data

  defined_methods <- dplyr::tbl(mass_spec_standards_con, dbplyr::sql("SELECT * FROM methods")) %>%
    dplyr::collect()
  if (!dots$chromatographic_method %in% defined_methods$chromatographicMethod) {
    stop(dots$chromatographic_method, " is not a defined value for chromatographic_method; existing methods are:\n", paste(defined_methods$chromatographicMethod, collapse = "\n"))
  }

  debugr::dwatch(msg = "determined defined methods. [mass_spec<pipeline_wrappers.R>::pipeline_direct_infusion]")

  matched_analytical_methods <- defined_methods %>%
    dplyr::filter(
      chromatographicMethod == dots$chromatographic_method,
      collisionEnergies == dots$collision_energies,
      mode == dots$mode
    )

  if (nrow(matched_analytical_methods) == 1) {
    matched_method <- paste0(gsub(" ", "_", matched_analytical_methods$analyticalMethod[1]), "-unnormalized-matched_single_energy.msp")
  } else if (nrow(matched_analytical_methods) == 0) {
    warning("no analytical methods were matched;")
    matched_method <- ""
  } else {
    warning(nrow(matched_analytical_methods), " analytical methods were matched")
    matched_method <- ""
  }

  standards_data <- get_standards_subset(mass_spec_standards_con,
    matched_method = matched_method,
    r_scripts_path = clamr_config$r_scripts_path,
    chemical_class = chemical_class,
    ms_mode = clamr_config$mode,
    chromatographic_method = clamr_config$chromatographic_method
  )

  # done with standards
  DBI::dbDisconnect(mass_spec_standards_con)

  # TODO: retrieve samples from clamr_config
  samples <- c("")
  direct_infusion_results <- mzkitcpp::direct_infusion_results(samples, standards_data, TRUE)

  # TODO: do something with results

  cat("Direct Infusion analysis DONE\n")
  invisible(0)
}

#' Pipeline QC
#'
#' Create a QC directory containing results from internal and external standards.
#'
#' @param ... all arguments are passed in from mzkit.R using \code{do.call(rwrapper, as.list(formatted_flags))}
#' @return 0 invisibly, function called for side effects.
#'
#' @export
pipeline_qc <- function(...) {
  dots <- list(...)

  test_output_folder(dots)
  output_folder <- normalizePath(dots$output_folder)

  # test for appropriate configuration syntax for the standards database
  validate_standard_search_environment(dots)

  # connect to standards database
  mass_spec_standards_con <- DBI::dbConnect(
    RMySQL::MySQL(),
    user = dots$standard_db_user,
    password = Sys.getenv(dots$standard_db_passwd_key),
    dbname = dots$dbname,
    host = Sys.getenv(dots$standard_db_host_key)
  )

  # validate mzrollDB
  test_mzroll_db_file(dots)
  mzroll_db_path <- normalizePath(file.path(dots$output_folder, dots$mzroll_db_file))
  mzroll_db_con <- clamr::mzroll_db_sqlite(mzroll_db_path)

  if (is.null(dots$referenceSample)) {
    warning("No reference sample specified! Ensure that the reference sample field exists
            in your configuration file ('referenceSample=x') in the list of
            'rwrapper_parameters' under 'pipeline_qc'.
            The basepeak plot will not be generated.")
    dots$referenceSample <- ""
  }

  if (!file.exists(dots$referenceSample)) {
    warning(
      paste("The provided reference sample file provided does not exist:'
                  ", dots$referenceSample, "'\n"),
      "The basepeak plot will not be generated."
    )
    dots$referenceSample <- ""
  }

  # format config
  clamr_config <- clamr::build_clamr_config(dots)

  # extract samples
  before <- Sys.time()

  samples_df <- dplyr::tbl(mzroll_db_con, dbplyr::sql("SELECT name as samplename, filename FROM samples")) %>%
    dplyr::collect() %>%
    # remove extension
    dplyr::mutate(samplename = stringr::str_remove(samplename, "\\.[a-zA-Z]+$")) %>%
    # dplyr::sample_n(10) %>% # debugging
    {
      .
    }

  DBI::dbDisconnect(mzroll_db_con)

  after <- Sys.time()
  debugr::dwatch(msg = paste("extracted samples for qc in", difftime(after, before, units = "secs"), "seconds. [mass_spec<pipeline_wrappers.R>::pipeline_qc]"))

  # verify that sample raw files can be found
  missing_files <- samples_df$filename[!file.exists(samples_df$filename)]
  if (length(missing_files) != 0) {
    stop(length(missing_files), " files could not be found: ", paste(missing_files, collapse = ", "))
  }

  # extract data for all internal and external standards
  before <- Sys.time()

  extracted_qc_standards <- clamqc::qc_extract_standards(samples_df,
    clamr_config,
    mzroll_db_path,
    mass_spec_standards_con,
    query_methodId = dots$methodId
  )

  after <- Sys.time()
  debugr::dwatch(msg = paste("extracted qc standards in", difftime(after, before, units = "secs"), "seconds. [mass_spec<pipeline_wrappers.R>::pipeline_qc]"))

  # summaries & plots to report

  # standard numerical summaries
  before <- Sys.time()

  standard_abundances <- clamqc::qc_standard_summaries(extracted_qc_standards)

  after <- Sys.time()
  debugr::dwatch(msg = paste("computed standard_abundances in", difftime(after, before, units = "secs"), "seconds. [mass_spec<pipeline_wrappers.R>::pipeline_qc]"))

  # standard summary plots
  before <- Sys.time()

  summary_plots <- extracted_qc_standards %>%
    # generate plots for each method separately
    dplyr::mutate(summary_plots = purrr::map(standard_data, .f = clamqc::qc_standard_plots)) %>%
    dplyr::select(-standard_data) %>%
    tidyr::unnest(summary_plots)

  after <- Sys.time()
  debugr::dwatch(msg = paste("created summary_plots in", difftime(after, before, units = "secs"), "seconds. [mass_spec<pipeline_wrappers.R>::pipeline_qc]"))

  # summarize basepeak

  if (dots$referenceSample != "") {
    before <- Sys.time()
    basepeak_plot <- clamqc::basepeak_comparison(samples_df, clamr_config$referenceSample)

    after <- Sys.time()
    debugr::dwatch(msg = paste("created basepeak_plot in", difftime(after, before, units = "secs"), "seconds. [mass_spec<pipeline_wrappers.R>::pipeline_qc]"))
  } else {
    warning("basepeak generation skipped")
  }

  # saving reports
  before <- Sys.time()

  # Create the QC directory if it does not already exist.
  qc_subfolder <- file.path(output_folder, "QC")
  if (!dir.exists(qc_subfolder)) {
    dir.create(qc_subfolder)
  }

  readr::write_tsv(standard_abundances, file.path(qc_subfolder, "standard_abundances.tsv"))

  debugr::dwatch(msg = "write_tsv completed. [mass_spec<pipeline_wrappers.R>::pipeline_qc]")

  summary_plots %>%
    dplyr::mutate(label = paste(plot_type, standard, sep = "_")) %>%
    {
      purrr::walk2(.$grob, .$label, function(grob, label, output_folder) {
        labelled_plot <- grob + ggplot2::ggtitle(stringr::str_replace_all(label, "[_-]", " "))
        ggplot2::ggsave(labelled_plot, filename = file.path(output_folder, "QC", paste0(label, ".pdf")), height = 8, width = 10)
      }, output_folder = output_folder)
    }

  if (dots$referenceSample != "") {
    ggplot2::ggsave(basepeak_plot, filename = file.path(output_folder, "QC", "basepeak.pdf"), height = 12, width = 12)
  }

  after <- Sys.time()
  debugr::dwatch(msg = paste("saved all qc data and plots in", difftime(after, before, units = "secs"), "seconds. [mass_spec<pipeline_wrappers.R>::pipeline_qc]"))

  cat("Standard QC is DONE - QC summaries added to QC folder\n")
  invisible(0)
}

#' Pipeline EDA
#'
#' Save heatmaps of all compounds and just identified compounds
#'
#' @param ... all arguments are passed in from mzkit.R using \code{do.call(rwrapper, as.list(formatted_flags))}
#' @return 0 invisibly, function called for side effects.
#'
#' @export
pipeline_eda <- function(...) {
  dots <- list(...)

  # output_folder, mzroll_db_file, MS1tol and MS2tol are required arguments which will be passed as ...
  # along with any optional arguments ms2_driven_aligner

  test_output_folder(dots)
  output_folder <- normalizePath(dots$output_folder)

  test_mzroll_db_file(dots)
  mzroll_db_path <- normalizePath(file.path(dots$output_folder, dots$mzroll_db_file))

  # heatmap of knowns
  mzroll_list_identified <- calicomics::process_mzroll(mzroll_db_path,
    standard_databases = NULL,
    only_identified = TRUE
  )

  debugr::dwatch(msg = "determined mzroll_list_identified. [mass_spec<pipeline_wrappers.R>::pipeline_eda]")

  if (is.null(mzroll_list_identified) || mzroll_list_identified$peakgroups %>% nrow() < 2 || mzroll_list_identified$peaks %>% nrow() < 2) {
    cat("Not enough data to generate heatmap.\n")
    return(invisible(0))
  }

  mzroll_list_identified$samples <- mzroll_list_identified$samples %>%
    dplyr::mutate(sampleName = calicomics::remove_constant_name(name))

  identified_heatmap <- calicomics::plot_heatmap(mzroll_list_identified, sample.var = "sampleName", feature.var = "peak_label", plot_type = "grob")

  # Create the reports subdirectory if it does not already exist.
  reports_subfolder <- file.path(output_folder, "reports")
  if (!dir.exists(reports_subfolder)) {
    dir.create(reports_subfolder)
  }

  ggplot2::ggsave(identified_heatmap, filename = file.path(reports_subfolder, paste0("heatmap_identified.pdf")), height = 25, width = 15)

  # heatmap of all peakgroups
  mzroll_list_all <- calicomics::process_mzroll(mzroll_db_path,
    standard_databases = NULL,
    only_identified = FALSE
  )
  mzroll_list_all$samples <- mzroll_list_all$samples %>%
    dplyr::mutate(sampleName = calicomics::remove_constant_name(name))

  all_heatmap <- calicomics::plot_heatmap(mzroll_list_all, sample.var = "sampleName", feature.var = "groupId", plot_type = "grob")
  ggplot2::ggsave(all_heatmap, filename = file.path(reports_subfolder, paste0("heatmap_all.pdf")), height = 20, width = 15)

  cat("Summary heatmaps added to reports\n")
  invisible(0)
}
