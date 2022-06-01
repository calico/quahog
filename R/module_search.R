#' Pipeline Standard Search
#'
#' Identify peakgroups which match known compounds
#'
#' Using a matching model compare experimental features of peakgroups to candidate compounds which they may match. The mzrollDB file is updated
#' by annotating peakgroups with their best match and adding all matches
#' as a separate table. Feasible matches are written out as a .msp file
#' to support manual validation in MAVEN.
#'
#' @param ... all arguments are passed in from mzkit.R using \code{do.call(rwrapper, as.list(formatted_flags))}
#'
#' @return 0 invisibly, function called for side effects.
#'
#' @export
pipeline_standard_search <- function(...) {
  debugr::dwatch(msg = "Entered pipeline_standard_search function. [mass_spec<pipeline_wrappers.R>::pipeline_standard_search]")

  # based on peak group features (mz, rt, fragmentation) and adduct type
  # search for standard ions which match the peak well

  dots <- list(...)

  # output_folder, mzroll_db_file, MS1tol and MS2tol are required arguments which will be passed as ...
  # along with any optional arguments for identify_peakgroups

  test_output_folder(dots)
  output_folder <- normalizePath(dots$output_folder)

  debugr::dwatch(msg = "Resolved output_folder. [mass_spec<pipeline_wrappers.R>::pipeline_standard_search]")
  debugr::dwatch(msg = paste("output_folder=", output_folder))

  test_mzroll_db_file(dots)
  mzroll_db_path <- normalizePath(file.path(dots$output_folder, dots$mzroll_db_file))

  debugr::dwatch(msg = "Just before test_searchStandards_dots. [mass_spec<pipeline_wrappers.R>::pipeline_standard_search]")

  test_searchStandards_dots(dots)

  debugr::dwatch(msg = "Just before matching model. [mass_spec<pipeline_wrappers.R>::pipeline_standard_search]")

  # determine which matching method to use and (if needed) find its scoring parameters
  matching_model <- setup_matching_model(dots)

  # load standard data for the mode/chemical class of interest

  # Hook to use cached .rds file instead of accessing external standards database. Used for CI testing.
  if (grepl("\\.rds$", dots$dbname)) {
    matched_method <- gsub(".rds$", ".msp", dots$dbname)
    mass_spec_standards_con <- "use-rds-file"
  } else {

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

    sql_result <- DBI::dbSendQuery(mass_spec_standards_con, "SET SESSION WAIT_TIMEOUT=9999999")
    sql_result <- DBI::dbSendQuery(mass_spec_standards_con, "SET SESSION MAX_EXECUTION_TIME=9999999")

    debugr::dwatch(msg = "Succesfully defined connection to mass_spec standards [mass_spec<pipeline_wrappers.R>::pipeline_standard_search]")

    # test chromatography method validity based on standards data

    matched_method <- clamdb::match_method(dots, mass_spec_standards_con)
  }

  ## Matching ##

  # connect to sql-lite database and setup config

  mzroll_db_con <- clamr::mzroll_db_sqlite(mzroll_db_path)

  clamr_config <- clamr::build_clamr_config(dots)

  before <- Sys.time()
  # align features of experimental and standards
  peakgroup_standard_features <- clamdb::align_features_to_standards(mzroll_db_con, clamr_config, peakgroups = NULL, mass_spec_standards_con, matched_method)

  after <- Sys.time()
  matched_time <- paste("Matched experimental data to standards db in", difftime(after, before, units = "secs"), "seconds\n")

  debugr::dwatch(msg = paste(matched_time, "[mass_spec<pipeline_wrappers.R>::pipeline_standard_search]"))

  cat(matched_time)

  # Option to exclude isotopes from identification
  if ("exclude_isotopes" %in% names(dots) && dots$exclude_isotopes == "T") {
    isotope_adducts <- clamr::standard_coelution_library %>%
      dplyr::filter(type == "isotopologue")

    before <- Sys.time()

    peakgroup_standard_features <- peakgroup_standard_features %>%
      dplyr::filter(!stringr::str_detect(adductName, isotope_adducts$name))

    after <- Sys.time()

    debugr::dwatch(msg = paste("Filtered out isotopologue peak groups in", difftime(after, before, units = "secs"), "seconds\n"))
  } else {
    debugr::dwatch(msg = "Peak groups designated as isotopologues will be searched against the standards database.")
  }

  if (any(class(mass_spec_standards_con) == "MySQLConnection")) {
    # done with standards
    DBI::dbDisconnect(mass_spec_standards_con)
  }

  if (matching_model$score_method == "lipids_fast_match") {
    debugr::dwatch(msg = "using lipids_fast_match model for analyte identification and scoring.")

    before <- Sys.time()

    scored_matches_idsOnly <- mzkitcpp::scored_matches_lipids_rcpp(peakgroup_standard_features, dots, FALSE, 1, 0.30)

    match_scores_complete <- scored_matches_idsOnly$match_scores %>%
      dplyr::inner_join(peakgroup_standard_features, by = c("groupId" = "groupId", "ionId" = "ionId")) %>%
      dplyr::mutate(originalCompoundName = compoundName) %>%
      dplyr::mutate(compoundName = summarizedCompoundName) %>%
      dplyr::select(-summarizedCompoundName) %>%
      dplyr::mutate(adductName = adductName.x) %>%
      dplyr::select(-adductName.x, -adductName.y)

    top1_matches_complete <- scored_matches_idsOnly$top1_matches %>%
      dplyr::inner_join(peakgroup_standard_features, by = c("groupId" = "groupId", "ionId" = "ionId")) %>%
      dplyr::mutate(originalCompoundName = compoundName) %>%
      dplyr::mutate(compoundName = summarizedCompoundName) %>%
      dplyr::select(-summarizedCompoundName) %>%
      dplyr::mutate(adductName = adductName.x) %>%
      dplyr::select(-adductName.x, -adductName.y)

    top1_discoveries_complete <- top1_matches_complete %>%
      dplyr::filter(is_match == TRUE)

    scored_matches <- list(
      "match_scores" = match_scores_complete,
      "top1_matches" = top1_matches_complete,
      "top1_discoveries" = top1_discoveries_complete
    )

    after <- Sys.time()

    debugr::dwatch(msg = paste("Determined scored matches for lipids in", difftime(after, before, units = "secs"), "s [mass_spec<pipeline_wrappers.R>::pipeline_standard_search]"))
  } else {
    before <- Sys.time()

    # generate features which compare similarity of experimental and standard full scan and fragmentation features.
    comparison_features <- clamdb::compare_aligned_features(peakgroup_standard_features,
      clamr_config$MS2,
      n_top_spectra_summed = 3L,
      quality_weights = c("purity" = 2, "quality" = 1),
      frag_similarity_methods = "all",
      minimum_ic_fraction = 0.005
    ) %>%
      tidyr::unnest(qc_features) %>%
      tidyr::unnest(fullscan_features) %>%
      tidyr::unnest(fragmentation_features)

    after <- Sys.time()
    debugr::dwatch(msg = paste("produced comparison features in", difftime(after, before, units = "secs"), "s [mass_spec<pipeline_wrappers.R>::pipeline_standard_search]"))

    # scoring comparison features to identify metabolite matches
    scored_matches <- do.call(clamdb::score_standard_matches, append(list(grouping_vars = "groupId", comparison_features = comparison_features), matching_model))

    debugr::dwatch(msg = "determined scored_matches. [mass_spec<pipeline_wrappers.R>::pipeline_standard_search]")

    # Add these columns for schema agreement with lipidomics samples
    scored_matches$match_scores <- scored_matches$match_scores %>%
      dplyr::mutate(originalCompoundName = compoundName) %>%
      dplyr::mutate(summarizationLevel = 4)

    # Use adduct name from standards instead of from coelution labeling
    scored_matches$match_scores <- scored_matches$match_scores %>%
      dplyr::mutate(adductName = sapply(stdCompounds, function(l) {
        l$adductName
      }))

    scored_matches$top1_matches <- scored_matches$top1_matches %>%
      dplyr::mutate(adductName = sapply(stdCompounds, function(l) {
        l$adductName
      }))

    scored_matches$top1_discoveries <- scored_matches$top1_discoveries %>%
      dplyr::mutate(adductName = sapply(stdCompounds, function(l) {
        l$adductName
      }))
  }

  # Pass along all matches as discoveries (mzkitchen 10 test)
  # This line can be commented out to reinstate the previous threshold.
  scored_matches$top1_discoveries <- scored_matches$top1_matches

  # Label entries above threshold
  scored_matches$top1_discoveries <- scored_matches$top1_discoveries %>%
    dplyr::mutate(label = ifelse(score >= 0.809, "l", ""))

  # matches performed
  match_scores <- scored_matches$match_scores %>%
    dplyr::mutate(
      compoundName = ifelse(is.na(smiles), compoundName, paste(compoundName, smiles, sep = ": ")),
      compoundDB = matched_method
    ) %>%
    dplyr::select(groupId, ionId, compoundName, adductName, summarizationLevel, originalCompoundName, score, is_match)

  cat(nrow(match_scores), "matches scored\n")

  # identifications made
  identity_calls <- scored_matches$top1_discoveries %>%
    dplyr::arrange(compoundName) %>%
    dplyr::mutate(
      compoundName_new = ifelse(is.na(smiles), compoundName, paste(compoundName, smiles, sep = ": ")),
      compoundDB_new = matched_method
    )

  cat(nrow(identity_calls), "peaks identified\n")

  # msp export for lipids is now handled by mzkitcpp::scored_matches_lipids_rcpp
  if (matching_model$score_method != "lipids_fast_match") {
    before <- Sys.time()
    save_top_ion_msp(peakgroup_standard_features, scored_matches, output_folder, max_reported_ions = Inf, normalize_intensities = FALSE)
    after <- Sys.time()
    debugr::dwatch(msg = paste("Exported msp files in ", difftime(after, before, units = "secs"), "s [mass_spec<pipeline_wrappers.R>::pipeline_standard_search]"))
  }

  # update peakgroups with new IDs

  # Use the adduct names associated with the identified compounds instead of coelution scoring-derived adduct assignments.
  updated_peakgroups <- dplyr::tbl(mzroll_db_con, dbplyr::sql("SELECT * FROM peakgroups")) %>%
    dplyr::collect() %>%
    dplyr::left_join(identity_calls %>%
      dplyr::select(groupId,
        ms2Score_new = score,
        compoundName_new,
        compoundId_new = ionId,
        compoundDB_new,
        adductName_new = adductName,
        label_new = label
      ) %>%
      dplyr::mutate(compoundId_new = as.character(compoundId_new)),
    by = "groupId"
    ) %>%
    dplyr::mutate(searchTableName = ifelse(!is.na(compoundName_new), "clamDB", searchTableName)) %>%
    dplyr::mutate(
      ms2Score = ms2Score_new,
      compoundName = compoundName_new,
      compoundId = compoundId_new,
      compoundDB = compoundDB_new,
      label = label_new
    ) %>%
    dplyr::mutate(adductName = ifelse(!is.na(adductName_new), adductName_new, adductName)) %>%
    dplyr::select(-ms2Score_new, -compoundName_new, -compoundId_new, -compoundDB_new, -adductName_new, -label_new)

  peakgroups_no_parent <- updated_peakgroups %>% dplyr::filter(parentGroupId == 0)
  child_peakgroups <- updated_peakgroups %>% dplyr::filter(parentGroupId != 0)

  # If no children, no need for any re-numbering.
  if (nrow(child_peakgroups) > 0) {
    parent_peakgroups <- peakgroups_no_parent %>%
      dplyr::filter(groupId %in% child_peakgroups$parentGroupId) %>%
      dplyr::mutate(parentGroupId = groupId) # assign the parent to itself temporarily

    parents_and_children <- rbind(child_peakgroups, parent_peakgroups) %>%
      dplyr::arrange(parentGroupId, -ms2Score)

    all_others <- updated_peakgroups %>% dplyr::filter(!groupId %in% parents_and_children$groupId)

    lastParentGroupId <- as.integer(-1)
    clusterGroupId <- as.integer(-1)
    for (i in 1:nrow(parents_and_children)) {
      currentParentGroupId <- parents_and_children$parentGroupId[i]

      if (lastParentGroupId != currentParentGroupId) {
        parents_and_children$parentGroupId[i] <- 0L
        clusterGroupId <- parents_and_children$groupId[i]
      } else {
        parents_and_children$parentGroupId[i] <- clusterGroupId
      }

      lastParentGroupId <- currentParentGroupId
    }

    # Keep metaGroupId to always be 0 to compare/contrast effect of MAVEN-based clustering
    corrected_peak_groups <- rbind(parents_and_children, all_others) %>%
      dplyr::mutate(metaGroupId = 0L)
  } else {
    corrected_peak_groups <- updated_peakgroups %>% dplyr::mutate(metaGroupId = 0L)
  }

  # overwrite peakgroups with new identifications
  clamshell::sqlite_write(mzroll_db_con, "peakgroups", corrected_peak_groups, overwrite = TRUE)

  # add all matches performed to .mzrollDB
  clamshell::sqlite_write(mzroll_db_con, "matches", match_scores, overwrite = TRUE)

  # done with mzrollDB
  DBI::dbDisconnect(mzroll_db_con)

  cat("Standard search DONE - compoundName, compoundId, compoundDB, ms2Score variables of the peakgroups table have been overwritten\n")
  invisible(0)
}

setup_matching_model <- function(dots) {
  stopifnot("r_scripts_path" %in% names(dots))

  if (!("matching_model" %in% names(dots))) {
    warning("matching_model not provided in config, baseline matching method will be used")
    matching_model <- list(score_method = "baseline")
  } else if (dots$matching_model == "baseline") {
    matching_model <- list(score_method = "baseline")
  } else if (dots$matching_model == "lipids_fast_match") {
    matching_model <- list(score_method = "lipids_fast_match")
  } else {
    asset_path <- do.call(file.path, as.list((c(stringr::str_split(dots$r_scripts_path, pattern = .Platform$file.sep)[[1]] %>%
      {
        .[-length(.)]
      }, "open_CLaM", "matching_models", dots$matching_model))))

    if (file.exists(dots$matching_model)) {
      matching_model <- readRDS(dots$matching_model)
    } else if (file.exists(asset_path)) {
      matching_model <- readRDS(asset_path)
    } else {
      warning("matching_model could not be found from config, baseline matching method will be used")
      matching_model <- list(score_method = "baseline")
    }
  }
  stopifnot("list" %in% class(matching_model))
  stopifnot("score_method" %in% names(matching_model))

  matching_model
}

validate_standard_search_environment <- function(dots) {

  # check for required dots

  required_dots <- c(
    "dbname",
    "standard_db_user",
    "standard_db_passwd_key",
    "standard_db_host_key"
  )

  missing_required_dots <- setdiff(required_dots, names(dots))

  if (length(missing_required_dots) != 0) {
    stop(
      length(missing_required_dots),
      " variables are required and must be specified within the config:\n",
      paste(missing_required_dots, collapse = "\n")
    )
  }

  # check for existence of key_variables
  error_message <- glue::glue()

  missing_standard_pwd <- Sys.getenv(dots$standard_db_passwd_key) == ""
  if (missing_standard_pwd) {
    error_message <- glue::glue(
      error_message,
      glue::glue("standard database password is missing,
                 {dots$standard_db_passwd_key} should be added to your .Renviron
                 with an appropriate password")
    )
  }

  missing_standard_host <- Sys.getenv(dots$standard_db_host_key) == ""
  if (missing_standard_host) {
    error_message <- glue::glue(
      error_message,
      glue::glue("standard database host is missing,
                 {dots$standard_db_host_key} should be added to your .Renviron
                 with an appropriate host url")
    )
  }

  if (missing_standard_pwd | missing_standard_host) {
    stop(error_message)
  }

  invisible(0)
}

save_top_ion_msp <- function(peakgroup_standard_features, scored_matches, output_folder, max_reported_ions = 10000, normalize_intensities = TRUE) {
  if (!is.infinite(max_reported_ions)) {
    top_ions <- scored_matches$match_scores %>%
      dplyr::group_by(ionId) %>%
      dplyr::arrange(desc(score)) %>%
      dplyr::slice(1) %>%
      dplyr::ungroup() %>%
      dplyr::slice(1:max_reported_ions) %>%
      dplyr::select(ionId) %>%
      {
        .$ionId
      }
  } else {
    top_ions <- unique(scored_matches$match_scores$ionId)
  }

  unique_top_ion_features <- peakgroup_standard_features %>%
    dplyr::filter(ionId %in% top_ions) %>%
    dplyr::group_by(ionId) %>%
    dplyr::slice(1) %>%
    dplyr::select(-groupId, -parentGroupId, -adductName, -peaks, -scans) %>%
    dplyr::ungroup()

  clamdb::generate_msp(unique_top_ion_features, normalize_intensities = normalize_intensities, minimum_ic_fraction = 0.001, save_name = "mzkit", dir_path = file.path(output_folder, "libraries"))
}

test_searchStandards_dots <- function(dots) {

  # verify that configuration contains appropriate (and appropriately) formatted fields for searching standards

  missing_instrument_tolerances <- setdiff(c("MS1tol", "MS2tol"), names(dots))
  if (length(missing_instrument_tolerances) != 0) {
    stop("Instrument tolerances must be specified with ", paste(missing_instrument_tolerances, collapse = ", "))
  }

  if (!"mode" %in% names(dots) || !dots$mode %in% c("negative", "positive")) {
    stop("mode is a required flag when calling searchStandards; valid options are negative and positive")
  }
  if (!"chemical_class" %in% names(dots) || !dots$chemical_class %in% c("polar", "lipid", "fatty acid")) {
    stop("chemical_class is a required flag when calling searchStandards; valid options are polar, lipid, fatty acid")
  }
  if (!"collision_energies" %in% names(dots) || all(class(dots$collision_energies) != "character") || length(dots$collision_energies) != 1) {
    stop("collision_energies is a required flag when calling searchStandards; format example: \"20,40,80\"")
  }
  if (!"chromatographic_method" %in% names(dots) || all(class(dots$collision_energies) != "character") || length(dots$collision_energies) != 1) {
    stop("chromatographic_method is a required flag when calling searchStandards; format example: \"MS-Chrom-001-A C18-TBA\"")
  }

  invisible(TRUE)
}
