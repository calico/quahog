#' Format Flagged Arguments
#'
#' Convert a command-line input / json input a named list.
#'
#' @param input_args an unnamed vector of arguments of the form "mode=negative"
#'
#' @return a named list
#'
#' @examples
#' format_flagged_arguments(c("MS1tol=10ppm", "MS2tol=20ppm", "mode=negative", "fu='b a r'"))
#'
#' @export
format_flagged_arguments <- function(input_args) {
  flagged_args <- input_args[grepl("=", input_args)]

  if (length(flagged_args) == 0) {
    stop("No arguments detected, pass arguments with variable=value")
  }

  formatted_flags <- lapply(flagged_args, function(x) {
    an_arg <- strsplit(x, split = "=")[[1]]
    if (length(an_arg) != 2) {
      stop(x, " is misformatted")
    }
    out <- an_arg[2]
    names(out) <- an_arg[1]
    out
  })

  # generate a named vector of flags
  formatted_flags <- do.call("c", formatted_flags) %>%
    # remove quotes
    {
      structure(stringr::str_replace_all(., '[\'\"]', ""), names = names(.))
    }

  # test for input duplication
  if (any(names(formatted_flags) == "")) {
    stop("some input flags did not contain a variable")
  }
  if (length(names(formatted_flags)) != length(unique(names(formatted_flags)))) {
    stop("some input flags were duplicated")
  }

  formatted_flags
}

#' Reclass commmand line args
#'
#' Using a list of character arguments (generated from parsing class-less command line inputs),
#' change the class of inputs based on default arguments, and class specifications, then pad the list with the
#' unspecified defaults.
#'
#' @param chr_list a named list of specified arguments of class character.
#' @param call_fxn an unquoted function which is meant to be called with the processed arguments.
#' @param default_args a list containing default values of arguments to be used if a value is not specified in \code{chr_list}.
#' @param other_args a named vector containing the classes of other arguments (besides \code{default_args}) whose class should be set.
#'
#' @return a named list of arguments to be supplied to \code{call_fxn}
#'
#' @details Currently vector arguments are not supported.
#'
#' @examples
#' reclass_cmd_args(
#'   chr_list = list(cosine_cutoff = "0.96", print_plots = "TRUE", show_specific_samples = "1", fu = "bar"),
#'   call_fxn = clamr::ms2_driven_aligner,
#'   default_args = list(group_by_charge = FALSE, maximum_mz_set_size = 200L, cosine_cutoff = 0.95, spline_ridge_penalty = 200, print_plots = TRUE),
#'   other_args = c(show_specific_samples = "integer")
#' )
#'
#' @export
reclass_cmd_args <- function(chr_list, call_fxn, default_args = NULL, other_args = NULL) {
  if (length(intersect(names(default_args), names(other_args))) != 0) {
    stop(paste(intersect(names(default_args), names(other_args)), collapse = " & "), " are specified in both
          default_args and other_args; either include a default value or class, not both")
  }

  chr_list_matching_args <- chr_list[intersect(names(chr_list), names(formals(call_fxn)))]

  # class changes based on class of matching default args
  default_arg_classes <- purrr::map(default_args, class)
  reclassed_w_default <- chr_list_matching_args[names(chr_list_matching_args) %in% names(default_args)]
  reclassed_w_default_classes <- default_arg_classes[names(chr_list_matching_args)[names(chr_list_matching_args) %in% names(default_args)]]
  reclassed_w_default <- purrr::map2(reclassed_w_default, reclassed_w_default_classes, function(x, y) {
    class(x) <- y
    x
  })

  # class changes based class specifications in other_args

  reclassed_w_other <- chr_list_matching_args[names(chr_list_matching_args) %in% names(other_args)]
  reclassed_w_other_classes <- other_args[names(chr_list_matching_args)[names(chr_list_matching_args) %in% names(other_args)]]
  reclassed_w_other <- purrr::map2(reclassed_w_other, reclassed_w_other_classes, function(x, y) {
    class(x) <- y
    x
  })

  # leave non-reclassed args as is

  nonreclassed_args <- chr_list_matching_args[setdiff(names(chr_list_matching_args), union(names(reclassed_w_default), names(reclassed_w_other)))]

  # combine all args

  chr_list_matching_args <- append(nonreclassed_args, append(reclassed_w_default, reclassed_w_other))
  chr_list_defaults <- default_args[setdiff(names(default_args), names(chr_list_matching_args))]

  append(chr_list_matching_args, chr_list_defaults)
}
