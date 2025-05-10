# configure_style_lint_repo.R
#
# This script configures lintr for line length and other common linters,
# then styles and lints R and Rmd files in the current directory and its subdirectories.
# Output can be logged to a file.
#
# USAGE:
# 1. Save this file in the root of your Git repository.
# 2. Open R.
# 3. Set your working directory to the root of your Git repository:
#    setwd("path/to/your/repo")
# 4. Run the script:
#    source("configure_style_lint_repo.R")
#
# ---

# --- Initial Setup ---
cat("--- Script Initializing ---\n")
LOG_OUTPUT_TO_FILE <- TRUE # SET TO FALSE if you only want console output
LOG_FILE_NAME <- "code_quality_log.txt" # Base name for the log file
log_con <- NULL
actual_log_file_name <- LOG_FILE_NAME

if (LOG_OUTPUT_TO_FILE) {
  if (file.exists(LOG_FILE_NAME)) {
    timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
    actual_log_file_name <- paste0(tools::file_path_sans_ext(LOG_FILE_NAME), "_", timestamp, ".txt")
  }
  cat(paste0("Output will be logged to: ", actual_log_file_name, "\n"))

  tryCatch(
    {
      log_con <- file(actual_log_file_name, open = "wt")
      sink(log_con, type = "output")
      sink(log_con, type = "message", append = TRUE)
    },
    error = function(e) {
      LOG_OUTPUT_TO_FILE <<- FALSE
      log_con <<- NULL
      base::cat(paste0("CRITICAL ERROR: Could not open log file '", actual_log_file_name, "'. Logging to file is disabled.\nError: ", e$message, "\nProceeding with console output only.\n"))
    }
  )
}

# --- IMPORTANT: `styler` vs `lintr` ---
# (Explanation kept the same)
cat("--- IMPORTANT: `styler` vs `lintr` ---\n")
cat("This script uses two tools for code quality:\n")
cat("1. `styler`: Reformats your code LAYOUT (spacing, indentation, line breaks around operators).\n")
cat("   It CHANGES your files directly if their layout doesn't match its style guide.\n")
cat("   If `styler` reports 'File unchanged', its layout was already compliant.\n")
cat("2. `lintr`: Checks for a WIDER range of style issues (e.g., variable names, commented code),\n")
cat("   potential bugs, and code complexity. These are often the 'underlines' in your editor.\n")
cat("   `lintr` DOES NOT automatically change your files. This script WILL let you CONFIGURE `lintr`\n")
cat("   to ignore specific checks project-wide (via a .lintr file), which can reduce 'underlines'.\n")
cat("   Otherwise, lintr issues require manual code changes.\n")
cat("------------------------------------------\n\n")


# Function to ensure packages are installed and loaded
ensure_packages <- function(packages) {
  for (pkg in packages) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      cat(paste0("Installing package: ", pkg, "...\n"))
      install.packages(pkg)
    }
    suppressPackageStartupMessages(library(pkg, character.only = TRUE))
  }
}

required_packages <- c("styler", "lintr")
cat("Ensuring required packages are installed and loaded...\n")
tryCatch(
  {
    ensure_packages(required_packages)
    cat("Packages checked and loaded.\n")
  },
  error = function(e) {
    cat("Error during package management. Please check your R environment.\n")
    print(e)
    if (LOG_OUTPUT_TO_FILE && !is.null(log_con)) {
      if (sink.number(type = "message") > 0) sink(type = "message")
      if (sink.number(type = "output") > 0) sink(type = "output")
      if (inherits(log_con, "connection") && isOpen(log_con)) close(log_con)
      base::cat(paste0("\nError during package setup. Log (if active and created) is at: ", file.path(getwd(), actual_log_file_name), "\n"))
    }
    stop("Halting script due to package management error.")
  }
)

cat("\n--- Step 0: Configure Linters for the Project (via .lintr file) ---\n")
dot_lintr_file <- file.path(".", ".lintr")

# Temporary storage for user inputs during interactive session
# This is because we will sink back to file *after* all readline prompts
# in configure_lintr_interactive are done.
temporary_user_inputs_for_log <- character(0)

# Main function for interactive linter configuration
configure_lintr_interactive <- function(current_log_con_arg, current_log_output_to_file_arg) {
  # --- SINK MANAGEMENT: Revert to console for all interactive prompts in this function ---
  output_was_sinking_before_interactive <- FALSE
  message_was_sinking_before_interactive <- FALSE

  if (current_log_output_to_file_arg && !is.null(current_log_con_arg)) {
    if (sink.number(type = "output") > 0) {
      output_was_sinking_before_interactive <- TRUE
      sink(type = "output") # Revert to console
    }
    if (sink.number(type = "message") > 0) {
      message_was_sinking_before_interactive <- TRUE
      sink(type = "message") # Revert to console
    }
  }
  # All cat() and readline() from here until end of prompts will go to CONSOLE

  cat("--- Lintr Configuration ---\n")
  cat("Let's configure specific linters for your project.\n")
  # (Keep introductory messages the same)
  cat("For each, you can set a value, disable it, or keep its default behavior.\n")
  cat("Answering with just [Enter] usually picks the default.\n")
  cat("NOTE: If your typed characters appear as '$' or other strange symbols,\n")
  cat("      and direct readline() tests work, try simplifying sink management further\n")
  cat("      or testing in a different R console interface.\n")


  linters_to_manage <- list(
    "line_length_linter" = list(
      description = "Maximum line length for code.",
      q_text_configure = "Enter max line length (e.g., 100, 120)",
      type = "length_or_disable", default_value = 120
    ),
    "commented_code_linter" = list(
      description = "Flags commented-out R code (e.g., # old_function_call()).",
      q_text = "Disable warnings for commented-out R code?",
      type = "yes_no_disable", default_string_if_active = "no (linter active)"
    ),
    "object_name_linter" = list(
      description = "Checks variable/function naming styles (default prefers 'snake_case' or 'dotted.case').",
      q_text = "Disable warnings for object naming conventions?",
      type = "yes_no_disable", default_string_if_active = "no (linter active)"
    ),
    "object_length_linter" = list(
      description = "Checks if variable/function names are too long.",
      q_text_configure = "Enter max object name length (e.g., 30, 40)",
      type = "length_or_disable", default_value = 30
    ),
    "T_and_F_symbol_linter" = list(
      description = "Ensures TRUE/FALSE are used instead of T/F.",
      q_text = "Disable warnings for using T/F instead of TRUE/FALSE?",
      type = "yes_no_disable", default_string_if_active = "no (linter active)"
    ),
    "seq_linter" = list(
      description = "Warns against 1:length(x), preferring seq_along(x) or seq_len(length(x)).",
      q_text = "Disable warnings for 1:length(x) style loops?",
      type = "yes_no_disable", default_string_if_active = "no (linter active)"
    )
  )

  configured_linter_settings <- character()
  temporary_user_inputs_for_log <<- character(0) # Clear any previous temp logs

  for (linter_name in names(linters_to_manage)) {
    linter_details <- linters_to_manage[[linter_name]]

    cat(paste0("\n--- Configuring: ", linter_name, " ---\n")) # To console
    cat(linter_details$description, "\n") # To console

    user_choice_raw <- ""
    action_taken_msg_console <- ""
    linter_setting_for_file <- NA_character_

    if (linter_details$type == "length_or_disable") {
      prompt_msg <- paste0(linter_details$q_text_configure, " [Default: ", linter_details$default_value, " if blank, or type 'disable']: ")
      user_choice_raw <- readline(prompt = prompt_msg) # Reads from console
      temporary_user_inputs_for_log <<- c(temporary_user_inputs_for_log, paste0("User input for '", linter_name, " config': '", user_choice_raw, "'"))

      if (tolower(user_choice_raw) == "disable") {
        linter_setting_for_file <- paste0(linter_name, " = NULL")
        action_taken_msg_console <- paste0("  -> ", linter_name, " will be DISABLED.\n")
      } else if (user_choice_raw == "") {
        linter_setting_for_file <- paste0(linter_name, "(", linter_details$default_value, ")")
        action_taken_msg_console <- paste0("  -> ", linter_name, " set to default: ", linter_details$default_value, ".\n")
      } else if (grepl("^[0-9]+$", user_choice_raw)) {
        val <- as.integer(user_choice_raw)
        if (val >= 10) {
          linter_setting_for_file <- paste0(linter_name, "(", val, ")")
          action_taken_msg_console <- paste0("  -> ", linter_name, " set to: ", val, ".\n")
        } else {
          linter_setting_for_file <- paste0(linter_name, "(", linter_details$default_value, ")")
          action_taken_msg_console <- paste0("  -> Invalid length '", user_choice_raw, "'. ", linter_name, " set to default: ", linter_details$default_value, ".\n")
        }
      } else {
        linter_setting_for_file <- paste0(linter_name, "(", linter_details$default_value, ")")
        action_taken_msg_console <- paste0("  -> Unrecognized input '", user_choice_raw, "'. ", linter_name, " set to default: ", linter_details$default_value, ".\n")
      }
    } else if (linter_details$type == "yes_no_disable") {
      prompt_msg <- paste0(linter_details$q_text, " (y/n) [Default: '", linter_details$default_string_if_active, "' if blank]: ")
      user_choice_raw <- readline(prompt = prompt_msg) # Reads from console
      temporary_user_inputs_for_log <<- c(temporary_user_inputs_for_log, paste0("User input for '", linter_name, " disable choice': '", user_choice_raw, "'"))


      choice_lower <- tolower(substr(user_choice_raw, 1, 1))
      if (choice_lower == "y") {
        linter_setting_for_file <- paste0(linter_name, " = NULL")
        action_taken_msg_console <- paste0("  -> ", linter_name, " will be DISABLED.\n")
      } else if (choice_lower == "n" || user_choice_raw == "") {
        action_taken_msg_console <- paste0("  -> ", linter_name, " will use its default behavior (ACTIVE).\n")
      } else {
        action_taken_msg_console <- paste0("  -> Unrecognized input '", user_choice_raw, "'. ", linter_name, " will use its default behavior (ACTIVE).\n")
      }
    }
    cat(action_taken_msg_console) # To console
    if (!is.na(linter_setting_for_file)) {
      configured_linter_settings <- c(configured_linter_settings, linter_setting_for_file)
    }
  }

  # --- SINK MANAGEMENT: Restore sinks to log file (if they were active before) ---
  if (current_log_output_to_file_arg && !is.null(current_log_con_arg)) {
    if (output_was_sinking_before_interactive) sink(current_log_con_arg, type = "output", append = TRUE)
    if (message_was_sinking_before_interactive) sink(current_log_con_arg, type = "message", append = TRUE)
  }
  # From here, all cat() and print() should go back to the log file if configured

  # Now write the temporarily stored user inputs to the log
  if (current_log_output_to_file_arg && !is.null(current_log_con_arg) && length(temporary_user_inputs_for_log) > 0) {
    for (log_entry in temporary_user_inputs_for_log) {
      cat(log_entry, "\n", file = current_log_con_arg, append = TRUE)
    }
  } else if (!current_log_output_to_file_arg && length(temporary_user_inputs_for_log) > 0) {
    # If logging to file was off, but we still want to see these (e.g. for debugging)
    cat("--- Captured User Inputs (would be in log if enabled) ---\n")
    for (log_entry in temporary_user_inputs_for_log) {
      cat(log_entry, "\n")
    }
    cat("----------------------------------------------------------\n")
  }


  lintr_config_content_to_write <- ""
  if (length(configured_linter_settings) > 0) {
    lintr_options_string <- paste0("  ", configured_linter_settings, collapse = ",\n")
    lintr_config_content_to_write <- paste0(
      "linters: linters_with_defaults(\n",
      lintr_options_string,
      "\n)" # The closing parenthesis for linters_with_defaults() itself is NOT indented.
      # Individual options *within* the parentheses are indented.
    )
  } else {
    lintr_config_content_to_write <- "linters: linters_with_defaults()"
    cat("No specific linter customizations made that override defaults or disable linters; using basic `linters_with_defaults()`.\n")
  }

  cat("\nGenerated .lintr content will be:\n---\n", lintr_config_content_to_write, "\n---\n") # Goes to log

  tryCatch(
    {
      writeLines(lintr_config_content_to_write, dot_lintr_file)
      cat(paste0("Successfully created/updated '", dot_lintr_file, "'.\n"))
    },
    error = function(e) {
      cat(paste0("Error: Could not write to '", dot_lintr_file, "'. Please check permissions.\nError details: ", e$message, "\n"))
      cat("Proceeding without writing/updating .lintr file. Default lintr behavior will apply.\n")
    }
  )
  cat("---------------------------------------------\n\n")
}


# Main logic for .lintr configuration
if (file.exists(dot_lintr_file)) {
  # Temporarily ensure console interaction
  output_sinking_for_overwrite_q <- LOG_OUTPUT_TO_FILE && !is.null(log_con) && sink.number(type = "output") > 0
  message_sinking_for_overwrite_q <- LOG_OUTPUT_TO_FILE && !is.null(log_con) && sink.number(type = "message") > 0

  if (output_sinking_for_overwrite_q) sink(type = "output")
  if (message_sinking_for_overwrite_q) sink(type = "message")

  cat("A '.lintr' file already exists.\n")
  overwrite_choice <- readline(prompt = "Do you want to reconfigure linters (potentially overwriting/updating it)? (yes/no/keep) [Default: keep]: ")

  # Resume logging (if they were active)
  if (output_sinking_for_overwrite_q && !is.null(log_con)) sink(log_con, type = "output", append = TRUE)
  if (message_sinking_for_overwrite_q && !is.null(log_con)) sink(log_con, type = "message", append = TRUE)

  # Log the choice immediately after resuming sinks (or to console if LOG_OUTPUT_TO_FILE is false)
  overwrite_log_msg <- paste0("User choice for existing .lintr: '", overwrite_choice, "'\n")
  if (LOG_OUTPUT_TO_FILE && !is.null(log_con)) {
    cat(overwrite_log_msg, file = log_con, append = TRUE)
  } else if (!LOG_OUTPUT_TO_FILE) {
    cat(overwrite_log_msg)
  }


  if (tolower(substr(overwrite_choice, 1, 1)) == "y") {
    cat("Reconfiguring .lintr as requested...\n") # This goes to log if active
    configure_lintr_interactive(log_con, LOG_OUTPUT_TO_FILE)
  } else {
    cat("Keeping existing '.lintr' file. Lintr will use its current settings.\n") # To log
    cat("---------------------------------------------\n\n")
  }
} else {
  cat("No '.lintr' file found. Proceeding to configure linters interactively.\n") # To log
  configure_lintr_interactive(log_con, LOG_OUTPUT_TO_FILE)
}


# --- Styling ---
cat("--- Step 1: Styling .R, .r, .Rmd, .rmd files with styler ---\n")
cat("This reREFORMATS CODE LAYOUT (indentation, spacing) and saves changes directly to files.\n")
cat("IMPORTANT: Commit any current work to Git BEFORE this script to easily review styler's changes.\n\n")

tryCatch(
  {
    cat("Styling .R and .r files...\n")
    styler::style_dir(path = ".", filetype = c("R", "r"), recursive = TRUE)
    cat(".R/.r files styling attempt complete.\n\n")

    cat("Styling .Rmd and .rmd files...\n")
    styler::style_dir(path = ".", filetype = c("Rmd", "rmd"), recursive = TRUE)
    cat(".Rmd/.rmd files styling attempt complete.\n\n")
  },
  error = function(e) {
    cat("Error during styling:\n")
    print(e)
    cat("\n")
  }
)
cat("Step 1: Styling complete.\n---------------------------------------------\n\n")

# --- Linting ---
cat("--- Step 2: Linting R/Rmd files with lintr ---\n")
cat("Using settings from '.lintr' (if it exists and is valid) or lintr defaults.\n")
cat("`lintr` IDENTIFIES issues; it does NOT change your files. Issues appear in the log below and in your editor.\n\n")

r_files <- list.files(path = ".", pattern = "\\.R$", recursive = TRUE, full.names = TRUE, ignore.case = TRUE)
rmd_files <- list.files(path = ".", pattern = "\\.Rmd$", recursive = TRUE, full.names = TRUE, ignore.case = TRUE)
all_files_to_lint <- c(r_files, rmd_files)

self_script_name_pattern <- "configure_style_lint_repo.R"
normalized_self_path <- normalizePath(file.path(".", self_script_name_pattern), mustWork = FALSE) # Make sure it's specific enough
all_files_to_lint <- all_files_to_lint[normalizePath(all_files_to_lint, mustWork = FALSE) != normalized_self_path]


if (length(all_files_to_lint) > 0) {
  cat(paste0("Found ", length(all_files_to_lint), " R/Rmd files to lint (excluding this script itself).\n\n"))
  lint_results_summary <- list()
  total_lint_issues <- 0

  for (file_path in all_files_to_lint) {
    cat("Linting file:", file_path, "...\n")
    tryCatch(
      {
        current_lints <- lintr::lint(file_path)
        num_current_issues <- length(current_lints)
        total_lint_issues <- total_lint_issues + num_current_issues
        lint_results_summary[[basename(file_path)]] <- num_current_issues

        if (num_current_issues > 0) {
          cat("Found", num_current_issues, "linting issue(s) in", basename(file_path), ".\n")
          cat(paste("  (See full details for these", num_current_issues, "issues below or in the log file if LOG_OUTPUT_TO_FILE is TRUE)\n"))
          print(current_lints)
        } else {
          cat("No linting issues found in", basename(file_path), ".\n")
        }
        cat("---\n")
      },
      error = function(e) {
        cat("Error linting file:", file_path, "\n")
        print(e)
        lint_results_summary[[basename(file_path)]] <- "Error during linting"
        cat("---\n")
      }
    )
  }

  cat("\n--- Linting Summary ---\n")
  if (length(lint_results_summary) > 0) {
    for (fname in names(lint_results_summary)) {
      cat(paste0("File: ", fname, " - Issues/Status: ", lint_results_summary[[fname]], "\n"))
    }
  }
  cat(paste0("Total linting issues found across all files (excluding errors): ", total_lint_issues, "\n"))
} else {
  cat("No .R or .Rmd files found to lint (excluding this script itself).\n")
}
cat("---------------------------------------------\n\n")
cat("--- Repository Processing Finished ---\n")

# Restore console output and close log file connection
if (LOG_OUTPUT_TO_FILE && !is.null(log_con)) {
  if (sink.number(type = "message") > 0) sink(type = "message")
  if (sink.number(type = "output") > 0) sink(type = "output")
  if (inherits(log_con, "connection") && isOpen(log_con)) {
    close(log_con)
  }
  base::cat(paste0("\nFull output logged to: ", file.path(getwd(), actual_log_file_name), "\n"))
  base::cat("Review the log file for detailed messages and linting reports.\n")
  base::cat("Remember: 'styler' reformats layout. 'lintr' finds other issues; use the .lintr config or fix code manually.\n")
} else if (!LOG_OUTPUT_TO_FILE) {
  base::cat("\nLogging to file was disabled. All output was to console.\n")
  base::cat("Remember: 'styler' reformats layout. 'lintr' finds other issues; use the .lintr config or fix code manually.\n")
}
