#!/usr/bin/env Rscript

# Jitter analysis driver for the snow crab ADMB model.
# Runs independent jittered starts in parallel, summarizes convergence, and
# plots final parameter values and reported time series.

parse_args <- function(args = commandArgs(trailingOnly = TRUE)) {
  out <- list(
    n = 20L,
    cores = max(1L, parallel::detectCores() - 1L),
    jitter_sd = 0.1,
    seed = 123,
    gradient_cutoff = 0.001,
    require_hessian = TRUE,
    include_failed_in_plots = FALSE,
    overwrite = FALSE,
    model_dir = normalizePath(".", winslash = "/", mustWork = FALSE),
    output_dir = NULL
  )

  for (arg in args) {
    if (!grepl("^--", arg)) {
      stop("Arguments must use --name=value format. Problem argument: ", arg)
    }
    parts <- strsplit(sub("^--", "", arg), "=", fixed = TRUE)[[1]]
    key <- gsub("-", "_", parts[1])
    value <- if (length(parts) > 1) paste(parts[-1], collapse = "=") else "TRUE"

    if (!key %in% names(out)) {
      stop("Unknown argument: --", parts[1])
    }

    if (is.logical(out[[key]])) {
      out[[key]] <- tolower(value) %in% c("true", "t", "1", "yes", "y")
    } else if (is.integer(out[[key]])) {
      out[[key]] <- as.integer(value)
    } else if (is.numeric(out[[key]])) {
      out[[key]] <- as.numeric(value)
    } else {
      out[[key]] <- value
    }
  }

  if (is.null(out$output_dir)) {
    out$output_dir <- file.path(out$model_dir, "jitter_runs")
  }
  out$model_dir <- normalizePath(out$model_dir, winslash = "/", mustWork = TRUE)
  out$output_dir <- normalizePath(out$output_dir, winslash = "/", mustWork = FALSE)
  out
}

num_scan <- function(x) {
  if (length(x) == 0) return(numeric(0))
  suppressWarnings(scan(text = paste(x, collapse = " "), quiet = TRUE))
}

extract_numbers <- function(x) {
  hits <- regmatches(
    x,
    gregexpr("[-+]?[0-9]*\\.?[0-9]+([eE][-+]?[0-9]+)?", x, perl = TRUE)
  )[[1]]
  as.numeric(hits)
}

parse_par <- function(path) {
  lines <- readLines(path, warn = FALSE)
  header <- lines[grepl("^# Number of parameters", lines)]
  meta <- list(n_parameters = NA_integer_, objective = NA_real_, max_gradient = NA_real_)
  if (length(header)) {
    nums <- extract_numbers(header[1])
    if (length(nums) >= 3) {
      meta$n_parameters <- as.integer(nums[1])
      meta$objective <- nums[2]
      meta$max_gradient <- nums[3]
    }
  }

  block_start <- grep("^# [A-Za-z0-9_]+:", lines)
  blocks <- list()
  if (length(block_start)) {
    for (i in seq_along(block_start)) {
      start <- block_start[i]
      end <- if (i < length(block_start)) block_start[i + 1] - 1L else length(lines)
      name <- sub("^# ([A-Za-z0-9_]+):.*$", "\\1", lines[start])
      values <- num_scan(lines[(start + 1L):end])
      blocks[[name]] <- values
    }
  }
  list(meta = meta, blocks = blocks)
}

write_pin <- function(blocks, path, digits = 12) {
  con <- file(path, open = "wt")
  on.exit(close(con), add = TRUE)
  for (nm in names(blocks)) {
    writeLines(paste0("# ", nm, ":"), con)
    vals <- formatC(blocks[[nm]], digits = digits, format = "fg", flag = "#")
    if (length(vals) == 1L) {
      writeLines(vals, con)
    } else {
      chunks <- split(vals, ceiling(seq_along(vals) / 8))
      for (chunk in chunks) writeLines(paste(chunk, collapse = " "), con)
    }
  }
}

parse_std_names <- function(path) {
  if (!file.exists(path)) return(character(0))
  lines <- readLines(path, warn = FALSE)
  lines <- lines[-1]
  fields <- strsplit(trimws(lines), "\\s+")
  names <- vapply(fields, function(x) if (length(x) >= 2) x[2] else NA_character_, character(1))
  unique(stats::na.omit(names))
}

parse_tpl_bounds <- function(path) {
  lines <- readLines(path, warn = FALSE)
  lines <- sub("//.*$", "", lines)
  lines <- trimws(lines)
  init_lines <- lines[grepl("^init_bounded_", lines)]
  bounds <- data.frame(
    name = character(), lower = numeric(), upper = numeric(),
    stringsAsFactors = FALSE
  )

  for (line in init_lines) {
    pieces <- strsplit(line, "\\s+", perl = TRUE)[[1]]
    if (length(pieces) < 2) next
    name <- sub("\\(.*$", "", pieces[2])
    inside <- regmatches(line, gregexpr("\\([^()]*\\)", line))[[1]]
    if (!length(inside)) next
    args <- trimws(strsplit(gsub("[()]", "", inside[length(inside)]), ",", fixed = TRUE)[[1]])
    if (length(args) == 2L) {
      lower <- suppressWarnings(as.numeric(args[1]))
      upper <- suppressWarnings(as.numeric(args[2]))
    } else if (length(args) == 4L) {
      lower <- suppressWarnings(as.numeric(args[3]))
      upper <- suppressWarnings(as.numeric(args[4]))
    } else if (length(args) >= 5L) {
      lower <- suppressWarnings(as.numeric(args[length(args) - 2L]))
      upper <- suppressWarnings(as.numeric(args[length(args) - 1L]))
    } else {
      lower <- NA_real_
      upper <- NA_real_
    }
    if (!is.na(lower) && !is.na(upper)) {
      bounds <- rbind(bounds, data.frame(
        name = name,
        lower = lower,
        upper = upper,
        stringsAsFactors = FALSE
      ))
    }
  }
  bounds[!duplicated(bounds$name), ]
}

clamp <- function(x, lower, upper, eps = 1e-8) {
  if (is.finite(lower)) x <- pmax(x, lower + eps * max(1, abs(lower)))
  if (is.finite(upper)) x <- pmin(x, upper - eps * max(1, abs(upper)))
  x
}

jitter_values <- function(values, lower = -Inf, upper = Inf, jitter_sd = 0.1) {
  values <- as.numeric(values)
  if (!length(values)) return(values)

  if (is.finite(lower) && is.finite(upper) && lower >= 0 && upper <= 1) {
    eps <- 1e-6
    p <- (values - lower) / (upper - lower)
    p <- pmin(pmax(p, eps), 1 - eps)
    z <- stats::qlogis(p) + stats::rnorm(length(values), 0, jitter_sd)
    out <- lower + (upper - lower) * stats::plogis(z)
  } else if (is.finite(lower) && lower >= 0 && all(values > 0)) {
    out <- values * exp(stats::rnorm(length(values), 0, jitter_sd))
  } else if (is.finite(lower) && is.finite(upper)) {
    scale <- pmax(abs(values), (upper - lower) / 4, 1e-6)
    out <- values + stats::rnorm(length(values), 0, jitter_sd * scale)
  } else {
    scale <- pmax(abs(values), 1)
    out <- values + stats::rnorm(length(values), 0, jitter_sd * scale)
  }

  clamp(out, lower, upper)
}

make_jittered_blocks <- function(base_blocks, active_names, bounds, jitter_sd) {
  out <- base_blocks
  jitter_names <- intersect(names(base_blocks), active_names)
  for (nm in jitter_names) {
    row <- bounds[bounds$name == nm, , drop = FALSE]
    lower <- if (nrow(row)) row$lower[1] else -Inf
    upper <- if (nrow(row)) row$upper[1] else Inf
    out[[nm]] <- jitter_values(base_blocks[[nm]], lower, upper, jitter_sd)
  }
  out
}

copy_model_inputs <- function(model_dir, run_dir) {
  needed <- c("snow_down.exe", "snow_down.DAT", "catch_dat.DAT")
  optional <- c("snow_down.tpl")
  files <- file.path(model_dir, c(needed, optional))
  missing <- needed[!file.exists(file.path(model_dir, needed))]
  if (length(missing)) stop("Missing required model files: ", paste(missing, collapse = ", "))
  dir.create(run_dir, recursive = TRUE, showWarnings = FALSE)
  file.copy(files[file.exists(files)], run_dir, overwrite = TRUE)
}

find_admb_output <- function(run_dir, ext) {
  expected <- file.path(run_dir, paste0("snow_down.", ext))
  if (file.exists(expected)) return(expected)
  hits <- list.files(
    run_dir,
    pattern = paste0("\\.", ext, "$"),
    ignore.case = TRUE,
    full.names = TRUE
  )
  if (length(hits)) hits[1] else expected
}

run_one <- function(i, config, base_blocks, active_names, bounds) {
  run_id <- sprintf("run_%03d", i)
  run_dir <- file.path(config$output_dir, run_id)
  if (dir.exists(run_dir) && !config$overwrite) {
    stop("Run directory exists. Use --overwrite=true to replace: ", run_dir)
  }
  if (dir.exists(run_dir) && config$overwrite) {
    unlink(run_dir, recursive = TRUE, force = TRUE)
  }
  dir.create(run_dir, recursive = TRUE, showWarnings = FALSE)
  copy_model_inputs(config$model_dir, run_dir)

  set.seed(config$seed + i)
  pin_blocks <- make_jittered_blocks(base_blocks, active_names, bounds, config$jitter_sd)
  write_pin(pin_blocks, file.path(run_dir, "snow_down.pin"))
  write.csv(
    data.frame(
      parameter = rep(names(pin_blocks), lengths(pin_blocks)),
      index = unlist(lapply(pin_blocks, seq_along)),
      start_value = unlist(pin_blocks, use.names = FALSE)
    ),
    file.path(run_dir, "jittered_start_values.csv"),
    row.names = FALSE
  )

  oldwd <- getwd()
  on.exit(setwd(oldwd), add = TRUE)
  setwd(run_dir)
  start_time <- Sys.time()
  status <- tryCatch(
    system2(file.path(".", "snow_down.exe"), stdout = "stdout.log", stderr = "stderr.log"),
    error = function(e) {
      writeLines(conditionMessage(e), "system_error.log")
      999L
    }
  )
  elapsed_seconds <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))

  par_path <- find_admb_output(run_dir, "par")
  rep_path <- find_admb_output(run_dir, "rep")
  par <- if (file.exists(par_path)) parse_par(par_path) else list(meta = list(
    n_parameters = NA_integer_, objective = NA_real_, max_gradient = NA_real_
  ), blocks = list())

  hessian_files <- c(
    list.files(run_dir, pattern = "\\.cor$", ignore.case = TRUE, full.names = TRUE),
    list.files(run_dir, pattern = "\\.hes$", ignore.case = TRUE, full.names = TRUE)
  )
  has_hessian <- length(hessian_files) > 0 &&
    any(file.info(hessian_files)$size > 0, na.rm = TRUE)

  data.frame(
    run_id = run_id,
    run_dir = normalizePath(run_dir, winslash = "/", mustWork = FALSE),
    exit_status = as.integer(status),
    elapsed_seconds = elapsed_seconds,
    objective = par$meta$objective,
    max_gradient = par$meta$max_gradient,
    has_hessian = has_hessian,
    par_exists = file.exists(par_path),
    rep_exists = file.exists(rep_path),
    stringsAsFactors = FALSE
  )
}

section_lines <- function(path) {
  lines <- readLines(path, warn = FALSE)
  idx <- grep("^\\$", lines)
  out <- list()
  if (!length(idx)) return(out)
  for (i in seq_along(idx)) {
    start <- idx[i]
    end <- if (i < length(idx)) idx[i + 1L] - 1L else length(lines)
    nm <- sub("^\\$", "", trimws(lines[start]))
    out[[nm]] <- lines[(start + 1L):end]
  }
  out
}

section_vector <- function(sections, name) {
  if (!name %in% names(sections)) return(numeric(0))
  num_scan(sections[[name]])
}

section_matrix <- function(sections, name) {
  if (!name %in% names(sections)) return(matrix(numeric(0), nrow = 0))
  rows <- lapply(sections[[name]], num_scan)
  rows <- rows[lengths(rows) > 0]
  if (!length(rows)) return(matrix(numeric(0), nrow = 0))
  width <- max(lengths(rows))
  mat <- matrix(NA_real_, nrow = length(rows), ncol = width)
  for (i in seq_along(rows)) mat[i, seq_along(rows[[i]])] <- rows[[i]]
  mat
}

years_from_sections <- function(sections) {
  styr <- section_vector(sections, "styr")[1]
  endyr <- section_vector(sections, "endyr")[1]
  if (is.na(styr) || is.na(endyr)) return(integer(0))
  seq.int(styr, endyr)
}

parse_run_outputs <- function(summary) {
  pars <- list()
  ts <- list()
  sel <- list()

  for (r in seq_len(nrow(summary))) {
    run_id <- summary$run_id[r]
    run_dir <- summary$run_dir[r]

    par_path <- find_admb_output(run_dir, "par")
    if (file.exists(par_path)) {
      par <- parse_par(par_path)
      scalar_names <- c(
        "log_avg_rec", "log_m_mu", "prop_rec", "sigma_m", "sigma_q", "log_f",
        "fish_ret_sel_50", "fish_ret_sel_50_post",
        "fish_ret_sel_slope", "fish_ret_sel_slope_post",
        "fish_tot_sel_offset", "fish_tot_sel_slope",
        "surv_omega", "surv_alpha1", "surv_beta1", "surv_alpha2", "surv_beta2"
      )
      for (nm in intersect(scalar_names, names(par$blocks))) {
        vals <- par$blocks[[nm]]
        pars[[length(pars) + 1L]] <- data.frame(
          run_id = run_id,
          parameter = if (length(vals) == 1L) nm else paste0(nm, "[", seq_along(vals), "]"),
          value = vals,
          objective = summary$objective[r],
          max_gradient = summary$max_gradient[r],
          accepted = summary$accepted[r],
          stringsAsFactors = FALSE
        )
      }
    }

    rep_path <- find_admb_output(run_dir, "rep")
    if (!file.exists(rep_path)) next
    sections <- section_lines(rep_path)
    years <- years_from_sections(sections)
    sizes <- section_vector(sections, "sizes")

    add_vector_ts <- function(label, values) {
      if (!length(values) || !length(years)) return(NULL)
      n <- min(length(values), length(years))
      data.frame(
        run_id = run_id,
        year = years[seq_len(n)],
        quantity = label,
        value = values[seq_len(n)],
        accepted = summary$accepted[r],
        stringsAsFactors = FALSE
      )
    }

    imm_pop <- section_matrix(sections, "pred_imm_pop_num")
    mat_pop <- section_matrix(sections, "pred_mat_pop_num")
    nat_m <- section_matrix(sections, "natural mortality")
    mat_m <- section_matrix(sections, "mature natural mortality")

    chunks <- list(
      add_vector_ts("recruitment", section_vector(sections, "recruits")),
      add_vector_ts("fishing_mortality", section_vector(sections, "est_fishing_mort")),
      add_vector_ts("survey_immature_abundance", section_vector(sections, "imm_numbers_pred")),
      add_vector_ts("survey_mature_abundance", section_vector(sections, "mat_numbers_pred"))
    )

    if (nrow(imm_pop) && nrow(mat_pop) && length(years)) {
      n <- min(nrow(imm_pop), nrow(mat_pop), length(years))
      chunks[[length(chunks) + 1L]] <- data.frame(
        run_id = run_id,
        year = years[seq_len(n)],
        quantity = "total_abundance",
        value = rowSums(imm_pop[seq_len(n), , drop = FALSE], na.rm = TRUE) +
          rowSums(mat_pop[seq_len(n), , drop = FALSE], na.rm = TRUE),
        accepted = summary$accepted[r],
        stringsAsFactors = FALSE
      )
      chunks[[length(chunks) + 1L]] <- data.frame(
        run_id = run_id,
        year = years[seq_len(n)],
        quantity = "immature_population",
        value = rowSums(imm_pop[seq_len(n), , drop = FALSE], na.rm = TRUE),
        accepted = summary$accepted[r],
        stringsAsFactors = FALSE
      )
      chunks[[length(chunks) + 1L]] <- data.frame(
        run_id = run_id,
        year = years[seq_len(n)],
        quantity = "mature_population",
        value = rowSums(mat_pop[seq_len(n), , drop = FALSE], na.rm = TRUE),
        accepted = summary$accepted[r],
        stringsAsFactors = FALSE
      )
    }

    if (nrow(nat_m) && length(years)) {
      n <- min(nrow(nat_m), length(years))
      chunks[[length(chunks) + 1L]] <- data.frame(
        run_id = run_id,
        year = years[seq_len(n)],
        quantity = "natural_mortality_mean",
        value = rowMeans(nat_m[seq_len(n), , drop = FALSE], na.rm = TRUE),
        accepted = summary$accepted[r],
        stringsAsFactors = FALSE
      )
    }

    if (nrow(mat_m) && length(years)) {
      n <- min(nrow(mat_m), length(years))
      chunks[[length(chunks) + 1L]] <- data.frame(
        run_id = run_id,
        year = years[seq_len(n)],
        quantity = "mature_natural_mortality_mean",
        value = rowMeans(mat_m[seq_len(n), , drop = FALSE], na.rm = TRUE),
        accepted = summary$accepted[r],
        stringsAsFactors = FALSE
      )
    }

    chunks <- Filter(Negate(is.null), chunks)
    if (length(chunks)) ts[[length(ts) + 1L]] <- do.call(rbind, chunks)

    selectivity_sections <- c(
      survey_selectivity = "survey selectivity",
      mature_survey_selectivity = "mature survey selectivity",
      retained_fish_selectivity = "ret_fish_sel",
      total_fish_selectivity = "total_fish_sel"
    )
    for (kind in names(selectivity_sections)) {
      mat <- section_matrix(sections, selectivity_sections[[kind]])
      if (!nrow(mat) || !length(years) || !length(sizes)) next
      nyr <- min(nrow(mat), length(years))
      nsz <- min(ncol(mat), length(sizes))
      grid <- expand.grid(
        year = years[seq_len(nyr)],
        size = sizes[seq_len(nsz)]
      )
      sel[[length(sel) + 1L]] <- data.frame(
        run_id = run_id,
        selectivity = kind,
        year = grid$year,
        size = grid$size,
        value = as.vector(t(mat[seq_len(nyr), seq_len(nsz), drop = FALSE])),
        accepted = summary$accepted[r],
        stringsAsFactors = FALSE
      )
    }
  }

  list(
    parameters = if (length(pars)) do.call(rbind, pars) else data.frame(),
    timeseries = if (length(ts)) do.call(rbind, ts) else data.frame(),
    selectivity = if (length(sel)) do.call(rbind, sel) else data.frame()
  )
}

save_plot <- function(plot, path, width = 9, height = 6) {
  ggplot2::ggsave(path, plot, width = width, height = height, units = "in", dpi = 150)
}

plot_outputs <- function(parsed, summary, plot_dir, include_failed_in_plots = FALSE) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required for plotting. Install it with install.packages('ggplot2').")
  }
  dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
  gg <- asNamespace("ggplot2")

  plot_summary <- if (include_failed_in_plots) summary else summary[summary$accepted, , drop = FALSE]
  plot_parameters <- if (include_failed_in_plots) parsed$parameters else parsed$parameters[parsed$parameters$accepted, , drop = FALSE]
  plot_ts <- if (include_failed_in_plots) parsed$timeseries else parsed$timeseries[parsed$timeseries$accepted, , drop = FALSE]
  plot_sel <- if (include_failed_in_plots) parsed$selectivity else parsed$selectivity[parsed$selectivity$accepted, , drop = FALSE]

  if (nrow(plot_summary)) {
    p <- gg$ggplot(plot_summary, gg$aes(x = objective)) +
      gg$geom_histogram(bins = 30, color = "white", fill = "#4C78A8") +
      gg$theme_bw() +
      gg$labs(x = "Objective function value", y = "Number of runs")
    save_plot(p, file.path(plot_dir, "objective_histogram.png"), 7, 5)
  }

  if (nrow(plot_parameters)) {
    p <- gg$ggplot(plot_parameters, gg$aes(x = objective, y = value)) +
      gg$geom_point(alpha = 0.75, size = 1.8, color = "#2F4B7C") +
      gg$facet_wrap(stats::as.formula("~ parameter"), scales = "free_y") +
      gg$theme_bw() +
      gg$labs(x = "Objective function value", y = "Final parameter value")
    save_plot(p, file.path(plot_dir, "important_parameters_vs_objective.png"), 12, 8)
  }

  if (nrow(plot_ts)) {
    for (quantity in unique(plot_ts$quantity)) {
      dat <- plot_ts[plot_ts$quantity == quantity, , drop = FALSE]
      p <- gg$ggplot(dat, gg$aes(x = year, y = value, group = run_id)) +
        gg$geom_line(alpha = 0.35, color = "#3B6EA8") +
        gg$theme_bw() +
        gg$labs(x = "Year", y = quantity)
      save_plot(p, file.path(plot_dir, paste0("timeseries_", quantity, ".png")), 9, 5)
    }
  }

  if (nrow(plot_sel)) {
    mean_sel <- stats::aggregate(
      value ~ run_id + selectivity + year,
      data = plot_sel,
      FUN = mean,
      na.rm = TRUE
    )
    p <- gg$ggplot(mean_sel, gg$aes(x = year, y = value, group = run_id)) +
      gg$geom_line(alpha = 0.35, color = "#4C78A8") +
      gg$facet_wrap(stats::as.formula("~ selectivity"), scales = "free_y") +
      gg$theme_bw() +
      gg$labs(x = "Year", y = "Mean selectivity across sizes")
    save_plot(p, file.path(plot_dir, "selectivity_mean_timeseries.png"), 12, 7)

    years <- sort(unique(plot_sel$year))
    selected_years <- unique(round(stats::quantile(years, probs = c(0, 0.5, 1), names = FALSE)))
    at_size <- plot_sel[plot_sel$year %in% selected_years, , drop = FALSE]
    p <- gg$ggplot(at_size, gg$aes(x = size, y = value, group = run_id)) +
      gg$geom_line(alpha = 0.25, color = "#2F4B7C") +
      gg$facet_grid(stats::as.formula("selectivity ~ year"), scales = "free_y") +
      gg$theme_bw() +
      gg$labs(x = "Size", y = "Selectivity")
    save_plot(p, file.path(plot_dir, "selectivity_at_size_selected_years.png"), 12, 8)
  }
}

run_snow_jitter <- function(config = parse_args(character())) {
  if (!requireNamespace("parallel", quietly = TRUE)) {
    stop("Package 'parallel' is required.")
  }

  dir.create(config$output_dir, recursive = TRUE, showWarnings = FALSE)
  summary_dir <- file.path(config$output_dir, "summary")
  dir.create(summary_dir, recursive = TRUE, showWarnings = FALSE)

  base_par_path <- file.path(config$model_dir, "snow_down.par")
  base_std_path <- file.path(config$model_dir, "snow_down.std")
  tpl_path <- file.path(config$model_dir, "snow_down.tpl")

  base_par <- parse_par(base_par_path)
  active_names <- intersect(parse_std_names(base_std_path), names(base_par$blocks))
  bounds <- parse_tpl_bounds(tpl_path)

  if (!length(active_names)) {
    stop("No active parameters found from snow_down.std that match snow_down.par blocks.")
  }

  writeLines(
    c(
      paste("model_dir:", config$model_dir),
      paste("output_dir:", config$output_dir),
      paste("n:", config$n),
      paste("cores:", config$cores),
      paste("jitter_sd:", config$jitter_sd),
      paste("seed:", config$seed),
      paste("gradient_cutoff:", config$gradient_cutoff),
      paste("require_hessian:", config$require_hessian),
      paste("include_failed_in_plots:", config$include_failed_in_plots),
      paste("active_parameters:", paste(active_names, collapse = ", "))
    ),
    file.path(summary_dir, "jitter_config.txt")
  )

  cores <- min(config$cores, config$n)
  if (.Platform$OS.type == "windows" && cores > 1L) {
    cl <- parallel::makeCluster(cores)
    on.exit(parallel::stopCluster(cl), add = TRUE)
    parallel::clusterExport(
      cl,
      varlist = c(
        "config", "base_par", "active_names", "bounds", "run_one",
        "copy_model_inputs", "make_jittered_blocks", "jitter_values",
        "clamp", "write_pin", "parse_par", "num_scan", "extract_numbers",
        "find_admb_output"
      ),
      envir = environment()
    )
    run_summary <- do.call(rbind, parallel::parLapply(
      cl, seq_len(config$n), run_one,
      config = config,
      base_blocks = base_par$blocks,
      active_names = active_names,
      bounds = bounds
    ))
  } else {
    run_summary <- do.call(rbind, parallel::mclapply(
      seq_len(config$n), run_one,
      config = config,
      base_blocks = base_par$blocks,
      active_names = active_names,
      bounds = bounds,
      mc.cores = cores
    ))
  }

  run_summary$accepted <- with(
    run_summary,
    exit_status == 0 &
      par_exists &
      rep_exists &
      (!config$require_hessian | has_hessian) &
      !is.na(max_gradient) &
      max_gradient <= config$gradient_cutoff
  )

  write.csv(run_summary, file.path(summary_dir, "run_summary.csv"), row.names = FALSE)
  write.csv(
    data.frame(
      n_runs = nrow(run_summary),
      n_accepted = sum(run_summary$accepted),
      mean_elapsed_seconds_all = mean(run_summary$elapsed_seconds, na.rm = TRUE),
      median_elapsed_seconds_all = stats::median(run_summary$elapsed_seconds, na.rm = TRUE),
      mean_elapsed_seconds_accepted = mean(run_summary$elapsed_seconds[run_summary$accepted], na.rm = TRUE),
      median_elapsed_seconds_accepted = stats::median(run_summary$elapsed_seconds[run_summary$accepted], na.rm = TRUE)
    ),
    file.path(summary_dir, "runtime_summary.csv"),
    row.names = FALSE
  )

  parsed <- parse_run_outputs(run_summary)
  write.csv(parsed$parameters, file.path(summary_dir, "parameter_values.csv"), row.names = FALSE)
  write.csv(parsed$timeseries, file.path(summary_dir, "timeseries_values.csv"), row.names = FALSE)
  write.csv(parsed$selectivity, file.path(summary_dir, "selectivity_values.csv"), row.names = FALSE)

  plot_outputs(
    parsed,
    run_summary,
    file.path(summary_dir, "plots"),
    include_failed_in_plots = config$include_failed_in_plots
  )

  invisible(list(config = config, summary = run_summary, parsed = parsed))
}

if (sys.nframe() == 0L) {
  cfg <- parse_args()
  result <- run_snow_jitter(cfg)
  cat("Completed snow crab jitter analysis\n")
  cat("Runs:", nrow(result$summary), "\n")
  cat("Accepted:", sum(result$summary$accepted), "\n")
  cat("Summary:", file.path(cfg$output_dir, "summary"), "\n")
}
