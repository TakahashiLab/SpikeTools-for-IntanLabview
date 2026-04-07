args <- commandArgs(trailingOnly = TRUE)

default_csv <- "J:/OneDrive - 同志社大学/共有データ(山中)/salmon/salmon2025-2/salmon2025-2_result/speed_heading_data2025_2_trout.csv"
fps <- 30
analysis_start_sec <- 15 * 60

csv_path <- if (length(args) >= 1) args[[1]] else default_csv
out_dir <- if (length(args) >= 2) args[[2]] else dirname(csv_path)

has_non_ascii <- function(x) {
  grepl("[^\x01-\x7F]", x)
}

normalize_windows_path <- function(path) {
  path <- gsub("/", "\\\\", path)
  normalizePath(path, winslash = "\\", mustWork = FALSE)
}

powershell_copy <- function(src, dst, overwrite = TRUE) {
  src <- enc2native(normalize_windows_path(src))
  dst <- enc2native(normalize_windows_path(dst))
  dst_dir <- dirname(dst)
  dir.create(dst_dir, recursive = TRUE, showWarnings = FALSE)
  overwrite_flag <- if (overwrite) " -Force" else ""
  cmd <- sprintf(
    "Copy-Item -LiteralPath '%s' -Destination '%s'%s",
    gsub("'", "''", src),
    gsub("'", "''", dst),
    overwrite_flag
  )
  status <- system2(
    "powershell",
    c("-NoProfile", "-ExecutionPolicy", "Bypass", "-Command", cmd),
    stdout = FALSE,
    stderr = FALSE
  )
  if (!identical(status, 0L)) {
    stop("Failed to copy file via PowerShell: ", src, " -> ", dst)
  }
}

stage_input_for_r <- function(path) {
  if (!has_non_ascii(path)) {
    return(path)
  }
  ext <- tools::file_ext(path)
  staged <- file.path(
    tempdir(),
    paste0("salmonBehav2025_input", if (nzchar(ext)) paste0(".", ext) else "")
  )
  powershell_copy(path, staged, overwrite = TRUE)
  staged
}

prepare_output_dir_for_r <- function(path) {
  if (!has_non_ascii(path)) {
    dir.create(path, recursive = TRUE, showWarnings = FALSE)
    return(list(write_dir = path, final_dir = path, staged = FALSE))
  }
  write_dir <- file.path(tempdir(), "salmonBehav2025_output")
  dir.create(write_dir, recursive = TRUE, showWarnings = FALSE)
  list(write_dir = write_dir, final_dir = path, staged = TRUE)
}

csv_path_r <- stage_input_for_r(csv_path)
out_info <- prepare_output_dir_for_r(out_dir)
out_dir_r <- out_info$write_dir

message("Reading CSV: ", csv_path)
data <- read.csv(csv_path_r, stringsAsFactors = FALSE)

required_cols <- c("Condition", "Cohort", "Fish", "Time", "Speed", "Heading")
missing_cols <- setdiff(required_cols, names(data))
if (length(missing_cols) > 0) {
  stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
}

if ("Level" %in% names(data)) {
  data <- data[data$Level == "sample", , drop = FALSE]
}

if ("PairKey_Fish" %in% names(data)) {
  data$PairKey_Fish[data$PairKey_Fish == ""] <- NA
} else {
  data$PairKey_Fish <- paste0("C", data$Cohort, "_F", data$Fish)
}

data <- data[is.finite(data$Time), , drop = FALSE]
data <- data[data$Time >= analysis_start_sec, , drop = FALSE]

if (nrow(data) == 0) {
  stop("No sample rows remained after filtering Time >= ", analysis_start_sec, " sec.")
}

circular_mean_deg <- function(theta_deg) {
  theta_deg <- theta_deg[is.finite(theta_deg)]
  if (length(theta_deg) == 0) {
    return(NA_real_)
  }
  theta_rad <- theta_deg * pi / 180
  c_mean <- mean(cos(theta_rad))
  s_mean <- mean(sin(theta_rad))
  atan2(s_mean, c_mean) * 180 / pi
}

finite_mean <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) == 0) {
    return(NA_real_)
  }
  mean(x)
}

finite_var <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) < 2) {
    return(NA_real_)
  }
  stats::var(x)
}

finite_diff_mean <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) < 2) {
    return(NA_real_)
  }
  mean(diff(x))
}

heading_revolutions_mean <- function(theta_deg) {
  theta_deg <- theta_deg[is.finite(theta_deg)]
  if (length(theta_deg) < 2) {
    return(NA_real_)
  }
  mean(cumsum(c(0, diff(theta_deg))))
}

momentary_activity_index <- function(speed) {
  speed <- speed[is.finite(speed)]
  if (length(speed) < 2) {
    return(NA_real_)
  }
  speed_sd <- stats::sd(speed)
  if (!is.finite(speed_sd)) {
    return(NA_real_)
  }
  threshold <- 3 * speed_sd
  mean(speed > threshold)
}

mean_distance_10s <- function(speed, fps) {
  speed <- speed[is.finite(speed)]
  if (length(speed) == 0) {
    return(NA_real_)
  }
  chunk_size <- 10 * fps
  starts <- seq.int(1, length(speed), by = chunk_size)
  distances <- numeric(length(starts))
  for (i in seq_along(starts)) {
    idx <- starts[i]:min(starts[i] + chunk_size - 1, length(speed))
    distances[i] <- sum(speed[idx]) / fps
  }
  mean(distances)
}

mean_speed_per_minute <- function(speed, fps) {
  speed <- speed[is.finite(speed)]
  if (length(speed) == 0) {
    return(NA_real_)
  }
  chunk_size <- 60 * fps
  starts <- seq.int(1, length(speed), by = chunk_size)
  minute_means <- numeric(length(starts))
  for (i in seq_along(starts)) {
    idx <- starts[i]:min(starts[i] + chunk_size - 1, length(speed))
    minute_means[i] <- mean(speed[idx])
  }
  mean(minute_means)
}

summarize_group <- function(df) {
  df <- df[order(df$Time), , drop = FALSE]
  speed <- df$Speed
  heading <- df$Heading
  valid_speed <- is.finite(speed)
  valid_heading <- is.finite(heading)
  speed_for_heading <- speed[is.finite(speed)]
  speed_threshold <- if (length(speed_for_heading) > 1) 3 * stats::sd(speed_for_heading) else Inf
  heading_keep <- valid_heading & valid_speed & speed <= speed_threshold

  data.frame(
    Session = if ("Session" %in% names(df)) df$Session[1] else NA_character_,
    Condition = df$Condition[1],
    Cohort = df$Cohort[1],
    Fish = df$Fish[1],
    PairKey_Fish = df$PairKey_Fish[1],
    MeanSpeed = finite_mean(speed),
    MeanAcceleration = finite_diff_mean(speed),
    ActivityFrequency = momentary_activity_index(speed),
    MeanHeading_deg = circular_mean_deg(heading[heading_keep]),
    HeadingRevolutions = heading_revolutions_mean(heading[heading_keep]),
    Distance_10s = mean_distance_10s(speed, fps),
    SpeedVariance = finite_var(speed),
    MeanSpeedPerMinute = mean_speed_per_minute(speed, fps),
    N = nrow(df),
    NValidSpeed = sum(valid_speed),
    NValidHeading = sum(valid_heading),
    stringsAsFactors = FALSE
  )
}

group_keys <- paste(data$Condition, data$Cohort, data$Fish, data$PairKey_Fish, sep = "||")
grouped <- split(data, group_keys, drop = TRUE)
averaged_data <- do.call(rbind, lapply(grouped, summarize_group))
rownames(averaged_data) <- NULL

safe_paired_t <- function(df, control_label, other_label, metric) {
  sub <- df[df$Condition %in% c(control_label, other_label), c("PairKey_Fish", "Condition", metric), drop = FALSE]
  sub <- sub[is.finite(sub[[metric]]), , drop = FALSE]
  if (nrow(sub) == 0) {
    return(data.frame(CondA = control_label, CondB = other_label, Metric = metric, N = 0, p = NA_real_, t = NA_real_))
  }

  wide <- reshape(sub, idvar = "PairKey_Fish", timevar = "Condition", direction = "wide")
  col_a <- paste0(metric, ".", control_label)
  col_b <- paste0(metric, ".", other_label)
  if (!(col_a %in% names(wide)) || !(col_b %in% names(wide))) {
    return(data.frame(CondA = control_label, CondB = other_label, Metric = metric, N = 0, p = NA_real_, t = NA_real_))
  }

  valid <- is.finite(wide[[col_a]]) & is.finite(wide[[col_b]])
  if (sum(valid) < 2) {
    return(data.frame(CondA = control_label, CondB = other_label, Metric = metric, N = sum(valid), p = NA_real_, t = NA_real_))
  }

  tt <- stats::t.test(wide[[col_a]][valid], wide[[col_b]][valid], paired = TRUE)
  data.frame(
    CondA = control_label,
    CondB = other_label,
    Metric = metric,
    N = sum(valid),
    MeanA = mean(wide[[col_a]][valid]),
    MeanB = mean(wide[[col_b]][valid]),
    p = unname(tt$p.value),
    t = unname(tt$statistic),
    stringsAsFactors = FALSE
  )
}

run_group_tests <- function(df, control_label, pattern, group_label) {
  conditions <- unique(df$Condition)
  candidates <- conditions[grepl(pattern, conditions)]
  candidates <- setdiff(candidates, control_label)
  metrics <- c("MeanSpeed", "Distance_10s")
  out <- list()
  for (cond in candidates) {
    for (metric in metrics) {
      out[[length(out) + 1]] <- safe_paired_t(df, control_label, cond, metric)
    }
  }
  if (length(out) == 0) {
    return(data.frame(Group = character(), CondA = character(), CondB = character(), Metric = character(), N = integer(), MeanA = numeric(), MeanB = numeric(), p = numeric(), t = numeric()))
  }
  res <- do.call(rbind, out)
  res$Group <- group_label
  res[, c("Group", "CondA", "CondB", "Metric", "N", "MeanA", "MeanB", "p", "t")]
}

sea_tests <- run_group_tests(averaged_data, "controlSea1", "Sea", "Sea")
rotation_tests <- run_group_tests(averaged_data, "controlT1", "rotation|controlT", "Rotation")
microt_tests <- run_group_tests(averaged_data, "controlT1", "microT|controlT", "microT")
pairwise_tests <- rbind(sea_tests, rotation_tests, microt_tests)

summary_csv <- file.path(out_dir_r, "salmonBehav2025_summary.csv")
pairwise_csv <- file.path(out_dir_r, "salmonBehav2025_pairwise_tests.csv")

write.csv(averaged_data, summary_csv, row.names = FALSE)
write.csv(pairwise_tests, pairwise_csv, row.names = FALSE)

if (isTRUE(out_info$staged)) {
  final_summary_csv <- file.path(out_info$final_dir, "salmonBehav2025_summary.csv")
  final_pairwise_csv <- file.path(out_info$final_dir, "salmonBehav2025_pairwise_tests.csv")
  powershell_copy(summary_csv, final_summary_csv, overwrite = TRUE)
  powershell_copy(pairwise_csv, final_pairwise_csv, overwrite = TRUE)
  summary_csv <- final_summary_csv
  pairwise_csv <- final_pairwise_csv
}

message("Rows read: ", nrow(data))
message("Unique conditions: ", paste(sort(unique(averaged_data$Condition)), collapse = ", "))
message("Summary written: ", summary_csv)
message("Pairwise tests written: ", pairwise_csv)

print(utils::head(averaged_data))
print(pairwise_tests)
