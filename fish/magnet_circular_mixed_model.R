args <- commandArgs(trailingOnly = TRUE)

`%||%` <- function(x, y) if (is.null(x) || length(x) == 0 || identical(x, "")) y else x

csv_path <- args[[1]] %||% "speed_heading_data2025.csv"
out_dir <- args[[2]] %||% dirname(csv_path)
window_start <- if (length(args) >= 3) as.numeric(args[[3]]) else 300
window_end <- if (length(args) >= 4) as.numeric(args[[4]]) else 600
condition_arg <- if (length(args) >= 5) args[[5]] else ""

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
set.seed(1)

message("Reading CSV: ", csv_path)
dat_all <- read.csv(csv_path, stringsAsFactors = FALSE)

required_cols <- c("Level", "Condition", "Cohort", "Fish", "Time", "HeadingMag", "cosMag", "sinMag")
missing_cols <- setdiff(required_cols, names(dat_all))
if (length(missing_cols) > 0) {
  stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
}

derive_order_in_session <- function(dat) {
  if ("OrderInSession" %in% names(dat) && any(is.finite(suppressWarnings(as.numeric(dat$OrderInSession))))) {
    dat$OrderInSession <- suppressWarnings(as.numeric(dat$OrderInSession))
    return(dat)
  }
  if (!("Session" %in% names(dat))) {
    stop("Missing required columns: OrderInSession or Session")
  }
  session_key <- data.frame(
    Cohort = dat$Cohort,
    Session = dat$Session,
    stringsAsFactors = FALSE
  )
  session_key <- unique(session_key)
  session_key <- session_key[!(is.na(session_key$Cohort) | is.na(session_key$Session) | session_key$Session == ""), , drop = FALSE]
  if (nrow(session_key) == 0) {
    stop("Could not derive OrderInSession because Session values are empty.")
  }
  session_key$SessionStamp <- suppressWarnings(as.POSIXct(substr(session_key$Session, 1, 19), format = "%Y-%m-%d %H-%M-%S", tz = "UTC"))
  session_key$SessionIndex <- seq_len(nrow(session_key))
  ord <- order(session_key$Cohort, is.na(session_key$SessionStamp), session_key$SessionStamp, session_key$SessionIndex)
  session_key <- session_key[ord, , drop = FALSE]
  session_key$OrderInSession <- ave(
    session_key$SessionIndex,
    session_key$Cohort,
    FUN = function(x) seq_along(x)
  )
  dat <- merge(dat, session_key[, c("Cohort", "Session", "OrderInSession")], by = c("Cohort", "Session"), all.x = TRUE, sort = FALSE)
  dat$OrderInSession <- suppressWarnings(as.numeric(dat$OrderInSession))
  dat
}

dat_all <- derive_order_in_session(dat_all)

add_distance10s <- function(dat) {
  if (!all(c("Cohort", "Fish", "Session", "Condition", "Time", "Speed") %in% names(dat))) {
    dat$Distance10s <- NA_real_
    return(dat)
  }
  dat$Distance10s <- NA_real_
  key <- paste(dat$Cohort, dat$Fish, dat$Session, dat$Condition, sep = "|")
  split_idx <- split(seq_len(nrow(dat)), key)
  for (idx in split_idx) {
    if (length(idx) < 2) next
    ord <- order(dat$Time[idx])
    ii <- idx[ord]
    t <- suppressWarnings(as.numeric(dat$Time[ii]))
    v <- suppressWarnings(as.numeric(dat$Speed[ii]))
    valid <- is.finite(t) & is.finite(v)
    if (sum(valid) < 2) next
    dt <- stats::median(diff(t[valid]), na.rm = TRUE)
    if (!is.finite(dt) || dt <= 0) next
    win <- max(1L, round(10 / max(dt, .Machine$double.eps)))
    D <- rep(NA_real_, length(ii))
    runs <- split(seq_along(ii), cumsum(c(1, diff(valid) != 0)))
    for (seg in runs) {
      if (!valid[seg[1]]) next
      if (length(seg) < win) next
      step <- v[seg] * dt
      dist10 <- stats::filter(step, rep(1, win), sides = 1)
      D[seg] <- as.numeric(dist10)
    }
    dat$Distance10s[ii] <- D
  }
  dat
}

dat_all <- add_distance10s(dat_all)

finite_mean <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) == 0) return(NA_real_)
  mean(x)
}

format_num <- function(x, digits = 3) {
  ifelse(is.finite(x), formatC(x, digits = digits, format = "fg", flag = "#"), "NA")
}

format_p <- function(p) {
  if (!is.finite(p)) return("NA")
  if (p < 1e-4) return(formatC(p, digits = 2, format = "e"))
  formatC(p, digits = 4, format = "fg", flag = "#")
}

wrap_deg <- function(x) ((x + 180) %% 360) - 180

circular_mean_deg <- function(theta_deg) {
  theta_deg <- theta_deg[is.finite(theta_deg)]
  if (length(theta_deg) == 0) return(NA_real_)
  theta_rad <- theta_deg * pi / 180
  atan2(mean(sin(theta_rad)), mean(cos(theta_rad))) * 180 / pi
}

circular_rbar <- function(theta_deg) {
  theta_deg <- theta_deg[is.finite(theta_deg)]
  if (length(theta_deg) == 0) return(NA_real_)
  theta_rad <- theta_deg * pi / 180
  sqrt(mean(cos(theta_rad))^2 + mean(sin(theta_rad))^2)
}

first_non_missing <- function(x) {
  idx <- which(!(is.na(x) | x == ""))
  if (length(idx) == 0) return(NA)
  x[idx[1]]
}

sign_test_p <- function(x) {
  x <- x[is.finite(x)]
  x <- x[abs(x) > 1e-8]
  if (length(x) == 0) return(NA_real_)
  pos <- sum(x > 0)
  stats::binom.test(pos, length(x), p = 0.5)$p.value
}

bootstrap_circular_mean_ci <- function(x, n_boot = 2000, conf = 0.95) {
  x <- x[is.finite(x)]
  if (length(x) < 2) return(c(lo = NA_real_, hi = NA_real_))
  alpha <- (1 - conf) / 2
  boot <- replicate(n_boot, circular_mean_deg(sample(x, length(x), replace = TRUE)))
  c(
    lo = as.numeric(stats::quantile(boot, alpha, na.rm = TRUE)),
    hi = as.numeric(stats::quantile(boot, 1 - alpha, na.rm = TRUE))
  )
}

shift_permutation_p <- function(x, n_perm = 2000) {
  x <- x[is.finite(x)]
  x <- wrap_deg(x)
  if (length(x) < 2) return(NA_real_)
  obs <- circular_mean_deg(x)
  if (!is.finite(obs)) return(NA_real_)
  perm_stats <- replicate(n_perm, {
    flipped <- x * sample(c(-1, 1), length(x), replace = TRUE)
    abs(circular_mean_deg(flipped))
  })
  mean(perm_stats >= abs(obs), na.rm = TRUE)
}

summarize_shift_distribution <- function(x, reference_condition, condition, cohort_group = NA_character_, n_boot = 2000, n_perm = 2000) {
  x <- wrap_deg(x[is.finite(x)])
  if (length(x) == 0) {
    return(data.frame(
      ReferenceCondition = reference_condition,
      Condition = condition,
      CohortGroup = cohort_group,
      MeanShiftDeg = NA_real_,
      MedianShiftDeg = NA_real_,
      CI_lo = NA_real_,
      CI_hi = NA_real_,
      SignTestP = NA_real_,
      PermutationP = NA_real_,
      PositiveFraction = NA_real_,
      N = 0,
      stringsAsFactors = FALSE
    ))
  }
  ci <- bootstrap_circular_mean_ci(x, n_boot = n_boot)
  data.frame(
    ReferenceCondition = reference_condition,
    Condition = condition,
    CohortGroup = cohort_group,
    MeanShiftDeg = circular_mean_deg(x),
    MedianShiftDeg = stats::median(x, na.rm = TRUE),
    CI_lo = unname(ci["lo"]),
    CI_hi = unname(ci["hi"]),
    SignTestP = sign_test_p(x),
    PermutationP = shift_permutation_p(x, n_perm = n_perm),
    PositiveFraction = mean(x > 0, na.rm = TRUE),
    N = length(x),
    stringsAsFactors = FALSE
  )
}

circular_correlation <- function(theta_deg_a, theta_deg_b, n_perm = 2000) {
  ok <- is.finite(theta_deg_a) & is.finite(theta_deg_b)
  theta_deg_a <- theta_deg_a[ok]
  theta_deg_b <- theta_deg_b[ok]
  if (length(theta_deg_a) < 5) return(list(r = NA_real_, p = NA_real_, n = length(theta_deg_a)))
  a <- theta_deg_a * pi / 180
  b <- theta_deg_b * pi / 180
  a_bar <- atan2(sum(sin(a)), sum(cos(a)))
  b_bar <- atan2(sum(sin(b)), sum(cos(b)))
  num <- sum(sin(a - a_bar) * sin(b - b_bar))
  den <- sqrt(sum(sin(a - a_bar)^2) * sum(sin(b - b_bar)^2))
  r_obs <- if (den > 0) num / den else NA_real_
  if (!is.finite(r_obs)) return(list(r = r_obs, p = NA_real_, n = length(theta_deg_a)))
  perm_r <- rep(NA_real_, n_perm)
  for (i in seq_len(n_perm)) {
    bp <- sample(b, length(b), replace = FALSE)
    b_bar_p <- atan2(sum(sin(bp)), sum(cos(bp)))
    num_p <- sum(sin(a - a_bar) * sin(bp - b_bar_p))
    den_p <- sqrt(sum(sin(a - a_bar)^2) * sum(sin(bp - b_bar_p)^2))
    perm_r[i] <- if (den_p > 0) num_p / den_p else NA_real_
  }
  list(r = r_obs, p = mean(abs(perm_r) >= abs(r_obs), na.rm = TRUE), n = length(theta_deg_a))
}

permutation_correlation_test <- function(x, y, n_perm = 2000, method = "pearson") {
  ok <- is.finite(x) & is.finite(y)
  x <- x[ok]
  y <- y[ok]
  if (length(x) < 5) return(list(r = NA_real_, p = NA_real_, n = length(x)))
  r_obs <- suppressWarnings(stats::cor(x, y, method = method))
  if (!is.finite(r_obs)) return(list(r = NA_real_, p = NA_real_, n = length(x)))
  perm_r <- replicate(n_perm, suppressWarnings(stats::cor(x, sample(y, length(y), replace = FALSE), method = method)))
  list(r = r_obs, p = mean(abs(perm_r) >= abs(r_obs), na.rm = TRUE), n = length(x))
}

extract_table_col <- function(tbl, patterns) {
  for (pat in patterns) {
    hits <- names(tbl)[grepl(pat, names(tbl), ignore.case = TRUE)]
    if (length(hits) > 0) return(hits[1])
  }
  NA_character_
}

prepare_sample_data <- function(dat, start_sec, end_sec, keep_conditions = NULL) {
  out <- dat[dat$Level == "sample", , drop = FALSE]
  out <- out[is.finite(out$Time), , drop = FALSE]
  out <- out[out$Time >= start_sec & out$Time <= end_sec, , drop = FALSE]
  if (!is.null(keep_conditions) && length(keep_conditions) > 0) {
    out <- out[tolower(out$Condition) %in% keep_conditions, , drop = FALSE]
  }
  if (nrow(out) == 0) return(out)
  out$Condition <- tolower(out$Condition)
  out$Pair <- paste0("C", out$Cohort, "_F", out$Fish)
  out
}

choose_activity_metrics <- function(dat) {
  candidate_metrics <- c("Speed", "Distance10s", "TurnRate", "Curvature", "Tortuosity")
  present <- candidate_metrics[candidate_metrics %in% names(dat)]
  present[vapply(present, function(metric) any(is.finite(suppressWarnings(as.numeric(dat[[metric]]))), na.rm = TRUE), logical(1))]
}

aggregate_fish_condition <- function(dat_window, activity_metrics = choose_activity_metrics(dat_window)) {
  if (nrow(dat_window) == 0) return(data.frame())
  agg_key <- interaction(dat_window$Pair, dat_window$Condition, dat_window$OrderInSession, drop = TRUE, lex.order = TRUE)
  split_dat <- split(dat_window, agg_key, drop = TRUE)
  rows <- lapply(split_dat, function(df) {
    base_row <- data.frame(
      Pair = first_non_missing(df$Pair),
      Condition = first_non_missing(df$Condition),
      Cohort = df$Cohort[1],
      AgeGroup = if ("AgeGroup" %in% names(df)) first_non_missing(as.character(df$AgeGroup)) else NA_character_,
      Fish = df$Fish[1],
      OrderInSession = df$OrderInSession[1],
      MeanHeadingDeg = circular_mean_deg(df$HeadingMag),
      cosMag = finite_mean(df$cosMag),
      sinMag = finite_mean(df$sinMag),
      N = nrow(df),
      stringsAsFactors = FALSE
    )
    for (metric in activity_metrics) {
      base_row[[metric]] <- finite_mean(df[[metric]])
    }
    base_row
  })
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

build_prev_condition <- function(df) {
  df$PrevCondition <- "none"
  for (cohort in unique(df$Cohort)) {
    idx <- which(df$Cohort == cohort)
    ord <- order(df$OrderInSession[idx], df$Condition[idx])
    ii <- idx[ord]
    if (length(ii) > 1) df$PrevCondition[ii[-1]] <- df$Condition[ii[-length(ii)]]
  }
  df
}

build_cohort_group <- function(df) {
  cohort_ids <- sort(unique(df$Cohort))
  seq_labels <- character(length(cohort_ids))
  for (i in seq_along(cohort_ids)) {
    cohort <- cohort_ids[i]
    sub <- df[df$Cohort == cohort, , drop = FALSE]
    sub <- sub[order(sub$OrderInSession, sub$Condition), , drop = FALSE]
    sub <- sub[!duplicated(sub$OrderInSession), , drop = FALSE]
    seq_labels[i] <- paste(sub$Condition, collapse = " -> ")
  }
  df$CohortGroup <- unname(setNames(seq_labels, cohort_ids)[as.character(df$Cohort)])
  df
}

finalize_agg <- function(agg) {
  agg <- build_prev_condition(agg)
  agg <- build_cohort_group(agg)
  agg$ConditionFamily <- condition_family(agg$Condition)
  if ("AgeGroup" %in% names(agg)) {
    agg$AgeGroup <- ifelse(is.na(agg$AgeGroup) | agg$AgeGroup == "", NA_character_, as.character(agg$AgeGroup))
    agg$AgeGroup <- factor(agg$AgeGroup)
  }
  agg$Condition <- factor(agg$Condition)
  agg$ConditionFamily <- factor(agg$ConditionFamily)
  agg$PrevCondition <- factor(agg$PrevCondition)
  agg$CohortGroup <- factor(agg$CohortGroup)
  agg$Cohort <- factor(agg$Cohort)
  agg$Fish <- factor(agg$Fish)
  agg
}

pairwise_circular_ci <- function(df, angle_col = "MeanHeadingDeg", group_col = "Condition", id_col = "Pair", n_boot = 2000) {
  groups <- sort(unique(as.character(df[[group_col]])))
  if (length(groups) < 2) return(data.frame())
  out <- list()
  k <- 1
  for (i in seq_len(length(groups) - 1)) {
    for (j in (i + 1):length(groups)) {
      a <- groups[i]
      b <- groups[j]
      xa <- df[as.character(df[[group_col]]) == a, c(id_col, angle_col), drop = FALSE]
      xb <- df[as.character(df[[group_col]]) == b, c(id_col, angle_col), drop = FALSE]
      names(xa) <- c(id_col, "AngleA")
      names(xb) <- c(id_col, "AngleB")
      merged <- merge(xa, xb, by = id_col)
      merged <- merged[is.finite(merged$AngleA) & is.finite(merged$AngleB), , drop = FALSE]
      if (nrow(merged) == 0) next
      delta <- wrap_deg(merged$AngleA - merged$AngleB)
      boot <- rep(NA_real_, n_boot)
      if (length(delta) >= 2) {
        for (bidx in seq_len(n_boot)) {
          samp <- sample(delta, length(delta), replace = TRUE)
          boot[bidx] <- circular_mean_deg(samp)
        }
      }
      out[[k]] <- data.frame(
        GroupA = a,
        GroupB = b,
        N = nrow(merged),
        MeanDeltaDeg = circular_mean_deg(delta),
        RbarDelta = circular_rbar(delta),
        CI_lo = if (all(!is.na(boot))) as.numeric(stats::quantile(boot, 0.025, na.rm = TRUE)) else NA_real_,
        CI_hi = if (all(!is.na(boot))) as.numeric(stats::quantile(boot, 0.975, na.rm = TRUE)) else NA_real_,
        stringsAsFactors = FALSE
      )
      k <- k + 1
    }
  }
  if (length(out) == 0) return(data.frame())
  do.call(rbind, out)
}

choose_reference_condition <- function(df) {
  preferred <- c("controlsea1", "controlt1", "control", "sham")
  conds <- unique(as.character(df$Condition))
  hit <- preferred[preferred %in% conds]
  if (length(hit) > 0) return(hit[1])
  ord <- aggregate(OrderInSession ~ Condition, data = df, FUN = function(x) min(as.numeric(x), na.rm = TRUE))
  ord <- ord[order(ord$OrderInSession, ord$Condition), , drop = FALSE]
  as.character(ord$Condition[1])
}

individual_alignment_analysis <- function(df, reference_condition = NULL) {
  ref_cond <- reference_condition %||% choose_reference_condition(df)
  ref_rows <- df[df$Condition == ref_cond, c("Pair", "MeanHeadingDeg"), drop = FALSE]
  names(ref_rows)[2] <- "ReferenceHeadingDeg"
  aligned <- merge(df, ref_rows, by = "Pair", all.x = TRUE)
  aligned$ReferenceSource <- ifelse(is.finite(aligned$ReferenceHeadingDeg), ref_cond, "pair_mean")
  need_pair_mean <- !is.finite(aligned$ReferenceHeadingDeg)
  if (any(need_pair_mean)) {
    pair_mean <- aggregate(MeanHeadingDeg ~ Pair, data = df, FUN = circular_mean_deg)
    names(pair_mean)[2] <- "PairMeanHeadingDeg"
    aligned <- merge(aligned, pair_mean, by = "Pair", all.x = TRUE)
    aligned$ReferenceHeadingDeg[need_pair_mean] <- aligned$PairMeanHeadingDeg[need_pair_mean]
  } else {
    aligned$PairMeanHeadingDeg <- NA_real_
  }
  aligned$AlignedHeadingDeg <- wrap_deg(aligned$MeanHeadingDeg - aligned$ReferenceHeadingDeg)
  cond_levels <- unique(as.character(aligned$Condition))
  by_condition <- do.call(rbind, lapply(cond_levels, function(cond) {
    x <- aligned$AlignedHeadingDeg[as.character(aligned$Condition) == cond]
    boot <- rep(NA_real_, 2000)
    xv <- x[is.finite(x)]
    if (length(xv) >= 2) {
      for (i in seq_len(2000)) boot[i] <- circular_mean_deg(sample(xv, length(xv), replace = TRUE))
    }
    data.frame(
      Condition = cond,
      MeanAlignedDeg = circular_mean_deg(x),
      RbarAligned = circular_rbar(x),
      CI_lo = if (all(!is.na(boot))) as.numeric(stats::quantile(boot, 0.025, na.rm = TRUE)) else NA_real_,
      CI_hi = if (all(!is.na(boot))) as.numeric(stats::quantile(boot, 0.975, na.rm = TRUE)) else NA_real_,
      SignTestP = sign_test_p(x),
      PosCount = sum(x > 0, na.rm = TRUE),
      NegCount = sum(x < 0, na.rm = TRUE),
      N = sum(is.finite(x)),
      stringsAsFactors = FALSE
    )
  }))
  rownames(by_condition) <- NULL
  list(reference_condition = ref_cond, aligned = aligned, by_condition = by_condition)
}

fit_lmer_side <- function(df, response) {
  if (!requireNamespace("lme4", quietly = TRUE)) return(list(ok = FALSE, reason = "Package 'lme4' is not installed."))
  formula_txt <- paste0(response, " ~ Condition + OrderInSession + PrevCondition + CohortGroup + Condition:CohortGroup + (1|Cohort) + (1|Fish)")
  engine <- if (requireNamespace("lmerTest", quietly = TRUE)) "lmerTest" else "lme4"
  fit <- tryCatch(
    if (engine == "lmerTest") lmerTest::lmer(stats::as.formula(formula_txt), data = df, REML = FALSE)
    else lme4::lmer(stats::as.formula(formula_txt), data = df, REML = FALSE),
    error = function(e) e
  )
  if (inherits(fit, "error")) return(list(ok = FALSE, reason = conditionMessage(fit), formula = formula_txt, engine = engine))
  anova_tbl <- tryCatch(as.data.frame(stats::anova(fit)), error = function(e) NULL)
  if (!is.null(anova_tbl) && length(rownames(anova_tbl)) > 0) {
    anova_tbl$Term <- rownames(anova_tbl)
    rownames(anova_tbl) <- NULL
  }
  coef_tbl <- tryCatch(as.data.frame(summary(fit)$coefficients), error = function(e) NULL)
  if (!is.null(coef_tbl) && length(rownames(coef_tbl)) > 0) {
    coef_tbl$Term <- rownames(coef_tbl)
    rownames(coef_tbl) <- NULL
  }
  list(ok = TRUE, formula = formula_txt, model = fit, anova = anova_tbl, coefficients = coef_tbl, engine = engine)
}

write_model_result <- function(result, prefix, out_dir) {
  txt_path <- file.path(out_dir, paste0(prefix, "_summary.txt"))
  if (!isTRUE(result$ok)) {
    writeLines(c("status: failed", paste("formula:", result$formula %||% "n/a"), paste("reason:", result$reason %||% "unknown")), txt_path)
    return(invisible(NULL))
  }
  capture.output({
    cat("engine:", result$engine %||% "unknown", "\n")
    cat("formula:", result$formula, "\n\n")
    cat("anova\n")
    print(result$anova)
    cat("\ncoefficients\n")
    print(result$coefficients)
  }, file = txt_path)
  if (!is.null(result$anova)) write.csv(result$anova, file.path(out_dir, paste0(prefix, "_anova.csv")), row.names = FALSE)
  if (!is.null(result$coefficients)) write.csv(result$coefficients, file.path(out_dir, paste0(prefix, "_coefficients.csv")), row.names = FALSE)
}

extract_mixed_model_terms <- function(result, response_label) {
  if (!isTRUE(result$ok) || is.null(result$anova)) return(data.frame())
  A <- result$anova
  term_col <- if ("Term" %in% names(A)) "Term" else names(A)[1]
  p_col <- extract_table_col(A, c("^Pr\\(", "^p$", "p.value", "p_value", "prob"))
  f_col <- extract_table_col(A, c("^F", "F.value"))
  data.frame(
    Response = response_label,
    Term = as.character(A[[term_col]]),
    F_value = if (!is.na(f_col)) suppressWarnings(as.numeric(A[[f_col]])) else NA_real_,
    p_value = if (!is.na(p_col)) suppressWarnings(as.numeric(A[[p_col]])) else NA_real_,
    stringsAsFactors = FALSE
  )
}

fit_bayesian_circular <- function(df) {
  if (!requireNamespace("bpnreg", quietly = TRUE)) return(list(ok = FALSE, reason = "Package 'bpnreg' is not installed."))
  theta_rad <- df$MeanHeadingDeg * pi / 180
  ok <- is.finite(theta_rad)
  if (sum(ok) < 10) return(list(ok = FALSE, reason = "Too few finite mean angles for Bayesian circular model."))
  work <- df[ok, , drop = FALSE]
  work$theta <- theta_rad[ok]
  fit <- tryCatch(
    bpnreg::bpnr(pred.I = theta ~ Condition + OrderInSession + PrevCondition + CohortGroup, data = work, its = 2000, burn = 1000, seed = 1),
    error = function(e) e
  )
  if (inherits(fit, "error")) return(list(ok = FALSE, reason = conditionMessage(fit)))
  list(ok = TRUE, model = fit)
}

cohort_pairwise_analysis <- function(df, n_boot = 1000) {
  groups <- unique(as.character(df$CohortGroup))
  out <- list()
  for (grp in groups) {
    sub <- df[as.character(df$CohortGroup) == grp, , drop = FALSE]
    pw <- pairwise_circular_ci(sub, n_boot = n_boot)
    if (nrow(pw) == 0) next
    pw$CohortGroup <- grp
    out[[length(out) + 1]] <- pw
  }
  if (length(out) == 0) return(data.frame())
  out <- do.call(rbind, out)
  rownames(out) <- NULL
  out
}

cohort_consistency_summary <- function(cohort_pairwise, pooled_pairwise) {
  if (nrow(cohort_pairwise) == 0 || nrow(pooled_pairwise) == 0) return(data.frame())
  pair_key <- paste(pooled_pairwise$GroupA, pooled_pairwise$GroupB, sep = " vs ")
  pooled_sign <- sign(pooled_pairwise$MeanDeltaDeg)
  names(pooled_sign) <- pair_key
  cohort_pairwise$Pair <- paste(cohort_pairwise$GroupA, cohort_pairwise$GroupB, sep = " vs ")
  out <- by(cohort_pairwise, cohort_pairwise$Pair, function(df) {
    key <- unique(df$Pair)[1]
    ps <- pooled_sign[[key]]
    same_sign <- sum(sign(df$MeanDeltaDeg) == ps, na.rm = TRUE)
    data.frame(
      Pair = key,
      PooledMeanDeltaDeg = pooled_pairwise$MeanDeltaDeg[pair_key == key][1],
      CohortGroups = nrow(df),
      SameDirectionGroups = same_sign,
      SameDirectionFraction = same_sign / nrow(df),
      MinGroupDeltaDeg = min(df$MeanDeltaDeg, na.rm = TRUE),
      MaxGroupDeltaDeg = max(df$MeanDeltaDeg, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, out)
  rownames(out) <- NULL
  out
}

build_window_specs <- function(dat, base_start, base_end) {
  t_min <- floor(min(dat$Time[is.finite(dat$Time)], na.rm = TRUE))
  t_max <- ceiling(max(dat$Time[is.finite(dat$Time)], na.rm = TRUE))
  width <- base_end - base_start
  specs <- list(
    list(label = "base", start = base_start, end = base_end),
    list(label = "early_shift", start = max(t_min, base_start - 60), end = max(t_min + width, base_end - 60)),
    list(label = "late_shift", start = min(t_max - width, base_start + 60), end = min(t_max, base_end + 60)),
    list(label = "wide_context", start = max(t_min, base_start - 60), end = min(t_max, base_end + 60)),
    list(label = "narrow_core", start = base_start + floor(width / 4), end = base_end - floor(width / 4))
  )
  out <- list()
  seen <- character(0)
  for (sp in specs) {
    if (!is.finite(sp$start) || !is.finite(sp$end) || sp$end <= sp$start) next
    key <- paste(sp$start, sp$end, sep = "_")
    if (key %in% seen) next
    seen <- c(seen, key)
    out[[length(out) + 1]] <- sp
  }
  out
}

window_sensitivity_analysis <- function(dat, base_start, base_end, keep_conditions) {
  specs <- build_window_specs(dat, base_start, base_end)
  out <- list()
  meta <- list()
  for (sp in specs) {
    dat_win <- prepare_sample_data(dat, sp$start, sp$end, keep_conditions)
    agg_win <- aggregate_fish_condition(dat_win)
    if (nrow(agg_win) == 0) next
    agg_win <- finalize_agg(agg_win)
    pw <- pairwise_circular_ci(agg_win, n_boot = 1000)
    if (nrow(pw) == 0) next
    pw$WindowLabel <- sp$label
    pw$WindowStart <- sp$start
    pw$WindowEnd <- sp$end
    pw$Rows <- nrow(agg_win)
    out[[length(out) + 1]] <- pw
    meta[[length(meta) + 1]] <- data.frame(WindowLabel = sp$label, WindowStart = sp$start, WindowEnd = sp$end, AggregatedRows = nrow(agg_win), stringsAsFactors = FALSE)
  }
  list(pairwise = if (length(out) > 0) do.call(rbind, out) else data.frame(), windows = if (length(meta) > 0) do.call(rbind, meta) else data.frame())
}

quality_sensitivity_analysis <- function(agg) {
  if (nrow(agg) == 0) return(data.frame())
  ref_n <- stats::median(agg$N[is.finite(agg$N)], na.rm = TRUE)
  thresholds <- unique(sort(round(c(0, 0.5, 0.7, 0.9) * ref_n)))
  out <- list()
  for (thr in thresholds) {
    sub <- agg[agg$N >= thr, , drop = FALSE]
    if (nrow(sub) < 6) next
    pw <- pairwise_circular_ci(sub, n_boot = 1000)
    if (nrow(pw) == 0) next
    pw$MinSamplesPerFishCondition <- thr
    pw$RetainedRows <- nrow(sub)
    pw$RetainedFraction <- nrow(sub) / nrow(agg)
    out[[length(out) + 1]] <- pw
  }
  if (length(out) == 0) return(data.frame())
  out <- do.call(rbind, out)
  rownames(out) <- NULL
  out
}

pairwise_continuous_tests_by_group <- function(df, metrics, group_col = "Condition", id_col = "Pair", adjust_within = c("Metric")) {
  if (length(metrics) == 0 || nrow(df) == 0 || !(group_col %in% names(df)) || !(id_col %in% names(df))) return(data.frame())
  groups <- sort(unique(as.character(df[[group_col]])))
  if (length(groups) < 2) return(data.frame())
  out <- list()
  for (metric in metrics) {
    for (i in seq_len(length(groups) - 1)) {
      for (j in (i + 1):length(groups)) {
        a <- groups[i]
        b <- groups[j]
        xa <- df[as.character(df[[group_col]]) == a, c(id_col, metric), drop = FALSE]
        xb <- df[as.character(df[[group_col]]) == b, c(id_col, metric), drop = FALSE]
        names(xa) <- c(id_col, "A")
        names(xb) <- c(id_col, "B")
        merged <- merge(xa, xb, by = id_col)
        merged <- merged[is.finite(merged$A) & is.finite(merged$B), , drop = FALSE]
        if (nrow(merged) < 2) next
        tt <- tryCatch(stats::t.test(merged$A, merged$B, paired = TRUE), error = function(e) NULL)
        out[[length(out) + 1]] <- data.frame(
          Metric = metric,
          GroupA = a,
          GroupB = b,
          N = nrow(merged),
          MeanA = mean(merged$A, na.rm = TRUE),
          MeanB = mean(merged$B, na.rm = TRUE),
          MeanDiff = mean(merged$A - merged$B, na.rm = TRUE),
          CI_lo = if (!is.null(tt)) unname(tt$conf.int[1]) else NA_real_,
          CI_hi = if (!is.null(tt)) unname(tt$conf.int[2]) else NA_real_,
          t_value = if (!is.null(tt)) unname(tt$statistic) else NA_real_,
          p_value = if (!is.null(tt)) unname(tt$p.value) else NA_real_,
          stringsAsFactors = FALSE
        )
      }
    }
  }
  if (length(out) == 0) return(data.frame())
  out <- do.call(rbind, out)
  out$p_adj <- NA_real_
  if ("p_value" %in% names(out)) {
    if (identical(adjust_within, c("Metric")) && "Metric" %in% names(out)) {
      for (metric in unique(as.character(out$Metric))) {
        idx <- which(as.character(out$Metric) == metric & is.finite(out$p_value))
        if (length(idx) > 0) out$p_adj[idx] <- stats::p.adjust(out$p_value[idx], method = "holm")
      }
    } else {
      split_key <- interaction(lapply(adjust_within, function(nm) if (nm %in% names(out)) out[[nm]] else rep(NA_character_, nrow(out))), drop = TRUE, lex.order = TRUE)
      for (key in unique(split_key)) {
        idx <- which(split_key == key & is.finite(out$p_value))
        if (length(idx) > 0) out$p_adj[idx] <- stats::p.adjust(out$p_value[idx], method = "holm")
      }
    }
  }
  rownames(out) <- NULL
  out
}

pairwise_continuous_tests <- function(df, metrics) {
  pairwise_continuous_tests_by_group(df, metrics, group_col = "Condition", id_col = "Pair", adjust_within = c("Metric"))
}

add_significance_flags <- function(tbl, lo_col = "CI_lo", hi_col = "CI_hi", p_col = "p_adj", alpha = 0.05) {
  if (!is.data.frame(tbl) || nrow(tbl) == 0) return(tbl)
  tbl$SignificantByCI <- NA
  if (all(c(lo_col, hi_col) %in% names(tbl))) {
    tbl$SignificantByCI <- is.finite(tbl[[lo_col]]) & is.finite(tbl[[hi_col]]) & (tbl[[lo_col]] > 0 | tbl[[hi_col]] < 0)
  }
  tbl$SignificantByP <- NA
  if (p_col %in% names(tbl)) {
    tbl$SignificantByP <- is.finite(tbl[[p_col]]) & tbl[[p_col]] < alpha
  }
  tbl
}

add_condition_p_adjustment <- function(tbl, method = "holm") {
  if (!is.data.frame(tbl) || nrow(tbl) == 0 || !all(c("Term", "p_value") %in% names(tbl))) return(tbl)
  tbl$p_adj <- NA_real_
  idx <- which(tbl$Term == "Condition" & is.finite(tbl$p_value))
  if (length(idx) > 0) tbl$p_adj[idx] <- stats::p.adjust(tbl$p_value[idx], method = method)
  tbl
}

cohort_pairwise_continuous_analysis <- function(df, metrics) {
  if (length(metrics) == 0 || nrow(df) == 0 || !("CohortGroup" %in% names(df))) return(data.frame())
  out <- list()
  groups <- unique(as.character(df$CohortGroup))
  for (grp in groups) {
    sub <- df[as.character(df$CohortGroup) == grp, , drop = FALSE]
    res <- pairwise_continuous_tests(sub, metrics)
    if (!is.data.frame(res) || nrow(res) == 0) next
    res$CohortGroup <- grp
    out[[length(out) + 1]] <- res
  }
  if (length(out) == 0) return(data.frame())
  out <- do.call(rbind, out)
  rownames(out) <- NULL
  out
}

age_group_pairwise_analysis <- function(df, n_boot = 1000) {
  if (nrow(df) == 0 || !("AgeGroup" %in% names(df))) return(data.frame())
  age_groups <- unique(as.character(df$AgeGroup))
  age_groups <- age_groups[!(is.na(age_groups) | age_groups == "")]
  out <- list()
  for (grp in age_groups) {
    sub <- df[as.character(df$AgeGroup) == grp, , drop = FALSE]
    pw <- pairwise_circular_ci(sub, n_boot = n_boot)
    if (nrow(pw) == 0) next
    pw$AgeGroup <- grp
    out[[length(out) + 1]] <- pw
  }
  if (length(out) == 0) return(data.frame())
  out <- do.call(rbind, out)
  rownames(out) <- NULL
  out
}

age_group_pairwise_continuous_analysis <- function(df, metrics) {
  if (length(metrics) == 0 || nrow(df) == 0 || !("AgeGroup" %in% names(df))) return(data.frame())
  out <- list()
  age_groups <- unique(as.character(df$AgeGroup))
  age_groups <- age_groups[!(is.na(age_groups) | age_groups == "")]
  for (grp in age_groups) {
    sub <- df[as.character(df$AgeGroup) == grp, , drop = FALSE]
    res <- pairwise_continuous_tests(sub, metrics)
    if (!is.data.frame(res) || nrow(res) == 0) next
    res$AgeGroup <- grp
    out[[length(out) + 1]] <- res
  }
  if (length(out) == 0) return(data.frame())
  out <- do.call(rbind, out)
  rownames(out) <- NULL
  out
}

fit_age_lmer_side <- function(df, response) {
  if (!("AgeGroup" %in% names(df))) return(list(ok = FALSE, reason = "AgeGroup column is not available.", formula = NA_character_))
  work <- df[!(is.na(df$AgeGroup) | as.character(df$AgeGroup) == ""), , drop = FALSE]
  if (nrow(work) == 0 || length(unique(as.character(work$AgeGroup))) < 2) {
    return(list(ok = FALSE, reason = "AgeGroup requires at least two non-empty groups.", formula = NA_character_))
  }
  if (!requireNamespace("lme4", quietly = TRUE)) return(list(ok = FALSE, reason = "Package 'lme4' is not installed."))
  formula_txt <- paste0(response, " ~ Condition * AgeGroup + OrderInSession + PrevCondition + (1|Cohort) + (1|Fish)")
  engine <- if (requireNamespace("lmerTest", quietly = TRUE)) "lmerTest" else "lme4"
  fit <- tryCatch(
    if (engine == "lmerTest") lmerTest::lmer(stats::as.formula(formula_txt), data = work, REML = FALSE)
    else lme4::lmer(stats::as.formula(formula_txt), data = work, REML = FALSE),
    error = function(e) e
  )
  if (inherits(fit, "error")) return(list(ok = FALSE, reason = conditionMessage(fit), formula = formula_txt, engine = engine))
  anova_tbl <- tryCatch(as.data.frame(stats::anova(fit)), error = function(e) NULL)
  if (!is.null(anova_tbl) && length(rownames(anova_tbl)) > 0) {
    anova_tbl$Term <- rownames(anova_tbl)
    rownames(anova_tbl) <- NULL
  }
  coef_tbl <- tryCatch(as.data.frame(summary(fit)$coefficients), error = function(e) NULL)
  if (!is.null(coef_tbl) && length(rownames(coef_tbl)) > 0) {
    coef_tbl$Term <- rownames(coef_tbl)
    rownames(coef_tbl) <- NULL
  }
  list(ok = TRUE, formula = formula_txt, model = fit, anova = anova_tbl, coefficients = coef_tbl, engine = engine)
}

extract_age_terms <- function(result, response_label) {
  out <- extract_mixed_model_terms(result, response_label)
  if (!is.data.frame(out) || nrow(out) == 0) return(out)
  keep <- grepl("^Condition$", out$Term) | grepl("^AgeGroup", out$Term) | grepl("^Condition:AgeGroup", out$Term)
  out[keep, , drop = FALSE]
}

order_pairwise_continuous_analysis <- function(df, metrics) {
  if (length(metrics) == 0 || nrow(df) == 0 || !("OrderInSession" %in% names(df))) return(data.frame())
  work <- df[is.finite(suppressWarnings(as.numeric(df$OrderInSession))), , drop = FALSE]
  if (nrow(work) == 0) return(data.frame())
  work$OrderLabel <- paste0("session", as.integer(work$OrderInSession))
  out <- pairwise_continuous_tests_by_group(work, metrics, group_col = "OrderLabel", id_col = "Pair", adjust_within = c("Metric"))
  if (!is.data.frame(out) || nrow(out) == 0) return(out)
  out$OrderA <- suppressWarnings(as.numeric(gsub("^session", "", out$GroupA)))
  out$OrderB <- suppressWarnings(as.numeric(gsub("^session", "", out$GroupB)))
  out
}

fit_order_lmer_side <- function(df, response) {
  if (!("OrderInSession" %in% names(df))) return(list(ok = FALSE, reason = "OrderInSession column is not available.", formula = NA_character_))
  work <- df[is.finite(suppressWarnings(as.numeric(df$OrderInSession))), , drop = FALSE]
  if (nrow(work) == 0 || length(unique(as.numeric(work$OrderInSession))) < 2) {
    return(list(ok = FALSE, reason = "OrderInSession requires at least two levels.", formula = NA_character_))
  }
  if (!requireNamespace("lme4", quietly = TRUE)) return(list(ok = FALSE, reason = "Package 'lme4' is not installed."))
  work$OrderInSession <- factor(as.integer(work$OrderInSession))
  formula_txt <- paste0(response, " ~ OrderInSession + (1|Cohort) + (1|Fish)")
  engine <- if (requireNamespace("lmerTest", quietly = TRUE)) "lmerTest" else "lme4"
  fit <- tryCatch(
    if (engine == "lmerTest") lmerTest::lmer(stats::as.formula(formula_txt), data = work, REML = FALSE)
    else lme4::lmer(stats::as.formula(formula_txt), data = work, REML = FALSE),
    error = function(e) e
  )
  if (inherits(fit, "error")) return(list(ok = FALSE, reason = conditionMessage(fit), formula = formula_txt, engine = engine))
  anova_tbl <- tryCatch(as.data.frame(stats::anova(fit)), error = function(e) NULL)
  if (!is.null(anova_tbl) && length(rownames(anova_tbl)) > 0) {
    anova_tbl$Term <- rownames(anova_tbl)
    rownames(anova_tbl) <- NULL
  }
  coef_tbl <- tryCatch(as.data.frame(summary(fit)$coefficients), error = function(e) NULL)
  if (!is.null(coef_tbl) && length(rownames(coef_tbl)) > 0) {
    coef_tbl$Term <- rownames(coef_tbl)
    rownames(coef_tbl) <- NULL
  }
  list(ok = TRUE, formula = formula_txt, model = fit, anova = anova_tbl, coefficients = coef_tbl, engine = engine)
}

extract_order_terms <- function(result, response_label) {
  out <- extract_mixed_model_terms(result, response_label)
  if (!is.data.frame(out) || nrow(out) == 0) return(out)
  out[grepl("^OrderInSession", out$Term), , drop = FALSE]
}

fit_subset_order_lmer_side <- function(df, response, family_label, age_group = NULL) {
  if (!("OrderInSession" %in% names(df)) || !("ConditionFamily" %in% names(df))) {
    return(list(ok = FALSE, reason = "OrderInSession or ConditionFamily column is not available.", formula = NA_character_))
  }
  work <- df[as.character(df$ConditionFamily) == family_label, , drop = FALSE]
  if (!is.null(age_group) && "AgeGroup" %in% names(work)) {
    work <- work[as.character(work$AgeGroup) == age_group, , drop = FALSE]
  }
  work <- work[is.finite(suppressWarnings(as.numeric(work$OrderInSession))), , drop = FALSE]
  if (nrow(work) == 0) {
    return(list(ok = FALSE, reason = "No rows remained after subset filtering.", formula = NA_character_))
  }
  if (length(unique(as.numeric(work$OrderInSession))) < 2) {
    return(list(ok = FALSE, reason = "OrderInSession requires at least two levels in the subset.", formula = NA_character_))
  }
  if (length(unique(as.character(work$Condition))) < 2) {
    return(list(ok = FALSE, reason = "Subset requires at least two conditions to separate order from condition structure.", formula = NA_character_))
  }
  if (!requireNamespace("lme4", quietly = TRUE)) return(list(ok = FALSE, reason = "Package 'lme4' is not installed."))
  work$OrderInSession <- factor(as.integer(work$OrderInSession))
  work$Condition <- factor(as.character(work$Condition))
  formula_txt <- paste0(response, " ~ OrderInSession + (1|Condition) + (1|Cohort) + (1|Fish)")
  engine <- if (requireNamespace("lmerTest", quietly = TRUE)) "lmerTest" else "lme4"
  fit <- tryCatch(
    if (engine == "lmerTest") lmerTest::lmer(stats::as.formula(formula_txt), data = work, REML = FALSE)
    else lme4::lmer(stats::as.formula(formula_txt), data = work, REML = FALSE),
    error = function(e) e
  )
  if (inherits(fit, "error")) return(list(ok = FALSE, reason = conditionMessage(fit), formula = formula_txt, engine = engine))
  anova_tbl <- tryCatch(as.data.frame(stats::anova(fit)), error = function(e) NULL)
  if (!is.null(anova_tbl) && length(rownames(anova_tbl)) > 0) {
    anova_tbl$Term <- rownames(anova_tbl)
    rownames(anova_tbl) <- NULL
  }
  coef_tbl <- tryCatch(as.data.frame(summary(fit)$coefficients), error = function(e) NULL)
  if (!is.null(coef_tbl) && length(rownames(coef_tbl)) > 0) {
    coef_tbl$Term <- rownames(coef_tbl)
    rownames(coef_tbl) <- NULL
  }
  list(ok = TRUE, formula = formula_txt, model = fit, anova = anova_tbl, coefficients = coef_tbl, engine = engine, subset_n = nrow(work))
}

subset_order_activity_analysis <- function(df, metrics, families = c("Sea", "Intensity")) {
  if (length(metrics) == 0 || nrow(df) == 0) {
    return(list(pairwise = data.frame(), model_terms = data.frame(), model_meta = data.frame(), model_results = list()))
  }
  work <- df
  work$ConditionFamily <- condition_family(work$Condition)
  age_groups <- if ("AgeGroup" %in% names(work)) unique(as.character(work$AgeGroup)) else character(0)
  age_groups <- age_groups[!(is.na(age_groups) | age_groups == "")]
  subset_specs <- list()
  for (fam in families) {
    for (age in age_groups) subset_specs[[length(subset_specs) + 1]] <- list(Family = fam, AgeGroup = age)
  }
  pairwise_rows <- list()
  model_term_rows <- list()
  meta_rows <- list()
  model_results <- list()
  for (spec in subset_specs) {
    fam <- spec$Family
    age <- spec$AgeGroup
    sub <- work[as.character(work$ConditionFamily) == fam & as.character(work$AgeGroup) == age, , drop = FALSE]
    sub <- sub[is.finite(suppressWarnings(as.numeric(sub$OrderInSession))), , drop = FALSE]
    if (nrow(sub) > 0) {
      sub$OrderLabel <- paste0("session", as.integer(sub$OrderInSession))
      pw <- pairwise_continuous_tests_by_group(sub, metrics, group_col = "OrderLabel", id_col = "Pair", adjust_within = c("Subset", "Metric"))
      if (is.data.frame(pw) && nrow(pw) > 0) {
        pw$Subset <- paste(fam, age, sep = "_")
        pw$ConditionFamily <- fam
        pw$AgeGroup <- age
        pairwise_rows[[length(pairwise_rows) + 1]] <- pw
      }
    }
    for (metric in metrics) {
      res <- fit_subset_order_lmer_side(work, metric, family_label = fam, age_group = age)
      key <- paste(fam, age, metric, sep = "__")
      model_results[[key]] <- res
      meta_rows[[length(meta_rows) + 1]] <- data.frame(
        ConditionFamily = fam,
        AgeGroup = age,
        Response = metric,
        OK = isTRUE(res$ok),
        N = res$subset_n %||% NA_real_,
        Formula = res$formula %||% NA_character_,
        Reason = if (isTRUE(res$ok)) "" else res$reason %||% "unknown",
        stringsAsFactors = FALSE
      )
      if (isTRUE(res$ok)) {
        terms <- extract_order_terms(res, metric)
        if (is.data.frame(terms) && nrow(terms) > 0) {
          terms$ConditionFamily <- fam
          terms$AgeGroup <- age
          model_term_rows[[length(model_term_rows) + 1]] <- terms
        }
      }
    }
  }
  list(
    pairwise = if (length(pairwise_rows) > 0) do.call(rbind, pairwise_rows) else data.frame(),
    model_terms = if (length(model_term_rows) > 0) do.call(rbind, model_term_rows) else data.frame(),
    model_meta = if (length(meta_rows) > 0) do.call(rbind, meta_rows) else data.frame(),
    model_results = model_results
  )
}

within_fish_consistency_analysis <- function(agg, alignment) {
  ref_cond <- alignment$reference_condition
  ref_df <- agg[as.character(agg$Condition) == ref_cond, c("Pair", "MeanHeadingDeg"), drop = FALSE]
  names(ref_df)[2] <- "ReferenceHeadingDeg"
  others <- setdiff(unique(as.character(agg$Condition)), ref_cond)
  corr_rows <- list()
  shift_rows <- list()
  pair_rows <- list()
  pair_cohort_rows <- list()
  pair_age_rows <- list()
  shift_distribution_rows <- list()
  shift_cohort_rows <- list()
  shift_age_rows <- list()
  pair_age_map <- if ("AgeGroup" %in% names(agg)) unique(agg[, c("Pair", "AgeGroup"), drop = FALSE]) else data.frame()
  for (cond in others) {
    cond_df <- agg[as.character(agg$Condition) == cond, c("Pair", "MeanHeadingDeg"), drop = FALSE]
    names(cond_df)[2] <- "ConditionHeadingDeg"
    merged <- merge(ref_df, cond_df, by = "Pair")
    cc <- circular_correlation(merged$ReferenceHeadingDeg, merged$ConditionHeadingDeg, n_perm = 2000)
    corr_rows[[length(corr_rows) + 1]] <- data.frame(ReferenceCondition = ref_cond, Condition = cond, CircularCorrelation = cc$r, PermutationP = cc$p, N = cc$n, stringsAsFactors = FALSE)
    aligned_sub <- alignment$aligned[as.character(alignment$aligned$Condition) == cond, c("Pair", "CohortGroup", "AlignedHeadingDeg"), drop = FALSE]
    if (nrow(pair_age_map) > 0) aligned_sub <- merge(aligned_sub, pair_age_map, by = "Pair", all.x = TRUE, sort = FALSE)
    names(aligned_sub)[names(aligned_sub) == "AlignedHeadingDeg"] <- "ShiftDeg"
    aligned_sub <- aligned_sub[is.finite(aligned_sub$ShiftDeg), , drop = FALSE]
    aligned_sub$ReferenceCondition <- ref_cond
    aligned_sub$Condition <- cond
    keep_cols <- intersect(c("ReferenceCondition", "Condition", "CohortGroup", "AgeGroup", "Pair", "ShiftDeg"), names(aligned_sub))
    shift_distribution_rows[[length(shift_distribution_rows) + 1]] <- aligned_sub[, keep_cols, drop = FALSE]
    x <- aligned_sub$ShiftDeg
    shift_rows[[length(shift_rows) + 1]] <- within(
      summarize_shift_distribution(x, reference_condition = ref_cond, condition = cond, cohort_group = NA_character_),
      RbarShift <- circular_rbar(x)
    )
    if (nrow(aligned_sub) > 0 && "CohortGroup" %in% names(aligned_sub)) {
      for (grp in unique(as.character(aligned_sub$CohortGroup))) {
        x_grp <- aligned_sub$ShiftDeg[as.character(aligned_sub$CohortGroup) == grp]
        shift_cohort_rows[[length(shift_cohort_rows) + 1]] <- summarize_shift_distribution(
          x_grp,
          reference_condition = ref_cond,
          condition = cond,
          cohort_group = grp
        )
      }
    }
    if (nrow(aligned_sub) > 0 && "AgeGroup" %in% names(aligned_sub)) {
      for (grp in unique(as.character(aligned_sub$AgeGroup))) {
        if (is.na(grp) || grp == "") next
        x_grp <- aligned_sub$ShiftDeg[as.character(aligned_sub$AgeGroup) == grp]
        age_row <- summarize_shift_distribution(
          x_grp,
          reference_condition = ref_cond,
          condition = cond,
          cohort_group = NA_character_
        )
        age_row$AgeGroup <- grp
        shift_age_rows[[length(shift_age_rows) + 1]] <- age_row
      }
    }
  }
  if (length(others) >= 2) {
    for (i in seq_len(length(others) - 1)) {
      for (j in (i + 1):length(others)) {
        a <- others[i]
        b <- others[j]
        sa <- alignment$aligned[as.character(alignment$aligned$Condition) == a, c("Pair", "AlignedHeadingDeg"), drop = FALSE]
        sb <- alignment$aligned[as.character(alignment$aligned$Condition) == b, c("Pair", "AlignedHeadingDeg"), drop = FALSE]
        names(sa)[2] <- "ShiftA"
        names(sb)[2] <- "ShiftB"
        merged <- merge(sa, sb, by = "Pair")
        merged <- merged[is.finite(merged$ShiftA) & is.finite(merged$ShiftB), , drop = FALSE]
        if (nrow(merged) < 5) next
        pair_rows[[length(pair_rows) + 1]] <- data.frame(ConditionA = a, ConditionB = b, PearsonR = suppressWarnings(stats::cor(merged$ShiftA, merged$ShiftB, method = "pearson")), SpearmanR = suppressWarnings(stats::cor(merged$ShiftA, merged$ShiftB, method = "spearman")), SameSignFraction = mean(sign(merged$ShiftA) == sign(merged$ShiftB)), N = nrow(merged), stringsAsFactors = FALSE)
        if ("CohortGroup" %in% names(alignment$aligned)) {
          sa2 <- alignment$aligned[as.character(alignment$aligned$Condition) == a, c("Pair", "CohortGroup", "AlignedHeadingDeg"), drop = FALSE]
          sb2 <- alignment$aligned[as.character(alignment$aligned$Condition) == b, c("Pair", "CohortGroup", "AlignedHeadingDeg"), drop = FALSE]
          names(sa2)[names(sa2) == "AlignedHeadingDeg"] <- "ShiftA"
          names(sb2)[names(sb2) == "AlignedHeadingDeg"] <- "ShiftB"
          merged2 <- merge(sa2, sb2[, c("Pair", "ShiftB"), drop = FALSE], by = "Pair")
          merged2 <- merged2[is.finite(merged2$ShiftA) & is.finite(merged2$ShiftB), , drop = FALSE]
          for (grp in unique(as.character(merged2$CohortGroup))) {
            subg <- merged2[as.character(merged2$CohortGroup) == grp, , drop = FALSE]
            if (nrow(subg) < 5) next
            pear <- permutation_correlation_test(subg$ShiftA, subg$ShiftB, n_perm = 2000, method = "pearson")
            spear <- permutation_correlation_test(subg$ShiftA, subg$ShiftB, n_perm = 2000, method = "spearman")
            pair_cohort_rows[[length(pair_cohort_rows) + 1]] <- data.frame(
              CohortGroup = grp,
              ConditionA = a,
              ConditionB = b,
              PearsonR = pear$r,
              PearsonP = pear$p,
              SpearmanR = spear$r,
              SpearmanP = spear$p,
              SameSignFraction = mean(sign(subg$ShiftA) == sign(subg$ShiftB)),
              N = nrow(subg),
              stringsAsFactors = FALSE
            )
          }
        }
        if (nrow(pair_age_map) > 0) {
          merged_age <- merge(merged, pair_age_map, by = "Pair", all.x = TRUE, sort = FALSE)
          for (grp in unique(as.character(merged_age$AgeGroup))) {
            if (is.na(grp) || grp == "") next
            subg <- merged_age[as.character(merged_age$AgeGroup) == grp, , drop = FALSE]
            if (nrow(subg) < 5) next
            pear <- permutation_correlation_test(subg$ShiftA, subg$ShiftB, n_perm = 2000, method = "pearson")
            spear <- permutation_correlation_test(subg$ShiftA, subg$ShiftB, n_perm = 2000, method = "spearman")
            pair_age_rows[[length(pair_age_rows) + 1]] <- data.frame(
              AgeGroup = grp,
              ConditionA = a,
              ConditionB = b,
              PearsonR = pear$r,
              PearsonP = pear$p,
              SpearmanR = spear$r,
              SpearmanP = spear$p,
              SameSignFraction = mean(sign(subg$ShiftA) == sign(subg$ShiftB)),
              N = nrow(subg),
              stringsAsFactors = FALSE
            )
          }
        }
      }
    }
  }
  list(
    reference_correlations = if (length(corr_rows) > 0) do.call(rbind, corr_rows) else data.frame(),
    shift_summary = if (length(shift_rows) > 0) do.call(rbind, shift_rows) else data.frame(),
    shift_pair_correlations = if (length(pair_rows) > 0) do.call(rbind, pair_rows) else data.frame(),
    shift_pair_cohort_correlations = if (length(pair_cohort_rows) > 0) do.call(rbind, pair_cohort_rows) else data.frame(),
    shift_distribution = if (length(shift_distribution_rows) > 0) do.call(rbind, shift_distribution_rows) else data.frame(),
    shift_cohort_summary = if (length(shift_cohort_rows) > 0) do.call(rbind, shift_cohort_rows) else data.frame(),
    shift_age_summary = if (length(shift_age_rows) > 0) do.call(rbind, shift_age_rows) else data.frame(),
    shift_pair_age_correlations = if (length(pair_age_rows) > 0) do.call(rbind, pair_age_rows) else data.frame()
  )
}

prepare_shift_cluster_data <- function(within_fish) {
  dat <- within_fish$shift_distribution
  if (!is.data.frame(dat) || nrow(dat) == 0) return(data.frame())
  conds <- sort(unique(as.character(dat$Condition)))
  if (length(conds) == 0) return(data.frame())
  pairs <- sort(unique(as.character(dat$Pair)))
  out <- data.frame(Pair = pairs, stringsAsFactors = FALSE)
  cohort_map <- unique(dat[, c("Pair", "CohortGroup"), drop = FALSE])
  out <- merge(out, cohort_map, by = "Pair", all.x = TRUE, sort = FALSE)
  ref_cond <- unique(as.character(dat$ReferenceCondition))
  out$ReferenceCondition <- if (length(ref_cond) > 0) ref_cond[1] else NA_character_
  for (cond in conds) {
    sub <- dat[as.character(dat$Condition) == cond, c("Pair", "ShiftDeg"), drop = FALSE]
    names(sub)[2] <- paste0("Shift_", cond)
    out <- merge(out, sub, by = "Pair", all.x = TRUE, sort = FALSE)
  }
  shift_cols <- grep("^Shift_", names(out), value = TRUE)
  out <- out[stats::complete.cases(out[, shift_cols, drop = FALSE]), , drop = FALSE]
  rownames(out) <- NULL
  out
}

evaluate_shift_clusters <- function(cluster_data, n_perm = 2000, nstart = 50) {
  shift_cols <- grep("^Shift_", names(cluster_data), value = TRUE)
  if (!is.data.frame(cluster_data) || length(shift_cols) == 0 || nrow(cluster_data) < 8) {
    return(list(assignments = data.frame(), summary = data.frame(), stats = data.frame(), centers = data.frame()))
  }
  X_raw <- as.matrix(cluster_data[, shift_cols, drop = FALSE])
  X <- scale(X_raw)
  if (any(!is.finite(X))) {
    return(list(assignments = data.frame(), summary = data.frame(), stats = data.frame(), centers = data.frame()))
  }
  set.seed(1)
  km2 <- stats::kmeans(X, centers = 2, nstart = nstart)
  obs_ratio <- km2$betweenss / km2$totss
  perm_ratio <- replicate(n_perm, {
    Xp <- X
    for (j in seq_len(ncol(Xp))) Xp[, j] <- sample(Xp[, j], length(Xp[, j]), replace = FALSE)
    km_perm <- stats::kmeans(Xp, centers = 2, nstart = 10)
    km_perm$betweenss / km_perm$totss
  })
  p_value <- mean(perm_ratio >= obs_ratio, na.rm = TRUE)

  centers_raw <- rowsum(X_raw, km2$cluster) / as.vector(table(km2$cluster))
  center_norm <- sqrt(rowSums(centers_raw^2))
  responder_id <- as.integer(names(which.max(center_norm))[1])
  cluster_role <- ifelse(seq_len(nrow(centers_raw)) == responder_id, "responder", "non_responder")

  assignments <- cluster_data
  assignments$Cluster <- paste0("cluster_", km2$cluster)
  assignments$ClusterRole <- ifelse(km2$cluster == responder_id, "responder", "non_responder")
  assignments$DistanceFromZero <- sqrt(rowSums(X_raw^2))

  summary_rows <- lapply(seq_len(nrow(centers_raw)), function(i) {
    idx <- which(km2$cluster == i)
    row <- data.frame(
      Cluster = paste0("cluster_", i),
      ClusterRole = cluster_role[i],
      N = length(idx),
      MeanDistanceFromZero = mean(sqrt(rowSums(X_raw[idx, , drop = FALSE]^2))),
      stringsAsFactors = FALSE
    )
    for (j in seq_along(shift_cols)) {
      row[[paste0("Mean_", shift_cols[j])]] <- mean(X_raw[idx, j], na.rm = TRUE)
      row[[paste0("Median_", shift_cols[j])]] <- stats::median(X_raw[idx, j], na.rm = TRUE)
    }
    row
  })
  summary_tbl <- do.call(rbind, summary_rows)

  centers_tbl <- data.frame(
    Cluster = paste0("cluster_", seq_len(nrow(centers_raw))),
    ClusterRole = cluster_role,
    stringsAsFactors = FALSE
  )
  for (j in seq_along(shift_cols)) centers_tbl[[shift_cols[j]]] <- centers_raw[, j]

  stats_tbl <- data.frame(
    Method = "kmeans_2clusters_on_shiftdeg",
    N = nrow(cluster_data),
    Dimensions = length(shift_cols),
    BetweenTotalSS = obs_ratio,
    PermutationP = p_value,
    SuggestedSubpopulation = is.finite(p_value) && p_value < 0.05,
    stringsAsFactors = FALSE
  )
  list(assignments = assignments, summary = summary_tbl, stats = stats_tbl, centers = centers_tbl)
}

plot_shift_clusters <- function(cluster_assignments, cluster_centers, out_dir) {
  if (!is.data.frame(cluster_assignments) || nrow(cluster_assignments) == 0) return(invisible(NULL))
  shift_cols <- grep("^Shift_", names(cluster_assignments), value = TRUE)
  if (length(shift_cols) < 2) return(invisible(NULL))
  xcol <- shift_cols[1]
  ycol <- shift_cols[2]
  png_path <- file.path(out_dir, "within_fish_shift_clusters.png")
  grDevices::png(png_path, width = 1800, height = 1400, res = 180)
  on.exit(grDevices::dev.off(), add = TRUE)
  roles <- as.character(cluster_assignments$ClusterRole)
  cols <- ifelse(roles == "responder", "#d95f02", "#1b9e77")
  cohorts <- as.character(cluster_assignments$CohortGroup)
  cohort_levels <- unique(cohorts)
  pchs <- setNames(seq_along(cohort_levels) + 14, cohort_levels)
  graphics::plot(
    cluster_assignments[[xcol]], cluster_assignments[[ycol]],
    xlab = paste0(sub("^Shift_", "", xcol), " shift (deg)"),
    ylab = paste0(sub("^Shift_", "", ycol), " shift (deg)"),
    main = "Within-fish shift clustering",
    pch = unname(pchs[cohorts]),
    col = cols,
    bg = cols
  )
  graphics::abline(h = 0, v = 0, lty = 2, col = "gray50")
  if (is.data.frame(cluster_centers) && nrow(cluster_centers) > 0) {
    graphics::points(cluster_centers[[xcol]], cluster_centers[[ycol]], pch = 8, cex = 2.2, lwd = 2)
    graphics::text(cluster_centers[[xcol]], cluster_centers[[ycol]], labels = cluster_centers$ClusterRole, pos = 3, cex = 0.9)
  }
  graphics::legend("topright",
    legend = c("responder", "non_responder"),
    col = c("#d95f02", "#1b9e77"),
    pch = 16, bty = "n"
  )
  graphics::legend("bottomright",
    legend = cohort_levels,
    pch = unname(pchs[cohort_levels]),
    title = "CohortGroup",
    bty = "n",
    cex = 0.8
  )
  invisible(png_path)
}

extract_bayes_table <- function(tbl) {
  if (is.null(tbl) || !is.data.frame(tbl) || nrow(tbl) == 0) return(data.frame())
  out <- tbl
  out$Term <- rownames(tbl)
  rownames(out) <- NULL
  mean_col <- extract_table_col(out, c("^mean$", "mean", "estimate"))
  lo_col <- extract_table_col(out, c("2.5", "l-95", "lower", "lo"))
  hi_col <- extract_table_col(out, c("97.5", "u-95", "upper", "hi"))
  if (!is.na(mean_col)) out$MeanValue <- suppressWarnings(as.numeric(out[[mean_col]]))
  if (!is.na(lo_col)) out$CI_lo <- suppressWarnings(as.numeric(out[[lo_col]]))
  if (!is.na(hi_col)) out$CI_hi <- suppressWarnings(as.numeric(out[[hi_col]]))
  out
}

write_bayesian_outputs <- function(bayes_result, out_dir) {
  bayes_txt <- file.path(out_dir, "bayesian_circular_summary.txt")
  if (!isTRUE(bayes_result$ok)) {
    writeLines(c("status: skipped", paste("reason:", bayes_result$reason %||% "unknown")), bayes_txt)
    return(list(model_fit = data.frame(), coef = data.frame(), coef_cat = data.frame(), coef_means = data.frame()))
  }
  bfit <- bayes_result$model
  circ_coef <- extract_bayes_table(if (!is.null(bfit$circ.coef)) as.data.frame(bfit$circ.coef) else NULL)
  circ_coef_cat <- extract_bayes_table(if (!is.null(bfit$circ.coef.cat)) as.data.frame(bfit$circ.coef.cat) else NULL)
  circ_coef_means <- extract_bayes_table(if (!is.null(bfit$circ.coef.means)) as.data.frame(bfit$circ.coef.means) else NULL)
  model_fit <- if (!is.null(bfit$model.fit)) as.data.frame(bfit$model.fit) else data.frame()
  if (nrow(circ_coef) > 0) write.csv(circ_coef, file.path(out_dir, "bayesian_circular_coefficients.csv"), row.names = FALSE)
  if (nrow(circ_coef_cat) > 0) write.csv(circ_coef_cat, file.path(out_dir, "bayesian_circular_coefficients_categorical.csv"), row.names = FALSE)
  if (nrow(circ_coef_means) > 0) write.csv(circ_coef_means, file.path(out_dir, "bayesian_circular_coefficients_means.csv"), row.names = FALSE)
  if (nrow(model_fit) > 0) write.csv(model_fit, file.path(out_dir, "bayesian_circular_model_fit.csv"), row.names = FALSE)
  capture.output({
    cat("status: fitted\n\n")
    if (nrow(model_fit) > 0) { cat("model_fit\n"); print(model_fit); cat("\n") }
    if (nrow(circ_coef_means) > 0) { cat("circ.coef.means\n"); print(circ_coef_means); cat("\n") }
    if (nrow(circ_coef_cat) > 0) { cat("circ.coef.cat\n"); print(circ_coef_cat); cat("\n") }
    if (nrow(circ_coef) > 0) { cat("circ.coef\n"); print(circ_coef); cat("\n") }
  }, file = bayes_txt)
  list(model_fit = model_fit, coef = circ_coef, coef_cat = circ_coef_cat, coef_means = circ_coef_means)
}

interpret_mixed_terms <- function(mixed_terms, response) {
  sub <- mixed_terms[mixed_terms$Response == response, , drop = FALSE]
  if (nrow(sub) == 0) return(sprintf("%s model could not be summarized.", response))
  wanted <- c("Condition", "OrderInSession", "PrevCondition", "CohortGroup", "Condition:CohortGroup")
  bits <- c()
  for (term in wanted) {
    idx <- which(sub$Term == term)
    if (length(idx) == 0) next
    bits <- c(bits, sprintf("%s p=%s", term, format_p(sub$p_value[idx[1]])))
  }
  if (length(bits) == 0) sprintf("%s model did not provide term-wise p-values.", response) else sprintf("%s: %s", response, paste(bits, collapse = " | "))
}

pretty_condition <- function(x) {
  key <- tolower(as.character(x))
  out <- ifelse(
    key == "controlsea1", "controlSea1",
    ifelse(
      key == "beringsea", "beringSea",
      ifelse(
        key == "okhotsksea", "okhotskSea",
        ifelse(key == "controlt1", "controlT1", ifelse(key == "microt10", "microT10", ifelse(key == "microt100", "microT100", as.character(x))))
      )
    )
  )
  unname(out)
}

condition_description <- function(x) {
  key <- tolower(as.character(x))
  out <- ifelse(
    key == "controlsea1", "controlsea1 (\u5bfe\u7167)",
    ifelse(
      key == "beringsea", "beringsea (\u30d9\u30fc\u30ea\u30f3\u30b0\u6d77\u78c1\u5834)",
      ifelse(
        key == "okhotsksea", "okhotsksea (\u30aa\u30db\u30fc\u30c4\u30af\u6d77\u78c1\u5834)",
        ifelse(key == "controlt1", "controlt1 (\u5bfe\u7167)", ifelse(key == "microt10", "microt10", ifelse(key == "microt100", "microt100", as.character(x))))
      )
    )
  )
  unname(out)
}

condition_family <- function(x) {
  key <- tolower(as.character(x))
  out <- ifelse(
    key %in% c("controlsea1", "beringsea", "okhotsksea"), "Sea",
    ifelse(key %in% c("controlt1", "microt10", "microt100"), "Intensity", "Other")
  )
  unname(out)
}

scalar_or_na <- function(x) {
  if (is.null(x) || length(x) == 0) return(NA_real_)
  suppressWarnings(as.numeric(x[[1]]))
}

format_deg_jp <- function(x, digits = 1, signed = FALSE) {
  x <- scalar_or_na(x)
  if (!is.finite(x)) return("NA")
  val <- formatC(x, digits = digits, format = "f")
  if (signed && x > 0) val <- paste0("+", val)
  paste0(val, "\u00b0")
}

format_ci_jp <- function(lo, hi, digits = 1, signed = TRUE) {
  paste0("[", format_deg_jp(lo, digits = digits, signed = signed), ", ", format_deg_jp(hi, digits = digits, signed = signed), "]")
}

pairwise_judgement <- function(lo, hi) {
  lo <- scalar_or_na(lo)
  hi <- scalar_or_na(hi)
  if (!is.finite(lo) || !is.finite(hi)) return("\u5224\u5b9a\u4e0d\u80fd")
  if (lo > 0 || hi < 0) return("\u6709\u610f\u5dee\u3042\u308a")
  if (abs(lo) <= 5 || abs(hi) <= 5) return("\u5883\u754c\u7684\uff08CI\u304c\u308f\u305a\u304b\u306b0\u3092\u542b\u3080\uff09")
  "\u6709\u610f\u5dee\u306a\u3057"
}

mixed_term_value <- function(mixed_terms, response, term, field = "p_value") {
  sub <- mixed_terms[mixed_terms$Response == response & mixed_terms$Term == term, , drop = FALSE]
  if (nrow(sub) == 0 || !(field %in% names(sub))) return(NA_real_)
  suppressWarnings(as.numeric(sub[[field]][1]))
}

window_label_jp <- function(lbl, start = NA_real_, end = NA_real_) {
  if (is.finite(scalar_or_na(start)) && is.finite(scalar_or_na(end))) {
    return(sprintf("%s (%s-%ss)", as.character(lbl), format_num(scalar_or_na(start), 0), format_num(scalar_or_na(end), 0)))
  }
  key <- tolower(as.character(lbl))
  out <- ifelse(
    key == "base", "base (300-600s)",
    ifelse(
      key == "early_shift", "early (240-540s)",
      ifelse(
        key == "late_shift", "late (360-660s)",
        ifelse(key == "narrow_core", "narrow (375-525s)", ifelse(key == "wide_context", "wide (240-660s)", as.character(lbl)))
      )
    )
  )
  unname(out)
}

pick_condition_bayes_terms <- function(bayes_out) {
  cand <- rbind(bayes_out$coef_means, bayes_out$coef_cat)
  if (!is.data.frame(cand) || nrow(cand) == 0 || !("Term" %in% names(cand))) return(data.frame())
  keep <- cand[grepl("^Condition", cand$Term), , drop = FALSE]
  if (nrow(keep) == 0) return(data.frame())
  keep
}

pretty_metric <- function(metric) {
  key <- tolower(as.character(metric))
  out <- ifelse(
    key == "speed", "Speed",
    ifelse(
      key == "distance10s", "Distance10s",
      ifelse(key == "turnrate", "TurnRate", ifelse(key == "curvature", "Curvature", ifelse(key == "tortuosity", "Tortuosity", as.character(metric))))
    )
  )
  unname(out)
}

build_results_interpretation <- function(agg, pairwise_ci, alignment, mixed_terms, bayes_out, cohort_pairwise, cohort_consistency, cohort_activity_pairwise, window_sensitivity, quality_sensitivity, within_fish, activity_pairwise, activity_mixed_terms, window_start, window_end) {
  conds <- unique(as.character(agg$Condition))
  conds_desc <- condition_description(conds)
  sample_n <- if (nrow(pairwise_ci) > 0 && "N" %in% names(pairwise_ci)) max(pairwise_ci$N, na.rm = TRUE) else sum(is.finite(agg$MeanHeadingDeg))
  lines <- c(
    "\u89e3\u6790\u6982\u8981",
    "\u89e3\u6790\u5bfe\u8c61: \u30b5\u30b1\u306e\u904a\u6cf3\u65b9\u4f4d\uff08heading\uff09\u306e\u5186\u74b0\u30c7\u30fc\u30bf",
    sprintf("\u6642\u9593\u7a93: %s-%s\u79d2\uff08sample\u671f\u9593\uff09", format_num(window_start, 0), format_num(window_end, 0)),
    sprintf("\u6761\u4ef6: %s", paste(conds_desc, collapse = "\u3001")),
    sprintf("\u30b5\u30f3\u30d7\u30eb\u30b5\u30a4\u30ba: \u5404\u6761\u4ef6 N=%d\uff08\u500b\u4f53\u00d7\u30b3\u30db\u30fc\u30c8\u306e\u7d44\u307f\u5408\u308f\u305b\uff09", sample_n),
    ""
  )

  lines <- c(lines, "1. \u30da\u30a2\u30ef\u30a4\u30ba\u5186\u74b0\u5e73\u5747\u5dee\uff08\u4e3b\u89e3\u6790\uff09")
  lines <- c(lines, "\u6bd4\u8f03 | \u5e73\u5747\u5dee | 95% bootstrap CI | Rbar | \u5224\u5b9a")
  near_zero_rows <- character()
  if (nrow(pairwise_ci) == 0) {
    lines <- c(lines, "\u30da\u30a2\u30ef\u30a4\u30ba\u6bd4\u8f03\u3092\u8a08\u7b97\u3067\u304d\u307e\u305b\u3093\u3067\u3057\u305f\u3002")
  } else {
    for (i in seq_len(nrow(pairwise_ci))) {
      row <- pairwise_ci[i, ]
      judge <- pairwise_judgement(row$CI_lo, row$CI_hi)
      lines <- c(lines, sprintf(
        "%s - %s | %s | %s | %s | %s",
        pretty_condition(row$GroupA),
        pretty_condition(row$GroupB),
        format_deg_jp(row$MeanDeltaDeg, digits = 1, signed = TRUE),
        format_ci_jp(row$CI_lo, row$CI_hi, digits = 1, signed = TRUE),
        format_num(row$RbarDelta, 3),
        judge
      ))
      if (is.finite(row$CI_lo) && is.finite(row$CI_hi) && row$CI_lo <= 0 && row$CI_hi >= 0 && min(abs(row$CI_lo), abs(row$CI_hi)) <= 5) {
        near_zero_rows <- c(near_zero_rows, sprintf(
          "%s vs %s \u306f CI \u304c 0 \u306b\u975e\u5e38\u306b\u8fd1\u304f\u3001%s \u306e\u30b7\u30d5\u30c8\u5099\u5019\u3068\u89e3\u91c8\u3067\u304d\u307e\u3059\u3002",
          pretty_condition(row$GroupA), pretty_condition(row$GroupB), format_deg_jp(row$MeanDeltaDeg, 1, TRUE)
        ))
      }
    }
    lines <- c(lines, "95% CI \u304c 0 \u3092\u8de8\u3050\u6bd4\u8f03\u306f\u3001\u5f93\u6765\u306e\u6709\u610f\u6c34\u6e96\u3067\u306f\u6709\u610f\u5dee\u306a\u3057\u3068\u89e3\u91c8\u3057\u307e\u3059\u3002")
    if (length(near_zero_rows) > 0) lines <- c(lines, near_zero_rows)
  }

  cos_f <- mixed_term_value(mixed_terms, "cosMag", "Condition", "F_value")
  cos_p <- mixed_term_value(mixed_terms, "cosMag", "Condition", "p_value")
  sin_f <- mixed_term_value(mixed_terms, "sinMag", "Condition", "F_value")
  sin_p <- mixed_term_value(mixed_terms, "sinMag", "Condition", "p_value")
  other_ps <- c(
    mixed_term_value(mixed_terms, "cosMag", "OrderInSession", "p_value"),
    mixed_term_value(mixed_terms, "sinMag", "OrderInSession", "p_value"),
    mixed_term_value(mixed_terms, "cosMag", "PrevCondition", "p_value"),
    mixed_term_value(mixed_terms, "sinMag", "PrevCondition", "p_value"),
    mixed_term_value(mixed_terms, "cosMag", "CohortGroup", "p_value"),
    mixed_term_value(mixed_terms, "sinMag", "CohortGroup", "p_value"),
    mixed_term_value(mixed_terms, "cosMag", "Condition:CohortGroup", "p_value"),
    mixed_term_value(mixed_terms, "sinMag", "Condition:CohortGroup", "p_value")
  )
  other_ps <- other_ps[is.finite(other_ps)]
  cos_p_adj <- mixed_term_value(mixed_terms, "cosMag", "Condition", "p_adj")
  sin_p_adj <- mixed_term_value(mixed_terms, "sinMag", "Condition", "p_adj")
  lines <- c(lines, "", "2. \u6df7\u5408\u30e2\u30c7\u30eb\uff08cos/sin\u6210\u5206\u306e\u7dda\u5f62\u6df7\u5408\u30e2\u30c7\u30eb\uff09", "\u5fdc\u7b54\u5909\u6570 | Condition \u306e F\u5024 | p\u5024 | p\u88dc\u6b63 | \u5224\u5b9a")
  lines <- c(lines, sprintf("cosMag | %s | %s | %s | %s", format_num(cos_f, 3), format_p(cos_p), format_p(cos_p_adj), ifelse(is.finite(cos_p_adj) && cos_p_adj < 0.05, "\u6709\u610f\u5dee\u3042\u308a", "\u6709\u610f\u5dee\u306a\u3057")))
  lines <- c(lines, sprintf("sinMag | %s | %s | %s | %s", format_num(sin_f, 3), format_p(sin_p), format_p(sin_p_adj), ifelse(is.finite(sin_p_adj) && sin_p_adj < 0.05, "\u6709\u610f\u5dee\u3042\u308a", "\u6709\u610f\u5dee\u306a\u3057")))
  if (length(other_ps) > 0) {
    lines <- c(lines, sprintf("OrderInSession\uff08\u9806\u5e8f\u52b9\u679c\uff09\u3001PrevCondition\uff08\u30ad\u30e3\u30ea\u30fc\u30aa\u30fc\u30d0\u30fc\u52b9\u679c\uff09\u3001CohortGroup\uff08\u30b3\u30db\u30fc\u30c8\u7fa4\u52b9\u679c\uff09\u3001\u4ea4\u4e92\u4f5c\u7528\u306f\u3059\u3079\u3066\u975e\u6709\u610f\u3067\u3057\u305f\uff08\u6700\u5c0f p = %s\uff09\u3002", format_p(min(other_ps, na.rm = TRUE))))
    lines <- c(lines, "\u78c1\u5834\u6761\u4ef6\u4ee5\u5916\u306e\u4ea4\u7d61\u8981\u56e0\u3082\u660e\u78ba\u306b\u306f\u691c\u51fa\u3055\u308c\u305a\u3001\u5b9f\u9a13\u30c7\u30b6\u30a4\u30f3\u7531\u6765\u306e\u5927\u304d\u306a\u30d0\u30a4\u30a2\u30b9\u306f\u793a\u3055\u308c\u307e\u305b\u3093\u3067\u3057\u305f\u3002")
  }

  lines <- c(lines, "", "2b. \u6d3b\u52d5\u91cf\u6307\u6a19\uff08Speed / Distance / TurnRate \u306a\u3069\uff09")
  if (!is.data.frame(activity_pairwise) || nrow(activity_pairwise) == 0) {
    lines <- c(lines, "\u6d3b\u52d5\u91cf\u6307\u6a19\u306e\u30da\u30a2\u30ef\u30a4\u30ba\u6bd4\u8f03\u306f\u8a08\u7b97\u3067\u304d\u307e\u305b\u3093\u3067\u3057\u305f\u3002")
  } else {
    lines <- c(lines, "\u6307\u6a19 | \u6bd4\u8f03 | \u5e73\u5747\u5dee | 95% CI | p\u5024 | p\u88dc\u6b63 | \u5224\u5b9a")
    lines <- c(lines, "\u5404\u884c\u52d5\u6307\u6a19\u306b\u3064\u3044\u3066 3 \u6761\u4ef6\u30da\u30a2\u9593\u306e\u5bfe\u5fdc\u3042\u308a t \u691c\u5b9a\u3092\u884c\u3044\u3001\u6307\u6a19\u3054\u3068\u306b Holm \u6cd5\u3067\u591a\u91cd\u6bd4\u8f03\u88dc\u6b63\u3057\u305f p \u5024\u3067\u5224\u5b9a\u3057\u3066\u3044\u307e\u3059\u3002")
    for (i in seq_len(nrow(activity_pairwise))) {
      row <- activity_pairwise[i, ]
      lines <- c(lines, sprintf(
        "%s | %s - %s | %s | [%s, %s] | %s | %s | %s",
        pretty_metric(row$Metric),
        pretty_condition(row$GroupA),
        pretty_condition(row$GroupB),
        format_num(row$MeanDiff, 3),
        format_num(row$CI_lo, 3),
        format_num(row$CI_hi, 3),
        format_p(row$p_value),
        format_p(row$p_adj),
        ifelse(is.finite(row$p_adj) && row$p_adj < 0.05, "\u6709\u610f\u5dee\u3042\u308a", "\u6709\u610f\u5dee\u306a\u3057")
      ))
    }
    sig_activity <- activity_pairwise[is.finite(activity_pairwise$p_adj) & activity_pairwise$p_adj < 0.05, , drop = FALSE]
    if (nrow(sig_activity) > 0) {
      lines <- c(lines, "\u6d3b\u52d5\u91cf\u3067\u306f\u3001\u4ee5\u4e0b\u306e\u30da\u30a2\u6bd4\u8f03\u3067\u6709\u610f\u5dee\u304c\u898b\u3089\u308c\u307e\u3057\u305f\u3002")
      for (i in seq_len(nrow(sig_activity))) {
        row <- sig_activity[i, ]
        direction_txt <- ifelse(is.finite(row$MeanDiff) && row$MeanDiff > 0,
          sprintf("%s \u306e\u65b9\u304c\u9ad8\u3044", pretty_condition(row$GroupA)),
          ifelse(is.finite(row$MeanDiff) && row$MeanDiff < 0,
            sprintf("%s \u306e\u65b9\u304c\u9ad8\u3044", pretty_condition(row$GroupB)),
            "\u65b9\u5411\u306f\u4e0d\u5b9a"
          )
        )
        lines <- c(lines, sprintf(
          "%s \u3067 %s - %s \u306b\u5dee\u304c\u3042\u308a\uff08\u5e73\u5747\u5dee %s, 95%% CI [%s, %s], p=%s\uff09\u3002%s\u3002",
          pretty_metric(row$Metric),
          pretty_condition(row$GroupA),
          pretty_condition(row$GroupB),
          format_num(row$MeanDiff, 3),
          format_num(row$CI_lo, 3),
          format_num(row$CI_hi, 3),
          format_p(row$p_adj),
          direction_txt
        ))
      }
    }
  }
  if (is.data.frame(activity_mixed_terms) && nrow(activity_mixed_terms) > 0) {
    lines <- c(lines, "\u6df7\u5408\u30e2\u30c7\u30eb\u306e Condition \u4e3b\u52b9\u679c")
    keep_metrics <- unique(as.character(activity_mixed_terms$Response))
    for (metric in keep_metrics) {
      sub <- activity_mixed_terms[activity_mixed_terms$Response == metric & activity_mixed_terms$Term == "Condition", , drop = FALSE]
      if (nrow(sub) == 0) next
      lines <- c(lines, sprintf(
        "%s: F=%s, p=%s, p\u88dc\u6b63=%s, %s",
        pretty_metric(metric),
        format_num(sub$F_value[1], 3),
        format_p(sub$p_value[1]),
        format_p(sub$p_adj[1]),
        ifelse(is.finite(sub$p_adj[1]) && sub$p_adj[1] < 0.05, "\u6761\u4ef6\u52b9\u679c\u3042\u308a", "\u6761\u4ef6\u52b9\u679c\u306f\u660e\u78ba\u3067\u306f\u3042\u308a\u307e\u305b\u3093")
      ))
    }
    sig_activity_model <- activity_mixed_terms[activity_mixed_terms$Term == "Condition" & is.finite(activity_mixed_terms$p_adj) & activity_mixed_terms$p_adj < 0.05, , drop = FALSE]
    if (nrow(sig_activity_model) == 0 && exists("sig_activity") && nrow(sig_activity) > 0) {
      lines <- c(lines, "\u305f\u3060\u3057\u3001\u3053\u308c\u3089\u306e activity \u306e\u6709\u610f\u5dee\u306f mixed model \u306e Condition \u4e3b\u52b9\u679c\u3068\u3057\u3066\u306f\u518d\u73fe\u3055\u308c\u3066\u304a\u3089\u305a\u3001\u9650\u5b9a\u7684\u306a\u30da\u30a2\u5dee\u3068\u3057\u3066\u89e3\u91c8\u3059\u308b\u306e\u304c\u9069\u5207\u3067\u3059\u3002")
    }
  }

  lines <- c(lines, "", "3. \u30d9\u30a4\u30ba\u5186\u74b0\u30e2\u30c7\u30eb\uff08bpnreg\uff09")
  cond_bayes <- pick_condition_bayes_terms(bayes_out)
  if (nrow(cond_bayes) == 0) {
    lines <- c(lines, "\u30d9\u30a4\u30ba\u5186\u74b0\u30e2\u30c7\u30eb\u306e\u4fc2\u6570\u8981\u7d04\u306f\u5229\u7528\u3067\u304d\u307e\u305b\u3093\u3067\u3057\u305f\u3002")
  } else {
    for (i in seq_len(nrow(cond_bayes))) {
      row <- cond_bayes[i, ]
      includes_zero <- is.finite(row$CI_lo) && is.finite(row$CI_hi) && row$CI_lo <= 0 && row$CI_hi >= 0
      lines <- c(lines, sprintf(
        "%s: mean=%s, 95%% HPD %s -> %s",
        row$Term,
        format_num(row$MeanValue, 3),
        format_ci_jp(row$CI_lo, row$CI_hi, digits = 2, signed = TRUE),
        ifelse(includes_zero, "0\u3092\u542b\u3080", "0\u3092\u542b\u307e\u306a\u3044")
      ))
    }
    lines <- c(lines, "\u30d9\u30a4\u30ba\u30e2\u30c7\u30eb\u3067\u3082\u6761\u4ef6\u52b9\u679c\u306e\u4e8b\u5f8c\u5206\u5e03\u304c\u660e\u78ba\u306b 0 \u304b\u3089\u96e2\u308c\u308b\u8a3c\u62e0\u306f\u5f37\u304f\u3042\u308a\u307e\u305b\u3093\u3002")
  }

  lines <- c(lines, "", "4. \u30b3\u30db\u30fc\u30c8\u65b9\u5411\u4e00\u81f4\u6027\uff08\u88dc\u52a9\u6240\u898b\uff09", "\u6bd4\u8f03 | \u540c\u65b9\u5411\u30b3\u30db\u30fc\u30c8\u7fa4 | \u4e00\u81f4\u7387 | \u7fa4\u9593Delta\u7bc4\u56f2")
  if (nrow(cohort_consistency) == 0) {
    lines <- c(lines, "\u30b3\u30db\u30fc\u30c8\u4e00\u81f4\u6027\u306f\u8a08\u7b97\u3067\u304d\u307e\u305b\u3093\u3067\u3057\u305f\u3002")
  } else {
    for (i in seq_len(nrow(cohort_consistency))) {
      row <- cohort_consistency[i, ]
      comps <- strsplit(as.character(row$Pair), " vs ", fixed = TRUE)[[1]]
      ga <- if (length(comps) >= 1) comps[1] else row$Pair
      gb <- if (length(comps) >= 2) comps[2] else ""
      lines <- c(lines, sprintf(
        "%s - %s | %d/%d | %s | %s \u301c %s",
        pretty_condition(ga),
        pretty_condition(gb),
        row$SameDirectionGroups,
        row$CohortGroups,
        paste0(format_num(row$SameDirectionFraction * 100, 0), "%"),
        format_deg_jp(row$MinGroupDeltaDeg, 1, TRUE),
        format_deg_jp(row$MaxGroupDeltaDeg, 1, TRUE)
      ))
    }
    strong_consistency <- cohort_consistency[cohort_consistency$SameDirectionGroups == cohort_consistency$CohortGroups, , drop = FALSE]
    if (nrow(strong_consistency) > 0) lines <- c(lines, "\u5c11\u306a\u304f\u3068\u3082\u4e00\u90e8\u306e\u6bd4\u8f03\u3067\u306f\u3001\u72ec\u7acb\u30b3\u30db\u30fc\u30c8\u7fa4\u3092\u307e\u305f\u3044\u3067\u52b9\u679c\u65b9\u5411\u304c\u4e00\u8cab\u3057\u3066\u3044\u307e\u3057\u305f\u3002")
  }
  sig_cohort_angle <- cohort_pairwise[isTRUE(cohort_pairwise$SignificantByCI), , drop = FALSE]
  if (nrow(sig_cohort_angle) > 0) {
    lines <- c(lines, "\u30b3\u30db\u30fc\u30c8 group \u5185\u306e\u89d2\u5ea6\u6bd4\u8f03\u3067\u6709\u610f\u5dee\u304c\u3042\u3063\u305f\u3082\u306e:")
    for (i in seq_len(nrow(sig_cohort_angle))) {
      row <- sig_cohort_angle[i, ]
      lines <- c(lines, sprintf(
        "%s | %s - %s | CI %s",
        row$CohortGroup,
        pretty_condition(row$GroupA),
        pretty_condition(row$GroupB),
        format_ci_jp(row$CI_lo, row$CI_hi, 1, TRUE)
      ))
    }
  } else {
    lines <- c(lines, "\u30b3\u30db\u30fc\u30c8 group \u5185\u306e\u89d2\u5ea6\u6bd4\u8f03\u3067\u306f\u3001CI \u57fa\u6e96\u3067\u6709\u610f\u5dee\u306f\u78ba\u8a8d\u3055\u308c\u307e\u305b\u3093\u3067\u3057\u305f\u3002")
  }
  if (is.data.frame(cohort_activity_pairwise) && nrow(cohort_activity_pairwise) > 0) {
    sig_cohort_activity <- cohort_activity_pairwise[isTRUE(cohort_activity_pairwise$SignificantByP), , drop = FALSE]
    if (nrow(sig_cohort_activity) > 0) {
      lines <- c(lines, "\u30b3\u30db\u30fc\u30c8 group \u5185\u306e\u6d3b\u52d5\u91cf\u6bd4\u8f03\u3067\u88dc\u6b63\u5f8c\u6709\u610f\u3060\u3063\u305f\u3082\u306e:")
      for (i in seq_len(nrow(sig_cohort_activity))) {
        row <- sig_cohort_activity[i, ]
        lines <- c(lines, sprintf(
          "%s | %s | %s - %s | p\u88dc\u6b63=%s",
          row$CohortGroup,
          pretty_metric(row$Metric),
          pretty_condition(row$GroupA),
          pretty_condition(row$GroupB),
          format_p(row$p_adj)
        ))
      }
    } else {
      lines <- c(lines, "\u30b3\u30db\u30fc\u30c8 group \u5185\u306e\u6d3b\u52d5\u91cf\u6bd4\u8f03\u3067\u306f\u3001\u88dc\u6b63\u5f8c\u306b\u6709\u610f\u5dee\u306f\u78ba\u8a8d\u3055\u308c\u307e\u305b\u3093\u3067\u3057\u305f\u3002")
    }
  }

  lines <- c(lines, "", "5. \u500b\u4f53\u5185\u4e00\u8cab\u6027\uff08Within-fish analysis\uff09", "\u53c2\u7167\u6761\u4ef6 -> \u5b9f\u9a13\u6761\u4ef6 | \u5186\u74b0\u76f8\u95a2 r | \u9806\u5217\u691c\u5b9a p | \u5224\u5b9a")
  if (nrow(within_fish$reference_correlations) == 0) {
    lines <- c(lines, "\u500b\u4f53\u5185\u4e00\u8cab\u6027\u306f\u8a08\u7b97\u3067\u304d\u307e\u305b\u3093\u3067\u3057\u305f\u3002")
  } else {
    for (i in seq_len(nrow(within_fish$reference_correlations))) {
      row <- within_fish$reference_correlations[i, ]
      lines <- c(lines, sprintf(
        "%s -> %s | %s | p = %s | %s",
        pretty_condition(row$ReferenceCondition),
        pretty_condition(row$Condition),
        format_num(row$CircularCorrelation, 3),
        format_p(row$PermutationP),
        ifelse(is.finite(row$PermutationP) && row$PermutationP < 0.05, "\u6709\u610f", "\u6709\u610f\u5dee\u306a\u3057")
      ))
    }
    sig_corr <- within_fish$reference_correlations[is.finite(within_fish$reference_correlations$PermutationP) & within_fish$reference_correlations$PermutationP < 0.05, , drop = FALSE]
    if (nrow(sig_corr) > 0) {
      lines <- c(lines, "\u500b\u4f53\u3054\u3068\u306e\u65b9\u5411\u9078\u597d\u306f\u6761\u4ef6\u3092\u8de8\u3044\u3067\u90e8\u5206\u7684\u306b\u4fdd\u5b58\u3055\u308c\u3066\u304a\u308a\u3001\u500b\u4f53\u5185\u76f8\u95a2\u304c\u4eca\u56de\u306e\u89e3\u6790\u3067\u6700\u3082\u660e\u77ad\u306a\u30b7\u30b0\u30ca\u30eb\u3067\u3059\u3002")
    }
  }
  if (nrow(within_fish$shift_summary) > 0) {
    lines <- c(lines, "ReferenceCondition -> Condition | MeanShift | MedianShift | 95% bootstrap CI | sign test p | permutation p | PositiveFraction")
    for (i in seq_len(nrow(within_fish$shift_summary))) {
      row <- within_fish$shift_summary[i, ]
      lines <- c(lines, sprintf(
        "%s -> %s | %s | %s | %s | %s | %s | %s",
        pretty_condition(row$ReferenceCondition),
        pretty_condition(row$Condition),
        format_deg_jp(row$MeanShiftDeg, 1, TRUE),
        format_deg_jp(row$MedianShiftDeg, 1, TRUE),
        format_ci_jp(row$CI_lo, row$CI_hi, 1, TRUE),
        format_p(row$SignTestP),
        format_p(row$PermutationP),
        paste0(format_num(row$PositiveFraction * 100, 0), "%")
      ))
    }
    non_sig_shift <- within_fish$shift_summary[is.finite(within_fish$shift_summary$SignTestP), , drop = FALSE]
    if (nrow(non_sig_shift) > 0) {
      lines <- c(lines, sprintf("\u4e00\u65b9\u3067 shift \u306e\u7b26\u53f7\u691c\u5b9a\u306f %s \u301c %s \u3067\u3001\u5168\u500b\u4f53\u304c\u540c\u65b9\u5411\u306b\u63c3\u3063\u3066\u30b7\u30d5\u30c8\u3059\u308b\u96c6\u56e3\u5fdc\u7b54\u306f\u652f\u6301\u3055\u308c\u307e\u305b\u3093\u3067\u3057\u305f\u3002",
        format_p(min(non_sig_shift$SignTestP, na.rm = TRUE)),
        format_p(max(non_sig_shift$SignTestP, na.rm = TRUE))
      ))
    }
    sig_perm_shift <- within_fish$shift_summary[is.finite(within_fish$shift_summary$PermutationP) & within_fish$shift_summary$PermutationP < 0.05, , drop = FALSE]
    if (nrow(sig_perm_shift) > 0) {
      for (i in seq_len(nrow(sig_perm_shift))) {
        row <- sig_perm_shift[i, ]
        lines <- c(lines, sprintf(
          "%s -> %s \u3067\u306f\u3001\u9b5a\u3054\u3068\u306e ShiftDeg \u5206\u5e03\u306e\u5186\u74b0\u5e73\u5747\u306f %s\u3001\u4e2d\u592e\u5024\u306f %s\u300195%% bootstrap CI \u306f %s \u3067\u3057\u305f\u3002permutation test \u3067\u306f p=%s \u3067\u30010 \u304b\u3089\u306e\u7cfb\u7d71\u7684\u306a\u56de\u8ee2\u30b7\u30d5\u30c8\u304c\u793a\u5506\u3055\u308c\u307e\u3059\u3002",
          pretty_condition(row$ReferenceCondition),
          pretty_condition(row$Condition),
          format_deg_jp(row$MeanShiftDeg, 1, TRUE),
          format_deg_jp(row$MedianShiftDeg, 1, TRUE),
          format_ci_jp(row$CI_lo, row$CI_hi, 1, TRUE),
          format_p(row$PermutationP)
        ))
      }
    }
  }
  if (is.data.frame(within_fish$shift_cohort_summary) && nrow(within_fish$shift_cohort_summary) > 0) {
    lines <- c(lines, "CohortGroup-specific ShiftDeg summary", "CohortGroup | ReferenceCondition -> Condition | MeanShift | MedianShift | 95% bootstrap CI | sign test p | permutation p | N")
    for (i in seq_len(nrow(within_fish$shift_cohort_summary))) {
      row <- within_fish$shift_cohort_summary[i, ]
      lines <- c(lines, sprintf(
        "%s | %s -> %s | %s | %s | %s | %s | %s | %d",
        row$CohortGroup,
        pretty_condition(row$ReferenceCondition),
        pretty_condition(row$Condition),
        format_deg_jp(row$MeanShiftDeg, 1, TRUE),
        format_deg_jp(row$MedianShiftDeg, 1, TRUE),
        format_ci_jp(row$CI_lo, row$CI_hi, 1, TRUE),
        format_p(row$SignTestP),
        format_p(row$PermutationP),
        row$N
      ))
    }
    sig_cohort_shift <- within_fish$shift_cohort_summary[is.finite(within_fish$shift_cohort_summary$PermutationP) & within_fish$shift_cohort_summary$PermutationP < 0.05, , drop = FALSE]
    if (nrow(sig_cohort_shift) > 0) {
      lines <- c(lines, "Some CohortGroups showed a non-zero bias in ShiftDeg.")
    } else {
      lines <- c(lines, "ShiftDeg bias was not reproduced significantly within individual CohortGroups.")
    }
  }

  lines <- c(lines, "", "6. \u6642\u9593\u7a93\u611f\u5ea6\u5206\u6790")
  if (nrow(window_sensitivity$pairwise) == 0) {
    lines <- c(lines, "\u6642\u9593\u7a93\u611f\u5ea6\u5206\u6790\u306f\u8a08\u7b97\u3067\u304d\u307e\u305b\u3093\u3067\u3057\u305f\u3002")
  } else {
    labels <- unique(as.character(window_sensitivity$pairwise$WindowLabel))
    label_display <- sapply(labels, function(lbl) {
      sub <- window_sensitivity$pairwise[as.character(window_sensitivity$pairwise$WindowLabel) == lbl, , drop = FALSE]
      window_label_jp(lbl, sub$WindowStart[1], sub$WindowEnd[1])
    })
    header <- paste(c("\u6bd4\u8f03", label_display), collapse = " | ")
    lines <- c(lines, header)
    pair_keys <- unique(paste(window_sensitivity$pairwise$GroupA, window_sensitivity$pairwise$GroupB, sep = " | "))
    for (pk in pair_keys) {
      sub <- window_sensitivity$pairwise[paste(window_sensitivity$pairwise$GroupA, window_sensitivity$pairwise$GroupB, sep = " | ") == pk, , drop = FALSE]
      sub <- sub[match(labels, as.character(sub$WindowLabel)), , drop = FALSE]
      comp_label <- sprintf("%s - %s", pretty_condition(sub$GroupA[1]), pretty_condition(sub$GroupB[1]))
      vals <- sapply(seq_len(nrow(sub)), function(i) format_deg_jp(sub$MeanDeltaDeg[i], 1, TRUE))
      lines <- c(lines, paste(c(comp_label, vals), collapse = " | "))
    }
    narrow_hit <- window_sensitivity$pairwise[tolower(as.character(window_sensitivity$pairwise$WindowLabel)) == "narrow_core" &
      ((tolower(as.character(window_sensitivity$pairwise$GroupA)) == "controlsea1" & tolower(as.character(window_sensitivity$pairwise$GroupB)) == "okhotsksea") |
       (tolower(as.character(window_sensitivity$pairwise$GroupA)) == "okhotsksea" & tolower(as.character(window_sensitivity$pairwise$GroupB)) == "controlsea1")), , drop = FALSE]
    if (nrow(narrow_hit) > 0) {
      row <- narrow_hit[1, ]
      lines <- c(lines, sprintf("narrow_core \u7a93\u3067\u306f %s - %s \u306e CI \u304c %s \u3067\u30010 \u3092\u308f\u305a\u304b\u306b\u5916\u308c\u308b\u53ef\u80fd\u6027\u304c\u3042\u308a\u307e\u3059\u3002\u591a\u91cd\u6bd4\u8f03\u3092\u8003\u3048\u308b\u3068\u614e\u91cd\u306a\u89e3\u91c8\u304c\u5fc5\u8981\u3067\u3059\u3002",
        pretty_condition(row$GroupA), pretty_condition(row$GroupB), format_ci_jp(row$CI_lo, row$CI_hi, 2, TRUE)
      ))
    }
    lines <- c(lines, "\u5168\u4f53\u3068\u3057\u3066\u306f\u3001\u52b9\u679c\u65b9\u5411\u306f\u6642\u9593\u7a93\u3092\u5909\u3048\u3066\u3082\u5927\u304d\u304f\u306f\u53cd\u8ee2\u305b\u305a\u3001\u6975\u7aef\u306a\u30a6\u30a3\u30f3\u30c9\u30a6\u4f9d\u5b58\u6027\u306f\u5f37\u304f\u3042\u308a\u307e\u305b\u3093\u3002")
  }

  lines <- c(lines, "", "\u7dcf\u5408\u89e3\u91c8", "\u7d50\u8ad6")
  all_cross_zero <- nrow(pairwise_ci) > 0 && all(pairwise_ci$CI_lo <= 0 & pairwise_ci$CI_hi >= 0, na.rm = TRUE)
  if (all_cross_zero) {
    lines <- c(lines, "\u96c6\u56e3\u30ec\u30d9\u30eb\u3067\u306e\u65b9\u4f4d\u306e\u6709\u610f\u306a\u6761\u4ef6\u5dee\u306f\u660e\u78ba\u306b\u306f\u691c\u51fa\u3055\u308c\u307e\u305b\u3093\u3067\u3057\u305f\u3002")
  } else {
    lines <- c(lines, "\u4e00\u90e8\u306e\u6bd4\u8f03\u3067\u306f\u96c6\u56e3\u30ec\u30d9\u30eb\u306e\u6761\u4ef6\u5dee\u304c\u793a\u5506\u3055\u308c\u307e\u3057\u305f\u304c\u3001\u5168\u4f53\u3068\u3057\u3066\u306f\u5f37\u3044\u8a3c\u62e0\u3067\u306f\u3042\u308a\u307e\u305b\u3093\u3002")
  }
  lines <- c(lines, "", "\u305f\u3060\u3057\u3001\u5f31\u3044\u304c\u4e00\u8cab\u3057\u305f\u30b7\u30b0\u30ca\u30eb\u3068\u3057\u3066\u4ee5\u4e0b\u304c\u6319\u3052\u3089\u308c\u307e\u3059\u3002")
  if (nrow(pairwise_ci) > 0) {
    ord_pw <- pairwise_ci[order(abs(pairwise_ci$MeanDeltaDeg), decreasing = TRUE), , drop = FALSE]
    top_n <- min(3, nrow(ord_pw))
    for (i in seq_len(top_n)) {
      row <- ord_pw[i, ]
      lines <- c(lines, sprintf("%s - %s: \u5e73\u5747\u30b7\u30d5\u30c8 %s\u3001CI %s\u3002",
        pretty_condition(row$GroupA), pretty_condition(row$GroupB), format_deg_jp(row$MeanDeltaDeg, 1, TRUE), format_ci_jp(row$CI_lo, row$CI_hi, 1, TRUE)
      ))
    }
  }
  if (is.data.frame(activity_pairwise) && nrow(activity_pairwise) > 0) {
    sig_activity <- activity_pairwise[is.finite(activity_pairwise$p_adj) & activity_pairwise$p_adj < 0.05, , drop = FALSE]
    if (nrow(sig_activity) > 0) {
      lines <- c(lines, "\u6d3b\u52d5\u91cf\u6307\u6a19\u3067\u306f\u4e00\u90e8\u306b\u6709\u610f\u306a\u30da\u30a2\u5dee\u304c\u3042\u308a\u3001\u3068\u304f\u306b Distance10s \u3084 Curvature \u306b\u6761\u4ef6\u9593\u5dee\u304c\u793a\u3055\u308c\u307e\u3057\u305f\u3002")
      if (is.data.frame(activity_mixed_terms) && nrow(activity_mixed_terms) > 0) {
        sig_activity_model <- activity_mixed_terms[activity_mixed_terms$Term == "Condition" & is.finite(activity_mixed_terms$p_adj) & activity_mixed_terms$p_adj < 0.05, , drop = FALSE]
        if (nrow(sig_activity_model) == 0) {
          lines <- c(lines, "\u305f\u3060\u3057\u3001activity \u6307\u6a19\u306e mixed model \u3067\u306f Condition \u4e3b\u52b9\u679c\u306f\u6709\u610f\u3067\u306f\u306a\u304f\u3001\u7d50\u679c\u306f\u9650\u5b9a\u7684\u306a\u30da\u30a2\u6bd4\u8f03\u3068\u3057\u3066\u614e\u91cd\u306b\u89e3\u91c8\u3059\u3079\u304d\u3067\u3059\u3002")
        }
      }
    }
  }
  if (nrow(cohort_consistency) > 0) lines <- c(lines, "\u4e00\u90e8\u306e\u6bd4\u8f03\u3067\u306f\u30b3\u30db\u30fc\u30c8\u7fa4\u9593\u3067\u52b9\u679c\u65b9\u5411\u304c\u9ad8\u3044\u4e00\u81f4\u7387\u3092\u793a\u3057\u307e\u3057\u305f\u3002")
  if (nrow(within_fish$reference_correlations) > 0 && any(within_fish$reference_correlations$PermutationP < 0.05, na.rm = TRUE)) {
    lines <- c(lines, "\u500b\u4f53\u5185\u306e\u65b9\u4f4d\u4fdd\u5b58\u6027\u306f\u6709\u610f\u3067\u3001\u78c1\u5834\u6761\u4ef6\u4e0b\u3067\u3082\u500b\u4f53\u3054\u3068\u306e\u65b9\u5411\u9078\u597d\u304c\u7dad\u6301\u3055\u308c\u308b\u53ef\u80fd\u6027\u304c\u3042\u308a\u307e\u3059\u3002")
  }
  if (length(other_ps) > 0) lines <- c(lines, "\u9806\u5e8f\u52b9\u679c\u30fb\u30ad\u30e3\u30ea\u30fc\u30aa\u30fc\u30d0\u30fc\u30fb\u30b3\u30db\u30fc\u30c8\u7fa4\u5dee\u306f\u5927\u304d\u304f\u306a\u304f\u3001\u89e3\u6790\u7d50\u679c\u3092\u5927\u304d\u304f\u6b6a\u3081\u308b\u4ea4\u7d61\u306e\u8a3c\u62e0\u306f\u4e4f\u3057\u3044\u3068\u8003\u3048\u3089\u308c\u307e\u3059\u3002")
  lines
}

build_results_interpretation_v2 <- function(agg, pairwise_ci, alignment, mixed_terms, bayes_out, cohort_pairwise, cohort_consistency, cohort_activity_pairwise, age_pairwise, age_activity_pairwise, age_mixed_terms, order_activity_pairwise, order_mixed_terms, subset_order_activity, window_sensitivity, quality_sensitivity, within_fish, shift_clusters, activity_pairwise, activity_mixed_terms, window_start, window_end) {
  conds <- unique(as.character(agg$Condition))
  conds_desc <- condition_description(conds)
  sample_n <- if (nrow(pairwise_ci) > 0 && "N" %in% names(pairwise_ci)) max(pairwise_ci$N, na.rm = TRUE) else sum(is.finite(agg$MeanHeadingDeg))
  lines <- c(
    "Analysis overview",
    "Target: circular heading data of salmon swimming direction",
    sprintf("Window: %s-%s s (sample period)", format_num(window_start, 0), format_num(window_end, 0)),
    sprintf("Conditions: %s", paste(conds_desc, collapse = ", ")),
    sprintf("Sample size: N=%d per condition (fish x cohort combinations)", sample_n),
    ""
  )
  lines <- c(lines, "1. Pairwise circular mean differences", "Comparison | Mean difference | 95% bootstrap CI | Rbar | Judgement")
  if (nrow(pairwise_ci) == 0) {
    lines <- c(lines, "Pairwise comparisons could not be computed.")
  } else {
    for (i in seq_len(nrow(pairwise_ci))) {
      row <- pairwise_ci[i, ]
      lines <- c(lines, sprintf("%s - %s | %s | %s | %s | %s",
        pretty_condition(row$GroupA), pretty_condition(row$GroupB),
        format_deg_jp(row$MeanDeltaDeg, 1, TRUE),
        format_ci_jp(row$CI_lo, row$CI_hi, 1, TRUE),
        format_num(row$RbarDelta, 3),
        pairwise_judgement(row$CI_lo, row$CI_hi)
      ))
    }
  }

  cos_f <- mixed_term_value(mixed_terms, "cosMag", "Condition", "F_value")
  cos_p <- mixed_term_value(mixed_terms, "cosMag", "Condition", "p_value")
  sin_f <- mixed_term_value(mixed_terms, "sinMag", "Condition", "F_value")
  sin_p <- mixed_term_value(mixed_terms, "sinMag", "Condition", "p_value")
  cos_p_adj <- mixed_term_value(mixed_terms, "cosMag", "Condition", "p_adj")
  sin_p_adj <- mixed_term_value(mixed_terms, "sinMag", "Condition", "p_adj")
  lines <- c(lines, "", "2. Mixed models for circular components", "Response | F for Condition | p | p_adj | judgement")
  lines <- c(lines, sprintf("cosMag | %s | %s | %s | %s", format_num(cos_f, 3), format_p(cos_p), format_p(cos_p_adj), ifelse(is.finite(cos_p_adj) && cos_p_adj < 0.05, "significant", "not significant")))
  lines <- c(lines, sprintf("sinMag | %s | %s | %s | %s", format_num(sin_f, 3), format_p(sin_p), format_p(sin_p_adj), ifelse(is.finite(sin_p_adj) && sin_p_adj < 0.05, "significant", "not significant")))

  if (is.data.frame(activity_pairwise) && nrow(activity_pairwise) > 0) {
    lines <- c(lines, "", "2b. Activity metrics", "Metric | Comparison | MeanDiff | 95% CI | p | p_adj | judgement")
    for (i in seq_len(nrow(activity_pairwise))) {
      row <- activity_pairwise[i, ]
      lines <- c(lines, sprintf("%s | %s - %s | %s | [%s, %s] | %s | %s | %s",
        pretty_metric(row$Metric),
        pretty_condition(row$GroupA), pretty_condition(row$GroupB),
        format_num(row$MeanDiff, 3),
        format_num(row$CI_lo, 3), format_num(row$CI_hi, 3),
        format_p(row$p_value), format_p(row$p_adj),
        ifelse(is.finite(row$p_adj) && row$p_adj < 0.05, "significant", "not significant")
      ))
    }
  }
  if (is.data.frame(order_activity_pairwise) && nrow(order_activity_pairwise) > 0) {
    lines <- c(lines, "", "2c. OrderInSession analyses", "Metric | Order comparison | MeanDiff | 95% CI | p | p_adj | judgement")
    lines <- c(lines, "This section tests whether activity differs by 1st / 2nd / 3rd session regardless of condition.")
    for (i in seq_len(nrow(order_activity_pairwise))) {
      row <- order_activity_pairwise[i, ]
      lines <- c(lines, sprintf("%s | %s vs %s | %s | [%s, %s] | %s | %s | %s",
        pretty_metric(row$Metric),
        row$GroupA, row$GroupB,
        format_num(row$MeanDiff, 3),
        format_num(row$CI_lo, 3), format_num(row$CI_hi, 3),
        format_p(row$p_value), format_p(row$p_adj),
        ifelse(is.finite(row$p_adj) && row$p_adj < 0.05, "significant", "not significant")
      ))
    }
    if (is.data.frame(order_mixed_terms) && nrow(order_mixed_terms) > 0) {
      lines <- c(lines, "Order-focused mixed models", "Response | F for OrderInSession | p | judgement")
      for (metric in unique(as.character(order_mixed_terms$Response))) {
        sub <- order_mixed_terms[order_mixed_terms$Response == metric & order_mixed_terms$Term == "OrderInSession", , drop = FALSE]
        if (nrow(sub) == 0) next
        lines <- c(lines, sprintf("%s | %s | %s | %s",
          pretty_metric(metric),
          format_num(sub$F_value[1], 3),
          format_p(sub$p_value[1]),
          ifelse(is.finite(sub$p_value[1]) && sub$p_value[1] < 0.05, "order effect suggested", "no clear order effect")
        ))
      }
    }
    sig_order <- order_activity_pairwise[is.finite(order_activity_pairwise$p_adj) & order_activity_pairwise$p_adj < 0.05, , drop = FALSE]
    if (nrow(sig_order) > 0) {
      lines <- c(lines, "If session1 is consistently lower than session2/session3 across metrics, that supports an order effect rather than a magnetic-condition effect.")
    } else {
      lines <- c(lines, "No robust activity difference was detected across 1st / 2nd / 3rd session order after correction.")
    }
  }
  if (is.list(subset_order_activity) && is.data.frame(subset_order_activity$model_terms) && nrow(subset_order_activity$model_terms) > 0) {
    lines <- c(lines, "", "2d. Subset-specific OrderInSession models", "ConditionFamily | AgeGroup | Response | F for OrderInSession | p | judgement")
    term_tbl <- subset_order_activity$model_terms
    term_tbl <- term_tbl[term_tbl$Term == "OrderInSession", , drop = FALSE]
    for (i in seq_len(nrow(term_tbl))) {
      row <- term_tbl[i, ]
      lines <- c(lines, sprintf("%s | %s | %s | %s | %s | %s",
        row$ConditionFamily,
        row$AgeGroup,
        pretty_metric(row$Response),
        format_num(row$F_value, 3),
        format_p(row$p_value),
        ifelse(is.finite(row$p_value) && row$p_value < 0.05, "order effect suggested", "not significant")
      ))
    }
    lines <- c(lines, "These subset models test order effects within Sea-only and Intensity-only datasets while treating Condition as a random intercept, so they are more direct tests of session-order artifacts.")
  }

  cond_bayes <- pick_condition_bayes_terms(bayes_out)
  lines <- c(lines, "", "3. Bayesian circular model")
  if (nrow(cond_bayes) == 0) {
    lines <- c(lines, "Bayesian coefficient summary was not available.")
  } else {
    for (i in seq_len(nrow(cond_bayes))) {
      row <- cond_bayes[i, ]
      lines <- c(lines, sprintf("%s | mean=%s | 95%% HPD %s", row$Term, format_num(row$MeanValue, 3), format_ci_jp(row$CI_lo, row$CI_hi, 2, TRUE)))
    }
  }

  lines <- c(lines, "", "4. Cohort consistency", "Comparison | same-direction groups | agreement | delta range")
  if (is.data.frame(cohort_consistency) && nrow(cohort_consistency) > 0) {
    for (i in seq_len(nrow(cohort_consistency))) {
      row <- cohort_consistency[i, ]
      comps <- strsplit(as.character(row$Pair), " vs ", fixed = TRUE)[[1]]
      ga <- if (length(comps) >= 1) comps[1] else row$Pair
      gb <- if (length(comps) >= 2) comps[2] else ""
      lines <- c(lines, sprintf("%s - %s | %d/%d | %s | %s to %s",
        pretty_condition(ga), pretty_condition(gb),
        row$SameDirectionGroups, row$CohortGroups,
        paste0(format_num(row$SameDirectionFraction * 100, 0), "%"),
        format_deg_jp(row$MinGroupDeltaDeg, 1, TRUE),
        format_deg_jp(row$MaxGroupDeltaDeg, 1, TRUE)
      ))
    }
  }
  if (is.data.frame(cohort_activity_pairwise) && nrow(cohort_activity_pairwise) > 0) {
    sig_cohort_activity <- cohort_activity_pairwise[isTRUE(cohort_activity_pairwise$SignificantByP), , drop = FALSE]
    if (nrow(sig_cohort_activity) == 0) lines <- c(lines, "No cohort-specific activity comparison survived adjusted significance.")
  }

  if ("AgeGroup" %in% names(agg) && any(!(is.na(agg$AgeGroup) | as.character(agg$AgeGroup) == ""))) {
    lines <- c(lines, "", "4b. AgeGroup analyses")
    age_levels <- unique(as.character(agg$AgeGroup))
    age_levels <- age_levels[!(is.na(age_levels) | age_levels == "")]
    if (length(age_levels) > 0) lines <- c(lines, sprintf("Age groups detected: %s", paste(age_levels, collapse = ", ")))
    if (is.data.frame(age_pairwise) && nrow(age_pairwise) > 0) {
      lines <- c(lines, "AgeGroup | Comparison | Mean difference | 95% bootstrap CI | Rbar | Judgement")
      for (i in seq_len(nrow(age_pairwise))) {
        row <- age_pairwise[i, ]
        lines <- c(lines, sprintf("%s | %s - %s | %s | %s | %s | %s",
          row$AgeGroup,
          pretty_condition(row$GroupA), pretty_condition(row$GroupB),
          format_deg_jp(row$MeanDeltaDeg, 1, TRUE),
          format_ci_jp(row$CI_lo, row$CI_hi, 1, TRUE),
          format_num(row$RbarDelta, 3),
          pairwise_judgement(row$CI_lo, row$CI_hi)
        ))
      }
    }
    if (is.data.frame(age_activity_pairwise) && nrow(age_activity_pairwise) > 0) {
      lines <- c(lines, "AgeGroup-specific activity tests", "AgeGroup | Metric | Comparison | MeanDiff | 95% CI | p_adj | judgement")
      for (i in seq_len(nrow(age_activity_pairwise))) {
        row <- age_activity_pairwise[i, ]
        lines <- c(lines, sprintf("%s | %s | %s - %s | %s | [%s, %s] | %s | %s",
          row$AgeGroup,
          pretty_metric(row$Metric),
          pretty_condition(row$GroupA), pretty_condition(row$GroupB),
          format_num(row$MeanDiff, 3),
          format_num(row$CI_lo, 3), format_num(row$CI_hi, 3),
          format_p(row$p_adj),
          ifelse(is.finite(row$p_adj) && row$p_adj < 0.05, "significant", "not significant")
        ))
      }
    }
    if (is.data.frame(age_mixed_terms) && nrow(age_mixed_terms) > 0) {
      lines <- c(lines, "AgeGroup interaction mixed models", "Response | Term | F | p")
      for (i in seq_len(nrow(age_mixed_terms))) {
        row <- age_mixed_terms[i, ]
        lines <- c(lines, sprintf("%s | %s | %s | %s",
          pretty_metric(row$Response),
          row$Term,
          format_num(row$F_value, 3),
          format_p(row$p_value)
        ))
      }
    }
  }

  lines <- c(lines, "", "5. Within-fish consistency", "ReferenceCondition -> Condition | Circular correlation r | permutation p | judgement")
  if (is.data.frame(within_fish$reference_correlations) && nrow(within_fish$reference_correlations) > 0) {
    for (i in seq_len(nrow(within_fish$reference_correlations))) {
      row <- within_fish$reference_correlations[i, ]
      lines <- c(lines, sprintf("%s -> %s | %s | p = %s | %s",
        pretty_condition(row$ReferenceCondition), pretty_condition(row$Condition),
        format_num(row$CircularCorrelation, 3), format_p(row$PermutationP),
        ifelse(is.finite(row$PermutationP) && row$PermutationP < 0.05, "significant", "not significant")
      ))
    }
  }
  if (is.data.frame(within_fish$shift_summary) && nrow(within_fish$shift_summary) > 0) {
    lines <- c(lines, "ReferenceCondition -> Condition | MeanShift | MedianShift | 95% bootstrap CI | sign test p | permutation p | PositiveFraction")
    for (i in seq_len(nrow(within_fish$shift_summary))) {
      row <- within_fish$shift_summary[i, ]
      lines <- c(lines, sprintf("%s -> %s | %s | %s | %s | %s | %s | %s",
        pretty_condition(row$ReferenceCondition), pretty_condition(row$Condition),
        format_deg_jp(row$MeanShiftDeg, 1, TRUE), format_deg_jp(row$MedianShiftDeg, 1, TRUE),
        format_ci_jp(row$CI_lo, row$CI_hi, 1, TRUE),
        format_p(row$SignTestP), format_p(row$PermutationP),
        paste0(format_num(row$PositiveFraction * 100, 0), "%")
      ))
    }
  }
  if (is.data.frame(within_fish$shift_cohort_summary) && nrow(within_fish$shift_cohort_summary) > 0) {
    lines <- c(lines, "CohortGroup-specific ShiftDeg summary", "CohortGroup | ReferenceCondition -> Condition | MeanShift | MedianShift | 95% bootstrap CI | sign test p | permutation p | N")
    for (i in seq_len(nrow(within_fish$shift_cohort_summary))) {
      row <- within_fish$shift_cohort_summary[i, ]
      lines <- c(lines, sprintf("%s | %s -> %s | %s | %s | %s | %s | %s | %d",
        row$CohortGroup, pretty_condition(row$ReferenceCondition), pretty_condition(row$Condition),
        format_deg_jp(row$MeanShiftDeg, 1, TRUE), format_deg_jp(row$MedianShiftDeg, 1, TRUE),
        format_ci_jp(row$CI_lo, row$CI_hi, 1, TRUE), format_p(row$SignTestP), format_p(row$PermutationP), row$N
      ))
    }
  }
  if (is.data.frame(within_fish$shift_pair_cohort_correlations) && nrow(within_fish$shift_pair_cohort_correlations) > 0) {
    lines <- c(lines, "CohortGroup-specific shift correlations", "CohortGroup | ConditionA vs ConditionB | Pearson r | permutation p | Spearman rho | permutation p | same-sign fraction | N")
    for (i in seq_len(nrow(within_fish$shift_pair_cohort_correlations))) {
      row <- within_fish$shift_pair_cohort_correlations[i, ]
      lines <- c(lines, sprintf("%s | %s vs %s | %s | %s | %s | %s | %s | %d",
        row$CohortGroup,
        pretty_condition(row$ConditionA), pretty_condition(row$ConditionB),
        format_num(row$PearsonR, 3), format_p(row$PearsonP),
        format_num(row$SpearmanR, 3), format_p(row$SpearmanP),
        paste0(format_num(row$SameSignFraction * 100, 0), "%"),
        row$N
      ))
    }
    lines <- c(lines, "Different CohortGroups can show different correlation signs because the relation between beringsea and okhotsksea shifts may depend on presentation sequence, random cohort composition, and limited per-group sample size.")
  }
  if (is.data.frame(within_fish$shift_age_summary) && nrow(within_fish$shift_age_summary) > 0) {
    lines <- c(lines, "AgeGroup-specific ShiftDeg summary", "AgeGroup | ReferenceCondition -> Condition | MeanShift | MedianShift | 95% bootstrap CI | sign test p | permutation p | N")
    for (i in seq_len(nrow(within_fish$shift_age_summary))) {
      row <- within_fish$shift_age_summary[i, ]
      lines <- c(lines, sprintf("%s | %s -> %s | %s | %s | %s | %s | %s | %d",
        row$AgeGroup, pretty_condition(row$ReferenceCondition), pretty_condition(row$Condition),
        format_deg_jp(row$MeanShiftDeg, 1, TRUE), format_deg_jp(row$MedianShiftDeg, 1, TRUE),
        format_ci_jp(row$CI_lo, row$CI_hi, 1, TRUE), format_p(row$SignTestP), format_p(row$PermutationP), row$N
      ))
    }
  }
  if (is.data.frame(within_fish$shift_pair_age_correlations) && nrow(within_fish$shift_pair_age_correlations) > 0) {
    lines <- c(lines, "AgeGroup-specific shift correlations", "AgeGroup | ConditionA vs ConditionB | Pearson r | permutation p | Spearman rho | permutation p | same-sign fraction | N")
    for (i in seq_len(nrow(within_fish$shift_pair_age_correlations))) {
      row <- within_fish$shift_pair_age_correlations[i, ]
      lines <- c(lines, sprintf("%s | %s vs %s | %s | %s | %s | %s | %s | %d",
        row$AgeGroup,
        pretty_condition(row$ConditionA), pretty_condition(row$ConditionB),
        format_num(row$PearsonR, 3), format_p(row$PearsonP),
        format_num(row$SpearmanR, 3), format_p(row$SpearmanP),
        paste0(format_num(row$SameSignFraction * 100, 0), "%"),
        row$N
      ))
    }
    lines <- c(lines, "If correlation signs differ by AgeGroup, developmental stage may be contributing to the heterogeneity in magnetic response patterns.")
  }
  if (is.list(shift_clusters) && is.data.frame(shift_clusters$stats) && nrow(shift_clusters$stats) > 0) {
    st <- shift_clusters$stats[1, ]
    lines <- c(lines, "ShiftDeg clustering", "Method | N | dimensions | between/total SS | permutation p | judgement")
    lines <- c(lines, sprintf("k-means (k=2) | %d | %d | %s | %s | %s",
      st$N, st$Dimensions, format_num(st$BetweenTotalSS, 3), format_p(st$PermutationP),
      ifelse(isTRUE(st$SuggestedSubpopulation), "possible subpopulation", "no clear evidence of subpopulation")
    ))
    if (is.data.frame(shift_clusters$summary) && nrow(shift_clusters$summary) > 0) {
      lines <- c(lines, "Cluster | role | N | mean distance")
      for (i in seq_len(nrow(shift_clusters$summary))) {
        row <- shift_clusters$summary[i, ]
        lines <- c(lines, sprintf("%s | %s | %d | %s", row$Cluster, row$ClusterRole, row$N, format_num(row$MeanDistanceFromZero, 2)))
      }
    }
    if (isTRUE(st$SuggestedSubpopulation)) {
      lines <- c(lines, "A two-group structure consistent with responder and non-responder patterns was exploratorily suggested in multivariate ShiftDeg space.")
    } else {
      lines <- c(lines, "There is currently no strong statistical evidence for a ShiftDeg subpopulation structure.")
    }
  }
  if (is.data.frame(window_sensitivity$pairwise) && nrow(window_sensitivity$pairwise) > 0) {
    lines <- c(lines, "", "6. Window sensitivity")
    labels <- unique(as.character(window_sensitivity$pairwise$WindowLabel))
    lines <- c(lines, paste(c("Comparison", labels), collapse = " | "))
    pair_keys <- unique(paste(window_sensitivity$pairwise$GroupA, window_sensitivity$pairwise$GroupB, sep = " | "))
    for (pk in pair_keys) {
      sub <- window_sensitivity$pairwise[paste(window_sensitivity$pairwise$GroupA, window_sensitivity$pairwise$GroupB, sep = " | ") == pk, , drop = FALSE]
      sub <- sub[match(labels, as.character(sub$WindowLabel)), , drop = FALSE]
      vals <- sapply(seq_len(nrow(sub)), function(i) format_deg_jp(sub$MeanDeltaDeg[i], 1, TRUE))
      lines <- c(lines, paste(c(sprintf("%s - %s", pretty_condition(sub$GroupA[1]), pretty_condition(sub$GroupB[1])), vals), collapse = " | "))
    }
  }
  lines <- c(lines, "", "Overall interpretation", "Within-fish consistency was the strongest signal, suggesting individual differences that may cancel at the population mean.")
  if (is.list(shift_clusters) && is.data.frame(shift_clusters$stats) && nrow(shift_clusters$stats) > 0 && isTRUE(shift_clusters$stats$SuggestedSubpopulation[1])) {
    lines <- c(lines, "ShiftDeg clustering suggested a candidate responder/non-responder split, consistent with fish-specific magnetic response patterns.")
  }
  lines
}

readme_lines_builder <- function(alignment) {
  c(
    "Outputs:",
    "1. circular_model_input.csv",
    "2. circular_mean_difference_ci.csv",
    "3. individual_alignment_rows.csv / individual_alignment_summary.csv",
    "4. individual_alignment_pairwise_ci.csv",
    "5. cohort_pairwise_effects.csv / cohort_direction_consistency.csv / cohort_activity_pairwise_tests.csv",
    "5b. age_group_pairwise_effects.csv / age_group_activity_pairwise_tests.csv / age_*_mixed_model_*.csv",
    "5c. order_activity_pairwise_tests.csv / order_activity_mixed_model_terms_summary.csv / order_activity_*_mixed_model_*.csv",
    "5d. subset_order_activity_pairwise_tests.csv / subset_order_activity_mixed_model_terms_summary.csv / subset_order_activity_model_meta.csv / subset_order_*_*.csv",
    "6. window_sensitivity_pairwise.csv / window_sensitivity_windows.csv",
    "7. quality_sensitivity_pairwise.csv",
    "8. within_fish_reference_correlations.csv / within_fish_shift_summary.csv / within_fish_shift_distribution.csv / within_fish_shift_cohort_summary.csv / within_fish_shift_pair_correlations.csv / within_fish_shift_pair_cohort_correlations.csv / within_fish_shift_age_summary.csv / within_fish_shift_pair_age_correlations.csv / within_fish_shift_cluster_*.csv / within_fish_shift_clusters.png",
    "9. cos_mixed_model_anova.csv / sin_mixed_model_anova.csv",
    "10. activity_pairwise_tests.csv / activity_*_anova.csv / activity_*_coefficients.csv",
    "11. bayesian_circular_summary.txt and coefficient csv files",
    "12. mixed_model_terms_summary.csv / activity_mixed_model_terms_summary.csv",
    "13. results_interpretation.txt",
    paste0("Alignment reference condition: ", alignment$reference_condition)
  )
}

if (nzchar(condition_arg)) keep_conditions <- tolower(trimws(strsplit(condition_arg, ",")[[1]])) else keep_conditions <- NULL
dat_window <- prepare_sample_data(dat_all, window_start, window_end, keep_conditions)
if (nrow(dat_window) == 0) stop("No sample rows available after filtering.")

activity_metrics <- choose_activity_metrics(dat_window)
agg <- finalize_agg(aggregate_fish_condition(dat_window, activity_metrics = activity_metrics))
write.csv(agg, file.path(out_dir, "circular_model_input.csv"), row.names = FALSE)

pairwise_ci <- pairwise_circular_ci(agg, n_boot = 2000)
pairwise_ci <- add_significance_flags(pairwise_ci, lo_col = "CI_lo", hi_col = "CI_hi", p_col = "")
write.csv(pairwise_ci, file.path(out_dir, "circular_mean_difference_ci.csv"), row.names = FALSE)

activity_pairwise <- pairwise_continuous_tests(agg, activity_metrics)
activity_pairwise <- add_significance_flags(activity_pairwise, lo_col = "CI_lo", hi_col = "CI_hi", p_col = "p_adj")
write.csv(activity_pairwise, file.path(out_dir, "activity_pairwise_tests.csv"), row.names = FALSE)

alignment <- individual_alignment_analysis(agg)
write.csv(alignment$aligned, file.path(out_dir, "individual_alignment_rows.csv"), row.names = FALSE)
write.csv(alignment$by_condition, file.path(out_dir, "individual_alignment_summary.csv"), row.names = FALSE)
write.csv(pairwise_circular_ci(alignment$aligned, angle_col = "AlignedHeadingDeg", group_col = "Condition", id_col = "Pair", n_boot = 2000), file.path(out_dir, "individual_alignment_pairwise_ci.csv"), row.names = FALSE)

cohort_pairwise <- cohort_pairwise_analysis(agg, n_boot = 1000)
cohort_pairwise <- add_significance_flags(cohort_pairwise, lo_col = "CI_lo", hi_col = "CI_hi", p_col = "")
write.csv(cohort_pairwise, file.path(out_dir, "cohort_pairwise_effects.csv"), row.names = FALSE)
cohort_consistency <- cohort_consistency_summary(cohort_pairwise, pairwise_ci)
write.csv(cohort_consistency, file.path(out_dir, "cohort_direction_consistency.csv"), row.names = FALSE)
cohort_activity_pairwise <- cohort_pairwise_continuous_analysis(agg, activity_metrics)
cohort_activity_pairwise <- add_significance_flags(cohort_activity_pairwise, lo_col = "CI_lo", hi_col = "CI_hi", p_col = "p_adj")
write.csv(cohort_activity_pairwise, file.path(out_dir, "cohort_activity_pairwise_tests.csv"), row.names = FALSE)
age_pairwise <- age_group_pairwise_analysis(agg, n_boot = 1000)
age_pairwise <- add_significance_flags(age_pairwise, lo_col = "CI_lo", hi_col = "CI_hi", p_col = "")
write.csv(age_pairwise, file.path(out_dir, "age_group_pairwise_effects.csv"), row.names = FALSE)
age_activity_pairwise <- age_group_pairwise_continuous_analysis(agg, activity_metrics)
age_activity_pairwise <- add_significance_flags(age_activity_pairwise, lo_col = "CI_lo", hi_col = "CI_hi", p_col = "p_adj")
write.csv(age_activity_pairwise, file.path(out_dir, "age_group_activity_pairwise_tests.csv"), row.names = FALSE)
order_activity_pairwise <- order_pairwise_continuous_analysis(agg, activity_metrics)
order_activity_pairwise <- add_significance_flags(order_activity_pairwise, lo_col = "CI_lo", hi_col = "CI_hi", p_col = "p_adj")
write.csv(order_activity_pairwise, file.path(out_dir, "order_activity_pairwise_tests.csv"), row.names = FALSE)
subset_order_activity <- subset_order_activity_analysis(agg, activity_metrics, families = c("Sea", "Intensity"))
if (is.data.frame(subset_order_activity$pairwise) && nrow(subset_order_activity$pairwise) > 0) {
  subset_order_activity$pairwise <- add_significance_flags(subset_order_activity$pairwise, lo_col = "CI_lo", hi_col = "CI_hi", p_col = "p_adj")
}
write.csv(subset_order_activity$pairwise, file.path(out_dir, "subset_order_activity_pairwise_tests.csv"), row.names = FALSE)
write.csv(subset_order_activity$model_terms, file.path(out_dir, "subset_order_activity_mixed_model_terms_summary.csv"), row.names = FALSE)
write.csv(subset_order_activity$model_meta, file.path(out_dir, "subset_order_activity_model_meta.csv"), row.names = FALSE)

window_sensitivity <- window_sensitivity_analysis(dat_all, window_start, window_end, keep_conditions)
write.csv(window_sensitivity$pairwise, file.path(out_dir, "window_sensitivity_pairwise.csv"), row.names = FALSE)
write.csv(window_sensitivity$windows, file.path(out_dir, "window_sensitivity_windows.csv"), row.names = FALSE)

quality_sensitivity <- quality_sensitivity_analysis(agg)
write.csv(quality_sensitivity, file.path(out_dir, "quality_sensitivity_pairwise.csv"), row.names = FALSE)

within_fish <- within_fish_consistency_analysis(agg, alignment)
write.csv(within_fish$reference_correlations, file.path(out_dir, "within_fish_reference_correlations.csv"), row.names = FALSE)
write.csv(within_fish$shift_summary, file.path(out_dir, "within_fish_shift_summary.csv"), row.names = FALSE)
write.csv(within_fish$shift_distribution, file.path(out_dir, "within_fish_shift_distribution.csv"), row.names = FALSE)
write.csv(within_fish$shift_cohort_summary, file.path(out_dir, "within_fish_shift_cohort_summary.csv"), row.names = FALSE)
write.csv(within_fish$shift_age_summary, file.path(out_dir, "within_fish_shift_age_summary.csv"), row.names = FALSE)
write.csv(within_fish$shift_pair_correlations, file.path(out_dir, "within_fish_shift_pair_correlations.csv"), row.names = FALSE)
write.csv(within_fish$shift_pair_cohort_correlations, file.path(out_dir, "within_fish_shift_pair_cohort_correlations.csv"), row.names = FALSE)
write.csv(within_fish$shift_pair_age_correlations, file.path(out_dir, "within_fish_shift_pair_age_correlations.csv"), row.names = FALSE)
shift_cluster_data <- prepare_shift_cluster_data(within_fish)
shift_clusters <- evaluate_shift_clusters(shift_cluster_data)
if (is.data.frame(shift_clusters$assignments) && nrow(shift_clusters$assignments) > 0) write.csv(shift_clusters$assignments, file.path(out_dir, "within_fish_shift_cluster_assignments.csv"), row.names = FALSE)
if (is.data.frame(shift_clusters$summary) && nrow(shift_clusters$summary) > 0) write.csv(shift_clusters$summary, file.path(out_dir, "within_fish_shift_cluster_summary.csv"), row.names = FALSE)
if (is.data.frame(shift_clusters$stats) && nrow(shift_clusters$stats) > 0) write.csv(shift_clusters$stats, file.path(out_dir, "within_fish_shift_cluster_stats.csv"), row.names = FALSE)
if (is.data.frame(shift_clusters$centers) && nrow(shift_clusters$centers) > 0) write.csv(shift_clusters$centers, file.path(out_dir, "within_fish_shift_cluster_centers.csv"), row.names = FALSE)
plot_shift_clusters(shift_clusters$assignments, shift_clusters$centers, out_dir)

cos_result <- fit_lmer_side(agg, "cosMag")
sin_result <- fit_lmer_side(agg, "sinMag")
write_model_result(cos_result, "cos_mixed_model", out_dir)
write_model_result(sin_result, "sin_mixed_model", out_dir)
mixed_terms <- rbind(extract_mixed_model_terms(cos_result, "cosMag"), extract_mixed_model_terms(sin_result, "sinMag"))
mixed_terms <- add_condition_p_adjustment(mixed_terms)
write.csv(mixed_terms, file.path(out_dir, "mixed_model_terms_summary.csv"), row.names = FALSE)

activity_results <- lapply(activity_metrics, function(metric) fit_lmer_side(agg, metric))
names(activity_results) <- activity_metrics
for (metric in activity_metrics) {
  prefix <- paste0("activity_", metric)
  write_model_result(activity_results[[metric]], prefix, out_dir)
}
activity_mixed_terms <- do.call(rbind, lapply(activity_metrics, function(metric) extract_mixed_model_terms(activity_results[[metric]], metric)))
activity_mixed_terms <- add_condition_p_adjustment(activity_mixed_terms)
write.csv(activity_mixed_terms, file.path(out_dir, "activity_mixed_model_terms_summary.csv"), row.names = FALSE)

age_cos_result <- fit_age_lmer_side(agg, "cosMag")
age_sin_result <- fit_age_lmer_side(agg, "sinMag")
write_model_result(age_cos_result, "age_cos_mixed_model", out_dir)
write_model_result(age_sin_result, "age_sin_mixed_model", out_dir)
age_activity_results <- lapply(activity_metrics, function(metric) fit_age_lmer_side(agg, metric))
names(age_activity_results) <- activity_metrics
for (metric in activity_metrics) {
  prefix <- paste0("age_activity_", metric)
  write_model_result(age_activity_results[[metric]], prefix, out_dir)
}
age_mixed_terms_parts <- c(
  list(extract_age_terms(age_cos_result, "cosMag"), extract_age_terms(age_sin_result, "sinMag")),
  lapply(activity_metrics, function(metric) extract_age_terms(age_activity_results[[metric]], metric))
)
age_mixed_terms_parts <- age_mixed_terms_parts[vapply(age_mixed_terms_parts, is.data.frame, logical(1))]
age_mixed_terms <- if (length(age_mixed_terms_parts) > 0) do.call(rbind, age_mixed_terms_parts) else data.frame()
if (is.data.frame(age_mixed_terms) && nrow(age_mixed_terms) > 0) {
  write.csv(age_mixed_terms, file.path(out_dir, "age_mixed_model_terms_summary.csv"), row.names = FALSE)
}
order_activity_results <- lapply(activity_metrics, function(metric) fit_order_lmer_side(agg, metric))
names(order_activity_results) <- activity_metrics
for (metric in activity_metrics) {
  prefix <- paste0("order_activity_", metric)
  write_model_result(order_activity_results[[metric]], prefix, out_dir)
}
order_mixed_terms_parts <- lapply(activity_metrics, function(metric) extract_order_terms(order_activity_results[[metric]], metric))
order_mixed_terms_parts <- order_mixed_terms_parts[vapply(order_mixed_terms_parts, is.data.frame, logical(1))]
order_mixed_terms <- if (length(order_mixed_terms_parts) > 0) do.call(rbind, order_mixed_terms_parts) else data.frame()
if (is.data.frame(order_mixed_terms) && nrow(order_mixed_terms) > 0) {
  write.csv(order_mixed_terms, file.path(out_dir, "order_activity_mixed_model_terms_summary.csv"), row.names = FALSE)
}
if (length(subset_order_activity$model_results) > 0) {
  for (nm in names(subset_order_activity$model_results)) {
    parts <- strsplit(nm, "__", fixed = TRUE)[[1]]
    if (length(parts) != 3) next
    fam <- tolower(parts[1])
    age <- tolower(parts[2])
    metric <- parts[3]
    prefix <- paste0("subset_order_", fam, "_", age, "_", metric)
    write_model_result(subset_order_activity$model_results[[nm]], prefix, out_dir)
  }
}

bayes_result <- fit_bayesian_circular(agg)
bayes_out <- write_bayesian_outputs(bayes_result, out_dir)

writeLines(build_results_interpretation_v2(agg, pairwise_ci, alignment, mixed_terms, bayes_out, cohort_pairwise, cohort_consistency, cohort_activity_pairwise, age_pairwise, age_activity_pairwise, age_mixed_terms, order_activity_pairwise, order_mixed_terms, subset_order_activity, window_sensitivity, quality_sensitivity, within_fish, shift_clusters, activity_pairwise, activity_mixed_terms, window_start, window_end), file.path(out_dir, "results_interpretation.txt"))
writeLines(readme_lines_builder(alignment), file.path(out_dir, "README_circular_models.txt"))

message("Wrote circular mixed-model outputs to: ", out_dir)
