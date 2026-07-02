suppressPackageStartupMessages({
  library(survival)
  library(glmnet)
  library(dplyr)
  library(purrr)
  library(tibble)
})

z_scale <- function(x) {
  as.numeric(scale(x))
}

complete_terms <- function(data, cols) {
  cols <- intersect(cols, names(data))
  stats::complete.cases(data[, cols, drop = FALSE])
}

bh_q <- function(p) {
  stats::p.adjust(p, method = "BH")
}

safe_ci <- function(beta, se) {
  c(exp(beta), exp(beta - 1.96 * se), exp(beta + 1.96 * se))
}

make_survival_formula <- function(time_col, event_col, exposure, covariates) {
  stats::as.formula(
    paste0("survival::Surv(", time_col, ", ", event_col, ") ~ ", paste(c(exposure, covariates), collapse = " + "))
  )
}

fit_cox_model <- function(data, time_col, event_col, exposure, covariates) {
  needed <- c(time_col, event_col, exposure, covariates)
  analytic <- data[complete_terms(data, needed), needed, drop = FALSE]
  fit <- survival::coxph(make_survival_formula(time_col, event_col, exposure, covariates), data = analytic, x = TRUE)
  s <- summary(fit)
  idx <- match(exposure, rownames(s$coefficients))
  beta <- s$coefficients[idx, "coef"]
  se <- s$coefficients[idx, "se(coef)"]
  ci <- safe_ci(beta, se)
  tibble(
    exposure = exposure,
    n = fit$n,
    events = fit$nevent,
    beta = beta,
    se = se,
    hr = ci[1],
    lo = ci[2],
    hi = ci[3],
    p = s$coefficients[idx, "Pr(>|z|)"]
  )
}

fit_outcome_exposure_grid <- function(data, outcome_specs, exposures, covariates) {
  purrr::map_dfr(names(outcome_specs), function(outcome_name) {
    purrr::map_dfr(exposures, function(exposure) {
      res <- fit_cox_model(
        data = data,
        time_col = outcome_specs[[outcome_name]]$time,
        event_col = outcome_specs[[outcome_name]]$event,
        exposure = exposure,
        covariates = covariates
      )
      res$outcome <- outcome_name
      res
    })
  }) |>
    group_by(outcome) |>
    mutate(q = bh_q(p)) |>
    ungroup()
}

sequential_covariate_sets <- function(base_covariates, nutrient_covariates = character(), adiposity_covariates = character(), biomarker_covariates = character()) {
  list(
    base = base_covariates,
    nutrient_adjusted = c(base_covariates, nutrient_covariates),
    adiposity_adjusted = c(base_covariates, nutrient_covariates, adiposity_covariates),
    biomarker_adjusted = c(base_covariates, nutrient_covariates, adiposity_covariates, biomarker_covariates)
  )
}

survival_at <- function(fit, data, time_point) {
  sf <- survival::survfit(fit, newdata = data)
  ss <- summary(sf, times = time_point, extend = TRUE)
  as.numeric(ss$surv)
}

standardized_risk_difference <- function(fit, data, exposure, high_value, low_value, time_point) {
  high_data <- data
  low_data <- data
  high_data[[exposure]] <- high_value
  low_data[[exposure]] <- low_value
  risk_high <- 1 - mean(survival_at(fit, high_data, time_point), na.rm = TRUE)
  risk_low <- 1 - mean(survival_at(fit, low_data, time_point), na.rm = TRUE)
  1000 * (risk_high - risk_low)
}

simulate_risk_difference <- function(fit, data, exposure, high_value, low_value, time_point, draws = 500, seed = 1) {
  set.seed(seed)
  beta_hat <- stats::coef(fit)
  vc <- stats::vcov(fit)
  coef_draws <- MASS::mvrnorm(draws, mu = beta_hat, Sigma = vc)
  original_coef <- beta_hat
  rd <- numeric(draws)
  for (i in seq_len(draws)) {
    fit$coefficients <- coef_draws[i, ]
    rd[i] <- standardized_risk_difference(fit, data, exposure, high_value, low_value, time_point)
  }
  fit$coefficients <- original_coef
  tibble(
    rd_per_1000 = stats::median(rd, na.rm = TRUE),
    rd_lo_per_1000 = stats::quantile(rd, 0.025, na.rm = TRUE),
    rd_hi_per_1000 = stats::quantile(rd, 0.975, na.rm = TRUE)
  )
}

exposure_contrast_values <- function(data, exposure) {
  qs <- stats::quantile(data[[exposure]], probs = c(0.2, 0.8), na.rm = TRUE)
  low <- mean(data[[exposure]][data[[exposure]] <= qs[[1]]], na.rm = TRUE)
  high <- mean(data[[exposure]][data[[exposure]] >= qs[[2]]], na.rm = TRUE)
  c(low = low, high = high)
}

scan_endpoints <- function(data, endpoint_table, exposure, covariates, minimum_events = 100) {
  rows <- purrr::pmap_dfr(endpoint_table, function(endpoint, time, event, ...) {
    events <- sum(data[[event]] == 1, na.rm = TRUE)
    if (events < minimum_events) {
      return(tibble(endpoint = endpoint, n = NA_integer_, events = events, hr = NA_real_, lo = NA_real_, hi = NA_real_, p = NA_real_))
    }
    res <- fit_cox_model(data, time, event, exposure, covariates)
    res$endpoint <- endpoint
    res
  })
  rows |>
    mutate(q = bh_q(p)) |>
    arrange(q, p)
}

count_positive_endpoints <- function(scan_results, alpha = 0.05) {
  scan_results |>
    filter(!is.na(q), q < alpha, hr > 1) |>
    summarise(n_positive = dplyr::n(), .groups = "drop")
}

make_folds <- function(n, k = 5, seed = 1) {
  set.seed(seed)
  sample(rep(seq_len(k), length.out = n))
}

residualize_exposure <- function(data, exposure, covariates) {
  form <- stats::as.formula(paste(exposure, "~", paste(covariates, collapse = " + ")))
  stats::resid(stats::lm(form, data = data))
}

fit_crossfit_biomarker_profile <- function(data, exposure, marker_cols, covariates, k = 5, alpha = 0.5, seed = 1) {
  marker_cols <- intersect(marker_cols, names(data))
  needed <- c(exposure, marker_cols, covariates)
  analytic <- data[complete_terms(data, needed), needed, drop = FALSE]
  y <- residualize_exposure(analytic, exposure, covariates)
  x <- as.matrix(analytic[, marker_cols, drop = FALSE])
  for (j in seq_len(ncol(x))) {
    x[is.na(x[, j]), j] <- stats::median(x[, j], na.rm = TRUE)
  }
  folds <- make_folds(nrow(analytic), k, seed)
  pred <- rep(NA_real_, nrow(analytic))
  for (fold in seq_len(k)) {
    train <- folds != fold
    test <- folds == fold
    x_train <- scale(x[train, , drop = FALSE])
    center <- attr(x_train, "scaled:center")
    scalev <- attr(x_train, "scaled:scale")
    x_test <- sweep(sweep(x[test, , drop = FALSE], 2, center, "-"), 2, scalev, "/")
    fit <- glmnet::cv.glmnet(x_train, y[train], alpha = alpha, family = "gaussian")
    pred[test] <- as.numeric(stats::predict(fit, newx = x_test, s = "lambda.min"))
  }
  out <- analytic
  out$biomarker_profile <- z_scale(pred)
  out
}

profile_performance <- function(profile_data, exposure, outcome_specs) {
  exposure_perf <- stats::cor(profile_data$biomarker_profile, profile_data[[exposure]], use = "complete.obs")
  tibble(metric = "exposure_correlation", value = exposure_perf) |>
    bind_rows(purrr::map_dfr(names(outcome_specs), function(outcome_name) {
      fit <- fit_cox_model(profile_data, outcome_specs[[outcome_name]]$time, outcome_specs[[outcome_name]]$event, "biomarker_profile", character())
      fit$metric <- paste0("outcome_", outcome_name)
      fit |> select(metric, hr, lo, hi, p)
    }))
}

fit_interaction <- function(data, time_col, event_col, exposure, modifier, covariates) {
  data[[exposure]] <- z_scale(data[[exposure]])
  data[[modifier]] <- z_scale(data[[modifier]])
  needed <- c(time_col, event_col, exposure, modifier, covariates)
  analytic <- data[complete_terms(data, needed), needed, drop = FALSE]
  form <- stats::as.formula(
    paste0("survival::Surv(", time_col, ", ", event_col, ") ~ ", exposure, " * ", modifier, " + ", paste(covariates, collapse = " + "))
  )
  fit <- survival::coxph(form, data = analytic)
  s <- summary(fit)
  term <- paste0(exposure, ":", modifier)
  idx <- match(term, rownames(s$coefficients))
  beta <- s$coefficients[idx, "coef"]
  se <- s$coefficients[idx, "se(coef)"]
  ci <- safe_ci(beta, se)
  tibble(exposure = exposure, modifier = modifier, n = fit$n, events = fit$nevent, hr_interaction = ci[1], lo = ci[2], hi = ci[3], p = s$coefficients[idx, "Pr(>|z|)"])
}

fit_interaction_grid <- function(data, outcome_specs, exposures, modifiers, covariates) {
  purrr::map_dfr(names(outcome_specs), function(outcome_name) {
    out <- outcome_specs[[outcome_name]]
    purrr::map_dfr(exposures, function(exposure) {
      purrr::map_dfr(modifiers, function(modifier) {
        res <- fit_interaction(data, out$time, out$event, exposure, modifier, covariates)
        res$outcome <- outcome_name
        res
      })
    })
  }) |>
    group_by(outcome) |>
    mutate(q = bh_q(p)) |>
    ungroup()
}

make_joint_strata <- function(data, exposure, modifier) {
  e <- cut(data[[exposure]], breaks = stats::quantile(data[[exposure]], c(0, 1 / 3, 2 / 3, 1), na.rm = TRUE), include.lowest = TRUE, labels = c("low", "middle", "high"))
  g <- cut(data[[modifier]], breaks = stats::quantile(data[[modifier]], c(0, 1 / 3, 2 / 3, 1), na.rm = TRUE), include.lowest = TRUE, labels = c("low", "middle", "high"))
  interaction(e, g, sep = "_")
}

check_proportional_hazards <- function(data, time_col, event_col, exposure, covariates, max_n = 80000, seed = 1) {
  needed <- c(time_col, event_col, exposure, covariates)
  analytic <- data[complete_terms(data, needed), needed, drop = FALSE]
  if (nrow(analytic) > max_n) {
    set.seed(seed)
    analytic <- analytic[sample(seq_len(nrow(analytic)), max_n), , drop = FALSE]
  }
  fit <- survival::coxph(make_survival_formula(time_col, event_col, exposure, covariates), data = analytic)
  z <- survival::cox.zph(fit)
  p <- z$table[exposure, "p"]
  tibble(n = fit$n, events = fit$nevent, exposure = exposure, p = p)
}

early_event_exclusion <- function(data, time_col, event_col, lag_time) {
  keep <- is.na(data[[time_col]]) | data[[event_col]] == 0 | data[[time_col]] > lag_time
  data[keep, , drop = FALSE]
}

fit_lagged_sensitivity <- function(data, outcome_specs, exposure, covariates, lag_time) {
  purrr::map_dfr(names(outcome_specs), function(outcome_name) {
    out <- outcome_specs[[outcome_name]]
    lagged <- early_event_exclusion(data, out$time, out$event, lag_time)
    res <- fit_cox_model(lagged, out$time, out$event, exposure, covariates)
    res$outcome <- outcome_name
    res$lag_time <- lag_time
    res
  })
}

fit_stratified_consistency <- function(data, outcome_specs, exposure, covariates, strata_col, minimum_events = 100) {
  purrr::map_dfr(sort(unique(data[[strata_col]])), function(level) {
    subset_data <- data[data[[strata_col]] == level, , drop = FALSE]
    purrr::map_dfr(names(outcome_specs), function(outcome_name) {
      out <- outcome_specs[[outcome_name]]
      if (sum(subset_data[[out$event]] == 1, na.rm = TRUE) < minimum_events) {
        return(tibble())
      }
      res <- fit_cox_model(subset_data, out$time, out$event, exposure, covariates)
      res$outcome <- outcome_name
      res$stratum <- level
      res$stratum_variable <- strata_col
      res
    })
  })
}
