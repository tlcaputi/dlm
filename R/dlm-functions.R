#' Reverse Cumulative Sum
#'
#' Gives the reverse cumulative sum of a vector.
#' @param x A vector
#' @return A vector 
#' @export
#' 
revcumsum <- function(x){
  x <- rev(cumsum(rev(x)))
}

#' Standard Error Cumulative Sum
#'
#' Gives the cumulative sum of a variance-covariance matrix.
#' @param cov A variance-covariance matrix
#' @return A vector of SEs
#' @export
#' 
secumsum <- function(cov) {
    cov <- as.matrix(cov)
    L <- dim(cov)[1]
    if (L == 1) {
        return(c(sqrt(cov)))
    }
    se <- numeric(L)
    for (i in seq_len(L)) {
        a <- matrix(rep(1, i), nrow = 1)
        se[i] <- sqrt(a %*% cov[1:i, 1:i, drop = FALSE] %*% t(a))
    }
    se
}




#' Reverse Standard Error Cumulative Sum
#'
#' Gives the cumulative sum of a variance-covariance matrix.
#' @param cov A variance-covariance matrix
#' @return A vector of SEs
#' @export
#' 
serevcumsum <- function(cov) {
    cov <- as.matrix(cov)
    L <- dim(cov)[1]
    if (L == 1) {
        return(c(sqrt(cov)))
    }
    se <- numeric(L)
    for (i in seq.int(L, 1)) {
        a <- matrix(rep(1, L - i + 1), nrow = 1)
        se[i] <- sqrt(a %*% cov[i:L, i:L, drop = FALSE] %*% t(a))
    }
    se
}

# ---- Internal helpers (not exported) ----

# Prepare exposure data: filter NAs, deduplicate, check/fix balance
.prepare_exposure <- function(data, exposure_data, exposure, unit, time) {
  if (exposure %in% names(data)) {
    data <- data %>% select(-!!sym(exposure))
  }
  data <- data %>% filter(!is.na(!!sym(unit)), !is.na(!!sym(time)))
  exposure_data <- exposure_data %>%
    filter(!is.na(!!sym(unit)), !is.na(!!sym(time)), !is.na(!!sym(exposure))) %>%
    select(!!sym(unit), !!sym(time), !!sym(exposure)) %>%
    unique()

  if (!plm::is.pbalanced(exposure_data, index = c(unit, time))) {
    warning("The exposure data is not balanced. Filling missing values with NAs.")
    min_t <- min(exposure_data[[time]], na.rm = TRUE)
    max_t <- max(exposure_data[[time]], na.rm = TRUE)
    exposure_data <- tidyr::complete(
      exposure_data,
      !!rlang::sym(unit),
      !!rlang::sym(time) := min_t:max_t
    )
    if (!plm::is.pbalanced(exposure_data, index = c(unit, time))) {
      stop("The exposure data is still not balanced after filling missing values.")
    }
  }
  list(data = data, exposure_data = exposure_data)
}

# Create all leads and lags of exposure using data.table::shift (single grouped op)
.create_leads_lags <- function(exposure_data, exposure, unit, time, from_rt, to_rt) {
  dt <- data.table::as.data.table(exposure_data)
  data.table::setorderv(dt, c(unit, time))

  leads <- character(0)
  lags <- character(0)

  # Leads: from_rt < 0 means abs(from_rt)-1 leads
  n_leads <- abs(from_rt) - 1
  if (n_leads > 0) {
    lead_ns <- n_leads:1
    leads <- paste0(exposure, "_lead", lead_ns)
    dt[, (leads) := data.table::shift(get(exposure), n = lead_ns, type = "lead"), by = c(unit)]
  }

  # Lags: 0 through to_rt
  lag_ns <- 0:to_rt
  lags <- paste0(exposure, "_lag", lag_ns)
  dt[, (lags) := data.table::shift(get(exposure), n = lag_ns, type = "lag"), by = c(unit)]

  list(exposure_data = as.data.frame(dt), leads = leads, lags = lags)
}

# Convert gamma coefficients to beta (event-study) coefficients
.gamma_to_beta <- function(gamma, vcov, from_rt, to_rt, ref_period) {
  num_vars <- abs(from_rt) + to_rt
  time_to_event <- sort(setdiff(from_rt:to_rt, ref_period))

  if (from_rt == ref_period) {
    after_periods <- 1:num_vars
    coefs <- cumsum(gamma[after_periods])
    ses <- secumsum(as.matrix(vcov)[after_periods, after_periods, drop = FALSE])
  } else {
    num_before <- length(time_to_event[time_to_event < ref_period])
    before_periods <- 1:num_before
    after_periods <- (num_before + 1):num_vars
    coefs <- c(
      -revcumsum(gamma[before_periods]),
      cumsum(gamma[after_periods])
    )
    ses <- c(
      serevcumsum(vcov[before_periods, before_periods, drop = FALSE]),
      secumsum(vcov[after_periods, after_periods, drop = FALSE])
    )
  }

  data.frame(time_to_event = time_to_event, coef = coefs, se = ses)
}


#' Generate Data
#'
#' Produces a dataset perfect for a standard TWFE model
#' @param seed Seed for the random data
#' @param n_groups Number of groups
#' @param n_times Number of time periods
#' @param treat_prob Probability of treatment
#' @return A dataset 
#' @export
#' 
generate_data = function(seed=1234, n_groups = 26^2, n_times = 20, treat_prob = 0.4){

    ## This generates a very simple dataset for a binary, standard TWFE model.
    ## There are 26^2 groups and 20 time periods. 40% of units are treated at a random time between 7:9.
    ## 60% are never treated. 0% are always treated.

    old_seed <- globalenv()$.Random.seed
    on.exit({
        if (is.null(old_seed)) {
            rm(".Random.seed", envir = globalenv())
        } else {
            assign(".Random.seed", old_seed, envir = globalenv())
        }
    })
    set.seed(seed)

    # Create a panel of n_times time periods and n_groups groups
    times = 1:n_times
    groups = glue("group{1:n_groups}")
    panel = expand.grid(group = groups, time = times)
    # Give treat_prob of units a random treatment time between 7:9
    treatments = data.frame(group = groups) %>%
        dplyr::sample_frac(treat_prob) %>%
        dplyr::mutate(
            treatment_time = sample(7:9, size = dplyr::n(), replace = TRUE)
        )
    df = merge(panel, treatments, by=c("group"), all.x = T)

    # Define treatment, years_to_treat, post, and outcome
    df = df %>% mutate(
        # whether the unit was ever treated (at any point)
        treat = case_when(
            is.na(treatment_time) ~ 0,
            !is.na(treatment_time) ~ 1
        ),
        # how many years before/after treatment (-1000 for untreated units)
        years_to_treatment = case_when(
            treat == 0 ~ -1000,
            treat == 1 ~ time - treatment_time
        ),
        # before or after treatment
        post = case_when(
            treat == 0 ~ 0,
            treat == 1 ~ as.numeric(time >= treatment_time)
        ),
        # random outcome with effect of -3 after treatment
        outcome = stats::rnorm(dplyr::n(), sd = 5) + (I(years_to_treatment >= 0) * -3)
    )

    # Make it pretty
    df = df %>% arrange(group, time)

    # The result will have columns: group, time, treat, years_to_treatment, post, outcome
    # The effect of the event should be -3.
    return(df)

}



#' Distributed Lags Model
#'
#' This is the distributed lags model / continuous event study. 
#' @param data A data frame containing the unit, time, outcome, covariates, and any additional fixed effects
#' @param exposure_data A data frame containing the unit, time, and exposure variables
#' @param from_rt The starting lag period.
#' @param to_rt The ending lag period.
#' @param outcome The outcome variable.
#' @param exposure The exposure variable.
#' @param unit The unit identifier.
#' @param time The time variable.
#' @param covariates Vector of covariates for the model.
#' @param addl_fes Vector of additional fixed effects for the model.
#' @param ref_period Reference period (default -1)
#' @param weights Weights to be included in the regression
#' @return A list containing model results, coefficients, and plots.
#' @export
#' 
distributed_lags_model = function(data, exposure_data, from_rt, to_rt, outcome, exposure, unit, time, covariates = NULL, addl_fes = NULL, ref_period = -1, weights = NULL, dd=F, n=2, dict = NULL){

  for(v in c(unit, time, outcome, covariates, addl_fes)){
    if(!v %in% names(data)){
      warning(glue("Variable {v} not found in outcome data"))
    }
  }
  for(v in c(unit, time, exposure)){
    if(!v %in% names(exposure_data)){
      warning(glue("Variable {v} not found in exposure data"))
    }
  }
  if(!(ref_period %in% (from_rt:to_rt))){
    stop("ref_period must be in from_rt:to_rt")
  }

  # Prepare data
  prep <- .prepare_exposure(data, exposure_data, exposure, unit, time)
  data <- prep$data
  exposure_data <- prep$exposure_data

  MINTIME <- min(exposure_data[[time]], na.rm = TRUE)
  MAXTIME <- max(exposure_data[[time]], na.rm = TRUE)
  logger::log_info("MINTIME: {MINTIME}, MAXTIME: {MAXTIME}")

  data_years_included <- (MINTIME + abs(to_rt)):(MAXTIME - abs(from_rt) + 1)
  data_years_included <- intersect(data_years_included, unique(data[[time]]))

  # Create leads/lags efficiently
  logger::log_info("Creating leads and lags")
  ll <- .create_leads_lags(exposure_data, exposure, unit, time, from_rt, to_rt)
  exposure_data <- ll$exposure_data
  leads <- ll$leads
  lags <- ll$lags
  logger::log_info("Done creating leads and lags")

  # Merge
  tmp <- merge(data, exposure_data, by = c(unit, time), all.x = TRUE)
  tmp <- tmp %>% mutate(unit := !!sym(unit), time := !!sym(time))

  # Build formula
  leads_lags_str <- paste0(c(leads, lags), collapse = " + ")
  fes <- c("unit", "time", addl_fes)
  fe_str <- paste0(fes, collapse = " + ")
  if (!is.null(covariates)) {
    fmla_str <- glue("{outcome} ~ {leads_lags_str} + {paste0(covariates, collapse = ' + ')} | {fe_str}")
  } else {
    fmla_str <- glue("{outcome} ~ {leads_lags_str} | {fe_str}")
  }
  logger::log_info("Formula: {fmla_str}")

  # Estimate
  fmla <- stats::as.formula(fmla_str)
  if (!is.null(weights)) {
    wt_fmla <- stats::as.formula(paste0("~", weights))
    model <- fixest::feols(fmla, data = tmp, cluster = ~unit, fixef.rm = 'none', weights = wt_fmla)
  } else {
    model <- fixest::feols(fmla, data = tmp, cluster = ~unit, fixef.rm = 'none')
  }

  num_vars <- abs(from_rt) + to_rt
  coefficients <- model$coefficients
  if (length(coefficients) != (length(lags) + length(leads) + length(covariates))) {
    logger::log_info("Not all coefficients were estimated")
    print(coefficients)
    stop("Not all coefficients were estimated")
  }

  gamma <- model$coefficients[1:num_vars]
  vcov <- stats::vcov(model, cluster = ~unit)[1:num_vars, 1:num_vars]

  # Transform gamma -> beta
  betas <- .gamma_to_beta(gamma, vcov, from_rt, to_rt, ref_period)
  

  outcome_name = "Coefficient"
  exposure_name = "Time to Treatment"
  if(!is.null(dict)){
    if(outcome %in% names(dict)){
      outcome_name = dict[outcome]
    } 
    if(exposure %in% names(dict)){
      exposure_name = dict[exposure]
    } 
  }
  
  # This just creates a plot (adding in reference period)
  # and is not really necessary for the results
  plotdf = rbind.data.frame(betas, data.frame(time_to_event = ref_period, coef = 0, se = 0))
  plotdf = plotdf %>% mutate(
    time_to_event_str = case_when(
      time_to_event == from_rt ~ glue("{from_rt}+"),
      time_to_event == to_rt ~ glue("{to_rt}+"),
      T ~ as.character(time_to_event)
    )
  )
  plotdf = plotdf %>% arrange(time_to_event)
  plotdf %>% as.data.frame() %>% print()
  p = ggplot(plotdf, aes(x = time_to_event, y = coef))
  p = p + geom_line(color = "darkblue")
  p = p + geom_point(color = "darkblue")
  p = p + geom_errorbar(aes(ymin = coef - 1.96*se, ymax = coef + 1.96*se), width = 0.2, color = "darkblue")
  p = p + geom_hline(yintercept = 0, linetype = "dashed")
  p = p + geom_vline(xintercept = ref_period+0.5, linetype = "dashed")
  min_included_year = min(data_years_included, na.rm = T)
  max_included_year = max(data_years_included, na.rm = T)
  p = p + labs(x = exposure_name, y = outcome_name, caption = glue("N={comma(nobs(model))} | From {min_included_year} To {max_included_year} | {Sys.time()}"))
  p = p + theme_bw()

  if(dd){
    out = twfe_companion(
      data = data, 
      exposure_data = exposure_data, 
      from_rt = from_rt, 
      to_rt = to_rt, 
      outcome = outcome, 
      exposure = exposure, 
      unit = unit, 
      time = time, 
      covariates = covariates, 
      addl_fes = addl_fes, 
      ref_period = ref_period, 
      weights = weights, 
      dd = dd, 
      n = n,
      remove_unit_FE = remove_unit_FE
    )
    # out = do.call(twfe_companion, list(...))
    p = add_caption_to_plot(p, out)
    # p = p + labs(caption = out)
  }

  return(list(
    "betas" = betas,
    "plot" = p,
    "model" = model,
    "vcov" = vcov,
    "data_years_included" = data_years_included,
    "fmla_str" = fmla_str,
    "from_rt" = from_rt,
    "to_rt" = to_rt,
    "exposure" = exposure,
    "outcome" = outcome
  ))

}



#' Standard Two-Way Fixed Effects Model for Comparison
#'
#' This produces the standard two-way fixed effects model for comparison to the distributed lags model
#' @param data A data frame containing the main dataset.
#' @param from_rt The starting lag period.
#' @param to_rt The ending lag period.
#' @param outcome The outcome variable.
#' @param time The time variable.
#' @param unit The unit identifier.
#' @param time_to_treatment The variable representing time to treatment.
#' @param treat The treatment indicator variable.
#' @param covariates Optional covariates for the model.
#' @param ref_period The reference period for summing coefficients (default is -1).
#' @return A list containing model results, coefficients, and plots.
#' @export
#' 
standard_twfe_for_comparison = function(data, from_rt, to_rt, outcome, time, unit, time_to_treatment, treat, covariates = NULL, ref_period = -1, weights = NULL){

    # Capture the minimum and maximum time
    MINTIME = min(data[[time]], na.rm = T)
    MAXTIME = max(data[[time]], na.rm = T)

    # It's easier to just assign names to the variables we want to use.
    # Although this means that the original data can't use these variable names...
    tmp = data %>% mutate(
        time_to_treatment := !!sym(time_to_treatment),
        treat := !!sym(treat),
        unit := !!sym(unit),
        time := !!sym(time),
        outcome := !!sym(outcome)
    ) %>% filter(time >= MINTIME + abs(to_rt), time <= MAXTIME - abs(from_rt) + 1)

    # This is the standard TWFE model with binning at from_rt and to_rt
    if(!is.null(covariates)){
        fmla = glue(
            "
               {outcome}
                ~
                i(time_to_treatment, treat, ref=c({ref_period}, -1000), bin=.('{from_rt}+' = ~x<={from_rt}, '{to_rt}+' = ~x>={to_rt}))
                + .[covariates]
                |
                unit
                + time
            "
        )
    } else {
        fmla = glue(
            "
                {outcome}
                ~
                i(time_to_treatment, treat, ref=c({ref_period}, -1000), bin=.('{from_rt}+' = ~x<={from_rt}, '{to_rt}+' = ~x>={to_rt}))
                |
                unit
                + time
            "
        )
    }

    arguments = c(
        "data = tmp",
        "cluster = ~unit"
    )
    if(!is.null(weights)){
        arguments = c(arguments, glue("weights = ~{weights}"))
    }
    cmd = glue("fixest::feols({fmla}, {paste0(arguments, collapse = ', ')})")

    # Actually estimates the model
    fixest_model = eval(parse(text=cmd))

    # We plot it and export the results
    fixest_plot = ggiplot(fixest_model)
    fixest_betas = as.data.frame(fixest_model$coeftable)
    names(fixest_betas) = c("coef", "se", "tval", "pval")

    return(list("betas" = fixest_betas, "plot" = fixest_plot, "model" = fixest_model, "cmd" = cmd))
}


#' TWFE Companion
#'
#' This is the companion TWFE, which adds the DD estimate to the bottom of the figure
#' @param data A data frame containing the unit, time, outcome, covariates, and any additional fixed effects
#' @param exposure_data A data frame containing the unit, time, and exposure variables
#' @param from_rt The starting lag period.
#' @param to_rt The ending lag period.
#' @param outcome The outcome variable.
#' @param exposure The exposure variable.
#' @param unit The unit identifier.
#' @param time The time variable.
#' @param covariates Vector of covariates for the model.
#' @param addl_fes Vector of additional fixed effects for the model.
#' @param ref_period Reference period (default -1)
#' @param weights Weights to be included in the regression
#' @return A list containing model results, coefficients, and plots.
#' @export
#' 
twfe_companion = function(data, exposure_data, from_rt, to_rt, outcome, exposure, unit, time, covariates = NULL, addl_fes = NULL, ref_period = -1, weights = NULL, dd = F, n = 2, remove_unit_FE = FALSE){
  
  # Capture the minimum and maximum time
  MINTIME = min(exposure_data[[time]], na.rm = T)
  MAXTIME = max(exposure_data[[time]], na.rm = T)

  df = merge(data, exposure_data, by=c(unit, time), all.x = T)
  logger::log_info("DATA N: {nrow(data)}")
  logger::log_info("DF N: {nrow(df)}")
  logger::log_info("EXPOSURE N: {nrow(exposure_data)}")

  # It's easier to just assign names to the variables we want to use.
  # Although this means that the original data can't use these variable names...
  tmp = df %>% mutate(
      unit := !!sym(unit),
      time := !!sym(time),
      outcome := !!sym(outcome)
  ) %>% filter(time >=( MINTIME + abs(to_rt)), time <= (MAXTIME - abs(from_rt) + 1))
  
  logger::log_info("DD FROM: {MINTIME + abs(to_rt)}")
  logger::log_info("DD TO: {MAXTIME - abs(from_rt) + 1}")
  logger::log_info("DD N: {nrow(tmp)}")
  logger::log_info("DD FROM2: {min(tmp[[time]], na.rm = T)}")
  logger::log_info("DD TO2: {max(tmp[[time]], na.rm = T)}")

  # Add in FEs
  fes = c("unit", "time", addl_fes)
  if(remove_unit_FE){
    fes = c("time", addl_fes)
  }
  fe_str = paste0(fes, collapse = ' + ')

  # Generate the formula
  if(!is.null(covariates)){
    covariate_str = paste0(covariates, collapse = ' + ')
    fmla_str = glue("{outcome} ~ {exposure} + {covariate_str} | {fe_str}")
  } else {
    fmla_str = glue("{outcome} ~ {exposure} | {fe_str}")
  }

  # if(!is.null(addl_fes)){
  #   addl_fe_str = paste0(addl_fes, collapse = ' + ')
  #   fmla_str = glue("{fmla_str} + {addl_fe_str}")
  # }

  # Estimate model with arguments
  fmla = stats::as.formula(fmla_str)
  if(!is.null(weights)){
    wt_fmla = stats::as.formula(paste0("~", weights))
    model = fixest::feols(fmla, data = tmp, cluster = ~unit, fixef.rm = 'none', weights = wt_fmla)
  } else {
    model = fixest::feols(fmla, data = tmp, cluster = ~unit, fixef.rm = 'none')
  }
  logger::log_info("DD formula: {fmla_str}")

  ct = as.data.frame(fixest::coeftable(model))
  ci = stats::confint(model)
  beta = ct[1, 1]
  se = ct[1, 2]
  tval = ct[1, 3]
  pval = ct[1, 4]
  lo95 = ci[1, 1]
  hi95 = ci[1, 2]
  dvm = fixest::fitstat(model, "my", simplify = T)
  nval = stats::nobs(model)
  logger::log_info("DD NVAL: {nval}")

  format0 = function(x, num_digits){
    return(format(x, digits = num_digits, nsmall = num_digits))
  }

  out = glue("DD beta: {format0(beta, n)} ({format0(lo95, n)} to {format0(hi95, n)}, p = {format0(pval, n)}), N={comma(nval)}, DVM = {format0(dvm, n)}")
  return(out)
  
}
  





#' Get Caption from Plot
#'
#' @param p A ggplot object
#' @return The caption for the ggplot
#' @export
#' 
get_caption_from_plot = function(p){
    plot_info <- ggplot_build(p)
    original_caption <- plot_info$plot$labels$caption
    return(original_caption)
}


#' Add Caption to Plot
#'
#' @param p A ggplot object
#' @param caption_addition The caption to add
#' @param sep The separator between the original caption and the addition
#' @return The caption for the ggplot
#' @export
#' 
add_caption_to_plot = function(p, caption_addition, sep="\n"){
    original_caption = get_caption_from_plot(p)
    if(is.null(original_caption)){
        caption = caption_addition
    } else {
        caption = glue("{original_caption}{sep}{caption_addition}")
    }
    p = p + labs(caption = caption)
    return(p)
}



#' Distributed Lags Models
#'
#' This is the distributed lags model / continuous event study. It allows for multiple outcomes.
#' @param data A data frame containing the unit, time, outcome, covariates, and any additional fixed effects
#' @param exposure_data A data frame containing the unit, time, and exposure variables
#' @param from_rt The starting lag period.
#' @param to_rt The ending lag period.
#' @param outcome The outcome variable.
#' @param exposure The exposure variable.
#' @param unit The unit identifier.
#' @param time The time variable.
#' @param covariates Vector of covariates for the model.
#' @param addl_fes Vector of additional fixed effects for the model.
#' @param ref_period Reference period (default -1)
#' @param weights Weights to be included in the regression
#' @param dd Whether to include the DD estimate in the plot
#' @param n Number of digits to round to (for the DD estimates)
#' @param dict A dictionary of variable names
#' @param remove_unit_FE Whether to remove the unit fixed effects
#' @param addl_arguments Additional arguments to be included in the model, as strings
#' @param model_type The type of model to be used (default is "feols")
#' @param plot_type The type of plot to be used (default is "ribbon")
#' @return A list containing model results, coefficients, and plots.
#' @export
#' 
distributed_lags_models = function(data, exposure_data, from_rt, to_rt, outcomes, exposure, unit, time, covariates = NULL, addl_fes = NULL, ref_period = -1, weights = NULL, dd=F, n=2, dict = NULL, remove_unit_FE = FALSE, addl_arguments = c(), model_type = "feols", plot_type = "ribbon"){

  outcomes <- c(outcomes)

  for (v in c(unit, time, outcomes, covariates, addl_fes)) {
    if (!v %in% names(data)) warning(glue("Variable {v} not found in outcome data"))
  }
  for (v in c(unit, time, exposure)) {
    if (!v %in% names(exposure_data)) warning(glue("Variable {v} not found in exposure data"))
  }
  if (!(ref_period %in% (from_rt:to_rt))) stop("ref_period must be in from_rt:to_rt")

  # Prepare data
  prep <- .prepare_exposure(data, exposure_data, exposure, unit, time)
  data <- prep$data
  exposure_data <- prep$exposure_data

  MINTIME <- min(exposure_data[[time]], na.rm = TRUE)
  MAXTIME <- max(exposure_data[[time]], na.rm = TRUE)
  logger::log_info("MINTIME: {MINTIME}, MAXTIME: {MAXTIME}")

  data_years_included <- (MINTIME + abs(to_rt)):(MAXTIME - abs(from_rt) + 1)
  data_years_included <- intersect(data_years_included, unique(data[[time]]))

  # Create leads/lags efficiently (once for all outcomes)
  logger::log_info("Creating leads and lags")
  ll <- .create_leads_lags(exposure_data, exposure, unit, time, from_rt, to_rt)
  exposure_data <- ll$exposure_data
  leads <- ll$leads
  lags <- ll$lags
  logger::log_info("Done creating leads and lags")

  # Merge (once for all outcomes)
  tmp <- merge(data, exposure_data, by = c(unit, time), all.x = TRUE)
  tmp <- tmp %>% mutate(unit := !!sym(unit), time := !!sym(time))


  .list = lapply(outcomes, function(outcome){

    res = tryCatch({
    logger::log_info("Processing {outcome}")

    # Generate the formula

    if(remove_unit_FE){
      fes = c("time", addl_fes)
    } else {
      fes = c("unit", "time", addl_fes)
    }
    fe_str = paste0(fes, collapse = ' + ')

    leads_lags_str = paste0(c(leads, lags), collapse = " + ")
    if(!is.null(covariates)){
      covariate_str = paste0(covariates, collapse = ' + ')
      fmla_str = glue("{outcome} ~ {leads_lags_str} + {covariate_str} | {fe_str}")
    } else {
      fmla_str = glue("{outcome} ~ {leads_lags_str} | {fe_str}")
    }
    
    logger::log_info("Formula:")
    print(fmla_str)
    
    # Estimate model with arguments
    arguments = c("data = tmp", "cluster = ~unit", "fixef.rm = 'none'", addl_arguments)
    if(!is.null(weights)){
      arguments = c(arguments, glue("weights = ~{weights}"))
    }
    cmd = glue("fixest::{model_type}({fmla_str}, {paste0(arguments, collapse = ', ')})")
    model = eval(parse(text=cmd))
    
    num_vars <- abs(from_rt) + to_rt
    coefficients <- model$coefficients
    if (length(coefficients) != (length(lags) + length(leads) + length(covariates))) {
      logger::log_info("Not all coefficients were estimated")
      stop("Not all coefficients were estimated")
    }
    nobs_model <- nobs(model)
    logger::log_info("DLM Model has N = {comma(nobs_model)}")

    gamma <- model$coefficients[1:num_vars]
    vcov <- vcov(model, cluster = ~unit)[1:num_vars, 1:num_vars]

    # Transform gamma -> beta
    betas <- .gamma_to_beta(gamma, vcov, from_rt, to_rt, ref_period)

    outcome_name <- "Coefficient"
    exposure_name <- "Time to Unit Change"
    if(!is.null(dict)){
      if(outcome %in% names(dict)){
        outcome_name = dict[outcome]
      } 
      if(exposure %in% names(dict)){
        exposure_name = glue("Time to Unit Change in {dict[exposure]}")
      } 
    }
    
    # This just creates a plot (adding in reference period)
    # and is not really necessary for the results
    plotdf = rbind.data.frame(betas, data.frame(time_to_event = ref_period, coef = 0, se = 0))
    plotdf = plotdf %>% mutate(
      time_to_event_str = case_when(
        time_to_event == from_rt ~ glue("{from_rt}+"),
        time_to_event == to_rt ~ glue("{to_rt}+"),
        T ~ as.character(time_to_event)
      )
    )
    plotdf = plotdf %>% arrange(time_to_event)
    p = ggplot(plotdf, aes(x = time_to_event, y = coef))
    if(tolower(plot_type) == "ribbon"){
      p = p + 
      geom_ribbon(aes(ymin = coef - 1.96 * se, ymax = coef + 1.96 * se), fill = "lightgrey", alpha = 0.6) +
      geom_line(aes(y = coef), color = "darkblue") +
      geom_point(aes(y = coef), color = "darkblue")
    } else {
      p = p + geom_line(color = "darkblue")
      p = p + geom_point(color = "darkblue")
      p = p + geom_errorbar(aes(ymin = coef - 1.96*se, ymax = coef + 1.96*se), width = 0.2, color = "darkblue")
    }
    p = p + geom_hline(yintercept = 0, linetype = "dashed")
    p = p + geom_vline(xintercept = ref_period+0.5, linetype = "dashed")
    min_included_year = min(data_years_included, na.rm = T)
    max_included_year = max(data_years_included, na.rm = T)
    cap = glue("N={comma(nobs_model)} | From {min_included_year} To {max_included_year} | {Sys.time()}")
    logger::log_info("CAPTION: {cap}")
    p = p + labs(x = exposure_name, y = outcome_name, caption = cap)
    p = p + theme_bw()

    if(dd){
      out = twfe_companion(
        data = data, 
        exposure_data = exposure_data, 
        from_rt = from_rt, 
        to_rt = to_rt, 
        outcome = outcome, 
        exposure = exposure, 
        unit = unit, 
        time = time, 
        covariates = covariates, 
        addl_fes = addl_fes, 
        ref_period = ref_period, 
        weights = weights, 
        dd = dd, 
        n = n
      )
      # out = do.call(twfe_companion, list(...))
      p = add_caption_to_plot(p, out)
      logger::log_info("ADDED TO CAPTION: {out}")
      # p = p + labs(caption = out)
    }

    list(
      "betas" = betas, 
      "plot" = p, 
      "model" = model, 
      "vcov" = vcov, 
      "data_years_included" = data_years_included, 
      "fmla_str" = fmla_str, 
      "from_rt" = from_rt, 
      "to_rt" = to_rt, 
      "cmd" = cmd,
      "exposure" = exposure,
      "outcome" = outcome
    )
    

    }, 
    error = function(cond){
      logger::log_info("Error in outcome {outcome}")
      logger::log_info("Returning null")
      print(cond)
      return(NULL)
    })

    return(res)
    
  })

  .list = .list[!sapply(.list, is.null)]
  
  # if(length(.list) == 0){
  #   logger::log_error("No outcomes were estimated")
  #   return(NULL)
  # } else if(length(.list) == 1){
  #   return(.list[[1]])
  # } else {
  #   return(.list)
  # }


  if(length(.list) == 0){
    logger::log_error("No outcomes were estimated")
    return(NULL)
  } else {
    return(.list)
  }


}





#' Distributed Lags Models v2
#'
#' This is experimental but may be faster for multiple outcomes.
#' @param data A data frame containing the unit, time, outcome, covariates, and any additional fixed effects
#' @param exposure_data A data frame containing the unit, time, and exposure variables
#' @param from_rt The starting lag period.
#' @param to_rt The ending lag period.
#' @param outcome The outcome variable.
#' @param exposure The exposure variable.
#' @param unit The unit identifier.
#' @param time The time variable.
#' @param covariates Vector of covariates for the model.
#' @param addl_fes Vector of additional fixed effects for the model.
#' @param ref_period Reference period (default -1)
#' @param weights Weights to be included in the regression
#' @param dd Whether to include the DD estimate in the plot
#' @param n Number of digits to round to (for the DD estimates)
#' @param dict A dictionary of variable names
#' @param remove_unit_FE Whether to remove the unit fixed effects
#' @param addl_arguments Additional arguments to be included in the model, as strings
#' @return A list containing model results, coefficients, and plots.
#' @export
#' 
distributed_lags_models2 = function(data, exposure_data, from_rt, to_rt, outcomes, exposure, unit, time, covariates = NULL, addl_fes = NULL, ref_period = -1, weights = NULL, dd=F, n=2, dict = NULL, remove_unit_FE = FALSE, addl_arguments = c(), model_type = "feols", plot_type = "ribbon"){

  outcomes <- c(outcomes)

  # Filter to valid outcomes
  outcomes <- outcomes[outcomes %in% names(data)]
  outcomes <- outcomes[sapply(outcomes, function(x) !all(is.na(data[[x]])))]
  if (length(outcomes) == 0) {
    logger::log_error("No outcomes were specified or all outcomes are NA")
    return(NULL)
  }

  for (v in c(unit, time, outcomes, covariates, addl_fes)) {
    if (!v %in% names(data)) warning(glue("Variable {v} not found in outcome data"))
  }
  for (v in c(unit, time, exposure)) {
    if (!v %in% names(exposure_data)) warning(glue("Variable {v} not found in exposure data"))
  }
  if (!(ref_period %in% (from_rt:to_rt))) stop("ref_period must be in from_rt:to_rt")

  # Prepare data
  prep <- .prepare_exposure(data, exposure_data, exposure, unit, time)
  data <- prep$data
  exposure_data <- prep$exposure_data
  data <- data %>% select(any_of(c(unit, time, outcomes, exposure, covariates, addl_fes, weights)))

  logger::log_info("========== RUNNING DISTRIBUTED LAGS MODELS ===========")
  logger::log_info("Data N: {nrow(data)}, Exposure N: {nrow(exposure_data)}")
  logger::log_info("from_rt: {from_rt}, to_rt: {to_rt}")
  logger::log_info("Outcomes: {paste(outcomes, collapse = ', ')}")

  MINTIME <- min(exposure_data[[time]], na.rm = TRUE)
  MAXTIME <- max(exposure_data[[time]], na.rm = TRUE)
  logger::log_info("MINTIME: {MINTIME}, MAXTIME: {MAXTIME}")

  data_years_included <- (MINTIME + abs(to_rt)):(MAXTIME - abs(from_rt) + 1)
  data_years_included <- intersect(data_years_included, unique(data[[time]]))

  # Create leads/lags efficiently (once for all outcomes)
  logger::log_info("Creating leads and lags")
  ll <- .create_leads_lags(exposure_data, exposure, unit, time, from_rt, to_rt)
  exposure_data <- ll$exposure_data
  leads <- ll$leads
  lags <- ll$lags
  logger::log_info("Done creating leads and lags")

  # Merge (once for all outcomes)
  tmp <- merge(data, exposure_data, by = c(unit, time), all.x = TRUE)
  tmp <- tmp %>% mutate(unit := !!sym(unit), time := !!sym(time))
  
  
  
  if(remove_unit_FE){
    fes = c("time", addl_fes)
  } else {
    fes = c("unit", "time", addl_fes)
  }
  fe_str = paste0(fes, collapse = ' + ')
  
  leads_lags_str = paste0(c(leads, lags), collapse = " + ")
  if(!is.null(covariates)){
    covariate_str = paste0(covariates, collapse = ' + ')
    fmla_str = glue(".[outcomes] ~ {leads_lags_str} + {covariate_str} | {fe_str}")
  } else {
    fmla_str = glue(".[outcomes] ~ {leads_lags_str} | {fe_str}")
  }

  logger::log_info("Exposure: {exposure}")
  logger::log_info("Formula:")
  print(fmla_str)
  
  # Estimate model with arguments
  arguments = c("data = tmp", "cluster = ~unit", "fixef.rm = 'none'", addl_arguments)
  if(!is.null(weights)){
    arguments = c(arguments, glue("weights = ~{weights}"))
  }
  cmd = glue("fixest::{model_type}({fmla_str}, {paste0(arguments, collapse = ', ')})")
  models = eval(parse(text=cmd))
  
  if(length(outcomes) == 1){
    models = list(models)
  } 

  logger::log_info("length(models): {length(models)}")
  
  .list = lapply(models, function(model){
    
    logger::log_info("Capturing outcome")
    outcome <- as.character(formula(model))[2]
    
    res = tryCatch({
      logger::log_info("Processing {outcome}")
      
      # Generate the formula
      
      if(remove_unit_FE){
        fes = c("time", addl_fes)
      } else {
        fes = c("unit", "time", addl_fes)
      }
      fe_str = paste0(fes, collapse = ' + ')
      
      leads_lags_str = paste0(c(leads, lags), collapse = " + ")
      if(!is.null(covariates)){
        covariate_str = paste0(covariates, collapse = ' + ')
        fmla_str = glue("{outcome} ~ {leads_lags_str} + {covariate_str} | {fe_str}")
      } else {
        fmla_str = glue("{outcome} ~ {leads_lags_str} | {fe_str}")
      }
      
      logger::log_info("Outcome: {outcome}")
      logger::log_info("Formula:")
      print(fmla_str)
      
      # Estimate model with arguments
      arguments = c("data = tmp", "cluster = ~unit", "fixef.rm = 'none'", addl_arguments)
      if(!is.null(weights)){
        arguments = c(arguments, glue("weights = ~{weights}"))
      }
      cmd = glue("fixest::{model_type}({fmla_str}, {paste0(arguments, collapse = ', ')})")
      # model = eval(parse(text=cmd))
      
      num_vars <- abs(from_rt) + to_rt
      coefficients <- model$coefficients
      if (length(coefficients) != (length(lags) + length(leads) + length(covariates))) {
        logger::log_info("Not all coefficients were estimated")
        stop("Not all coefficients were estimated")
      }
      nobs_model <- nobs(model)
      logger::log_info("DLM Model has N = {comma(nobs_model)}")

      gamma <- model$coefficients[1:num_vars]
      vcov <- vcov(model, cluster = ~unit)[1:num_vars, 1:num_vars]

      # Transform gamma -> beta
      betas <- .gamma_to_beta(gamma, vcov, from_rt, to_rt, ref_period)
      
      
      outcome_name = "Coefficient"
      exposure_name = "Time to Unit Change"
      if(!is.null(dict)){
        if(outcome %in% names(dict)){
          outcome_name = dict[outcome]
        } 
        if(exposure %in% names(dict)){
          exposure_name = glue("Time to Unit Change in {dict[exposure]}")
        } 
      }
      
      # This just creates a plot (adding in reference period)
      # and is not really necessary for the results
      plotdf = rbind.data.frame(betas, data.frame(time_to_event = ref_period, coef = 0, se = 0))
      plotdf = plotdf %>% mutate(
        time_to_event_str = case_when(
          time_to_event == from_rt ~ glue("{from_rt}+"),
          time_to_event == to_rt ~ glue("{to_rt}+"),
          T ~ as.character(time_to_event)
        )
      )
      plotdf = plotdf %>% arrange(time_to_event)
      p = ggplot(plotdf, aes(x = time_to_event, y = coef))
      if(tolower(plot_type) == "ribbon"){
        p = p + 
        geom_ribbon(aes(ymin = coef - 1.96 * se, ymax = coef + 1.96 * se), fill = "lightgrey", alpha = 0.6) +
        geom_line(aes(y = coef), color = "darkblue") +
        geom_point(aes(y = coef), color = "darkblue")
      } else {
        p = p + geom_line(color = "darkblue")
        p = p + geom_point(color = "darkblue")
        p = p + geom_errorbar(aes(ymin = coef - 1.96*se, ymax = coef + 1.96*se), width = 0.2, color = "darkblue")
      }
      p = p + geom_hline(yintercept = 0, linetype = "dashed")
      p = p + geom_vline(xintercept = ref_period+0.5, linetype = "dashed")
      min_included_year = min(data_years_included, na.rm = T)
      max_included_year = max(data_years_included, na.rm = T)
      cap = glue("N={comma(nobs_model)} | From {min_included_year} To {max_included_year} | {Sys.time()}")
      logger::log_info("CAPTION: {cap}")
      p = p + labs(x = exposure_name, y = outcome_name, caption = cap)
      p = p + theme_bw()
      
      if(dd){
        out = twfe_companion(
          data = data, 
          exposure_data = exposure_data, 
          from_rt = from_rt, 
          to_rt = to_rt, 
          outcome = outcome, 
          exposure = exposure, 
          unit = unit, 
          time = time, 
          covariates = covariates, 
          addl_fes = addl_fes, 
          ref_period = ref_period, 
          weights = weights, 
          dd = dd, 
          n = n
        )
        # out = do.call(twfe_companion, list(...))
        p = add_caption_to_plot(p, out)
        logger::log_info("ADDED TO CAPTION: {out}")
        # p = p + labs(caption = out)
      }
      
      list(
        "betas" = betas, 
        "plot" = p, 
        "model" = model, 
        "vcov" = vcov, 
        "data_years_included" = data_years_included, 
        "fmla_str" = fmla_str, 
        "from_rt" = from_rt, 
        "to_rt" = to_rt, 
        "cmd" = cmd,
        "exposure" = exposure,
        "outcome" = outcome
      )
      
      
    }, 
    error = function(cond){
      logger::log_info("Error in outcome {outcome}")
      logger::log_info("Returning null")
      print(cond)
      return(NULL)
    })
    
    gc()
    return(res)
    
  })
  
  .list = .list[!sapply(.list, is.null)]
  
  # if(length(.list) == 0){
  #   logger::log_error("No outcomes were estimated")
  #   return(NULL)
  # } else if(length(.list) == 1){
  #   return(.list[[1]])
  # } else {
  #   return(.list)
  # }
  
  
  if(length(.list) == 0){
    logger::log_error("No outcomes were estimated")
    return(NULL)
  } else {
    return(.list)
  }
  
  
}

