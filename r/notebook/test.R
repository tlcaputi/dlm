# roxygen2::roxygenise()

# Load package


# rm(distributed_lags_models)
# devtools::install_github("tlcaputi/dlm")

# devtools::install_github("tlcaputi/dlm@backup-master-before-rewind", force = TRUE)
devtools::install_github("tlcaputi/dlm", force = TRUE)


rm(list = ls())


distributed_lags_models2 = function(data, exposure_data, from_rt, to_rt, outcomes, exposure, unit, time, covariates = NULL, addl_fes = NULL, ref_period = -1, weights = NULL, dd=F, n=2, dict = NULL, remove_unit_FE = FALSE, addl_arguments = c(), model_type = "feols"){
  
  # if outcomes is a string, make it a vector of that string
  outcomes = c(outcomes)
  
  for(v in c(unit, time, outcomes, covariates, addl_fes)){
    if(!v %in% names(data)){
      warning(glue("Variable {v} not found in outcome data"))
    }
  }
  
  for(v in c(unit, time, exposure)){
    if(!v %in% names(exposure_data)){
      warning(glue("Variable {v} not found in outcome data"))
    }
  }
  
  if(!(ref_period %in% (from_rt:to_rt))){
    stop("ref_period must be in from_rt:to_rt")
  }
  
  
  # Make sure that the outcome data doesn't have the exposure in it, which would cause problems in the merges
  try({
    data = data %>% select(-!!sym(exposure))
  })
  
  try({
    data = data %>% filter(!is.na(!!sym(unit)), !is.na(!!sym(time)))
  })
  
  try({
    exposure_data = exposure_data %>% filter(!is.na(!!sym(unit)), !is.na(!!sym(time)), !is.na(!!sym(exposure)))
  })
  
  try({
    exposure_data = exposure_data %>% select(!!sym(unit), !!sym(time), !!sym(exposure)) %>% unique()
  })
  
  
  if(!plm::is.pbalanced(exposure_data, index = c(unit, time))){
    warning("

    ****** WARNING!!!! ******
    The exposure data is not balanced after taking unique values of unit, time, and exposure.
    We are going to fill in the missing values with NAs. You should make sure this is okay.
    
    ")
    
    min_time = min(exposure_data[[time]], na.rm = T)
    max_time = max(exposure_data[[time]], na.rm = T)
    cmd = glue("
      exposure_data = exposure_data %>% tidyr::complete({unit}, {time}={min_time}:{max_time})
    ")
    eval(parse(text=cmd))
    # exposure_data = exposure_data %>% tidyr::complete(!!sym(unit), !!sym(time))
    
    if(!plm::is.pbalanced(exposure_data, index = c(unit, time))){
      stop("The exposure data is still not balanced after taking unique values of unit, time, and exposure.")
    }
    
  }
  
  
  
  # Capture the minimum and maximum time
  MINTIME = min(exposure_data[[time]], na.rm = T)
  MAXTIME = max(exposure_data[[time]], na.rm = T)
  log_info("MINTIME: {MINTIME}")
  log_info("MAXTIME: {MAXTIME}")
  
  # This tells us the data years that are included
  data_years_included = (MINTIME + abs(to_rt)):(MAXTIME - abs(from_rt) + 1)
  log_info("Data years that could be included, based upon the exposure data")
  print(data_years_included)
  data_years_included = intersect(data_years_included, unique(data[[time]]))
  log_info("Data years that actually are included, based upon the exposure and outcome data")
  print(data_years_included)
  
  # Create lag and lead variables in the distance data
  ## TODO: we only really need to create these once, not every time. Not sure how much time that saves.
  log_info("Creating leads and lags")
  leads = c()
  
  if(abs(from_rt) > 1){
    for(i in (abs(from_rt) - 1):1){
      exposure_data = exposure_data %>% 
        group_by(!!sym(unit)) %>% 
        mutate(
          "{exposure}_lead{i}" := dplyr::lead(!!sym(exposure), i)
        )
      leads = c(leads, glue("{exposure}_lead{i}"))
    }
  } 
  leads_str = paste0(leads, collapse = " + ")
  
  lags = c()
  for(i in 0:to_rt){
    exposure_data = exposure_data %>% 
      group_by(!!sym(unit)) %>% 
      mutate(
        "{exposure}_lag{i}" := dplyr::lag(!!sym(exposure), i)
      )
    lags = c(lags, glue("{exposure}_lag{i}"))
  }
  lags_str = paste0(lags, collapse = " + ")
  log_info("Done creating leads and lags")
  
  # Merge together the data and distance data
  # data = setDT(data); exposure_data = setDT(exposure_data)
  tmp = merge(data, exposure_data, by=c(unit, time), all.x = T)
  
  # Define unit and time
  tmp = tmp %>% mutate(unit := !!sym(unit), time := !!sym(time))
  
  
  
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
  
  log_info("Formula:")
  print(fmla_str)
  
  # Estimate model with arguments
  arguments = c("data = tmp", "cluster = ~unit", "fixef.rm = 'none'", addl_arguments)
  if(!is.null(weights)){
    arguments = c(arguments, glue("weights = ~{weights}"))
  }
  cmd = glue("fixest::{model_type}({fmla_str}, {paste0(arguments, collapse = ', ')})")
  log_info("CMD:")
  print(cmd)
  models = eval(parse(text=cmd))
  print(str(models))
  
  log_info("length(models): {length(models)}")
  
  .list = lapply(models, function(model){
    
    log_info("Capturing outcome")
    outcome <- as.character(formula(model))[2]
    
    res = tryCatch({
      log_info("Processing {outcome}")
      
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
      
      log_info("Formula:")
      print(fmla_str)
      
      # Estimate model with arguments
      arguments = c("data = tmp", "cluster = ~unit", "fixef.rm = 'none'", addl_arguments)
      if(!is.null(weights)){
        arguments = c(arguments, glue("weights = ~{weights}"))
      }
      cmd = glue("fixest::{model_type}({fmla_str}, {paste0(arguments, collapse = ', ')})")
      # model = eval(parse(text=cmd))
      
      log_info("Coefficients:")
      num_vars = abs(from_rt)+to_rt
      coefficients = model$coefficients
      if(length(coefficients) != (length(lags) + length(leads) + length(covariates))){
        log_info("Not all coefficients were estimated")
        print(coefficients)
        stop("Not all coefficients were estimated")
      }
      print(model$coefficients[1:num_vars])
      nobs_model = nobs(model)
      log_info("DLM Model has N = {comma(nobs_model)}")
      # Extract coefficients and regressions from the model
      gamma = model$coefficients[1:num_vars]
      vcov = vcov(model, cluster= ~unit)[1:num_vars, 1:num_vars]
      
      # Sum them up to the reference period
      
      if(from_rt == ref_period){
        
        time_to_event = sort(setdiff(from_rt:to_rt, c(ref_period)))
        after_periods = 1:num_vars
        log_info("After periods:")
        print(after_periods)
        log_info("vcov")
        print(vcov)
        
        coefs = c(
          cumsum(gamma[after_periods])
        )
        ses = c(
          secumsum(as.matrix(vcov)[after_periods, after_periods])
        )
        betas = data.frame(
          time_to_event = time_to_event,
          coef = coefs,
          se = ses
        )
        
      } else if(from_rt < ref_period) {
        
        time_to_event = sort(setdiff(from_rt:to_rt, c(ref_period)))
        num_before_periods = length(time_to_event[time_to_event < ref_period])
        before_periods = 1:num_before_periods
        after_periods = (num_before_periods + 1):num_vars
        
        coefs = c(
          -revcumsum(gamma[before_periods]),
          cumsum(gamma[after_periods])
        )
        ses = c(
          serevcumsum(vcov[before_periods, before_periods]),
          secumsum(vcov[after_periods, after_periods])
        )
        betas = data.frame(
          time_to_event = time_to_event,
          coef = coefs,
          se = ses
        )
      }
      
      
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
      p = p + geom_line(color = "darkblue")
      p = p + geom_point(color = "darkblue")
      p = p + geom_errorbar(aes(ymin = coef - 1.96*se, ymax = coef + 1.96*se), width = 0.2, color = "darkblue")
      p = p + geom_hline(yintercept = 0, linetype = "dashed")
      p = p + geom_vline(xintercept = ref_period+0.5, linetype = "dashed")
      min_included_year = min(data_years_included, na.rm = T)
      max_included_year = max(data_years_included, na.rm = T)
      cap = glue("N={comma(nobs_model)} | From {min_included_year} To {max_included_year} | {Sys.time()}")
      log_info("CAPTION: {cap}")
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
        log_info("ADDED TO CAPTION: {out}")
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
      log_info("Error in outcome {outcome}")
      log_info("Returning null")
      print(cond)
      return(NULL)
    })
    
    return(res)
    
  })
  
  .list = .list[!sapply(.list, is.null)]
  
  # if(length(.list) == 0){
  #   log_error("No outcomes were estimated")
  #   return(NULL)
  # } else if(length(.list) == 1){
  #   return(.list[[1]])
  # } else {
  #   return(.list)
  # }
  
  
  if(length(.list) == 0){
    log_error("No outcomes were estimated")
    return(NULL)
  } else {
    return(.list)
  }
  
  
}


library(dlm)
library(dplyr)
# devtools::load_all(".")

# Generate Data
raw = generate_data()
raw = generate_data(seed = 12345, n_groups = 40^2, n_times = 400, treat_prob = 0.1)
raw = raw %>% mutate(outcome2 = outcome - 1 + rnorm(nrow(raw), 0, 1))

# Separate data
data = raw %>% select(-post)
exposure_data = raw %>% select(group, time, post) %>% unique()
outcomes = c("outcome")


mods2 = distributed_lags_models2(
  data = data, 
  exposure_data = exposure_data, 
  from_rt = -3, 
  to_rt = 3, 
  remove_unit_FE = F,
  outcome = outcomes, 
  exposure = "post", 
  unit = "group", 
  time = "time",
  dd = T,
  n = 4
)


# Distributed Lags Model
mods1 = distributed_lags_models(
    data = data, 
    exposure_data = exposure_data, 
    from_rt = -3, 
    to_rt = 3, 
    remove_unit_FE = F,
    outcome = outcomes, 
    exposure = "post", 
    unit = "group", 
    time = "time",
    dd = T,
    n = 4
)



mods1[[1]]$betas
mods2[[1]]$betas


mods1[[2]]$betas
mods2[[2]]$betas

rm(list = ls())

library(dlm)
library(dplyr)
# devtools::load_all(".")

# Generate Data
raw = generate_data()
raw = generate_data(seed = 12345, n_groups = 40^2, n_times = 400, treat_prob = 0.1)
raw = raw %>% mutate(outcome2 = outcome - 1 + rnorm(nrow(raw), 0, 1))

# Separate data
data = raw %>% select(-post)
exposure_data = raw %>% select(group, time, post) %>% unique()
outcomes = c("outcome")

# Distributed Lags Model
mods1 = distributed_lags_models(
    data = data, 
    exposure_data = exposure_data, 
    from_rt = -3, 
    to_rt = 3, 
    remove_unit_FE = F,
    outcome = outcomes, 
    exposure = "post", 
    unit = "group", 
    time = "time",
    dd = T,
    n = 4
)


mods2 = distributed_lags_models2(
  data = data, 
  exposure_data = exposure_data, 
  from_rt = -3, 
  to_rt = 3, 
  remove_unit_FE = F,
  outcome = outcomes, 
  exposure = "post", 
  unit = "group", 
  time = "time",
  dd = T,
  n = 4
)


mods1[[1]]$betas
mods2[[1]]$betas


mods1[[2]]$betas
mods2[[2]]$betas


length(mods1)
mods1$betas
mods2$betas

# 
# mod$betas
# # mod$plot
# # str(mod)
# 
# 
# # TWFE model for comparison
# comp = standard_twfe_for_comparison(
#     data = raw,
#     from_rt = -3,
#     to_rt = 3,
#     outcome = "outcome",
#     unit = "group",
#     time = "time",
#     time_to_treatment = "years_to_treatment",
#     treat = "treat"
# )
# 
# 
# # Should be the same
# print(mod$betas)
# print(comp$betas)
# 
# all(abs(mod$betas$coef - comp$betas$coef) < 0.0001)
# 
# 
# 
# mod$plot
