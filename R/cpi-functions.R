#' Inflation Adjust
#'
#' Gives the reverse cumulative sum of a vector.
#' @param df Data frame
#' @param old_price_var Old price variable name
#' @param new_price_var New price variable name
#' @param year_var Year variable name
#' @param to_year Year to which to adjust
#' @return Dataframe 
#' @export
#' 
inflation_adjust = function(df, old_price_var, new_price_var, year_var, to_year){
  data(cpi)
  baseline_val = cpi$cpi[cpi$year == to_year]
  cpi = cpi %>% 
      mutate(adjustment_factor = cpi / baseline_val) %>%
      select(
          "{year_var}" := year, 
          adjustment_factor
      ) %>%
      unique()
  df = setDT(df); cpi = setDT(cpi)
  df = merge(df, cpi, by = year_var, all.x = T)
  df = as.data.frame(df)
  df = df %>% 
    mutate(
        "{new_price_var}" := !!sym(old_price_var) / adjustment_factor
    ) %>% 
    select(-adjustment_factor)

  return(df)

}

