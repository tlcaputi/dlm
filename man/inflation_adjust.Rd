% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cpi-functions.R
\name{inflation_adjust}
\alias{inflation_adjust}
\title{Inflation Adjust}
\usage{
inflation_adjust(
  df,
  old_price_var,
  new_price_var,
  year_var = "year",
  to_year = 2010
)
}
\arguments{
\item{df}{Data frame}

\item{old_price_var}{Old price variable name}

\item{new_price_var}{New price variable name}

\item{year_var}{Year variable name}

\item{to_year}{Year to which to adjust}
}
\value{
Dataframe
}
\description{
Gives the reverse cumulative sum of a vector.
}
