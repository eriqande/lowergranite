\name{check_subgroup_level_completeness}
\alias{check_subgroup_level_completeness}
\title{Test that all STOCK..CHARACTERISTIC are present.}
\usage{
check_subgroup_level_completeness(Onames)
}
\arguments{
  \item{Onames}{the names of the different origin groups.
  This will be Originnames, typically.}
}
\description{
When using the ".." notation to define subroups of stocks (for example stock-by-sex or
stock-by-age) it is required that any subgroup designations be applied to all stocks.
This function just checks that this is the case.  If not, it throws an error.  Of course, if none of the
stock group designations have a ".." in them, then there is no problem.
}

