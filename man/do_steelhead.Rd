\name{do_steelhead}
\alias{do_steelhead}
\title{Do the steelhead stock analysis in the current working directory}
\usage{
do_steelhead(STOCK.DATA.XLSX, stock_group_start_col, nsim, B, DO_GSI_ON_PROP)
}
\arguments{
  \item{STOCK.DATA.XLSX}{The basename of the file with the
  data in it for running the simulation The function
  assumes that the file is in the "data_files" directory of
  the package.}

  \item{stock_group_start_col}{The column at which the
  stock groups start in the xlsx file (or, equivalently, in
  the data frame W).  That data frame, (whether it was
  passed in as W or read from an xlsx file) must have
  column names that correspond to the groupings of stocks
  that you want to do the bootstrapping for.  (For example,
  "UPSALM", "MFSALM", "SFSALM", ... ).  The column at which
  those names start must be passed to this function as
  \code{stock_group_start_col}. You must specify a value
  for this.  There is no default.  \emph{Note:} there must
  be no other columns in the data frame after the stock
  group columns.}

  \item{nsim}{number of simulation replicates to do. The
  default is 2 --- much lower than it should be (should be
  more like 500) because it takes a long time and this is
  better for testing}

  \item{B}{number of bootstrap replicates. The default is 5
  --- much lower than it should be (should be more like
  500) because it takes a long time and this is better for
  testing.}

  \item{DO_GSI_ON_PROP}{if set to TRUE then gsi_sim is used
  to create assignments that replace the assignments in the
  variable Prop.  If FALSE then the true origins are used.}
}
\description{
Note that you have to be sure that stock_group_start_col is appropriate for the
input file.  The input file can be modified to do stock-sex or stock-age if desired.
}

