\name{run_boot_gsi_analysis}
\alias{run_boot_gsi_analysis}
\title{Run the bootstrap analysis, with gsi module if specified in parameters}
\usage{
run_boot_gsi_analysis(W = NULL, stock_group_start_col,
  DAT.DIR = system.file("data_files", package = "lowergranite", mustWork = T),
  WORK.DIR = getwd(), STOCK.DATA.XLSX = file.path(DAT.DIR,
  "SH11SIMPOP_StockSex.xlsx"), drop.these.groups = NULL, collaps = c(1, 1,
  1, 1, 1, 2, 2, 2, 3, 4, 5, 6, 7, 8, 9, 10, 10, 10, 10, 11, 11, 11, 11, 11, 11,
  11, 11), DO_GSI_ON_PROP = FALSE, GSISIM = gsi_simBinaryPath(),
  GSI_SEEDS = c(NA, NA), BLFILE = file.path(DAT.DIR,
  "sthd_base_v3_187.txt"), RUFILE = file.path(DAT.DIR, "sthd_base_v3_rg.txt"),
  alph = 0.1, B = 5, nsim = 2, console_messages_to = "",
  reset_booty_seed = 0, GroupMin = 0, Run = "2011 Steelhead Stock")
}
\arguments{
\item{W}{a data frame that holds all the stock data in it used to drive the simulation.  This is the data
frame that we typically get by reading in the .xlsx file named in STOCK.DATA.XLSX.  However, we can also
just pass a data frame in directly.  By default W is NULL, in which case the stock data will be obtained
by reading in the STOCK.DATA.XLSX file.  If W is not NULL, then W will be used instead
of STOCK.DATA.XLSX.  If both W and STOCK.DATA.XLSX are NULL, that is an error.}

\item{stock_group_start_col}{The column at which the stock groups start in the xlsx file (or, equivalently,
in the data frame W).  That data frame, (whether it was passed in as W or read from an xlsx file) must have column
names that correspond to the groupings of stocks that you want to do the bootstrapping for.  (For example,
"UPSALM", "MFSALM", "SFSALM", ... ).  The column at which those names start must be passed to this function as
\code{stock_group_start_col}. You must specify a value for this.  There is no default.  \emph{Note:} there must be no other columns in the data frame after the stock group
columns.}

\item{DAT.DIR}{The directory where all the data files are.  Defaults to the directory
"data_files" in the installed package}

\item{WORK.DIR}{The working directory to do this in.  Default = current working directory.
Note that gsi_sim will also be run in this directory.}

\item{STOCK.DATA.XLSX}{the path of the file that has the stock data in it used to drive the simulations. This can be NULL, in which
case parameter \code{W} must be specified.}

\item{drop.these.groups}{A character vector of the names of the stocks or stock-by-age or stock-by-sex groups that you want
to drop from the analysis (typically because they are at such low numbers that there are bootstrap reps when none of
them occur).  Must follow the convention of column naming in the file.  For example \code{c("UPSALM..BY04", "MFSALM..BY08")}.
If this is non-null then only columns not matching any of the entries in this character vector will be retained.}

\item{collaps}{A vector of numbers in 1,...,N telling which weeks should be lumped together into "statistical weeks" from Kirk's code. It looks
like this gets used in a lot of the bootstrapping functions, but is not a formal parameter of the bootstrapping functions.}

\item{DO_GSI_ON_PROP}{if set to TRUE then gsi_sim is used to create assignments that replace the assignments in the variable Prop.  If FALSE
then the true origins are used.}

\item{GSISIM}{path to the gsi_sim executable.}

\item{GSI_SEEDS}{vector of two positive integers that will be written to the
gsi_sim_seeds file to make reproducible results. If any elements of the vector
are NA, then gsi_sim_seeds is not modified.}

\item{BLFILE}{path to the gsi_sim baseline file.}

\item{RUFILE}{path to the gsi_sim reporting group file.}

\item{alph}{One minus the size of the desired confidence interval to be calculated}

\item{B}{number of bootstrap replicates. The default is 5 --- much lower than it should be
(should be more like 500)
because it takes a long time and this is better for testing.}

\item{nsim}{number of simulation replicates to do. The default is 2 --- much lower than it should be
(should be more like 500)
because it takes a long time and this is better for testing}

\item{console_messages_to}{path to a file you want the console messages written to.  Note that it will
always append these to a file.  Default is "" which means send it to the console.}

\item{reset_booty_seed}{for some reason, calling gsi_sim seems to get the random number generator out
of state in a way that cannot be restored by saving .Random.seed and then setting that value back to
itself.  It is odd and vexing.  Anyway, in order to test that comparable results are obtained with
the super-informo gsi data, we have this.  It should be an integer.  If >0 then it will be passed
to \code{\link{set.seed()}} after the gsi code has been run (before entering the bootstrap loop.)
For an example of its use, see the test files.}
}
\value{
returns the \code{sumrys} object from Kirk's script.
}
\description{
This is basically Kirk's script that eric has wrapped up into a function, and made
some significant changes to.  Notably, the interface to gsi_sim has been rewritten
to make it much more time efficient.  The simulated "true" data sets are all simulated
first and stored in a list, and then, if doing the gsi_sim part, all of that is passed
to gsi_sim in one fell swoop and the Groop of each fish is replaced by its "gsi-inferred"
Groop.
}
\details{
Note that the default values here are set up to do an analysis of the steelhead data,
and, by default, to not use gsi_sim assigments.
}
\note{
Running gsi_sim this produced warnings that look like this:
\code{Warning message:}
\code{In system2(GSISIM, args = gsi.args, stdout = T) :}
\code{line 98 may be truncated in call to system(, intern = TRUE)}.
This was \emph{not} a problem.
The cause of this is that the output file from gsi_sim has a line it
that shows that the
command line looked like after all the --multi-fix-mix commands
were stuck onto it using
gsi_sim's --command-file option.  That line is so large that system()
chops it into two
(or many, if you are doing large nsim).  This line is
far away from any of the important
lines we will grep out, however, so it is not of any consequence.
I wrapped it in a \code{\link{suppressWarnings}} to not bark too much.

Also note that this function still writes stuff out to an xlsx file of a name
that is not specified in the parameters yet.  It will probably be better
to just return it as a data frame or matrix anyway.
}
\examples{
# Do a very short run with known stock of origin:
set.seed(5)
known_stock_result1 <- run_boot_gsi_analysis(stock_group_start_col = 9, nsim = 10, B = 50, DO_GSI_ON_PROP = F)

# Do a short run using the gsi assignments
gsi_result1 <- run_boot_gsi_analysis(stock_group_start_col = 9, nsim = 5, B = 10, DO_GSI_ON_PROP = T)
}

