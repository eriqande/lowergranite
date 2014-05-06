\name{slurp_all_results}
\alias{slurp_all_results}
\title{Read all the results from all the runs with output in Dirs into a named list}
\usage{
slurp_all_results(Dirs)
}
\arguments{
\item{Dirs}{a character vector of directories to slurp results from. Ideally
you will just give the basenae of them all as these will become the names
of the output list.}
}
\description{
Read all the results from all the runs with output in Dirs into a named list
}
\examples{
\dontrun{
setwd("~/prj/LowerGranite/runs/lowergranite_big_run_CincoDeMayo_2014/")
big_run <- slurp_all_results(dir(pattern = "^[CS].*"))
}
}

