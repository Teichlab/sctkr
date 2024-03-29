% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/CellTypeCompositionAnalysis.R
\name{plot_ranef}
\alias{plot_ranef}
\title{Plot estimated random effects of variables}
\usage{
plot_ranef(ranef_tbl, vars = NULL, celltypes = NULL, maxFC = 3)
}
\arguments{
\item{ranef_tbl}{Table of random effect estimates (mean and sd) extracted from fitted model
(data.frame-like object)}

\item{vars}{A list specifying required categorical sample metadata variables and their order,
e.g. list(size=c('small', 'medium', 'big'), price=c('cheap', 'expensive'))}

\item{celltypes}{If specified, only keep those cell types (vector of str)}

\item{maxFC}{Cap fold change by this value (numeric, default 3)}
}
\value{
ggplot object of a dot plot
}
\description{
This function plots random effects of variables on cell type composition as a dot plot
}
