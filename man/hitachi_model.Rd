% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/hitachi_model.R
\name{hitachi_model}
\alias{hitachi_model}
\title{hitachi_model}
\usage{
hitachi_model(data)
}
\arguments{
\item{data}{tibble with column S (Substrate concentration), v (initial velocity) and std.error (of the slope of the lm)}
}
\description{
This function fits models model.nls1 and model.nls2 to initial velocity in dependence of substrate concentration and provides a grid for fitting.
}
\examples{
hitachi_model()
}
