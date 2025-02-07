% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rfhelps.R
\name{build_data}
\alias{build_data}
\title{Build model data}
\usage{
build_data(
  years,
  ages,
  lengths,
  weight_array,
  maturity_array,
  survey_df,
  landings_df,
  bin_adjust = 0.5,
  min_size = 22,
  use_discards = FALSE,
  inf_length = 100,
  Q_prior_max = inf_length
)
}
\arguments{
\item{years}{vector of years to include}

\item{ages}{vector of ages to use in the model}

\item{lengths}{vector of lengths to include}

\item{weight_array}{array of weights with dimnames for male/female, year, length, etc.}

\item{maturity_array}{array of maturity ogives to use, with appropriate dimnames}

\item{survey_df}{data.frame of survey data to use}

\item{landings_df}{data.frame of landings data, should include landings bounds}

\item{min_size}{minimum required landing size}

\item{use_discards}{use discards or not}

\item{inf_length}{How far out to push the length matrix for the length plus group?}

\item{Q_prior_max}{Q prior max length}

\item{bin_ajdust}{midpoint of the length bins}
}
\description{
Setup the model data
}
