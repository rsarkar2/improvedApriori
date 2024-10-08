% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/improvedApriori.R
\docType{_PACKAGE}
\name{improvedApriori}
\alias{improvedApriori}
\title{R Package for implementing improvedApriori to identify frequently co-occurring drug/vaccine and AE pair}
\usage{
improvedApriori(
  data = data,
  numdrug = 5,
  numae = 13,
  numreports = 13,
  second_par = "ROR",
  minsup = 3,
  minthresh = 1
)
}
\arguments{
\item{data}{The binarized drug/vaccine safety data. Data should be saved as a data frame.}

\item{numdrug}{The number of unique drugs/vaccines. Default value is 5.}

\item{numae}{The number of unique adverse events. Default value is 13.}

\item{numreports}{The number of unique report IDs. Default value is 13.}

\item{second_par}{The choice for second hyper parameter, i.e., ROR, PRR, RR, Confidence.}

\item{minsup}{The threshold for the first hyper parameter - Support. Default value is 3.}

\item{minthresh}{The threshold for the second hyper parameter. Default value is 1.}
}
\value{
Data frame containing frequently co-occurring drugs/vaccines and AEs.
}
\description{
This is an R package to identify frequently co-occurring drug/vaccine and adverse event (AE) pairs from safety datasets using either the traditional Apriori method ("Confidence" as second hyper parameter) or an improved Apriori method (using "ROR", "RR", or "PRR" as the second hyper parameter)
}
\examples{
improvedApriori(data=kuo_data, numdrug=5, numae=13, numreports=13, second_par="ROR", minsup=3, minthresh=1)

}
