\name{f_circular_bloc_bootstrap}
\alias{f_circular_bloc_bootstrap}
\title{
Function computing a circular block bootstrap
}
\description{
This function computes a circular block-resampling bootstrap of a matrix of returns.
}
\usage{
f_circular_bloc_bootstrap(m_input_data_series, input_c, input_b, input_prob)
}
\arguments{
  \item{m_input_data_series}{A matrix of assets or portfolios returns (one per column)}
  \item{input_c}{A numerical value (number of wrapping the data around in a circle)}
  \item{input_b}{A numerical value (length of block size - time dimension)}
  \item{input_prob}{A numerical value (probability)}
}
\value{
\item{RSRL}{A numerical value (bootstrapped RSRL)}
\item{mRSRL}{A numerical value (bootstrapped mRSRL)}
\item{bootstapped_series}{A matrix of numerical values (bootstrapped returns)}
}
\references{
Efron, B. "Bootstrap methods: another look at the jackknife." The Annals of Statistics 7 (1979): 1-26.

Hall, Peter, Joel L. Horowitz, and Bing-Yi Jing. "On blocking rules for the bootstrap with dependent data." Biometrika 82.3 (1995): 561-574.

Politis, Dimitris N., and Joseph P. Romano. "A circular block-resampling procedure for stationary data." Exploring the limits of bootstrap 2635270 (1992).
}
\author{
Jean-Baptiste Hasse
}
\examples{
# NOT RUN {

  # Load data
  data("data_efficient_portfolios_returns")
  m_example_returns <- data_efficient_portfolios_returns[,1:2]

  # Compute Circular bootstap
  f_circular_bloc_bootstrap(m_example_returns, 10, 2, 0.95)

# }
}


\keyword{Econometrics}
\keyword{Bootstap}
