\name{f_test_RSRL}
\alias{f_test_RSRL}
\title{
Function computing coefficients and significance levels of the RSRL and mRSRL
}
\description{
This function computes coefficients and significance levels of the RSRL and mRSRL. It performs the (under)diversification test of a given portfolio compared to its benchmark.
}
\usage{
f_test_RSRL(v_input_p_r, v_input_b_r, input_c, input_b, input_sim, b_input_s, input_prob)
}
\arguments{
  \item{v_input_p_r}{A vector of portfolio returns}
  \item{v_input_b_r}{A vector of portfolio returns}
  \item{input_c}{A numerical value (number of data repetitions)}
  \item{input_b}{A numerical value (size of the block - time dimension)}
  \item{input_sim}{A numerical value (number of simulations)}
  \item{b_input_s}{A boolean value (percentile or studentized bootstrap)}
  \item{input_prob}{A numerical value (probability)}
}
\value{

\item{RSRL}{A numerical value (RSRL coefficient)}
\item{Signif_level_RSRL}{Numerical value (RSRL significance level)}
\item{mRSRL}{A numerical value (RSRL coefficient)}
\item{Signif_level_mRSRL}{Numerical value (mRSRL significance level)}

}
\references{
Candelon, Bertrand, Franz Fuerst, and Jean-Baptiste Hasse. "Diversification Potential in Real Estate Portfolios." (2020) Cambridge Working Paper.
}
\author{
Jean-Baptiste Hasse
}
\examples{
# NOT RUN {

  # Load data
  data("data_efficient_portfolios_returns")

  # Prepare data
  v_port <- data_efficient_portfolios_returns[,2]
  v_bench <- data_efficient_portfolios_returns[,1]

  # Test RSRL and mRSRL
  f_test_RSRL(v_port, v_bench, 10, 2, 1000, TRUE, 0.95)

# }
}


\keyword{Diversification}
\keyword{Finance}
\keyword{Econometrics}
