\name{f_RSRL}
\alias{f_RSRL}
\title{
Function computing the RSRL or the RSRL
}
\description{
This function computes the relative Sharpe ratio loss (RSRL) or its modified version (mRSRL) from two vectors of financial returns (a given portfolio and its benchmark). RSRL and mRSRL are both (under)diversification measures. Compared to RSRL, mRSRL is robust to the non-normality of returns.
}
\usage{
f_RSRL(v_input_data_portfolio, v_input_data_benchmark, b_input_RSRL_modified, input_prob)
}
\arguments{
  \item{v_input_data_portfolio}{A vector of returns}
  \item{v_input_data_benchmark}{A vector of returns}
  \item{b_input_RSRL_modified}{A boolean value}
  \item{input_prob}{A numerical value}
}

\value{
\item{result}{A numeric value}
}
\references{
Calvet, Laurent E., John Y. Campbell, and Paolo Sodini. "Down or out: Assessing the welfare costs of household investment mistakes." Journal of Political Economy 115.5 (2007): 707-747.

Candelon, Bertrand, Franz Fuerst, and Jean-Baptiste Hasse. "Diversification Potential in Real Estate Portfolios." (2020) Cambridge Working Paper.
}
\author{
Jean-Baptiste Hasse
}
\examples{
# NOT RUN {

  # Load data
  data("data_efficient_portfolios_returns")

  # Prepare variables
  v_port <- data_efficient_portfolios_returns[,2]
  v_bench <- data_efficient_portfolios_returns[,1]

  # Compute RSRL
  f_RSRL(v_port, v_bench, FALSE, 0.95)

  # Compute mRSRL
  f_RSRL(v_port, v_bench, TRUE, 0.95)

# }
}


\keyword{Risk}
\keyword{Econometrics}
