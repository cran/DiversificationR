\name{f_SR}
\alias{f_SR}
\title{
Function computing the Sharpe ratio or one of its modified version
}
\description{
This function computes the Sharpe ratio (SR) or one of its modified version (mSR) from two vectors of financial returns (a given portfolios and its benchmark).}
\usage{
f_SR(v_input_data_portfolio, v_input_data_benchmark, c_input_method, input_prob)
}
\arguments{
  \item{v_input_data_portfolio}{A vector of numerical values (returns)}
  \item{v_input_data_benchmark}{A vector of numerical values (returns)}
  \item{c_input_method}{A vector of characters (method)}
  \item{input_prob}{A numerical value (probability)}
}
\value{
\item{result}{A numeric value}
}
\references{
Bali, Turan G., Stephen J. Brown, and K. Ozgur Demirtas. "Do hedge funds outperform stocks and bonds?." Management Science 59.8 (2013): 1887-1903.

Favre, Laurent, and José-Antonio Galeano. "Mean-modified value-at-risk optimization with hedge funds." The journal of alternative investments 5.2 (2002): 21-25.

Gregoriou, Greg N., and Jean-Pierre Gueyie. "Risk-adjusted performance of funds of hedge funds using a modified Sharpe ratio." The Journal of wealth management 6.3 (2003): 77-83.

Sharpe, William F. "The sharpe ratio." Journal of Portfolio Management 21.1 (1994): 49-58.

Sharpe, William F. "Mutual fund performance." The Journal of business 39.1 (1966): 119-138.
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
  v_rf <- v_bench


  # Compute the Reward-to-Variablity Ratio as in Sharpe (1966)
  f_SR(v_port, v_rf, "", 0.95)

  # Compute the Sharpe ratio as in Sharpe (1994)
  f_SR(v_port, v_bench, "S", 0.95)

  # Compute the modified Sharpe ratio as in Favre and Galeano (2002) and Gregoriou and Gueyie (2003)
  f_SR(v_port, v_bench, "FG-GG", 0.95)

  # Compute the modified Sharpe ratio as in Bali et al. (2013)
  f_SR(v_port, v_bench, "BBD", 0.95)

# }
}


\keyword{Risk}
\keyword{Finance}
\keyword{Econometrics}
\keyword{Shape ratio}
