\name{f_VaR}
\alias{f_VaR}
\title{
Function computing Value-at-Risk and modified Value-at-Risk
}
\description{
This function computes the Value-at-Risk (VaR) or the modified Value-at-Risk (mVaR) from a vector of financial returns. mVaR is also called the Cornish-Fisher expansion of Value-at-Risk. Compared to classic VaR, mVaR adequately accounts for the non-normality of returns.
}
\usage{
f_VaR(v_input_data, b_input_var_modified, input_prob)}
\arguments{
  \item{v_input_data}{A vector including an asset or portfolio returns}
  \item{b_input_var_modified}{A boolean to compute VaR or mVaR}
  \item{input_prob}{A numerical value (probability)}
}

\value{
\item{result}{A numeric value}
}
\references{
Cornish, Edmund A., and Ronald A. Fisher. "Moments and cumulants in the specification of distributions." Revue de l'Institut international de Statistique (1938): 307-320.

Jorion, Philippe. "Risk2: Measuring the risk in value at risk." Financial analysts journal 52.6 (1996): 47-56.
}
\author{
Jean-Baptiste Hasse
}
\examples{
# NOT RUN {

  # Load data
  data("data_efficient_portfolios_returns")

  # Prepare variables
  v_port <- data_efficient_portfolios_returns[,1]

  # Compute VaR
  f_VaR(v_port, FALSE, 0.95)

  # Compute modified VaR
  f_VaR(v_port, TRUE, 0.95)

# }
}


\keyword{Risk}
\keyword{Finance}
\keyword{Econometrics}
\keyword{VaR}

