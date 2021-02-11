\name{f_diversification_measurement}
\alias{f_diversification_measurement}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Function computing portfolio diversification measures
}
\description{
This function computes several portfolio diversification measures: Portfolio Diversification Index (PDI), Diversification Ratio (DR), Diversification Delta (DD) and Diversification Delta Star (DD*).
}
\usage{
f_diversification_measurement(v_input_weights, m_input_returns, c_input_method)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{v_input_weights}{A vector of numerical values (asset weights)}
  \item{m_input_returns}{A matrix of numerical values (asset returns)}
  \item{c_input_method}{A character value (name of the diversification measure)}
}
\value{
\item{result}{A numeric value}
}
\references{
Rudin, Alexander M. "A portfolio diversification index." The Journal of Portfolio Management 32.2 (2006): 81-89.

Choueifaty, Yves, and Yves Coignard. "Toward maximum diversification." The Journal of Portfolio Management 35.1 (2008): 40-51.

Vermorken, Maximilian A., Francesca R. Medda, and Thomas Schroder. "The diversification delta: A higher-moment measure for portfolio diversification." The Journal of Portfolio Management 39.1 (2012): 67-74.

Flores, Yuri Salazar, et al. "The diversification delta: A different perspective." The Journal of Portfolio Management 43.4 (2017): 112-124.
}
\author{
Jean-Baptiste Hasse
}
\examples{
# NOT RUN {

  # Load data
  data("data_efficient_portfolios_returns")
  m_assets_returns <- data_efficient_portfolios_returns
  number_assets <- length(m_assets_returns[1,])
  v_weights <- rep(1/number_assets, number_assets)

  # Portfolio Diversification Index (PDI)
  f_diversification_measurement(v_weights, m_assets_returns, "Portfolio_Diversification_Index")

  # Diversification Ratio (DR)
  f_diversification_measurement(v_weights, m_assets_returns, "Diversification_Ratio")

  # Diversification Delta (DD)
  f_diversification_measurement(v_weights, m_assets_returns, "Diversification_Delta")

  # Diversification Delta Star (DD*)
  f_diversification_measurement(v_weights, m_assets_returns, "Diversification_Delta_Star")

# }
}

\keyword{Diversification}
\keyword{Finance}
\keyword{Econometrics}
