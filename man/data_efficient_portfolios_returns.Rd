\name{data_efficient_portfolios_returns}
\alias{data_efficient_portfolios_returns}
\docType{data}
\title{
Efficient portfolios returns
}
\description{
This dataset includes efficient real estate portfolios returns from 1999 to 2018 (annual frequency). Overall, country- and -sector level portfolios are computed in both Markowitz and Black-Litterman frameworks.
}
\usage{data("data_efficient_portfolios_returns")}
\format{
  The format is:
 num [1:19, 1:6] 7.87 6.93 6.32 6.92 7.1 ...
 - attr(*, "dimnames")=List of 2
  ..$ : NULL
  ..$ : chr [1:6] "v_Overall_M" "v_Countries_M" "v_Sectors_M" "v_Overall_BL" ...
}

\source{
Author's own calculations based on MSCI dataset.
}
\references{
Candelon, Bertrand, Franz Fuerst, and Jean-Baptiste Hasse. "Diversification Potential in Real Estate Portfolios." (2020) Cambridge Working Paper.
}
\examples{
data(data_efficient_portfolios_returns)
head(data_efficient_portfolios_returns)
}
\keyword{datasets}
