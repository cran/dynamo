\name{dynamo-package}
\alias{dynamo-package}
\alias{dynamo}
\docType{package}
\title{Dynamic Models}
\description{
DynaMo (Dynamic Models) is a package for time series analysis with a special focus on the 
time series models used in Financial Econometrics. The package provides methods for 
simulation, estimation, inference, regularization and prediction of a set of 
univariate dynamic models including: ARMA, ARMA-GARCH, ACD, MEM, Spline-MEM.
}
\details{
\tabular{ll}{
Package: \tab dynamo\cr
Type: \tab Package\cr
Version: \tab 0.1.3\cr
Date: \tab 2007-12-20\cr
License: \tab GPL 2\cr
}
rdm 
dm
}
\author{
Christian T. Brownlees

Maintainer: Christian T. Brownlees <ctb@ds.unifi.it>
}
\references{
}
\keyword{ package }
\seealso{}
\examples{
#y ~ rdm( ~arma(1,1)+garch(1,1) , 1000 )
#a11g11 <- dm( y~arma(1,1)+garch(1,1) )
1+1
}
