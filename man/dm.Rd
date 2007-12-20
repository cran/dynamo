\name{dm}
\alias{dm}
\alias{print.dm}
\alias{dm.spec}
\title{ Fitting Dynamic Models }
\description{
  'dm' is used to fit dynamic models. It can be used to estimate 
  several univariate dynamic models including ARMA, GARCH, ARMA-GARCH, ACD and MEM.
}
\usage{
dm(formula,innovations=NULL,data=parent.frame(),est='mle',maxiter=NULL,param=NULL,ftol=NULL,log='none')
}
\arguments{
	\item{formula}{a dynamo symbolic description of the model to be simulated.}
	\item{innovations}{the family of distributions of model innovations.}
	\item{data}{data frame with exogenous data.}
	\item{est}{estimation method.}
	\item{maxiter}{max iter.}
	\item{param}{dynamic model parameters.}
	\item{ftol}{tollerance.}
	\item{log}{log level.}
}
\details{
}
\value{
  'dm' returns an object of 'class' '"dm"'.

  An object of class '"dm"' is a list containg:

  \item{coefficients}{fitted parameters.}
}
\references{}
\author{Christian T. Brownlees}
\note{}
\seealso{ \code{\link{summary.dm}} \code{\link{print.dm}} \code{\link{rdm}} }
\examples{
# AR(1) GARCH(1,1) 
# simulate data
#y <- rdm( y~0+ar(1)+garch(1,1) , innovations='st' , 1000 , param=c(0.2,0.2,0.1,0.8,0.1) )
# estimate
#ar1garch11 <- dm( y~0+ar(1)+garch(1,1) , innovations='st' )
# view estimation results
#summary( ar1garch11 )

# ACD(1,1)
# simulate data
#y <- rdm( y~acd(1,1), 1000 , param=c(0.1,0.8,0.1,20.0) )
# estimate
#acd11 <- dm( y~acd(1,1)  )
# view estimation results
#summary( acd11 )
1+1
}
\keyword{ models }
