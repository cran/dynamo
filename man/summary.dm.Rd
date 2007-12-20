\name{summary.dm}
\alias{summary.dm}
\alias{print.summary.dm}
\title{Summarizing Dynamic Model Fits}
\description{
   summary method for class "dm"
}
\usage{
summary.dm( object , diag.resid.lag=5 , ... )
}
\arguments{
	\item{object}{a dynamo symbolic description of the model to be simulated.}
	\item{diag.resid.lag}{the family of distributions of model innovations.}
	\item{...}{some desc}
}
\details{}
\value{
	The function summary.dm computes and returns a list of summary statistics of the
	fitted dynamic model given in obejct.
}
\references{}
\author{Christian T. Brownlees}
\note{}
\seealso{\code{\link{rdm}}}
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
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ models }
% \keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
