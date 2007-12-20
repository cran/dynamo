\name{plot.dm}
\alias{plot.dm}
\title{Plot Diagnostics for a dm Object}
\description{
  Plots of the data and fitted dynamic model.
}
\usage{
plot.dm( x , ... )
}
\arguments{
	\item{x}{a dynamic model object.}
	\item{...}{other parameters to be passed through to plotting functions.}
}
\details{}
\author{Christian T. Brownlees}

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
\keyword{ models }
