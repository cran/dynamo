\name{predict.dm}
\alias{predict.dm}
\title{Predict Method for Dynamic Model Fits}
\description{
    Predicted values based on dynamic model object.
}
\usage{
predict.dm( object, type=stop('Forecast type not specified!',call.=FALSE), iterate=NULL, nfor=1, ql=0.05, qu=0.95, yo=NULL, Xo=matrix(0), Zo=matrix(0), ... )
}
\arguments{
	\item{object}{a dynamic model object.}
	\item{type}{number of digits.}
	\item{iterate}{some desc}
	\item{nfor}{some desc}
	\item{ql}{some desc}
	\item{qu}{some desc}
	\item{yo}{some desc}
	\item{Xo}{some desc}
	\item{Zo}{some desc}
	\item{...}{some desc}
}
\details{}
\references{}
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
