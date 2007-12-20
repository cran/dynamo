\name{rdm}
\alias{rdm}
\title{Dynamic Model Simulation}
\description{
	'rdm' is used to simulate dynamic models.
}
\usage{
rdm(formula,innovations=NULL,n,param,data=parent.frame(),seed=as.integer(Sys.time()),lst=FALSE)
}
\arguments{
	\item{formula}{a dynamo symbolic description of the model to be simulated.}
	\item{innovations}{the family of distributions of model innovations.}
	\item{n}{number of observations to be simulated}
	\item{param}{dynamic model parameters}
	\item{data}{data frame with exogenous data}
	\item{seed}{simulation seed}
	\item{lst}{data frame option}
}
\details{}
\value{
	If the 'lst' argument (which stands for 'data.frame') is FALSE (the default)
	the function returns the simulated series.
	If 'lst' is TRUE the function returns a data frame containing the simulated series and
	the simulated conditional mean and variance of the series.
}
\references{ ~put references to the literature/web site here ~ }
\author{Christian T. Brownlees}
\note{}
\seealso{\code{\link{dm}}}
\examples{
# ACD(1,1)
# simulate data
y <- rdm( ~acd(1,1), n=1000, param=c(0.2,0.1,0.8))

# plot
plot( y )
}
\keyword{ models }
