#     Description of this R script:
#     R function that constructs the design and runs the br gd pg algorithm using the glamlasso procedure.
#     Copyright (C) 2018 Adam Lund
#
#     This program is free software: you can redistribute it and/or modify
#     it under the terms of the GNU General Public License as published by
#     the Free Software Foundation, either version 3 of the License, or
#     (at your option) any later version.
#
#     This program is distributed in the hope that it will be useful,
#     but WITHOUT ANY WARRANTY; without even the implied warranty of
#     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#     GNU General Public License for more details.
#
#     You should have received a copy of the GNU General Public License
#     along with this program.  If not, see <http://www.gnu.org/licenses/>
#
#
#' @name fitmodel
#'
#' @title Sparse Estimation of Dynamical Array Models
#'
#' @description  An implementation of the method proposed in \cite{Lund and Hansen, 2018} for fitting 3-dimensional dynamical array models. 
#' The implementation is based on the glamlasso package, see  \cite{Lund et. al, 2017}, for efficient design matrix free lasso 
#' regularized estimation in a generalized linear array model. 
#' The implementation uses a block relaxation scheme to fit each individual component in the model using functions from the glamlasso package.
#'
#' @usage  fitmodel(V,
#'          phix, phiy, phil, phit,
#'          penaltyfactor,
#'          nlambda = 15,
#'          lambdaminratio = 0.0001,
#'          reltolinner = 10^-4,
#'          reltola = 10^-4,
#'          maxalt = 10)
#'
#' @param V is the \eqn{N_x\times Ny\times Nt\times G} array containing the data
#' @param phix is a \eqn{N_x\times px} matrix containing spatial basis functions
#' @param phiy is a \eqn{N_y\times py} matrix containing spatial basis functions
#' @param phil is a \eqn{Lp1\times pl} matrix containing temporal basis functions
#' @param phit is a \eqn{Nt\times pt} matrix containing temporal basis functions
#' @param penaltyfactor An array of size \eqn{p_1 \times \cdots \times p_d}. Is multiplied  with each element in \code{lambda} to allow differential shrinkage on the coefficients. 
#' @param nlambda the number of  \code{lambda} values to fit
#' @param lambdaminratio the smallest value for \code{lambda}, given as a fraction of 
#' @param reltolinner the convergence tolerance used with glamlasso 
#' @param reltola the convergence tolerance used for the alternation loop 
#' @param maxalt maximum nuber of alternations
#' 
#' @details This package contains an implementation of the method proposed in \cite{Lund and Hansen, 2018} 
#' for fitting a (partial) 3-dimensional 3-component dynamical array model to a \eqn{N_x\times N_y \times N_t \times G} data array \eqn{V}, where \eqn{N_x,N_y,N_t, G\in \{1,2,\ldots\}}. 
#' Note that \eqn{N_x, N_y, N_t} gives the number of observations in the two spatial dimensions and the temporal dimension respectively and \eqn{G} gives the number of trials. 
#' Let \eqn{L} be positive integer giving the length of the modelled delay, \eqn{M:=N_t - L - 1} and
#'  \eqn{\phi^{i}} be a \eqn{N_i\times p_i} matrix for \eqn{i\in \{x,y,l\}}. Then define a \eqn{M\times p_xp_yp_l} matrix 
#' \deqn{\Phi^{xyt}_{g} =  \left(  
#'  \begin{array}{ccc}    vec(\Phi^{xyl}_{1,g})
#'    \\
#'    \vdots
#'    \\
#'    vec(\Phi^{xyl}_{M,g})
#'    \end{array}
#'    \right),
#' }
#' for each \eqn{g\in \{1,\ldots,G\}}.  Here  using   the so called rotated \eqn{H}-transform \eqn{\rho} from \cite{Currie et al., 2006}, 
#'  \eqn{\Phi^{xyt}_{k, g}} is a  \eqn{p_x\times p_y\times p_l} array that, for each \eqn{k\in \{1,\ldots,M\}},  can  be computed  as 
#'\deqn{\Phi^{xyl}_{k, g} :=  \rho((\phi^l)^\top, \rho((\phi^y)^\top, \rho((\phi^x)^\top, V_{,,(k - L):(k - 1), g}))).}
#' Then we can write the model equation as for each trial \eqn{g} as
#' \deqn{V_{,,(L + 1):N_t, g}  = vec(\rho(\phi^{t}, \rho(\phi^{y}, \rho(\phi^{x}, A))) + \rho(\Phi^{xyt}, \rho(\phi^{y}, \rho(\phi^{x}, B))) + V_{,,-1:(M - 1), g} \odot C + E}
#' where \eqn{A} and \eqn{B} are  resp. \eqn{p_x\times p_y \times p_t} and \eqn{p_x\times p_y \times p_xp_yp_l}  coefficient arrays that are estimated and \eqn{C} 
#' is a \eqn{N_x\times N_y\times M} array containing \eqn{M} copies of the \eqn{N_x\times N_y} matrix \eqn{\rho(\phi^{y}, \rho(\phi^{x}, \Gamma))}, 
#' where \eqn{\Gamma} is a \eqn{p_x\times p_y} coefficient matrix  that is estimated.  
#' Finally \eqn{E} is a \eqn{N_x\times N_y\times M} array with Gaussian noise. See \cite{Lund and Hansen, 2018} for more details.
#' 
#' @return An object with S3 Class "dynamo". 
#' \item{Out}{A list where the first G items are the individual fits to the trials containing:}
#' \item{Out$V}{The \eqn{N_x\times N_y \times N_t} data array for trial \eqn{g}.}
#' \item{Out$Phixyl}{a \eqn{M\times p_xp_yp_l} matrix,  the convolution tensor.}  
#' \item{Out$BetaS}{\eqn{p_xp_yp_t\times}\code{nlambda} matrix containing the estimated parameters for the stimulus component for each \eqn{\lambda} value}
#' \item{Out$BetaF}{\eqn{p_xp_yp_xp_yp_l\times}\code{nlambda} matrix containing the estimated parameters for the filter component for each \eqn{\lambda} value}
#' \item{Out$BetaH}{\eqn{p_xp_y\times}\code{nlambda} matrix containing the estimated parameters for the instantaneous memory component for each \eqn{\lambda} value}
#' \item{Out$lambda}{vector of length \code{nlambda} containing the penalty parameters.}
#' \item{Out$Obj}{\code{maxalt}\eqn{\times}\code{nlambda} matrix containing the obective values for each alternation and each \eqn{\lambda}}
#' \item{phix}{is a \eqn{N_x\times p_x} matrix containing the basis functions over the spatial \eqn{x} domain.}
#' \item{phiy}{is a \eqn{N_y\times p_y} matrix containing the basis functions over the spatial \eqn{y} domain.}
#' \item{phil}{is a \eqn{L+1\times p_l} matrix containing the basis functions over the temporal domain for the delay.}
#' \item{phit}{is a \eqn{M\times p_t} matrix containing the basis functions over the temporal domain for the stimulus.}
#' 
#' @author  Adam Lund
#' 
#' Maintainer: Adam Lund, \email{adam.lund@@math.ku.dk}
#' 
#' @references 
#' Lund, A. and N. R. Hansen (2018). Sparse Network  Estimation for  Dynamical Spatio-temporal Array Models. 
#' \emph{In preparation},  url = {https://}.
#' 
#' Lund, A., M. Vincent, and N. R. Hansen (2017). Penalized estimation in 
#' large-scale generalized linear array models. 
#' \emph{Journal of Computational and Graphical Statistics}, 26, 3, 709-724.  url = {https://doi.org/10.1080/10618600.2017.1279548}.
#'
#' Currie, I. D., M. Durban, and P. H. C. Eilers (2006). Generalized linear
#' array models with applications to multidimensional smoothing. 
#' \emph{Journal of the Royal Statistical Society. Series B}. 68, 259-280. url = {http://dx.doi.org/10.1111/j.1467-9868.2006.00543.x}.
#' 
#' @keywords package 
#'
#' @examples 
#' \dontrun{
#' # Example showcasing the application from Lund and Hansen (2018).
#' ############################### data  
#' data(V)
#' ############################### constants 
#' Nx <- dim(V)[1]
#' Ny <- dim(V)[2]
#' Nt <- dim(V)[3]
#' L <- 50          #lag length in steps
#' Lp1 <- L + 1     #number of lag time points (= initial points)
#' t0 <- 0
#' M <- Nt - Lp1      #number of modelled time points
#' sl <- floor(200 / 0.6136) - 0 + 1   #stim start counted from -tau
#' sr <- sl + floor(250 / 0.6136)  #stim end counted from -tau
#' ##no. of basis func.
#' px <- 8
#' py <- 8
#' pl <- max(ceiling(Lp1 / 5), 4)
#' pt <- max(ceiling((Nt - sl) / 25), 4)
#' degx <- 2
#' degy <- 2
#' degl <- 3
#' degt <- 3
#' ############################### basis functions
#' library(MortalitySmooth)
#' phix <- round(MortSmooth_bbase(x = 1:Nx, xl = 1, xr = Nx, ndx = px - degx, deg = degx), 10)
#' phiy <- round(MortSmooth_bbase(x = 1:Ny, xl = 1, xr = Ny, ndx = py - degy, deg = degy), 10)
#' phil <- round(MortSmooth_bbase(x = -tau:(t0 - 1), xl = -tau, xr = (t0 - 1), 
#' ndx = pl - degl, deg = degl), 10)
#' phit <- round(MortSmooth_bbase(x = sl:Nt, xl = sl, xr = Nt, ndx = pt - degt, deg = degt), 10)
#' phit <- rbind(matrix(0,  (sl - 1) - Lp1, dim(phit)[2]), phit)
#' ############################### penalty weights  
#' wt <- c(1, 1, 2, 2, 3, 3, 3, 3, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 3, 3, 3, 3, 2, 1, 1)
#' penSlist <- list(matrix(1,  px,  py), matrix(1 / wt, dim(phit)[2], 1))
#' penF <- array(1, c(px, py, px * py * pl))
#' penH <- matrix(1, px, py)
#' penaltyfactor <- list(penSlist, penF, penH)
#' ############################### run algorithm
#' system.time({Fit <- fitmodel(V,
#'                              phix, phiy, phil, phit,
#'                              penaltyfactor, 
#'                              nlambda = 10,
#'                              lambdaminratio = 10^-1,
#'                              reltolinner = 10^-4,
#'                              reltola = 10^-4,
#'                              maxalt = 10)})
#'                               
#' ###############################  get one fit
#' modelno <- 6
#' fit <- Fit[[1]]
#' A <- array(fit$BetaS[, modelno], c(px, py, pt))
#' B <- array(fit$BetaF[, modelno], c(px, py, px * py * pl))
#' C <- array(fit$BetaH[, modelno], c(px, py))
#' shat <- RH(phit, RH(phiy, RH(phix, A)))
#' beta <- array(B, c(px, py, px, py, pl))
#' what <- RH(phil, RH(phiy, RH(phix, RH(phiy, RH(phix, beta)))))
#' ############################### compute summary network quantities
#' wbar <- apply(abs(what), c(1, 2, 3, 4), sum)
#' win <- apply(wbar, c(1, 2), sum)  
#' wout <- apply(wbar, c(3, 4), sum) 
#' indeg <- apply((what != 0), c(1, 2), sum)
#' outdeg <- apply((what != 0), c(3, 4), sum)
#' winnorm <- ifelse(indeg > 0, win / indeg, win)
#' woutnorm <- ifelse(outdeg > 0, wout / outdeg, wout)
#' ############################### plot summary network quantities
#' par(mfrow = c(2, 2), oma = c(0, 0 ,1, 0), mar = c(0, 0, 1, 0))
#' image(winnorm, main = paste("Time aggregated in effects"), axes = FALSE)
#' image(woutnorm, main = paste("Time aggregated out effects"), axes = FALSE)
#' timepoint <- which(shat[9, 9, ] == min(shat[9, 9, ]))
#' image(shat[,, timepoint], axes = FALSE, main = "Stimulus function")
#' plot(shat[1, 1, ], ylim = range(shat), type="l", main = "Stimulus function")
#' for(i in 1:Nx){for(j in 1:Ny){lines(shat[i, j, ])}}
#' abline(v = sl - Lp1, lty = 2)
#' abline(v = sr - Lp1, lty = 2)
#' }

library(glamlasso)
library(abind)
fitmodel <- function(V,
                     phix, phiy, phil, phit,
                     penaltyfactor,
                     nlambda = 15,
                     lambdaminratio = 0.0001,
                     reltolinner = 10^-4,
                     reltola = 10^-4,
                     maxalt = 10){

Out <- list()
class(Out) <- "dynamo"

if(length(dim(V)) == 3){V <- array(V, c(dim(V), 1))}
Nx <- dim(V)[1]
Ny <- dim(V)[2]
Nt <- dim(V)[3]
G <- dim(V)[4]
Lp1 <- dim(phil)[1] + 1
M <- Nt - Lp1
t0 = 0

px <- dim(phix)[2]
py <- dim(phiy)[2]
pl <- dim(phil)[2]
pt <- dim(phit)[2]
pxyl <- px * py * pl

penSlist<-penaltyfactor[[1]]

penS <- outer(penaltyfactor[[1]][[1]], c(penaltyfactor[[1]][[2]]))
penF <- penaltyfactor[[2]]
penH <- penaltyfactor[[3]]

for(g in 1:G){#fit G trials

############################### get design  ###################################
Ytr <- V[, , (Lp1 + 1):Nt, g]
n <- length(Ytr)
Vm1 <- V[ , , Lp1:(Nt - 1), g]

out <- design(V[, , 1:Nt, g], phix, phiy, phil, Nt, Lp1, t0 = 0)

Phixyl <- out$Phixyl
HtY <- out$HtY

BetasS <- matrix(NA, px * py * pt, nlambda)
BetasF <- matrix(NA, px * py * pxyl, nlambda)
BetasH <- matrix(NA, px * py, nlambda)
penaltyfactorall <- abind(penS, penaltyfactor[[2]], penaltyfactor[[3]])

############################### lambda ################################
absgradzeroall = abs(abind(RH(t(phit), RH(t(phiy), RH(t(phix), Ytr))), RH(t(Phixyl), RH(t(phiy), RH(t(phix), Ytr))), matrix(HtY, px, py)))  / n
absgradzeropencoef = absgradzeroall * (penaltyfactorall > 0)
penaltyfactorpencoef = (penaltyfactorall == 0) * 1 + penaltyfactorall
lambdamax = max(absgradzeropencoef / penaltyfactorpencoef)
lambda <- lambdamax * exp(seq( 0, log(lambdaminratio), length.out = nlambda))

############################### initialize####################
BetaS12 <- matrix(0.01, px, py)
BetaS3 <- matrix(0.01, pt, 1)
BetaS <- outer(BetaS12, c(BetaS3))
BetaF <- array(0, c(px, py, pxyl))
BetaH <- matrix(0, px, py)

YhatS <- RH(phit, RH(phiy, RH(phix,  BetaS)))
YhatF <- RH(Phixyl, RH(phiy, RH(phix, BetaF)))
YhatH <- Vm1 * array(RH(phiy, RH(phix, BetaH)), c(Nx, Ny, M))

Obj <- matrix(NA, maxalt + 1, nlambda)

############################### run alg####################

for(j in 1:nlambda){

obj <- rep(NA, maxalt + 1)
obj[1] <- sum((Ytr - YhatS - YhatF - YhatH)^2) / n + lambda[j] * (sum(abs(penS * BetaS)) + sum(abs(penF * BetaF)) + sum(abs(penH * BetaH)))

for(a in 1:maxalt){#alternations over components

for(b in 1:3){#alternate over components

if(b == 1){#stim

Y <- Ytr - YhatF - YhatH
fit <- glamlassoRR(list(phix, phiy, phit),
                   Y,
                   thetainit = list(BetaS12, BetaS3),
                   lambda = lambda[j],
                   penaltyfactor = penSlist,
                   reltolinner = reltolinner)

BetaS12 <- matrix(fit$coef12, px ,py)
BetaS3 <- matrix(fit$coef3, pt, 1)
BetaS <- array(outer(fit$coef12, c(fit$coef3)), c(px ,py ,pt))
YhatS <- RH(phit, RH(phiy, RH(phix, BetaS)))

}

if(b == 2){#filter

Y <- Ytr - YhatS - YhatH

fit <- glamlasso(list(phix, phiy, Phixyl),
                 Y,
                 thetainit = BetaF,
                 lambda = lambda[j],
                 penaltyfactor = penF,
                 reltolinner = reltolinner)

BetaFrot <- array(fit$coef, c(pxyl, py, px))
YhatF <- aperm(RH(phix, RH(phiy, RH(Phixyl, BetaFrot)))	, c(3, 2, 1))
BetaF <- aperm(BetaFrot, c(3, 2, 1))

BetaF <- array(fit$coef, c(px, py, pxyl))
YhatF <- RH(Phixyl, RH(phiy, RH(phix, BetaF)))

}

if(b == 3){#h

Y <- Ytr - YhatS - YhatF

fit <- glamlassoS(list(phix, phiy),
                  Y,
                  Vm1,
                  thetainit = BetaH,
                  lambda = lambda[j],
                  penaltyfactor = penH,
                  reltolinner = reltolinner)

BetaH <- matrix(fit$coef, px, py)
YhatH <- Vm1 * array(phix %*%  BetaH %*% t(phiy), c(Nx, Ny, M))

}

} #blocks

obj[a + 1] <- sum((Ytr - YhatS - YhatF - YhatH)^2) / n + lambda[j] * (sum(abs(penS * BetaS)) + sum(abs(penF * BetaF)) + sum(abs(penH * BetaH)))
if((obj[a + 1] - obj[a]) / obj[a] < reltola){break}

} #alt

Obj[, j] <- obj
BetasS[, j] <- c(BetaS)
BetasF[, j] <- c(BetaF)
BetasH[, j] <- c(BetaH)

} #lambda

fit <- list()
fit$V <- V
fit$Phixyl <- Phixyl
fit$BetaS <- BetasS
fit$BetaF <- BetasF
fit$BetaH <- BetasH
fit$lambda <- lambda
fit$Obj <- Obj

Out[[g]] <- fit
print(g)

}

Out$phix <- phix
Out$phiy <- phiy
Out$phil <- phil
Out$phit <- phit

Out

}

