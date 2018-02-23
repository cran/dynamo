/*
    Intended for use with R.
    Copyright (C) 2018 Adam Lund

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>
*/

//// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include "auxfunc.h"
using namespace std;
using namespace arma;

//////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////// construct design ///////////////////////////////////
//[[Rcpp::export]]
Rcpp::List getdes(arma::mat Ymat, //the repsonse V_0,,,V_M
                  arma::mat Vm1mat, // V_-1,,,V_M-1
                  arma::mat Vmat, //all  !!! data V_-L...V_M
                  arma::mat phix,arma::mat phiy,arma::mat Phil,
                  int t0//, //=1 include lag or in filter and =0 dont
                  ){
  //arma::mat Vmatrot, //all (!) rotated !!! data V_-L...V_M

//it seems that as Phixyl is computed using rotated data Vmat we
//actually compute Philyx and use that in the estimation...

/////////////////////////// Dimensions of entire data //////////////////////////////////
//L = tau is always = #lags in the model and Phit.n_rows = L when t0 = 0
//and = L + 1 when t0 = 1  so in other words Phit.n_rows = Lfilt = L + t0!.

int L = Phil.n_rows - t0,
    //Nt = Vmatrot.n_rows,
    Nx = phix.n_rows, Ny = phiy.n_rows, Nt = Vmat.n_elem / (Nx * Ny),
    px = phix.n_cols, py = phiy.n_cols, pl = Phil.n_cols;

int M = Nt - L - 1;

arma::mat HtY(px, py), Htmp, //Phixylrot(M, px * py * pl),
          Phixyl(M, px * py * pl),
//Vlagmatrot = Vmatrot(span(0, L - 1), span::all),
Vlagmat = Vmat.cols(0, Ny * L - 1);

/////////////////////////// Get Phixyl ///////////////////////////////////////////
for(int k = 0; k < M; k++){ //time k submatrix

//if(k > 0){Vlagmatrot = join_cols(Vlagmatrot(span(1, L - 1), span::all), Vmatrot(L - 1 + k, span::all));}
if(k > 0){Vlagmat = join_rows(Vlagmat.cols(Ny, Ny * L - 1), Vmat.cols(Ny * (L + k - 1), Ny * (L + k) - 1));}

//Phixylrot.row(k) = vectorise(RHmat(phix.t(), RHmat(phiy.t(), RHmat(Phil.t(), Vlagmatrot, Ny, Nx), Nx, pl), pl, py)).t(); //correct order?????
Phixyl.row(k) = vectorise(RHmat(Phil.t(), RHmat(phiy.t(), RHmat(phix.t(), Vlagmat, Ny, L), L, px), px, py)).t(); //correct order?????

}

////////////////////////// compute Hty
HtY.fill(0);
for(int i = 0; i < M; i++){HtY = HtY + phix.t() * (Vm1mat(span::all, span(i * Ny, (i + 1) * Ny - 1)) % Ymat(span::all, span(i * Ny, (i + 1) * Ny - 1))) * phiy;}

Rcpp::List output = Rcpp::List::create(Rcpp::Named("Phixyl") = Phixyl,
                                    //   Rcpp::Named("Phixylrot") = Phixylrot,
                                       Rcpp::Named("HtY") = HtY);

return output;

}

//////////////////////////////////////////////////////////////////////////////////////////
// ////////////////////////////// simulate from fitted dynamical model  //////////////////////////////
// //[[Rcpp::export]]
// Rcpp::List SimV(Rcpp::NumericVector shat,
//                 Rcpp::NumericVector fhat,
//                 arma::mat hhat, //should be zero when t0=1
//                 Rcpp::NumericVector v,
//                 arma::mat E,
//                 int t0, //include t0 (=0 or 1) in the filter
//                 int M, //# sim steps
//                 int L, //# lag steps
//                 double Deltat){
// 
//   Rcpp::NumericVector vecs(shat);
//   Rcpp::IntegerVector SDim = vecs.attr("dim");
//   arma::cube Shat(vecs.begin(), SDim[0], SDim[1], SDim[2], false);
// 
//   Rcpp::NumericVector vecf(fhat);
//   Rcpp::IntegerVector FDim = vecf.attr("dim");
//   arma::cube Fhat(vecf.begin(), FDim[0], FDim[1], FDim[2], false);
// 
//   Rcpp::NumericVector vecv(v);
//   Rcpp::IntegerVector VDim = vecv.attr("dim");
//   arma::cube V(vecv.begin(), VDim[0], VDim[1], VDim[2], false);
// 
//   int mn_index, N1 = Fhat.n_rows, N2 = Fhat.n_cols;
// 
//   double filter = 0;
// 
//   ////////////////////////// Simulate the M time steps ////////////////////////
//   for(int k = 0; k < M; k++){//text time index 0,...,M-1
// 
//     for(int n = 0; n < N2; n++){//c++ index
//       for(int m = 0; m < N1; m++){//c++ index
//         mn_index = m + N1 * n;
// 
//         for(int j = 0; j < N2; j++){//c++ index
//           for(int i = 0; i < N1; i++){//c++ index
//             for(int l = -L; l < t0; l++){//text time index -L,...,-1, t0 = 0 means no lag1 in filter
// 
//               filter += V(m, n, k + l + L) * Fhat(m, n, (l + L) + i * L + j * L * N1);
// 
//             }
//           }
//         }
// 
// V(m, n, k + 1 + L) =   Shat(m, n, k) * Deltat + pow(Deltat, 2) * filter + hhat(m, n) * V(m, n, k + L) + E(k, mn_index);
// filter = 0;
// 
//       }
//     }
//   }
// 
//   Rcpp::List output = Rcpp::List::create(Rcpp::Named("V") = V);
// 
//   return output;
// 
// }
// 
// //////////////////////////////////////////////////////////////////////////////////////////
// ////////////////////////////// get fitted drift //////////////////////////////
// //[[Rcpp::export]]
// Rcpp::List getMu(Rcpp::NumericVector shat,
//                  Rcpp::NumericVector fhat,
//                  arma::mat hhat, //should be zeroarma::mat when t0=1
//                  Rcpp::NumericVector v,
//                  int L, //# lag steps not # lag points!!!
//                  int t0,
//                  double Deltat){
// 
//   Rcpp::NumericVector vecs(shat);
//   Rcpp::IntegerVector sDim = vecs.attr("dim");
//   arma::cube Shat(vecs.begin(), sDim[0], sDim[1], sDim[2], false);
// 
//   Rcpp::NumericVector vecf(fhat);
//   Rcpp::IntegerVector fDim = vecf.attr("dim");
//   arma::cube Fhat(vecf.begin(), fDim[0], fDim[1], fDim[2], false);
// 
//   Rcpp::NumericVector vecv(v);
//   Rcpp::IntegerVector VDim = vecv.attr("dim");
//   arma::cube V(vecv.begin(), VDim[0], VDim[1], VDim[2], false);
// 
//   int N = V.n_slices, N1 = V.n_rows, N2 = V.n_cols, M = N - L - 1;
//   arma::cube F(N1, N2, M), Mu(N1, N2, M);
// 
//   double sumlagmn = 0;
//   for(int k = 0; k < M; k++){//text time index 0,...,M-1
// 
//     for(int n = 0; n < N2; n++){//c++ index
//       for(int m = 0; m < N1; m++){//c++ index
// 
//         for(int j = 0; j < N2; j++){//c++ index
//           for(int i = 0; i < N1; i++){//c++ index
//             for(int l = -L; l < t0; l++){//text time index -L,...,t0 (t0=0 means no lag 1 in filter)
// 
//               sumlagmn += V(i, j, k + l + L) * Fhat(m, n, l + L + i * L + j * L * N1);
// 
//             }
//           }
//         }
// 
//         F(m, n, k) = Deltat * sumlagmn;
//         sumlagmn = 0;
// 
//       }
//     }
// 
//   }
// 
//   //get the drift...
//   for(int k = 0; k < M; k++){//text time index 0,...,M-1
// 
//     for(int m = 0; m < N1; m++){//c++ index
//       for(int n = 0; n < N2; n++){//c++ index
// 
//         Mu(m, n, k) = (Shat(m, n, k) + F(m, n, k)) * Deltat + hhat(m, n) * V(m, n, k + L);
// 
//         //h is not * by deltat as h = (1 + h_old*delta) when lag1 not in filter ie t0=0!!!
//       }
//     }
// 
//   }
// 
//   Rcpp::List output = Rcpp::List::create(Rcpp::Named("F") = F,
//                                          Rcpp::Named("Mu") = Mu);
// 
//   return output;
// 
// }
//
// //////////////////////////////////////////////////////////////////////////////////////////
// /////////////////////////////proccees raw data to obtaian a film V//////////////////////////////////////////
// //[[Rcpp::export]]
// Rcpp::List getV(arma::mat const DeltaF,
//                 arma::mat Diode,
//                 double noise,
//                 double Deltat,
//                 double Deltax,
//                 double Deltay,
//                 double tstart,
//                 double tend,
//                 double xostart, //spatial original
//                 double xoend,
//                 double yostart,
//                 double yoend,
//                 double xcstart, //spatial cropping
//                 double xcend,
//                 double ycstart,
//                 double ycend){
// 
//   /////////////////////////// Dimensions of entire data ///////////////////////////////////
//   const size_t N = (tend - tstart) / Deltat + 1;  //#time points (indices) as in text
//   int d1o = (xoend - xostart) / Deltax + 1;       //#x points (indices) as in text
//   int d2o = (yoend - yostart) / Deltay + 1;       //#y points (indices) as in text
// 
//   /////////////////////////// Restructure data as 3d array ////////////////////////////////
//   arma::cube Vo(d1o, d2o, N);
//   for(int nt = 0; nt < N; nt++){
//     for(int nx = 0; nx < d1o; nx++){
//       for(int ny = 0; ny < d2o; ny++){
// 
//         if((Diode(nx, ny) - 1) != -1){
// 
//           Vo(nx, ny, nt) = DeltaF((Diode(nx, ny) - 1), nt);
// 
//         }else{
// 
//           //Vo(nx, ny, nt) = sqrt(noise)/100 * Rcpp::as<double>(Rcpp::rnorm(1, 0, 1));
//           Vo(nx, ny, nt) = 0;
// 
//         }
// 
//       }
//     }
//   }
// 
//   /////////////////////////// Coordinates and dimensions of cropped data  /////////////////
//   double xstart = xcstart;                    //top cut
//   double xend = xcend;                        //bottom cut
//   double ystart = ycstart;                    //left cut
//   double yend = ycend;                        //right cut
//   int d1 = (xend - xstart) / Deltax + 1;      //#x points (indices) to model
//   int d2 = (yend - ystart) / Deltay + 1;      //#y points (indices) to model
// 
//   //crop INDICES, everyarma::vector starts at 0 and ends at length(vector) - 1
//   int xistart =  xstart / Deltax - 1;                                     //top index
//   int xiend = (xend - xostart) / Deltax + 1 - 1;                          //bottom index
//   int yistart =  ystart / Deltay - 1;                                     //left index
//   int yiend = (yend - yostart) / Deltay + 1 - 1;                          //right index
// 
//   /////////////////////////// Crop spatial dims in 3d data array /////////////////////////
//   arma::cube V(d1, d2, N);
//   V.fill(0.42);
//   for(int nt = 0; nt < N; nt++){
//     for(int nx = xistart; nx <= xiend; nx++){
//       for(int ny = yistart; ny <= yiend; ny++){
// 
//         V(nx - xistart, ny - yistart, nt) = Vo(nx, ny, nt);
// 
//       }
//     }
//   }
// 
//   Rcpp::List output = Rcpp::List::create(Rcpp::Named("V") = V);
//   return output;
// }
//
// //////////////////////////////////////////////////////////////////////////////////////////
// ////////////////////////MSE///////////////////////////////////
// //[[Rcpp::export]]
// Rcpp::List MSE(Rcpp::NumericVector Beta, Rcpp::NumericVector Xs, Rcpp::NumericVector Ys){
// 
//   Rcpp::NumericVector vecB(Beta);
//   Rcpp::IntegerVector BDim = vecB.attr("dim");
//   arma::cube B(vecB.begin(), BDim[0], BDim[1], BDim[2], false);
//   Rcpp::NumericVector vecX(Xs);
//   Rcpp::IntegerVector XDim = vecX.attr("dim");
//   arma::cube X(vecX.begin(), XDim[0], XDim[1], XDim[2], false);
//   Rcpp::NumericVector vecY(Ys);
//   Rcpp::IntegerVector YDim = vecY.attr("dim");
//   arma::cube Y(vecY.begin(), YDim[0], YDim[1], YDim[2], false);
// 
//   /////////////////////////// MSE ////////////////////////////////
// 
//   int NoM = B.n_cols;      //#models
//   int P = Y.n_cols;        //#pixels
//   int NoP = P * X.n_cols;  //#par
//   int NoG = X.n_slices;    //#groups
// 
//   arma::cube Yfit(NoG, NoM, Y.n_rows * Y.n_cols);
//   arma::mat MSE(NoG, NoM);
//   arma::mat yfit;
//   arma::vec ydiff;
// 
//   arma::mat Bhat(X.n_cols, P);
//   arma::vec tmp(Bhat.memptr(), Bhat.n_elem, false, false);
//   for(int m = 0; m < NoM; m++){
//     for(int l = 0; l < NoG; l++){
// 
//       tmp = B.subcube(arma::span(l), arma::span(m), arma::span::all);
//       yfit =  X.slice(l) * Bhat;
//       ydiff =arma::vectorise(Y.slice(l) - yfit);
// 
//       Yfit(arma::span(l), arma::span(m), arma::span::all) =arma::vectorise(yfit);
//       MSE(l, m) = accu(pow(ydiff, 2)) / (Y.n_rows * Y.n_cols);
// 
//     }
//   }
// 
//   Rcpp::List output = Rcpp::List::create(Rcpp::Named("MSE") = MSE, Rcpp::Named("Yfit") = Yfit);
//   return output;
// 
// }
// 

