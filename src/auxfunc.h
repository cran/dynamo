/*
    Functions utilized by rcppfunc.cpp. 

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
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#include <RcppArmadillo.h>
#include <math.h>
using namespace std;
using namespace arma;

////////////////////////////////// Auxiliary functions
//////////////////// Direct RH-transform of a flat 3d array (matrix) M by a matrix X
mat RHmat(mat const& X, mat const& M,int col, int sli){

int rowx = X.n_rows;

////matrix multiply
mat XM = X * M;

////make matrix into rotated (!) cube (on matrix form)
mat Mnew(col, sli * rowx);
for (int s = 0; s < sli; s++) {

for (int c = 0; c < col; c++) {

for (int r = 0; r < rowx; r++) {

Mnew(c, s + r * sli) = XM(r, c + s * col);

}

}

}

return Mnew;

}

