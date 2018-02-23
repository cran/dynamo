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
#' @name dynamo_internal
#' @aliases design
#' @title Function that constructs design components.
#' @description  Internal design function.
#' @details This function is not intended for use by users.
#' @keywords internal
#' @usage ...
#' @author Adam Lund
#' 
#' ################# construct design components
design <- function(V,
                   phix,  phiy, Phil,
                   Nt, Lp1,
                   t0 = 0) {

Vmat <- matrix(V, dim(V)[1], dim(V)[2] * dim(V)[3])
Ymat <- matrix(V[, , (Lp1 + 1):Nt], dim(V[, , (Lp1 + 1):Nt])[1], dim(V[, , (Lp1 + 1):Nt])[2] * dim(V[, , (Lp1 + 1):Nt])[3])
Vm1mat <- matrix(V[, , Lp1:(Nt - 1)], dim(V[, , Lp1:(Nt - 1)])[1], dim(V[, , Lp1:(Nt - 1)])[2] * dim(V[, , Lp1:(Nt - 1)])[3])

out <- getdes(Ymat,
              Vm1mat,
              Vmat,
              phix,  phiy,  Phil,
              t0)

out

}

