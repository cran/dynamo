/*!
 *  dynamo
 *  Copyright (C) 2007, 2008 Christian T. Brownlees
 *
 *  This file is part of dynamo.
 * 
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 * 
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 * 
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef DYNAMO_BASE_H
#define DYNAMO_BASE_H

#include <math.h>
#include <string.h>

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_matrix_int.h>
#include <gsl/gsl_specfunc.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_statistics.h>
#include <gsl/gsl_multimin.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_sf.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_multifit.h>

//! \brief Lagged values of the v vector.
#define LAGV(v,i) (v->size-i-1)
//! \brief Lagged values of the M matrix.
#define LAGM(M,i) (M->size1-i-1)

/*! 
 * \brief Data check type enumeration.
 *
 * \sa dynamo_data_check
 */
enum DYNAMO_DATA_CHECK { 
	DYNAMO_DATA_CHECK_FULL=0 , /*!< Check all variables. */
	DYNAMO_DATA_CHECK_EXPL /*!< Check explanatory variables only. */
};

// typedefs
/*!
 * \brief Random deviates generator function type definition.
 * 
 * The typedef provides the function type definition for the functions used to simulated random deviates from the various random variables.
 *
 * \param y_t (unit size) vector of current value to be simulated
 * \param spec specification vector
 * \param param param vector
 * \param cond_t cond
 * \param r GSL random number generator
 * \param flag simulation status flag
 *
 * On exit the flag variable will be set to 1 if the simulated deviate is finite, 0 otherwise.
 *
 * \sa dynamo_sim
 */
typedef void (*dynamo_rand_t)(gsl_vector *y_t, const gsl_vector_int *spec, const gsl_vector *param, const gsl_vector *cond_t, const gsl_rng *r, int *flag);
typedef double (*dynamo_quant_t)(double p, const gsl_vector_int *spec, const gsl_vector *param, const gsl_vector *cond_t);

typedef double (*dynamo_cond_loglik_t)( const gsl_vector *y_t, const gsl_vector_int *spec, const gsl_vector *param, const gsl_vector *cond_t, int *flag);

/*!
 * \brief Conditional mean function type definition.
 * 
 * \param cond_t 
 * \param spec
 */
typedef void (*dynamo_cond_mean_t)(gsl_vector *cond_t, const gsl_vector_int *spec, const gsl_vector *param, const gsl_vector *y_l, const gsl_matrix *cond_l, const gsl_vector *x_l, int *flag);
typedef void (*dynamo_cond_var_t)(gsl_vector *cond_t, const gsl_vector_int *spec, const gsl_vector *param, const gsl_vector *y_l, const gsl_matrix *cond_l, const gsl_vector *z_l, int *flag);

/*!
 * \brief Random deviates functions array.
 * \sa dynamo_rand_t, DYNAMO_NEWS
 */
extern dynamo_rand_t dynamo_rand_vector[];
extern dynamo_quant_t dynamo_quant_vector[];
/*!
 * \brief Conditional log density function array. 
 * \sa dynamo_cond_loglik_t
 */
extern dynamo_cond_loglik_t dynamo_cond_loglik_vector[];
extern dynamo_cond_mean_t dynamo_cond_mean_vector[];
extern dynamo_cond_var_t dynamo_cond_var_vector[];

/*!
 * \brief Correct model specification check.
 *
 * \return 1 if the check is successfully passed, 0 otherwise.
 */
int dynamo_check_spec(char *msg, const gsl_vector_int *spec);

/*!
 * \brief Correct model input dimension check.
 *
 * The function checks for the correct array dimensions of the dynamo routines input vectors conditionally on the user defined dynamo specificaton.
 *
 * \return 1 if the check is successfully passed, 0 otherwise.
 */
int dynamo_check_dim(char *msg, const gsl_vector_int *spec, const gsl_vector *param, const gsl_vector *y, const gsl_matrix *cond, const gsl_matrix *X, const gsl_matrix *Z);

/*!
 * \brief Correct data check.
 *
 * \return 1 if the check is successfully passed, 0 otherwise.
 */
int dynamo_check_data(char *msg, int checktype, const gsl_vector_int *spec, const gsl_vector *param, const gsl_vector *y, const gsl_matrix *X, const gsl_matrix *Z);

/*!
 * \brief Max number of lags in the model.
 */ 
size_t dynamo_spec_maxlag(const gsl_vector_int *spec);

size_t dynamo_spec_mean_pmul(const gsl_vector_int *spec);
size_t dynamo_spec_var_pmul(const gsl_vector_int *spec);

size_t dynamo_spec_ndim_mean_x(const gsl_vector_int *spec); /**** ???? */

size_t dynamo_spec_mean_const(const gsl_vector_int *spec);
size_t dynamo_spec_mean_plag(const gsl_vector_int *spec, size_t i);
size_t dynamo_spec_mean_qlag(const gsl_vector_int *spec, size_t j);
size_t dynamo_spec_mean_rreg(const gsl_vector_int *spec, size_t k);
size_t dynamo_spec_mean_xreg(const gsl_vector_int *spec, size_t l);

size_t dynamo_spec_var_const(const gsl_vector_int *spec);
size_t dynamo_spec_var_plag(const gsl_vector_int *spec, size_t i);
size_t dynamo_spec_var_qlag(const gsl_vector_int *spec, size_t j);
size_t dynamo_spec_var_rreg(const gsl_vector_int *spec, size_t k);
size_t dynamo_spec_var_xreg(const gsl_vector_int *spec, size_t l);

size_t dynamo_spec_news(const gsl_vector_int *spec, size_t i);

/*!
 * \brief Normal random deviates.
 * \sa dynamo_rand_t
 */
void dynamo_rand_norm(gsl_vector *y_t, const gsl_vector_int *spec, const gsl_vector *param, const gsl_vector *cond_t, const gsl_rng *r, int *flag);
/*!
 * \brief Standardised Student t random deviates.
 * \sa dynamo_rand_t
 */
void dynamo_rand_st(gsl_vector *y_t, const gsl_vector_int *spec, const gsl_vector *param, const gsl_vector *cond_t, const gsl_rng *r, int *flag);
/*!
 * \brief Exponential random deviates.
 * \sa dynamo_rand_t
 */
void dynamo_rand_exp(gsl_vector *y_t, const gsl_vector_int *spec, const gsl_vector *param, const gsl_vector *cond_t, const gsl_rng *r, int *flag);
/*!
 * \brief Gamma random deviates.
 * \sa dynamo_rand_t
 */
void dynamo_rand_gamma(gsl_vector *y_t, const gsl_vector_int *spec, const gsl_vector *param, const gsl_vector *cond_t, const gsl_rng *r, int *flag);
/*!
 * \brief Weibull random deviates.
 * \sa dynamo_rand_t
 */
void dynamo_rand_weibull(gsl_vector *y_t, const gsl_vector_int *spec, const gsl_vector *param, const gsl_vector *cond_t, const gsl_rng *r, int *flag);

double dynamo_quant_norm(double p, const gsl_vector_int *spec, const gsl_vector *param, const gsl_vector *cond_t);
double dynamo_quant_st(double p, const gsl_vector_int *spec, const gsl_vector *param, const gsl_vector *cond_t);
double dynamo_quant_exp(double p, const gsl_vector_int *spec, const gsl_vector *param, const gsl_vector *cond_t);
double dynamo_quant_gamma(double p, const gsl_vector_int *spec, const gsl_vector *param, const gsl_vector *cond_t);
double dynamo_quant_weibull(double p, const gsl_vector_int *spec, const gsl_vector *param, const gsl_vector *cond_t);

/*!
 * \brief Conditionally normal loglikelihood.
 * \sa dynamo_cond_loglik_t
 * \todo Code a more efficient implementation!
 */
double dynamo_cond_loglik_norm( const gsl_vector *y_t, const gsl_vector_int *spec, const gsl_vector *param, const gsl_vector *cond_t, int *flag);
double dynamo_cond_loglik_exp( const gsl_vector *y_t, const gsl_vector_int *spec, const gsl_vector *param, const gsl_vector *cond_t, int *flag);
double dynamo_cond_loglik_gamma( const gsl_vector *y_t, const gsl_vector_int *spec, const gsl_vector *param, const gsl_vector *cond_t, int *flag);
/*!
 * \f[
 * Weibull\left(\frac{\Gamma\left(1+\frac{1}{\phi}\right)}{\mu_t},\phi \right)=
 *  \phi\left(\frac{\Gamma\left(1+\frac{1}{\phi}\right)}{\mu_t}\right)^{\phi} x^{\phi-1} \exp\left\{- \left(\frac{x_t\Gamma\left(1+\frac{1}{\phi}\right)}{\mu_t}\right)^{\phi} \right\}
 * \f]
 */
double dynamo_cond_loglik_weibull( const gsl_vector *y_t, const gsl_vector_int *spec, const gsl_vector *param, const gsl_vector *cond_t, int *flag);

void dynamo_uncond( gsl_vector *uncond, const gsl_vector_int *spec, const gsl_vector *param, const gsl_matrix *X, const gsl_matrix *Z);

void dynamo_cond_mean_none(gsl_vector *cond_t, const gsl_vector_int *spec, const gsl_vector *param, const gsl_vector *y_l, const gsl_matrix *cond_l, const gsl_vector *x_l, int *flag);
void dynamo_cond_mean_const(gsl_vector *cond_t, const gsl_vector_int *spec, const gsl_vector *param, const gsl_vector *y_l, const gsl_matrix *cond_l, const gsl_vector *x_l, int *flag);
void dynamo_cond_mean_arma(gsl_vector *cond_t, const gsl_vector_int *spec, const gsl_vector *param, const gsl_vector *y_l, const gsl_matrix *cond_l, const gsl_vector *x_l, int *flag);
void dynamo_cond_mean_mem(gsl_vector *cond_t, const gsl_vector_int *spec, const gsl_vector *param, const gsl_vector *y_l, const gsl_matrix *cond_l, const gsl_vector *x_l, int *flag);
void dynamo_cond_mean_amem(gsl_vector *cond_t, const gsl_vector_int *spec, const gsl_vector *param, const gsl_vector *y_l, const gsl_matrix *cond_l, const gsl_vector *x_l, int *flag);
void dynamo_cond_mean_cmem(gsl_vector *cond_t, const gsl_vector_int *spec, const gsl_vector *param, const gsl_vector *y_l, const gsl_matrix *cond_l, const gsl_vector *x_l, int *flag);
void dynamo_cond_mean_acmem(gsl_vector *cond_t, const gsl_vector_int *spec, const gsl_vector *param, const gsl_vector *y_l, const gsl_matrix *cond_l, const gsl_vector *x_l, int *flag);

void dynamo_cond_var_none(gsl_vector *cond_t, const gsl_vector_int *spec, const gsl_vector *param, const gsl_vector *y_l, const gsl_matrix *cond_l, const gsl_vector *z_l, int *flag);
void dynamo_cond_var_const(gsl_vector *cond_t, const gsl_vector_int *spec, const gsl_vector *param, const gsl_vector *y_l, const gsl_matrix *cond_l, const gsl_vector *z_l, int *flag);
void dynamo_cond_var_garch(gsl_vector *cond_t, const gsl_vector_int *spec, const gsl_vector *param, const gsl_vector *y_l, const gsl_matrix *cond_l, const gsl_vector *z_l, int *flag);
void dynamo_cond_var_gjrgarch(gsl_vector *cond_t, const gsl_vector_int *spec, const gsl_vector *param, const gsl_vector *y_l, const gsl_matrix *cond_l, const gsl_vector *z_l, int *flag);
void dynamo_cond_var_egarch(gsl_vector *cond_t, const gsl_vector_int *spec, const gsl_vector *param, const gsl_vector *y_l, const gsl_matrix *cond_l, const gsl_vector *z_l, int *flag);
void dynamo_cond_var_aparch(gsl_vector *cond_t, const gsl_vector_int *spec, const gsl_vector *param, const gsl_vector *y_l, const gsl_matrix *cond_l, const gsl_vector *z_l, int *flag);


// utils
double dynamo_cond_mean_pqsum(const gsl_vector_int *spec, const gsl_vector *param);
double dynamo_cond_var_pqsum(const gsl_vector_int *spec, const gsl_vector *param);


// explanatory variables management


// spline management functions
size_t dynamo_expl_spline_ndim(const gsl_vector_int *spec);
int dynamo_expl_spline_dim_check(const gsl_vector_int *spec, gsl_matrix *X, const gsl_vector *x, const gsl_vector *knots);
int dynamo_expl_spline_data_check(const gsl_vector_int *spec, gsl_matrix *X, const gsl_vector *x, const gsl_vector *knots);

/*!
 * \brief Set equidistant knots.
 */
void dynamo_expl_spline_knots_set(gsl_vector *knots, const gsl_vector_int *spec, double xl, double xr);

/*!
 * \brief Linear basis expansion.
 */
void dynamo_expl_spline_X_set(gsl_matrix *X, const gsl_vector_int *spec, const gsl_vector *x, const gsl_vector *knots);


// output
typedef int (*dynamo_printf_t)(const char *,...);

extern dynamo_printf_t dynamo_printf;

void dynamo_printf_set( dynamo_printf_t f );

#endif // DYNAMO_BASE_H

