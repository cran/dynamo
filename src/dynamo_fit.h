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

#ifndef DYNAMO_FIT_H
#define DYNAMO_FIT_H

#include "dynamo.h"
#include "dynamo_base.h"

#define DYNAMO_FIT_NMET_MAX 10
#define DYNAMO_FIT_MAXIT_DEFAULT 100
#define DYNAMO_FIT_FTOL_DEFAULT 4

#define DYNAMO_STEP(i) (pow(10.0,-(((double)i)+0.5)/1.5))

enum DYNAMO_FIT_METHOD { 
	DYNAMO_FIT_METHOD_GARCH_INIT=0, 
	DYNAMO_FIT_METHOD_MEM_INIT, 
	DYNAMO_FIT_METHOD_AMOEBA, 
	DYNAMO_FIT_METHOD_BFGS
};

enum DYNAMO_FIT_STATUS { 
	DYNAMO_FIT_STATUS_OK=0, 
	DYNAMO_FIT_STATUS_FAILURE 
};

/*!
 * \brief Optimization fmin arguments structure definition.
 */
typedef struct dynamo_fmin_args{

	const gsl_vector_int *spec;

	gsl_matrix *cond;
	gsl_matrix *dcond;
	gsl_matrix *dmean; 
	gsl_matrix *dvar;

	gsl_vector *h;
	gsl_matrix *G;

	const gsl_vector *y;
	const gsl_matrix *X;	
	const gsl_matrix *Z;

	const gsl_vector *pentune;
	const gsl_matrix *penpar;

} dynamo_fmin_args_t;


typedef void (*dynamo_cond_dloglik_t)(gsl_vector *g_mean, gsl_vector *g_var, gsl_vector *g_news, const gsl_vector_int *spec, const gsl_vector *param, const gsl_vector *y_t, const gsl_vector *cond_t, int *flag);

void dynamo_cond_dloglik_norm(gsl_vector *g_mean, gsl_vector *g_var, gsl_vector *g_news, const gsl_vector_int *spec, const gsl_vector *param, const gsl_vector *y_t, const gsl_vector *cond_t, int *flag);
void dynamo_cond_dloglik_exp(gsl_vector *g_mean, gsl_vector *g_var, gsl_vector *g_news, const gsl_vector_int *spec, const gsl_vector *param, const gsl_vector *y_t, const gsl_vector *cond_t, int *flag);
void dynamo_cond_dloglik_gamma(gsl_vector *g_mean, gsl_vector *g_var, gsl_vector *g_news, const gsl_vector_int *spec, const gsl_vector *param, const gsl_vector *y_t, const gsl_vector *cond_t, int *flag);

void dynamo_cond_dloglik_wei(gsl_vector *g_mean, gsl_vector *g_var, gsl_vector *g_news, const gsl_vector_int *spec, const gsl_vector *param, const gsl_vector *y_t, const gsl_vector *cond_t, int *flag);

typedef void (*dynamo_cond_dmean_t)( gsl_matrix *dmean, const gsl_vector_int *spec, const gsl_vector *param, const gsl_vector *y_l, const gsl_matrix *cond_l, const gsl_vector *x_t, const gsl_vector *z_t, int *flag);

void dynamo_cond_dmean_none( gsl_matrix *dmean, const gsl_vector_int *spec, const gsl_vector *param, const gsl_vector *y_l, const gsl_matrix *cond_l, const gsl_vector *x_t, const gsl_vector *z_t, int *flag);
void dynamo_cond_dmean_const( gsl_matrix *dmean, const gsl_vector_int *spec, const gsl_vector *param, const gsl_vector *y_l, const gsl_matrix *cond_l, const gsl_vector *x_t, const gsl_vector *z_t, int *flag);
void dynamo_cond_dmean_mem( gsl_matrix *dmean, const gsl_vector_int *spec, const gsl_vector *param, const gsl_vector *y_l, const gsl_matrix *cond_l, const gsl_vector *x_t, const gsl_vector *z_t, int *flag);

typedef void (*dynamo_cond_dvar_t)( gsl_matrix *dvar, const gsl_vector_int *spec, const gsl_vector *param, const gsl_vector *y_l, const gsl_matrix *cond_l, const gsl_vector *x_t, const gsl_vector *z_t, int *flag);
void dynamo_cond_dvar_none( gsl_matrix *dvar, const gsl_vector_int *spec, const gsl_vector *param, const gsl_vector *y_l, const gsl_matrix *cond_l, const gsl_vector *x_l, const gsl_vector *z_l, int *flag);
void dynamo_cond_dvar_const( gsl_matrix *dvar, const gsl_vector_int *spec, const gsl_vector *param, const gsl_vector *y_l, const gsl_matrix *cond_l, const gsl_vector *x_l, const gsl_vector *z_l, int *flag);
void dynamo_cond_dvar_garch( gsl_matrix *dvar, const gsl_vector_int *spec, const gsl_vector *param, const gsl_vector *y_l, const gsl_matrix *cond_l, const gsl_vector *x_l, const gsl_vector *z_l, int *flag);

/*!
 * \brief Fit method type definition.
 */
typedef int (*dynamo_fit_method_t)(gsl_vector *param, const dynamo_fmin_args_t *args, const gsl_vector_int *fit); 

typedef double (*dynamo_fmin_t)(const gsl_vector *param, void *args);
typedef void (*dynamo_dfmin_t)(const gsl_vector *param, void *args, gsl_vector * g);
typedef void (*dynamo_fdfmin_t)(const gsl_vector *param, void *args, double * f, gsl_vector * g);

double dynamo_fmin_aloglik(const gsl_vector *param, void *args);

void dynamo_dfmin_aloglik_ana(const gsl_vector *param, void *args, gsl_vector *g);
void dynamo_dfmin_aloglik_num(const gsl_vector *param, void *args, gsl_vector *g);

void dynamo_fdfmin_aloglik_ana(const gsl_vector *param, void *args, double *f, gsl_vector *g);
void dynamo_fdfmin_aloglik_num(const gsl_vector *param, void *args, double *f, gsl_vector *g);

void dynamo_fit_method_set(gsl_vector_int *met, size_t *nmet, const gsl_vector_int *spec, const gsl_vector_int *fit);

dynamo_fdfmin_t dynamo_fdfmin_set(const gsl_vector_int *spec, const gsl_vector_int *opts);
dynamo_dfmin_t dynamo_dfmin_set(const gsl_vector_int *spec, const gsl_vector_int *opts);

extern dynamo_cond_dmean_t dynamo_cond_dmean_vector[];
extern dynamo_cond_dvar_t dynamo_cond_dvar_vector[];

extern char *dynamo_fit_method_str[];
extern dynamo_fit_method_t dynamo_fit_method_vector[];
extern dynamo_fmin_t dynamo_fmin_vector[];

/*!
 * \brief Conditional log likelihood function.
 *
 * The function computes the conditional log likelihood function, that is the log likelihood conditional on the first 
 * \f$ n_\mathsf{lag} \f$ observations, where \f$ n_\mathsf{lag} \f$ is determined as the
 *
 * \return the log likelihood if no error is encountered in the computation, -HUGE_REAL otherwise
 */
double dynamo_loglik(gsl_matrix *cond, const gsl_vector_int *spec, const gsl_vector *param, const gsl_vector *y, const gsl_matrix *X, const gsl_matrix *Z);

/*!
 * \brief Analytic conditional log densities first derivatives.
 */
void dynamo_dlogden_ana(gsl_matrix *G, double *l, gsl_matrix *dcond, gsl_matrix *cond, const gsl_vector_int *spec, const gsl_vector *param, const gsl_vector *y, const gsl_matrix *X, const gsl_matrix *Z);

/*!
 * \brief Numerical conditional log densities first derivatives.
 */
void dynamo_dlogden_num(gsl_matrix *G, double *l, gsl_matrix *dcond, gsl_matrix *cond, const gsl_vector_int *spec, const gsl_vector *param, const gsl_vector *y, const gsl_matrix *X, const gsl_matrix *Z, const gsl_vector *h);

/*! 
 * \brief Numerical first derivatives of the log likelihood
 */
void dynamo_dloglik_num(gsl_vector *g, double *l, gsl_matrix *G, gsl_matrix *dcond, gsl_matrix *cond, const gsl_vector_int *spec, const gsl_vector *param, const gsl_vector *y, const gsl_matrix *X, const gsl_matrix *Z, const gsl_vector *h);

void dynamo_dloglik_ana(gsl_vector *g, double *l, gsl_matrix *G, gsl_matrix *dcond, gsl_matrix *cond, const gsl_vector_int *spec, const gsl_vector *param, const gsl_vector *y, const gsl_matrix *X, const gsl_matrix *Z);

// fit methods
int dynamo_fit_method_mem_init(gsl_vector *param, const dynamo_fmin_args_t *args, const gsl_vector_int *opts);
int dynamo_fit_method_garch_init(gsl_vector *param, const dynamo_fmin_args_t *args, const gsl_vector_int *setup);
int dynamo_fit_method_amoeba(gsl_vector *param, const dynamo_fmin_args_t *args, const gsl_vector_int *opts);
int dynamo_fit_method_bfgs(gsl_vector *param, const dynamo_fmin_args_t *args, const gsl_vector_int *opts);

// vc
void dynamo_fit_mle_vc_ana(gsl_matrix *vc, const gsl_vector *param, void *args);

void dynamo_fit_res( const gsl_vector_int *spec, gsl_vector *res, const gsl_vector *y, const gsl_matrix *cond);

#endif // DYNAMO_FIT_H 

