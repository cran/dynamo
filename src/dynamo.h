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

//! \file dynamo.h Dynamo library main header file.
#ifndef DYNAMO_H
#define DYNAMO_H

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_matrix_int.h>

/** @defgroup DM_MAIN DynaMo 
 *
 * \brief Library users' API.
 *
 * @{
 */

#define DYNAMO_INFO_MSG_MAXLENGTH 1000

//! \brief Maximum size (number of fields) of the specification vector.
#define DYNAMO_SPEC_FIELDS 20

//! \brief Maximum size mean ``P'' order
#define DYNAMO_SPEC_MEAN_P_MAX 5
//! \brief Maximum size mean ``Q'' order
#define DYNAMO_SPEC_MEAN_Q_MAX 5
#define DYNAMO_SPEC_MEAN_R_MAX 30
//! \brief Maximum size mean ``P'' order
#define DYNAMO_SPEC_VAR_P_MAX 5
//! \brief Maximum size mean ``Q'' order
#define DYNAMO_SPEC_VAR_Q_MAX 5
#define DYNAMO_SPEC_VAR_R_MAX 30

//! \brief Maximum size (number of fields) of the setup vector.
#define DYNAMO_FIT_FIELDS 10



//! \brief Default maximum number of iterations.
#define DYNAMO_MAXIT_DEFAULT 100
//! \brief Default floating point tollerance.
#define DYNAMO_FTOL_DEFAULT 4


//! \brief Specification shortcut for mean field vaue.
#define DMS_M(spec) ( (size_t) gsl_vector_int_get(spec,DYNAMO_SPEC_MEAN) ) 
//! \brief Specification shortcut for variance field vaue.
#define DMS_V(spec) ( (size_t) gsl_vector_int_get(spec,DYNAMO_SPEC_VAR) )
//! \brief Specification shortcut for news field vaue.
#define DMS_N(spec)  ( (size_t) gsl_vector_int_get(spec,DYNAMO_SPEC_NEWS) )

//! \brief Specification shortcut for mean equation ``P'' order field value.
#define DMS_MP(spec) ( (size_t) gsl_vector_int_get(spec,DYNAMO_SPEC_MEAN_P) )
//! \brief Specification shortcut for mean equation ``Q'' order.
#define DMS_MQ(spec) ( (size_t) gsl_vector_int_get(spec,DYNAMO_SPEC_MEAN_Q) )
//! \brief Specification shortcut for mean equation number of regressors.
#define DMS_MR(spec) ( (size_t) gsl_vector_int_get(spec,DYNAMO_SPEC_MEAN_R) )
#define DMS_MX(spec) ( (size_t) gsl_vector_int_get(spec,DYNAMO_SPEC_MEAN_X) )
#define DMS_VP(spec) ( (size_t) gsl_vector_int_get(spec,DYNAMO_SPEC_VAR_P) )
#define DMS_VQ(spec) ( (size_t) gsl_vector_int_get(spec,DYNAMO_SPEC_VAR_Q) )
#define DMS_VR(spec) ( (size_t) gsl_vector_int_get(spec,DYNAMO_SPEC_VAR_R) )
#define DMS_VX(spec) ( (size_t) gsl_vector_int_get(spec,DYNAMO_SPEC_VAR_X) )

//! \brief Parameter shortcut for mean constant.
#define DMP_MC(spec,param) ( gsl_vector_get(param,0) )
#define DMP_MP(spec,param,i) ( gsl_vector_get(param,dynamo_spec_mean_plag(spec,i)) )
#define DMP_MQ(spec,param,i) ( gsl_vector_get(param,dynamo_spec_mean_qlag(spec,i)) )
#define DMP_MR(spec,param,i) ( gsl_vector_get(param,dynamo_spec_mean_rreg(spec,i)) )
#define DMP_MX(spec,param,i) ( gsl_vector_get(param,dynamo_spec_mean_xreg(spec,i)) )
#define DMP_VC(spec,param) ( gsl_vector_get(param,dynamo_spec_var_const(spec)) )
#define DMP_VP(spec,param,i) ( gsl_vector_get(param,dynamo_spec_var_plag(spec,i)) )
#define DMP_VQ(spec,param,i) ( gsl_vector_get(param,dynamo_spec_var_qlag(spec,i)) )
#define DMP_VR(spec,param,i) ( gsl_vector_get(param,dynamo_spec_var_rreg(spec,i)) )
#define DMP_VX(spec,param,i) ( gsl_vector_get(param,dynamo_spec_var_xreg(spec,i)) )
#define DMP_N(spec,param,i) ( gsl_vector_get(param,dynamo_spec_news(spec,i)) )

//! \brief Conditional parameter shortcut for mean.
#define DMC_MU(spec,cond) ( gsl_vector_get(cond,0) )
//! \brief Conditional parameter shortcut for variance.
#define DMC_SIG2(spec,cond) ( gsl_vector_get(cond,1) )

//! \brief Conditional parameter shortcut for mean index.
#define DMC_MUI(spec) ( 0 )
//! \brief Conditional parameter shortcut for variance index.
#define DMC_SIG2I(spec) ( 1 )

/*!
 * \brief Conditional parameter shortcut for tau (mean component) index.
 *
 * \sa DYNAMO_MEAN_CMEM, DYNAMO_MEAN_SMEM, DYNAMO_MEAN_BSMEM, DYNAMO_MEAN_ACMEM, DYNAMO_MEAN_ASMEM, DYNAMO_MEAN_ABSMEM
 */ 
#define DMC_MU_TAUI(spec) ( 2 )

/*!
 * \brief Conditional parameter shortcut for g (mean component) index.
 *
 * \sa DYNAMO_MEAN_CMEM, DYNAMO_MEAN_SMEM, DYNAMO_MEAN_BSMEM, DYNAMO_MEAN_ACMEM, DYNAMO_MEAN_ASMEM, DYNAMO_MEAN_ABSMEM
 */ 
#define DMC_MU_GI(spec) ( 3 )


//! \brief Specification vector map.
enum DYNAMO_SPEC {
	DYNAMO_SPEC_MEAN=0, /*!< Mean equation field */
	DYNAMO_SPEC_VAR, /*!< Variance equation field */
	DYNAMO_SPEC_NEWS, /*!< News equation field */
	DYNAMO_SPEC_MEAN_P, /*!< Mean equation ``P'' order field */
	DYNAMO_SPEC_MEAN_Q, /*!< Mean equation ``Q'' order field */
	DYNAMO_SPEC_MEAN_R, /*!< Mean equation field */
	DYNAMO_SPEC_MEAN_X, /*!< Mean equation extra component field */
	DYNAMO_SPEC_MEAN_K, /*!< Mean equation knots field */
	DYNAMO_SPEC_VAR_P, /*!< Variance equation ``P'' order field */
	DYNAMO_SPEC_VAR_Q, /*!< Variance equation ``Q'' order field */
	DYNAMO_SPEC_VAR_R,  /*!< Variance equation ``R'' order field */
	DYNAMO_SPEC_VAR_X,  /*!< Variance equation ``R'' order field */
	DYNAMO_SPEC_ATTR /*!< Attributes */
};

//! \brief Mean equation specifications.
enum DYNAMO_MEAN { 
	DYNAMO_MEAN_NONE=0, /*!< No mean equation. */
	DYNAMO_MEAN_CONST, /*!< Constant mean equation. */
	DYNAMO_MEAN_ARMA, /*!< ARMA equation with no constant. */
	DYNAMO_MEAN_CARMA, /*!< ARMA equation with constant. */
	DYNAMO_MEAN_MEM, /*!< MEM equation. */
	DYNAMO_MEAN_AMEM, /*!< Asymmetric MEM equation. */
	DYNAMO_MEAN_CMEM, /*!< Component MEM. */
	DYNAMO_MEAN_SMEM, /*!< Spline Component MEM. */
	DYNAMO_MEAN_BSMEM, /*!< B-Spline Component MEM. */
	DYNAMO_MEAN_ACMEM, /*!< Asymmetric Component MEM.*/
	DYNAMO_MEAN_ASMEM, /*!< Asymmetric Spline Component MEM.*/
	DYNAMO_MEAN_ABSMEM, /*!< Asymmetric B-SPline Component MEM.*/
};

//! \brief Variance equation specifications.
enum DYNAMO_VAR { 
	DYNAMO_VAR_NONE=0, /*!< No variance equation. */
	DYNAMO_VAR_CONST, /*!< Constant variance equation. */
	DYNAMO_VAR_GARCH, /*!< GARCH variance equation. */
	DYNAMO_VAR_GJRGARCH, /*!< GJR GARCH variance equation. */
	DYNAMO_VAR_EGARCH, /*!< EGARCH variance equation. */
	DYNAMO_VAR_APARCH /*!< APARCH variance equation. */
};

//! \brief Innovations distributions.
enum DYNAMO_NEWS { 
	DYNAMO_NEWS_NORM=0, /*!< Normal innovations. */
	DYNAMO_NEWS_ST, /*!< Standardised Student t innovations. */
	DYNAMO_NEWS_SSKT, /*!< Standardised skewed Student t innovations */ 
	DYNAMO_NEWS_GED, /*!< GED innovations */
	DYNAMO_NEWS_EXP, /*!< Unit Exponential Innovations */
	DYNAMO_NEWS_GAMMA, /*!< Unit Gamma Innovations */ 
	DYNAMO_NEWS_WEIBULL /*!< Unit Weibull Innovations */
};

//! \brief Fit setup.
enum DYNAMO_FIT {
	DYNAMO_FIT_OBJ=0, 
	DYNAMO_FIT_MAXIT, 
	DYNAMO_FIT_STARTVAL, 
	DYNAMO_FIT_FTOL, 
	DYNAMO_FIT_LOGLEV
};

//! \brief Objective Function
enum DYNAMO_OBJ { 
	DYNAMO_OBJ_LOGLIK=0, /*!< Log Likelihood */ 
	DYNAMO_OBJ_PLOGLIK_RIDGE, /*!< Penalised Log Likelihood with Ridge Penalty */ 
	DYNAMO_OBJ_PLOGLIK_GRIDGE /*!< Penalised Log Likelihood with Generalised Ridge Penalty */ 
};

enum DYNAMO_STARTVAL { 
	DYNAMO_STARTVAL_DEFAULT=0, /*!< Use default starting values. */
	DYNAMO_STARTVAL_USERDEFINED  /*!< Use user defined starting values. */
};

//! \brief Log messages depth.
enum DYNAMO_LOGLEV { 
	DYNAMO_LOGLEV_NONE=0,/*!< No log messages. */
	DYNAMO_LOGLEV_BASIC, /*!< Basic log messages. */
	DYNAMO_LOGLEV_DETAIL /*!< Detailed log messages. */
};

//! \brief Routines return info code.
enum DYNAMO_INFO { 
	DYNAMO_INFO_OK=0, /*!< No problem during routine execution. */
	DYNAMO_INFO_BADSPEC, /*!< The model specification is not valid. */
	DYNAMO_INFO_NONCONF, /*!< The model arguments are not conformable. */
	DYNAMO_INFO_BADDATA, /*!< The input data contains illegal values. */
	DYNAMO_INFO_BADPARAM, /*!< The model parameters are illegal. */
	DYNAMO_INFO_BADCOMP /*!< A numerical error occured during the computations. */
};

//! \brief Prediction type.
enum DYNAMO_PRED { DYNAMO_PRED_STATIC=0 , DYNAMO_PRED_DYNAMIC };

//! \brief Iterated prediction type.
enum DYNAMO_IPRED { DYNAMO_IPRED_FIXED=0 , DYNAMO_IPRED_RECURSIVE , DYNAMO_IPRED_ROLLING };

//! \brief Mean equation labels.
extern char *dynamo_mean_str[];
//! \brief Variance equation labels.
extern char *dynamo_var_str[];
//! \brief Innovations distributions labels.
extern char *dynamo_news_str[];
//! \brief Info codes labels.
extern char *dynamo_info_str[];

/*!
 * \brief Number of parameters.
 */ 
size_t dynamo_spec_ndim(const gsl_vector_int *spec);
size_t dynamo_spec_ndim_mean(const gsl_vector_int *spec);
size_t dynamo_spec_ndim_var(const gsl_vector_int *spec);
size_t dynamo_spec_ndim_news(const gsl_vector_int *spec);

/*!
 * \brief Number of variables (columns) of the X mean explanatory variables matrix.
 */
size_t dynamo_spec_nvar_x(const gsl_vector_int *spec);
/*!
 * \brief Number of variables (columns) of the Z var explanatory variables matrix.
 */
size_t dynamo_spec_nvar_z(const gsl_vector_int *spec);

int dynamo_spec_has_mean_const(const gsl_vector_int *spec);
int dynamo_spec_has_var_const(const gsl_vector_int *spec);
int dynamo_spec_has_mean_reg(const gsl_vector_int *spec);
int dynamo_spec_has_var_reg(const gsl_vector_int *spec);
/*!
 * \brief Number of conditional parameters series.
 */ 
size_t dynamo_spec_ncnd(const gsl_vector_int *spec);

/*!
 * \brief Dynamic model simulation.
 * 
 * \param y_sim 
 * \param cond_sim 
 * \param spec specification vector
 * \param I binary indicator vector for asymetric indicator models  
 * \param X matrix of mean equation regressors
 * \param Z matrix of variance equation regressors
 * \param seed random number generator seed number
 * \param info_msg 
 *
 * \return Dynamo info code.
 *
 * \sa DYNAMO_INFO
 */
int dynamo_sim(gsl_vector *y_sim, gsl_matrix *cond_sim, const gsl_vector_int *spec, const gsl_vector *param, const gsl_matrix *X, const gsl_matrix *Z, long *seed, char *info_msg);

/*!
 * \brief Dynamic model fit.
 *
 * \return Dynamo info code.
 *
 * \sa DYNAMO_INFO
 */
int dynamo_fit(gsl_matrix *cond_fit, gsl_vector *param_fit, 
	gsl_vector *res, gsl_matrix *param_fit_avc, gsl_vector *gradient, double *obj,
	const gsl_vector_int *spec, const gsl_vector *y, const gsl_matrix *X, const gsl_matrix *Z, const gsl_vector_int *opts, char *info_msg);

/*!
 * \brief Dynamic model regularization.
 *
 * \return Dynamo info code.
 *
 * \sa DYNAMO_INFO
 */
int dynamo_reg(gsl_matrix *cond_fit, gsl_vector *param_fit, const gsl_vector_int *spec, const gsl_vector *y, char *info_msg);

/*!
 * \brief Dynamic model prediction.
 *
 * \return Dynamo info code.
 *
 * \sa DYNAMO_INFO
 */
int dynamo_pred(gsl_matrix *cond_pred, gsl_matrix *quant,
	int prediction, double ql, double qu,
	const gsl_vector_int *spec, const gsl_vector *param, const gsl_vector *y, const gsl_matrix *X, const gsl_matrix *Z,
	const gsl_vector *yo, const gsl_matrix *Xo, const gsl_matrix *Zo,
	char *info_msg);

/*!
 * \brief Dynamic model iterative prediction.
 *
 * \return Dynamo info code.
 *
 * \sa DYNAMO_INFO
 */
int dynamo_ipred(gsl_matrix *cond_ipred, const gsl_vector_int *spec, const gsl_vector *param, char *info_msg);

/*!
 * \brief Mean explanatory variables generation.
 *
 * \sa DYNAMO_INFO
 */
int dynamo_expl_X_gen( gsl_matrix *X, const gsl_vector_int *spec, const gsl_matrix *Xin);

/*!
 * \brief Variance explanatory variables generation.
 *
 * \sa DYNAMO_INFO
 */
int dynamo_expl_Z_gen( gsl_matrix *Z, const gsl_vector_int *spec, const gsl_matrix *Zin);

// utils
int dynamo_string2const(char *str);

/** @} */

#endif // DYNAMO_H

