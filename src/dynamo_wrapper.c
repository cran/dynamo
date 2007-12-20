#include <R.h>
#include "dynamo.h"

void dynamo_sim_wrapper(int *info, char **info_msg, double *yR, double *condR, int *nobsR, int *ncndR, char **specattrR, int *specnumR, double *paramR, int *nparamR, double *XR, double *ZR, int *seedR)
{
	size_t nobs, ncnd;
	long seed;
	gsl_vector_int *spec;
	gsl_vector_view param;
	gsl_vector_view y;
	gsl_matrix_view cond, X, Z;

	// init vars
	spec = gsl_vector_int_calloc(DYNAMO_SPEC_FIELDS);
	gsl_vector_int_set( spec, DYNAMO_SPEC_MEAN, dynamo_string2const(specattrR[0])); 
	gsl_vector_int_set( spec, DYNAMO_SPEC_VAR, dynamo_string2const(specattrR[1])); 
	gsl_vector_int_set( spec, DYNAMO_SPEC_NEWS, dynamo_string2const(specattrR[2]));
	gsl_vector_int_set( spec, DYNAMO_SPEC_MEAN_P, specnumR[0]);
	gsl_vector_int_set( spec, DYNAMO_SPEC_MEAN_Q, specnumR[1]);
	gsl_vector_int_set( spec, DYNAMO_SPEC_MEAN_R, specnumR[2]);
	gsl_vector_int_set( spec, DYNAMO_SPEC_MEAN_X, specnumR[3]);
	gsl_vector_int_set( spec, DYNAMO_SPEC_VAR_P, specnumR[4]);
	gsl_vector_int_set( spec, DYNAMO_SPEC_VAR_Q, specnumR[5]);
	gsl_vector_int_set( spec, DYNAMO_SPEC_VAR_R, specnumR[6]);
	gsl_vector_int_set( spec, DYNAMO_SPEC_VAR_X, specnumR[7]);

	nobs = (size_t) *nobsR;	
	ncnd = dynamo_spec_ncnd(spec); *ncndR = ncnd; // this function call should be safe
	seed = (size_t) *seedR;

	param = gsl_vector_view_array( paramR, *nparamR);

	y = gsl_vector_view_array( yR, nobs);
	cond = gsl_matrix_view_array( condR, nobs, ncnd);

	if( dynamo_spec_has_mean_reg(spec) ) X = gsl_matrix_view_array( XR, nobs, DMS_MR(spec)); // TODO: fix API
	if( dynamo_spec_has_var_reg(spec) ) Z = gsl_matrix_view_array( ZR, nobs, DMS_VR(spec) ); // TODO: fix API

	// sim
	*info = dynamo_sim( &y.vector, &cond.matrix, spec, &param.vector, &X.matrix, &Z.matrix, &seed, *info_msg);

	// cleaning up
	gsl_vector_int_free( spec );
}

void dynamo_fit_wrapper(int *info, char **info_msg, double *yR, double *condR, double *residR, int *nobsR, int *ncndR, char **specattrR, int *specnumR, double *paramR, double *avcR, double *gradientR, double *obj, int *nparamR, double *XR, double *ZR, char **fitattrR, int *fitnumR)
{
	size_t nobs, ncnd;
	gsl_vector_int *spec, *fit;
	gsl_vector_view param, gradient;
	gsl_matrix_view avc;
	gsl_vector_view y, resid;
	gsl_matrix_view cond, X, Z;

	// init vars
	spec = gsl_vector_int_calloc(DYNAMO_SPEC_FIELDS);
	gsl_vector_int_set( spec, DYNAMO_SPEC_MEAN, dynamo_string2const(specattrR[0])); 
	gsl_vector_int_set( spec, DYNAMO_SPEC_VAR,  dynamo_string2const(specattrR[1])); 
	gsl_vector_int_set( spec, DYNAMO_SPEC_NEWS, dynamo_string2const(specattrR[2]));
	gsl_vector_int_set( spec, DYNAMO_SPEC_MEAN_P, specnumR[0]);
	gsl_vector_int_set( spec, DYNAMO_SPEC_MEAN_Q, specnumR[1]);
	gsl_vector_int_set( spec, DYNAMO_SPEC_MEAN_R, specnumR[2]);
	gsl_vector_int_set( spec, DYNAMO_SPEC_MEAN_X, specnumR[3]);
	gsl_vector_int_set( spec, DYNAMO_SPEC_VAR_P, specnumR[4]);
	gsl_vector_int_set( spec, DYNAMO_SPEC_VAR_Q, specnumR[5]);
	gsl_vector_int_set( spec, DYNAMO_SPEC_VAR_R, specnumR[6]);
	gsl_vector_int_set( spec, DYNAMO_SPEC_VAR_X, specnumR[7]);

	// fit
	fit = gsl_vector_int_calloc(DYNAMO_FIT_FIELDS);
	gsl_vector_int_set( fit, DYNAMO_FIT_OBJ, dynamo_string2const(fitattrR[0]));
	gsl_vector_int_set( fit, DYNAMO_FIT_MAXIT, fitnumR[0]);
	gsl_vector_int_set( fit, DYNAMO_FIT_STARTVAL, fitnumR[1]);
	gsl_vector_int_set( fit, DYNAMO_FIT_FTOL, fitnumR[2]);
	gsl_vector_int_set( fit, DYNAMO_FIT_LOGLEV, dynamo_string2const(fitattrR[1]));

	nobs = (size_t) *nobsR;	
	ncnd = dynamo_spec_ncnd(spec); *ncndR = ncnd; // this function call should be safe

	*nparamR = dynamo_spec_ndim(spec);
	param = gsl_vector_view_array( paramR, *nparamR);
	avc = gsl_matrix_view_array( avcR, *nparamR, *nparamR);
	gradient = gsl_vector_view_array( gradientR, *nparamR);

	y = gsl_vector_view_array( yR, nobs);
	cond = gsl_matrix_view_array( condR, nobs, ncnd);
	resid = gsl_vector_view_array( residR, nobs); 

	if( dynamo_spec_has_mean_reg(spec) ) X = gsl_matrix_view_array( XR, nobs, DMS_MR(spec)); // TODO: fix API
	if( dynamo_spec_has_var_reg(spec) ) Z = gsl_matrix_view_array( ZR, nobs, DMS_VR(spec) ); // TODO: fix API

	*info = dynamo_fit( &cond.matrix, &param.vector, &resid.vector, &avc.matrix, &gradient.vector, obj, spec, &y.vector, &X.matrix, &Z.matrix, fit, *info_msg);

	// cleaning up
	gsl_vector_int_free( spec );
	gsl_vector_int_free( fit );
}

void dynamo_pred_wrapper(int *info, char **info_msg, 
	char **predictionR, double *quantR,
	double *cond_predR, double *quant_predR,
	char **specattrR, int *specnumR, double *paramR, 
	double *yR, double *XR, double *ZR, int *nobsR, 
	double *yoR, double *XoR, double *ZoR, int *nforR)
{
	int prediction;
	double ql, qu;
	size_t nobs, nfor, ncnd;
	gsl_vector_int *spec;
	gsl_vector_view param;
	gsl_matrix_view cond_pred, quant_pred;
	gsl_vector_view y, yo;
	gsl_matrix_view X, Z, Xo, Zo;

	// init vars
	spec = gsl_vector_int_calloc(DYNAMO_SPEC_FIELDS);
	gsl_vector_int_set( spec, DYNAMO_SPEC_MEAN, dynamo_string2const(specattrR[0])); 
	gsl_vector_int_set( spec, DYNAMO_SPEC_VAR, dynamo_string2const(specattrR[1])); 
	gsl_vector_int_set( spec, DYNAMO_SPEC_NEWS, dynamo_string2const(specattrR[2]));
	gsl_vector_int_set( spec, DYNAMO_SPEC_MEAN_P, specnumR[0]);
	gsl_vector_int_set( spec, DYNAMO_SPEC_MEAN_Q, specnumR[1]);
	gsl_vector_int_set( spec, DYNAMO_SPEC_MEAN_R, specnumR[2]);
	gsl_vector_int_set( spec, DYNAMO_SPEC_MEAN_X, specnumR[3]);
	gsl_vector_int_set( spec, DYNAMO_SPEC_VAR_P,  specnumR[4]);
	gsl_vector_int_set( spec, DYNAMO_SPEC_VAR_Q,  specnumR[5]);
	gsl_vector_int_set( spec, DYNAMO_SPEC_VAR_R,  specnumR[6]);
	gsl_vector_int_set( spec, DYNAMO_SPEC_VAR_X,  specnumR[7]);

	prediction = dynamo_string2const(predictionR[0]);
	ql = quantR[0]; qu = quantR[1];

	nobs = (size_t) *nobsR;
	nfor = (size_t) *nforR;
	ncnd = dynamo_spec_ncnd(spec);

	cond_pred = gsl_matrix_view_array( cond_predR, nfor, ncnd);
	quant_pred = gsl_matrix_view_array( quant_predR, nfor, 2);

	param = gsl_vector_view_array( paramR, dynamo_spec_ndim(spec));

	y = gsl_vector_view_array( yR, nobs);
	if( prediction==DYNAMO_PRED_STATIC && nfor>1 ) yo = gsl_vector_view_array(yoR,nfor);

	if( dynamo_spec_has_mean_reg(spec) )
	{
		X  = gsl_matrix_view_array(  XR, nobs, DMS_MR(spec)); //! \todo TODO: fix API
		Xo = gsl_matrix_view_array( XoR, nfor, DMS_MR(spec)); //! \todo TODO: fix API
	}
	if( dynamo_spec_has_var_reg(spec) )
	{
		Z  = gsl_matrix_view_array(  ZR, nobs, DMS_VR(spec)); //! \todo TODO: fix API
		Zo = gsl_matrix_view_array( ZoR, nfor, DMS_VR(spec)); //! \todo TODO: fix API
	}
	
	*info = dynamo_pred( &cond_pred.matrix, &quant_pred.matrix, prediction, ql, qu,
	spec, &param.vector, &y.vector, &X.matrix, &Z.matrix, &yo.vector, &Xo.matrix, &Zo.matrix, *info_msg);

	// cleaning up
	gsl_vector_int_free(spec);
}

