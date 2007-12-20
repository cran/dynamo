
#include "dynamo_base.h"
#include "dynamo.h"

dynamo_quant_t dynamo_quant_vector[] = {
	dynamo_quant_norm,
	dynamo_quant_st,
	dynamo_quant_norm,
	dynamo_quant_norm,
	dynamo_quant_exp,
	dynamo_quant_gamma,
	dynamo_quant_weibull
};

int dynamo_pred(gsl_matrix *cond_pred, gsl_matrix *quant,
	int prediction, double ql, double qu, 
	const gsl_vector_int *spec, const gsl_vector *param, const gsl_vector *y, const gsl_matrix *X, const gsl_matrix *Z,
	const gsl_vector *yo, const gsl_matrix *Xo, const gsl_matrix *Zo,
	char *info_msg)
{
	dynamo_cond_mean_t dynamo_cond_mean;
	dynamo_cond_var_t dynamo_cond_var;
	dynamo_quant_t dynamo_quant;

	size_t t;
	size_t nobs, nlag, nfor, ncnd;
	int goodcomp=1;

	gsl_matrix *cond_ext;
	gsl_matrix_view cond_l, cond_final;
	gsl_vector_view cond_0, cond_t;

	gsl_vector *y_ext;
	gsl_vector_view y_t, y_l;

	gsl_matrix *X_ext = NULL, *Z_ext = NULL;
	gsl_matrix_view X_ext_in, X_ext_out, Z_ext_in, Z_ext_out;
	gsl_vector_view x_l, z_l;
	
	// check inputs
	//if( !dynamo_check_spec( info_msg, spec) ) return DYNAMO_INFO_BADSPEC;
	//if( !dynamo_check_dim( info_msg, spec, param, y_sim, cond_sim, X, Z) ) return DYNAMO_INFO_NONCONF;
	//if( !dynamo_check_data( info_msg, DYNAMO_DATA_CHECK_EXPL, spec, param, y_sim, X, Z) ) return DYNAMO_INFO_BADDATA;

	nobs = y->size;
	nlag = dynamo_spec_maxlag(spec);
	nfor = cond_pred->size1;
	ncnd = dynamo_spec_ncnd(spec);

	// set the process functions
	dynamo_cond_mean = dynamo_cond_mean_vector[DMS_M(spec)];
	dynamo_cond_var = dynamo_cond_var_vector[DMS_V(spec)];
	dynamo_quant = dynamo_quant_vector[DMS_N(spec)];

	// allocate extended series and initialise elements to zero
	y_ext = gsl_vector_calloc( nobs+nfor+nlag);
	cond_ext = gsl_matrix_calloc( nobs+nfor+nlag, dynamo_spec_ncnd(spec));

	// TODO: BUG: number of rows is equal to what?
	if( dynamo_spec_has_mean_reg(spec) ) X_ext = gsl_matrix_calloc( nobs+nfor, gsl_vector_int_get(spec,DYNAMO_SPEC_MEAN_R));
	if( dynamo_spec_has_var_reg(spec) ) Z_ext = gsl_matrix_calloc( nobs+nfor, gsl_vector_int_get(spec,DYNAMO_SPEC_VAR_R));

	// set data
	if( dynamo_spec_has_mean_reg(spec) )
	{
		X_ext_in = gsl_matrix_submatrix( X_ext, 0, 0, nobs, X_ext->size2);
		gsl_matrix_memcpy( &X_ext_in.matrix, X);
		X_ext_out = gsl_matrix_submatrix( X_ext, nobs, 0, nfor, X_ext->size2);
		gsl_matrix_memcpy( &X_ext_out.matrix, Xo);
	}

	if( dynamo_spec_has_var_reg(spec) )
	{
		Z_ext_in = gsl_matrix_submatrix( Z_ext, 0, 0, nobs, Z_ext->size2);
		gsl_matrix_memcpy( &Z_ext_in.matrix, Z);
		Z_ext_out = gsl_matrix_submatrix( Z_ext, nobs, 0, nfor, Z_ext->size2);
		gsl_matrix_memcpy( &Z_ext_out.matrix, Zo);
	}

	// init past of the series
	if(nlag)
	{
		cond_0 = gsl_matrix_row(cond_ext,0);
		dynamo_uncond( &cond_0.vector, spec, param, X, Z);
		gsl_vector_set( y_ext, 0, DMC_MU(spec,&cond_0.vector));

		for(t=1;t<nlag;++t)
		{
			gsl_vector_set( y_ext, t, DMC_MU(spec,&cond_0.vector));

			cond_t = gsl_matrix_row( cond_ext, t);
			gsl_vector_memcpy( &cond_t.vector, &cond_0.vector);
		}
	}

	for(t=nlag;t<nobs+nlag;++t)
	{
		gsl_vector_set( y_ext, t, gsl_vector_get(y,t-nlag) );
	}

	switch( prediction )
	{
		case DYNAMO_PRED_STATIC:
			for(t=nobs+nlag;t<nfor+nobs+nlag-1;++t) gsl_vector_set( y_ext, t, gsl_vector_get(yo,t-nobs-nlag));
			break;
		case DYNAMO_PRED_DYNAMIC:
			break;
	}

	// predict
	for(t=nlag;t<nlag+nobs+nfor;++t)
	{
		// organize variables
		y_t = gsl_vector_subvector( y_ext, t, 1);
		cond_t = gsl_matrix_row( cond_ext, t);

		if( nlag )
		{
			y_l = gsl_vector_subvector( y_ext, t-nlag, nlag);
			cond_l = gsl_matrix_submatrix( cond_ext, t-nlag, 0, nlag, ncnd);
		}

		if( dynamo_spec_has_mean_reg(spec) ) x_l = gsl_matrix_row( (gsl_matrix *) X_ext, t-nlag);
		if( dynamo_spec_has_var_reg(spec) )  z_l = gsl_matrix_row( (gsl_matrix *) Z_ext, t-nlag);
	
		// start the action
		dynamo_cond_mean( &cond_t.vector, spec, param, &y_l.vector, &cond_l.matrix, &x_l.vector, &goodcomp);
		if( !goodcomp ) break;

		dynamo_cond_var( &cond_t.vector, spec, param, &y_l.vector, &cond_l.matrix, &z_l.vector, &goodcomp);
		if( !goodcomp ) break;

		if( prediction == DYNAMO_PRED_DYNAMIC && t>=(nobs+nlag) ) gsl_vector_set( y_ext, t, DMC_MU(spec,&cond_t.vector));
	}

	// badcomp info management
	if( !goodcomp ) sprintf( info_msg, "an error occured in the computations at t=%d", t-nlag+1);

	cond_final = gsl_matrix_submatrix( cond_ext, nlag+nobs, 0, nfor, dynamo_spec_ncnd(spec));
	gsl_matrix_memcpy( cond_pred, &cond_final.matrix);

	for( t=0; t<nfor; ++t)
	{
		cond_t = gsl_matrix_row( cond_pred, t);
		gsl_matrix_set( quant, t, 0, dynamo_quant(ql,spec,param,&cond_t.vector));
		gsl_matrix_set( quant, t, 1, dynamo_quant(qu,spec,param,&cond_t.vector));
	}

	// cleaning up
	if(dynamo_spec_has_mean_reg(spec)) gsl_matrix_free(X_ext);
	if(dynamo_spec_has_var_reg(spec)) gsl_matrix_free(Z_ext);

	return goodcomp ? DYNAMO_INFO_OK : DYNAMO_INFO_BADCOMP;
}


double dynamo_quant_norm(double p, const gsl_vector_int *spec, const gsl_vector *param, const gsl_vector *cond_t)
{
	return DMC_MU(spec,cond_t)+gsl_cdf_gaussian_Pinv(p,sqrt(DMC_SIG2(spec,cond_t)));
}

double dynamo_quant_st(double p, const gsl_vector_int *spec, const gsl_vector *param, const gsl_vector *cond_t)
{
	//return DMC_MU(spec,cond_t)+sqrt(DMC_SIG2(spec,cond_t))*chdm_cdf_st_Pinv(p,NEWS(theta,spec,0));
	return 1;
}

double dynamo_quant_exp(double p, const gsl_vector_int *spec, const gsl_vector *param, const gsl_vector *cond_t)
{
	return gsl_cdf_exponential_Pinv( p, DMC_MU(spec,cond_t) );
}

double dynamo_quant_gamma(double p, const gsl_vector_int *spec, const gsl_vector *param, const gsl_vector *cond_t)
{
	return gsl_cdf_gamma_Pinv( p, DMP_N(spec,param,0), DMC_MU(spec,cond_t)/DMP_N(spec,param,0) );
}

double dynamo_quant_weibull(double p, const gsl_vector_int *spec, const gsl_vector *param, const gsl_vector *cond_t)
{
	return gsl_cdf_weibull_Pinv( p, pow(gsl_sf_gamma(1.0+1.0/DMP_N(spec,param,0))/DMC_MU(spec,cond_t),-1.0), DMP_N(spec,param,0));
}

int
dynamo_ipred(gsl_matrix *ipred, const gsl_vector_int *spec, const gsl_vector *theta, char *info_msg)
{
	return 0;
}

