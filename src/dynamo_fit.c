
#include "dynamo_fit.h"

dynamo_cond_loglik_t dynamo_cond_loglik_vector[] = { 
	dynamo_cond_loglik_norm, 
	dynamo_cond_loglik_norm, 
	dynamo_cond_loglik_norm, 
	dynamo_cond_loglik_norm, 
	dynamo_cond_loglik_exp, 
	dynamo_cond_loglik_gamma, 
	dynamo_cond_loglik_weibull
};

dynamo_cond_dloglik_t dynamo_cond_dloglik_vector[] = { 
	dynamo_cond_dloglik_norm, 
	dynamo_cond_dloglik_exp, 
	dynamo_cond_dloglik_exp, 
	dynamo_cond_dloglik_exp, 
	dynamo_cond_dloglik_exp, 
	dynamo_cond_dloglik_gamma, 
	dynamo_cond_dloglik_wei
};

dynamo_cond_dmean_t dynamo_cond_dmean_vector[] = {
	dynamo_cond_dmean_none,
	dynamo_cond_dmean_const,
	dynamo_cond_dmean_const,
	dynamo_cond_dmean_const,
	dynamo_cond_dmean_mem,
};

dynamo_cond_dvar_t dynamo_cond_dvar_vector[] = {
	dynamo_cond_dvar_none,
	dynamo_cond_dvar_const,
	dynamo_cond_dvar_garch
};

char *dynamo_fit_method_str[] = {
	"GARCH parameters initialisation",
	"MEM parameters initialisation",
	"Derivative Free Search",
	"BFGS"
};

dynamo_fit_method_t dynamo_fit_method_vector[] = {
	dynamo_fit_method_garch_init,
	dynamo_fit_method_mem_init,
	dynamo_fit_method_amoeba,
	dynamo_fit_method_bfgs
};

dynamo_fmin_t dynamo_fmin_vector[] = {
	dynamo_fmin_aloglik
};

double dynamo_fmin_aloglik(const gsl_vector *param, void *args)
{
	dynamo_fmin_args_t *a=args;
	return	-dynamo_loglik( a->cond, a->spec, param, a->y, a->X, a->Z);
}

void dynamo_dfmin_aloglik_ana(const gsl_vector *param, void *args, gsl_vector *g)
{
	double f;
	dynamo_fdfmin_aloglik_ana(param,args,&f,g);
}

void dynamo_fdfmin_aloglik_ana(const gsl_vector *param, void *args, double *f, gsl_vector *g)
{
	dynamo_fmin_args_t *a=args;
	dynamo_dloglik_ana( g, f, a->G, a->dcond, a->cond, a->spec, param, a->y, a->X, a->Z);
	*f = -(*f);
	gsl_vector_scale( g , -1.0 );
}

void dynamo_dfmin_aloglik_num(const gsl_vector *param, void *args, gsl_vector *g)
{
	double f;
	dynamo_fdfmin_aloglik_num(param,args,&f,g);
}

void dynamo_fdfmin_aloglik_num(const gsl_vector *param, void *args, double *f, gsl_vector *g)
{
	dynamo_fmin_args_t *a=args;
	dynamo_dloglik_num( g, f, a->G, a->dcond, a->cond, a->spec, param, a->y, a->X, a->Z, a->h);
	*f = -(*f);
	gsl_vector_scale( g , -1.0 );
}

void dynamo_fit_method_set(gsl_vector_int *met, size_t *nmet, const gsl_vector_int *spec, const gsl_vector_int *setup)
{
	*nmet = 0;

	if( gsl_vector_int_get(setup,DYNAMO_FIT_STARTVAL) == DYNAMO_STARTVAL_DEFAULT )
	{
		// initial values mean equation
		switch( gsl_vector_int_get(spec,DYNAMO_SPEC_MEAN) )
		{
			case DYNAMO_MEAN_MEM:
				gsl_vector_int_set( met, *nmet, DYNAMO_FIT_METHOD_MEM_INIT);
				*nmet += 1;
				break;
		}
	
		// initial values variance equation
		switch( gsl_vector_int_get(spec,DYNAMO_SPEC_VAR) )
		{
			case DYNAMO_VAR_GARCH:
				gsl_vector_int_set( met , *nmet , DYNAMO_FIT_METHOD_GARCH_INIT );
				*nmet += 1;
				break;
		}
	}

	// derivative free search
	gsl_vector_int_set( met , *nmet , DYNAMO_FIT_METHOD_AMOEBA );
	*nmet += 1;

	// optimisation
	switch( gsl_vector_int_get(setup,DYNAMO_FIT_OBJ) )
	{
		case DYNAMO_OBJ_LOGLIK:
		case DYNAMO_OBJ_PLOGLIK_RIDGE:
		case DYNAMO_OBJ_PLOGLIK_GRIDGE:
			gsl_vector_int_set( met , *nmet , DYNAMO_FIT_METHOD_BFGS );
			*nmet += 1;
			break;
	}
}

dynamo_dfmin_t dynamo_dfmin_set(const gsl_vector_int *spec, const gsl_vector_int *opts)
{
	return dynamo_dfmin_aloglik_ana;
}

dynamo_fdfmin_t dynamo_fdfmin_set(const gsl_vector_int *spec, const gsl_vector_int *opts)
{
	return dynamo_fdfmin_aloglik_ana;
}

int dynamo_fit(gsl_matrix *cond_fit, gsl_vector *param_fit, 
	gsl_vector *res, gsl_matrix *param_fit_avc, gsl_vector *gradient, double *obj,
	const gsl_vector_int *spec, const gsl_vector *y, const gsl_matrix *X, const gsl_matrix *Z, const gsl_vector_int *opts, char *info_msg)
{
	size_t i;
	size_t nobs, ndim;
	size_t m, nmet;
	int status=DYNAMO_FIT_STATUS_OK;
	dynamo_fmin_args_t args;
	gsl_vector_int *met;
	gsl_vector *param_cur;

	dynamo_fdfmin_t fdfmin;

	// check inputs
	if( !dynamo_check_spec( info_msg, spec) ) return DYNAMO_INFO_BADSPEC;
	if( !dynamo_check_dim( info_msg, spec, param_fit, y, cond_fit, X, Z) ) return DYNAMO_INFO_NONCONF;
	if( !dynamo_check_data( info_msg, DYNAMO_DATA_CHECK_FULL, spec, param_fit, y, X, Z) ) return DYNAMO_INFO_BADDATA;

	nobs = y->size;
	ndim = dynamo_spec_ndim(spec);

	// obj function arg initialisation
	args.spec = spec;
	args.cond = cond_fit;
	args.y = y;
	args.X = X;
	args.Z = Z;

	args.G = gsl_matrix_calloc( nobs, ndim);
	args.h = gsl_vector_calloc( ndim);

	// tmp
	gsl_vector_set_all(args.h,0.05);

	param_cur = gsl_vector_alloc( ndim);
	gsl_vector_memcpy( param_cur, param_fit);

	met = gsl_vector_int_alloc( DYNAMO_FIT_NMET_MAX ); 
	dynamo_fit_method_set( met, &nmet, spec, opts);

	for( m=0; m<nmet; ++m)
	{
		if( gsl_vector_int_get(opts,DYNAMO_FIT_LOGLEV) ) dynamo_printf(" - Step %d of %d: %s\n",m+1,nmet, dynamo_fit_method_str[ gsl_vector_int_get(met,m) ]);

		status = dynamo_fit_method_vector[ gsl_vector_int_get(met,m) ]( param_cur,  &args, opts);

		if( status == DYNAMO_FIT_STATUS_OK ) gsl_vector_memcpy( param_fit , param_cur );

		if( gsl_vector_int_get(opts,DYNAMO_FIT_LOGLEV) )
		{
			dynamo_printf("   Current estimate = ");
			for( i=0; i<dynamo_spec_ndim(spec); ++i) dynamo_printf("%.2g ", gsl_vector_get(param_fit,i));
			dynamo_printf("\n");
		}
	}

	// compute objetive function and gradient at param estimate
	fdfmin = dynamo_fdfmin_set( spec, opts);
	fdfmin( param_fit, &args, obj, gradient);

	// asy cov. matrix
	dynamo_fit_mle_vc_ana( param_fit_avc, param_fit, &args);
	dynamo_fit_res( spec, res, y, cond_fit);

	// cleaning up
	gsl_vector_free( param_cur );
	gsl_vector_int_free( met );

	gsl_matrix_free( args.G );
	gsl_vector_free( args.h );

	return DYNAMO_INFO_OK;
}

double dynamo_loglik(gsl_matrix *cond, const gsl_vector_int *spec, const gsl_vector *param, const gsl_vector *y, const gsl_matrix *X, const gsl_matrix *Z)
{
	dynamo_cond_mean_t dynamo_cond_mean;
	dynamo_cond_var_t dynamo_cond_var;
	dynamo_cond_loglik_t dynamo_cond_loglik;

	double loglik=0.0;
	size_t t;
	size_t nobs, nlag, ncnd;
	int goodcomp=1;
	
	gsl_matrix_view cond_l;
	gsl_vector_view cond_0, cond_t;

	gsl_vector_view y_t, y_l;
	gsl_vector_view x_l, z_l;

	//if( !dynamo_check_const_loglik(spec,param,y,X,Z) ) return -HUGE_VAL;
	
	nobs = y->size;
	nlag = dynamo_spec_maxlag(spec);
	ncnd = dynamo_spec_ncnd(spec);

	// set the process functions
	dynamo_cond_mean = dynamo_cond_mean_vector[DMS_M(spec)];
	dynamo_cond_var = dynamo_cond_var_vector[DMS_V(spec)];
	dynamo_cond_loglik = dynamo_cond_loglik_vector[DMS_N(spec)];

	// TODO: init using no params!
	if(nlag)
	{
		// non dip da param 
		cond_0 = gsl_matrix_row(cond,0);
		//dynamo_uncond( &cond_0.vector, spec, param, X, Z);

		// TODO: not very nice
		gsl_vector_set( &cond_0.vector, DMC_MUI(spec), gsl_stats_mean(y->data,1,nobs));
		gsl_vector_set( &cond_0.vector, DMC_SIG2I(spec), gsl_stats_variance(y->data,1,nobs));

		for(t=1;t<nlag;++t)
		{
			cond_t = gsl_matrix_row( cond, t);
			gsl_vector_memcpy( &cond_t.vector, &cond_0.vector);
		}
	}

	// compute loglik
	for(t=nlag; t<nobs; ++t)
	{
		// organize variables
		// TODO: not nice!
		y_t = gsl_vector_subvector( (gsl_vector *) y, t, 1);
		cond_t = gsl_matrix_row( cond, t);

		if( nlag )
		{
			// TODO: not nice!
			y_l = gsl_vector_subvector( (gsl_vector *) y, t-nlag, nlag);
			cond_l = gsl_matrix_submatrix( cond, t-nlag, 0, nlag, ncnd);
		}

		if( dynamo_spec_has_mean_reg(spec) ) x_l = gsl_matrix_row( (gsl_matrix *) X, t);
		if( dynamo_spec_has_var_reg(spec) )  z_l = gsl_matrix_row( (gsl_matrix *) Z, t);
	
		// start the action
		dynamo_cond_mean( &cond_t.vector, spec, param, &y_l.vector, &cond_l.matrix, &x_l.vector, &goodcomp);
		if( !goodcomp ) break;

		dynamo_cond_var( &cond_t.vector, spec, param, &y_l.vector, &cond_l.matrix, &z_l.vector, &goodcomp);
		if( !goodcomp ) break;

		// compute loglikelihood
		loglik += dynamo_cond_loglik( &y_t.vector, spec, param, &cond_t.vector, &goodcomp);
		if( !goodcomp ) break;
	}

	//printf("ll %e\n",loglik);

	return goodcomp?loglik:-HUGE_VAL;
}

void dynamo_dlogden_ana(gsl_matrix *G, double *l, gsl_matrix *dcond, gsl_matrix *cond, const gsl_vector_int *spec, const gsl_vector *param, const gsl_vector *y, const gsl_matrix *X, const gsl_matrix *Z)
{
	dynamo_cond_mean_t dynamo_cond_mean;
	dynamo_cond_var_t dynamo_cond_var;
	dynamo_cond_loglik_t dynamo_cond_loglik;
	dynamo_cond_dloglik_t dynamo_cond_dloglik;
	dynamo_cond_dmean_t dynamo_cond_dmean;
	dynamo_cond_dvar_t dynamo_cond_dvar;

	double loglik=0.0;
	size_t t;
	size_t nobs, nlag, ncnd;
	int goodcomp=1;
	
	gsl_matrix_view cond_l;
	gsl_vector_view cond_0, cond_t;

	gsl_vector_view y_t, y_l;
	gsl_vector_view x_l, z_l;

	gsl_vector_view g, g_mean, g_var, g_news;
	gsl_vector_view dmean_t, dvar_t;

	// ACTHUNG
	gsl_matrix *dmean=NULL, *dvar=NULL;
	if( dynamo_spec_ndim_mean(spec) ) dmean = gsl_matrix_calloc( dynamo_spec_maxlag(spec)+1, dynamo_spec_ndim_mean(spec));
	if( dynamo_spec_ndim_var(spec) ) dvar = gsl_matrix_calloc( dynamo_spec_maxlag(spec)+1, dynamo_spec_ndim_var(spec));

	//if( !dynamo_check_const_loglik(spec,param,y,X,Z) ) return -HUGE_VAL;
	
	nobs = y->size;
	nlag = dynamo_spec_maxlag(spec);
	ncnd = dynamo_spec_ncnd(spec);

	// set the process functions
	dynamo_cond_mean = dynamo_cond_mean_vector[DMS_M(spec)];
	dynamo_cond_var = dynamo_cond_var_vector[DMS_V(spec)];
	dynamo_cond_loglik = dynamo_cond_loglik_vector[DMS_N(spec)];

	dynamo_cond_dmean = dynamo_cond_dmean_vector[DMS_M(spec)];
	dynamo_cond_dvar = dynamo_cond_dvar_vector[DMS_V(spec)];
	dynamo_cond_dloglik = dynamo_cond_dloglik_vector[DMS_N(spec)];

	// TODO: init using no params!
	if(nlag)
	{
		// non dip da param 
		cond_0 = gsl_matrix_row(cond,0);
		//dynamo_uncond( &cond_0.vector, spec, param, X, Z);

		for(t=1;t<nlag;++t)
		{
			cond_t = gsl_matrix_row( cond, t);
			gsl_vector_memcpy( &cond_t.vector, &cond_0.vector);
		}
	}

	// compute loglik
	for(t=nlag; t<nobs; ++t)
	{
		// organize variables
		// TODO: not nice!
		y_t = gsl_vector_subvector( (gsl_vector *) y, t, 1);
		cond_t = gsl_matrix_row( cond, t);

		if( nlag )
		{
			// TODO: not nice!
			y_l = gsl_vector_subvector( (gsl_vector *) y, t-nlag, nlag);
			cond_l = gsl_matrix_submatrix( cond, t-nlag, 0, nlag, ncnd);
		}

		if( dynamo_spec_has_mean_reg(spec) ) x_l = gsl_matrix_row( (gsl_matrix *) X, t);
		if( dynamo_spec_has_var_reg(spec) )  z_l = gsl_matrix_row( (gsl_matrix *) Z, t);
	
		// start the action
		dynamo_cond_mean( &cond_t.vector, spec, param, &y_l.vector, &cond_l.matrix, &x_l.vector, &goodcomp);
		if( !goodcomp ) break;

		dynamo_cond_var( &cond_t.vector, spec, param, &y_l.vector, &cond_l.matrix, &z_l.vector, &goodcomp);
		if( !goodcomp ) break;

		// compute loglikelihood
		loglik += dynamo_cond_loglik( &y_t.vector, spec, param, &cond_t.vector, &goodcomp);
		if( !goodcomp ) break;

		g = gsl_matrix_row( G , t-nlag );
		if( dynamo_spec_ndim_mean(spec) ) g_mean = gsl_vector_subvector( &g.vector, 0, dynamo_spec_ndim_mean(spec));
		if( dynamo_spec_ndim_var(spec) ) g_var = gsl_vector_subvector( &g.vector, dynamo_spec_ndim_mean(spec), dynamo_spec_ndim_var(spec));
		if( dynamo_spec_ndim_news(spec) ) g_news = gsl_vector_subvector( &g.vector, dynamo_spec_ndim_mean(spec)+dynamo_spec_ndim_var(spec), dynamo_spec_ndim_news(spec));

		dynamo_cond_dloglik( &g_mean.vector, &g_var.vector, &g_news.vector, spec, param, &y_t.vector, &cond_t.vector, &goodcomp);

		dynamo_cond_dmean( dmean, spec, param, &y_l.vector, &cond_l.matrix, &x_l.vector, &z_l.vector, &goodcomp);
		dynamo_cond_dvar( dvar,  spec, param, &y_l.vector, &cond_l.matrix, &x_l.vector, &z_l.vector, &goodcomp);
	
		if( dynamo_spec_ndim_mean(spec) )
		{
			dmean_t = gsl_matrix_row( dmean , nlag );
			gsl_vector_mul( &g_mean.vector , &dmean_t.vector );
		}

		if( dynamo_spec_ndim_var(spec) )
		{
			dvar_t = gsl_matrix_row( dvar , nlag );
			gsl_vector_mul( &g_var.vector , &dvar_t.vector );
		}
	}

	*l = loglik;

	// ACTHUNG
	if( dynamo_spec_ndim_mean(spec) ) gsl_matrix_free( dmean );
	if( dynamo_spec_ndim_var(spec) ) gsl_matrix_free( dvar );

}

void dynamo_dlogden_num(gsl_matrix *G, double *l, gsl_matrix *dcond, gsl_matrix *cond, const gsl_vector_int *spec, const gsl_vector *param, const gsl_vector *y, const gsl_matrix *X, const gsl_matrix *Z, const gsl_vector *h)
{
	dynamo_cond_mean_t dynamo_cond_mean;
	dynamo_cond_var_t dynamo_cond_var;
	dynamo_cond_loglik_t dynamo_cond_loglik;

	double loglik=0.0, logden_plus, logden_minus;
	size_t i,t;
	size_t nobs, nlag, ncnd, ndim;
	int goodcomp=1;
	
	gsl_matrix_view cond_l;
	gsl_vector_view cond_0, cond_t;

	gsl_vector_view y_t, y_l;
	gsl_vector_view x_l, z_l;

	gsl_vector *param_plus, *param_minus;

	//if( !dynamo_check_const_loglik(spec,param,y,X,Z) ) return -HUGE_VAL;
	
	nobs = y->size;
	nlag = dynamo_spec_maxlag(spec);
	ncnd = dynamo_spec_ncnd(spec);
	ndim = dynamo_spec_ndim(spec);

	// TODO: remove from here!
	param_plus = gsl_vector_calloc(ndim);
	param_minus = gsl_vector_calloc(ndim);

	gsl_vector_memcpy(param_plus,param);
	gsl_vector_memcpy(param_minus,param);

	// set the process functions
	dynamo_cond_mean = dynamo_cond_mean_vector[DMS_M(spec)];
	dynamo_cond_var = dynamo_cond_var_vector[DMS_V(spec)];
	dynamo_cond_loglik = dynamo_cond_loglik_vector[DMS_N(spec)];

	// TODO: init using no params!
	if(nlag)
	{
		// non dip da param 
		cond_0 = gsl_matrix_row(cond,0);
		// TODO: not very nice
		gsl_vector_set( &cond_0.vector, DMC_MUI(spec), gsl_stats_mean(y->data,1,nobs));

		for(t=1;t<nlag;++t)
		{
			cond_t = gsl_matrix_row( cond, t);
			gsl_vector_memcpy( &cond_t.vector, &cond_0.vector);
		}
	}

	// compute loglik
	for(t=nlag; t<nobs; ++t)
	{
		// organize variables
		// TODO: not nice!
		y_t = gsl_vector_subvector( (gsl_vector *) y, t, 1);
		cond_t = gsl_matrix_row( cond, t);

		if( nlag )
		{
			// TODO: not nice!
			y_l = gsl_vector_subvector( (gsl_vector *) y, t-nlag, nlag);
			cond_l = gsl_matrix_submatrix( cond, t-nlag, 0, nlag, ncnd);
		}

		if( dynamo_spec_has_mean_reg(spec) ) x_l = gsl_matrix_row( (gsl_matrix *) X, t);
		if( dynamo_spec_has_var_reg(spec) )  z_l = gsl_matrix_row( (gsl_matrix *) Z, t);
	
		// start the action
		dynamo_cond_mean( &cond_t.vector, spec, param, &y_l.vector, &cond_l.matrix, &x_l.vector, &goodcomp);
		if( !goodcomp ) break;

		dynamo_cond_var( &cond_t.vector, spec, param, &y_l.vector, &cond_l.matrix, &z_l.vector, &goodcomp);
		if( !goodcomp ) break;

		// compute loglikelihood
		loglik += dynamo_cond_loglik( &y_t.vector, spec, param, &cond_t.vector, &goodcomp);
		if( !goodcomp ) break;

		for(i=0;i<ndim;++i)
		{	
			// set increments
			gsl_vector_set( param_plus, i, gsl_vector_get(param,i)+gsl_vector_get(h,i));
			gsl_vector_set( param_minus, i, gsl_vector_get(param,i)-gsl_vector_get(h,i));

			dynamo_cond_mean( &cond_t.vector, spec, param_plus, &y_l.vector, &cond_l.matrix, &x_l.vector, &goodcomp);
			dynamo_cond_mean( &cond_t.vector, spec, param_minus, &y_l.vector, &cond_l.matrix, &x_l.vector, &goodcomp);
			if( !goodcomp ) break;
			
			dynamo_cond_var( &cond_t.vector, spec, param_plus, &y_l.vector, &cond_l.matrix, &z_l.vector, &goodcomp);
			dynamo_cond_var( &cond_t.vector, spec, param_minus, &y_l.vector, &cond_l.matrix, &z_l.vector, &goodcomp);
			if( !goodcomp ) break;

			logden_plus = dynamo_cond_loglik( &y_t.vector, spec, param_plus, &cond_t.vector, &goodcomp);
			logden_minus = dynamo_cond_loglik( &y_t.vector, spec, param_minus, &cond_t.vector, &goodcomp);
			if( !goodcomp ) break;
			
			gsl_matrix_set( G, t, i, (logden_plus-logden_minus)/(2.0*gsl_vector_get(h,i)) );

			// reset increments
			gsl_vector_set( param_plus, i, gsl_vector_get(param,i)-gsl_vector_get(h,i));
			gsl_vector_set( param_minus, i, gsl_vector_get(param,i)+gsl_vector_get(h,i));
		}

		if( !goodcomp ) break;
	}

	*l = loglik;

	gsl_vector_free(param_plus);
	gsl_vector_free(param_minus);
}

void dynamo_dloglik_ana(gsl_vector *g, double *l, gsl_matrix *G, gsl_matrix *dcond, gsl_matrix *cond, const gsl_vector_int *spec, const gsl_vector *param, const gsl_vector *y, const gsl_matrix *X, const gsl_matrix *Z)
{
	gsl_vector *one = NULL; // bad!
	
	one = gsl_vector_alloc( y->size );
	gsl_vector_set_all(one,1.0);

	dynamo_dlogden_ana(G,l,dcond,cond,spec,param,y,X,Z);

	gsl_blas_dgemv(CblasTrans,1.0,G,one,0.0,g);

	gsl_vector_free(one);
}

void dynamo_dloglik_num(gsl_vector *g, double *l, gsl_matrix *G, gsl_matrix *dcond, gsl_matrix *cond, const gsl_vector_int *spec, const gsl_vector *param, const gsl_vector *y, const gsl_matrix *X, const gsl_matrix *Z, const gsl_vector *h)
{
	gsl_vector *one; // bad!
	
	one = gsl_vector_alloc( y->size );
	gsl_vector_set_all(one,1.0);

	dynamo_dlogden_num(G,l,dcond,cond,spec,param,y,X,Z,h);

	gsl_blas_dgemv(CblasTrans,1.0,G,one,0.0,g);

	gsl_vector_free(one);
}

double dynamo_cond_loglik_norm( const gsl_vector *y, const gsl_vector_int *spec, const gsl_vector *param, const gsl_vector *cond_t, int *flag)
{
	double ld = log( gsl_ran_gaussian_pdf( gsl_vector_get(y,0)-DMC_MU(spec,cond_t), sqrt(DMC_SIG2(spec,cond_t)) ) );
	*flag = gsl_finite(ld);
	return ld;
}

double dynamo_cond_loglik_exp( const gsl_vector *y, const gsl_vector_int *spec, const gsl_vector *param, const gsl_vector *cond_t, int *flag)
{
	//printf("mu %f\n",DMC_MC(spec,cond_t));

	double ld = log( gsl_ran_exponential_pdf( gsl_vector_get(y,0), DMC_MU(spec,cond_t) ) );
	*flag = gsl_finite(ld);
	//printf("%f > %f\n",log(gsl_ran_exponential_pdf( gsl_vector_get(y,0), DMC_MU(spec,cond_t) )),DMC_MU(spec,cond_t));
	return ld;
}

double dynamo_cond_loglik_gamma( const gsl_vector *y, const gsl_vector_int *spec, const gsl_vector *param, const gsl_vector *cond_t, int *flag)
{
	double ld = log( gsl_ran_gamma_pdf( gsl_vector_get(y,0), DMP_N(spec,param,0), DMC_MU(spec,cond_t)/DMP_N(spec,param,0)  ) );
	*flag = gsl_finite(ld);
	return ld;
}

double dynamo_cond_loglik_weibull(const gsl_vector *y, const gsl_vector_int *spec, const gsl_vector *param, const gsl_vector *cond_t, int *flag)
{
        double ld = log( gsl_ran_weibull_pdf( gsl_vector_get(y,0), pow(gsl_sf_gamma(1.0+1.0/DMP_N(spec,param,0))/DMC_MU(spec,cond_t),-1.0), DMP_N(spec,param,0)));
	*flag = gsl_finite(ld);
	return ld;
}

void dynamo_cond_dloglik_norm(gsl_vector *g_mean, gsl_vector *g_var, gsl_vector *g_news, const gsl_vector_int *spec, const gsl_vector *param, const gsl_vector *y_t, const gsl_vector *cond_t, int *flag)
{
	//gsl_vector_set_all( g_mean , (y_t-mu_t)/sigma2_t );
	gsl_vector_set_all( g_var , -0.5*(1.0/DMC_SIG2(spec,cond_t)) + (0.5*gsl_pow_2(gsl_vector_get(y_t,0)-DMC_MU(spec,cond_t)))/gsl_pow_2(DMC_SIG2(spec,cond_t)) );
}

/*
void dynamo_cond_dloglik_st(gsl_vector *g_mean, gsl_vector *g_var, gsl_vector *g_news, const gsl_vector_int *spec, const gsl_vector *param, double y_t, double mu_t, double sigma2_t)
{
	double z, eta, dnews, dvar;

	z = (y_t-mu_t)/sqrt(sigma2_t);
	eta = gsl_vector_get( param , dynamo_spec_news(spec,0) );

	dvar = ((eta+1.0)/(2.0*eta)) * ( 1.0 / (1.0 + (eta/(1.0-2.0*eta))*z*z ) ) * (eta*gsl_pow_2(y_t-mu_t))/(1.0-2.0*eta) * (1.0/(sigma2_t*sigma2_t)) - 1.0/(2.0*sigma2_t);

	gsl_vector_set_all( g_var , dvar );

	dnews = 1.0/(2*eta*(1.0-2.0*eta)) - (1.0/(2.0*eta*eta)) * ( gsl_sf_psi((eta+1.0)/(2.0*eta)) - gsl_sf_psi(1.0/(2.0*eta)) ) +
		- (eta+1.0)/(2*eta*(1.0-2.0*eta))*( (z*z) / (1.0-2*eta+eta*z*z) ) + 1.0/(2*eta*eta)*log( 1.0 + (eta/(1.0-2.0*eta))*z*z );

	gsl_vector_set_all( g_news , dnews ); 
}

void dynamo_cond_dloglik_sskt(gsl_vector *g_mean, gsl_vector *g_var, gsl_vector *g_news, const gsl_vector_int *spec, const gsl_vector *param, double y_t, double mu_t, double sigma2_t)
{
}
*/

void dynamo_cond_dloglik_exp(gsl_vector *g_mean, gsl_vector *g_var, gsl_vector *g_news, const gsl_vector_int *spec, const gsl_vector *param, const gsl_vector *y_t, const gsl_vector *cond_t, int *flag)
{
	gsl_vector_set_all( g_mean , ((gsl_vector_get(y_t,0)-DMC_MU(spec,cond_t))/gsl_pow_2(DMC_MU(spec,cond_t))) );
}

void dynamo_cond_dloglik_gamma(gsl_vector *g_mean, gsl_vector *g_var, gsl_vector *g_news, const gsl_vector_int *spec, const gsl_vector *param, const gsl_vector *y_t, const gsl_vector *cond_t, int *flag)
{
	double phi = DMP_N(spec,param,0); 
	double mu_t= DMC_MU(spec,cond_t);

	gsl_vector_set_all( g_mean , phi*((gsl_vector_get(y_t,0)-mu_t)/(mu_t*mu_t)) );
	gsl_vector_set_all( g_news , log(phi)+1-gsl_sf_psi(phi)+log(gsl_vector_get(y_t,0))-log(mu_t)-gsl_vector_get(y_t,0)/mu_t );
}

void dynamo_cond_dloglik_wei(gsl_vector *g_mean, gsl_vector *g_var, gsl_vector *g_news, const gsl_vector_int *spec, const gsl_vector *param, const gsl_vector *y_t, const gsl_vector *cond_t, int *flag)
{
	double phi = DMP_N(spec,param,0); 
	double mu_t= DMC_MU(spec,cond_t);

	gsl_vector_set_all( g_mean , (phi/mu_t)* (gsl_sf_exp(phi*(log(gsl_sf_gamma(1.0+1.0/phi))+log(gsl_vector_get(y_t,0)/mu_t)))-1.0) );
	gsl_vector_set_all( g_news , (1.0/phi)+(log(gsl_sf_gamma(1.0+1.0/phi))-(1.0/phi)*gsl_sf_psi(1.0+1.0/phi)+log(gsl_vector_get(y_t,0)/mu_t))*(1.0-gsl_sf_exp(phi*(log(gsl_sf_gamma(1.0+1.0/phi))+log(gsl_vector_get(y_t,0)/mu_t)))) );
}

void dynamo_cond_dmean_none( gsl_matrix *dmean, const gsl_vector_int *spec, const gsl_vector *param, const gsl_vector *y_l, const gsl_matrix *cond_l, const gsl_vector *x_l, const gsl_vector *z_l, int *flag)
{
}

void dynamo_cond_dmean_const( gsl_matrix *dmean, const gsl_vector_int *spec, const gsl_vector *param, const gsl_vector *y_l, const gsl_matrix *cond_l, const gsl_vector *x_l, const gsl_vector *z_l, int *flag)
{
}

void dynamo_cond_dmean_mem( gsl_matrix *dmean, const gsl_vector_int *spec, const gsl_vector *param, const gsl_vector *y_l, const gsl_matrix *cond_l, const gsl_vector *x_l, const gsl_vector *z_l, int *flag)
{
	size_t i,j;
	size_t p,q,r;

	p = DMS_MP(spec);
	q = DMS_MQ(spec);
	r = DMS_MR(spec);
	gsl_vector_view dmean_t, dmean_t_q;

	// house keeping
	for(i=0;i<dynamo_spec_maxlag(spec);++i)
	{
		for(j=0;j<dynamo_spec_ndim_mean(spec);++j)
		{
			gsl_matrix_set( dmean , i , j , gsl_matrix_get(dmean,i+1,j) );
		}
	}

	dmean_t = gsl_matrix_row( dmean , dynamo_spec_maxlag(spec) );

	// d omega
	gsl_vector_set(&dmean_t.vector,0,1.0);
	// d alpha
	for(i=0;i<p;++i) gsl_vector_set( &dmean_t.vector , dynamo_spec_mean_plag(spec,i) , gsl_vector_get( y_l, LAGV(y_l,i)) );
	// d beta
	for(i=0;i<q;++i) gsl_vector_set( &dmean_t.vector , dynamo_spec_mean_qlag(spec,i) , gsl_matrix_get( cond_l, LAGM(cond_l,i), DMC_MUI(spec)));
	// d gamma
	for(i=0;i<r;++i) gsl_vector_set( &dmean_t.vector  , dynamo_spec_mean_rreg(spec,i) , gsl_vector_get(x_l,i) );

	for(i=0;i<q;++i)
	{
		dmean_t_q = gsl_matrix_row( dmean , dynamo_spec_maxlag(spec)-1-i );

		gsl_blas_daxpy( DMP_MQ(spec,param,i) , &dmean_t_q.vector , &dmean_t.vector );
	}

}

void dynamo_cond_dvar_none( gsl_matrix *dvar, const gsl_vector_int *spec, const gsl_vector *param, const gsl_vector *y_l, const gsl_matrix *cond_l, const gsl_vector *x_t, const gsl_vector *z_t, int *flag)
{
}

void dynamo_cond_dvar_const( gsl_matrix *dvar, const gsl_vector_int *spec, const gsl_vector *param, const gsl_vector *y_l, const gsl_matrix *cond_l, const gsl_vector *x_t, const gsl_vector *z_t, int *flag)
{
}

void dynamo_cond_dvar_garch( gsl_matrix *dvar, const gsl_vector_int *spec, const gsl_vector *param, const gsl_vector *y_l, const gsl_matrix *cond_l, const gsl_vector *x_l, const gsl_vector *z_l, int *flag)
{
	size_t i,j;
	size_t p,q,r;

	p = DMS_VP(spec);
	q = DMS_VQ(spec);
	r = DMS_VR(spec);
	gsl_vector_view dvar_t, dvar_t_q;

	// house keeping
	for(i=0;i<dynamo_spec_maxlag(spec);++i)
	{
		for(j=0;j<dynamo_spec_ndim_var(spec);++j)
		{
			gsl_matrix_set( dvar , i , j , gsl_matrix_get(dvar,i+1,j) );
		}
	}

	dvar_t = gsl_matrix_row( dvar , dynamo_spec_maxlag(spec) );

	// d omega
	gsl_vector_set(&dvar_t.vector,0,1.0);

	// d alpha
	for(i=0;i<p;++i) gsl_vector_set( &dvar_t.vector , dynamo_spec_var_plag(spec,i) , gsl_pow_2( gsl_vector_get( y_l, LAGV(y_l,i)) - gsl_matrix_get( cond_l, LAGM(cond_l,i), DMC_MUI(spec))) );
	// d beta
	for(i=0;i<q;++i) gsl_vector_set( &dvar_t.vector , dynamo_spec_var_qlag(spec,i) , gsl_matrix_get( cond_l, LAGM(cond_l,i), DMC_SIG2I(spec)) );
	// d gamma
	for(i=0;i<r;++i) gsl_vector_set( &dvar_t.vector  , dynamo_spec_var_rreg(spec,i) , gsl_vector_get(z_l,i) );

	for(i=0;i<q;++i)
	{
		dvar_t_q = gsl_matrix_row( dvar , dynamo_spec_maxlag(spec)-1-i );

		gsl_blas_daxpy( DMP_VQ(spec,param,i) , &dvar_t_q.vector , &dvar_t.vector );
	}
}

// fit methods
int dynamo_fit_method_garch_init(gsl_vector *param, const dynamo_fmin_args_t *args, const gsl_vector_int *setup)
{
	int status;
	size_t i,j;
	size_t na=11, nb=11, nc=21;
	double alpha[] = { 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9 , 0.93 };
	double beta[] = { 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9 , 0.93 };
	double gamma[] = { -1.0 , -0.9 , -0.8 , -0.7 , -0.6 , -0.5 , -0.4 , -0.3 , -0.2 , -0.1 , 0.0 , 0.1 , 0.2 , 0.3 , 0.4 , 0.5 , 0.6 , 0.7 , 0.8 , 0.9 , 1.0 };
	size_t ai, bi, ci; 
	size_t ai_star=0, bi_star=0, ci_star=0;
	double loglik, loglik_star=0.0;
	double var;
	gsl_vector *z_mean=NULL, *ones_n=NULL;
	size_t p,q,r;
	double omega; 

	p = DMS_VP(args->spec);
	q = DMS_VQ(args->spec);
	r = DMS_VR(args->spec);
	var = gsl_stats_variance( gsl_vector_const_ptr(args->y,0), 1, args->y->size);

	if( r )
	{
		z_mean = gsl_vector_calloc( r);
		ones_n = gsl_vector_alloc( args->y->size);
		gsl_vector_set_all( ones_n, 1.0);
		gsl_blas_dgemv( CblasTrans, 1.0/((double) args->y->size), args->Z, ones_n, 0.0, z_mean);
	}

	if( p+q > 0 )
	{

		// initial omega
		gsl_vector_set( param , dynamo_spec_var_const(args->spec) , var*(1.0-(alpha[0]+beta[0])) );
		// initial alpha
		for( i=0; i<p ; ++i ) gsl_vector_set( param , dynamo_spec_var_plag(args->spec,i) , alpha[0]/( (double) p ) );
		// initial beta
		for( i=0; i<q ; ++i ) gsl_vector_set( param , dynamo_spec_var_qlag(args->spec,i) , beta[0]/( (double) q ) );
		// exp
		for( i=0; i<r ; ++i ) gsl_vector_set( param , dynamo_spec_var_rreg(args->spec,i) , 0.0 );

		// initial news
		switch( gsl_vector_int_get(args->spec,DYNAMO_SPEC_NEWS) )
		{
			case DYNAMO_NEWS_ST:
				break;
		}	

		// initial loglik star
		loglik_star = -dynamo_fmin_aloglik(param, (void *) args);

		for( ai=0; ai<na; ++ai)
		{
			for( bi=0; bi<nb && (alpha[ai]+beta[bi])<0.999; ++bi)
			{
				// set omega
				gsl_vector_set( param , dynamo_spec_var_const(args->spec) , var*(1.0-(alpha[ai]+beta[bi])) );
				// set alpha
				for( i=0; i<p ; ++i ) gsl_vector_set( param , dynamo_spec_var_plag(args->spec,i) , alpha[ai]/( (double) p ) );
				// set beta
				for( i=0; i<q ; ++i ) gsl_vector_set( param , dynamo_spec_var_qlag(args->spec,i) , beta[bi]/( (double) q ) );

				// compute loglik
				loglik = -dynamo_fmin_aloglik(param,(void *) args);
				
				// cmp
				if( loglik > loglik_star )
				{
					ai_star = ai;
					bi_star = bi;
					loglik_star = loglik;
				}

				/*
				printf("%f %f %f: %f\n",
					gsl_vector_get(param,dynamo_spec_var_const(args->spec)),
					gsl_vector_get(param,dynamo_spec_var_plag(args->spec,0)),
					gsl_vector_get(param,dynamo_spec_var_qlag(args->spec,0)),
					loglik
					);
				*/
			}
		}
		
		// final omega
		gsl_vector_set( param , dynamo_spec_var_const(args->spec) , var*(1.0-(alpha[ai_star]+beta[bi_star])) );
		// final alpha
		for( i=0; i<p ; ++i ) gsl_vector_set( param , dynamo_spec_var_plag(args->spec,i) , alpha[ai_star]/( (double) p ) );
		// final beta
		for( i=0; i<q ; ++i ) gsl_vector_set( param , dynamo_spec_var_qlag(args->spec,i) , beta[bi_star]/( (double) q ) );
	}
	else
	{
		gsl_vector_set( param, dynamo_spec_mean_const(args->spec) , var );
	}


	switch( DMS_N(args->spec) )
	{
		case DYNAMO_NEWS_ST:
			{
				double ekurt;
				dynamo_fmin_aloglik(param,(void *) args);
				for(ekurt=0,i=0;i<args->y->size;++i)
				{
					ekurt += gsl_pow_2( gsl_pow_2( gsl_vector_get(args->y,i)-gsl_matrix_get(args->cond,i,DMC_MUI(args->spec)))/gsl_matrix_get(args->cond,i,DMC_SIG2I(args->spec)) );
				}
				ekurt = ekurt/( (double) args->y->size ) - 3.0;
				gsl_vector_set( param , dynamo_spec_news(args->spec,0) , 1.0/(6.0/ekurt + 4.0) );
			}
	} 

	for( i=0; i<r ; ++i )
	{
		ci_star = 10;
		gsl_vector_set( param , dynamo_spec_var_rreg(args->spec,i) , gamma[ci_star] ); 
		omega=var*(1.0-alpha[ai_star]-beta[ai_star]);
		for(j=0;j<r;++j)
		{
			omega -= gsl_vector_get( param , dynamo_spec_var_rreg(args->spec,j) )*gsl_vector_get(z_mean,j);
		}
		gsl_vector_set( param , dynamo_spec_var_const(args->spec) , omega );
		loglik_star = -dynamo_fmin_aloglik(param, (void *) args);

		for( ci=0; ci<nc; ++ci )
		{
			gsl_vector_set( param , dynamo_spec_var_rreg(args->spec,i) , gamma[ci] ); 

			omega=var*(1.0-alpha[ai_star]-beta[ai_star]);
			for(j=0;j<r;++j)
			{
				omega -= gsl_vector_get( param , dynamo_spec_var_rreg(args->spec,j) )*gsl_vector_get(z_mean,j);
			}
			gsl_vector_set( param , dynamo_spec_var_const(args->spec) , omega );

			loglik = -dynamo_fmin_aloglik(param,(void *) args);

			if( loglik > loglik_star && gsl_finite(loglik)==1 )
			{
				ci_star = ci;
				loglik_star = loglik;
			}

			//printf("%f: o %f g %f m %f v %f\n",loglik,omega,gamma[ci],gsl_vector_get(z_mean,i),var);
		}

		gsl_vector_set( param , dynamo_spec_var_rreg(args->spec,i) , gamma[ci_star] );
		omega=var*(1.0-alpha[ai_star]-beta[ai_star]);
		for(j=0;j<r;++j)
		{
			omega -= gsl_vector_get( param , dynamo_spec_var_rreg(args->spec,j) )*gsl_vector_get(z_mean,j);
		}
		gsl_vector_set( param , dynamo_spec_var_const(args->spec) , omega );
	}

	if( r )
	{
		gsl_vector_free( z_mean );
		gsl_vector_free( ones_n );
	}

	loglik = -dynamo_fmin_aloglik(param,(void *) args);
	if( gsl_finite(loglik) ) status = DYNAMO_FIT_STATUS_OK;
	else status = DYNAMO_FIT_STATUS_FAILURE;

	return status;
}

int dynamo_fit_method_mem_init(gsl_vector *param, const dynamo_fmin_args_t *args, const gsl_vector_int *setup)
{
	size_t i,j;
	size_t nobs;
	size_t na=11, nb=11, nc=21;
	double alpha[] = { 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.95};
	double beta[] = { 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.95};
	double gamma[] = { -1.0 , -0.9, -0.8, -0.7, -0.6, -0.5, -0.4, -0.3, -0.2, -0.1, 0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0};
	double gamma_fac = 2.0;
	size_t ai, bi, ci; 
	size_t ai_star=0, bi_star=0, ci_star=0;
	double loglik, loglik_star=0.0;
	double mean;
	gsl_vector *x_mean=NULL, *ones_n=NULL;
	size_t p,q,r;
	double omega, ab;

	// 
	double gamma_phi;
	double weibull_upper_limit = 40;
	double weibull_step = 0.5;
	double weibull_guess, weibull_guess_star;


	p = DMS_MP(args->spec);
	q = DMS_MQ(args->spec);
	r = DMS_MR(args->spec);
	nobs = args->y->size;
	mean = gsl_stats_mean( gsl_vector_const_ptr(args->y,0), 1, args->y->size); 

	if( r )
	{
		x_mean = gsl_vector_calloc( r);
		ones_n = gsl_vector_alloc( args->y->size);
		gsl_vector_set_all( ones_n, 1.0);
		gsl_blas_dgemv( CblasTrans, 1.0/((double) args->y->size), args->X, ones_n, 0.0, x_mean);
	}

	// initial omega
	gsl_vector_set( param, dynamo_spec_mean_const(args->spec), mean);
	// initial alpha
	for( i=0; i<p; ++i) gsl_vector_set( param, dynamo_spec_mean_plag(args->spec,i), 0.0);
	// initial beta
	for( i=0; i<q; ++i) gsl_vector_set( param, dynamo_spec_mean_qlag(args->spec,i), 0.0);
	// initial gamma
	for( i=0; i<r; ++i) gsl_vector_set( param, dynamo_spec_mean_rreg(args->spec,i), 0.0);

	// init news
	switch( gsl_vector_int_get(args->spec,DYNAMO_SPEC_NEWS) )
	{
		case DYNAMO_NEWS_GAMMA:
		case DYNAMO_NEWS_WEIBULL:
			gsl_vector_set( param, dynamo_spec_news(args->spec,0), 1.0);
			break;
	}	

	// initial loglik star
	loglik_star = -dynamo_fmin_aloglik(param, (void *) args);

	// DEBUG
	//fprintf(stderr,"loglik initial: %f\n",loglik_star);

	for( ai=0; ai<na; ++ai)
	{
		for( bi=0; bi<nb && (alpha[ai]+beta[bi])<0.999; ++bi)
		{
			// set omega
			gsl_vector_set( param, dynamo_spec_mean_const(args->spec), mean*(1.0-(alpha[ai]+beta[bi])));
			// set alpha
			for( i=0; i<p; ++i) gsl_vector_set( param, dynamo_spec_mean_plag(args->spec,i), alpha[ai]/( (double) p ));
			// set beta
			for( i=0; i<q; ++i) gsl_vector_set( param, dynamo_spec_mean_qlag(args->spec,i), beta[bi]/( (double) q ));

			// compute loglik
			loglik = -dynamo_fmin_aloglik(param, (void *) args);
			
			// cmp
			if( loglik > loglik_star )
			{
				ai_star = ai;
				bi_star = bi;
				loglik_star = loglik;
			}
		}
	}

	// final omega
	gsl_vector_set( param, dynamo_spec_mean_const(args->spec), mean*(1.0-(alpha[ai_star]+beta[bi_star])));
	// final alpha
	for( i=0; i<p; ++i) gsl_vector_set( param , dynamo_spec_mean_plag(args->spec,i) , alpha[ai_star]/( (double) p ) );
	// final beta
	for( i=0; i<q; ++i) gsl_vector_set( param , dynamo_spec_mean_qlag(args->spec,i) , beta[bi_star]/( (double) q ) );
	// final pers.
	ab = alpha[ai_star]+beta[bi_star];


	// init innovations
	switch( DMS_N(args->spec) )
	{
		// init gamma with neil shephard's MoM estimator
		case DYNAMO_NEWS_GAMMA:
		{
			dynamo_fmin_aloglik(param,(void *) args);
			
			gamma_phi = 0;
			for( i=0; i<nobs; ++i)
			{
				gamma_phi += gsl_pow_2((gsl_vector_get(args->y,i)/(1e-4+gsl_matrix_get(args->cond,i,DMC_MUI(args->spec))))-1.0);
			}
			gamma_phi /= nobs;

			// in case anything went wrong
			if( !gsl_finite(gamma_phi) ) gamma_phi = 1.0;

			gsl_vector_set( param, dynamo_spec_news(args->spec,0), 1.0/gamma_phi);
		}
		break;

		// init weibull with grid search between weibull_step and weibull_upper_limit
		case DYNAMO_NEWS_WEIBULL:
		{
			weibull_guess_star = 1;
			loglik_star = -dynamo_fmin_aloglik(param, (void *) args);

			for( weibull_guess=weibull_step; weibull_guess<weibull_upper_limit; weibull_guess += weibull_step)
			{
				gsl_vector_set( param, dynamo_spec_news(args->spec,0), weibull_guess);
				loglik = -dynamo_fmin_aloglik(param, (void *) args);

				// cmp
				if( loglik > loglik_star )
				{
					weibull_guess_star = weibull_guess;
					loglik_star = loglik;
				}
			}

			gsl_vector_set( param, dynamo_spec_news(args->spec,0), weibull_guess_star);
		}
		break;
	}

	for( i=0; i<r ; ++i )
	{
		ci_star = 10;

		gsl_vector_set( param , dynamo_spec_mean_rreg(args->spec,i) , gamma_fac*gamma[ci_star] ); 
		omega=mean*(1.0-ab);
		for(j=0;j<r;++j)
		{
			omega -= gsl_vector_get( param , dynamo_spec_mean_rreg(args->spec,j) )*gsl_vector_get(x_mean,j);
		}
		gsl_vector_set( param , dynamo_spec_mean_const(args->spec) , omega );

		loglik_star = -dynamo_fmin_aloglik(param, (void *) args);

		for( ci=0; ci<nc; ++ci )
		{
			gsl_vector_set( param , dynamo_spec_mean_rreg(args->spec,i) , gamma_fac*gamma[ci] ); 

			omega=mean*(1.0-ab);
			for(j=0;j<r;++j)
			{
				omega -= gsl_vector_get( param , dynamo_spec_mean_rreg(args->spec,j) )*gsl_vector_get(x_mean,j);
			}
			gsl_vector_set( param , dynamo_spec_mean_const(args->spec) , omega );

			loglik = -dynamo_fmin_aloglik(param,(void *) args);

			if( loglik > loglik_star && gsl_finite(loglik)==1 )
			{
				ci_star = ci;
				loglik_star = loglik;
			}
		}

		gsl_vector_set( param , dynamo_spec_mean_rreg(args->spec,i) , gamma_fac*gamma[ci_star] );
		omega=mean*(1.0-ab);
		for(j=0;j<r;++j)
		{
			omega -= gsl_vector_get( param , dynamo_spec_mean_rreg(args->spec,j) )*gsl_vector_get(x_mean,j);
		}
		gsl_vector_set( param , dynamo_spec_mean_const(args->spec) , omega );
	}


	if(r)
	{
		gsl_vector_free( x_mean );
		gsl_vector_free( ones_n );
	}

	// compute loglik to check if everything is ok
	loglik = -dynamo_fmin_aloglik(param,(void *) args);

	return gsl_finite(loglik)? DYNAMO_FIT_STATUS_OK: DYNAMO_FIT_STATUS_FAILURE;
}


int dynamo_fit_method_amoeba(gsl_vector *param, const dynamo_fmin_args_t *args, const gsl_vector_int *fit)
{
	int status=GSL_CONTINUE;
	double fval;
	size_t ndim = param->size;  
	const gsl_multimin_fminimizer_type *T = gsl_multimin_fminimizer_nmsimplex;
	gsl_multimin_fminimizer *s = NULL;
	gsl_vector *ss;
	gsl_multimin_function minex_func;

	int loglev, maxiter;
	double ftol;

	size_t iter = 0, i;
	double size;

	loglev = gsl_vector_int_get(fit,DYNAMO_FIT_LOGLEV);
	maxiter = gsl_vector_int_get(fit,DYNAMO_FIT_MAXIT)==0? DYNAMO_FIT_MAXIT_DEFAULT : gsl_vector_int_get(fit,DYNAMO_FIT_MAXIT);
	ftol = gsl_vector_int_get(fit,DYNAMO_FIT_FTOL)==0? pow(10,-DYNAMO_FIT_FTOL_DEFAULT) : pow(10,-gsl_vector_int_get(fit,DYNAMO_FIT_FTOL));

	if( maxiter < 0 ) return DYNAMO_FIT_STATUS_OK;

	// step size 
	ss = gsl_vector_alloc(ndim);
	gsl_vector_set_all(ss,0.01);

	// Initialize method and iterate
	minex_func.f = dynamo_fmin_vector[ gsl_vector_int_get(fit,DYNAMO_FIT_OBJ)];
	minex_func.n = ndim;
	minex_func.params = (void *) args;

	s = gsl_multimin_fminimizer_alloc( T, ndim);
	gsl_multimin_fminimizer_set(s, &minex_func, param, ss);

	do{
		iter++;
		status = gsl_multimin_fminimizer_iterate(s);

		if (status) break;

		if( !gsl_finite(s->fval) ) break;

		size = gsl_multimin_fminimizer_size(s);
		status = gsl_multimin_test_size( size, ftol);

		switch( loglev )
		{
			case DYNAMO_LOGLEV_BASIC:
				dynamo_printf ("    - %4d)\n", iter);
				break;
			case DYNAMO_LOGLEV_DETAIL:
				dynamo_printf ("    - %4d) ", iter);
				for (i=0;i<ndim;i++) dynamo_printf ("%10.3g ", gsl_vector_get(s->x,i));
				dynamo_printf (" - f = %7.3f   size = %.3g\n",s->fval,size);
				break;
		}	
	}
	while ( status == GSL_CONTINUE && iter < maxiter );

	if( loglev )
	{
		if (status == GSL_SUCCESS) dynamo_printf ("    - Converged in %d iterations! :)\n",iter);
		else dynamo_printf ("    - Not converged after %d iterations! :(\n",maxiter);
	}
	
	// copy results
	gsl_vector_memcpy( param, s->x);

	// cleaning up
	gsl_vector_free(ss);
	gsl_multimin_fminimizer_free(s);

	// check
	fval = -dynamo_fmin_aloglik(param,(void *) args);
	if( gsl_finite(fval) ) status = DYNAMO_FIT_STATUS_OK;
	else status = DYNAMO_FIT_STATUS_FAILURE;

	return status;
}

int dynamo_fit_method_bfgs(gsl_vector *param, const dynamo_fmin_args_t *args, const gsl_vector_int *opts)
{
	int status;

	size_t ndim = param->size;  
	const gsl_multimin_fdfminimizer_type *T = gsl_multimin_fdfminimizer_vector_bfgs;
	gsl_multimin_fdfminimizer *s;
	gsl_multimin_function_fdf minex_func;

	int loglev;
	double ftol;

	size_t i;
	int iter = 0, maxiter;
	size_t tries = 0,maxtries=9;
	int status_opt, status_grd;

	double dxn;
	
	loglev = gsl_vector_int_get(opts,DYNAMO_FIT_LOGLEV);
	maxiter = gsl_vector_int_get(opts,DYNAMO_FIT_MAXIT)==0? DYNAMO_MAXIT_DEFAULT : gsl_vector_int_get(opts,DYNAMO_FIT_MAXIT);
	ftol = gsl_vector_int_get(opts,DYNAMO_FIT_FTOL)==0? pow(10,-DYNAMO_FTOL_DEFAULT) : pow(10,-gsl_vector_int_get(opts,DYNAMO_FIT_FTOL));

	if( maxiter < 0 ) return DYNAMO_FIT_STATUS_OK;

	// Initialize method and iterate
	minex_func.f = dynamo_fmin_vector[ gsl_vector_int_get(opts,DYNAMO_FIT_OBJ)];
	minex_func.df = dynamo_dfmin_set(args->spec,opts);
	minex_func.fdf = dynamo_fdfmin_set(args->spec,opts);
	minex_func.n = ndim;
	minex_func.params = (void *) args;

	s = gsl_multimin_fdfminimizer_alloc( T ,ndim);
	status_grd = GSL_CONTINUE;

	do{
		tries++;

		gsl_multimin_fdfminimizer_set(s, &minex_func, param, DYNAMO_STEP(tries), 1e-4);
		gsl_multimin_fdfminimizer_restart(s);

		do{
			status_opt = gsl_multimin_fdfminimizer_iterate(s);

			if ( status_opt ) break;

			iter++;
		
			if( !gsl_finite(s->f) ) break;

			status_grd = gsl_multimin_test_gradient(s->gradient,ftol);

			switch( loglev )
			{

				case DYNAMO_LOGLEV_BASIC:
					dynamo_printf("    - %4d)\n", iter);
					break;

				case DYNAMO_LOGLEV_DETAIL:
					dynamo_printf("    - %4d) ", iter);
					for (i=0;i<ndim;i++) dynamo_printf ("%10.3g ", gsl_vector_get(s->x,i));
					dynamo_printf(" - f = %7.3f   s = %.1e   |g| = %.3g\n", s->f, DYNAMO_STEP(tries), sqrt(gsl_blas_dnrm2(s->gradient)));
					break;
			}	

		}while ( status_grd == GSL_CONTINUE && iter < maxiter );

		//! \todo TODO: sure about this?
		dxn=0;
		for(i=0;i<ndim;++i) dxn = gsl_pow_2( gsl_vector_get(param,i)-gsl_vector_get(s->x,i));
		if( dxn < ftol ) status_grd = GSL_SUCCESS;

		gsl_vector_memcpy( param, s->x);

	}while( status_grd && tries < maxtries && iter < maxiter );

	if( loglev )
	{
		if (status_opt == GSL_SUCCESS) dynamo_printf ("    - Converged in %d iterations! :)\n",iter);
		else
		{
			if(iter==maxiter) dynamo_printf ("    - Not converged after %d iterations! :(\n",maxiter);
			else dynamo_printf ("    - Unable to find a minimum using a step size of %.1e! :(\n", DYNAMO_STEP(tries));
		}
	}
	
	// copy results
	gsl_vector_memcpy( param, gsl_multimin_fdfminimizer_x(s));

	// cleaing nup
	gsl_multimin_fdfminimizer_free(s);

	// check
	if( gsl_finite(s->f) ) status = DYNAMO_FIT_STATUS_OK;
	else status = DYNAMO_FIT_STATUS_FAILURE;

	return status;
}

void dynamo_fit_mle_vc_ana(gsl_matrix *vc, const gsl_vector *param, void *args)
{
	double ll;
	int s;
	dynamo_fmin_args_t *a=args;
	size_t ndim = dynamo_spec_ndim(a->spec);
	gsl_vector *g = gsl_vector_calloc( ndim);
	gsl_matrix *info = gsl_matrix_calloc( ndim, ndim);
	gsl_permutation *p = gsl_permutation_alloc(ndim);

	dynamo_dloglik_ana( g, &ll, a->G, a->dcond, a->cond, a->spec, param, a->y, a->X, a->Z);

	gsl_blas_dgemm(CblasTrans,CblasNoTrans,1.0,a->G,a->G,0.0,info);
	gsl_linalg_LU_decomp(info,p,&s);
	gsl_linalg_LU_invert(info,p,vc);

	gsl_vector_free( g );
	gsl_matrix_free( info );
	gsl_permutation_free( p );
}

void dynamo_fit_res( const gsl_vector_int *spec, gsl_vector *res, const gsl_vector *y, const gsl_matrix *cond)
{
	size_t i;
	size_t nobs = res->size;

	switch( DMS_N(spec) )
	{
		case DYNAMO_NEWS_NORM: 
		case DYNAMO_NEWS_ST:
			for(i=0;i<nobs;++i) gsl_vector_set( res, i, (gsl_vector_get(y,i)-gsl_matrix_get(cond,i,DMC_MUI(spec)))/sqrt(gsl_matrix_get(cond,i,DMC_SIG2I(spec))) );
			break;

		case DYNAMO_NEWS_EXP:
		case DYNAMO_NEWS_GAMMA:
		case DYNAMO_NEWS_WEIBULL:
			for(i=0;i<nobs;++i) gsl_vector_set( res, i, (gsl_vector_get(y,i)/gsl_matrix_get(cond,i,DMC_MUI(spec))));
			break;
	}
}

