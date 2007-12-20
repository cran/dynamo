
#include "dynamo_base.h"
#include "dynamo.h"

dynamo_rand_t dynamo_rand_vector[] = { dynamo_rand_norm, dynamo_rand_st, dynamo_rand_norm, dynamo_rand_norm, dynamo_rand_exp, dynamo_rand_gamma, dynamo_rand_weibull};

double dynamo_ran_st(const gsl_rng *r, double eta)
{
	return eta<1e-2? gsl_ran_gaussian(r,1.0) : (1.0/sqrt( 1.0/(1.0-2.0*eta) ))*gsl_ran_tdist(r,1.0/eta);
}

int
dynamo_sim(gsl_vector *y_sim, gsl_matrix *cond_sim, const gsl_vector_int *spec, const gsl_vector *param, const gsl_matrix *X, const gsl_matrix *Z, long *seed, char *info_msg)
{
	dynamo_cond_mean_t dynamo_cond_mean;
	dynamo_cond_var_t dynamo_cond_var;
	dynamo_rand_t dynamo_rand;

	size_t t;
	size_t nobs, nlag, ncnd;
	int goodcomp=1;

	gsl_matrix *cond_sim_ext;
	gsl_matrix_view cond_l, cond_sim_final;
	gsl_vector_view cond_0, cond_t;

	gsl_vector *y_sim_ext;
	gsl_vector_view y_t, y_l, y_sim_final;

	gsl_vector_view x_l, z_l;

	const gsl_rng_type * RT;
	gsl_rng *r;

	// check inputs
	if( !dynamo_check_spec( info_msg, spec) ) return DYNAMO_INFO_BADSPEC;
	if( !dynamo_check_dim( info_msg, spec, param, y_sim, cond_sim, X, Z) ) return DYNAMO_INFO_NONCONF;
	if( !dynamo_check_data( info_msg, DYNAMO_DATA_CHECK_EXPL, spec, param, y_sim, X, Z) ) return DYNAMO_INFO_BADDATA;

	nobs = y_sim->size;
	nlag = dynamo_spec_maxlag(spec);
	ncnd = dynamo_spec_ncnd(spec);

	// set the process functions
	dynamo_cond_mean = dynamo_cond_mean_vector[DMS_M(spec)];
	dynamo_cond_var = dynamo_cond_var_vector[DMS_V(spec)];
	dynamo_rand = dynamo_rand_vector[DMS_N(spec)];

	// random number initialisation
	gsl_rng_env_setup();
	RT = gsl_rng_default;
	r = gsl_rng_alloc(RT);
	gsl_rng_set(r,*seed);

	// allocate extended series and initialise elements to zero
	y_sim_ext = gsl_vector_calloc( nobs+nlag);
	cond_sim_ext = gsl_matrix_calloc( nobs+nlag, dynamo_spec_ncnd(spec));

	// init past of the series
	if(nlag)
	{
		cond_0 = gsl_matrix_row(cond_sim_ext,0);
		dynamo_uncond( &cond_0.vector, spec, param, X, Z);
		gsl_vector_set( y_sim_ext, 0, DMC_MU(spec,&cond_0.vector));

		for(t=1;t<nlag;++t)
		{
			gsl_vector_set( y_sim_ext, t, DMC_MU(spec,&cond_0.vector));

			cond_t = gsl_matrix_row( cond_sim_ext, t);
			gsl_vector_memcpy( &cond_t.vector, &cond_0.vector);
		}
	}
		
	// simulate
	for(t=nlag;t<nlag+nobs;++t)
	{
		// organize variables
		y_t = gsl_vector_subvector( y_sim_ext, t, 1);
		cond_t = gsl_matrix_row( cond_sim_ext, t);

		if( nlag )
		{
			y_l = gsl_vector_subvector( y_sim_ext, t-nlag, nlag);
			cond_l = gsl_matrix_submatrix( cond_sim_ext, t-nlag, 0, nlag, ncnd);
		}

		if( dynamo_spec_has_mean_reg(spec) ) x_l = gsl_matrix_row( (gsl_matrix *) X, t-nlag);
		if( dynamo_spec_has_var_reg(spec) )  z_l = gsl_matrix_row( (gsl_matrix *) Z, t-nlag);
	
		// start the action
		dynamo_cond_mean( &cond_t.vector, spec, param, &y_l.vector, &cond_l.matrix, &x_l.vector, &goodcomp);
		if( !goodcomp ) break;

		dynamo_cond_var( &cond_t.vector, spec, param, &y_l.vector, &cond_l.matrix, &z_l.vector, &goodcomp);
		if( !goodcomp ) break;

		dynamo_rand( &y_t.vector, spec, param, &cond_t.vector, r, &goodcomp);
		if( !goodcomp ) break;
	}

	// badcomp info management
	if( !goodcomp ) sprintf( info_msg, "an error occured in the computations at t=%d", t-nlag+1);

	// copy results
	y_sim_final = gsl_vector_subvector( y_sim_ext, nlag, nobs);
	gsl_vector_memcpy( y_sim, &y_sim_final.vector);

	cond_sim_final = gsl_matrix_submatrix( cond_sim_ext, nlag, 0, nobs, dynamo_spec_ncnd(spec));
	gsl_matrix_memcpy( cond_sim, &cond_sim_final.matrix);

	// free extended series and random number
	gsl_rng_free(r);
	gsl_vector_free(y_sim_ext);
	gsl_matrix_free(cond_sim_ext);

	return goodcomp ? DYNAMO_INFO_OK : DYNAMO_INFO_BADCOMP;
}

void dynamo_rand_norm(gsl_vector *y_t, const gsl_vector_int *spec, const gsl_vector *param, const gsl_vector *cond_t, const gsl_rng *r, int *flag)
{
	double y = DMC_MU(spec,cond_t) + gsl_ran_gaussian( r, sqrt(DMC_SIG2(spec,cond_t)));
	*flag = gsl_finite(y);
	if( *flag ) gsl_vector_set( y_t, 0, y);
}

void dynamo_rand_st(gsl_vector *y_t, const gsl_vector_int *spec, const gsl_vector *param, const gsl_vector *cond_t, const gsl_rng *r, int *flag)
{
	double y = DMC_MU(spec,cond_t) + sqrt(DMC_SIG2(spec,cond_t))*dynamo_ran_st( r, DMP_N(spec,param,0));
	*flag = gsl_finite(y);
	if( *flag ) gsl_vector_set( y_t, 0, y);
}

void dynamo_rand_exp(gsl_vector *y_t, const gsl_vector_int *spec, const gsl_vector *param, const gsl_vector *cond_t, const gsl_rng *r, int *flag)
{
	double y = gsl_ran_exponential( r, DMC_MU(spec,cond_t));
	*flag = gsl_finite(y);
	if( *flag ) gsl_vector_set( y_t, 0, y);
}

void dynamo_rand_gamma(gsl_vector *y_t, const gsl_vector_int *spec, const gsl_vector *param, const gsl_vector *cond_t, const gsl_rng *r, int *flag)
{
	double y = DMC_MU(spec,cond_t) * gsl_ran_gamma( r, DMP_N(spec,param,0), 1.0/DMP_N(spec,param,0));
	*flag = gsl_finite(y);
	if( *flag ) gsl_vector_set( y_t, 0, y);
}

void dynamo_rand_weibull(gsl_vector *y_t, const gsl_vector_int *spec, const gsl_vector *param, const gsl_vector *cond_t, const gsl_rng *r, int *flag)
{
        double y = gsl_ran_weibull( r, pow(gsl_sf_gamma(1.0+1.0/DMP_N(spec,param,0))/DMC_MU(spec,cond_t),-1.0), DMP_N(spec,param,0));
	*flag = gsl_finite(y);
	if( *flag ) gsl_vector_set( y_t, 0, y);
}
