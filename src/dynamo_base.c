
#include "dynamo.h"
#include "dynamo_base.h"

char *dynamo_mean_str[] = 
{
	"None",
	"Cons",
	"ARMA",
	"ARMA with Constant",
	"MEM",
	"Asymmetric MEM",
	"Component MEM",
	"Spline MEM",
	"B-Spline MEM",
	"Asymmetric Component MEM"
};

char *dynamo_var_str[] = 
{
	"None",
	"Constant",
	"GARCH",
	"GJRGARCH",
	"EGARCH",
	"APARCH"
};

char *dynamo_news_str[] = 
{
	"Normal",
	"Standardised Student t",
	"Skewed Standardised Student t",
	"Exponential",
	"Gamma",
	"Weibull"
};

char *dynamo_info_str[] = 
{
	"OK",
	"Bad Specification",
	"Non Conformable Arguments",
	"Bad Data",
	"Bad Parameters",
	"Bad Computation"
};

dynamo_cond_mean_t dynamo_cond_mean_vector[] = 
{ 
	dynamo_cond_mean_none, 	
	dynamo_cond_mean_const, 
	dynamo_cond_mean_arma, 
	dynamo_cond_mean_arma, 
	dynamo_cond_mean_mem, 
	dynamo_cond_mean_amem,
	dynamo_cond_mean_cmem,
	dynamo_cond_mean_cmem, 
	dynamo_cond_mean_cmem 
};

dynamo_cond_var_t dynamo_cond_var_vector[] = 
{
	dynamo_cond_var_none, 
	dynamo_cond_var_const,
	dynamo_cond_var_garch, 
	dynamo_cond_var_gjrgarch 
};

int dynamo_check_spec(char *info_msg, const gsl_vector_int *spec)
{
	int goodspec=0;
	int mean, var, news;

	mean = gsl_vector_int_get( spec , DYNAMO_SPEC_MEAN );
	var = gsl_vector_int_get( spec , DYNAMO_SPEC_VAR );
	news = gsl_vector_int_get( spec , DYNAMO_SPEC_NEWS );

	// 1) check spec size
	if( spec->size < DYNAMO_SPEC_FIELDS )
	{
		sprintf(info_msg,"illegal specification vector size (%d instead of %d)",spec->size,DYNAMO_SPEC_FIELDS);
		return goodspec;
	}

	// 2) check if value are in range

	// 3) check model specification
	// ARMA GARCH like
	if( 	( mean==DYNAMO_MEAN_NONE || mean==DYNAMO_MEAN_CONST || mean==DYNAMO_MEAN_ARMA || mean==DYNAMO_MEAN_CARMA ) && 
		( var==DYNAMO_VAR_CONST  || var==DYNAMO_VAR_GARCH || var==DYNAMO_VAR_GJRGARCH || var==DYNAMO_VAR_EGARCH || var==DYNAMO_VAR_APARCH ) &&
		( news==DYNAMO_NEWS_NORM || news==DYNAMO_NEWS_ST || news==DYNAMO_NEWS_SSKT || news==DYNAMO_NEWS_GED )	)
	{

		// check if X field is not used

		goodspec = 1;
	}

	// MEM like
	if( 	( 
			mean==DYNAMO_MEAN_CONST || mean==DYNAMO_MEAN_MEM || mean==DYNAMO_MEAN_CMEM || mean==DYNAMO_MEAN_SMEM || mean==DYNAMO_MEAN_BSMEM 
			|| mean==DYNAMO_MEAN_AMEM || mean==DYNAMO_MEAN_ACMEM || mean==DYNAMO_MEAN_ASMEM || mean==DYNAMO_MEAN_ABSMEM 
		) && 
		( var==DYNAMO_VAR_NONE ) &&
		( news==DYNAMO_NEWS_EXP || news==DYNAMO_NEWS_GAMMA || news==DYNAMO_NEWS_WEIBULL )	)
	{
		goodspec = 1;
	}

	// componente MEM like
	

	if( !goodspec )
	{
		sprintf(info_msg,"illegal triple: mean equation = %s  variance equation = %s innovation term = %s",
			dynamo_mean_str[mean], dynamo_var_str[var], dynamo_news_str[news]);
		return goodspec;
	}

	// 4) check model orders and number of explanatory variables
	if( DMS_MP(spec)<0 && DYNAMO_SPEC_MEAN_P_MAX>DMS_MP(spec) ) goodspec=0;
	if( DMS_MP(spec)<0 && DYNAMO_SPEC_MEAN_Q_MAX>DMS_MP(spec) ) goodspec=0;

	if( !goodspec )
	{ 
		sprintf(info_msg,"illegal mean p-q orders (%d,%d)",DMS_MP(spec),DMS_MP(spec)); 
		return goodspec;
	}

	if( DMS_MR(spec)<0 && DYNAMO_SPEC_MEAN_R_MAX>DMS_MR(spec) ) goodspec=0;

	if( !goodspec )
	{ 
		sprintf(info_msg,"illegal number of mean regressors (%d)",DMS_MR(spec)); 
		return goodspec;
	}

	//if( DMS_MX(spec)<0 && DYNAMO_SPEC_MEAN_X_MAX>DMS_MX(spec) ) goodspec=0;

	/*
	if( gsl_vector_int_get(spec,DYNAMO_SPEC_MEAN_X)<0 && DYNAMO_SPEC_MEAN_R_MAX>gsl_vector_int_get(spec,DYNAMO_SPEC_MEAN_X) ) goodspec=0;
	if( gsl_vector_int_get(spec,DYNAMO_SPEC_VAR_P)<0 && DYNAMO_SPEC_VAR_P_MAX>gsl_vector_int_get(spec,DYNAMO_SPEC_VAR_P) ) goodspec=0;
	if( gsl_vector_int_get(spec,DYNAMO_SPEC_VAR_Q)<0 && DYNAMO_SPEC_VAR_Q_MAX>gsl_vector_int_get(spec,DYNAMO_SPEC_VAR_Q) ) goodspec=0;
	if( gsl_vector_int_get(spec,DYNAMO_SPEC_VAR_R)<0 && DYNAMO_SPEC_VAR_R_MAX>gsl_vector_int_get(spec,DYNAMO_SPEC_VAR_R) ) goodspec=0;
	*/

	return goodspec;
}

int dynamo_check_dim(char *info_msg, const gsl_vector_int *spec, const gsl_vector *param, const gsl_vector *y, const gsl_matrix *cond, const gsl_matrix *X, const gsl_matrix *Z)
{
	int mean;
	int gooddim = 1;

	if( param->size != dynamo_spec_ndim(spec) ) 
	{
		sprintf( info_msg,"Illegal param vector size (%d instead of %d).",param->size,dynamo_spec_ndim(spec));

		gooddim=0;
		return gooddim;
	}

	if( cond->size1 != y->size || cond->size2 != dynamo_spec_ncnd(spec) )
	{
		sprintf( info_msg,"Illegal cond matrix size (%d,%d instead of %d %d).",cond->size1,cond->size2,y->size,dynamo_spec_ncnd(spec));

		gooddim=0;
		return gooddim;
	}

	if( dynamo_spec_has_mean_reg(spec) )
	{
		mean = DMS_M(spec);
		if( X->size1 != y->size ) gooddim = 0;

		if( mean==DYNAMO_MEAN_AMEM || mean==DYNAMO_MEAN_ACMEM || mean==DYNAMO_MEAN_ASMEM || mean==DYNAMO_MEAN_ABSMEM )
		{
			if( X->size2 != (DMS_MP(spec)+DMS_MR(spec)+dynamo_spec_ndim_mean_x(spec)) ) gooddim = 0;
		}
		else
		{
			if( X->size2 != (DMS_MR(spec)+dynamo_spec_ndim_mean_x(spec)) ) gooddim = 0;
		}

		if( !gooddim )
		{
			sprintf( info_msg,"Illegal X matrix dimensions for the given mean specification (%d,%d).",X->size1,X->size2);
			return gooddim;
		}
	}
	
	if( dynamo_spec_has_var_reg(spec) )
	{
		if( Z->size1 != y->size ) gooddim = 0;
		if( Z->size2 != (DMS_VR(spec)+DMS_VX(spec)) ) gooddim = 0;

		if( !gooddim )
		{
			sprintf( info_msg,"Illegal Z matrix dimensions for the given variance specification (%d,%d).",Z->size1,Z->size2);
			return gooddim;
		}
	}

	return gooddim;
}

int dynamo_check_data(char *msg, int checktype, const gsl_vector_int *spec, const gsl_vector *param, const gsl_vector *y, const gsl_matrix *X, const gsl_matrix *Z)
{
	int gooddim = 1;
	return gooddim;
}

size_t dynamo_spec_maxlag(const gsl_vector_int *spec)
{
	size_t maxlag = 0;
	
	maxlag = DMS_MP(spec)>maxlag? DMS_MP(spec): maxlag;
	maxlag = DMS_MQ(spec)>maxlag? DMS_MQ(spec): maxlag;
	maxlag = DMS_VP(spec)>maxlag? DMS_VP(spec): maxlag;
	maxlag = DMS_VQ(spec)>maxlag? DMS_VQ(spec): maxlag;

	return maxlag;
}

size_t dynamo_spec_ncnd(const gsl_vector_int *spec)
{
	size_t ncnd = 2;

	switch( DMS_M(spec) )
	{
		case DYNAMO_MEAN_CMEM:
		case DYNAMO_MEAN_SMEM:
		case DYNAMO_MEAN_BSMEM:
		case DYNAMO_MEAN_ACMEM:
		case DYNAMO_MEAN_ASMEM:
		case DYNAMO_MEAN_ABSMEM:
			ncnd += 2;
		break;
	}

	return ncnd;
}

size_t dynamo_spec_ndim(const gsl_vector_int *spec)
{
	return dynamo_spec_ndim_mean(spec) + dynamo_spec_ndim_var(spec) + dynamo_spec_ndim_news(spec);
}

size_t dynamo_spec_ndim_mean(const gsl_vector_int *spec)
{
	return dynamo_spec_has_mean_const(spec)+dynamo_spec_mean_pmul(spec)*DMS_MP(spec)+DMS_MQ(spec)+DMS_MR(spec)+dynamo_spec_ndim_mean_x(spec);
}

size_t dynamo_spec_ndim_var(const gsl_vector_int *spec)
{
	return dynamo_spec_has_var_const(spec)+dynamo_spec_var_pmul(spec)*DMS_VP(spec)+DMS_VQ(spec)+DMS_VR(spec)+DMS_VX(spec);
}

size_t dynamo_spec_ndim_news(const gsl_vector_int *spec)
{
	size_t ndim=0;

	switch( DMS_N(spec) )
	{
		case DYNAMO_NEWS_ST: 
		case DYNAMO_NEWS_GAMMA:
		case DYNAMO_NEWS_WEIBULL:
			ndim+=1; 
			break; 

		case DYNAMO_NEWS_SSKT:
		case DYNAMO_NEWS_GED: 
			ndim+=2; 
			break; 
	}

	return ndim;
}

int dynamo_spec_has_mean_const(const gsl_vector_int *spec)
{
	return ( DMS_M(spec)==DYNAMO_MEAN_CARMA || DMS_M(spec)==DYNAMO_MEAN_MEM || DMS_M(spec)==DYNAMO_MEAN_AMEM );
}

int dynamo_spec_has_var_const(const gsl_vector_int *spec)
{
	return (DMS_V(spec)>DYNAMO_VAR_NONE)*1;
}

int dynamo_spec_has_mean_reg(const gsl_vector_int *spec)
{
	return (DMS_MR(spec)+DMS_MX(spec))>0;
}

int dynamo_spec_has_var_reg(const gsl_vector_int *spec)
{
	return DMS_VR(spec)>0;
}

size_t dynamo_spec_mean_pmul(const gsl_vector_int *spec)
{
	if( DMS_M(spec) == DYNAMO_MEAN_AMEM || DMS_M(spec) == DYNAMO_MEAN_ACMEM ) return 1;
	else return 1;
}

size_t dynamo_spec_ndim_mean_x(const gsl_vector_int *spec)
{
	size_t ndim;

	switch( DMS_M(spec) )
	{
		case DYNAMO_MEAN_SMEM: 
		case DYNAMO_MEAN_BSMEM:
			//ndim = dynamo_spline_ndim(spec);
			ndim = 1;
			break; 
		case DYNAMO_MEAN_CMEM:
			ndim = DMS_MX(spec);
			break;
		default: 
			ndim = 0;
			break;
	}

	return ndim;
}

size_t dynamo_spec_var_pmul(const gsl_vector_int *spec)
{
	if( DMS_V(spec) == DYNAMO_VAR_GJRGARCH ) return 2;
	else return 1;
}

size_t dynamo_spec_mean_const(const gsl_vector_int *spec)
{
	return 0;
}
size_t dynamo_spec_mean_plag(const gsl_vector_int *spec, size_t i)
{
	return dynamo_spec_has_mean_const(spec) + i;
}
size_t dynamo_spec_mean_qlag(const gsl_vector_int *spec, size_t j)
{
	return dynamo_spec_has_mean_const(spec) + dynamo_spec_mean_pmul(spec)*DMS_MP(spec) + j;
}

size_t dynamo_spec_mean_rreg(const gsl_vector_int *spec, size_t k)
{
	return dynamo_spec_has_mean_const(spec) + dynamo_spec_mean_pmul(spec)*DMS_MP(spec) + DMS_MQ(spec) + k;
}

size_t dynamo_spec_mean_xreg(const gsl_vector_int *spec, size_t l)
{
	return dynamo_spec_has_mean_const(spec) + dynamo_spec_mean_pmul(spec)*DMS_MP(spec) + DMS_MQ(spec) + DMS_MR(spec) + l;
}

size_t dynamo_spec_var_const(const gsl_vector_int *spec)
{
	return dynamo_spec_ndim_mean(spec);
}

size_t dynamo_spec_var_plag(const gsl_vector_int *spec, size_t i)
{
	return dynamo_spec_ndim_mean(spec) + dynamo_spec_has_var_const(spec) + i;
}

size_t dynamo_spec_var_qlag(const gsl_vector_int *spec, size_t j)
{
	return dynamo_spec_ndim_mean(spec) + dynamo_spec_has_var_const(spec) + dynamo_spec_var_pmul(spec)*DMS_VP(spec) + j;
}

size_t dynamo_spec_var_rreg(const gsl_vector_int *spec, size_t k)
{
	return dynamo_spec_ndim_mean(spec) + dynamo_spec_has_var_const(spec) + dynamo_spec_var_pmul(spec)*DMS_VP(spec) + DMS_VQ(spec) + k;
}

size_t dynamo_spec_var_xreg(const gsl_vector_int *spec, size_t l)
{
	return dynamo_spec_ndim_mean(spec) + dynamo_spec_has_var_const(spec) + dynamo_spec_var_pmul(spec)*DMS_VP(spec) + DMS_VQ(spec) + DMS_VR(spec) + l;
}

size_t dynamo_spec_news(const gsl_vector_int *spec, size_t i)
{
	return dynamo_spec_ndim_mean(spec) + dynamo_spec_ndim_var(spec) + i;
}

void dynamo_uncond( gsl_vector *uncond, const gsl_vector_int *spec, const gsl_vector *param, const gsl_matrix *X, const gsl_matrix *Z)
{
	switch( DMS_M(spec) )
	{
		case DYNAMO_MEAN_NONE:
		case DYNAMO_MEAN_CONST:
		break;

		case DYNAMO_MEAN_MEM:
			gsl_vector_set( uncond , DMC_MUI(spec) , DMP_MC(spec,param)/(1.0-dynamo_cond_mean_pqsum(spec,param)) );
		break;
	}

	switch( DMS_V(spec) )
	{
		case DYNAMO_VAR_NONE:
		case DYNAMO_VAR_CONST:
		break;

		case DYNAMO_VAR_GARCH:
			gsl_vector_set( uncond , DMC_SIG2I(spec), DMP_VC(spec,param)/(1.0-dynamo_cond_var_pqsum(spec,param)) );
		break;
	}
}

void dynamo_cond_mean_none(gsl_vector *cond_t, const gsl_vector_int *spec, const gsl_vector *param, const gsl_vector *y_l, const gsl_matrix *cond_l, const gsl_vector *x_l, int *flag)
{
}

void dynamo_cond_mean_const(gsl_vector *cond_t, const gsl_vector_int *spec, const gsl_vector *param, const gsl_vector *y_l, const gsl_matrix *cond_l, const gsl_vector *x_l, int *flag)
{
	gsl_vector_set( cond_t , DMC_MUI(spec) , DMP_MC(spec,param) );
}

void dynamo_cond_mean_arma(gsl_vector *cond_t, const gsl_vector_int *spec, const gsl_vector *param, const gsl_vector *y_l, const gsl_matrix *cond_l, const gsl_vector *x_l, int *flag)
{
	size_t i;
	double mu;

	if( dynamo_spec_has_mean_const(spec) )	mu = DMP_MC( spec, param);
	else mu = 0.0;

	for(i=0; i<DMS_MP(spec); ++i) 
	{
		mu += DMP_MP(spec,param,i) * gsl_vector_get( y_l, LAGV(y_l,i));
	}
	for(i=0; i<DMS_MQ(spec); ++i) 
	{
		mu += DMP_MQ(spec,param,i) * gsl_matrix_get( cond_l, LAGM(cond_l,i), DMC_MUI(spec));
	}
	for(i=0; i<DMS_MR(spec); ++i)
	{
		mu += DMP_MR(spec,param,i) * gsl_vector_get( x_l, i); 
	}

	*flag = gsl_finite(mu);
	if( *flag ) gsl_vector_set( cond_t, DMC_MUI(spec), mu);
}

void dynamo_cond_mean_mem(gsl_vector *cond_t, const gsl_vector_int *spec, const gsl_vector *param, const gsl_vector *y_l, const gsl_matrix *cond_l, const gsl_vector *x_l, int *flag)
{
	size_t i;
	double mu;

	mu = DMP_MC( spec, param);
	//printf("%f >",mu);
	for(i=0; i<DMS_MP(spec); ++i) 
	{
		mu += DMP_MP(spec,param,i) * gsl_vector_get( y_l, LAGV(y_l,i));
	}
	//printf("%f >",mu);
	for(i=0; i<DMS_MQ(spec); ++i) 
	{
		mu += DMP_MQ(spec,param,i) * gsl_matrix_get( cond_l, LAGM(cond_l,i), DMC_MUI(spec));
	}
	//printf("%f >",mu);
	for(i=0; i<DMS_MR(spec); ++i)
	{
		mu += DMP_MR(spec,param,i) * gsl_vector_get( x_l, i); 
	}
	//printf("\n");

	*flag = (mu>0.0 && gsl_finite(mu));
	if( *flag ) gsl_vector_set( cond_t, DMC_MUI(spec), mu);
}

void dynamo_cond_mean_amem(gsl_vector *cond_t, const gsl_vector_int *spec, const gsl_vector *param, const gsl_vector *y_l, const gsl_matrix *cond_l, const gsl_vector *x_l, int *flag)
{
	size_t i;
	double mu;

	mu = DMP_MC( spec, param);
	for(i=0; i<DMS_MP(spec); ++i) 
	{
		mu += DMP_MP(spec,param,2*i) * gsl_vector_get( y_l, LAGV(y_l,i));
		mu += DMP_MP(spec,param,2*i+1) * gsl_vector_get( y_l, LAGV(y_l,i))*gsl_vector_get( x_l, LAGV(x_l,i)); 
	}
	for(i=0; i<DMS_MQ(spec); ++i) 
	{
		mu += DMP_MQ(spec,param,i) * gsl_matrix_get( cond_l, LAGM(cond_l,i), DMC_MUI(spec));
	}
	for(i=0; i<DMS_MR(spec); ++i)
	{
		mu += DMP_MR(spec,param,i) * gsl_vector_get( x_l, DMS_MP(spec) + i); 
	}

	*flag = (mu>0.0 && gsl_finite(mu));
	if( *flag ) gsl_vector_set( cond_t, DMC_MUI(spec), mu);
}

void dynamo_cond_mean_cmem(gsl_vector *cond_t, const gsl_vector_int *spec, const gsl_vector *param, const gsl_vector *y_l, const gsl_matrix *cond_l, const gsl_vector *x_l, int *flag)
{
	size_t i;
	double tau, g;

	tau=0;
	for(i=0; i<DMS_MX(spec); ++i)
	{
		tau += DMP_MX(spec,param,i) * gsl_vector_get( x_l, DMS_MR(spec)+i);
	}
	tau = exp( tau );

	g = 1.0-dynamo_cond_mean_pqsum(spec,param);
	for(i=0; i<DMS_MP(spec); ++i) 
	{
		g += DMP_MP(spec,param,2*i) * ( gsl_vector_get( y_l, LAGV(y_l,i)) / gsl_matrix_get( cond_l, LAGM(cond_l,i), DMC_MU_TAUI(spec)) );
	}
	for(i=0; i<DMS_MQ(spec); ++i) 
	{
		g += DMP_MQ(spec,param,i) * gsl_matrix_get( cond_l, LAGM(cond_l,i), DMC_MU_GI(spec) );
	}
	for(i=0; i<DMS_MR(spec); ++i)
	{
		g += DMP_MR(spec,param,i) * gsl_vector_get( x_l, i); 
	}

	*flag = (tau*g>0.0 && gsl_finite(tau*g));
	if( *flag )
	{
		gsl_vector_set( cond_t, DMC_MUI(spec), tau*g);
		gsl_vector_set( cond_t, DMC_MU_TAUI(spec), tau);
		gsl_vector_set( cond_t, DMC_MU_GI(spec), g);
	}
}

void dynamo_cond_var_none(gsl_vector *cond_t, const gsl_vector_int *spec, const gsl_vector *param, const gsl_vector *y_l, const gsl_matrix *cond_l, const gsl_vector *z_l, int *flag)
{
	gsl_vector_set( cond_t, DMC_SIG2I(spec), 1.0);
}

void dynamo_cond_var_const(gsl_vector *cond_t, const gsl_vector_int *spec, const gsl_vector *param, const gsl_vector *y_l, const gsl_matrix *cond_l, const gsl_vector *z_l, int *flag)
{
	gsl_vector_set( cond_t , DMC_SIG2I(spec) , DMP_VC(spec,param) );
}

void dynamo_cond_var_garch(gsl_vector *cond_t, const gsl_vector_int *spec, const gsl_vector *param, const gsl_vector *y_l, const gsl_matrix *cond_l, const gsl_vector *z_l, int *flag)
{
	size_t i;
	double sig2;

	sig2 = DMP_VC( spec, param);
	for( i=0; i<DMS_VP(spec); ++i) 
	{
		sig2 += DMP_VP(spec,param,i) * gsl_pow_2( gsl_vector_get( y_l, LAGV(y_l,i)) - gsl_matrix_get( cond_l, LAGM(cond_l,i), DMC_MUI(spec)));
	}
	for( i=0; i<DMS_VQ(spec); ++i) 
	{
		sig2 += DMP_VQ(spec,param,i) * gsl_matrix_get( cond_l, LAGM(cond_l,i), DMC_SIG2I(spec));
	}
	for( i=0; i<DMS_VR(spec); ++i)
	{
		sig2 += DMP_VR(spec,param,i) * gsl_vector_get( z_l, i); 
	}

	*flag = (sig2>0.0 && gsl_finite(sig2));
	if( *flag ) gsl_vector_set( cond_t, DMC_SIG2I(spec), sig2);
}

void dynamo_cond_var_gjrgarch(gsl_vector *cond_t, const gsl_vector_int *spec, const gsl_vector *param, const gsl_vector *y_l, const gsl_matrix *cond_l, const gsl_vector *z_l, int *flag)
{
	size_t i;
	double sig2;

	sig2 = DMP_VC( spec, param);
	for( i=0; i<DMS_VP(spec); ++i) 
	{
		double negind = ((gsl_vector_get( y_l, LAGV(y_l,i)) - gsl_matrix_get( cond_l, LAGM(cond_l,i), DMC_MUI(spec)))<0);
		sig2 += DMP_VP(spec,param,2*i) * gsl_pow_2( gsl_vector_get( y_l, LAGV(y_l,i)) - gsl_matrix_get( cond_l, LAGM(cond_l,i), DMC_MUI(spec)));
		sig2 += DMP_VP(spec,param,2*i+1) * gsl_pow_2( gsl_vector_get( y_l, LAGV(y_l,i)) - gsl_matrix_get( cond_l, LAGM(cond_l,i), DMC_MUI(spec)))*negind;
	}
	for( i=0; i<DMS_VQ(spec); ++i) 
	{
		sig2 += DMP_VQ(spec,param,i) * gsl_matrix_get( cond_l, LAGM(cond_l,i), DMC_SIG2I(spec));
	}
	for( i=0; i<DMS_VR(spec); ++i)
	{
		sig2 += DMP_VR(spec,param,i) * gsl_vector_get( z_l, i); 
	}

	*flag = (sig2>0.0 && gsl_finite(sig2));
	if( *flag  ) gsl_vector_set( cond_t, DMC_SIG2I(spec), sig2);
}

// utils
int dynamo_string2const(char *str)
{
	     if( strcmp(str,"DYNAMO_SPEC_MEAN") == 0 ) return DYNAMO_SPEC_MEAN;
	else if( strcmp(str,"DYNAMO_SPEC_VAR") == 0 ) return DYNAMO_SPEC_VAR;
	else if( strcmp(str,"DYNAMO_SPEC_NEWS") == 0 ) return DYNAMO_SPEC_NEWS;
	else if( strcmp(str,"DYNAMO_SPEC_MEAN_P") == 0 ) return DYNAMO_SPEC_MEAN_P;
	else if( strcmp(str,"DYNAMO_SPEC_MEAN_Q") == 0 ) return DYNAMO_SPEC_MEAN_Q;
	else if( strcmp(str,"DYNAMO_SPEC_MEAN_R") == 0 ) return DYNAMO_SPEC_MEAN_R;
	else if( strcmp(str,"DYNAMO_SPEC_VAR_P") == 0 ) return DYNAMO_SPEC_VAR_P;
	else if( strcmp(str,"DYNAMO_SPEC_VAR_Q") == 0 ) return DYNAMO_SPEC_VAR_Q;
	else if( strcmp(str,"DYNAMO_SPEC_VAR_R") == 0 ) return DYNAMO_SPEC_VAR_R;

	else if( strcmp(str,"DYNAMO_MEAN_NONE") == 0 ) return DYNAMO_MEAN_NONE;
	else if( strcmp(str,"DYNAMO_MEAN_CONST") == 0 ) return DYNAMO_MEAN_CONST;
	else if( strcmp(str,"DYNAMO_MEAN_ARMA") == 0 ) return DYNAMO_MEAN_ARMA;
	else if( strcmp(str,"DYNAMO_MEAN_CARMA") == 0 ) return DYNAMO_MEAN_CARMA;
	else if( strcmp(str,"DYNAMO_MEAN_MEM") == 0 ) return DYNAMO_MEAN_MEM;
	else if( strcmp(str,"DYNAMO_MEAN_CMEM") == 0 ) return DYNAMO_MEAN_CMEM;

	else if( strcmp(str,"DYNAMO_VAR_NONE") == 0 ) return DYNAMO_VAR_NONE;
	else if( strcmp(str,"DYNAMO_VAR_CONST") == 0 ) return DYNAMO_VAR_CONST;
	else if( strcmp(str,"DYNAMO_VAR_GARCH") == 0 ) return DYNAMO_VAR_GARCH;
	else if( strcmp(str,"DYNAMO_VAR_GJRGARCH") == 0 ) return DYNAMO_VAR_GJRGARCH;

	else if( strcmp(str,"DYNAMO_NEWS_NORM") == 0 ) return DYNAMO_NEWS_NORM;
	else if( strcmp(str,"DYNAMO_NEWS_ST") == 0 ) return DYNAMO_NEWS_ST;
	else if( strcmp(str,"DYNAMO_NEWS_SSKT") == 0 ) return DYNAMO_NEWS_SSKT;
	else if( strcmp(str,"DYNAMO_NEWS_EXP") == 0 ) return DYNAMO_NEWS_EXP;
	else if( strcmp(str,"DYNAMO_NEWS_GAMMA") == 0 ) return DYNAMO_NEWS_GAMMA;
	else if( strcmp(str,"DYNAMO_NEWS_WEIBULL") == 0 ) return DYNAMO_NEWS_WEIBULL;

	else if( strcmp(str,"DYNAMO_STARTVAL_DEFAULT") == 0 ) return DYNAMO_STARTVAL_DEFAULT;
	else if( strcmp(str,"DYNAMO_STARTVAL_USERDEFINED") == 0 ) return DYNAMO_STARTVAL_USERDEFINED;

	else if( strcmp(str,"DYNAMO_OBJ_LOGLIK") == 0 ) return DYNAMO_OBJ_LOGLIK;
	else if( strcmp(str,"DYNAMO_OBJ_PLOGLIK_RIDGE") == 0 ) return DYNAMO_OBJ_PLOGLIK_RIDGE;
	else if( strcmp(str,"DYNAMO_OBJ_PLOGLIK_GRIDGE") == 0 ) return DYNAMO_OBJ_PLOGLIK_GRIDGE;

	else if( strcmp(str,"DYNAMO_LOGLEV_NONE") == 0 ) return DYNAMO_LOGLEV_NONE;
	else if( strcmp(str,"DYNAMO_LOGLEV_BASIC") == 0 ) return DYNAMO_LOGLEV_BASIC;
	else if( strcmp(str,"DYNAMO_LOGLEV_DETAIL") == 0 ) return DYNAMO_LOGLEV_DETAIL;

	else if( strcmp(str,"DYNAMO_PRED_STATIC") == 0 ) return DYNAMO_PRED_STATIC;
	else if( strcmp(str,"DYNAMO_PRED_DYNAMIC") == 0 ) return DYNAMO_PRED_DYNAMIC;
	
	else if( strcmp(str,"DYNAMO_IPRED_FIXED") == 0 ) return DYNAMO_IPRED_FIXED;
	else if( strcmp(str,"DYNAMO_IPRED_RECURSIVE") == 0 ) return DYNAMO_IPRED_RECURSIVE;
	else if( strcmp(str,"DYNAMO_IPRED_ROLLING") == 0 ) return DYNAMO_IPRED_ROLLING;

	return -1;
}

double dynamo_cond_mean_pqsum(const gsl_vector_int *spec, const gsl_vector *param)
{
	size_t i;
	double sum=0.0;

	// TODO: fix for asymmetric models
	for(i=0;i<DMS_MP(spec);++i) sum += DMP_MP(spec,param,i);
	for(i=0;i<DMS_MQ(spec);++i) sum += DMP_MQ(spec,param,i);;

	return sum;
}

double dynamo_cond_var_pqsum(const gsl_vector_int *spec, const gsl_vector *param)
{
	size_t i;
	double sum=0.0;
	
	// TODO: fix for asymmetric models
	for(i=0;i<DMS_VP(spec);++i) sum += DMP_VP(spec,param,i);
	for(i=0;i<DMS_VQ(spec);++i) sum += DMP_VQ(spec,param,i);;

	return sum;
}

int dynamo_expl_X_gen( gsl_matrix *X, const gsl_vector_int *spec, const gsl_matrix *Xin)
{
	return 0;
}

int dynamo_expl_Z_gen( gsl_matrix *Z, const gsl_vector_int *spec, const gsl_matrix *Zin)
{
	return 0;
}


// output
dynamo_printf_t dynamo_printf = printf;

void dynamo_printf_set( dynamo_printf_t f)
{
	dynamo_printf = f;
}

