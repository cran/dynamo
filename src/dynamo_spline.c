
#include "dynamo_base.h"
#include "dynamo.h"

size_t dynamo_spline_ndim(const gsl_vector_int *spec)
{
	size_t ndim = 0;

	switch( DMS_M(spec) )
	{
		case DYNAMO_MEAN_SMEM: ndim = 4+DMS_MX(spec); break;
		case DYNAMO_MEAN_BSMEM: ndim = 3+DMS_MX(spec); break;
	} 

	return ndim;
}

int dynamo_spline_dim_check(const gsl_vector_int *spec, gsl_matrix *X, const gsl_vector *x, const gsl_vector *knots)
{
	int goodspline = 1;
	return goodspline;
}

int dynamo_spline_data_check(const gsl_vector_int *spec, gsl_matrix *X, const gsl_vector *x, const gsl_vector *knots)
{
	int goodspline = 1;
	return goodspline;
}

void dynamo_spline_knots_set(gsl_vector *knots, const gsl_vector_int *spec, double xl, double xr)
{
	size_t i;
	double dx;

	dx = (xr-xl)/( (double) (DMS_MX(spec)+1));

	switch( DMS_M(spec) ) 
	{
		case DYNAMO_MEAN_SMEM:
		{	
			for(i=0; i<4; ++i) gsl_vector_set( knots, i, 0);
			for(i=4; i<DMS_MX(spec); ++i) gsl_vector_set( knots, 4+i, xl+dx+(i-4)*dx);
		}
		break;

		case DYNAMO_MEAN_BSMEM:
		{
			for(i=0; i<DMS_MX(spec); ++i) gsl_vector_set( knots, 4+i, xl+dx+(i-4)*dx);
		}
		break;
	}
}

void dynamo_spline_X_set(gsl_matrix *X, const gsl_vector_int *spec, const gsl_vector *x, const gsl_vector *knots)
{
	size_t i, t;

	switch( DMS_M(spec) ) 
	{
		case DYNAMO_MEAN_SMEM:
		{	

			for(t=0;t<X->size1;t++)	
			{
				gsl_matrix_set( X, t, 0, 1.0);
				gsl_matrix_set( X, t, 1, (double) gsl_vector_get(x,t));
				gsl_matrix_set( X, t, 2, (double) gsl_pow_2(gsl_vector_get(x,t)));
				gsl_matrix_set( X, t, 3, (double) gsl_pow_3(gsl_vector_get(x,t)));

				for(i=4;i<knots->size;++i)
				{
					gsl_matrix_set( X, t, i, 
						gsl_pow_3( (gsl_vector_get(x,t)-gsl_vector_get(knots,i))*((gsl_vector_get(x,t)-gsl_vector_get(knots,i))>0) ) );
				}
			}

		}
		break;

		case DYNAMO_MEAN_BSMEM:
		{
		/*
			size_t l;
			xl = (double) -1;
			xr = (double) (T->size1);
			dx = (xr-xl)/((double) gsl_vector_int_get(spec,CHDM_SPEC_TREND_NKNOTS));

			x = gsl_vector_alloc(T->size1);
			k = gsl_vector_alloc(T->size2);
//			X = gsl_matrix_alloc(T->size1,T->size2);
			K = gsl_matrix_alloc(T->size1,T->size2);
			P = gsl_matrix_alloc(T->size1,T->size2);
			T_up = gsl_matrix_alloc(T->size1,T->size2);

			for(t=0;t<T->size1;++t) gsl_vector_set( x , t , (double) t );
			for(i=0;i<T->size2;++i) gsl_matrix_set_col( X , i , x );

			for(i=0;i<T->size2;++i) gsl_vector_set( knots , i , xl+dx*( i - 3.0 ) );
			for(t=0;t<T->size1;++t) gsl_matrix_set_row( K , t , knots );
			
			for(t=0;t<T->size1;++t)
			{
				for(i=0;i<T->size2;++i)
				{
					gsl_matrix_set( P , t , i , (gsl_matrix_get(X,t,i)-gsl_matrix_get(K,t,i))/dx ); 
					gsl_matrix_set( T , t , i , (gsl_matrix_get(K,t,i)<=gsl_matrix_get(X,t,i)) && (gsl_matrix_get(X,t,i)<(gsl_matrix_get(K,t,i)+dx)) );
				}
			}

			//printf("xl xr dx %f %f %f\n",xl,xr,dx);
			//printf("t\n");
			//gsl_vector_fprintf( stdout , k , "%f" );

			for(dg=1;dg<=3;++dg)
			{
				for(t=0;t<T->size1;++t)
				{
					for(i=0;i<T->size2;++i)
					{
						gsl_matrix_set( T_up , t , i ,
						gsl_matrix_get( P , t , i )*gsl_matrix_get( T , t , i ) 
						+ ( (double) dg + 1.0 - gsl_matrix_get( P , t , i ) )*gsl_matrix_get( T , t , (i+1)<T->size2?(i+1):0 ) );
						gsl_matrix_set( T_up , t , i , gsl_matrix_get( T_up , t , i )/((double) dg) );
					}
				}

				gsl_matrix_memcpy( T , T_up );
			}
	
			//gsl_vector_free( x );
			//gsl_vector_free( k );
			//gsl_matrix_free( X );
			gsl_matrix_free( P );
			gsl_matrix_free( K );
			gsl_matrix_free( T_up );
		*/
		}
		break;

	}
}

