/* horizontal reduce by a float factor with a kernel
 *
 * 29/1/16
 * 	- from shrinkh.c
 * 10/3/16
 * 	- add other kernels
 * 15/8/16
 * 	- rename xshrink as hshrink for consistency
 * 9/9/16
 * 	- add @centre option
 */

/*

    This file is part of VIPS.
    
    VIPS is free software; you can redistribute it and/or modify
    it under the terms of the GNU Lesser General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU Lesser General Public License for more details.

    You should have received a copy of the GNU Lesser General Public License
    along with this program; if not, write to the Free Software
    Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
    02110-1301  USA

 */

/*

    These files are distributed with VIPS - http://www.vips.ecs.soton.ac.uk

 */

/*
#define DEBUG
 */

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif /*HAVE_CONFIG_H*/
#include <vips/intl.h>

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <vips/vips.h>
#include <vips/debug.h>
#include <vips/internal.h>

#include "presample.h"
#include "templates.h"

typedef struct _VipsReduceh {
	VipsResample parent_instance;

	double hshrink;		/* Reduce factor */

	/* The thing we use to make the kernel.
	 */
	VipsKernel kernel;

	/* Use centre rather than corner sampling convention.
	 */
	gboolean centre;

	/* Number of points in kernel.
	 */
	int n_point;

	/* Precalculated interpolation matrices. int (used for pel
	 * sizes up to short), and double (for all others). We go to
	 * scale + 1 so we can round-to-nearest safely.
	 */
	int *matrixi[VIPS_TRANSFORM_SCALE + 1];
	double *matrixf[VIPS_TRANSFORM_SCALE + 1];

} VipsReduceh;

typedef VipsResampleClass VipsReducehClass;

/* We need C linkage for this.
 */
extern "C" {
G_DEFINE_TYPE( VipsReduceh, vips_reduceh, VIPS_TYPE_RESAMPLE );
}

/* Get n points. @shrink is the shrink factor, so 2 for a 50% reduction. 
 */
int
vips_reduce_get_points( VipsKernel kernel, double shrink )
{
	switch( kernel ) {
	case VIPS_KERNEL_NEAREST:
		return( 1 );

	case VIPS_KERNEL_LINEAR:
		return( rint( 2 * shrink ) + 1 );

	case VIPS_KERNEL_CUBIC:
	case VIPS_KERNEL_MITCHELL:
		return( rint( 4 * shrink ) + 1 );

	case VIPS_KERNEL_LANCZOS2:
		/* Needs to be in sync with calculate_coefficients_lanczos().
		 */
		return( rint( 2 * 2 * shrink ) + 1 );

	case VIPS_KERNEL_LANCZOS3:
	case VIPS_KERNEL_APPROX_LANCZOS3:
		return( rint( 2 * 3 * shrink ) + 1 );

	default:
		g_assert_not_reached();
		return( 0 );
	}
}

/* Calculate a mask element. 
 */
void
vips_reduce_make_mask( double *c, VipsKernel kernel, double shrink, double x )
{
	switch( kernel ) {
	case VIPS_KERNEL_NEAREST:
		c[0] = 1.0;
		break;

	case VIPS_KERNEL_LINEAR:
		calculate_coefficients_triangle( c, shrink, x );
		break;

	case VIPS_KERNEL_CUBIC:
		/* Catmull-Rom.
		 */
		calculate_coefficients_cubic( c, shrink, x, 0.0, 0.5 );
		break;

	case VIPS_KERNEL_MITCHELL:
		calculate_coefficients_cubic( c, shrink, x,
			1.0 / 3.0, 1.0 / 3.0 );
		break;

	case VIPS_KERNEL_LANCZOS2:
		calculate_coefficients_lanczos( c, 2, shrink, x );
		break;

	case VIPS_KERNEL_LANCZOS3:
		calculate_coefficients_lanczos( c, 3, shrink, x );
		break;
	case VIPS_KERNEL_APPROX_LANCZOS3:
		calculate_coefficients_approx_lanczos( c, 3, shrink, x );
		break;
	default:
		g_assert_not_reached();
		break;
	}
}

#ifndef restrict
#define restrict __restrict__
#endif

template <typename T, int max_value>
static void inline
reduceh_unsigned_int_tab( VipsReduceh *reduceh,
	VipsPel *pout, const VipsPel *pin,
	const int bands, const int * restrict cx )
{
	T* restrict out = (T *) pout;
	const T* restrict in = (T *) pin;
	const int n = reduceh->n_point;

	for( int z = 0; z < bands; z++ ) {
		int sum;

		sum = reduce_sum<T, int>( in + z, bands, cx, n );
		sum = unsigned_fixed_round( sum );
		sum = VIPS_CLIP( 0, sum, max_value );

		out[z] = sum;
	}
}

template <typename T, int min_value, int max_value>
static void inline
reduceh_signed_int_tab( VipsReduceh *reduceh,
	VipsPel *pout, const VipsPel *pin,
	const int bands, const int * restrict cx )
{
	T* restrict out = (T *) pout;
	const T* restrict in = (T *) pin;
	const int n = reduceh->n_point;

	for( int z = 0; z < bands; z++ ) {
		int sum;

		sum = reduce_sum<T, int>( in, bands, cx, n );
		sum = signed_fixed_round( sum );
		sum = VIPS_CLIP( min_value, sum, max_value );

		out[z] = sum;

		in += 1;
	}
}

/* Floating-point version.
 */
template <typename T>
static void inline
reduceh_float_tab( VipsReduceh *reduceh,
	VipsPel *pout, const VipsPel *pin,
	const int bands, const double *cx )
{
	T* restrict out = (T *) pout;
	const T* restrict in = (T *) pin;
	const int n = reduceh->n_point;

	for( int z = 0; z < bands; z++ ) {
		out[z] = reduce_sum<T, double>( in, bands, cx, n );
		in += 1;
	}
}

/* 32-bit int output needs a double intermediate.
 */

template <typename T, int max_value>
static void inline
reduceh_unsigned_int32_tab( VipsReduceh *reduceh,
	VipsPel *pout, const VipsPel *pin,
	const int bands, const double * restrict cx )
{
	T* restrict out = (T *) pout;
	const T* restrict in = (T *) pin;
	const int n = reduceh->n_point;

	for( int z = 0; z < bands; z++ ) {
		double sum;

		sum = reduce_sum<T, double>( in, bands, cx, n );
		out[z] = VIPS_CLIP( 0, sum, max_value );

		in += 1;
	}
}

template <typename T, int min_value, int max_value>
static void inline
reduceh_signed_int32_tab( VipsReduceh *reduceh,
	VipsPel *pout, const VipsPel *pin,
	const int bands, const double * restrict cx )
{
	T* restrict out = (T *) pout;
	const T* restrict in = (T *) pin;
	const int n = reduceh->n_point;

	for( int z = 0; z < bands; z++ ) {
		double sum;

		sum = reduce_sum<T, double>( in, bands, cx, n );
		sum = VIPS_CLIP( min_value, sum, max_value );
		out[z] = sum;

		in += 1;
	}
}


/* Ultra-high-quality version for double images.
 */
template <typename T>
static void inline
reduceh_notab( VipsReduceh *reduceh,
               VipsPel *pout, const VipsPel *pin,
               const int bands, double x )
{
	T* restrict out = (T *) pout;
	const T* restrict in = (T *) pin;
	const int n = reduceh->n_point;

	double cx[MAX_POINT];

	vips_reduce_make_mask( cx, reduceh->kernel, reduceh->hshrink, x );

	for( int z = 0; z < bands; z++ ) {
		out[z] = reduce_sum<T, double>( in, bands, cx, n );

		in += 1;
	}
}

#define EPSILON  (1.0e-12)

static inline double reciprocal(const double x)
{
	double sign;

	/*
	  Return 1/x where x is perceptible (not unlimited or infinitesimal).
	*/
	sign=x < 0.0 ? -1.0 : 1.0;
	if ((sign*x) >= EPSILON)
		return (1.0/x);
	return (sign/EPSILON);
}

template <typename T, int max_value>
static void inline
reduceh_notab_blend( VipsReduceh *reduceh,
               VipsPel *pout, const VipsPel *pin,
               const int bands, double x )
{
	T* restrict out = (T *) pout;
	const T* restrict in = (T *) pin;
	const int n = reduceh->n_point;
	int band;
	const int stride = bands;

	double cx[MAX_POINT];

	vips_reduce_make_mask( cx, reduceh->kernel, reduceh->hshrink, x );

	int alpha_band = bands - 1;

	for( band = 0; band < bands - 1; band++ )
	{
		double sum = 0;
		double normalize = 0;

		for( int i = 0; i < n; i++ )
		{
			double alpha = (double)in[i * stride + alpha_band] * (1.0 / max_value);

			sum += alpha * (double)cx[i] * (double)in[i * stride + band];
			normalize += alpha;
		}
//		normalize = 1;
		out[band] = (T)VIPS_CLIP(0, sum * reciprocal( normalize), max_value);
	}

	double sum = 0;

	for( int i = 0; i < n - 1; i++ ) {
		sum += cx[i] * in[i * stride + alpha_band];
	}

	out[band] = (T)VIPS_CLIP(0, sum, max_value);
}


static int
vips_reduceh_gen( VipsRegion *out_region, void *seq,
                  void *a, void *b, gboolean *stop )
{
	VipsImage *in = (VipsImage *) a;
	VipsReduceh *reduceh = (VipsReduceh *) b;
	VipsRegion *ir = (VipsRegion *) seq;
	VipsRect *r = &out_region->valid;

	/* Double bands for complex.
	 */
	const int bands = in->Bands *
	                  (vips_band_format_iscomplex( in->BandFmt ) ? 2 : 1);

	VipsRect s;

#ifdef DEBUG
	printf( "vips_reduceh_gen: generating %d x %d at %d x %d\n",
		r->width, r->height, r->left, r->top );
#endif /*DEBUG*/
	printf( "vips_reduceh_gen: generating %d x %d at %d x %d\n",
	        r->width, r->height, r->left, r->top );

	s.left = r->left * reduceh->hshrink;
	s.top = r->top;
	s.width = r->width * reduceh->hshrink;
	s.height = r->height;
	if( vips_region_prepare( ir, &s ) )
		return (-1);

	VIPS_GATE_START( "vips_reduceh_gen: work" );
	int resize_filter_support = 3;
	double resize_filter_scale = (1.0/3.0);
	double scale = reduceh->hshrink;
	double support = scale * resize_filter_support;

	typedef unsigned short T;
	const int max_value = USHRT_MAX;
	scale = reciprocal(scale);


	for( int x = 0; x < r->width; x++ ) {
		double bisect = (double) (x + 0.5) / scale + EPSILON;
		size_t start = (ssize_t) VIPS_MAX( bisect - support + 0.5, 0.0 );
		size_t stop = (ssize_t) VIPS_MIN( bisect + support + 0.5, in->Xsize );
		double weight[1000];
		double density = 0;
		int n;

		for( n = 0; n < (stop - start); n++ ) {
			double wx = VIPS_ABS(scale * ((double) (start + n) - bisect + 0.5));
			weight[n] = sinc_fast( wx * resize_filter_scale) * sinc_fast( wx );
			density += weight[n];
		}

		if( n == 0 )
			continue;

		if( (density != 0.0) && (density != 1.0) ) {
			/*
			  Normalize.
			*/
			density = reciprocal( density );
			for( int i = 0; i < n; i++ ) {
				weight[i] *= density;
			}
		}

		const T *p = (const T *) VIPS_REGION_ADDR(
			ir, ir->valid.left + start, r->top );

		for( int y = 0; y < r->height; y++ ) {
			for( int i = 0; i < bands; i++ ) {
				T *q = (T *) VIPS_REGION_ADDR( out_region, r->left + x,
				                               r->top + y );
				double pixel = 0;

				if( i == bands - 1 ) {
					/*
					  No alpha blending.
					*/
					for( int j = 0; j < n; j++ ) {
						int k = y * n + j;
						pixel += weight[j] * p[k * bands + i];
					}
					q[i] = (T) VIPS_CLIP( 0, pixel, max_value );
//					printf( "%f,%f,%f,%f,%d\n",
//					        weight[0],
//					        weight[1],
//					        weight[2],
//					        weight[3],
//					        q[i] );

					continue;
				}

				/*
		          Alpha blending.
		        */
				double gamma = 0.0;
				for( int j = 0; j < n; j++ ) {
					int k = y * n + j;

					T alpha_value = p[k * bands + bands - 1];
					T pixel_value = p[k * bands + i];
					double alpha = (1.0 / max_value) * alpha_value;
					pixel += alpha * weight[j] * pixel_value;
					gamma += alpha;
				}
				gamma = reciprocal( gamma );
				q[i] = VIPS_CLIP( 0, gamma * pixel, max_value );

//				printf( "%f,%f,%f,%f,%d\n",
//				        weight[0],
//				        weight[1],
//				        weight[2],
//				        weight[3],
//				        q[i] );

			}
		}
	}

	VIPS_GATE_STOP( "vips_reduceh_gen: work" );

	VIPS_COUNT_PIXELS( out_region, "vips_reduceh_gen" );

	return (0);
}

static int
vips_reduceh_build( VipsObject *object )
{
	VipsObjectClass *object_class = VIPS_OBJECT_GET_CLASS( object );
	VipsResample *resample = VIPS_RESAMPLE( object );
	VipsReduceh *reduceh = (VipsReduceh *) object;
	VipsImage **t = (VipsImage **)
		vips_object_local_array( object, 2 );

	VipsImage *in;
	int width;

	if( VIPS_OBJECT_CLASS( vips_reduceh_parent_class )->build( object ) )
		return( -1 );

	in = resample->in;

	if( reduceh->hshrink < 1 ) {
		vips_error( object_class->nickname,
			"%s", _( "reduce factors should be >= 1" ) );
		return( -1 );
	}

	if( reduceh->hshrink == 1 )
		return( vips_image_write( in, resample->out ) );

	/* Build the tables of pre-computed coefficients.
	 */
	reduceh->n_point =
		vips_reduce_get_points( reduceh->kernel, reduceh->hshrink );
	g_info( "reduceh: %d point mask", reduceh->n_point );
	if( reduceh->n_point > MAX_POINT ) {
		vips_error( object_class->nickname,
			"%s", _( "reduce factor too large" ) );
		return( -1 );
	}
	for( int x = 0; x < VIPS_TRANSFORM_SCALE + 1; x++ ) {
		reduceh->matrixf[x] =
			VIPS_ARRAY( object, reduceh->n_point, double );
		reduceh->matrixi[x] =
			VIPS_ARRAY( object, reduceh->n_point, int );
		if( !reduceh->matrixf[x] ||
			!reduceh->matrixi[x] )
			return( -1 );

		vips_reduce_make_mask( reduceh->matrixf[x],
			reduceh->kernel, reduceh->hshrink,
			(float) x / VIPS_TRANSFORM_SCALE );

		for( int i = 0; i < reduceh->n_point; i++ )
			reduceh->matrixi[x][i] = reduceh->matrixf[x][i] *
				VIPS_INTERPOLATE_SCALE;

#ifdef DEBUG
		printf( "vips_reduceh_build: mask %d\n    ", x );
		for( int i = 0; i < reduceh->n_point; i++ )
			printf( "%d ", reduceh->matrixi[x][i] );
		printf( "\n" );
#endif /*DEBUG*/
	}

	/* Unpack for processing.
	 */
	if( vips_image_decode( in, &t[0] ) )
		return( -1 );
	in = t[0];

	/* Add new pixels around the input so we can interpolate at the edges.
	 * In centre mode, we read 0.5 pixels more to the right, so we must
	 * enlarge a little further.
	 */
//	width = in->Xsize + reduceh->n_point - 1;
//	if( reduceh->centre )
//		width += 1;
//	if( vips_embed( in, &t[1],
//		reduceh->n_point / 2 - 1, 0,
//		width, in->Ysize,
//		"extend", VIPS_EXTEND_COPY,
//		(void *) NULL ) )
//		return( -1 );
//	in = t[1];

	if( vips_image_pipelinev( resample->out,
		VIPS_DEMAND_STYLE_THINSTRIP, in, (void *) NULL ) )
		return( -1 );

	/* Size output. We need to always round to nearest, so round(), not
	 * rint().
	 *
	 * Don't change xres/yres, leave that to the application layer. For
	 * example, vipsthumbnail knows the true reduce factor (including the
	 * fractional part), we just see the integer part here.
	 */
	resample->out->Xsize = VIPS_ROUND_UINT(
		resample->in->Xsize / reduceh->hshrink );
	if( resample->out->Xsize <= 0 ) {
		vips_error( object_class->nickname,
			"%s", _( "image has shrunk to nothing" ) );
		return( -1 );
	}

#ifdef DEBUG
	printf( "vips_reduceh_build: reducing %d x %d image to %d x %d\n",
		in->Xsize, in->Ysize, 
		resample->out->Xsize, resample->out->Ysize );
#endif /*DEBUG*/
	printf( "vips_reduceh_build: reducing %d x %d image to %d x %d\n",
	        in->Xsize, in->Ysize,
	        resample->out->Xsize, resample->out->Ysize );

	if( vips_image_generate( resample->out,
		vips_start_one, vips_reduceh_gen, vips_stop_one,
		in, reduceh ) )
		return( -1 );

	vips_reorder_margin_hint( resample->out, reduceh->n_point );

	return( 0 );
}

static void
vips_reduceh_class_init( VipsReducehClass *reduceh_class )
{
	GObjectClass *gobject_class = G_OBJECT_CLASS( reduceh_class );
	VipsObjectClass *vobject_class = VIPS_OBJECT_CLASS( reduceh_class );
	VipsOperationClass *operation_class =
		VIPS_OPERATION_CLASS( reduceh_class );

	VIPS_DEBUG_MSG( "vips_reduceh_class_init\n" );

	gobject_class->set_property = vips_object_set_property;
	gobject_class->get_property = vips_object_get_property;

	vobject_class->nickname = "reduceh";
	vobject_class->description = _( "shrink an image horizontally" );
	vobject_class->build = vips_reduceh_build;

	operation_class->flags = VIPS_OPERATION_SEQUENTIAL;

	VIPS_ARG_DOUBLE( reduceh_class, "hshrink", 3,
		_( "Hshrink" ),
		_( "Horizontal shrink factor" ),
		VIPS_ARGUMENT_REQUIRED_INPUT,
		G_STRUCT_OFFSET( VipsReduceh, hshrink ),
		1, 1000000, 1 );

	VIPS_ARG_ENUM( reduceh_class, "kernel", 3,
		_( "Kernel" ),
		_( "Resampling kernel" ),
		VIPS_ARGUMENT_OPTIONAL_INPUT,
		G_STRUCT_OFFSET( VipsReduceh, kernel ),
		VIPS_TYPE_KERNEL, VIPS_KERNEL_LANCZOS3 );

	VIPS_ARG_BOOL( reduceh_class, "centre", 7,
		_( "Centre" ),
		_( "Use centre sampling convention" ),
		VIPS_ARGUMENT_OPTIONAL_INPUT,
		G_STRUCT_OFFSET( VipsReduceh, centre ),
		FALSE );

	/* Old name.
	 */
	VIPS_ARG_DOUBLE( reduceh_class, "xshrink", 3,
		_( "Xshrink" ),
		_( "Horizontal shrink factor" ),
		VIPS_ARGUMENT_REQUIRED_INPUT | VIPS_ARGUMENT_DEPRECATED,
		G_STRUCT_OFFSET( VipsReduceh, hshrink ),
		1, 1000000, 1 );

}

static void
vips_reduceh_init( VipsReduceh *reduceh )
{
	reduceh->kernel = VIPS_KERNEL_LANCZOS3;
}

/* See reduce.c for the doc comment.
 */

int
vips_reduceh( VipsImage *in, VipsImage **out, double hshrink, ... )
{
	va_list ap;
	int result;

	va_start( ap, hshrink );
	result = vips_call_split( "reduceh", ap, in, out, hshrink );
	va_end( ap );

	return( result );
}
