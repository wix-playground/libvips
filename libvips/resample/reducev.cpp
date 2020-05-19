/* vertical reduce by a float factor with a kernel
 *
 * 29/1/16
 * 	- from shrinkv.c
 * 10/3/16
 * 	- add other kernels
 * 21/3/16
 * 	- add vector path
 * 2/4/16
 * 	- better int mask creation ... we now adjust the scale to keep the sum
 * 	  equal to the target scale
 * 15/6/16
 * 	- better accuracy with smarter multiplication
 * 15/8/16
 * 	- rename yshrink as vshrink for consistency
 * 9/9/16
 * 	- add @centre option
 * 7/3/17
 * 	- add a seq line cache
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
#define DEBUG_PIXELS
#define DEBUG_COMPILE
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
#include <vips/vector.h>

#include "presample.h"
#include "templates.h"

/* We can't run more than this many passes. Larger than this and we
 * fall back to C.
 */
#define MAX_PASS (10)

/* The number of params we pass for coeffs. Orc limits this rather. 
 */
#define MAX_PARAM (8)

/* A pass with a vector. 
 */
typedef struct {
	int first;		/* The index of the first mask coff we use */
	int last;		/* The index of the last mask coff we use */

	int r;			/* Set previous result in this var */
	int d2;			/* Write new temp result here */

	int p[MAX_PARAM];	/* Mask coeffs passed in these */
	int n_param;

        /* The code we generate for this section of this mask. 
	 */
        VipsVector *vector;
} Pass;

typedef struct _VipsReducev {
	VipsResample parent_instance;

	double vshrink;		/* Shrink factor */

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

	/* And another set for orc: we want 2.6 precision.
	 */
	int *matrixo[VIPS_TRANSFORM_SCALE + 1];

	/* The passes we generate for this mask.
	 */
	int n_pass;	
	Pass pass[MAX_PASS];

} VipsReducev;

typedef VipsResampleClass VipsReducevClass;

/* We need C linkage for this.
 */
extern "C" {
G_DEFINE_TYPE( VipsReducev, vips_reducev, VIPS_TYPE_RESAMPLE );
}

static void
vips_reducev_finalize( GObject *gobject )
{
	VipsReducev *reducev = (VipsReducev *) gobject; 

	for( int i = 0; i < reducev->n_pass; i++ )
		VIPS_FREEF( vips_vector_free, reducev->pass[i].vector );
	reducev->n_pass = 0;
	for( int i = 0; i < VIPS_TRANSFORM_SCALE + 1; i++ ) {
		VIPS_FREE( reducev->matrixf[i] );
		VIPS_FREE( reducev->matrixi[i] );
		VIPS_FREE( reducev->matrixo[i] );
	}

	G_OBJECT_CLASS( vips_reducev_parent_class )->finalize( gobject );
}

#define TEMP( N, S ) vips_vector_temporary( v, (char *) N, S )
#define PARAM( N, S ) vips_vector_parameter( v, (char *) N, S )
#define SCANLINE( N, P, S ) vips_vector_source_scanline( v, (char *) N, P, S )
#define CONST( N, V, S ) vips_vector_constant( v, (char *) N, V, S )
#define ASM2( OP, A, B ) vips_vector_asm2( v, (char *) OP, A, B )
#define ASM3( OP, A, B, C ) vips_vector_asm3( v, (char *) OP, A, B, C )

/* Generate code for a section of the mask. first is the index we start
 * at, we set last to the index of the last one we use before we run 
 * out of intermediates / constants / parameters / sources or mask
 * coefficients.
 *
 * 0 for success, -1 on error.
 */
static int
vips_reducev_compile_section( VipsReducev *reducev, Pass *pass, gboolean first )
{
	VipsVector *v;
	int i;

#ifdef DEBUG_COMPILE
	printf( "starting pass %d\n", pass->first ); 
#endif /*DEBUG_COMPILE*/

	pass->vector = v = vips_vector_new( "reducev", 1 );

	/* We have two destinations: the final output image (8-bit) and the
	 * intermediate buffer if this is not the final pass (16-bit).
	 */
	pass->d2 = vips_vector_destination( v, "d2", 2 );

	/* "r" is the array of sums from the previous pass (if any).
	 */
	pass->r = vips_vector_source_name( v, "r", 2 );

	/* The value we fetch from the image, the accumulated sum.
	 */
	TEMP( "value", 2 );
	TEMP( "sum", 2 );

	/* Init the sum. If this is the first pass, it's a constant. If this
	 * is a later pass, we have to init the sum from the result 
	 * of the previous pass. 
	 */
	if( first ) {
		char c0[256];

		CONST( c0, 0, 2 );
		ASM2( "loadpw", "sum", c0 );
	}
	else 
		ASM2( "loadw", "sum", "r" );

	for( i = pass->first; i < reducev->n_point; i++ ) {
		char source[256];
		char coeff[256];

		SCANLINE( source, i, 1 );

		/* This mask coefficient.
		 */
		vips_snprintf( coeff, 256, "p%d", i );
		pass->p[pass->n_param] = PARAM( coeff, 2 );
		pass->n_param += 1;
		if( pass->n_param >= MAX_PARAM )
			return( -1 );

		/* Mask coefficients are 2.6 bits fixed point. We need to hold
		 * about -0.5 to 1.0, so -2 to +1.999 is as close as we can
		 * get. 
		 *
		 * We need a signed multiply, so the image pixel needs to
		 * become a signed 16-bit value. We know only the bottom 8 bits
		 * of the image and coefficient are interesting, so we can take
		 * the bottom bits of a 16x16->32 multiply. 
		 *
		 * We accumulate the signed 16-bit result in sum.
		 */
		ASM2( "convubw", "value", source );
		ASM3( "mullw", "value", "value", coeff );
		ASM3( "addssw", "sum", "sum", "value" );

		/* We've used this coeff.
		 */
		pass->last = i;

		if( vips_vector_full( v ) )
			break;

		/* orc 0.4.24 and earlier hate more than about five lines at
		 * once :( 
		 */
		if( i - pass->first > 3 )
			break;
	}

	/* If this is the end of the mask, we write the 8-bit result to the
	 * image, otherwise write the 16-bit intermediate to our temp buffer. 
	 */
	if( pass->last >= reducev->n_point - 1 ) {
		char c32[256];
		char c6[256];
		char c0[256];
		char c255[256];

		CONST( c32, 32, 2 );
		ASM3( "addw", "sum", "sum", c32 );
		CONST( c6, 6, 2 );
		ASM3( "shrsw", "sum", "sum", c6 );

		/* You'd think "convsuswb", convert signed 16-bit to unsigned
		 * 8-bit with saturation, would be quicker, but it's a lot
		 * slower.
		 */
		CONST( c0, 0, 2 );
		ASM3( "maxsw", "sum", c0, "sum" ); 
		CONST( c255, 255, 2 );
		ASM3( "minsw", "sum", c255, "sum" ); 

		ASM2( "convwb", "d1", "sum" );
	}
	else 
		ASM2( "copyw", "d2", "sum" );

	if( !vips_vector_compile( v ) ) 
		return( -1 );

#ifdef DEBUG_COMPILE
	printf( "done coeffs %d to %d\n", pass->first, pass->last );
	vips_vector_print( v );
#endif /*DEBUG_COMPILE*/

	return( 0 );
}

static int
vips_reducev_compile( VipsReducev *reducev )
{
	Pass *pass;

	/* Generate passes until we've used up the whole mask.
	 */
	for( int i = 0;; ) {
		/* Allocate space for another pass.
		 */
		if( reducev->n_pass == MAX_PASS ) 
			return( -1 );
		pass = &reducev->pass[reducev->n_pass];
		reducev->n_pass += 1;

		pass->first = i;
		pass->r = -1;
		pass->d2 = -1;
		pass->n_param = 0;

		if( vips_reducev_compile_section( reducev,
			pass, reducev->n_pass == 1 ) )
			return( -1 );
		i = pass->last + 1;

		if( i >= reducev->n_point )
			break;
	}

	return( 0 );
}

/* Our sequence value.
 */
typedef struct {
	VipsReducev *reducev;
	VipsRegion *ir;		/* Input region */

	/* In vector mode we need a pair of intermediate buffers to keep the 
	 * results of each pass in.
	 */
	signed short *t1;
	signed short *t2;
} Sequence;

static int
vips_reducev_stop( void *vseq, void *a, void *b )
{
	Sequence *seq = (Sequence *) vseq;

	VIPS_UNREF( seq->ir );
	VIPS_FREE( seq->t1 );
	VIPS_FREE( seq->t2 );

	return( 0 );
}

static void *
vips_reducev_start( VipsImage *out, void *a, void *b )
{
	VipsImage *in = (VipsImage *) a;
	VipsReducev *reducev = (VipsReducev *) b;
	int sz = VIPS_IMAGE_N_ELEMENTS( in );

	Sequence *seq;

	if( !(seq = VIPS_NEW( out, Sequence )) )
		return( NULL );

	/* Init!
	 */
	seq->reducev = reducev;
	seq->ir = NULL;
	seq->t1 = NULL;
	seq->t2 = NULL;

	/* Attach region and arrays.
	 */
	seq->ir = vips_region_new( in );
	seq->t1 = VIPS_ARRAY( NULL, sz, signed short );
	seq->t2 = VIPS_ARRAY( NULL, sz, signed short );
	if( !seq->ir || 
		!seq->t1 || 
		!seq->t2  ) {
		vips_reducev_stop( seq, NULL, NULL );
		return( NULL );
	}

	return( seq );
}

static inline double reciprocal(const double x)
{
	const double EPSILON = 1.0e-12;

	double
		sign;

	/*
	  Return 1/x where x is perceptible (not unlimited or infinitesimal).
	*/
	sign=x < 0.0 ? -1.0 : 1.0;
	if ((sign*x) >= EPSILON)
		return(1.0/x);
	return(sign/EPSILON);
}

/* You'd think this would vectorise, but gcc hates mixed types in nested loops
 * :-(
 */
template <typename T, int max_value>
static void inline
reducev_unsigned_int_tab( VipsReducev *reducev, VipsPel *pout,
                          const VipsPel *pin, const int width_times_bands, const int lskip,
                          const int *restrict coefficients, const int bands )
{
	const gboolean has_alpha = TRUE; //TODO: accept as parameter or something
	T* restrict out = (T *) pout;
	const T* restrict in = (T *) pin;
	const int stride = lskip / sizeof( T );

	for( int x = 0; x < width_times_bands; x++ ) {

		int band = x % bands;

		if( !has_alpha || band == bands - 1 ) {
			//No alpha blending
			int sum = 0;
			const T *in_element = in + x;

			for( int y = 0; y < reducev->n_point; y++ ) {
				sum += coefficients[y] * in_element[0];
				in_element += stride;
			}

			sum = signed_fixed_round( sum );
			sum = VIPS_CLIP( 0, sum, max_value );

			out[x] = sum;
			continue;
		}

		//Alpha blending
		const T *in_element = in + x;
		const T *in_alpha_element = in + x - band + bands - 1;

		double sum = 0;
		double gamma = 0;
		for( int y = 0; y < reducev->n_point; y++ ) {
			double alpha = (double)in_alpha_element[0] / max_value;
			double multiplied_alpha = alpha * coefficients[y];

			sum += multiplied_alpha * in_element[0];
			gamma += multiplied_alpha;
			
			in_element += stride;
			in_alpha_element += stride;
		}
		gamma = reciprocal( gamma );

		sum = VIPS_CLIP( 0, gamma * sum, max_value );
		out[x] = sum;
	}
}

template <typename T, int min_value, int max_value>
static void inline
reducev_signed_int_tab( VipsReducev *reducev,
	VipsPel *pout, const VipsPel *pin,
	const int ne, const int lskip, const int * restrict cy )
{
	T* restrict out = (T *) pout;
	const T* restrict in = (T *) pin;
	const int n = reducev->n_point;
	const int l1 = lskip / sizeof( T );

	for( int z = 0; z < ne; z++ ) {
		int sum;

		sum = reduce_sum<T, int>( in + z, l1, cy, n );
		sum = signed_fixed_round( sum ); 
		sum = VIPS_CLIP( min_value, sum, max_value ); 

		out[z] = sum;
	}
}

/* Floating-point version.
 */
template <typename T>
static void inline
reducev_float_tab( VipsReducev *reducev,
	VipsPel *pout, const VipsPel *pin,
	const int ne, const int lskip, const double * restrict cy )
{
	T* restrict out = (T *) pout;
	const T* restrict in = (T *) pin;
	const int n = reducev->n_point;
	const int l1 = lskip / sizeof( T );

	for( int z = 0; z < ne; z++ ) 
		out[z] = reduce_sum<T, double>( in + z, l1, cy, n );
}

/* 32-bit int output needs a double intermediate.
 */

template <typename T, int max_value>
static void inline
reducev_unsigned_int32_tab( VipsReducev *reducev,
	VipsPel *pout, const VipsPel *pin,
	const int ne, const int lskip, const double * restrict cy )
{
	T* restrict out = (T *) pout;
	const T* restrict in = (T *) pin;
	const int n = reducev->n_point;
	const int l1 = lskip / sizeof( T );

	for( int z = 0; z < ne; z++ ) {
		double sum;

		sum = reduce_sum<T, double>( in + z, l1, cy, n );
		out[z] = VIPS_CLIP( 0, sum, max_value ); 
	}
}

template <typename T, int min_value, int max_value>
static void inline
reducev_signed_int32_tab( VipsReducev *reducev,
	VipsPel *pout, const VipsPel *pin,
	const int ne, const int lskip, const double * restrict cy )
{
	T* restrict out = (T *) pout;
	const T* restrict in = (T *) pin;
	const int n = reducev->n_point;
	const int l1 = lskip / sizeof( T );

	for( int z = 0; z < ne; z++ ) {
		double sum;

		sum = reduce_sum<T, double>( in + z, l1, cy, n );
		out[z] = VIPS_CLIP( min_value, sum, max_value ); 
	}
}

/* Ultra-high-quality version for double images.
 */
template <typename T, int max_value>
static void inline
reducev_notab( VipsReducev *reducev,
	VipsPel *pout, const VipsPel *pin,
	const int width_times_bands, const int lskip, double y, int bands )
{
	double coefficients[MAX_POINT];
	vips_reduce_make_mask( coefficients, reducev->kernel, reducev->vshrink, y );

	const gboolean has_alpha = TRUE; //TODO: accept as parameter or something
	T* restrict out = (T *) pout;
	const T* restrict in = (T *) pin;
	const int stride = lskip / sizeof( T );

	for( int x = 0; x < width_times_bands; x++ ) {

		int band = x % bands;

		if( !has_alpha || band == bands - 1 ) {
			//No alpha blending
			double sum = 0;
			const T *in_element = in + x;

			for( int i = 0; i < reducev->n_point; i++ ) {
				sum += coefficients[i] * in_element[0];
				in_element += stride;
			}

			out[x] = sum;
			continue;
		}

		//Alpha blending
		const T *in_element = in + x;
		const T *in_alpha_element = in + x - band + bands - 1;

		double sum = 0;
		double gamma = 0;
		for( int i = 0; i < reducev->n_point; i++ ) {
			double alpha = (double)in_alpha_element[0] / max_value;
			double multiplied_alpha = alpha * coefficients[i];

			sum += multiplied_alpha * in_element[0];
			gamma += multiplied_alpha;

			in_element += stride;
			in_alpha_element += stride;
		}
		gamma = reciprocal( gamma );

		sum = VIPS_CLIP( 0, gamma * sum, max_value );
		out[x] = sum;
	}
}

static int
vips_reducev_gen( VipsRegion *out_region, void *vseq, 
	void *a, void *b, gboolean *stop )
{
	VipsImage *in = (VipsImage *) a;
	VipsReducev *reducev = (VipsReducev *) b;
	Sequence *seq = (Sequence *) vseq;
	VipsRegion *ir = seq->ir;
	VipsRect *r = &out_region->valid;

	/* Double bands for complex.
	 */
	const int bands = in->Bands * 
		(vips_band_format_iscomplex( in->BandFmt ) ?  2 : 1);
	int num_elements = r->width * bands;

	VipsRect s;

#ifdef DEBUG
	printf( "vips_reducev_gen: generating %d x %d at %d x %d\n",
		r->width, r->height, r->left, r->top ); 
#endif /*DEBUG*/

	s.left = r->left;
	s.top = r->top * reducev->vshrink;
	s.width = r->width;
	s.height = r->height * reducev->vshrink + reducev->n_point;
	if( reducev->centre )
		s.height += 1;
	if( vips_region_prepare( ir, &s ) )
		return( -1 );

	VIPS_GATE_START( "vips_reducev_gen: work" ); 

	for( int y = 0; y < r->height; y ++ ) { 
		VipsPel *q = 
			VIPS_REGION_ADDR( out_region, r->left, r->top + y );
		const double Y = (r->top + y) * reducev->vshrink + 
			(reducev->centre ? 0.5 : 0.0);
		VipsPel *p = VIPS_REGION_ADDR( ir, r->left, (int) Y ); 
		const int scaled_y = Y * VIPS_TRANSFORM_SCALE * 2;
		const int siy = scaled_y & (VIPS_TRANSFORM_SCALE * 2 - 1);
		const int coefficients_index = (siy + 1) >> 1;
		const int *coefficients_y_int = reducev->matrixi[coefficients_index];
		const double *coefficients_y_float = reducev->matrixf[coefficients_index];
		const int lskip = VIPS_REGION_LSKIP( ir );

		switch( in->BandFmt ) {
		case VIPS_FORMAT_UCHAR:
			reducev_unsigned_int_tab
				<unsigned char, UCHAR_MAX>(
				reducev,
				q, p, num_elements, lskip, coefficients_y_int, bands );
			break;

		case VIPS_FORMAT_CHAR:
			reducev_signed_int_tab
				<signed char, SCHAR_MIN, SCHAR_MAX>(
				reducev,
				q, p, num_elements, lskip, coefficients_y_int );
			break;

		case VIPS_FORMAT_USHORT:
			reducev_unsigned_int_tab
				<unsigned short, USHRT_MAX>(
				reducev,
				q, p, num_elements, lskip, coefficients_y_int, bands);
			break;

		case VIPS_FORMAT_SHORT:
			reducev_signed_int_tab
				<signed short, SHRT_MIN, SHRT_MAX>(
				reducev,
				q, p, num_elements, lskip, coefficients_y_int );
			break;

		case VIPS_FORMAT_UINT:
			reducev_unsigned_int32_tab
				<unsigned int, INT_MAX>(
				reducev,
				q, p, num_elements, lskip, coefficients_y_float );
			break;

		case VIPS_FORMAT_INT:
			reducev_signed_int32_tab
				<signed int, INT_MIN, INT_MAX>(
				reducev,
				q, p, num_elements, lskip, coefficients_y_float );
			break;

		case VIPS_FORMAT_FLOAT:
		case VIPS_FORMAT_COMPLEX:
			reducev_float_tab<float>( reducev,
			                          q, p, num_elements, lskip, coefficients_y_float );
			break;

		case VIPS_FORMAT_DPCOMPLEX:
		case VIPS_FORMAT_DOUBLE:
			reducev_notab<double, USHRT_MAX>( reducev,
			                       q, p, num_elements, lskip, Y - (int) Y, bands );
			break;

		default:
			g_assert_not_reached();
			break;
		}
	}

	VIPS_GATE_STOP( "vips_reducev_gen: work" ); 

	VIPS_COUNT_PIXELS( out_region, "vips_reducev_gen" ); 

	return( 0 );
}

/* Process uchar images with a vector path.
 */
static int
vips_reducev_vector_gen( VipsRegion *out_region, void *vseq, 
	void *a, void *b, gboolean *stop )
{
	VipsImage *in = (VipsImage *) a;
	VipsReducev *reducev = (VipsReducev *) b;
	Sequence *seq = (Sequence *) vseq;
	VipsRegion *ir = seq->ir;
	VipsRect *r = &out_region->valid;
	int ne = r->width * in->Bands;

	VipsExecutor executor[MAX_PASS];
	VipsRect s;

#ifdef DEBUG_PIXELS
	printf( "vips_reducev_vector_gen: generating %d x %d at %d x %d\n",
		r->width, r->height, r->left, r->top ); 
#endif /*DEBUG_PIXELS*/

	s.left = r->left;
	s.top = r->top * reducev->vshrink;
	s.width = r->width;
	s.height = r->height * reducev->vshrink + reducev->n_point;
	if( reducev->centre )
		s.height += 1;
	if( vips_region_prepare( ir, &s ) )
		return( -1 );

#ifdef DEBUG_PIXELS
	printf( "vips_reducev_vector_gen: preparing %d x %d at %d x %d\n",
		s.width, s.height, s.left, s.top ); 
#endif /*DEBUG_PIXELS*/

	for( int i = 0; i < reducev->n_pass; i++ ) 
		vips_executor_set_program( &executor[i], 
			reducev->pass[i].vector, ne );

	VIPS_GATE_START( "vips_reducev_vector_gen: work" ); 

	for( int y = 0; y < r->height; y ++ ) { 
		VipsPel *q = 
			VIPS_REGION_ADDR( out_region, r->left, r->top + y );
		const double Y = (r->top + y) * reducev->vshrink + 
			(reducev->centre ? 0.5 : 0.0); 
		const int py = (int) Y; 
		const int sy = Y * VIPS_TRANSFORM_SCALE * 2;
		const int siy = sy & (VIPS_TRANSFORM_SCALE * 2 - 1);
		const int ty = (siy + 1) >> 1;
		const int *cyo = reducev->matrixo[ty];

#ifdef DEBUG_PIXELS
		printf( "starting row %d\n", y + r->top ); 
		printf( "coefficients:\n" );
		for( int i = 0; i < reducev->n_point; i++ ) 
			printf( "\t%d - %d\n", i, cyo[i] );
		printf( "first column of pixel values:\n" ); 
		for( int i = 0; i < reducev->n_point; i++ ) 
			printf( "\t%d - %d\n", i, 
				*VIPS_REGION_ADDR( ir, r->left, r->top + y + i ) ); 
#endif /*DEBUG_PIXELS*/

		/* We run our n passes to generate this scanline.
		 */
		for( int i = 0; i < reducev->n_pass; i++ ) {
			Pass *pass = &reducev->pass[i]; 

			vips_executor_set_scanline( &executor[i], 
				ir, r->left, py );
			vips_executor_set_array( &executor[i],
				pass->r, seq->t1 );
			vips_executor_set_array( &executor[i],
				pass->d2, seq->t2 );
			for( int j = 0; j < pass->n_param; j++ ) 
				vips_executor_set_parameter( &executor[i],
					pass->p[j], cyo[j + pass->first] ); 
			vips_executor_set_destination( &executor[i], q );
			vips_executor_run( &executor[i] );

			VIPS_SWAP( signed short *, seq->t1, seq->t2 );
		}

#ifdef DEBUG_PIXELS
		printf( "pixel result:\n" );
		printf( "\t%d\n", *q ); 
#endif /*DEBUG_PIXELS*/
	}

	VIPS_GATE_STOP( "vips_reducev_vector_gen: work" ); 

	VIPS_COUNT_PIXELS( out_region, "vips_reducev_vector_gen" ); 

	return( 0 );
}

static int
vips_reducev_raw( VipsReducev *reducev, VipsImage *in, VipsImage **out ) 
{
	VipsObjectClass *object_class = VIPS_OBJECT_GET_CLASS( reducev );
	VipsResample *resample = VIPS_RESAMPLE( reducev );

	VipsGenerateFn generate;

	/* Build masks.
	 */
	for( int y = 0; y < VIPS_TRANSFORM_SCALE + 1; y++ ) {
		reducev->matrixf[y] = 
			VIPS_ARRAY( NULL, reducev->n_point, double ); 
		if( !reducev->matrixf[y] )
			return( -1 ); 

		vips_reduce_make_mask( reducev->matrixf[y],
			reducev->kernel, reducev->vshrink, 
			(float) y / VIPS_TRANSFORM_SCALE ); 

#ifdef DEBUG
		printf( "%6.2g", (double) y / VIPS_TRANSFORM_SCALE ); 
		for( int i = 0; i < reducev->n_point; i++ ) 
			printf( ", %6.2g", reducev->matrixf[y][i] );
		printf( "\n" ); 
#endif /*DEBUG*/
	}

	/* uchar and ushort need an int version of the masks.
	 */
	if( VIPS_IMAGE_SIZEOF_ELEMENT( in ) <= 2 ) 
		for( int y = 0; y < VIPS_TRANSFORM_SCALE + 1; y++ ) {
			reducev->matrixi[y] = 
				VIPS_ARRAY( NULL, reducev->n_point, int ); 
			if( !reducev->matrixi[y] )
				return( -1 ); 

			vips_vector_to_fixed_point( 
				reducev->matrixf[y], reducev->matrixi[y], 
				reducev->n_point, VIPS_INTERPOLATE_SCALE );
		}

	/* And we need an 2.6 version if we will use the vector path.
	 */
	if( in->BandFmt == VIPS_FORMAT_UCHAR &&
		vips_vector_isenabled() ) 
		for( int y = 0; y < VIPS_TRANSFORM_SCALE + 1; y++ ) {
			reducev->matrixo[y] = 
				VIPS_ARRAY( NULL, reducev->n_point, int ); 
			if( !reducev->matrixo[y] )
				return( -1 ); 

			vips_vector_to_fixed_point( 
				reducev->matrixf[y], reducev->matrixo[y], 
				reducev->n_point, 64 );
		}

	/* Try to build a vector version, if we can.
	 */
	generate = vips_reducev_gen;
//	if( in->BandFmt == VIPS_FORMAT_UCHAR &&
//		vips_vector_isenabled() &&
//		!vips_reducev_compile( reducev ) ) {
//		g_info( "reducev: using vector path" );
//		generate = vips_reducev_vector_gen;
//	}

	*out = vips_image_new();
	if( vips_image_pipelinev( *out, 
		VIPS_DEMAND_STYLE_FATSTRIP, in, (void *) NULL ) )
		return( -1 );

	/* Size output. We need to always round to nearest, so round(), not
	 * rint().
	 *
	 * Don't change xres/yres, leave that to the application layer. For
	 * example, vipsthumbnail knows the true reduce factor (including the
	 * fractional part), we just see the integer part here.
	 */
	(*out)->Ysize = VIPS_ROUND_UINT( 
		resample->in->Ysize / reducev->vshrink );
	if( (*out)->Ysize <= 0 ) { 
		vips_error( object_class->nickname, 
			"%s", _( "image has shrunk to nothing" ) );
		return( -1 );
	}

#ifdef DEBUG
	printf( "vips_reducev_build: reducing %d x %d image to %d x %d\n", 
		in->Xsize, in->Ysize, 
		(*out)->Xsize, (*out)->Ysize );  
#endif /*DEBUG*/

	if( vips_image_generate( *out,
		vips_reducev_start, generate, vips_reducev_stop, 
		in, reducev ) )
		return( -1 );

	vips_reorder_margin_hint( *out, reducev->n_point ); 

	return( 0 );
}

static int
vips_reducev_build( VipsObject *object )
{
	VipsObjectClass *object_class = VIPS_OBJECT_GET_CLASS( object );
	VipsResample *resample = VIPS_RESAMPLE( object );
	VipsReducev *reducev = (VipsReducev *) object;
	VipsImage **t = (VipsImage **) vips_object_local_array( object, 6 );

	VipsImage *in;
	int height;

	if( VIPS_OBJECT_CLASS( vips_reducev_parent_class )->build( object ) )
		return( -1 );

	in = resample->in; 

	if( reducev->vshrink < 1 ) { 
		vips_error( object_class->nickname, 
			"%s", _( "reduce factor should be >= 1" ) );
		return( -1 );
	}

	if( reducev->vshrink == 1 ) 
		return( vips_image_write( in, resample->out ) );

	reducev->n_point = 
		vips_reduce_get_points( reducev->kernel, reducev->vshrink ); 
	g_info( "reducev: %d point mask", reducev->n_point );
	if( reducev->n_point > MAX_POINT ) {
		vips_error( object_class->nickname, 
			"%s", _( "reduce factor too large" ) );
		return( -1 );
	}

	/* Unpack for processing.
	 */
	if( vips_image_decode( in, &t[0] ) )
		return( -1 );
	in = t[0];

	//TEMP TEMP TEMP
	if( vips_colourspace(in, &t[4], VIPS_INTERPRETATION_RGB16, NULL))
		return( -1 );
	in = t[4];

	if( vips_cast_double(in, &t[5], NULL))
		return( -1 );
	in = t[5];
	//TEMP TEMP TEMP

	/* Add new pixels around the input so we can interpolate at the edges.
	 */
	height = in->Ysize + reducev->n_point - 1;
	if( reducev->centre )
		height += 1;
	if( vips_embed( in, &t[1], 
		0, reducev->n_point / 2 - 1, 
		in->Xsize, height, 
		"extend", VIPS_EXTEND_COPY,
		(void *) NULL ) )
		return( -1 );
	in = t[1];

	if( vips_reducev_raw( reducev, in, &t[2] ) )
		return( -1 );
	in = t[2];

	/* Large reducev will throw off sequential mode. Suppose thread1 is
	 * generating tile (0, 0), but stalls. thread2 generates tile
	 * (0, 1), 128 lines further down the output. After it has done,
	 * thread1 tries to generate (0, 0), but by then the pixels it needs
	 * have gone from the input image line cache if the reducev is large.
	 *
	 * To fix this, put another seq on the output of reducev. Now we'll
	 * always have the previous XX lines of the shrunk image, and we won't
	 * fetch out of order. 
	 */
	if( vips_image_get_typeof( in, VIPS_META_SEQUENTIAL ) ) { 
		g_info( "reducev sequential line cache" ); 

		if( vips_sequential( in, &t[3], 
			"tile_height", 10,
			// "trace", TRUE,
			(void *) NULL ) )
			return( -1 );
		in = t[3];
	}

	if( vips_image_write( in, resample->out ) )
		return( -1 ); 

	return( 0 );
}

static void
vips_reducev_class_init( VipsReducevClass *reducev_class )
{
	GObjectClass *gobject_class = G_OBJECT_CLASS( reducev_class );
	VipsObjectClass *vobject_class = VIPS_OBJECT_CLASS( reducev_class );
	VipsOperationClass *operation_class = 
		VIPS_OPERATION_CLASS( reducev_class );

	VIPS_DEBUG_MSG( "vips_reducev_class_init\n" );

	gobject_class->finalize = vips_reducev_finalize;
	gobject_class->set_property = vips_object_set_property;
	gobject_class->get_property = vips_object_get_property;

	vobject_class->nickname = "reducev";
	vobject_class->description = _( "shrink an image vertically" );
	vobject_class->build = vips_reducev_build;

	operation_class->flags = VIPS_OPERATION_SEQUENTIAL;

	VIPS_ARG_DOUBLE( reducev_class, "vshrink", 3, 
		_( "Vshrink" ), 
		_( "Vertical shrink factor" ),
		VIPS_ARGUMENT_REQUIRED_INPUT,
		G_STRUCT_OFFSET( VipsReducev, vshrink ),
		1, 1000000, 1 );

	VIPS_ARG_ENUM( reducev_class, "kernel", 4, 
		_( "Kernel" ), 
		_( "Resampling kernel" ),
		VIPS_ARGUMENT_OPTIONAL_INPUT,
		G_STRUCT_OFFSET( VipsReducev, kernel ),
		VIPS_TYPE_KERNEL, VIPS_KERNEL_LANCZOS3 );

	VIPS_ARG_BOOL( reducev_class, "centre", 7, 
		_( "Centre" ), 
		_( "Use centre sampling convention" ),
		VIPS_ARGUMENT_OPTIONAL_INPUT,
		G_STRUCT_OFFSET( VipsReducev, centre ),
		FALSE );

	/* Old name.
	 */
	VIPS_ARG_DOUBLE( reducev_class, "yshrink", 3, 
		_( "Yshrink" ), 
		_( "Vertical shrink factor" ),
		VIPS_ARGUMENT_REQUIRED_INPUT | VIPS_ARGUMENT_DEPRECATED,
		G_STRUCT_OFFSET( VipsReducev, vshrink ),
		1, 1000000, 1 );

}

static void
vips_reducev_init( VipsReducev *reducev )
{
	reducev->kernel = VIPS_KERNEL_LANCZOS3;
}

/* See reduce.c for the doc comment.
 */

int
vips_reducev( VipsImage *in, VipsImage **out, double vshrink, ... )
{
	va_list ap;
	int result;

	va_start( ap, vshrink );
	result = vips_call_split( "reducev", ap, in, out, vshrink );
	va_end( ap );

	return( result );
}
