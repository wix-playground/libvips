/* various interpolation templates
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


#ifndef restrict
#define restrict __restrict
#endif

/*
 * Various casts which assume that the data is already in range. (That
 * is, they are to be used with monotone samplers.)
 */
template <typename T> static T inline
to_fptypes( const double val )
{
	const T newval = val;

	return( newval );
}

template <typename T> static T inline
to_withsign( const double val )
{
	const int sign_of_val = 2 * ( val >= 0. ) - 1;
	const int rounded_abs_val = .5 + sign_of_val * val;
	const T newval = sign_of_val * rounded_abs_val;

	return( newval );
}

template <typename T> static T inline
to_nosign( const double val )
{
	const T newval = .5 + val;

	return( newval );
}

/*
 * Various bilinear implementation templates. Note that no clampling
 * is used: There is an assumption that the data is such that
 * over/underflow is not an issue:
 */

/*
 * Bilinear interpolation for float and double types. The first four
 * inputs are weights, the last four are the corresponding pixel
 * values:
 */
template <typename T> static T inline
bilinear_fptypes(
	const double w_times_z,
	const double x_times_z,
	const double w_times_y,
	const double x_times_y,
	const double tre_thr,
	const double tre_thrfou,
	const double trequa_thr,
	const double trequa_thrfou )
{
	const T newval =
		w_times_z * tre_thr +
		x_times_z * tre_thrfou +
		w_times_y * trequa_thr +
		x_times_y * trequa_thrfou;

	return( newval );
}

/*
 * Bilinear interpolation for signed integer types:
 */
template <typename T> static T inline
bilinear_withsign(
	const double w_times_z,
	const double x_times_z,
	const double w_times_y,
	const double x_times_y,
	const double tre_thr,
	const double tre_thrfou,
	const double trequa_thr,
	const double trequa_thrfou )
{
	const double val =
		w_times_z * tre_thr +
		x_times_z * tre_thrfou +
		w_times_y * trequa_thr +
		x_times_y * trequa_thrfou;

	const int sign_of_val = 2 * ( val >= 0. ) - 1;

	const int rounded_abs_val = .5 + sign_of_val * val;

	const T newval = sign_of_val * rounded_abs_val;

	return( newval );
}

/*
 * Bilinear Interpolation for unsigned integer types:
 */
template <typename T> static T inline
bilinear_nosign(
	const double w_times_z,
	const double x_times_z,
	const double w_times_y,
	const double x_times_y,
	const double tre_thr,
	const double tre_thrfou,
	const double trequa_thr,
	const double trequa_thrfou )
{
	const T newval =
		w_times_z * tre_thr +
		x_times_z * tre_thrfou +
		w_times_y * trequa_thr +
		x_times_y * trequa_thrfou +
		0.5;

	return( newval );
}

/*
 * Bicubic (Catmull-Rom) interpolation templates:
 */

static int inline
unsigned_fixed_round( int v )
{
	const int round_by = VIPS_INTERPOLATE_SCALE >> 1;

	return( (v + round_by) >> VIPS_INTERPOLATE_SHIFT );
}

/* Fixed-point integer bicubic, used for 8-bit types.
 */
template <typename T> static int inline
bicubic_unsigned_int(
	const T uno_one, const T uno_two, const T uno_thr, const T uno_fou,
	const T dos_one, const T dos_two, const T dos_thr, const T dos_fou,
	const T tre_one, const T tre_two, const T tre_thr, const T tre_fou,
	const T qua_one, const T qua_two, const T qua_thr, const T qua_fou,
	const int* restrict cx, const int* restrict cy )
{
	const int c0 = cx[0];
	const int c1 = cx[1];
	const int c2 = cx[2];
	const int c3 = cx[3];

	const int r0 = unsigned_fixed_round( 
		c0 * uno_one +
		c1 * uno_two +
		c2 * uno_thr +
		c3 * uno_fou ); 
	const int r1 = unsigned_fixed_round( 
		c0 * dos_one +
		c1 * dos_two +
		c2 * dos_thr +
		c3 * dos_fou ); 
	const int r2 = unsigned_fixed_round( 
		c0 * tre_one +
		c1 * tre_two +
		c2 * tre_thr +
		c3 * tre_fou ); 
	const int r3 = unsigned_fixed_round( 
		c0 * qua_one +
		c1 * qua_two +
		c2 * qua_thr +
		c3 * qua_fou ); 

	return( unsigned_fixed_round( 
		cy[0] * r0 +
		cy[1] * r1 +
		cy[2] * r2 +
		cy[3] * r3 ) ); 
}

static int inline
signed_fixed_round( int v )
{
	const int sign_of_v = 2 * (v > 0) - 1;
	const int round_by = sign_of_v * (VIPS_INTERPOLATE_SCALE >> 1);

	return( (v + round_by) >> VIPS_INTERPOLATE_SHIFT );
}

/* Fixed-point integer bicubic, used for 8-bit types.
 */
template <typename T> static int inline
bicubic_signed_int(
	const T uno_one, const T uno_two, const T uno_thr, const T uno_fou,
	const T dos_one, const T dos_two, const T dos_thr, const T dos_fou,
	const T tre_one, const T tre_two, const T tre_thr, const T tre_fou,
	const T qua_one, const T qua_two, const T qua_thr, const T qua_fou,
	const int* restrict cx, const int* restrict cy )
{
	const int c0 = cx[0];
	const int c1 = cx[1];
	const int c2 = cx[2];
	const int c3 = cx[3];

	const int r0 = signed_fixed_round( 
		c0 * uno_one +
		c1 * uno_two +
		c2 * uno_thr +
		c3 * uno_fou ); 
	const int r1 = signed_fixed_round( 
		c0 * dos_one +
		c1 * dos_two +
		c2 * dos_thr +
		c3 * dos_fou ); 
	const int r2 = signed_fixed_round( 
		c0 * tre_one +
		c1 * tre_two +
		c2 * tre_thr +
		c3 * tre_fou ); 
	const int r3 = signed_fixed_round( 
		c0 * qua_one +
		c1 * qua_two +
		c2 * qua_thr +
		c3 * qua_fou ); 

	return( signed_fixed_round( 
		cy[0] * r0 +
		cy[1] * r1 +
		cy[2] * r2 +
		cy[3] * r3 ) ); 
}

template <typename T> static T inline
cubic_float(
	const T one, const T two, const T thr, const T fou,
	const double* restrict cx )
{
	return( cx[0] * one +
		 cx[1] * two +
		 cx[2] * thr +
		 cx[3] * fou );
}

/* Floating-point bicubic, used for int/float/double types.
 */
template <typename T> static T inline
bicubic_float(
	const T uno_one, const T uno_two, const T uno_thr, const T uno_fou,
	const T dos_one, const T dos_two, const T dos_thr, const T dos_fou,
	const T tre_one, const T tre_two, const T tre_thr, const T tre_fou,
	const T qua_one, const T qua_two, const T qua_thr, const T qua_fou,
	const double* restrict cx, const double* restrict cy )
{
	const double r0 = cubic_float<T>( 
		uno_one, uno_two, uno_thr, uno_fou, cx ); 
	const double r1 = cubic_float<T>( 
		dos_one, dos_two, dos_thr, dos_fou, cx ); 
	const double r2 = cubic_float<T>( 
		tre_one, tre_two, tre_thr, tre_fou, cx ); 
	const double r3 = cubic_float<T>( 
		qua_one, qua_two, qua_thr, qua_fou, cx ); 

	return( cubic_float<T>( r0, r1, r2, r3, cy ) ); 
}

/* Given an offset in [0,1] (we can have x == 1 when building tables),
 * calculate c0, c1, c2, c3, the catmull-rom coefficients. This is called
 * from the interpolator as well as from the table builder.
 */
static void inline
calculate_coefficients_catmull( double c[4], const double x )
{
	/* Nicolas believes that the following is an hitherto unknown
	 * hyper-efficient method of computing Catmull-Rom coefficients. It
	 * only uses 4* & 1+ & 5- for a total of only 10 flops to compute
	 * four coefficients.
	 */
	const double cr1  = 1. - x;
	const double cr2  = -.5 * x;
	const double cr3  = cr1 * cr2;
	const double cone = cr1 * cr3;
	const double cfou = x * cr3;
	const double cr4  = cfou - cone;
	const double ctwo = cr1 - cone + cr4;
	const double cthr = x - cfou - cr4;

	g_assert( x >= 0. && x <= 1. );

	c[0] = cone;
	c[3] = cfou;
	c[1] = ctwo;
	c[2] = cthr;
}

/* Given an x in [0,1] (we can have x == 1 when building tables),
 * calculate c0 .. c(@shrink + 1), the triangle coefficients. This is called
 * from the interpolator as well as from the table builder.
 */
static void inline
calculate_coefficients_triangle( double *c,
                                 const double shrink, const double x )
{
	/* Needs to be in sync with vips_reduce_get_points().
	 */
	const int n_points = rint( 2 * shrink ) + 1;

	int i;
	double sum;

	sum = 0;
	for( i = 0; i < n_points; i++ ) {
		double xp = (i - (shrink - 0.5) - x) / shrink;

		double l;

		l = 1.0 - VIPS_FABS( xp );
		l = VIPS_MAX( 0.0, l );

		c[i] = l;
		sum += l;
	}

	for( i = 0; i < n_points; i++ )
		c[i] /= sum;
}

/* Generate a cubic filter. See:
 *
 * Mitchell and Netravali, Reconstruction Filters in Computer Graphics
 * Computer Graphics, Volume 22, Number 4, August 1988.
 *
 * B = 1,   C = 0   - cubic B-spline
 * B = 1/3, C = 1/3 - Mitchell
 * B = 0,   C = 1/2 - Catmull-Rom spline
 */
static void inline
calculate_coefficients_cubic( double *c,
                              const double shrink, const double x, double B, double C )
{
	/* Needs to be in sync with vips_reduce_get_points().
	 */
	const int n_points = rint( 4 * shrink ) + 1; 

	int i;
	double sum;

	sum = 0;
	for( i = 0; i < n_points; i++ ) {
		const double xp = (i - (2 * shrink - 1) - x) / shrink;
		const double axp = VIPS_FABS( xp ); 
		const double axp2 = axp * axp;
		const double axp3 = axp2 * axp;

		double l;

		if( axp <= 1 )
			l = ((12 - 9 * B - 6 * C) * axp3 +
			     (-18 + 12 * B + 6 * C) * axp2 +
			     (6 - 2 * B)) / 6;
		else if( axp <= 2 )
			l = ((-B - 6 * C) * axp3 +
			     (6 * B + 30 * C) * axp2 +
			     (-12 * B - 48 * C) * axp +
			     (8 * B + 24 * C)) / 6;
		else
			l = 0.0;

		c[i] = l;
		sum += l;
	}

	for( i = 0; i < n_points; i++ )
		c[i] /= sum;
}

/* Given an x in [0,1] (we can have x == 1 when building tables),
 * calculate c0 .. c(@a * @shrink + 1), the lanczos coefficients. This is called
 * from the interpolator as well as from the table builder.
 *
 * @a is the number of lobes, so usually 2 or 3. @shrink is the reduction
 * factor, so 1 for interpolation, 2 for a x2 reduction, etc. We need more
 * points for large decimations to avoid aliasing.
 */
static void inline
calculate_coefficients_lanczos( double *c,
                                const int a, const double shrink, const double x )
{
	/* Needs to be in sync with vips_reduce_get_points().
	 */
	const int n_points = rint( 2 * a * shrink ) + 1; 

	int i;
	double sum;

	sum = 0;
	for( i = 0; i < n_points; i++ ) {
		double xp = (i - (n_points - 2) / 2 - x) / shrink;

		double l;

		if( xp == 0.0 )
			l = 1.0;
		else if( xp < -a )
			l = 0.0;
		else if( xp > a )
			l = 0.0;
		else
			l = (double) a * sin( VIPS_PI * xp ) *
			    sin( VIPS_PI * xp / (double) a ) /
			    (VIPS_PI * VIPS_PI * xp * xp);

		c[i] = l;
		sum += l;
	}

	for( i = 0; i < n_points; i++ )
		c[i] /= sum;
}

static double sinc_fast(double x)
{
	/*
	  Approximations of the sinc function sin(pi x)/(pi x) over the interval
	  [-4,4] constructed by Nicolas Robidoux and Chantal Racette with funding
	  from the Natural Sciences and Engineering Research Council of Canada.

	  Although the approximations are polynomials (for low order of
	  approximation) and quotients of polynomials (for higher order of
	  approximation) and consequently are similar in form to Taylor polynomials /
	  Pade approximants, the approximations are computed with a completely
	  different technique.

	  Summary: These approximations are "the best" in terms of bang (accuracy)
	  for the buck (flops). More specifically: Among the polynomial quotients
	  that can be computed using a fixed number of flops (with a given "+ - * /
	  budget"), the chosen polynomial quotient is the one closest to the
	  approximated function with respect to maximum absolute relative error over
	  the given interval.

	  The Remez algorithm, as implemented in the boost library's minimax package,
	  is the key to the construction: http://www.boost.org/doc/libs/1_36_0/libs/
	  math/doc/sf_and_dist/html/math_toolkit/backgrounders/remez.html

	  If outside of the interval of approximation, use the standard trig formula.
	*/
	if (x > 4.0)
	{
		return (sin( VIPS_PI * x ) / (VIPS_PI * x));
	}

	/*
	  The approximations only depend on x^2 (sinc is an even function).
	*/
	const double xx = x*x;

	/*
      Max. abs. rel. error 2.2e-8 < 1/2^25.
    */
	const double c0 = 0.173611107357320220183368594093166520811e-2L;
	const double c1 = -0.384240921114946632192116762889211361285e-3L;
	const double c2 = 0.394201182359318128221229891724947048771e-4L;
	const double c3 = -0.250963301609117217660068889165550534856e-5L;
	const double c4 = 0.111902032818095784414237782071368805120e-6L;
	const double c5 = -0.372895101408779549368465614321137048875e-8L;
	const double c6 = 0.957694196677572570319816780188718518330e-10L;
	const double c7 = -0.187208577776590710853865174371617338991e-11L;
	const double c8 = 0.253524321426864752676094495396308636823e-13L;
	const double c9 = -0.177084805010701112639035485248501049364e-15L;
	const double p =
		c0+xx*(c1+xx*(c2+xx*(c3+xx*(c4+xx*(c5+xx*(c6+xx*(c7+xx*(c8+xx*c9))))))));
	return((xx-1.0)*(xx-4.0)*(xx-9.0)*(xx-16.0)*p);

}

static void inline
calculate_coefficients_approx_lanczos( double *c,
                                       const int a, const double shrink, const double x )
{
	/* Needs to be in sync with vips_reduce_get_points().
	*/
	const int n_points = rint( 2 * a * shrink ) + 1;


	int i;
	double sum;

	sum = 0;
	for( i = 0; i < n_points; i++ ) {
		double xp = (i - (n_points - 2) / 2 - x) / shrink;

		double l;

		if( xp == 0.0 )
			l = 1.0;
		else if( xp < -a )
			l = 0.0;
		else if( xp > a )
			l = 0.0;
		else
//			l = (double) a * sin( VIPS_PI * xp ) *
//			    sin( VIPS_PI * xp / (double) a ) /
//			    (VIPS_PI * VIPS_PI * xp * xp);
			l = sinc_fast(xp) * sinc_fast(xp / a);

		c[i] = l;
		sum += l;
	}

	for( i = 0; i < n_points; i++ )
		c[i] /= sum;
}

/* Our inner loop for resampling with a convolution. Operate on elements of 
 * type T, gather results in an intermediate of type IT.
 */
template <typename T, typename IT>
static IT
reduce_sum( const T * restrict in, int stride, const IT * restrict coefficients, int n )
{
	IT sum;

	sum = 0; 
	for( int i = 0; i < n; i++ ) {
		sum += coefficients[i] * in[0];
		in += stride;
	}

	return( sum ); 
}
