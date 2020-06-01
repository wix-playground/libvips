//
// Created by Alon Neubach on 01/06/2020.
//

#ifndef LIBVIPS_REDUCE_H
#define LIBVIPS_REDUCE_H

#define EPSILON  (1.0e-12)

static inline double reciprocal( const double x )
{
	double sign;

	/*
	  Return 1/x where x is perceptible (not unlimited or infinitesimal).
	*/
	sign = x < 0.0 ? -1.0 : 1.0;
	if( sign * x >= EPSILON )
		return( 1.0 / x );

	return( sign / EPSILON );
}

static inline void
calculate_weights( double factor, double bisect, int start,
                   double *weights, int n )
{
	const double resize_filter_scale = (1.0 / 3.0);
	double density = 0;
	factor = reciprocal( factor );

	for( int i = 0; i < n; i++ ) {
		double wx = VIPS_ABS(
			((double) (start + i) - bisect + 0.5) * factor );
		weights[i] = sinc_fast( wx * resize_filter_scale ) * sinc_fast( wx );
		density += weights[i];
	}

	if( (density != 0.0) && (density != 1.0) ) {
		/*
		  Normalize.
		*/
		density = 1 / density;
		for( int i = 0; i < n; i++ ) {
			weights[i] *= density;
		}
	}
}

template <typename T, int max_value>
static double
calculate_pixel_with_alpha_blend( int stride, int alpha_index,
                                  const double *weights, int n, int band_index,
                                  const unsigned short *source_bands )
{
	double alpha_sum = 0.0;
	double destination_pixel = 0;

	for( int i = 0; i < n; i++ ) {
		unsigned short source_alpha = source_bands[alpha_index];
		unsigned short source_pixel = source_bands[band_index];
		double alpha = weights[i] * source_alpha / max_value;

		destination_pixel += alpha * source_pixel;
		alpha_sum += alpha;

		source_bands += stride;
	}

	return VIPS_CLIP( 0, destination_pixel / alpha_sum, max_value );
}

template <typename T,int max_value>
static double
calculate_pixel_no_alpha_blend( int stride,
                                const double *weights, int n, int band_index,
                                const unsigned short *source_bands )
{
	double destination_pixel = 0;
	for( int i = 0; i < n; i++ ) {
		T source_pixel = source_bands[band_index];
		destination_pixel += weights[i] * source_pixel;

		source_bands += stride;
	}

	return VIPS_CLIP( 0, destination_pixel, max_value );
}


#endif //LIBVIPS_REDUCE_H
