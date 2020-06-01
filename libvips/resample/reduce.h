//
// Created by Alon Neubach on 01/06/2020.
//

#ifndef LIBVIPS_REDUCE_H
#define LIBVIPS_REDUCE_H

static inline void
calculate_filter( double factor, double bisect, int filter_start,
                  double *weights, int n )
{
	const double resize_filter_scale = (1.0 / 3.0);
	double density = 0;

	for( int i = 0; i < n; i++ ) {
		double wx = VIPS_ABS(
			((double) (filter_start + i) - bisect + 0.5) / factor );
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

template <typename T,int max_value>
static double
apply_filter_no_alpha( int stride,
                       const double *weights, int n, int band_index,
                       const VipsPel* p )
{
	double destination_pixel = 0;
	for( int i = 0; i < n; i++ ) {
		T source_pixel = ((T*)p)[band_index];
		destination_pixel += weights[i] * source_pixel;

		p += stride;
	}

	return VIPS_CLIP( 0, destination_pixel, max_value );
}

template <typename T, int max_value>
static double
apply_filter_with_alpha( int stride, int alpha_index,
                         const double *weights, int n, int band_index,
                         const VipsPel* p )
{
	double alpha_sum = 0.0;
	double destination_pixel = 0;

	for( int i = 0; i < n; i++ ) {
		T source_alpha = ((T*)p)[alpha_index];
		T source_pixel = ((T*)p)[band_index];
		double alpha = weights[i] * source_alpha / max_value;

		destination_pixel += alpha * source_pixel;
		alpha_sum += alpha;

		p += stride;
	}

	return VIPS_CLIP( 0, destination_pixel / alpha_sum, max_value );
}


template <typename T, int max_value>
static void reduce_inner_dimension( const VipsImage *in, const double *filter,
                                    int filter_size, const int filter_stride,
                                    int inner_dimension_size, const VipsPel *p,
                                    VipsPel *q, int source_inner_stride,
                                    int destination_inner_stride )
{
	gboolean has_alpha = vips_image_hasalpha( (VipsImage *) in );
	const int num_bands = in->Bands *
	                      (vips_band_format_iscomplex( in->BandFmt ) ?  2 : 1);
	int alpha_index = num_bands - 1;

	for( int i = 0; i < inner_dimension_size; i++ ) {
		for( int band_index = 0; band_index < num_bands; band_index++ ) {
			T pixel;

			if( !has_alpha || band_index == alpha_index ) {
				pixel = apply_filter_no_alpha<T, max_value>(
					filter_stride, filter, filter_size, band_index, p );
			} else {
				pixel = apply_filter_with_alpha<T, max_value>(
					filter_stride, alpha_index, filter, filter_size, band_index, p );
			}

			((T *) q)[band_index] = pixel;
		} // for band_index

		p += source_inner_stride;
		q += destination_inner_stride;
	} // for i
}


#endif //LIBVIPS_REDUCE_H
