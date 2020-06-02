//
// Created by Alon Neubach on 01/06/2020.
//

#ifndef LIBVIPS_REDUCE_H
#define LIBVIPS_REDUCE_H

#include <float.h>
#define EPSILON  (1.0e-12)

static inline void
calculate_filter( double factor, int destination_index, int max_source_size,
                  double *filter, int *out_filter_size, int *out_filter_start)
{
	const int filter_support = 3;
	double support = factor * filter_support;

	double bisect = (double) (destination_index + 0.5) * factor + EPSILON;
	int filter_stop = (int) VIPS_MIN( bisect + support + 0.5, max_source_size );
	int filter_start = (int) VIPS_MAX( bisect - support + 0.5, 0.0 );
	int filter_size = filter_stop - filter_start;

	*out_filter_size = filter_size;
	*out_filter_start = filter_start;

	if (filter == NULL) {
		return;
	}

	const double resize_filter_scale = (1.0 / 3.0);
	double density = 0;

	for( int i = 0; i < filter_size; i++ ) {
		double wx = VIPS_ABS(
			((double) (filter_start + i) - bisect + 0.5) / factor );
		filter[i] = sinc_fast( wx * resize_filter_scale ) * sinc_fast( wx );
		density += filter[i];
	}

	if( (density != 0.0) && (density != 1.0) ) {
		/*
		  Normalize.
		*/
		density = 1 / density;
		for( int i = 0; i < filter_size; i++ ) {
			filter[i] *= density;
		}
	}
}

template <typename T>
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

	return destination_pixel;
}

template <typename T, int max_value>
static double
apply_filter_with_alpha( int stride, const double *weights, int n,
                         int band_index, int alpha_index,
                         const VipsPel *p )
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

	return destination_pixel / alpha_sum;
}


template <typename T, bool clip, int min_value, int max_value>
static void
reduce_inner_dimension( const VipsImage *in, const double *filter,
                        int filter_size, const int filter_stride,
                        int inner_dimension_size, const VipsPel *p,
                        VipsPel *q, int source_stride,
                        int destination_stride )
{
	if (filter_size == 0) {
		return;
	}

	gboolean has_alpha = vips_image_hasalpha( (VipsImage *) in );
	const int num_bands = in->Bands *
	                      (vips_band_format_iscomplex( in->BandFmt ) ?  2 : 1);
	int alpha_index = num_bands - 1;

	for( int i = 0; i < inner_dimension_size; i++ ) {
		for( int band_index = 0; band_index < num_bands; band_index++ ) {
			double pixel;

			if( !has_alpha || band_index == alpha_index ) {
				pixel = apply_filter_no_alpha<T>(
					filter_stride, filter, filter_size, band_index, p );
			} else {
				pixel = apply_filter_with_alpha<T, max_value>(
					filter_stride, filter, filter_size, band_index, alpha_index, p );
			}

			((T *) q)[band_index] = clip ? VIPS_CLIP(min_value, pixel, max_value) : pixel;
		} // for band_index

		p += source_stride;
		q += destination_stride;
	} // for i
}

static void reduce_inner_dimension_band_fmt(const VipsImage *in, const double *filter,
                                            int filter_size, const int filter_stride,
                                            int inner_dimension_size, const VipsPel *p,
                                            VipsPel *q, int source_inner_stride,
                                            int destination_inner_stride)
{
	switch( in->BandFmt ) {
		case VIPS_FORMAT_UCHAR:
			reduce_inner_dimension
				<unsigned char, true, 0, UCHAR_MAX>(in, filter, filter_size, filter_stride, inner_dimension_size, p, q, source_inner_stride, destination_inner_stride);
			break;

		case VIPS_FORMAT_CHAR:
			reduce_inner_dimension
				<signed char, true, SCHAR_MIN, SCHAR_MAX>(in, filter, filter_size, filter_stride, inner_dimension_size, p, q, source_inner_stride, destination_inner_stride);
			break;

		case VIPS_FORMAT_USHORT:
			reduce_inner_dimension
				<unsigned short, true, 0, USHRT_MAX>(in, filter, filter_size, filter_stride, inner_dimension_size, p, q, source_inner_stride, destination_inner_stride);
			break;

		case VIPS_FORMAT_SHORT:
			reduce_inner_dimension
				<signed short, true, SHRT_MIN, SHRT_MAX>(in, filter, filter_size, filter_stride, inner_dimension_size, p, q, source_inner_stride, destination_inner_stride);
			break;

		case VIPS_FORMAT_UINT:
			reduce_inner_dimension
				<unsigned int, true, 0, INT_MAX>(in, filter, filter_size, filter_stride, inner_dimension_size, p, q, source_inner_stride, destination_inner_stride);
			break;

		case VIPS_FORMAT_INT:
			reduce_inner_dimension
				<signed int, true, INT_MIN, INT_MAX>(in, filter, filter_size, filter_stride, inner_dimension_size, p, q, source_inner_stride, destination_inner_stride);
			break;

		case VIPS_FORMAT_FLOAT:
		case VIPS_FORMAT_COMPLEX:
			reduce_inner_dimension<float, false, 0, 1>(in, filter, filter_size, filter_stride, inner_dimension_size, p, q, source_inner_stride, destination_inner_stride);
			break;

		case VIPS_FORMAT_DOUBLE:
		case VIPS_FORMAT_DPCOMPLEX:
			reduce_inner_dimension<double, false, 0, 1>(in, filter, filter_size, filter_stride, inner_dimension_size, p, q, source_inner_stride, destination_inner_stride);
			break;

		default:
			g_assert_not_reached();
			break;
	}
}
#endif //LIBVIPS_REDUCE_H
