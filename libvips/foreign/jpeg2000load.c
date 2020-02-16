/* load jpeg2000 images with OpenJPEG
 *
 * 10/2/20
 * 	- from heifload.c
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
#define DEBUG_VERBOSE
#define VIPS_DEBUG
#define DEBUG
 */

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif /*HAVE_CONFIG_H*/

#include <vips/intl.h>

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <vips/vips.h>
#include <vips/debug.h>
#include <vips/internal.h>

#ifdef HAVE_JPEG2000_DECODER

#include <openjpeg.h>

#include "pforeign.h"

#define JP2_RFC3745_MAGIC "\x00\x00\x00\x0c\x6a\x50\x20\x20\x0d\x0a\x87\x0a"
#define JP2_MAGIC "\x0d\x0a\x87\x0a"
#define J2K_CODESTREAM_MAGIC "\xff\x4f\xff\x51"
#define J2K_CFMT 1
#define JP2_CFMT 2

#define VIPS_TYPE_FOREIGN_LOAD_JPEG2000 (vips_foreign_load_jpeg2000_get_type())
#define VIPS_FOREIGN_LOAD_JPEG2000(obj) \
    (G_TYPE_CHECK_INSTANCE_CAST( (obj), \
    VIPS_TYPE_FOREIGN_LOAD_JPEG2000, VipsForeignLoadJPEG2000 ))
#define VIPS_FOREIGN_LOAD_JPEG2000_CLASS(klass) \
    (G_TYPE_CHECK_CLASS_CAST( (klass), \
    VIPS_TYPE_FOREIGN_LOAD_JPEG2000, VipsForeignLoadJPEG2000Class))
#define VIPS_IS_FOREIGN_LOAD_JPEG2000(obj) \
    (G_TYPE_CHECK_INSTANCE_TYPE( (obj), VIPS_TYPE_FOREIGN_LOAD_JPEG2000 ))
#define VIPS_IS_FOREIGN_LOAD_JPEG2000_CLASS(klass) \
    (G_TYPE_CHECK_CLASS_TYPE( (klass), VIPS_TYPE_FOREIGN_LOAD_JPEG2000 ))
#define VIPS_FOREIGN_LOAD_JPEG2000_GET_CLASS(obj) \
    (G_TYPE_INSTANCE_GET_CLASS( (obj), \
    VIPS_TYPE_FOREIGN_LOAD_JPEG2000, VipsForeignLoadJPEG2000Class ))

typedef struct _VipsForeignLoadJPEG2000 {
    VipsForeignLoad parent_object;

    gboolean has_alpha;

    OPJ_COLOR_SPACE color_space;

    /* Size of image.
     */
    int width;
    int height;

    int decode_format;
    opj_stream_t *read_stream;
    opj_codec_t *codec;
    opj_image_t *jpeg_image;

} VipsForeignLoadJPEG2000;

typedef struct _VipsForeignLoadJPEG2000Class {
    VipsForeignLoadClass parent_class;

    /* Open the reader, eg. call opj_read_header() etc. This
     * has to be a vfunc so generate can restart after minimise.
     */
    int (*open)(VipsForeignLoadJPEG2000 *jpeg2000);

} VipsForeignLoadJPEG2000Class;

G_DEFINE_ABSTRACT_TYPE(VipsForeignLoadJPEG2000, vips_foreign_load_jpeg2000,
                       VIPS_TYPE_FOREIGN_LOAD);

static void
vips_foreign_load_jpeg2000_close(VipsForeignLoadJPEG2000 *jpeg2000) {
    VIPS_FREEF(opj_stream_destroy, jpeg2000->read_stream);
    VIPS_FREEF(opj_destroy_codec, jpeg2000->codec);
    VIPS_FREEF(opj_image_destroy, jpeg2000->jpeg_image);
}

static void
vips_foreign_load_jpeg2000_dispose(GObject *gobject) {
    VipsForeignLoadJPEG2000 *jpeg2000 = (VipsForeignLoadJPEG2000 *) gobject;

    vips_foreign_load_jpeg2000_close(jpeg2000);

    G_OBJECT_CLASS(vips_foreign_load_jpeg2000_parent_class)->
            dispose(gobject);
}

void
vips__jpeg2000_error(const char *fmt, const char *message) {
    vips_error("jpeg2000", fmt, message);
}

static int
vips_foreign_load_jpeg2000_format(const char *buf, int len) {
    if (len >= 12) {
        if (memcmp(buf, JP2_RFC3745_MAGIC, 12) == 0 || memcmp(buf, JP2_MAGIC, 4) == 0) {
            return JP2_CFMT;
        } else if (memcmp(buf, J2K_CODESTREAM_MAGIC, 4) == 0) {
            return J2K_CFMT;
        } else
            return (-1);
    }

    return (0);
}

static int
vips_foreign_load_jpeg2000_is_a(const char *buf, int len) {
    return vips_foreign_load_jpeg2000_format(buf, len) > 0;
}

static VipsForeignFlags
vips_foreign_load_jpeg2000_get_flags(VipsForeignLoad *load) {
    return (VIPS_FOREIGN_SEQUENTIAL);
}

static int
vips_foreign_load_jpeg2000_set_header(VipsForeignLoadJPEG2000 *jpeg2000, VipsImage *out) {
    int bands;

    // todo: jpeg2000->has_alpha =

    bands = jpeg2000->has_alpha ? 4 : 3;

    // todo .. EXIF, XMP and IPCT.

    /* FIXME should probably check the profile type ... lcms seems to be
     * able to load at least rICC and prof.
     */
    if(jpeg2000->jpeg_image->icc_profile_len) {
        vips_image_set_blob(out, VIPS_META_ICC_NAME,
            (VipsCallbackFn) NULL,
            jpeg2000->jpeg_image->icc_profile_buf,
            jpeg2000->jpeg_image->icc_profile_len);
    }

    /* todo: VIPS_INTERPRETATION_sRGB is not always the right answer...
     * all grey images, perhaps.
     */
    vips_image_pipelinev(out, VIPS_DEMAND_STYLE_SMALLTILE, NULL);
    vips_image_init_fields(out, jpeg2000->width, jpeg2000->height, bands,
                           VIPS_FORMAT_UCHAR, VIPS_CODING_NONE, VIPS_INTERPRETATION_sRGB,
                           1.0, 1.0);

    return (0);
}

static int
vips_foreign_load_jpeg2000_header(VipsForeignLoad *load) {
    VipsObjectClass *class = VIPS_OBJECT_GET_CLASS(load);
    VipsForeignLoadJPEG2000 *jpeg2000 = (VipsForeignLoadJPEG2000 *) load;
    VipsForeignLoadJPEG2000Class *jpeg2000_class =
            VIPS_FOREIGN_LOAD_JPEG2000_GET_CLASS(jpeg2000);

    if (jpeg2000_class->open(jpeg2000))
        return (-1);

    /* Read the main header of the codestream and if necessary the JP2 boxes*/
    if (! opj_read_header(jpeg2000->read_stream, jpeg2000->codec, &jpeg2000->jpeg_image)) {
        vips__jpeg2000_error("failed to read the header", NULL);
        vips_foreign_load_jpeg2000_close(jpeg2000);
        return (-1);
    }

    jpeg2000->width = jpeg2000->jpeg_image->x1;
    jpeg2000->height = jpeg2000->jpeg_image->y1;
    jpeg2000->color_space = jpeg2000->jpeg_image->color_space;

    //todo: has_alpha
    //todo: colorspace

    if (vips_foreign_load_jpeg2000_set_header(jpeg2000, load->out)) {
        vips_foreign_load_jpeg2000_close(jpeg2000);
        return (-1);
    }

    vips_foreign_load_jpeg2000_close(jpeg2000);

    return (0);
}

static int
vips_foreign_load_jpeg2000_generate(VipsRegion *or,
                                    void *seq, void *a, void *b, gboolean *stop) {
    VipsForeignLoadJPEG2000 *jpeg2000 = (VipsForeignLoadJPEG2000 *) a;
    VipsObjectClass *class = VIPS_OBJECT_GET_CLASS(jpeg2000);
    VipsForeignLoadJPEG2000Class *jpeg2000_class =
            VIPS_FOREIGN_LOAD_JPEG2000_GET_CLASS(jpeg2000);
    VipsRect *r = &or->valid;

    opj_codec_t* codec = NULL;

    g_assert(r->height == 1);

    if (jpeg2000_class->open(jpeg2000))
        return (-1);

    switch (jpeg2000->decode_format) {
        case J2K_CFMT: { /* JPEG-2000 codestream */
            /* Get a decoder handle */
            codec = opj_create_decompress(OPJ_CODEC_J2K);
            break;
        }
        case JP2_CFMT: { /* JPEG 2000 compressed image data */
            /* Get a decoder handle */
            codec = opj_create_decompress(OPJ_CODEC_JP2);
            break;
        }
        default:
            vips__jpeg2000_error(&"failed to create codec for format " [ (jpeg2000->decode_format + '0')], NULL);
            vips_foreign_load_jpeg2000_close(jpeg2000);
            return (-1);
    }

    if (!jpeg2000->img) {
        struct jpeg2000_decoding_options *options;
        enum jpeg2000_chroma chroma = jpeg2000->has_alpha ?
                                      jpeg2000_chroma_interleaved_RGBA :
                                      jpeg2000_chroma_interleaved_RGB;

        /* Only disable transforms if we have been able to fetch the
         * untransformed dimensions.
         */
        options = jpeg2000_decoding_options_alloc();

        error = jpeg2000_decode_image(jpeg2000->handle, &jpeg2000->img,
                                      jpeg2000_colorspace_RGB, chroma,
                                      options);
        jpeg2000_decoding_options_free(options);
        if (error.code) {
            vips__jpeg2000_error(&error);
            return (-1);
        }
    }

    memcpy(VIPS_REGION_ADDR(or, 0, r->top),
           jpeg2000->data + jpeg2000->stride * line,
           VIPS_IMAGE_SIZEOF_LINE(or->im));

    return (0);
}

static void
vips_foreign_load_jpeg2000_minimise(VipsObject *object, VipsForeignLoadJPEG2000 *jpeg2000) {
    vips_foreign_load_jpeg2000_close(jpeg2000);
}

static int
vips_foreign_load_jpeg2000_load(VipsForeignLoad *load) {
    VipsForeignLoadJPEG2000 *jpeg2000 = (VipsForeignLoadJPEG2000 *) load;
    VipsForeignLoadJPEG2000Class *class =
            VIPS_FOREIGN_LOAD_JPEG2000_GET_CLASS(jpeg2000);

    VipsImage **t = (VipsImage **)
            vips_object_local_array(VIPS_OBJECT(load), 3);

    if (class->open(jpeg2000))
        return (-1);

    t[0] = vips_image_new();
    if (vips_foreign_load_jpeg2000_set_header(jpeg2000, t[0]))
        return (-1);

    /* Close input immediately at end of read.
     */
    g_signal_connect(t[0], "minimise",
                     G_CALLBACK(vips_foreign_load_jpeg2000_minimise), jpeg2000);

    if (vips_image_generate(t[0],
                            NULL, vips_foreign_load_jpeg2000_generate, NULL, jpeg2000, NULL) ||
        vips_sequential(t[0], &t[1], NULL) ||
        vips_image_write(t[1], load->real))
        return (-1);

    return (0);
}

static int
vips_foreign_load_jpeg2000_open(VipsForeignLoadJPEG2000 *jpeg2000)
{
    return (0);
}

static void
vips_foreign_load_jpeg2000_class_init(VipsForeignLoadJPEG2000Class *class)
{
    GObjectClass *gobject_class = G_OBJECT_CLASS(class);
    VipsObjectClass *object_class = (VipsObjectClass *) class;
    VipsForeignLoadClass *load_class = (VipsForeignLoadClass *) class;
    VipsForeignLoadJPEG2000Class *jpeg2000_class =
            (VipsForeignLoadJPEG2000Class *) class;

    gobject_class->dispose = vips_foreign_load_jpeg2000_dispose;
    gobject_class->set_property = vips_object_set_property;
    gobject_class->get_property = vips_object_get_property;

    object_class->nickname = "jpeg2000load_base";
    object_class->description = _("load a JPEG2000 image");

    load_class->get_flags = vips_foreign_load_jpeg2000_get_flags;
    load_class->header = vips_foreign_load_jpeg2000_header;
    load_class->load = vips_foreign_load_jpeg2000_load;

    jpeg2000_class->open = vips_foreign_load_jpeg2000_open;
}

static void
vips_foreign_load_jpeg2000_init(VipsForeignLoadJPEG2000 *jpeg2000) {

}

typedef struct _VipsForeignLoadJPEG2000File {
    VipsForeignLoadJPEG2000 parent_object;

    /* Filename for load.
     */
    char *filename;

} VipsForeignLoadJPEG2000File;

typedef VipsForeignLoadJPEG2000Class VipsForeignLoadJPEG2000FileClass;

G_DEFINE_TYPE(VipsForeignLoadJPEG2000File, vips_foreign_load_jpeg2000_file,
              vips_foreign_load_jpeg2000_get_type());

static int
vips_foreign_load_jpeg2000_file_is_a(const char *filename) {
    char buf[12];

    if (vips__get_bytes(filename, (unsigned char *) buf, 12) != 12)
        return (0);

    return (vips_foreign_load_jpeg2000_is_a(buf, 12));
}

static int
vips_foreign_load_jpeg2000_file_header(VipsForeignLoad *load) {
    VipsForeignLoadJPEG2000 *jpeg2000 = (VipsForeignLoadJPEG2000 *) load;
    VipsForeignLoadJPEG2000File *file = (VipsForeignLoadJPEG2000File *) load;

    if (VIPS_FOREIGN_LOAD_CLASS(
            vips_foreign_load_jpeg2000_file_parent_class)->header(load)) {
        /* Close early if our base class fails to read.
         */
        vips_foreign_load_jpeg2000_close(jpeg2000);
        return (-1);
    }

    VIPS_SETSTR(load->out->filename, file->filename);

    return (0);
}

const char *vips__jpeg2000_suffs[] = {
        ".jp2",
        ".jpg2",
        ".j2k",
        ".jpx",
        ".jpf",
        ".jpc",
        ".j2c",
//        ".jpt",
        NULL
};

static int
vips_foreign_load_jpeg2000_file_open(VipsForeignLoadJPEG2000 *jpeg2000)
{
    VipsForeignLoadJPEG2000File *file = (VipsForeignLoadJPEG2000File *) jpeg2000;
    char buf[12];

    if (vips__get_bytes(file->filename, (unsigned char *) buf, 12) != 12)
        return (0);

    jpeg2000->read_stream = opj_stream_create_default_file_stream(file->filename, 1);
    if (!jpeg2000->read_stream) {
        vips__jpeg2000_error("failed to create the stream from the file %s", file->filename);
        return (-1);
    }

    jpeg2000->decode_format = vips_foreign_load_jpeg2000_format(buf, 12);

    return (VIPS_FOREIGN_LOAD_JPEG2000_CLASS(
            vips_foreign_load_jpeg2000_file_parent_class)->open(jpeg2000));
}

static void
vips_foreign_load_jpeg2000_file_class_init(VipsForeignLoadJPEG2000FileClass *class)
{
    GObjectClass *gobject_class = G_OBJECT_CLASS(class);
    VipsObjectClass *object_class = (VipsObjectClass *) class;
    VipsForeignClass *foreign_class = (VipsForeignClass *) class;
    VipsForeignLoadClass *load_class = (VipsForeignLoadClass *) class;
    VipsForeignLoadJPEG2000Class *jpeg2000_class =
            (VipsForeignLoadJPEG2000Class *) class;

    gobject_class->set_property = vips_object_set_property;
    gobject_class->get_property = vips_object_get_property;

    object_class->nickname = "jpeg2000load";

    foreign_class->suffs = vips__jpeg2000_suffs;

    load_class->is_a = vips_foreign_load_jpeg2000_file_is_a;
    load_class->header = vips_foreign_load_jpeg2000_file_header;

    jpeg2000_class->open = vips_foreign_load_jpeg2000_file_open;

    VIPS_ARG_STRING(class, "filename", 1,
                    _("Filename"),
                    _("Filename to load from"),
                    VIPS_ARGUMENT_REQUIRED_INPUT,
                    G_STRUCT_OFFSET(VipsForeignLoadJPEG2000File, filename),
                    NULL);
}

static void
vips_foreign_load_jpeg2000_file_init(VipsForeignLoadJPEG2000File *file)
{
}

typedef struct _VipsForeignLoadJPEG2000Buffer
{
    VipsForeignLoadJPEG2000 parent_object;

    /* Load from a buffer.
     */
    VipsArea *buf;

} VipsForeignLoadJPEG2000Buffer;

typedef VipsForeignLoadJPEG2000Class VipsForeignLoadJPEG2000BufferClass;

G_DEFINE_TYPE(VipsForeignLoadJPEG2000Buffer, vips_foreign_load_jpeg2000_buffer,
              vips_foreign_load_jpeg2000_get_type());

static gboolean
vips_foreign_load_jpeg2000_buffer_is_a(const void *buf, size_t len)
{
    return (vips_foreign_load_jpeg2000_is_a(buf, len));
}

static int
vips_foreign_load_jpeg2000_buffer_open(VipsForeignLoadJPEG2000 *jpeg2000)
{
    VipsForeignLoadJPEG2000Buffer *buffer = (VipsForeignLoadJPEG2000Buffer *) jpeg2000;

    jpeg2000->read_stream = opj_stream_create(buffer->buf->length ,1);
    if (!jpeg2000->read_stream) {
        vips__jpeg2000_error("failed to create the stream from buffer", NULL);
        return (-1);
    }

    // todo: free user data function?
    opj_stream_set_user_data(jpeg2000->read_stream, buffer->buf->data, NULL);

    jpeg2000->decode_format = vips_foreign_load_jpeg2000_format(buffer->buf->data, buffer->buf->length);

    return (VIPS_FOREIGN_LOAD_JPEG2000_CLASS(
            vips_foreign_load_jpeg2000_buffer_parent_class)->open(jpeg2000));
}

static void
vips_foreign_load_jpeg2000_buffer_class_init(
        VipsForeignLoadJPEG2000BufferClass *class)
{
    GObjectClass *gobject_class = G_OBJECT_CLASS(class);
    VipsObjectClass *object_class = (VipsObjectClass *) class;
    VipsForeignLoadClass *load_class = (VipsForeignLoadClass *) class;
    VipsForeignLoadJPEG2000Class *jpeg2000_class =
            (VipsForeignLoadJPEG2000Class *) class;

    gobject_class->set_property = vips_object_set_property;
    gobject_class->get_property = vips_object_get_property;

    object_class->nickname = "jpeg2000load_buffer";

    load_class->is_a_buffer = vips_foreign_load_jpeg2000_buffer_is_a;

    jpeg2000_class->open = vips_foreign_load_jpeg2000_buffer_open;

    VIPS_ARG_BOXED(class, "buffer", 1,
                   _("Buffer"),
                   _("Buffer to load from"),
                   VIPS_ARGUMENT_REQUIRED_INPUT,
                   G_STRUCT_OFFSET(VipsForeignLoadJPEG2000Buffer, buf),
                   VIPS_TYPE_BLOB);
}

static void
vips_foreign_load_jpeg2000_buffer_init(VipsForeignLoadJPEG2000Buffer *buffer)
{
}

#endif /*HAVE_JPEG2000_DECODER*/

/**
 * vips_jpeg2000load:
 * @filename: file to load
 * @out: (out): decompressed image
 * @...: %NULL-terminated list of optional named arguments
 *
 * Optional arguments:
 *
 * Read a JPEG2000 image file into a VIPS image.
 *
 * See also: vips_image_new_from_file().
 *
 * Returns: 0 on success, -1 on error.
 */
int
vips_jpeg2000load(const char *filename, VipsImage **out, ...)
{
    va_list ap;
    int result;

    va_start(ap, out);
    result = vips_call_split("jpeg2000load", ap, filename, out);
    va_end(ap);

    return (result);
}

/**
 * vips_jpeg2000load_buffer:
 * @buf: (array length=len) (element-type guint8): memory area to load
 * @len: (type gsize): size of memory area
 * @out: (out): image to write
 * @...: %NULL-terminated list of optional named arguments
 *
 * Optional arguments:
 *
 *
 * Read a JPEG2000 image file into a VIPS image. 
 * Exactly as vips_jpeg2000load(), but read from a memory buffer. 
 *
 * You must not free the buffer while @out is active. The 
 * #VipsObject::postclose signal on @out is a good place to free. 
 *
 * See also: vips_jpeg2000load().
 *
 * Returns: 0 on success, -1 on error.
 */
int
vips_jpeg2000load_buffer(void *buf, size_t len, VipsImage **out, ...)
{
    va_list ap;
    VipsBlob *blob;
    int result;

    /* We don't take a copy of the data or free it.
     */
    blob = vips_blob_new(NULL, buf, len);

    va_start(ap, out);
    result = vips_call_split("jpeg2000load_buffer", ap, blob, out);
    va_end(ap);

    vips_area_unref(VIPS_AREA(blob));

    return (result);
}
