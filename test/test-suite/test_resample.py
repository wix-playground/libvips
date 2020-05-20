# vim: set fileencoding=utf-8 :
import math

import pytest

import pyvips
from helpers import JPEG_FILE, HEIC_FILE, all_formats, have, LOGO3_PNG_FILE, IMAGES


# Run a function expecting a complex image on a two-band image
def run_cmplx(fn, image):
    if image.format == pyvips.BandFormat.FLOAT:
        new_format = pyvips.BandFormat.COMPLEX
    elif image.format == pyvips.BandFormat.DOUBLE:
        new_format = pyvips.BandFormat.DPCOMPLEX
    else:
        raise pyvips.Error("run_cmplx: not float or double")

    # tag as complex, run, revert tagging
    cmplx = image.copy(bands=1, format=new_format)
    cmplx_result = fn(cmplx)

    return cmplx_result.copy(bands=2, format=image.format)


def to_polar(image):
    """Transform image coordinates to polar.

    The image is transformed so that it is wrapped around a point in the
    centre. Vertical straight lines become circles or segments of circles,
    horizontal straight lines become radial spokes.
    """
    # xy image, zero in the centre, scaled to fit image to a circle
    xy = pyvips.Image.xyz(image.width, image.height)
    xy -= [image.width / 2.0, image.height / 2.0]
    scale = min(image.width, image.height) / float(image.width)
    xy *= 2.0 / scale

    # to polar, scale vertical axis to 360 degrees
    index = run_cmplx(lambda x: x.polar(), xy)
    index *= [1, image.height / 360.0]

    return image.mapim(index)


def to_rectangular(image):
    """Transform image coordinates to rectangular.

    The image is transformed so that it is unwrapped from a point in the
    centre. Circles or segments of circles become vertical straight lines,
    radial lines become horizontal lines.
    """
    # xy image, vertical scaled to 360 degrees
    xy = pyvips.Image.xyz(image.width, image.height)
    xy *= [1, 360.0 / image.height]

    # to rect, scale to image rect
    index = run_cmplx(lambda x: x.rect(), xy)
    scale = min(image.width, image.height) / float(image.width)
    index *= scale / 2.0
    index += [image.width / 2.0, image.height / 2.0]

    return image.mapim(index)


class TestResample:
    def test_affine(self):
        im = pyvips.Image.new_from_file(JPEG_FILE)

        # vsqbs is non-interpolatory, don't test this way
        for name in ["nearest", "bicubic", "bilinear", "nohalo", "lbb"]:
            x = im
            interpolate = pyvips.Interpolate.new(name)
            for i in range(4):
                x = x.affine([0, 1, 1, 0], interpolate=interpolate)

            assert (x - im).abs().max() == 0

    def test_affine_nohalo(self):
        im = pyvips.Image.new_from_file(LOGO3_PNG_FILE)
        new_width = 328.0
        new_height = 98.0

        in_width = im.width
        factor = in_width / new_width
        shrink = math.floor(factor)
        if shrink < 1:
            shrink = 1

        residual = shrink / factor
        print('factor=%s' % factor)
        print('shrink=%s' % shrink)
        print('residual=%s' % residual)

        if shrink > 1:
            print('shrinking')
            im = im.shrink(shrink, shrink)

            shrunk_width = int(im.width)
            shrunk_height = int(im.height)

            residualx = new_width / shrunk_width
            residualy = new_height / shrunk_height

            residual = min(residualx, residualy)
        print('after shrink residual=%s' % residual)

        for interp in ["nearest", "bicubic", "bilinear", "nohalo", "lbb"]:
            if residual:
                interpolate = pyvips.Interpolate.new(interp)
                im = im.affine([residual, 0, 0, residual], interpolate=interpolate)
                im.write_to_file('%s-affined-%s.png' % (LOGO3_PNG_FILE, interp))

            sharp = im.sharpen(mode='rgb', sigma=0.66, m2=1.0, x1=1.0)
            sharp.write_to_file('%s-affined-%s-sharpened-rgb.png' % (LOGO3_PNG_FILE, interp))

    def test_reduce(self):
        im = pyvips.Image.new_from_file(JPEG_FILE)
        # cast down to 0-127, the smallest range, so we aren't messed up by
        # clipping
        im = im.cast(pyvips.BandFormat.CHAR)

        for fac in [1, 1.1, 1.5, 1.999]:
            for fmt in all_formats:
                for kernel in ["nearest", "linear",
                               "cubic", "lanczos2", "lanczos3"]:
                    x = im.cast(fmt)
                    r = x.reduce(fac, fac, kernel=kernel)
                    d = abs(r.avg() - im.avg())
                    assert d < 2

        # try constant images ... should not change the constant
        for const in [0, 1, 2, 254, 255]:
            im = (pyvips.Image.black(10, 10) + const).cast("uchar")
            for kernel in ["nearest", "linear",
                           "cubic", "lanczos2", "lanczos3"]:
                # print "testing kernel =", kernel
                # print "testing const =", const
                shr = im.reduce(2, 2, kernel=kernel)
                d = abs(shr.avg() - im.avg())
                assert d == 0

    def test_resize(self):
        im = pyvips.Image.new_from_file(JPEG_FILE)
        im2 = im.resize(0.25)
        # in py3, round() does not round to nearest in the obvious way, so we
        # have to do it by hand
        assert im2.width == int(im.width / 4.0 + 0.5)
        assert im2.height == int(im.height / 4.0 + 0.5)

        # test geometry rounding corner case
        im = pyvips.Image.black(100, 1)
        x = im.resize(0.5)
        assert x.width == 50
        assert x.height == 1

    def test_resize_logo3__lanczos3(self):
        self.resize_and_sharpen(IMAGES + '/logo3.png', 328.0)

    def test_resize_and_sharpen_zetta(self):
        self.resize_and_sharpen(IMAGES + '/zetta.png', 436.0)

    def test_resize_and_sharpen_tiny(self):
        self.resize_and_sharpen(IMAGES + '/4x4.png', 3.0)

    def test_resize_and_sharpen_two_strip(self):
        self.resize_and_sharpen(IMAGES + '/two-strip.png', 16.0)

    @staticmethod
    def resize_and_sharpen(filename, new_width):
        im = pyvips.Image.new_from_file(filename)
        # im = im.colourspace('scrgb')
        # im = im.premultiply()
        # im = im.unpremultiply()

        kernel = 'lanczos3'
        # kernel = 'mitchell'
        # im = im.reduce(1 / (328.0 / 2382.0), 1 / (328.0 / 2382.0), kernel=kernel)
        print('new_width / im.width=', new_width / im.width)
        im = im.resize(new_width / im.width, kernel=kernel)
        im.write_to_file('%s.resized-lanczos.png' % filename)
        # im = im.thumbnail_image(328, linear=True)
        # im.write_to_file('%s.thumbnail-linear.png' % filename)
        im = im.sharpen(mode='rgb', sigma=0.66, m2=1.0, x1=1.0)
        im.write_to_file('%s.resized-lanczos-sharpened.png' % filename)

    def test_thumbnail_logo3(self):
        filename = IMAGES + '/logo3.png'
        im = pyvips.Image.new_from_file(filename)
        im = im.thumbnail_image(328, linear=True)
        im.write_to_file('%s.thumbnail-linear.png' % filename)
        im = im.sharpen(mode='rgb', sigma=0.66, m2=1.0, x1=1.0)
        im.write_to_file('%s.thumbnail-linear-sharpened.png' % filename)

    def test_resize_logo3_solid_bg__lanczos3(self):
        im = pyvips.Image.new_from_file(IMAGES + '/logo3-solid-bg.png')
        im = im.premultiply()
        kernel = 'lanczos3'
        # kernel = 'mitchell'
        im = im.resize(328.0 / 2382.0, kernel=kernel)
        im = im.unpremultiply()
        im = im.sharpen(mode='rgb', sigma=1 + 0.66 / 2, m2=1.0, x1=1.0)
        im.write_to_file('%s.resized-lanczos.png' % LOGO3_PNG_FILE)

    def test_shrink(self):
        im = pyvips.Image.new_from_file(JPEG_FILE)
        im2 = im.shrink(4, 4)
        # in py3, round() does not round to nearest in the obvious way, so we
        # have to do it by hand
        assert im2.width == int(im.width / 4.0 + 0.5)
        assert im2.height == int(im.height / 4.0 + 0.5)
        assert abs(im.avg() - im2.avg()) < 1

        im2 = im.shrink(2.5, 2.5)
        assert im2.width == int(im.width / 2.5 + 0.5)
        assert im2.height == int(im.height / 2.5 + 0.5)
        assert abs(im.avg() - im2.avg()) < 1

    @pytest.mark.skipif(not pyvips.at_least_libvips(8, 5),
                        reason="requires libvips >= 8.5")
    def test_thumbnail(self):
        im = pyvips.Image.thumbnail(JPEG_FILE, 100)

        assert im.height == 100
        assert im.bands == 3
        assert im.bands == 3

        # the average shouldn't move too much
        im_orig = pyvips.Image.new_from_file(JPEG_FILE)
        assert abs(im_orig.avg() - im.avg()) < 1

        # make sure we always get the right width
        for height in range(440, 1, -13):
            im = pyvips.Image.thumbnail(JPEG_FILE, height)
            assert im.height == height

        # should fit one of width or height
        im = pyvips.Image.thumbnail(JPEG_FILE, 100, height=300)
        assert im.width == 100
        assert im.height != 300
        im = pyvips.Image.thumbnail(JPEG_FILE, 300, height=100)
        assert im.width != 300
        assert im.height == 100

        # with @crop, should fit both width and height
        im = pyvips.Image.thumbnail(JPEG_FILE, 100,
                                    height=300, crop=True)
        assert im.width == 100
        assert im.height == 300

        im1 = pyvips.Image.thumbnail(JPEG_FILE, 100)
        with open(JPEG_FILE, 'rb') as f:
            buf = f.read()
        im2 = pyvips.Image.thumbnail_buffer(buf, 100)
        assert abs(im1.avg() - im2.avg()) < 1

        if have("heifload"):
            # this image is orientation 6 ... thumbnail should flip it
            im = pyvips.Image.new_from_file(HEIC_FILE)
            thumb = pyvips.Image.thumbnail(HEIC_FILE, 100)

            # original is landscape
            assert im.width > im.height

            # thumb should be portrait
            assert thumb.width < thumb.height
            assert thumb.height == 100

    def test_similarity(self):
        im = pyvips.Image.new_from_file(JPEG_FILE)
        im2 = im.similarity(angle=90)
        im3 = im.affine([0, -1, 1, 0])
        # rounding in calculating the affine transform from the angle stops
        # this being exactly true
        assert (im2 - im3).abs().max() < 50

    def test_similarity_scale(self):
        im = pyvips.Image.new_from_file(JPEG_FILE)
        im2 = im.similarity(scale=2)
        im3 = im.affine([2, 0, 0, 2])
        assert (im2 - im3).abs().max() == 0

    # added in 8.7
    def test_rotate(self):
        if have("rotate"):
            im = pyvips.Image.new_from_file(JPEG_FILE)
            im2 = im.rotate(90)
            im3 = im.affine([0, -1, 1, 0])
            # rounding in calculating the affine transform from the angle stops
            # this being exactly true
            assert (im2 - im3).abs().max() < 50

    def test_mapim(self):
        im = pyvips.Image.new_from_file(JPEG_FILE)

        p = to_polar(im)
        r = to_rectangular(p)

        # the left edge (which is squashed to the origin) will be badly
        # distorted, but the rest should not be too bad
        a = r.crop(50, 0, im.width - 50, im.height).gaussblur(2)
        b = im.crop(50, 0, im.width - 50, im.height).gaussblur(2)
        assert (a - b).abs().max() < 40

        # this was a bug at one point, strangely, if executed with debug
        # enabled
        mp = pyvips.Image.xyz(im.width, im.height)
        interp = pyvips.Interpolate.new('bicubic')
        assert im.mapim(mp, interpolate=interp).avg() == im.avg()


if __name__ == '__main__':
    pytest.main()
