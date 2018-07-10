/* Copyright (c) 2013, Robert Wang, email: robertwgh (at) gmail.com
    All rights reserved. https://github.com/robertwgh/ezSIFT

    Revision history:
      September, 15, 2013: initial version.
*/

#ifndef IMAGE_UTILITY_H
#define IMAGE_UTILITY_H

#include <list>

namespace ezsift {

struct SiftKeypoint;
struct MatchPair;

template <typename T>
class Image;

struct ImagePPM {
    int w;
    int h;
    unsigned char *img_r;
    unsigned char *img_g;
    unsigned char *img_b;
};

int read_pgm(const char *filename, unsigned char *&data, int &w, int &h);
void write_pgm(const char *filename, unsigned char *data, int w, int h);

int read_ppm(const char *filename, unsigned char *&data, int &w, int &h);
void write_ppm(const char *filename, unsigned char *data, int w, int h);
void write_rgb2ppm(const char *filename, unsigned char *r, unsigned char *g,
                   unsigned char *b, int w, int h);

void write_float_pgm(const char *filename, float *data, int w, int h, int mode);

void rasterCircle(ImagePPM *imgPPM, int r, int c, int radius);
void draw_red_circle(ImagePPM *imgPPM, int r, int c, int cR);
void draw_circle(ImagePPM *imgPPM, int r, int c, int cR, float thickness);
void draw_red_orientation(ImagePPM *imgPPM, int x, int y, float ori, int cR);

// Combine two images for interest points matching.
// Images are combined horizontally.
int combine_image(Image<unsigned char> &out_image,
                  const Image<unsigned char> &image1,
                  const Image<unsigned char> &image2);

// Draw circles to incidate the keypoints.
// Bars in the circle incidate the orientation of the keypoints.
void draw_keypoints_to_ppm_file(const char *out_filename,
                                const Image<unsigned char> &image,
                                std::list<SiftKeypoint> kpt_list);

// Draw lines between matched keypoints.
// The result image is stored in a ppm file.
int draw_match_lines_to_ppm_file(const char *filename,
                                 Image<unsigned char> &image1,
                                 Image<unsigned char> &image2,
                                 std::list<MatchPair> &match_list);

// Draw matched lines on a color RGB image.
int draw_line_to_rgb_image(const unsigned char *&data, int w, int h,
                           MatchPair &mp);

// Draw matched lines on an Image object.
int draw_line_to_image(Image<unsigned char> &image, MatchPair &mp);

/****************************************
 * Utility functions
 ***************************************/
// Image operations
// Get pixel from an image with unsigned char datatype.
inline unsigned char get_pixel(unsigned char *imageData, int w, int h, int r,
                               int c)
{
    unsigned char val;
    if (c >= 0 && c < w && r >= 0 && r < h) {
        val = imageData[r * w + c];
    }
    else if (c < 0) {
        val = imageData[r * w];
    }
    else if (c >= w) {
        val = imageData[r * w + w - 1];
    }
    else if (r < 0) {
        val = imageData[c];
    }
    else if (r >= h) {
        val = imageData[(h - 1) * w + c];
    }
    else {
        val = 0;
    }
    return val;
}

// Get pixel value from an image with float data type.
inline float get_pixel_f(float *imageData, int w, int h, int r, int c)
{
    float val;
    if (c >= 0 && c < w && r >= 0 && r < h) {
        val = imageData[r * w + c];
    }
    else if (c < 0) {
        val = imageData[r * w];
    }
    else if (c >= w) {
        val = imageData[r * w + w - 1];
    }
    else if (r < 0) {
        val = imageData[c];
    }
    else if (r >= h) {
        val = imageData[(h - 1) * w + c];
    }
    else {
        val = 0.0f;
    }
    return val;
}

} // end namespace ezsift

#endif // IMAGE_UTILITY_H
