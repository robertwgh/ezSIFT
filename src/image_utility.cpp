/*  Copyright (c) 2013, Robert Wang, email: robertwgh (at) gmail.com
   All rights reserved. https://github.com/robertwgh/ezSIFT

   Part of "img_utility.cpp" code referred to David Lowe's code
   and code from here https://sites.google.com/site/5kk73gpu2011/.

   Revision history:
      September, 15, 2013: initial version.
      July 2nd, 2018: code refactor.
*/

#include "image_utility.h"
#include "common.h"
#include "ezsift.h"
#include "image.h"

#include <cassert>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <ctype.h>
#include <limits>

namespace ezsift {

int read_pgm(const char *filename, unsigned char *&data, int &w, int &h)
{
    unsigned char *_data;
    FILE *in_file;
    char ch, type;
    int i;

    in_file = fopen(filename, "rb");
    if (!in_file) {
        fprintf(stderr, "ERROR(0): Fail to open file %s\n", filename);
        return -1;
    }
    // Determine pgm image type (only type three can be used)
    ch = getc(in_file);
    if (ch != 'P') {
        printf("ERROR(1): Not valid pgm/ppm file type\n");
        return -1;
    }
    ch = getc(in_file);
    // Convert the one digit integer currently represented as a character to
    // an integer(48 == '0')
    type = ch - 48;
    if (type != 5) {
        printf("ERROR(2): this file type (P%d) is not supported!\n", type);
        return -1;
    }
    while (getc(in_file) != '\n')
        ;                          // Skip to end of line
    while (getc(in_file) == '#') { // Skip comment lines
        while (getc(in_file) != '\n')
            ;
    }
    fseek(in_file, -1, SEEK_CUR); // Backup one character

    fscanf(in_file, "%d", &w);
    fscanf(in_file, "%d", &h);
    fscanf(in_file, "%d", &i); // Skipped here
    while (getc(in_file) != '\n')
        ;
    _data = (unsigned char *)malloc((w) * (h) * sizeof(unsigned char));

    fread(_data, sizeof(unsigned char), (w) * (h), in_file);
    data = _data;

    return 0;
}

void write_pgm(const char *filename, unsigned char *data, int w, int h)
{
    FILE *out_file;
    assert(w > 0);
    assert(h > 0);

    out_file = fopen(filename, "wb");
    if (!out_file) {
        fprintf(stderr, "Fail to open file: %s\n", filename);
        return;
    }

    fprintf(out_file, "P5\n");
    fprintf(out_file, "%d %d\n255\n", w, h);
    fwrite(data, sizeof(unsigned char), w * h, out_file);
    fclose(out_file);
}

void write_float_pgm(const char *filename, float *data, int w, int h, int mode)
{
    int i, j;
    unsigned char *charImg;
    int tmpInt;
    charImg = (unsigned char *)malloc(w * h * sizeof(unsigned char));
    for (i = 0; i < h; i++) {
        for (j = 0; j < w; j++) {
            if (mode == 1) { // clop
                if (data[i * w + j] >= 255.0) {
                    charImg[i * w + j] = 255;
                }
                else if (data[i * w + j] <= 0.0) {
                    charImg[i * w + j] = 0;
                }
                else {
                    charImg[i * w + j] = (int)data[i * w + j];
                }
            }
            else if (mode == 2) { // abs, x10, clop
                tmpInt = (int)(fabs(data[i * w + j]) * 10.0);
                if (fabs(data[i * w + j]) >= 255) {
                    charImg[i * w + j] = 255;
                }
                else if (tmpInt <= 0) {
                    charImg[i * w + j] = 0;
                }
                else {
                    charImg[i * w + j] = tmpInt;
                }
            }
            else {
                return;
            }
        }
    }
    write_pgm(filename, charImg, w, h);
    free(charImg);
}

void setPixelRed(ImagePPM *img, int r, int c)
{
    if ((r >= 0) && (r < img->h) && (c >= 0) && (c < img->w)) {
        img->img_r[r * img->w + c] = 0;
        img->img_g[r * img->w + c] = 0;
        img->img_b[r * img->w + c] = 255;
    }
}

void draw_red_circle(ImagePPM *imgPPM, int r, int c, int cR)
{
    int cx = -cR, cy = 0, err = 2 - 2 * cR; // II. Quadrant
    do {
        setPixelRed(imgPPM, r - cx, c + cy); //   I. Quadrant
        setPixelRed(imgPPM, r - cy, c - cx); //  II. Quadrant
        setPixelRed(imgPPM, r + cx, c - cy); // III. Quadrant
        setPixelRed(imgPPM, r + cy, c + cx); //  IV. Quadrant
        cR = err;
        if (cR > cx)
            err += ++cx * 2 + 1; // e_xy+e_x > 0
        if (cR <= cy)
            err += ++cy * 2 + 1; // e_xy+e_y < 0
    } while (cx < 0);
}

void draw_circle(ImagePPM *imgPPM, int r, int c, int cR, float thickness)
{
    int x, y;
    float f = thickness;
    for (x = -cR; x <= +cR; x++) // column
    {
        for (y = -cR; y <= +cR; y++) // row
        {
            if ((((x * x) + (y * y)) > (cR * cR) - (f / 2)) &&
                (((x * x) + (y * y)) < (cR * cR) + (f / 2)))
                setPixelRed(imgPPM, y + r, x + c);
        }
    }
}

// http://en.wikipedia.org/wiki/Midpoint_circle_algorithm
void rasterCircle(ImagePPM *imgPPM, int r, int c, int radius)
{
    int f = 1 - radius;
    int ddF_x = 1;
    int ddF_y = -2 * radius;
    int x = 0;
    int y = radius;

    int x0 = r;
    int y0 = c;

    setPixelRed(imgPPM, x0, y0 + radius);
    setPixelRed(imgPPM, x0, y0 - radius);
    setPixelRed(imgPPM, x0 + radius, y0);
    setPixelRed(imgPPM, x0 - radius, y0);

    while (x < y) {
        // ddF_x == 2 * x + 1;
        // ddF_y == -2 * y;
        // f == x*x + y*y - radius*radius + 2*x - y + 1;
        if (f >= 0) {
            y--;
            ddF_y += 2;
            f += ddF_y;
        }
        x++;
        ddF_x += 2;
        f += ddF_x;
        setPixelRed(imgPPM, x0 + x, y0 + y);
        setPixelRed(imgPPM, x0 - x, y0 + y);
        setPixelRed(imgPPM, x0 + x, y0 - y);
        setPixelRed(imgPPM, x0 - x, y0 - y);
        setPixelRed(imgPPM, x0 + y, y0 + x);
        setPixelRed(imgPPM, x0 - y, y0 + x);
        setPixelRed(imgPPM, x0 + y, y0 - x);
        setPixelRed(imgPPM, x0 - y, y0 - x);
    }
}

void draw_red_orientation(ImagePPM *imgPPM, int x, int y, float ori, int cR)
{
    int xe = (int)(x + cos(ori) * cR), ye = (int)(y + sin(ori) * cR);
    // Bresenham's line algorithm
    int dx = abs(xe - x), sx = x < xe ? 1 : -1;
    int dy = -abs(ye - y), sy = y < ye ? 1 : -1;
    int err = dx + dy, e2; /* error value e_xy */

    for (;;) { /* loop */
        setPixelRed(imgPPM, y, x);
        if (x == xe && y == ye)
            break;
        e2 = 2 * err;
        if (e2 >= dy) {
            err += dy;
            x += sx;
        } /* e_xy+e_x > 0 */
        if (e2 <= dx) {
            err += dx;
            y += sy;
        } /* e_xy+e_y < 0 */
    }
}

void skip_comment(FILE *fp)
{
    int c;
    while (isspace(c = getc(fp)))
        ;
    if (c != '#') {
        ungetc(c, fp);
        return;
    }

    do {
        c = getc(fp);
    } while (c != '\n' && c != EOF);
}

int read_ppm(const char *filename, unsigned char *&data, int &w, int &h)
{
    FILE *fp;

    if ((fp = fopen(filename, "rb")) == NULL) {
        printf("Could not read file: %s\n", filename);
        return -1;
    }

    int width, height, maxComp;
    char cookie[3];
    fscanf(fp, "%2s", cookie);
    if (strcmp("P6", cookie)) {
        printf("Wrong file type\n");
        fclose(fp);
        return -1;
    }
    skip_comment(fp);

    fscanf(fp, "%4d", &width);
    fscanf(fp, "%4d", &height);
    fscanf(fp, "%3d", &maxComp);
    fread(cookie, 1, 1, fp); // Read newline which follows maxval

    if (maxComp != 255) {
        printf("Data error: %d\n", maxComp);
        fclose(fp);
        return -1;
    }

    if (data == NULL) {
        data = new unsigned char[3 * width * height];
    }

    size_t res = fread(data, sizeof(unsigned char), 3 * width * height, fp);
    assert((int)res == 3 * width * height);
    fclose(fp);

    w = width;
    h = height;

    return 0;
}

void write_ppm(const char *filename, unsigned char *data, int w, int h)
{
    FILE *fp;
    if ((fp = fopen(filename, "wb")) == NULL) {
        printf("Cannot write to file %s\n", filename);
        return;
    }

    /* Write header */
    fprintf(fp, "P6\n");
    fprintf(fp, "%d %d\n", w, h);
    fprintf(fp, "255\n");

    fwrite(data, sizeof(unsigned char), w * h * 3, fp);
    fclose(fp);
}

void write_rgb2ppm(const char *filename, unsigned char *r, unsigned char *g,
                   unsigned char *b, int w, int h)
{
    FILE *out_file;
    int i;

    unsigned char *obuf =
        (unsigned char *)malloc(3 * w * h * sizeof(unsigned char));

    for (i = 0; i < w * h; i++) {
        obuf[3 * i + 0] = r[i];
        obuf[3 * i + 1] = g[i];
        obuf[3 * i + 2] = b[i];
    }
    out_file = fopen(filename, "wb");
    fprintf(out_file, "P6\n");
    fprintf(out_file, "%d %d\n255\n", w, h);
    fwrite(obuf, sizeof(unsigned char), 3 * w * h, out_file);
    fclose(out_file);
    free(obuf);
}

//////////////////////
// Helper Functions //
//////////////////////

// Combine two images horizontally
int combine_image(Image<unsigned char> &out_image,
                  const Image<unsigned char> &image1,
                  const Image<unsigned char> &image2)
{
    int w1 = image1.w;
    int h1 = image1.h;
    int w2 = image2.w;
    int h2 = image2.h;
    int dstW = w1 + w2;
    int dstH = MAX(h1, h2);

    out_image.init(dstW, dstH);

    unsigned char *srcData1 = image1.data;
    unsigned char *srcData2 = image2.data;
    unsigned char *dstData = out_image.data;

    for (int r = 0; r < dstH; r++) {
        if (r < h1) {
            memcpy(dstData, srcData1, w1 * sizeof(unsigned char));
        }
        else {
            memset(dstData, 0, w1 * sizeof(unsigned char));
        }
        dstData += w1;

        if (r < h2) {
            memcpy(dstData, srcData2, w2 * sizeof(unsigned char));
        }
        else {
            memset(dstData, 0, w2 * sizeof(unsigned char));
        }
        dstData += w2;
        srcData1 += w1;
        srcData2 += w2;
    }

    return 0;
}

// Helper callback function for merge match list.
bool same_match_pair(const MatchPair &first, const MatchPair &second)
{
    if (first.c1 == second.c1 && first.r1 == second.r1 &&
        first.c2 == second.c2 && first.r2 == second.r2)
        return true;
    else
        return false;
}

// Match keypoints from two images, using brutal force method.
// Use Euclidean distance as matching score.
int match_keypoints(std::list<SiftKeypoint> &kpt_list1,
                    std::list<SiftKeypoint> &kpt_list2,
                    std::list<MatchPair> &match_list)
{
    std::list<SiftKeypoint>::iterator kpt1;
    std::list<SiftKeypoint>::iterator kpt2;

    for (kpt1 = kpt_list1.begin(); kpt1 != kpt_list1.end(); kpt1++) {
        // Position of the matched feature.
        int r1 = (int)kpt1->r;
        int c1 = (int)kpt1->c;

        float *descr1 = kpt1->descriptors;
        float score1 = (std::numeric_limits<float>::max)(); // highest score
        float score2 = (std::numeric_limits<float>::max)(); // 2nd highest score

        // Position of the matched feature.
        int r2 = 0, c2 = 0;
        for (kpt2 = kpt_list2.begin(); kpt2 != kpt_list2.end(); kpt2++) {
            float score = 0;
            float *descr2 = kpt2->descriptors;
            float dif;
            for (int i = 0; i < DEGREE_OF_DESCRIPTORS; i++) {
                dif = descr1[i] - descr2[i];
                score += dif * dif;
            }

            if (score < score1) {
                score2 = score1;
                score1 = score;
                r2 = (int)kpt2->r;
                c2 = (int)kpt2->c;
            }
            else if (score < score2) {
                score2 = score;
            }
        }

#if (USE_FAST_FUNC == 1)
        if (fast_sqrt_f(score1 / score2) < SIFT_MATCH_NNDR_THR)
#else
        if (sqrtf(score1 / score2) < SIFT_MATCH_NNDR_THR)
#endif
        {
            MatchPair mp;
            mp.r1 = r1;
            mp.c1 = c1;
            mp.r2 = r2;
            mp.c2 = c2;

            match_list.push_back(mp);
        }
    }

    match_list.unique(same_match_pair);

#if PRINT_MATCH_KEYPOINTS
    std::list<MatchPair>::iterator p;
    int match_idx = 0;
    for (p = match_list.begin(); p != match_list.end(); p++) {
        printf("\tMatch %3d: (%4d, %4d) -> (%4d, %4d)\n", match_idx, p->r1,
               p->c1, p->r2, p->c2);
        match_idx++;
    }
#endif

    return 0;
}

void draw_keypoints_to_ppm_file(const char *out_filename,
                                const Image<unsigned char> &image,
                                std::list<SiftKeypoint> kpt_list)
{
    std::list<SiftKeypoint>::iterator it;
    ImagePPM imgPPM;
    int w = image.w;
    int h = image.h;
    int r, c;

    /*******************************
     * cR:
     * radius of the circle
     * cR = sigma * 4 * (2^O)
     *******************************/
    int cR;

    // initialize the imgPPM
    imgPPM.w = w;
    imgPPM.h = h;
    imgPPM.img_r = new unsigned char[w * h];
    imgPPM.img_g = new unsigned char[w * h];
    imgPPM.img_b = new unsigned char[w * h];

    int i, j;
    unsigned char *data = image.data;
    // Copy gray PGM images to color PPM images
    for (i = 0; i < h; i++) {
        for (j = 0; j < w; j++) {
            imgPPM.img_r[i * w + j] = data[i * w + j];
            imgPPM.img_g[i * w + j] = data[i * w + j];
            imgPPM.img_b[i * w + j] = data[i * w + j];
        }
    }

    for (it = kpt_list.begin(); it != kpt_list.end(); it++) {
        // derive circle radius cR
        cR = (int)it->scale;
        if (cR <= 1) { // avoid zero radius
            cR = 1;
        }
        r = (int)it->r;
        c = (int)it->c;
        //  draw_red_circle(&imgPPM, r, c, cR);
        rasterCircle(&imgPPM, r, c, cR);
        rasterCircle(&imgPPM, r, c, cR + 1);
        float ori = it->ori;
        draw_red_orientation(&imgPPM, c, r, ori, cR);
    }

    // write rendered image to output
    write_rgb2ppm(out_filename, imgPPM.img_r, imgPPM.img_g, imgPPM.img_b, w, h);

    // free allocated memory
    delete[] imgPPM.img_r;
    delete[] imgPPM.img_g;
    delete[] imgPPM.img_b;
    imgPPM.img_r = imgPPM.img_g = imgPPM.img_b = nullptr;
}

int export_kpt_list_to_file(const char *filename,
                            std::list<SiftKeypoint> &kpt_list,
                            bool bIncludeDescpritor)
{
    FILE *fp;
    fp = fopen(filename, "wb");
    if (!fp) {
        printf("Fail to open file: %s\n", filename);
        return -1;
    }

    fprintf(fp, "%u\t%d\n", static_cast<unsigned int>(kpt_list.size()), 128);

    std::list<SiftKeypoint>::iterator it;
    for (it = kpt_list.begin(); it != kpt_list.end(); it++) {
        fprintf(fp, "%d\t%d\t%.2f\t%.2f\t%.2f\t%.2f\t", it->octave, it->layer,
                it->r, it->c, it->scale, it->ori);
        if (bIncludeDescpritor) {
            for (int i = 0; i < 128; i++) {
                fprintf(fp, "%d\t", (int)(it->descriptors[i]));
            }
        }
        fprintf(fp, "\n");
    }

    fclose(fp);
    return 0;
}

// Draw a while line on a gray-scale image.
// MatchPair mp indicates the coordinates of the start point
// and the end point.
int draw_line_to_image(Image<unsigned char> &image, MatchPair &mp)
{
    int w = image.w;
    int r1 = mp.r1;
    int c1 = mp.c1;
    int r2 = mp.r2;
    int c2 = mp.c2 + w / 2;

    float k = (float)(r2 - r1) / (float)(c2 - c1);
    for (int c = c1; c < c2; c++) {
        // Line equation
        int r = (int)(k * (c - c1) + r1);
        image.data[r * w + c] = 190; // Draw a gray line.
    }
    return 0;
}

// Draw a line on the RGB color image.
int draw_line_to_rgb_image(unsigned char *&data, int w, int h, MatchPair &mp)
{
    int r1 = mp.r1;
    int c1 = mp.c1;
    int r2 = mp.r2;
    int c2 = mp.c2;

    float k = (float)(r2 - r1) / (float)(c2 - c1);
    for (int c = c1; c < c2; c++) {
        // Line equation
        int r = (int)(k * (c - c1) + r1);

        // Draw a blue line
        data[r * w * 3 + 3 * c] = 0;
        data[r * w * 3 + 3 * c + 1] = 0;
        data[r * w * 3 + 3 * c + 2] = 255;
    }

    return 0;
}

// Draw match lines between matched keypoints between two images.
int draw_match_lines_to_ppm_file(const char *filename,
                                 Image<unsigned char> &image1,
                                 Image<unsigned char> &image2,
                                 std::list<MatchPair> &match_list)
{
    Image<unsigned char> tmpImage;
    combine_image(tmpImage, image1, image2);

    int w = tmpImage.w;
    int h = tmpImage.h;
    unsigned char *srcData = tmpImage.data;
    unsigned char *dstData = new unsigned char[w * h * 3];

    for (int i = 0; i < w * h; i++) {
        dstData[i * 3 + 0] = srcData[i];
        dstData[i * 3 + 1] = srcData[i];
        dstData[i * 3 + 2] = srcData[i];
    }

    std::list<MatchPair>::iterator p;
    for (p = match_list.begin(); p != match_list.end(); p++) {
        MatchPair mp;
        mp.r1 = p->r1;
        mp.c1 = p->c1;
        mp.r2 = p->r2;
        mp.c2 = p->c2 + image1.w;
        draw_line_to_rgb_image(dstData, w, h, mp);
    }

    write_ppm(filename, dstData, w, h);

    delete[] dstData;
    dstData = nullptr;
    return 0;
}

} // end namespace ezsift
