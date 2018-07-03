/*  Copyright (c) 2013, Robert Wang, email: robertwgh (at) gmail.com
   All rights reserved. https://github.com/robertwgh/ezSIFT

   Part of "img_utility.cpp" code referred to David Lowe's code
   and code from here https://sites.google.com/site/5kk73gpu2011/.

   Revision history:
      September, 15, 2013: initial version.
      July 2nd, 2018: code refactor.
*/

#include "image_utility.h"

#include <cassert>
#include <cmath>
#include <cstdio>
#include <ctype.h>
#include <string>

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
    /* Determine pgm image type (only type three can be used)*/
    ch = getc(in_file);
    if (ch != 'P') {
        printf("ERROR(1): Not valid pgm/ppm file type\n");
        return -1;
    }
    ch = getc(in_file);
    /* Convert the one digit integer currently represented as a character to
     * an integer(48 == '0') */
    type = ch - 48;
    if (type != 5) {
        printf("ERROR(2): this file type (P%d) is not supported!\n", type);
        return -1;
    }
    while (getc(in_file) != '\n')
        ;                          /* Skip to end of line*/
    while (getc(in_file) == '#') { /* Skip comment lines */
        while (getc(in_file) != '\n')
            ;
    }
    fseek(in_file, -1, SEEK_CUR); /* Backup one character*/

    fscanf(in_file, "%d", &w);
    fscanf(in_file, "%d", &h);
    fscanf(in_file, "%d", &i); /* Skipped here */
    while (getc(in_file) != '\n')
        ;
    _data = (unsigned char *)malloc((w) * (h) * sizeof(unsigned char));

    fread(_data, sizeof(unsigned char), (w) * (h), in_file);
    data = _data;

    return 0;
} // read_pgm()

void write_pgm(const char *filename, unsigned char *data, int w, int h)
{
    FILE *out_file;
    assert(w > 0);
    assert(h > 0);

    out_file = fopen(filename, "wb");
    if (!out_file) {
        fprintf(stderr, "Fail to open file: %s\n", filename);
        exit(1);
    }

    fprintf(out_file, "P5\n");
    fprintf(out_file, "%d %d\n255\n", w, h);
    fwrite(data, sizeof(unsigned char), w * h, out_file);
    fclose(out_file);
} // write_pgm ()

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
} // write_float_pgm()

void setPixelRed(ImagePPM *img, int r, int c)
{
    if ((r >= 0) && (r < img->h) && (c >= 0) && (c < img->w)) {
        img->img_r[r * img->w + c] = 0;
        img->img_g[r * img->w + c] = 0;
        img->img_b[r * img->w + c] = 255;
    }
} // setPixelRed()

void draw_red_circle(ImagePPM *imgPPM, int r, int c, int cR)
{
    int cx = -cR, cy = 0, err = 2 - 2 * cR; /* II. Quadrant */
    do {
        setPixelRed(imgPPM, r - cx, c + cy); /*   I. Quadrant */
        setPixelRed(imgPPM, r - cy, c - cx); /*  II. Quadrant */
        setPixelRed(imgPPM, r + cx, c - cy); /* III. Quadrant */
        setPixelRed(imgPPM, r + cy, c + cx); /*  IV. Quadrant */
        cR = err;
        if (cR > cx)
            err += ++cx * 2 + 1; /* e_xy+e_x > 0 */
        if (cR <= cy)
            err += ++cy * 2 + 1; /* e_xy+e_y < 0 */
    } while (cx < 0);
} // draw_red_circle()

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
        // setPixelRed(imgPPM, x, y);
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
} // draw_red_orientation()

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
    } while (c != '\n' && c != EOF); /* gobble comment */
} // skip_comment()

int read_ppm(const char *filename, unsigned char *&data, int &w, int &h)
{
    // printf("Reading PPM file: %s\n", aFilename);
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
    fread(cookie, 1, 1, fp); /* Read newline which follows maxval */

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

} // end namespace ezsift
