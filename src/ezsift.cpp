/*  Copyright (c) 2013, Robert Wang, email: robertwgh (at) gmail.com
   All rights reserved. https://github.com/robertwgh/ezSIFT

   Some algorithms used in this code referred to:
   1. OpenCV: http://opencv.org/
   2. VLFeat: http://www.vlfeat.org/

   The SIFT algorithm was developed by David Lowe. More information can be found
   from: David G. Lowe, "Distinctive image features from scale-invariant
   keypoints," International Journal of Computer Vision, 60, 2 (2004), pp.
   91-110.

   Pay attention that the SIFT algorithm is patented. It is your responsibility
   to use the code in a legal way. Patent information: Method and apparatus for
   identifying scale invariant features in an image and use of same for locating
   an object in an image  David G. Lowe, US Patent 6,711,293 (March 23, 2004).
   Provisional application filed March 8, 1999. Asignee: The University of
   British Columbia.

   Revision history:
      September 15th, 2013: initial version.
      July 8th, 2014: fixed a bug in sample_2x in image.h. The bug only happened
   for image with odd width or height. July 2nd, 2018: code refactor.
*/

#include "ezsift.h"
#include "common.h"
#include "image.h"
#include "timer.h"
#include "vvector.h"

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <list>

namespace ezsift {

// Init sift parameters
void init_sift_parameters(bool doubleFirstOctave, float contrast_threshold,
                          float edge_threshold, float match_NDDR_threshold)
{
    SIFT_IMG_DBL = doubleFirstOctave;
    SIFT_CONTR_THR = contrast_threshold;
    SIFT_CURV_THR = edge_threshold;
    SIFT_MATCH_NNDR_THR = match_NDDR_threshold;
}

// Set up first Octave
// doubleFirstOctave = true, firstOcative=-1;
// doubleFirstOctave = false, firstOcative=0;
void double_original_image(bool doubleFirstOctave)
{
    SIFT_IMG_DBL = doubleFirstOctave;
    return;
}

// Compute octaves to build Gaussian Pyramid.
int build_octaves(const Image<unsigned char> &image,
                  std::vector<Image<unsigned char>> &octaves, int firstOctave,
                  int nOctaves)
{
    int w = image.w;
    int h = image.h;
    if (firstOctave == -1) {
        w = image.w * 2;
        h = image.h * 2;
    }

    for (int i = 0; i < nOctaves; i++) {
        if (i == 0 && firstOctave == -1) {
            octaves[i] = image.upsample_2x();
        }
        else if ((i == 0 && firstOctave != -1) ||
                 (i == 1 && firstOctave == -1)) {
            octaves[i] = image;
        }
        else {
            octaves[i] = octaves[(i - 1)].downsample_2x();
        }
        w = w / 2;
        h = h / 2;
    }
    return 0;
}

// Improved Gaussian Blurring Function
int gaussian_blur(const Image<float> &in_image, Image<float> &out_image,
                  std::vector<float> coef1d)
{
    int w = in_image.w;
    int h = in_image.h;
    int gR = static_cast<int>(coef1d.size()) / 2;

    Image<float> img_t(h, w);
    row_filter_transpose(in_image.data, img_t.data, w, h, &coef1d[0], gR);
    row_filter_transpose(img_t.data, out_image.data, h, w, &coef1d[0], gR);

    return 0;
}

// Apply Gaussian row filter to image, then transpose the image.
int row_filter_transpose(float *src, float *dst, int w, int h, float *coef1d,
                         int gR)
{
    float *row_buf = new float[w + gR * 2];
    float *row_start;
    int elemSize = sizeof(float);

    float *srcData = src;
    float *dstData = dst + w * h - 1;
    float partialSum = 0.0f;
    float *coef = coef1d;
    float *prow;

    float firstData, lastData;
    for (int r = 0; r < h; r++) {
        row_start = srcData + r * w;
        memcpy(row_buf + gR, row_start, elemSize * w);
        firstData = *(row_start);
        lastData = *(row_start + w - 1);
        for (int i = 0; i < gR; i++) {
            row_buf[i] = firstData;
            row_buf[i + w + gR] = lastData;
        }

        prow = row_buf;
        dstData = dstData - w * h + 1;
        for (int c = 0; c < w; c++) {
            partialSum = 0.0f;
            coef = coef1d;

            for (int i = -gR; i <= gR; i++) {
                partialSum += (*coef++) * (*prow++);
            }

            prow -= 2 * gR;
            *dstData = partialSum;
            dstData += h;
        }
    }
    delete[] row_buf;
    row_buf = nullptr;

    return 0;
}

// Build Gaussian pyramid using recursive method.
// The first layer is downsampled from last octave, layer=nLayers.
// All n-th layer is Gaussian blur from (n-1)-th layer.
int build_gaussian_pyramid(std::vector<Image<unsigned char>> &octaves,
                           std::vector<Image<float>> &gpyr, int nOctaves,
                           int nGpyrLayers)
{
    int nLayers = nGpyrLayers - 3;
    std::vector<std::vector<float>> gaussian_coefs =
        compute_gaussian_coefs(nOctaves, nGpyrLayers);

    int w, h;
    for (int i = 0; i < nOctaves; i++) {
        w = octaves[i].w;
        h = octaves[i].h;
        for (int j = 0; j < nGpyrLayers; j++) {
            if (i == 0 && j == 0) {
                gpyr[0].init(w, h);
                gaussian_blur(octaves[0].to_float(), gpyr[0],
                              gaussian_coefs[j]);
            }
            else if (i > 0 && j == 0) {
                gpyr[i * nGpyrLayers] =
                    gpyr[(i - 1) * nGpyrLayers + nLayers].downsample_2x();
            }
            else {
                gpyr[i * nGpyrLayers + j].init(w, h);
                gaussian_blur(gpyr[i * nGpyrLayers + j - 1],
                              gpyr[i * nGpyrLayers + j], gaussian_coefs[j]);
            }
        }
    }
    // Release octaves memory.
    octaves.clear();
    return 0;
}

// For build_gaussian_pyramid()
std::vector<std::vector<float>> compute_gaussian_coefs(int nOctaves,
                                                       int nGpyrLayers)
{
    // Compute all sigmas for different layers
    int nLayers = nGpyrLayers - 3;
    float sigma, sigma_pre;
    float sigma0 = SIFT_SIGMA;
    float k = powf(2.0f, 1.0f / nLayers);

    std::vector<float> sig(nGpyrLayers);
    sigma_pre = SIFT_IMG_DBL ? 2.0f * SIFT_INIT_SIGMA : SIFT_INIT_SIGMA;
    sig[0] = sqrtf(sigma0 * sigma0 - sigma_pre * sigma_pre);
    for (int i = 1; i < nGpyrLayers; i++) {
        sigma_pre = powf(k, (float)(i - 1)) * sigma0;
        sigma = sigma_pre * k;
        sig[i] = sqrtf(sigma * sigma - sigma_pre * sigma_pre);
    }

    std::vector<std::vector<float>> gaussian_coefs(nGpyrLayers);
    for (int i = 0; i < nGpyrLayers; i++) {
        // Compute Gaussian filter coefficients
        float factor = SIFT_GAUSSIAN_FILTER_RADIUS;
        int gR = (sig[i] * factor > 1.0f) ? (int)ceilf(sig[i] * factor) : 1;
        int gW = gR * 2 + 1;

        gaussian_coefs[i].resize(gW);
        float accu = 0.0f;
        float tmp;
        for (int j = 0; j < gW; j++) {
            tmp = (float)((j - gR) / sig[i]);
            gaussian_coefs[i][j] = expf(tmp * tmp * -0.5f) * (1 + j / 1000.0f);
            accu += gaussian_coefs[i][j];
        }
        for (int j = 0; j < gW; j++) {
            gaussian_coefs[i][j] = gaussian_coefs[i][j] / accu;
        } // End compute Gaussian filter coefs
    }
    return gaussian_coefs;
}

// Build difference of Gaussian pyramids.
int build_dog_pyr(std::vector<Image<float>> &gpyr,
                  std::vector<Image<float>> &dogPyr, int nOctaves,
                  int nDogLayers)
{
    int nGpyrLayers = nDogLayers + 1;

    int w, h;
    float *srcData1; // always data2-data1.
    float *srcData2;
    float *dstData;
    int index = 0;

    for (int i = 0; i < nOctaves; i++) {
        int row_start = i * nGpyrLayers;
        w = gpyr[row_start].w;
        h = gpyr[row_start].h;

        for (int j = 0; j < nDogLayers; j++) {
            dogPyr[i * nDogLayers + j].init(w, h);
            dstData = dogPyr[i * nDogLayers + j].data;

            srcData1 = gpyr[row_start + j].data;
            srcData2 = gpyr[row_start + j + 1].data;

            index = 0;
            while (index++ < w * h)
                *(dstData++) = *(srcData2++) - *(srcData1++);
        }
    }

    return 0;
}

// Build gradient pyramids.
int build_grd_rot_pyr(std::vector<Image<float>> &gpyr,
                      std::vector<Image<float>> &grdPyr,
                      std::vector<Image<float>> &rotPyr, int nOctaves,
                      int nLayers)
{
    int nGpyrLayers = nLayers + 3;
    int w, h;
    float dr, dc;
    float angle;

    float *srcData;
    float *grdData;
    float *rotData;

    for (int i = 0; i < nOctaves; i++) {
        // We only use gradient information from layers 1~n Layers.
        // Since keypoints only occur at these layers.

        w = gpyr[i * nGpyrLayers].w;
        h = gpyr[i * nGpyrLayers].h;
        for (int j = 1; j <= nLayers; j++) {
            int layer_index = i * nGpyrLayers + j;
            grdPyr[layer_index].init(w, h);
            rotPyr[layer_index].init(w, h);

            srcData = gpyr[layer_index].data;
            grdData = grdPyr[layer_index].data;
            rotData = rotPyr[layer_index].data;

            for (int r = 0; r < h; r++) {
                for (int c = 0; c < w; c++) {
                    dr = get_pixel_f(srcData, w, h, r + 1, c) -
                         get_pixel_f(srcData, w, h, r - 1, c);
                    dc = get_pixel_f(srcData, w, h, r, c + 1) -
                         get_pixel_f(srcData, w, h, r, c - 1);

#if (USE_FAST_FUNC == 1)
                    grdData[r * w + c] = fast_sqrt_f(dr * dr + dc * dc);
                    angle = fast_atan2_f(dr, dc); // atan2f(dr, dc + FLT_MIN);
#else
                    grdData[r * w + c] = sqrtf(dr * dr + dc * dc);
                    angle = atan2f(dr, dc + FLT_MIN);
                    angle = angle < 0 ? angle + _2PI : angle;
#endif
                    rotData[r * w + c] = angle;
                }
            }
        }
    }
    return 0;
}

// Compute orientation histogram for keypoint detection.
// Gradient information is computed in this function.
float compute_orientation_hist(const Image<float> &image, SiftKeypoint &kpt,
                               float *&hist)
{
    int nBins = SIFT_ORI_HIST_BINS;

    float kptr = kpt.ri;
    float kptc = kpt.ci;
    float kpt_scale = kpt.layer_scale;

    int kptr_i = (int)(kptr + 0.5f);
    int kptc_i = (int)(kptc + 0.5f);
    float d_kptr = kptr - kptr_i;
    float d_kptc = kptc - kptc_i;

    float sigma = SIFT_ORI_SIG_FCTR * kpt_scale;
    int win_radius = (int)(SIFT_ORI_RADIUS * kpt_scale);

    float *data = image.data;
    int w = image.w;
    int h = image.h;

    int r, c;
    float dr, dc;
    float magni, angle, weight;
    int bin;
    float fbin; // float point bin

    float *tmpHist = new float[nBins];
    memset(tmpHist, 0, nBins * sizeof(float));

    for (int i = -win_radius; i <= win_radius; i++) // rows
    {
        r = kptr_i + i;
        if (r <= 0 || r >= h - 1) // Cannot calculate dy
            continue;
        for (int j = -win_radius; j <= win_radius; j++) // columns
        {
            c = kptc_i + j;
            if (c <= 0 || c >= w - 1)
                continue;

            dr = get_pixel_f(data, w, h, r + 1, c) -
                 get_pixel_f(data, w, h, r - 1, c);
            dc = get_pixel_f(data, w, h, r, c + 1) -
                 get_pixel_f(data, w, h, r, c - 1);

#if (USE_FAST_FUNC == 1)
            magni = fast_sqrt_f(dr * dr + dc * dc);
            angle = fast_atan2_f(dr, dc); // Unit: degree, range: [-PI, PI]
#else
            magni = sqrtf(dr * dr + dc * dc);
            angle = atan2f(dr, dc); // Unit: degree, range: [-PI, PI]
            angle = angle < 0 ? angle + _2PI : angle;
#endif
            fbin = angle * nBins / _2PI;
            weight = expf(
                -1.0f *
                ((i - d_kptr) * (i - d_kptr) + (j - d_kptc) * (j - d_kptc)) /
                (2.0f * sigma * sigma));

#define SIFT_ORI_BILINEAR
#ifdef SIFT_ORI_BILINEAR
            bin = (int)(fbin - 0.5f);
            float d_fbin = fbin - 0.5f - bin;
            tmpHist[(bin + nBins) % nBins] += (1 - d_fbin) * magni * weight;
            tmpHist[(bin + 1) % nBins] += d_fbin * magni * weight;
#else
            bin = (int)(fbin);
            tmpHist[bin] += magni * weight;
#endif
        }
    }

#define TMPHIST(idx)                                                           \
    (idx < 0 ? tmpHist[0] : (idx >= nBins ? tmpHist[nBins - 1] : tmpHist[idx]))

#define USE_SMOOTH1 1
#if USE_SMOOTH1
    // Smooth the histogram. Algorithm comes from OpenCV.
    for (int i = 0; i < nBins; i++) {
        hist[i] = (TMPHIST(i - 2) + TMPHIST(i + 2)) * 1.0f / 16.0f +
                  (TMPHIST(i - 1) + TMPHIST(i + 1)) * 4.0f / 16.0f +
                  TMPHIST(i) * 6.0f / 16.0f;
    }
#else
    // Yet another smooth function
    // Smoothing algorithm comes from vl_feat implementation.
    for (int iter = 0; iter < 6; iter++) {
        float prev = TMPHIST(nBins - 1);
        float first = TMPHIST(0);
        int i;
        for (i = 0; i < nBins - 1; i++) {
            float newh = (prev + TMPHIST(i) + TMPHIST(i + 1)) / 3.0f;
            prev = hist[i];
            hist[i] = newh;
        }
        hist[i] = (prev + hist[i] + first) / 3.0f;
    }
#endif

    // Find the maximum item of the histogram
    float maxitem = hist[0];
    int max_i = 0;
    for (int i = 0; i < nBins; i++) {
        if (maxitem < hist[i]) {
            maxitem = hist[i];
            max_i = i;
        }
    }

    kpt.ori = max_i * _2PI / nBins;

    delete[] tmpHist;
    tmpHist = nullptr;
    return maxitem;
}

// Compute orientation histogram for keypoint detection.
// using pre-computed gradient information.
float compute_orientation_hist_with_gradient(const Image<float> &grdImage,
                                             const Image<float> &rotImage,
                                             SiftKeypoint &kpt, float *&hist)
{
    int nBins = SIFT_ORI_HIST_BINS;

    float kptr = kpt.ri;
    float kptc = kpt.ci;
    float kpt_scale = kpt.layer_scale;

    int kptr_i = (int)(kptr + 0.5f);
    int kptc_i = (int)(kptc + 0.5f);
    float d_kptr = kptr - kptr_i;
    float d_kptc = kptc - kptc_i;

    float sigma = SIFT_ORI_SIG_FCTR * kpt_scale;
    int win_radius = (int)(SIFT_ORI_RADIUS * kpt_scale);
    float exp_factor = -1.0f / (2.0f * sigma * sigma);

    float *grdData = grdImage.data;
    float *rotData = rotImage.data;
    int w = grdImage.w;
    int h = grdImage.h;

    int r, c;
    float magni, angle, weight;
    int bin;
    float fbin; // float point bin

    float *tmpHist = new float[nBins];
    memset(tmpHist, 0, nBins * sizeof(float));

    for (int i = -win_radius; i <= win_radius; i++) // rows
    {
        r = kptr_i + i;
        if (r <= 0 || r >= h - 1) // Cannot calculate dy
            continue;
        for (int j = -win_radius; j <= win_radius; j++) // columns
        {
            c = kptc_i + j;
            if (c <= 0 || c >= w - 1)
                continue;

            magni = grdData[r * w + c];
            angle = rotData[r * w + c];

            fbin = angle * nBins / _2PI;
            weight = expf(
                ((i - d_kptr) * (i - d_kptr) + (j - d_kptc) * (j - d_kptc)) *
                exp_factor);

#define SIFT_ORI_BILINEAR
#ifdef SIFT_ORI_BILINEAR
            bin = (int)(fbin - 0.5f);
            float d_fbin = fbin - 0.5f - bin;

            float mw = weight * magni;
            float dmw = d_fbin * mw;
            tmpHist[(bin + nBins) % nBins] += mw - dmw;
            tmpHist[(bin + 1) % nBins] += dmw;
#else
            bin = (int)(fbin);
            tmpHist[bin] += magni * weight;
#endif
        }
    }

#define TMPHIST(idx)                                                           \
    (idx < 0 ? tmpHist[0] : (idx >= nBins ? tmpHist[nBins - 1] : tmpHist[idx]))

#define USE_SMOOTH1 1
#if USE_SMOOTH1

    // Smooth the histogram. Algorithm comes from OpenCV.
    hist[0] = (tmpHist[0] + tmpHist[2]) * 1.0f / 16.0f +
              (tmpHist[0] + tmpHist[1]) * 4.0f / 16.0f +
              tmpHist[0] * 6.0f / 16.0f;
    hist[1] = (tmpHist[0] + tmpHist[3]) * 1.0f / 16.0f +
              (tmpHist[0] + tmpHist[2]) * 4.0f / 16.0f +
              tmpHist[1] * 6.0f / 16.0f;
    hist[nBins - 2] = (tmpHist[nBins - 4] + tmpHist[nBins - 1]) * 1.0f / 16.0f +
                      (tmpHist[nBins - 3] + tmpHist[nBins - 1]) * 4.0f / 16.0f +
                      tmpHist[nBins - 2] * 6.0f / 16.0f;
    hist[nBins - 1] = (tmpHist[nBins - 3] + tmpHist[nBins - 1]) * 1.0f / 16.0f +
                      (tmpHist[nBins - 2] + tmpHist[nBins - 1]) * 4.0f / 16.0f +
                      tmpHist[nBins - 1] * 6.0f / 16.0f;

    for (int i = 2; i < nBins - 2; i++) {
        hist[i] = (tmpHist[i - 2] + tmpHist[i + 2]) * 1.0f / 16.0f +
                  (tmpHist[i - 1] + tmpHist[i + 1]) * 4.0f / 16.0f +
                  tmpHist[i] * 6.0f / 16.0f;
    }

#else
    // Yet another smooth function
    // Algorithm comes from the vl_feat implementation.
    for (int iter = 0; iter < 6; iter++) {
        float prev = TMPHIST(nBins - 1);
        float first = TMPHIST(0);
        int i;
        for (i = 0; i < nBins - 1; i++) {
            float newh = (prev + TMPHIST(i) + TMPHIST(i + 1)) / 3.0f;
            prev = hist[i];
            hist[i] = newh;
        }
        hist[i] = (prev + hist[i] + first) / 3.0f;
    }
#endif

    // Find the maximum item of the histogram
    float maxitem = hist[0];
    int max_i = 0;
    for (int i = 0; i < nBins; i++) {
        if (maxitem < hist[i]) {
            maxitem = hist[i];
            max_i = i;
        }
    }

    kpt.ori = max_i * _2PI / nBins;

    delete[] tmpHist;
    tmpHist = nullptr;
    return maxitem;
}

// Keypoint detection.
int detect_keypoints(std::vector<Image<float>> &dogPyr,
                     std::vector<Image<float>> &grdPyr,
                     std::vector<Image<float>> &rotPyr, int nOctaves,
                     int nDogLayers, std::list<SiftKeypoint> &kpt_list)
{
    float *currData;
    float *lowData;
    float *highData;

    SiftKeypoint kpt;

    int w, h;
    int layer_index;
    int index;
    float val;

    int nBins = SIFT_ORI_HIST_BINS;
    float *hist = new float[nBins];
    int nGpyrLayers = nDogLayers + 1;

    // Some paper uses other thresholds, for example 3.0f for all cases
    // In Lowe's paper, |D(x)|<0.03 will be rejected.
    float contr_thr = 0.8f * SIFT_CONTR_THR;

    for (int i = 0; i < nOctaves; i++) {
        w = dogPyr[i * nDogLayers].w;
        h = dogPyr[i * nDogLayers].h;

        for (int j = 1; j < nDogLayers - 1; j++) {
            layer_index = i * nDogLayers + j;

            highData = dogPyr[layer_index + 1].data;
            currData = dogPyr[layer_index].data;
            lowData = dogPyr[layer_index - 1].data;

            for (int r = SIFT_IMG_BORDER; r < h - SIFT_IMG_BORDER; r++) {
                for (int c = SIFT_IMG_BORDER; c < w - SIFT_IMG_BORDER; c++) {
                    index = r * w + c;
                    val = currData[index];

                    bool bExtrema =
                        (val >= contr_thr && val > highData[index - w - 1] &&
                         val > highData[index - w] &&
                         val > highData[index - w + 1] &&
                         val > highData[index - 1] && val > highData[index] &&
                         val > highData[index + 1] &&
                         val > highData[index + w - 1] &&
                         val > highData[index + w] &&
                         val > highData[index + w + 1] &&
                         val > currData[index - w - 1] &&
                         val > currData[index - w] &&
                         val > currData[index - w + 1] &&
                         val > currData[index - 1] &&
                         val > currData[index + 1] &&
                         val > currData[index + w - 1] &&
                         val > currData[index + w] &&
                         val > currData[index + w + 1] &&
                         val > lowData[index - w - 1] &&
                         val > lowData[index - w] &&
                         val > lowData[index - w + 1] &&
                         val > lowData[index - 1] && val > lowData[index] &&
                         val > lowData[index + 1] &&
                         val > lowData[index + w - 1] &&
                         val > lowData[index + w] &&
                         val > lowData[index + w + 1]) || // Local min
                        (val <= -contr_thr && val < highData[index - w - 1] &&
                         val < highData[index - w] &&
                         val < highData[index - w + 1] &&
                         val < highData[index - 1] && val < highData[index] &&
                         val < highData[index + 1] &&
                         val < highData[index + w - 1] &&
                         val < highData[index + w] &&
                         val < highData[index + w + 1] &&
                         val < currData[index - w - 1] &&
                         val < currData[index - w] &&
                         val < currData[index - w + 1] &&
                         val < currData[index - 1] &&
                         val < currData[index + 1] &&
                         val < currData[index + w - 1] &&
                         val < currData[index + w] &&
                         val < currData[index + w + 1] &&
                         val < lowData[index - w - 1] &&
                         val < lowData[index - w] &&
                         val < lowData[index - w + 1] &&
                         val < lowData[index - 1] && val < lowData[index] &&
                         val < lowData[index + 1] &&
                         val < lowData[index + w - 1] &&
                         val < lowData[index + w] &&
                         val < lowData[index + w + 1]);

                    if (bExtrema) {
                        kpt.octave = i;
                        kpt.layer = j;
                        kpt.ri = (float)r;
                        kpt.ci = (float)c;

                        bool bGoodKeypoint = refine_local_extrema(
                            dogPyr, nOctaves, nDogLayers, kpt);

                        if (!bGoodKeypoint)
                            continue;

                        float max_mag = compute_orientation_hist_with_gradient(
                            grdPyr[i * nGpyrLayers + kpt.layer],
                            rotPyr[i * nGpyrLayers + kpt.layer], kpt, hist);

                        float threshold = max_mag * SIFT_ORI_PEAK_RATIO;

                        for (int ii = 0; ii < nBins; ii++) {
#define INTERPOLATE_ORI_HIST

#ifndef INTERPOLATE_ORI_HIST
                            if (hist[ii] >= threshold) {
                                kpt.mag = hist[ii];
                                kpt.ori = ii * _2PI / nBins;
                                kpt_list.push_back(kpt);
                            }
#else
                            // Use 3 points to fit a curve and find the accurate
                            // location of a keypoints
                            int left = ii > 0 ? ii - 1 : nBins - 1;
                            int right = ii < (nBins - 1) ? ii + 1 : 0;
                            float currHist = hist[ii];
                            float lhist = hist[left];
                            float rhist = hist[right];
                            if (currHist > lhist && currHist > rhist &&
                                currHist > threshold) {
                                // Refer to here:
                                // http://stackoverflow.com/questions/717762/how-to-calculate-the-vertex-of-a-parabola-given-three-points
                                float accu_ii =
                                    ii + 0.5f * (lhist - rhist) /
                                             (lhist - 2.0f * currHist + rhist);

                                // Since bin index means the starting point of a
                                // bin, so the real orientation should be bin
                                // index plus 0.5. for example, angles in bin 0
                                // should have a mean value of 5 instead of 0;
                                accu_ii += 0.5f;
                                accu_ii = accu_ii < 0 ? (accu_ii + nBins)
                                                      : accu_ii >= nBins
                                                            ? (accu_ii - nBins)
                                                            : accu_ii;
                                // The magnitude should also calculate the max
                                // number based on fitting But since we didn't
                                // actually use it in image matching, we just
                                // lazily use the histogram value.
                                kpt.mag = currHist;
                                kpt.ori = accu_ii * _2PI / nBins;
                                kpt_list.push_back(kpt);
                            }
#endif
                        }
                    }
                }
            }
        }
    }

    delete[] hist;
    hist = nullptr;
    return 0;
}

// Refine local keypoint extrema.
bool refine_local_extrema(std::vector<Image<float>> &dogPyr, int nOctaves,
                          int nDogLayers, SiftKeypoint &kpt)
{
    int nGpyrLayers = nDogLayers + 1;

    int w = 0, h = 0;
    int layer_idx = 0;
    int octave = kpt.octave;
    int layer = kpt.layer;
    int r = (int)kpt.ri;
    int c = (int)kpt.ci;

    float *currData = nullptr;
    float *lowData = nullptr;
    float *highData = nullptr;

    int xs_i = 0, xr_i = 0, xc_i = 0;
    float tmp_r = 0.0f, tmp_c = 0.0f, tmp_layer = 0.0f;
    float xr = 0.0f, xc = 0.0f, xs = 0.0f;
    float x_hat[3] = {xc, xr, xs};
    float dx = 0.0f, dy = 0.0f, ds = 0.0f;
    float dxx = 0.0f, dyy = 0.0f, dss = 0.0f, dxs = 0.0f, dys = 0.0f,
          dxy = 0.0f;

    tmp_r = (float)r;
    tmp_c = (float)c;
    tmp_layer = (float)layer;

    // Interpolation (x,y,sigma) 3D space to find sub-pixel accurate
    // location of keypoints.
    int i = 0;
    for (; i < SIFT_MAX_INTERP_STEPS; i++) {
        c += xc_i;
        r += xr_i;

        layer_idx = octave * nDogLayers + layer;
        w = dogPyr[layer_idx].w;
        h = dogPyr[layer_idx].h;
        currData = dogPyr[layer_idx].data;
        lowData = dogPyr[layer_idx - 1].data;
        highData = dogPyr[layer_idx + 1].data;

        dx = (get_pixel_f(currData, w, h, r, c + 1) -
              get_pixel_f(currData, w, h, r, c - 1)) *
             0.5f;
        dy = (get_pixel_f(currData, w, h, r + 1, c) -
              get_pixel_f(currData, w, h, r - 1, c)) *
             0.5f;
        ds = (get_pixel_f(highData, w, h, r, c) -
              get_pixel_f(lowData, w, h, r, c)) *
             0.5f;
        float dD[3] = {-dx, -dy, -ds};

        float v2 = 2.0f * get_pixel_f(currData, w, h, r, c);
        dxx = (get_pixel_f(currData, w, h, r, c + 1) +
               get_pixel_f(currData, w, h, r, c - 1) - v2);
        dyy = (get_pixel_f(currData, w, h, r + 1, c) +
               get_pixel_f(currData, w, h, r - 1, c) - v2);
        dss = (get_pixel_f(highData, w, h, r, c) +
               get_pixel_f(lowData, w, h, r, c) - v2);
        dxy = (get_pixel_f(currData, w, h, r + 1, c + 1) -
               get_pixel_f(currData, w, h, r + 1, c - 1) -
               get_pixel_f(currData, w, h, r - 1, c + 1) +
               get_pixel_f(currData, w, h, r - 1, c - 1)) *
              0.25f;
        dxs = (get_pixel_f(highData, w, h, r, c + 1) -
               get_pixel_f(highData, w, h, r, c - 1) -
               get_pixel_f(lowData, w, h, r, c + 1) +
               get_pixel_f(lowData, w, h, r, c - 1)) *
              0.25f;
        dys = (get_pixel_f(highData, w, h, r + 1, c) -
               get_pixel_f(highData, w, h, r - 1, c) -
               get_pixel_f(lowData, w, h, r + 1, c) +
               get_pixel_f(lowData, w, h, r - 1, c)) *
              0.25f;

        // The scale in two sides of the equation should cancel each other.
        float H[3][3] = {{dxx, dxy, dxs}, {dxy, dyy, dys}, {dxs, dys, dss}};
        float Hinvert[3][3];
        float det;

        // Matrix inversion
        // INVERT_3X3 = DETERMINANT_3X3, then SCALE_ADJOINT_3X3;
        // Using INVERT_3X3(Hinvert, det, H) is more convenient;
        // but using separate ones, we can check det==0 easily.
        float tmp;
        DETERMINANT_3X3(det, H);
        if (fabsf(det) < (std::numeric_limits<float>::min)())
            break;
        tmp = 1.0f / (det);
        // INVERT_3X3(Hinvert, det, H);
        SCALE_ADJOINT_3X3(Hinvert, tmp, H);
        MAT_DOT_VEC_3X3(x_hat, Hinvert, dD);

        xs = x_hat[2];
        xr = x_hat[1];
        xc = x_hat[0];

        // Update tmp data for keypoint update.
        tmp_r = r + xr;
        tmp_c = c + xc;
        tmp_layer = layer + xs;

        // Make sure there is room to move for next iteration.
        xc_i = ((xc >= SIFT_KEYPOINT_SUBPiXEL_THR && c < w - 2) ? 1 : 0) +
               ((xc <= -SIFT_KEYPOINT_SUBPiXEL_THR && c > 1) ? -1 : 0);

        xr_i = ((xr >= SIFT_KEYPOINT_SUBPiXEL_THR && r < h - 2) ? 1 : 0) +
               ((xr <= -SIFT_KEYPOINT_SUBPiXEL_THR && r > 1) ? -1 : 0);

        if (xc_i == 0 && xr_i == 0 && xs_i == 0)
            break;
    }

    // We MIGHT be able to remove the following two checking conditions.
    // Condition 1
    if (i >= SIFT_MAX_INTERP_STEPS)
        return false;
    // Condition 2.
    if (fabsf(xc) >= 1.5 || fabsf(xr) >= 1.5 || fabsf(xs) >= 1.5)
        return false;

    // If (r, c, layer) is out of range, return false.
    if (tmp_layer < 0 || tmp_layer > (nGpyrLayers - 1) || tmp_r < 0 ||
        tmp_r > h - 1 || tmp_c < 0 || tmp_c > w - 1)
        return false;

    {
        float value = get_pixel_f(currData, w, h, r, c) +
                      0.5f * (dx * xc + dy * xr + ds * xs);
        if (fabsf(value) < SIFT_CONTR_THR)
            return false;

        float trH = dxx + dyy;
        float detH = dxx * dyy - dxy * dxy;
        float response =
            (SIFT_CURV_THR + 1) * (SIFT_CURV_THR + 1) / (SIFT_CURV_THR);

        if (detH <= 0 || (trH * trH / detH) >= response)
            return false;
    }

    // Coordinates in the current layer.
    kpt.ci = tmp_c;
    kpt.ri = tmp_r;
    kpt.layer_scale = SIFT_SIGMA * powf(2.0f, tmp_layer / SIFT_INTVLS);

    int firstOctave = SIFT_IMG_DBL ? -1 : 0;
    float norm = powf(2.0f, (float)(octave + firstOctave));
    // Coordinates in the normalized format (compared to the original image).
    kpt.c = tmp_c * norm;
    kpt.r = tmp_r * norm;
    kpt.rlayer = tmp_layer;
    kpt.layer = layer;

    // Formula: Scale = sigma0 * 2^octave * 2^(layer/S);
    kpt.scale = kpt.layer_scale * norm;

    return true;
}

// Extract descriptor
// 1. Unroll the tri-linear part.
int extract_descriptor(std::vector<Image<float>> &grdPyr,
                       std::vector<Image<float>> &rotPyr, int nOctaves,
                       int nGpyrLayers, std::list<SiftKeypoint> &kpt_list)
{
    // Number of subregions, default 4x4 subregions.
    // The width of subregion is determined by the scale of the keypoint.
    // Or, in Lowe's SIFT paper[2004], width of subregion is 16x16.
    int nSubregion = SIFT_DESCR_WIDTH;
    int nHalfSubregion = nSubregion >> 1;

    // Number of histogram bins for each descriptor subregion.
    int nBinsPerSubregion = SIFT_DESCR_HIST_BINS;
    float nBinsPerSubregionPerDegree = (float)nBinsPerSubregion / _2PI;

    // 3-D structure for histogram bins (rbin, cbin, obin);
    // (rbin, cbin, obin) means (row of hist bin, column of hist bin,
    // orientation bin) In Lowe's paper, 4x4 histogram, each has 8 bins. that
    // means for each (rbin, cbin), there are 8 bins in the histogram.

    // In this implementation, histBin is a circular buffer.
    // we expand the cube by 1 for each direction.
    int nBins = nSubregion * nSubregion * nBinsPerSubregion;
    int nHistBins =
        (nSubregion + 2) * (nSubregion + 2) * (nBinsPerSubregion + 2);
    int nSliceStep = (nSubregion + 2) * (nBinsPerSubregion + 2);
    int nRowStep = (nBinsPerSubregion + 2);
    float *histBin = new float[nHistBins];

    float exp_scale = -2.0f / (nSubregion * nSubregion);

    for (std::list<SiftKeypoint>::iterator kpt = kpt_list.begin();
         kpt != kpt_list.end(); kpt++) {
        // Keypoint information
        int octave = kpt->octave;
        int layer = kpt->layer;

        float kpt_ori = kpt->ori;
        float kptr = kpt->ri;
        float kptc = kpt->ci;
        float kpt_scale = kpt->layer_scale;

        // Nearest coordinate of keypoints
        int kptr_i = (int)(kptr + 0.5f);
        int kptc_i = (int)(kptc + 0.5f);
        float d_kptr = kptr_i - kptr;
        float d_kptc = kptc_i - kptc;

        int layer_index = octave * nGpyrLayers + layer;
        int w = grdPyr[layer_index].w;
        int h = grdPyr[layer_index].h;
        float *grdData = grdPyr[layer_index].data;
        float *rotData = rotPyr[layer_index].data;

        // Note for Gaussian weighting.
        // OpenCV and vl_feat uses non-fixed size of subregion.
        // But they all use (0.5 * 4) as the Gaussian weighting sigma.
        // In Lowe's paper, he uses 16x16 sample region,
        // partition 16x16 region into 16 4x4 subregion.
        float subregion_width = SIFT_DESCR_SCL_FCTR * kpt_scale;
        int win_size =
            (int)(SQRT2 * subregion_width * (nSubregion + 1) * 0.5f + 0.5f);

        // Normalized cos() and sin() value.
        float sin_t = sinf(kpt_ori) / (float)subregion_width;
        float cos_t = cosf(kpt_ori) / (float)subregion_width;

        // Re-init histBin
        memset(histBin, 0, nHistBins * sizeof(float));

        // Start to calculate the histogram in the sample region.
        float rr, cc;
        float mag, angle, gaussian_weight;

        // Used for tri-linear interpolation.
        // int rbin_i, cbin_i, obin_i;
        float rrotate, crotate;
        float rbin, cbin, obin;
        float d_rbin, d_cbin, d_obin;

        // Boundary of sample region.
        int r, c;
        int left = MAX(-win_size, 1 - kptc_i);
        int right = MIN(win_size, w - 2 - kptc_i);
        int top = MAX(-win_size, 1 - kptr_i);
        int bottom = MIN(win_size, h - 2 - kptr_i);

        for (int i = top; i <= bottom; i++) // rows
        {
            for (int j = left; j <= right; j++) // columns
            {
                // Accurate position relative to (kptr, kptc)
                rr = i + d_kptr;
                cc = j + d_kptc;

                // Rotate the coordinate of (i, j)
                rrotate = (cos_t * cc + sin_t * rr);
                crotate = (-sin_t * cc + cos_t * rr);

                // Since for a bin array with 4x4 bins, the center is actually
                // at (1.5, 1.5)
                rbin = rrotate + nHalfSubregion - 0.5f;
                cbin = crotate + nHalfSubregion - 0.5f;

                // rbin, cbin range is (-1, d); if outside this range, then the
                // pixel is counted.
                if (rbin <= -1 || rbin >= nSubregion || cbin <= -1 ||
                    cbin >= nSubregion)
                    continue;

                // All the data need for gradient computation are valid, no
                // border issues.
                r = kptr_i + i;
                c = kptc_i + j;
                mag = grdData[r * w + c];
                angle = rotData[r * w + c] - kpt_ori;
                float angle1 = (angle < 0) ? (_2PI + angle)
                                           : angle; // Adjust angle to [0, 2PI)
                obin = angle1 * nBinsPerSubregionPerDegree;

                int x0, y0, z0;
                int x1, y1, z1;
                y0 = (int)floor(rbin);
                x0 = (int)floor(cbin);
                z0 = (int)floor(obin);
                d_rbin = rbin - y0;
                d_cbin = cbin - x0;
                d_obin = obin - z0;
                x1 = x0 + 1;
                y1 = y0 + 1;
                z1 = z0 + 1;

                // Gaussian weight relative to the center of sample region.
                gaussian_weight =
                    expf((rrotate * rrotate + crotate * crotate) * exp_scale);

                // Gaussian-weighted magnitude
                float gm = mag * gaussian_weight;
                // Tri-linear interpolation

                float vr1, vr0;
                float vrc11, vrc10, vrc01, vrc00;
                float vrco110, vrco111, vrco100, vrco101, vrco010, vrco011,
                    vrco000, vrco001;

                vr1 = gm * d_rbin;
                vr0 = gm - vr1;
                vrc11 = vr1 * d_cbin;
                vrc10 = vr1 - vrc11;
                vrc01 = vr0 * d_cbin;
                vrc00 = vr0 - vrc01;
                vrco111 = vrc11 * d_obin;
                vrco110 = vrc11 - vrco111;
                vrco101 = vrc10 * d_obin;
                vrco100 = vrc10 - vrco101;
                vrco011 = vrc01 * d_obin;
                vrco010 = vrc01 - vrco011;
                vrco001 = vrc00 * d_obin;
                vrco000 = vrc00 - vrco001;

                // int idx =  y0  * nSliceStep + x0  * nRowStep + z0;
                // All coords are offseted by 1. so x=[1, 4], y=[1, 4];
                // data for -1 coord is stored at position 0;
                // data for 8 coord is stored at position 9.
                // z doesn't need to move.
                int idx = y1 * nSliceStep + x1 * nRowStep + z0;
                histBin[idx] += vrco000;

                idx++;
                histBin[idx] += vrco001;

                idx += nRowStep - 1;
                histBin[idx] += vrco010;

                idx++;
                histBin[idx] += vrco011;

                idx += nSliceStep - nRowStep - 1;
                histBin[idx] += vrco100;

                idx++;
                histBin[idx] += vrco101;

                idx += nRowStep - 1;
                histBin[idx] += vrco110;

                idx++;
                histBin[idx] += vrco111;
            }
        }

        // Discard all the edges for row and column.
        // Only retrive edges for orientation bins.
        float *dstBins = new float[nBins];
        for (int i = 1; i <= nSubregion; i++) // slice
        {
            for (int j = 1; j <= nSubregion; j++) // row
            {
                int idx = i * nSliceStep + j * nRowStep;
                // comments: how this line works.
                // Suppose you want to write w=width, y=1, due to circular
                // buffer, we should write it to w=0, y=1; since we use a
                // circular buffer, it is written into w=width, y=1. Now, we
                // fectch the data back.
                histBin[idx] = histBin[idx + nBinsPerSubregion];

                // comments: how this line works.
                // Suppose you want to write x=-1 y=1, due to circular, it
                // should be at y=1, x=width-1; since we use circular buffer,
                // the value goes to y=0, x=width, now, we need to get it back.
                if (idx != 0)
                    histBin[idx + nBinsPerSubregion + 1] = histBin[idx - 1];

                int idx1 = ((i - 1) * nSubregion + j - 1) * nBinsPerSubregion;
                for (int k = 0; k < nBinsPerSubregion; k++) {
                    dstBins[idx1 + k] = histBin[idx + k];
                }
            }
        }

        // Normalize the histogram
        float sum_square = 0.0f;
        for (int i = 0; i < nBins; i++)
            sum_square += dstBins[i] * dstBins[i];

#if (USE_FAST_FUNC == 1)
        float thr = fast_sqrt_f(sum_square) * SIFT_DESCR_MAG_THR;
#else
        float thr = sqrtf(sum_square) * SIFT_DESCR_MAG_THR;
#endif

        float tmp = 0.0;
        sum_square = 0.0;
        // Cut off the numbers bigger than 0.2 after normalized.
        for (int i = 0; i < nBins; i++) {
            tmp = fmin(thr, dstBins[i]);
            dstBins[i] = tmp;
            sum_square += tmp * tmp;
        }

// Re-normalize
// The numbers are usually too small to store, so we use
// a constant factor to scale up the numbers.
#if (USE_FAST_FUNC == 1)
        float norm_factor = SIFT_INT_DESCR_FCTR / fast_sqrt_f(sum_square);
#else
        float norm_factor = SIFT_INT_DESCR_FCTR / sqrtf(sum_square);
#endif
        for (int i = 0; i < nBins; i++)
            dstBins[i] = dstBins[i] * norm_factor;

        memcpy(kpt->descriptors, dstBins, nBins * sizeof(float));

        if (dstBins) {
            delete[] dstBins;
            dstBins = nullptr;
        }
    }

    if (histBin) {
        delete[] histBin;
        histBin = nullptr;
    }

    return 0;
}

int sift_cpu(const Image<unsigned char> &image,
             std::list<SiftKeypoint> &kpt_list, bool bExtractDescriptors)
{
    // Index of the first octave.
    int firstOctave = (SIFT_IMG_DBL) ? -1 : 0;
    // Number of layers in one octave; same as s in the paper.
    int nLayers = SIFT_INTVLS;
    // Number of Gaussian images in one octave.
    int nGpyrLayers = nLayers + 3;
    // Number of DoG images in one octave.
    int nDogLayers = nLayers + 2;
    // Number of octaves according to the size of image.
    int nOctaves = (int)my_log2((float)fmin(image.w, image.h)) - 3 -
                   firstOctave; // 2 or 3, need further research

    // Build image octaves
    std::vector<Image<unsigned char>> octaves(nOctaves);
    build_octaves(image, octaves, firstOctave, nOctaves);

#if (DUMP_OCTAVE_IMAGE == 1)
    char foctave[256];
    for (int i = 0; i < nOctaves; i++) {
        sprintf(foctave, "octave_Octave-%d.pgm", i);
        write_pgm(foctave, octaves[i].data, octaves[i].w, octaves[i].h);
    }
#endif

    // Build Gaussian pyramid
    std::vector<Image<float>> gpyr(nOctaves * nGpyrLayers);
    build_gaussian_pyramid(octaves, gpyr, nOctaves, nGpyrLayers);

#if (DUMP_GAUSSIAN_PYRAMID_IMAGE == 1)
    char fgpyr[256];
    for (int i = 0; i < nOctaves; i++) {
        for (int j = 0; j < nGpyrLayers; j++) {
            sprintf(fgpyr, "gpyr-%d-%d.pgm", i, j);
            write_float_pgm(fgpyr, gpyr[i * nGpyrLayers + j].data,
                            gpyr[i * nGpyrLayers + j].w,
                            gpyr[i * nGpyrLayers + j].h, 1);
        }
    }
#endif

    // Build DoG pyramid
    std::vector<Image<float>> dogPyr(nOctaves * nDogLayers);
    build_dog_pyr(gpyr, dogPyr, nOctaves, nDogLayers);

#if (DUMP_DOG_IMAGE == 1)
    char fdog[256];
    Image<unsigned char> img_dog_t;
    for (int i = 0; i < nOctaves; i++) {
        for (int j = 0; j < nDogLayers; j++) {
            sprintf(fdog, "dog_Octave-%d_Layer-%d.pgm", i, j);
            img_dog_t = dogPyr[i * nDogLayers + j].to_unsigned char();
            write_pgm(fdog, img_dog_t.data, img_dog_t.w, img_dog_t.h);
        }
    }
#endif

    // Build gradient and rotation pyramids
    std::vector<Image<float>> grdPyr(nOctaves * nGpyrLayers);
    std::vector<Image<float>> rotPyr(nOctaves * nGpyrLayers);
    build_grd_rot_pyr(gpyr, grdPyr, rotPyr, nOctaves, nLayers);

    // Detect keypoints
    detect_keypoints(dogPyr, grdPyr, rotPyr, nOctaves, nDogLayers, kpt_list);

    // Extract descriptor
    if (bExtractDescriptors)
        extract_descriptor(grdPyr, rotPyr, nOctaves, nGpyrLayers, kpt_list);

    return 0;
}

} // end namespace ezsift
