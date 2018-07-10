/* Copyright (c) 2013, Robert Wang, email: robertwgh (at) gmail.com
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
   an object in an image David G. Lowe, US Patent 6,711,293 (March 23, 2004).
   Provisional application filed March 8, 1999. Asignee: The University of
   British Columbia.

   Revision history:
      September, 15, 2013: initial version.
      July 8th, 2014, re-organized source code.
*/

#ifndef EZSIFT_H
#define EZSIFT_H

#include "image.h"
#include <list>
#include <vector>

namespace ezsift {

/****************************************
 * Constant parameters
 ***************************************/

// default number of sampled intervals per octave
static int SIFT_INTVLS = 3;

// default sigma for initial gaussian smoothing
static float SIFT_SIGMA = 1.6f;

// the radius of Gaussian filter kernel;
// Gaussian filter mask will be (2*radius+1)x(2*radius+1).
// People use 2 or 3 most.
static float SIFT_GAUSSIAN_FILTER_RADIUS = 3.0f;

// default threshold on keypoint contrast |D(x)|
static float SIFT_CONTR_THR = 8.0f; // 8.0f;

// default threshold on keypoint ratio of principle curvatures
static float SIFT_CURV_THR = 10.0f;

// The keypoint refinement smaller than this threshold will be discarded.
static float SIFT_KEYPOINT_SUBPiXEL_THR = 0.6f;

// double image size before pyramid construction?
static bool SIFT_IMG_DBL = false; // true;

// assumed gaussian blur for input image
static float SIFT_INIT_SIGMA = 0.5f;

// width of border in which to ignore keypoints
static int SIFT_IMG_BORDER = 5;

// maximum steps of keypoint interpolation before failure
static int SIFT_MAX_INTERP_STEPS = 5;

// default number of bins in histogram for orientation assignment
static int SIFT_ORI_HIST_BINS = 36;

// determines gaussian sigma for orientation assignment
static float SIFT_ORI_SIG_FCTR =
    1.5f; // Can affect the orientation computation.

// determines the radius of the region used in orientation assignment
static float SIFT_ORI_RADIUS =
    3 * SIFT_ORI_SIG_FCTR; // Can affect the orientation computation.

// orientation magnitude relative to max that results in new feature
static float SIFT_ORI_PEAK_RATIO = 0.8f;

// maximum number of orientations for each keypoint location
// static const float SIFT_ORI_MAX_ORI = 4;

// determines the size of a single descriptor orientation histogram
static float SIFT_DESCR_SCL_FCTR = 3.f;

// threshold on magnitude of elements of descriptor vector
static float SIFT_DESCR_MAG_THR = 0.2f;

// factor used to convert floating-point descriptor to unsigned char
static float SIFT_INT_DESCR_FCTR = 512.f;

// default width of descriptor histogram array
static int SIFT_DESCR_WIDTH = 4;

// default number of bins per histogram in descriptor array
static int SIFT_DESCR_HIST_BINS = 8;

// default value of the nearest-neighbour distance ratio threshold
// |DR_nearest|/|DR_2nd_nearest|<SIFT_MATCH_NNDR_THR is considered as a match.
static float SIFT_MATCH_NNDR_THR = 0.65f;

#if 0
// intermediate type used for DoG pyramids
typedef short sift_wt;
static const int SIFT_FIXPT_SCALE = 48;
#else
// intermediate type used for DoG pyramids
typedef float sift_wt;
static const int SIFT_FIXPT_SCALE = 1;
#endif

/****************************************
 * Definitions
 ***************************************/
#define DEGREE_OF_DESCRIPTORS (128)
struct SiftKeypoint {
    int octave;   // octave number
    int layer;    // layer number
    float rlayer; // real number of layer number

    float r;     // normalized row coordinate
    float c;     // normalized col coordinate
    float scale; // normalized scale

    float ri;          // row coordinate in that layer.
    float ci;          // column coordinate in that layer.
    float layer_scale; // the scale of that layer

    float ori; // orientation in degrees.
    float mag; // magnitude

    float descriptors[DEGREE_OF_DESCRIPTORS];
};

// Match pair structure. Use for interest point matching.
struct MatchPair {
    int r1;
    int c1;
    int r2;
    int c2;
};

/****************************************
 *  SIFT Processing Functions
 ***************************************/
// Initialize SIFT parameters.
void init_sift_parameters(bool doubleFirstOctave = true,
                          float contrast_threshold = 8.0f,
                          float edge_threshold = 10.0f,
                          float match_NDDR_threshold = 0.6f);

// Enable doubling of original image.
void double_original_image(bool doubleFirstOctave);

// Efficient Gaussian Blur function.
// 1. Use row buf to handle border pixel.
// 2. hori processing and transpose
int gaussian_blur(const Image<float> &in_image, Image<float> &out_image,
                  std::vector<float> coef1d);

// Row filter and then transpose
int row_filter_transpose(float *src, float *dst, int w, int h, float *coef1d,
                         int gR);

// Build image octaves during the initialization.
int build_octaves(const Image<unsigned char> &image,
                  std::vector<Image<unsigned char>> &octaves, int firstOctave,
                  int nOctaves);

// Compute Gaussian filter coefficients for Gaussian Blur.
std::vector<std::vector<float>> compute_gaussian_coefs(int nOctaves,
                                                       int nGpyrLayers);

// Build Gaussian pyramid.
int build_gaussian_pyramid(std::vector<Image<unsigned char>> &octaves,
                           std::vector<Image<float>> &gpyr, int nOctaves,
                           int nGpyrLayers);

// Build DoG pyramid.
int build_dog_pyr(std::vector<Image<float>> &gpyr,
                  std::vector<Image<float>> &dogPyr, int nOctaves,
                  int nDogLayers);

// Build gradient and rotation pyramid.
int build_grd_rot_pyr(std::vector<Image<float>> &gpyr,
                      std::vector<Image<float>> &grdPyr,
                      std::vector<Image<float>> &rotPyr, int nOctaves,
                      int nLayers);

// Refine local extrema.
bool refine_local_extrema(std::vector<Image<float>> &dogPyr, int nOctaves,
                          int nDogLayers, SiftKeypoint &kpt);

// Export keypoint list to a file.
int export_kpt_list_to_file(const char *filename,
                            std::list<SiftKeypoint> &kpt_list,
                            bool bIncludeDescpritor);

// Compute orientation histogram.
float compute_orientation_hist(const Image<float> &image, SiftKeypoint &kpt,
                               float *&hist);

/****************************************
 *  SIFT Core Functions
 ***************************************/
// Detect keypoints.
int detect_keypoints(std::vector<Image<float>> &dogPyr,
                     std::vector<Image<float>> &grdPyr,
                     std::vector<Image<float>> &rotPyr, int nOctaves,
                     int nDogLayers, std::list<SiftKeypoint> &kpt_list);

// Extract descriptor.
int extract_descriptor(std::vector<Image<float>> &grdPyr,
                       std::vector<Image<float>> &rotPyr, int nOctaves,
                       int nGpyrLayers, std::list<SiftKeypoint> &kpt_list);

/****************************************
 *  SIFT Interface Functions
 ***************************************/
// Detect keypoints and extract descriptor.
int sift_cpu(const Image<unsigned char> &image,
             std::list<SiftKeypoint> &kpt_list, bool bExtractDescriptors);

// Match keypoints from two keypoint lists.
int match_keypoints(std::list<SiftKeypoint> &kpt_list1,
                    std::list<SiftKeypoint> &kpt_list2,
                    std::list<MatchPair> &match_list);

} // end namespace ezsift

#endif
