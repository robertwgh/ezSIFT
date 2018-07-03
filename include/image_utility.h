/*	Copyright (c) 2013, Robert Wang, email: robertwgh (at) gmail.com
	All rights reserved. https://sourceforge.net/p/ezsift

	Revision history:
		September, 15, 2013: initial version.
*/

#ifndef _IMG_IO_H_
#define _IMG_IO_H_

typedef struct {
  int w;
  int h;
  unsigned char * img_r;
  unsigned char * img_g;
  unsigned char * img_b;
} PPM_IMG;

int read_pgm(const char* filename, unsigned char * & data, int & w, int & h);
void write_pgm(const char* filename, unsigned char * data, int w, int h);

int read_ppm(const char * filename, unsigned char * & data, int & w, int & h);
void write_ppm(const char * filename, unsigned char * data, int w, int h);
void write_rgb2ppm(const char * filename, unsigned char* r, unsigned char* g, 
					unsigned char* b, int w, int h) ;

void write_float_pgm(const char* filename, float* data, int w, int h,  int mode);

void rasterCircle(PPM_IMG* imgPPM, int r, int c, int radius);
void draw_red_circle(PPM_IMG* imgPPM, int r, int c, int cR);
void draw_circle(PPM_IMG* imgPPM, int r, int c, int cR, float thickness);
void draw_red_orientation(PPM_IMG* imgPPM, int x, int y, float ori, int cR);

/****************************************
 * Utility functions
 ***************************************/
// Image operations
// Get pixel from an image with unsigned char datatype.
inline unsigned char get_pixel(
	unsigned char * imageData, 
	int w, int h, 
	int r, int c)
{	
	unsigned char val; 
	if ( c >= 0 && c < w && r >= 0 && r < h)
	{
		val = imageData[r * w + c];
	}else if (c < 0){
		val = imageData[r * w];
	}else if (c >= w){
		val = imageData[r * w + w - 1];
	}else if (r < 0){
		val = imageData[c];
	}else if (r >= h){
		val = imageData[(h-1) * w + c];
	}else{
		val = 0;
	}
	return val;
}

// Get pixel value from an image with float data type.
inline float get_pixel_f(
	float * imageData, 
	int w, int h, 
	int r, int c)
{	
	float val; 
	if ( c >= 0 && c < w && r >= 0 && r < h)
	{
		val = imageData[r * w + c];
	}else if (c < 0){
		val = imageData[r * w];
	}else if (c >= w){
		val = imageData[r * w + w - 1];
	}else if (r < 0){
		val = imageData[c];
	}else if (r >= h){
		val = imageData[(h-1) * w + c];
	}else{
		val = 0.0f;
	}
	return val;
}

#endif//_IMG_IO_H_
