/* colormap.h
 * 
 */

#define GL_GLEXT_PROTOTYPES
#define GL_SILENCE_DEPRECATION

#ifdef __linux__
#include <GL/glut.h>
#endif

#if defined(__APPLE__) && defined(__MACH__) 
#include <GLUT/glut.h>
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

/* rgb color in float format */
struct rgbf {
  float r; 
  float g; 
  float b; 
};

#ifndef SRGBACOLOR_T
struct srgbacolor {
  GLfloat c[4];
};
typedef struct srgbacolor rgbacolor;
#define SRGBACOLOR_T
#endif

/* convert rgb color from hexadecimal to rgb float rgbf */
void hex2rgbf(unsigned long c, struct rgbf *rgb); 

/* make the colormap */
int make_colormap_array(const unsigned long *color_in, const int nb_color_in, rgbacolor **color_out, const int nb_color_out );

