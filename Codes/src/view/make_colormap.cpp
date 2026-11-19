// ****************************************************************************
//
//              View - Part of SiMuScale - Multi-scale simulation framework
//
// ****************************************************************************
//
// Copyright: See the AUTHORS file provided with the package
// E-mail: simuscale-contact@lists.gforge.inria.fr
// Original Authors : Samuel Bernard, Carole Knibbe, David Parsons
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU Affero General Public License as
// published by the Free Software Foundation, either version 3 of the
// License, or (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Affero General Public License for more details.
//
// You should have received a copy of the GNU Affero General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.
//
// ****************************************************************************

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "make_colormap.h"

int make_colormap_array(const unsigned long *color_in, const int nb_color_in, rgbacolor **color_out, const int nb_color_out )
{
  struct rgbf *rgbf_in, rgbf;
  int i, j;
  int segment_size = (nb_color_out - 1)/(nb_color_in - 1);
  int nb_segments = nb_color_in - 1;

  if ( ( nb_color_out - 1 ) % (nb_color_in - 1 ) )
  {
    fprintf(stderr,"warning: (nb_color_out - 1 = %d) not divisible by (nb_color_in - 1 = %d)\n",nb_color_out-1,nb_color_in-1);
  }

  *color_out = (struct srgbacolor*)malloc(nb_color_out*sizeof(struct srgbacolor));

  rgbf_in = (struct rgbf *)malloc(nb_color_in*sizeof(struct rgbf));
  for ( i = 0; i < nb_color_in; i++ )
  {
    hex2rgbf(color_in[i],rgbf_in + i);
  }

  int row = 0;
  for ( i = 0; i < nb_segments; i++ )
  {
    for ( j = 0; j < segment_size; j++)
    {
      rgbf.r = (float)j/segment_size*rgbf_in[i+1].r + (1 - (float)j/segment_size)*rgbf_in[i].r;
      rgbf.g = (float)j/segment_size*rgbf_in[i+1].g + (1 - (float)j/segment_size)*rgbf_in[i].g;
      rgbf.b = (float)j/segment_size*rgbf_in[i+1].b + (1 - (float)j/segment_size)*rgbf_in[i].b;
      /* fprintf(stderr,"%d: %g %g %g\n", row, rgbf.r, rgbf.g, rgbf.b);   */
      (*color_out)[row].c[0] = rgbf.r;
      (*color_out)[row].c[1] = rgbf.g;
      (*color_out)[row].c[2] = rgbf.b;
      (*color_out)[row].c[3] = 1.0f;
      row++;
    }

  }
  rgbf.r = rgbf_in[nb_segments].r;
  rgbf.g = rgbf_in[nb_segments].g;
  rgbf.b = rgbf_in[nb_segments].b;
  /* fprintf(fmap,"%g %g %g\n", rgbf.r, rgbf.g, rgbf.b); */
  (*color_out)[row].c[0] = rgbf.r;
  (*color_out)[row].c[1] = rgbf.g;
  (*color_out)[row].c[2] = rgbf.b;
  (*color_out)[row].c[3] = 1.0f;
  

  free(rgbf_in);

  return 0;

}

void hex2rgbf(unsigned long c, struct rgbf *rgb)
{
  unsigned char crgb[3];
  crgb[0] = c >> 16;
  crgb[1] = (c >> 8) % 0x100;
  crgb[2] = c % 0x100;
  rgb->r = (float)crgb[0]/255;
  rgb->g = (float)crgb[1]/255;
  rgb->b = (float)crgb[2]/255;
}

