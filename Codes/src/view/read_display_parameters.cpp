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
#include "read_display_parameters.h"

enum parvals { PV_NULL, PV_DEFAULT, PV_RANGE, PV_PAUSE, PV_EXIT, PV_RANGE_CONST, PV_RANGE_LIN };

void display_par2(FILE *fid, const unsigned int step, struct display_pars *displ ) 
{
  int retval;
  char c;
  char keyword[32]; /* one of default or from */
  char par[32];     /* azimuth, inclination, ... */
  float val, val2;
  char sval[128];
  int ival;
  unsigned int step0, step1;
  unsigned int numiter;
  int readparval;
  static unsigned int iter = 0;

  /* default <par> <val>
   * from <step0> to <step1> : <par> <val0> [<val1>] <interp_method:not yet implemented>
   * ...
   *
   * 1:azimuth 2:inclination 3:twist 4:distance 5:pan_x 6:pan_y 
   * [7: signal_name [8: guides [9: hide[3] ] ] ]
   */

  rewind(fid);

  /* read file line by line */
  while ( ( retval = fscanf(fid, "%s", keyword) ) != EOF )
  {
    readparval = PV_NULL;

    if ( strncmp(keyword, "exit", 31) == 0 || 
         strncmp(keyword, "quit", 31) == 0 )
    {
      fscanf(fid, " at %u", &step0);
      if ( step == step0 )
      {
        exit( EXIT_SUCCESS );
      } 
    }
    else if ( strncmp(keyword, "pause", 31) == 0 ) /* pause time
                                               * pause at <step> for <numiter> <par> <val> <val2>
                                               * linear range over _iter_ 
                                               * */
    {
      fscanf(fid, " at %u for %u", &step0, &numiter);
      if ( step == step0 && iter < numiter )
      {
        readparval = PV_PAUSE;
        displ->pause = 1;
      }
    }
    else if ( strncmp(keyword, "from", 31) == 0 ) /* then scan _step0_ to _step1_ steps */
    {
      fscanf(fid, "%u to %u", &step0, &step1);
      if ( step >= step0 && step <= step1 )
      {
        readparval = PV_RANGE;
      }
      /* printf("--step %d, step0 %d, step1 %d, read %d\n", step, step0, step1, readparval);  */
    }
    else if ( strncmp(keyword, "default", 31) == 0 ) /* case default */
    {
      readparval = PV_DEFAULT;
    }

    /* read <par> <val> */
    if ( readparval != PV_NULL )
    {
      fscanf(fid, "%s", par); /* scan parameter name */
      int index;
      for ( index = 0; index < NB_DISP_PARS; index++ ) /* find index of parameter */
      {
        if ( strncmp(par,par_map.names[index],31) == 0 ) break;
      }
      if ( strncmp(par_map.types[index],"float",31) == 0 ) /* parameter is a float */
      {
        if ( readparval == PV_DEFAULT ) /* scan a single _val_ */
        {
          fscanf(fid, "%f", &val);
          displ->fval[index] = val;
        }
        else if ( readparval == PV_RANGE ) /* scan 1 or 2 _vals_ */
        {
          retval =  fscanf(fid, "%f%*[ ]%f", &val, &val2);
          if ( retval == 1 )
          {
            displ->fval[index] = val;
          }
          else if ( retval == 2 )
          {
            displ->fval[index] = val + (step - step0)*(val2 - val)/(step1 - step0); 
          }
        }
        else if ( readparval == PV_PAUSE ) /* scan two _vals_ */
        {
          printf("--PAUSE\n");
          retval =  fscanf(fid, "%f%*[ ]%f", &val, &val2);
          displ->fval[index] = val + iter*(val2 - val)/numiter;
          iter++;
        }
        /* printf("--step %d, update %s to %f\n", step, par_map.names[index], displ->fval[index]); */
      }
      else if ( strncmp(par_map.types[index],"string",31) == 0 ) /* parameter is a string */
      {
        fscanf(fid, "%s", sval); /* scan one value */
        strncpy(displ->sval,sval,127); 
      }
      else /* parameter is an int */
      {
        fscanf(fid, "%d", &ival); /* scan one value */
        displ->ival = ival;
      }

    }

    /* discard the rest of the line */
    do {retval = fscanf(fid, "%c", &c); /* printf("--discard '%c'\n",c); */ }  while ((c != '\n') && (retval != EOF));

  }

#if 0
  printf("--%d %f %f %f %f %f %f %f %f %f %s %d %d\n", step, displ->fval[0], displ->fval[1], displ->fval[2], 
      displ->fval[3], displ->fval[4], displ->fval[5], 
      displ->fval[6], displ->fval[7], displ->fval[8], 
      displ->sval, displ->ival, displ->pause);
#endif

}

void display_par1(FILE *fid, const unsigned int step, struct display_pars *displ ) 
{
  int retval;
  char c;
  char str[128], strx[32], stry[32], strz[32];
  char guides[2];

  (void)step; /* unused parameter */

  /* display_version 1 - each line specifies all display parameters */
  /* 1:azimuth 2:inclination 3:twist 4:distance 5:pan_x 6:pan_y 
   * [7: signal_name [8: guides [9: hide[3] ] ] ]
   */
  fscanf(fid, "%f %f %f %f %f %f", 
    displ->fval, displ->fval + 1, displ->fval + 2, displ->fval + 3, displ->fval + 4, displ->fval + 5);

  /* expects either whitespace, \n or EOF. whitespace may be followed one or more identifiers */
  retval = fscanf(fid, 
      "%*[^a-zA-Z_\n]%[a-zA-Z0-9_]"         /* str: signal name */
      "%*[^0-7\n]%[0-7]"                    /* guides: 0-7 */
      "%*[^0-9.\n-]%[0-9.eE-]"              /* strx: _x_plane float */
      "%*[^0-9.\n-]%[0-9.eE-]"              /* stry: _y_plane float */
      "%*[^0-9.\n-]%[0-9.eE-]",             /* strz: _z_plane float */
      str, guides, strx, stry, strz);
  /* printf("-- retval %d, %s %s %s %s %s\n",retval, str, guides, strx, stry, strz); */
  /* printf("-- retval = %d, %f %s '%s'\n", retval, _azim_angle, str, guides); */
  if ( retval >= 1 ) /*  found identifier, try to interpret it as a signal */
  {
    strncpy(displ->sval,str,127);
    displ->sval[127] = 0; /* always terminate the string */
  }
  if ( retval >= 2 ) /*  found a number,  try to interpret it as guides */
  {
    displ->ival = strtol(guides,NULL,10);
  }
  if ( retval == 5 ) /* found hide plane position */
  {
    displ->fval[6] = strtod(strx,NULL);
    displ->fval[7] = strtod(stry,NULL);
    displ->fval[8] = strtod(strz,NULL);
  }
  do {retval = fscanf(fid, "%c", &c); /* printf("%c", c); */}  while ((c != '\n') && (retval != EOF));
  if (retval == EOF) /* loop over the file */ 
  {
    rewind(fid);
    fscanf(fid, "version %*d\n"); /* scan first line */ 
  }

}
