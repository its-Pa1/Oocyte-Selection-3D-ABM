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

#define NB_DISP_PARS 11

/* 1:azimuth 2:inclination 3:twist 4:distance 5:pan_x 6:pan_y 
 * [7: signal_name [8: guides [9: hide[3] ] ] ] */
struct display_pars {
  float fval[9];
  char sval[128];
  unsigned int ival;
  int pause; 
};

const struct spar_map {
  char names[NB_DISP_PARS][32];
  char types[NB_DISP_PARS][32];
} par_map {
  .names = {
    "azimuth",
    "inclination",
    "twist",
    "distance",
    "pan_x",
    "pan_y",
    "hide_x",
    "hide_y",
    "hide_z",
    "signal",
    "guides"},
  .types = {
    "float",
    "float",
    "float",
    "float",
    "float",
    "float",
    "float",
    "float",
    "float",
    "string",
    "int"}
};

void display_par1(FILE *fid, const unsigned int step, struct display_pars *displ );
void display_par2(FILE *fid, const unsigned int step, struct display_pars *displ );
