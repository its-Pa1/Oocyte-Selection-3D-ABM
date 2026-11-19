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
// Additional contributions: Pawan Kumar
// GRN, Cell_Germ, Cell_Somatic, updated view files, Python scripts,
// and universal CMake / cross-platform support.
// These implementations were used to produce results reported in the
// article cited in the main README.
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

#define GL_GLEXT_PROTOTYPES
#define GL_SILENCE_DEPRECATION

#ifdef __linux__
#include <GL/glut.h>
#endif

#if defined(__APPLE__) && defined(__MACH__) 
#include <GLUT/glut.h>
#endif

#define ansi_esc "\033"
#define bb "1"
#define ul "4"
#define plain ""
#define fmt(style) ansi_esc "[" style "m"
#define num_opt fmt(plain)fmt(ul) "num" fmt(plain)
#define name_opt fmt(plain)fmt(ul) "name" fmt(plain)
#define file_opt fmt(plain)fmt(ul) "file" fmt(plain)
#define options_opt fmt(plain)fmt(ul) "options" fmt(plain)

#define DEFAULT_SIGNAL -2 
#define SIGNAL_NOT_FOUND -1
#define MAX_SIG_STRLEN 16
#define MAX_NBR_SIGNAL 1024
#define INITIAL_TS_RANGE 10.0f

#define PR_BS_CLEAR printf("\033[D\033[J")

#if 0
#define S_SLICES 10
#define S_STACKS 10
#endif

#define NB_COLORMAPS 13
#define NB_COLORS    65

/******* structures *******/

struct scell_sphere
{
  int id;
  float position[3];
  float orientation[3];
  float force[3];
  float velocity[3];
  float color[4];
  char  cell_tag[32];
  char  living_status;
  float signal;
  float radius;
  int   lineage;
  int   highlight;
  int   clip;
};
typedef struct scell_sphere cell_sphere;

struct worldsize {
  float x;
  float y;
  float z;
};

struct striangle
{
  GLfloat v1[3];
  GLfloat v2[3];
  GLfloat v3[3];
  GLfloat normal[3];
};
typedef struct striangle triangle;

struct stimeseries {
  GLfloat *t;
  GLfloat *y;
  GLfloat *y2;
  GLuint  *step;
  GLfloat tmin;
  GLfloat tmax;
  GLfloat ymin;
  GLfloat ymax;
  size_t  size;
  size_t  maxsize;
  size_t  id;
};
typedef struct stimeseries timeseries;

#ifndef SRGBACOLOR_T
struct srgbacolor {
  GLfloat c[4];
};
typedef struct srgbacolor rgbacolor;
#define SRGBACOLOR_T
#endif

/******* constants *******/
static const char menu_txt[5] = "menu";

/******* enums *******/

enum actions { MOVE_EYE, TWIST_EYE, ZOOM, MOVE_NONE, MOVE_PLANE, MOVE_PLANE_FINE };
enum color_models { COLOR_MODEL_NORMAL, COLOR_MODEL_INVERSE, COLOR_MODEL_FLAT, COLOR_MODEL_GLOW };
enum stencils { STENCIL_BACKGROUND, STENCIL_BORDER, STENCIL_SCALEBAR, STENCIL_LEGEND, STENCIL_CELL, STENCIL_TRACK, STENCIL_HOOVER, STENCIL_TAG };
enum update_signal_from { SIGNAL_INDEX, SIGNAL_NAME };
enum keyboard_input_mode { KEY_OFF, KEY_TRACKING, KEY_COLORMAP };
enum menu_items { MENU_COLORMAP, SUBMENU_COLORMAP_JET, SUBMENU_COLORMAP_SPRING, SUBMENU_COLORMAP_NEON,
  SUBMENU_COLORMAP_CUSTOM, MENU_HIDE_PLANE, MENU_CLIP_CELLS, MENU_NEXT_SIGNAL, MENU_PREVIOUS_SIGNAL,
  MENU_DISPLAY_SIGNAL, MENU_RESET, MENU_HELP, MENU_QUIT, SUBMENU_DISPLAY_SIGNAL_0, SUBMENU_COLORMAP_LIST_0 };
enum cell_geometry { GEOM_NORMAL, GEOM_DEAD, GEOM_SECTION, 
  GEOM_SECTION_X_HALF, GEOM_SECTION_X_14, GEOM_SECTION_X_34,
  GEOM_SECTION_Y_HALF, GEOM_SECTION_Y_14, GEOM_SECTION_Y_34, 
  GEOM_SECTION_Z_HALF, GEOM_SECTION_Z_14, GEOM_SECTION_Z_34,
  GEOM_ARROW, GEOM_SCALEBAR };
enum cell_hoovering { HOOV_DISABLED, HOOV_OFF, HOOV_ON };
enum shade_model { SHADE_SMOOTH, SHADE_FLAT };
enum extra_option { OPT_NONE, OPT_TAG, OPT_HIDE, OPT_GUIDE, OPT_SNAPSHOT, OPT_STEP, 
                    OPT_BACKGROUND, OPT_WINDOW_SIZE, OPT_ANIMATE, OPT_LINEAGE,
                    OPT_SCALEBARPOS };

/*** cell drawing masks ***/ 
#define DRAW_CELLS       0x4
#define DRAW_VELOCITIES  0x2
#define DRAW_FORCES      0x1

/*** guide drawing masks ***/
#define GUIDE_BORDER     0x1
#define GUIDE_VIEWCUT    0x2
#define GUIDE_SCALEBAR   0x4

/******* functions *******/

void init_color_map_spring();
void init_color_map_jet();
void init_color_map_neon();
void init_color_map_custom(char *filename);
void read_min_max_signals();
void read_header();
void read_file_line();
void read_display_parameters();
void animate(int value);
void init();
void display();
void draw_cell(GLint i);
int  locate_cell(GLint i);
void draw_force(GLint i);
void draw_velocity(GLint i);
void hoover_cell(GLint id, int color);
size_t get_cell_index(GLint id);
void auxdisplay();
void statsdisplay();
void create_menu();
void menu(int value);
void colormap_menu(int value);
void cell_signal_menu(int value);
void draw_ts(timeseries *ts, size_t current);
void reshape (int w, int h);
void reset_view();
void polar_view( GLfloat distance, GLfloat azimuth, GLfloat incidence, GLfloat twist);
void keyboard(unsigned char key, int x, int y);
void special_keyboard(int key, int x, int y);
void get_keyboard_input(unsigned char key);
void update_signal(int from_what);
void next_signal();
void previous_signal();
void set_signal(int sig);
void standstill();
void backtrack( unsigned int nbrframes);
void darken( float *color, float amount );
void lighten( float *color, float amount );
GLfloat distsq3v( GLfloat *p1, GLfloat *p2);
GLfloat normsq3v( GLfloat *v);
void assign3v(const GLfloat *source, size_t size_source, GLfloat *target, size_t size_target);
void copy_triangle(const triangle source, triangle *target);
void rotate3v(const GLfloat *axis, const GLfloat *angle, GLfloat *points, const size_t nbr_points);
void translate3v(const GLfloat *shift, GLfloat *points, size_t nbr_points);
void normal3v( GLfloat *p1, GLfloat *p2, GLfloat *p3, GLfloat *nv );
void object_material();
void dead_cell_geometry();
void arrow_geometry();
void scalebar_geometry();
void sphere_section(GLfloat s);
void triangularize(triangle **trgl, size_t nbr_triangles);
void draw_plane();
void move_plane(GLint x, GLint y);
void save_image(const char* filename, const char* file_type = "bmp");
void save_image();  // Overloaded version with defaults
void init_opencv_video(const char* filename, int width, int height);
void add_frame_to_video();
void cleanup_opencv_video();
void print_snapshot();
int get_lineage(int daughter);
void read_lineage_file();
int select_cell(GLint x, GLint y);
void hoover_main(GLint x, GLint y);
GLvoid mouse( GLint button, GLint state, GLint x, GLint y );
GLvoid motion( GLint x, GLint y );
GLvoid passive_motion( GLint x, GLint y );
void print_help(char* prog_name);
void on_exit();

/* constants */
const GLfloat _white_light[] = {1.0, 1.0, 1.0, 1.0};

/* base colors for viridis colormaps */
const struct { 
  char names[NB_COLORMAPS][64];
  unsigned long colors[NB_COLORMAPS][5];
}  _list_colormaps = { 
  .names = { 
    "jet", 
    "spring",
    "neon",
    "gray",
    "viridis",
    "inferno",
    "magma",
    "plasma",
    "cividis",
    "rocket",
    "mako",
    "turbo",
    "old_spring" },
  .colors = {
    {0x00007F, 0x0080FF, 0x7CFF79, 0xFF9400, 0x7F0000}, /* jet        */
    {0xFF00FF, 0xFF40BF, 0xFF807F, 0xFFC03F, 0xFFFF00}, /* spring     */
    {0x1E3553, 0x64BA6A, 0x77C5CD, 0xC876CA, 0xFAF351}, /* neon       */
    {0x020202, 0x222222, 0x555555, 0x858585, 0xFDFDFD}, /* gray       */
    {0x440154, 0x3B528B, 0x21908C, 0x5DC863, 0xFDE725}, /* viridis    */
    // {0x000004, 0x56106E, 0xBB3754, 0xF98C0A, 0xFCFFA4}, /* inferno    */
    {0x000004, 0x56106E, 0xBB3754, 0xE8743B, 0xFF9500}, /* cool_inferno - stops before yellow */
    {0x000004, 0x51127C, 0xB63679, 0xFB8861, 0xFCFDBF}, /* magma      */
    {0x0D0887, 0x7E03A8, 0xCC4678, 0xF89441, 0xF0F921}, /* plasma     */
    {0x00204D, 0x414D6B, 0x7C7B78, 0xBCAF6F, 0xFFEA46}, /* cividis    */
    {0x03051A, 0x611F53, 0xCB1B4F, 0xF4875E, 0xFAEBDD}, /* rocket     */
    {0x0B0405, 0x3E356B, 0x357BA2, 0x49C1AD, 0xDEF5E5}, /* mako       */
    {0x30123B, 0x28BBEC, 0xA2FC3C, 0xFB8022, 0x7A0403}, /* turbo      */
    // {0x116688, 0x84A05F, 0xFFDD33, 0xEE806E, 0xDD22AA}} /* old_spring */ 
    {0x0072B2, 0x56B4E9, 0xF0E442, 0xE69F00, 0xD55E00}} /* better old spring*/
    // {0xD55E00, 0xCC78BC, 0xCA9161, 0xFBAFE4, 0x949494}} /* blue-orange, no yellow */
    // {0x2D5F3F, 0x5A9F6A, 0x8FD699, 0xC4E8B4, 0xE6F5D0}} /* old_spring - fresh spring greens */ 

};
