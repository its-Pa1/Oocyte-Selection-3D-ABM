// ****************************************************************************
//
//              Populate - Part of SiMuScale - Multi-scale simulation framework
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


#define GL_GLEXT_PROTOTYPES
#define GL_SILENCE_DEPRECATION

#ifdef __linux__
#include <GL/glut.h>
#endif

#if defined(__APPLE__) && defined(__MACH__) 
#include <GLUT/glut.h>
#endif

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <getopt.h>

#define DEFAULT_SIGNAL -1 
#define SIGNAL_NOT_FOUND -1
#define MAX_SIG_STRLEN 16
#define MAX_NBR_SIGNAL 1024

#define ansi_esc "\033"
#define bb "1"
#define ul "4"
#define plain ""
#define fmt(style) ansi_esc "[" style "m"

#define S_SLICES 12
#define S_STACKS 12

/******* inline functions *******/

inline GLfloat sign(GLfloat x) 
{
  if ( x > 0.0f ) return 1.0f;
  if ( x < 0.0f ) return -1.0f;
  return 0.0f;
}

inline void clamp(GLfloat *v, GLfloat vmin, GLfloat vmax)
{
  if (*v > vmax)
  {
    *v = vmax;
  }
  else if ( *v < vmin )
  {
    *v = vmin;
  }
}

inline GLfloat normsq3v(GLfloat v[3])
{
  return v[0]*v[0] + v[1]*v[1] + v[2]*v[2];
}

/******* structures *******/

struct scell_sphere
{
  int id;
  float position[3];
  float color[3];
  float radius;
  int   colorindex;
};
typedef struct scell_sphere cell_sphere;

struct worldsize {
  float x;
  float y;
  float z;
};

GLfloat _colormap[7][3];
const int _colormap_size = 7;
GLfloat _min_signal = 0.0, _max_signal = 100.0;
GLfloat _min_clip_signal = -INFINITY, _max_clip_signal = INFINITY;
GLfloat _font_color[3] = {1.0, 1.0, 1.0}; /* default font color is white */
GLfloat _x_plane = 1.0 , _y_plane = 1.0 , _z_plane = 1.0;
int _pos_plane = 2; /* 0: x, 1: y, 2: z */
GLfloat _x_bulk, _y_bulk, _z_bulk;
GLfloat _radius = 1.0;
GLdouble _mouse_x, _mouse_y, _mouse_z;
int _within_world = 0;
#define MODE_ADD_BULK   0x01
#define MODE_REMOVE     0x02
#define MODE_PLANE      0x04
#define MODE_WORLD_BORDER 0x08
#define MODE_ADD_BULK_SELECT 0x10
int _mode = MODE_WORLD_BORDER;
GLfloat _mouse_shadow_color[4] = {0.0,0.0,0.0,0.3};
GLfloat _bulk_area_color[4] = {0.6,0.8,0.8,0.3};
GLfloat _mouse_overlay_color[4] = {0.8,0.8,0.8,0.8};
GLfloat _positioning_plane_color[3][4] = {{1.0, 0.9, 0.3, 0.3},
                                          {0.3, 1.0, 0.9, 0.3},
                                          {0.9, 0.3, 1.0, 0.3}}; 
int _current_color = 0;
FILE * _in_fid = NULL;
FILE * _flog  = NULL;
int _editflag = 0;
struct worldsize _worldsize = {40.0, 40.0, 40.0}; /* default world size */
cell_sphere * _mycells = NULL;
int _popsize = 0;
int _mycellscap = 0;
int _line = 0;
int _background_white = 0;
int _clipcell = 0;
char _inputstr[16];
char * _out_filename = NULL;
char * _in_filename = NULL;
char * _spheroid_pars[8] = {NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL};
int _nbr_spheroid = 0;
unsigned int _seed = 1;
GLuint _SphereDList;
enum actions { MOVE_EYE, TWIST_EYE, ZOOM, MOVE_NONE, MOVE_PLANE };
GLint _action;
GLdouble _xStart = 0.0, _yStart = 0.0;
GLfloat _near_clip, _far_clip, _distance, _twistAngle, _incAngle, _azimAngle;
GLfloat _pan_x, _pan_y, _pan_z;
GLfloat _fovy = 45.0;

void init();
void read_file();
void display();
void reshape (int w, int h);
void reset_view();
void polar_view( GLfloat distance, GLfloat azimuth, GLfloat incidence, GLfloat twist);
void keyboard(unsigned char key, int x, int y);
void special_keyboard(int key, int x, int y);
void mouse_on_plane( GLint x, GLint y);
void saveimage();
int add_cell(GLint x, GLint y);
void add_bulk();
void add_spheroid(int i);
void draw_bulk_area();
void draw_mouse_shadow();
int remove_cell(GLint x, GLint y);
void remove_cell_index(int to_remove);
int move_plane(GLint x, GLint y);
GLvoid mouse( GLint button, GLint state, GLint x, GLint y );
GLvoid motion( GLint x, GLint y );
GLvoid passive_motion( GLint x, GLint y );
void print_help(char* prog_name);
void on_exit();

void init_colormap()
{
  _colormap[0][0] = 1.000000; _colormap[0][1] = 0.200000; _colormap[0][2] = 0.000000;
  _colormap[1][0] = 1.000000; _colormap[1][1] = 0.666667; _colormap[1][2] = 0.000000;
  _colormap[2][0] = 0.019608; _colormap[2][1] = 0.800000; _colormap[2][2] = 0.200000;
  _colormap[3][0] = 0.000000; _colormap[3][1] = 0.666667; _colormap[3][2] = 1.000000;
  _colormap[4][0] = 0.666667; _colormap[4][1] = 0.200000; _colormap[4][2] = 1.000000;
  _colormap[5][0] = 0.200000; _colormap[5][1] = 0.200000; _colormap[5][2] = 0.200000;
  _colormap[6][0] = 1.000000; _colormap[6][1] = 1.000000; _colormap[6][2] = 1.000000;
}

int main(int argc, char** argv)
{

  /* Define allowed command-line options */
  const char * options_list = "he:s:o:r:";
  static struct option long_options_list[] = {
    { "help",     no_argument, NULL, 'h' },
    { "edit",     required_argument, NULL, 'e' },
    { "worlsize", required_argument, NULL, 's'},
    { "spheroid", required_argument, NULL, 'o'},
    { "seed", required_argument, NULL, 'r'},
    { 0, 0, 0, 0 }
  };
  char *endptr = NULL;


  /* Get actual values of the command-line options */
  int option;

  while ((option = getopt_long(argc, argv,
                               options_list, long_options_list, NULL)) != -1)
  {
    switch ( option )
    {
      case 'h' :
        print_help( argv[0] );
        exit( EXIT_SUCCESS );
      case 'e' :
        _editflag = 1; 
        _in_filename = optarg;
        printf("editing %s\n",_in_filename);
        break;
      case 's':
        _worldsize.x = strtof(optarg,&endptr);
        _worldsize.y = strtof(endptr,&endptr);
        _worldsize.z = strtof(endptr,&endptr);
        fprintf(_flog,"option to worldsize: %s\n",optarg);
        fprintf(_flog,"_worldsize: %f %f %f\n", _worldsize.x, _worldsize.y, _worldsize.z);
        break;
      case 'o':
        _spheroid_pars[_nbr_spheroid] = optarg;
        _nbr_spheroid++;
        break;
      case 'r':
        _seed = strtol(optarg,&endptr,10);
        break;
      default:
        print_help( argv[0] );
        exit( EXIT_FAILURE );
    }
  }

  argc -= optind;
  argv += optind;

  if ( argc )
  {
    _out_filename = argv[0];
  }



  glutInit(&argc, argv);
  glutInitDisplayMode (GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH | GLUT_STENCIL | GLUT_MULTISAMPLE);
  glutInitWindowSize (500, 500);
  glutInitWindowPosition (100, 100);
  glutCreateWindow ( "Simuscale - Populate cells" );
  init ();

  glutDisplayFunc(display);
  glutReshapeFunc(reshape);
  glutKeyboardFunc(keyboard);
  glutSpecialFunc(special_keyboard);
  glutMouseFunc(mouse);
  glutMotionFunc(motion);
  glutPassiveMotionFunc(passive_motion);

  atexit(on_exit);

  glutPostRedisplay();
  glutMainLoop();

  return 0;
}


void init(void)
{

  init_colormap();

  srand(_seed);
  
  _mycells = (cell_sphere*)malloc(sizeof(cell_sphere));
  _mycellscap = 1;

  glPixelStorei(GL_UNPACK_ALIGNMENT, 1);

   if ( _background_white )
     glClearColor (1.0, 1.0, 1.0, 1.0);
   else
     glClearColor (0.0, 0.0, 0.0, 0.0);

   glShadeModel (GL_SMOOTH);
   /* glShadeModel (GL_FLAT); */

   /* color of the light */
   GLfloat white_light[] = {1.0, 1.0, 1.0, 1.0};
   glLightfv(GL_LIGHT0, GL_AMBIENT, white_light);
   glLightfv(GL_LIGHT0, GL_DIFFUSE, white_light);
   glLightfv(GL_LIGHT0, GL_SPECULAR, white_light);

   /* position of the light */
   GLfloat light_position[] = {0.0, 0.0, 0.0, 1.0}; /* will move with the viewpoint */
   /* w=0 thus parallel rays like the sun, called directional light source */
   /* The direction (along z axis here) is transformed by the current modelview matrix */
   /* There is no attenuation for a directional light source */
   glLightfv(GL_LIGHT0, GL_POSITION, light_position);

   /* lighting model */
   GLfloat global_ambient_light[] = {0.3, 0.3, 0.3, 1.0}; /* to see objects even if no light source */
   glLightModelfv(GL_LIGHT_MODEL_AMBIENT, global_ambient_light);
   glLightModeli(GL_LIGHT_MODEL_LOCAL_VIEWER, GL_FALSE); /* infinite viewpoint */
   glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, GL_FALSE); /* back faces are inside the spheres, never seen */

   /* material for the objects */
   GLfloat mat_specular[] = {0.0, 0.0, 0.0, 1.0};
   GLfloat mat_ambient_refl[] = {0.2, 0.2, 0.2, 1.0};
   GLfloat mat_shininess[] = {0.0};
   GLfloat mat_diffuse_color[] = {0.5, 0.5, 0.5, 1.0};
   // GLfloat mat_emission[] = {0.05, 0.05, 0.05, 0.0};
   glMaterialfv(GL_FRONT, GL_AMBIENT, mat_ambient_refl);
   glMaterialfv(GL_FRONT, GL_SPECULAR, mat_specular);
   glMaterialfv(GL_FRONT, GL_SHININESS, mat_shininess);
   glMaterialfv(GL_FRONT, GL_DIFFUSE, mat_diffuse_color);
   // glMaterialfv(GL_FRONT, GL_EMISSION, mat_emission);
   glColorMaterial(GL_FRONT, GL_DIFFUSE); /* now glColor changes diffusion color */

   glBlendFunc(GL_SRC_ALPHA,GL_ONE_MINUS_SRC_ALPHA);

   glEnable(GL_LIGHT0);
   glEnable(GL_BLEND);
   glEnable(GL_RESCALE_NORMAL);
   glEnable(GL_DEPTH_TEST);
   glEnable(GL_COLOR_MATERIAL);
   glEnable(GL_MULTISAMPLE);

   /* Enable stencil operations */
   glEnable(GL_STENCIL_TEST);
   glStencilOp(GL_KEEP, GL_KEEP, GL_REPLACE);

   /* Set up near_clip and far_clip so that ( far_clip - near_clip ) > maxObjectSize, */
   /* and determine the viewing distance (adjust for zooming) */
   _near_clip = 0.1*_worldsize.y;
   _far_clip = _near_clip + 4.0*_worldsize.y;
   reset_view();

   /* display list to store the geometry of a sphere */
   _SphereDList = glGenLists(1);
   glNewList(_SphereDList, GL_COMPILE);
   glutSolidSphere(1.0, S_SLICES, S_STACKS);
   glEndList();

   if ( ( _flog = fopen("populate.log","w") ) == NULL )
   {
     fprintf(stderr, "Error, could not open populate.log.\n");
     exit(EXIT_FAILURE);
   }

   if ( _editflag ) 
   { /* try reading existing data */ 
     read_file();
   }

   if ( _nbr_spheroid ) 
   {
     for ( int i = 0; i < _nbr_spheroid; i++ )
     {
      add_spheroid(i); 
     }
   }

}

void read_file()
{
  int id;
  float x,y,z;
  float r;

  if ( ( _in_fid = fopen(_in_filename,"r") ) == NULL )
  {
    /* file does not exist, nothing to read */
    return;
  }
  while ( fscanf(_in_fid,"%d %f %f %f %f %d\n",&id,&x,&y,&z,&r,&_current_color) != EOF )
  {
    if ( _mycellscap <= _popsize )
    {
      _mycells = (cell_sphere*)realloc(_mycells,2*_mycellscap*sizeof(cell_sphere));
      _mycellscap *= 2;    
    }

    _mycells[_popsize].id = id;
    _mycells[_popsize].position[0] = x; 
    _mycells[_popsize].position[1] = y; 
    _mycells[_popsize].position[2] = z; 
    _mycells[_popsize].radius = r;
    _mycells[_popsize].color[0] = _colormap[_current_color][0]; 
    _mycells[_popsize].color[1] = _colormap[_current_color][1]; 
    _mycells[_popsize].color[2] = _colormap[_current_color][2];
    _mycells[_popsize].colorindex = _current_color;
    _popsize++;
    
  }
  fclose(_in_fid);

}

void display(void)
{

  char legend_string[128];
  GLubyte raster_rect[32] = {0x00, 0x00, 0x00, 0x00, 0x3f, 0xfc, 0x3f, 0xfc, 0x3f, 0xfc, 0x3f, 0xfc, 0x3f, 0xfc, 0x3f, 0xfc, 0x3f, 0xfc, 0x3f, 0xfc, 0x3f, 0xfc, 0x3f, 0xfc, 0x3f, 0xfc, 0x3f, 0xfc, 0x00, 0x00, 0x00, 0x00};
  GLubyte raster_border[32] = {0x1f, 0xf8, 0x7f, 0xfe, 0xc0, 0x03, 0xc0, 0x03, 0xc0, 0x03, 0xc0, 0x03, 0xc0, 0x03, 0xc0, 0x03, 0xc0, 0x03, 0xc0, 0x03, 0xc0, 0x03, 0xc0, 0x03, 0xc0, 0x03, 0xc0, 0x03, 0x7f, 0xfe, 0x1f, 0xf8};

  /***** Clear color buffer, depth buffer, stencil buffer *****/
  glClearStencil(0); /* default value */
  glClear (GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT | GL_STENCIL_BUFFER_BIT );

  /*************** Draw the text and legends ******************/
  snprintf(legend_string, 128, "r=%.2f, x=%.2f, y=%.2f, z=%.2f %s %s", 
      _radius,_mouse_x,_mouse_y,_mouse_z,_mode & MODE_PLANE ? "move plane" : "", 
      _mode & MODE_ADD_BULK ? "bulk" : ( _mode & MODE_REMOVE ? "remove" : "" ) ); 
  glColor3fv (_font_color);
  glWindowPos2i(15, 15); /* also sets current raster color to white */
  for (size_t i = 0; i < strlen(legend_string); i++) glutBitmapCharacter(GLUT_BITMAP_HELVETICA_12, legend_string[i]);

  snprintf(legend_string, 128, "| d: %s | m: %s | left: %s | %s", 
      _mode & MODE_REMOVE ? "add mode" : "remove mode",
      _mode & MODE_PLANE ? "rotate mode" : "plane mode",
      _mode & MODE_PLANE ? "move plane" : "rotate view",
      _within_world ? ( _mode & MODE_REMOVE ? "shift + left: remove cell" : "shift + left: add cell" ) : "" );
  glColor3fv (_font_color);
  glWindowPos2i(15, glutGet(GLUT_WINDOW_HEIGHT) - 15); /* also sets current raster color to white */
  for (size_t i = 0; i < strlen(legend_string); i++) glutBitmapCharacter(GLUT_BITMAP_HELVETICA_12, legend_string[i]);

  glColor3fv(_colormap[_current_color]);
  glWindowPos2i(350, 10);
  glBitmap(16, 16, 0.0, 0.0, 0, 0, raster_rect);
  glColor3f(1.0,1.0,1.0);
  glWindowPos2i(350, 10);
  glBitmap(16, 16, 0.0, 0.0, 0, 0, raster_border);
  snprintf(legend_string, 128, "color index=%d",_current_color); 
  glColor3fv (_font_color);
  glWindowPos2i(370, 15); /* also sets current raster color to white */
  for (size_t i = 0; i < strlen(legend_string); i++) glutBitmapCharacter(GLUT_BITMAP_HELVETICA_12, legend_string[i]);

  /*********************** Set the view **************************/
  glLoadIdentity ();             /* clear the matrix */
  /* viewing transformation  */
  /* gluLookAt (-1.5*_worldsize, -1.5*_worldsize, 1.5*_worldsize, _worldsize, _worldsize, 0.0, _worldsize, _worldsize, 1.5*_worldsize); */
  polar_view( _distance, _azimAngle, _incAngle, _twistAngle );

  /******************** Draw the world borders ********************/
  glColor3f (1.0, 1.0, 1.0); // White
  glLineWidth(1.0);
  glPolygonMode(GL_FRONT_AND_BACK, GL_LINE); // Mesh

  float world_size[3] = {_worldsize.x, _worldsize.y, _worldsize.z};
  float world_margins[3] = {0.05f*_worldsize.x, 0.05f*_worldsize.y, 0.05f*_worldsize.z};
  float total_world_size[3] = {world_size[0] + 2 * world_margins[0],
                               world_size[1] + 2 * world_margins[1],
                               world_size[2] + 2 * world_margins[2]};

  if ( _mode & MODE_WORLD_BORDER )
  {
    // Draw outer world borders
    glBegin(GL_QUAD_STRIP);
    glVertex3d(-world_margins[0], -world_margins[1], -world_margins[2]);
    glVertex3d(-world_margins[0], -world_margins[1] + total_world_size[1], -world_margins[2]);
    glVertex3d(-world_margins[0] + total_world_size[0], -world_margins[1], -world_margins[2]);
    glVertex3d(-world_margins[0] + total_world_size[0], -world_margins[1] + total_world_size[1], -world_margins[2]);
    glVertex3d(-world_margins[0] + total_world_size[0], -world_margins[1], -world_margins[2] + total_world_size[2]);
    glVertex3d(-world_margins[0] + total_world_size[0], -world_margins[1] + total_world_size[1], -world_margins[2] + total_world_size[2]);
    glVertex3d(-world_margins[0], -world_margins[1], -world_margins[2] + total_world_size[2]);
    glVertex3d(-world_margins[0], -world_margins[1] + total_world_size[1], -world_margins[2] + total_world_size[2]);
    glVertex3d(-world_margins[0], -world_margins[1], -world_margins[2]);
    glVertex3d(-world_margins[0], -world_margins[1] + total_world_size[1], -world_margins[2]);
    glEnd();

    // Draw inner world borders
    glBegin(GL_QUAD_STRIP);
    glVertex3d(0, 0, 0);
    glVertex3d(0, world_size[1], 0);
    glVertex3d(world_size[0], 0, 0);
    glVertex3d(world_size[0], world_size[1], 0);
    glVertex3d(world_size[0], 0, world_size[2]);
    glVertex3d(world_size[0], world_size[1], world_size[2]);
    glVertex3d(0, 0, world_size[2]);
    glVertex3d(0, world_size[1], world_size[2]);
    glVertex3d(0, 0, 0);
    glVertex3d(0, world_size[1], 0);
    glEnd();
  }
  glPolygonMode(GL_FRONT_AND_BACK, GL_FILL); // Back to fill mode

  /******************* Draw the cells ***************************/
  glEnable(GL_LIGHTING);
  for (int i = 0; i < _popsize; i++)
  {

    /**** draw cell i *****/
    glPushMatrix(); /* save V */
    glColor3fv(_mycells[i].color);
    glTranslatef(_mycells[i].position[0], _mycells[i].position[1], _mycells[i].position[2]); /* current matrix is VT  */
    glScalef(_mycells[i].radius, _mycells[i].radius, _mycells[i].radius ); /* current matrix is VTS */
    /* glCallList(_SphereDList); */
    glutSolidSphere(1.0, S_SLICES, S_STACKS);
    glPopMatrix(); /* current matrix is restored to V */

  }
  glDisable(GL_LIGHTING);

  /******************* Draw positioning plane *******************/
  glBegin(GL_QUADS);
  if ( _pos_plane == 0 ) /* x-plane x = _x_plane */ 
  {
    glColor4fv (_positioning_plane_color[0]); 
    glVertex3d(_x_plane, -world_margins[1], -world_margins[2]);
    glVertex3d(_x_plane, -world_margins[1], -world_margins[2] + total_world_size[2]);
    glVertex3d(_x_plane, -world_margins[1] + total_world_size[1], -world_margins[2] + total_world_size[2]);
    glVertex3d(_x_plane, -world_margins[1] + total_world_size[1], -world_margins[2]);
  }
  else if ( _pos_plane == 1 ) /* y-plane y = _y_plane */
  {
    glColor4fv (_positioning_plane_color[1]); 
    glVertex3d(-world_margins[0], _y_plane, -world_margins[2]);
    glVertex3d(-world_margins[0], _y_plane, -world_margins[2] + total_world_size[2]);
    glVertex3d(-world_margins[0] + total_world_size[0], _y_plane, -world_margins[2] + total_world_size[2]);
    glVertex3d(-world_margins[0] + total_world_size[0], _y_plane, -world_margins[2]);
  }
  else if ( _pos_plane == 2 ) /* z-plane z = _z_plane */
  {
    glColor4fv (_positioning_plane_color[2]); 
    glVertex3d(-world_margins[0], -world_margins[1], _z_plane);
    glVertex3d(-world_margins[0] + total_world_size[0], -world_margins[1], _z_plane);
    glVertex3d(-world_margins[0] + total_world_size[0], -world_margins[1] + total_world_size[1], _z_plane);
    glVertex3d(-world_margins[0], -world_margins[1] + total_world_size[1], _z_plane);


  }
  glEnd();


  /* guide cell at mouse position + square shadow at position plane */
  if ( _within_world ) /* if mouse points to coordinate inside world */ 
  {

    if ( _mode & MODE_REMOVE )
    {
      glutSetCursor(GLUT_CURSOR_DESTROY);
      draw_mouse_shadow();
    }
    else if ( _mode & MODE_ADD_BULK )
    {
      glutSetCursor(GLUT_CURSOR_LEFT_ARROW);
      if ( _mode & MODE_ADD_BULK_SELECT )
        draw_bulk_area();
    }
    else
    {
      glutSetCursor(GLUT_CURSOR_NONE);
      draw_mouse_shadow();
      /******* draw guide cell *******/ 
      glEnable(GL_LIGHTING);
      glDepthMask(GL_FALSE);
      glPushMatrix(); /* save V */
      // glColor4fv(_mouse_overlay_color);
      glColor4f(_colormap[_current_color][0],
                _colormap[_current_color][1],
                _colormap[_current_color][2],
                0.8);
      glTranslatef(_mouse_x, _mouse_y, _mouse_z); /* current matrix is VT  */
      glScalef(_radius,_radius,_radius); /* current matrix is VTS */
      /* glCallList(_SphereDList); */
      glutSolidSphere(1.0, S_SLICES, S_STACKS);
      glPopMatrix(); /* current matrix is restored to V */
      glDepthMask(GL_TRUE);
      glDisable(GL_LIGHTING);
    }
  }
  else
  {
    glutSetCursor(GLUT_CURSOR_LEFT_ARROW);
  }

  glutSwapBuffers();
}

void reshape (int w, int h)
{
  glViewport (0, 0, (GLsizei) w, (GLsizei) h);
  glMatrixMode (GL_PROJECTION);
  glLoadIdentity ();
  gluPerspective( _fovy, (1.0*w)/h, _near_clip, _far_clip );
  glMatrixMode (GL_MODELVIEW);
}


void reset_view()
{
  _distance = _near_clip + (_far_clip - _near_clip) / 2.0;
  _pan_x = 0.0; 
  _pan_y = 0.0;
  _pan_z = 0.0;
  _twistAngle = 0.0;	/* rotation of viewing volume (camera) */
  _incAngle = 85.0;
  _azimAngle = 30.0;
}


void polar_view( GLfloat distance, GLfloat azimuth, GLfloat incidence, GLfloat twist)
{
  glTranslatef( -_worldsize.x/2.0, -_worldsize.y/2.0, -distance);
  glTranslatef( _pan_x, _pan_y, _pan_z);
  glRotatef( -twist, 0.0f, 0.0f, 1.0);
  glRotatef( -incidence, 1.0f, 0.0f, 0.0);
  glRotatef( -azimuth, 0.0f, 0.0f, 1.0);
}


void keyboard(unsigned char key, int x, int y)
{

  switch (key)
  {
    case 'r': /* r or R : reset viewpoint */
      reset_view();
      break;
    case 'R':
      reset_view();
      break;
    case 'w': /* toggle show/hide world border display */
      _mode ^= MODE_WORLD_BORDER;
      break;
    case 'x': /* toggle x-plane positioning */
      _pos_plane = 0;
      break;
    case 'y': /* toggle y-plane positioning */
      _pos_plane = 1;
      break;
    case 'z': /* toggle z-plane positioning */
      _pos_plane = 2;
      break;
    case 'm':
      _mode ^= MODE_PLANE;
      break;
    case 'd':
      _mode ^= MODE_REMOVE;
      break;
    case '+':
      _radius += 0.1;
      break;
    case '-':
      _radius -= 0.1;
      break;
    case ']':
      _current_color++;
      _current_color %= _colormap_size;
      break;
    case '[':
      _current_color += _colormap_size - 1;
      _current_color %= _colormap_size;
      break;
    case 'u':
      remove_cell_index(_popsize - 1);
      break;
    case 'a':
      add_cell(x,y);
      break;
    case 'X':
      remove_cell(x,y);
      break;
    case 'b':
      _mode &= ~MODE_ADD_BULK_SELECT;
      _mode ^= MODE_ADD_BULK;
      break;
    case 's':
      saveimage();
      break;
    case 27:  /* escape to quit */
    case 'q': /* q to quit */
    case 'Q': /* Q to quit */
      exit(0);
  }
  glutPostRedisplay();
}

void special_keyboard(int key, int, int)
{
  int modifiers = glutGetModifiers();
  switch (key)
  {
    case GLUT_KEY_UP:
      /* printf("Up key pressed\n"); */
      if ( modifiers & GLUT_ACTIVE_SHIFT )
      {
        _pan_z += 0.025*_worldsize.x;
      }
      else
      {
        _pan_y -= 0.025*_worldsize.z;
      }
      glutPostRedisplay();
      break;
    case GLUT_KEY_DOWN:
      /* printf("Down key pressed\n"); */
      if ( modifiers & GLUT_ACTIVE_SHIFT )
      {
        _pan_z -= 0.025*_worldsize.x;
      }
      else
      {
        _pan_y += 0.025*_worldsize.z;
      }
      glutPostRedisplay();
      break;
    case GLUT_KEY_LEFT:
      /* printf("Left key pressed\n"); */
      _pan_x += 0.025*_worldsize.x;
      glutPostRedisplay();
      break;
    case GLUT_KEY_RIGHT:
      /* printf("Right key pressed\n"); */
      _pan_x -= 0.025*_worldsize.x;
      glutPostRedisplay();
      break;
  }
}

GLvoid mouse( GLint button, GLint state, GLint x, GLint y )
{
  int modifiers = glutGetModifiers();
  _action = MOVE_NONE;
  if (state == GLUT_DOWN)
  {
    switch (button)
    {
      case GLUT_LEFT_BUTTON:
        if ( modifiers & GLUT_ACTIVE_SHIFT )
        {
          if ( _mode & MODE_ADD_BULK  )
          {
            if ( _mode & MODE_ADD_BULK_SELECT ) 
            {
              add_bulk();
            }
            else
            {
              _x_bulk = _mouse_x;
              _y_bulk = _mouse_y;
              _z_bulk = _mouse_z;
            }
            _mode ^= MODE_ADD_BULK_SELECT;
          }
          else if ( _mode & MODE_REMOVE ) 
          {
            remove_cell(x,y);
          }
          else
          {
            add_cell(x,y); 
          }
        }
        else if ( _mode & MODE_PLANE )
        {
          _action = MOVE_PLANE;
        }
        else 
        {
          _action = MOVE_EYE;
        }
        break;
      case GLUT_MIDDLE_BUTTON:
        _action = TWIST_EYE;
        break;
      case GLUT_RIGHT_BUTTON:
        _action = ZOOM;
        break;
    }

    /* Update the saved mouse position */
    _xStart = x;
    _yStart = y;

  }

}

int add_cell(GLint x, GLint y)
{
  int window_height = glutGet(GLUT_WINDOW_HEIGHT);
  GLfloat depth;
  GLint stencil;
  GLdouble model[16];
  GLdouble projection[16];
  GLint view[4];
  GLdouble x_coord,y_coord,z_coord;

  glReadPixels(x, window_height - y - 1, 1, 1, GL_STENCIL_INDEX, GL_UNSIGNED_INT, &stencil);
  glReadPixels(x, window_height - y - 1, 1, 1, GL_DEPTH_COMPONENT, GL_FLOAT, &depth);

  glGetDoublev(GL_MODELVIEW_MATRIX,model);
  glGetDoublev(GL_PROJECTION_MATRIX,projection);
  glGetIntegerv(GL_VIEWPORT,view);

  /* get approximate nearest object coordinates */
  gluUnProject((GLdouble)x, (GLdouble)(window_height - y - 1),
      depth,model,projection,view,&x_coord,&y_coord,&z_coord);

  /* check that coords are within world boundaries */
  if ( x_coord < 0.0 || x_coord > _worldsize.x ||
       y_coord < 0.0 || y_coord > _worldsize.x || 
       z_coord < 0.0 || z_coord > _worldsize.x )
  {
    fprintf(stderr, "\aCoordinates (%f,%f,%f) are out of bound\n", x_coord,y_coord,z_coord);
    return 1;
  }
       
  

  if ( _mycellscap <= _popsize )
  {
    _mycells = (cell_sphere*)realloc(_mycells,2*_mycellscap*sizeof(cell_sphere));
    _mycellscap *= 2;    
  }

  _mycells[_popsize].id = _popsize+1;
  _mycells[_popsize].position[0] = x_coord; 
  _mycells[_popsize].position[1] = y_coord; 
  _mycells[_popsize].position[2] = z_coord; 
  _mycells[_popsize].radius = _radius;
  _mycells[_popsize].color[0] = _colormap[_current_color][0]; 
  _mycells[_popsize].color[1] = _colormap[_current_color][1]; 
  _mycells[_popsize].color[2] = _colormap[_current_color][2];
  _mycells[_popsize].colorindex = _current_color;

  _popsize++;


  glutPostRedisplay();
  
  return 0;
}

void add_bulk()
{
  GLfloat x0, x1, y0, y1, z0, z1;

  switch(_pos_plane)
  {
    case 0: /* x-plane */
      x0 = _x_plane;
      x1 = _x_plane + 0.1;
      y0 = _y_bulk;
      y1 = _mouse_y;
      z0 = _z_bulk;
      z1 = _mouse_z;
      break;
    case 1: /* y-plane */
      x0 = _x_bulk;
      x1 = _mouse_x;
      y0 = _y_plane;
      y1 = _y_plane - 0.1;
      z0 = _z_bulk;
      z1 = _mouse_z;
      break;
    case 2: /* z-plane */
      x0 = _x_bulk;
      x1 = _mouse_x;
      y0 = _y_bulk;
      y1 = _mouse_y;
      z0 = _z_plane;
      z1 = _z_plane + 0.1;
  }


  for ( int i = 0; i < 10; i++ )
  {

    if ( _mycellscap <= _popsize )
    {
      _mycells = (cell_sphere*)realloc(_mycells,2*_mycellscap*sizeof(cell_sphere));
      _mycellscap *= 2;    
    }

    _mycells[_popsize].id = _popsize+1;
    _mycells[_popsize].position[0] = x0 + (GLfloat)rand()/RAND_MAX*(x1 - x0); 
    _mycells[_popsize].position[1] = y0 + (GLfloat)rand()/RAND_MAX*(y1 - y0);
    _mycells[_popsize].position[2] = z0 + (GLfloat)rand()/RAND_MAX*(z1 - z0);
    _mycells[_popsize].radius = _radius;
    _mycells[_popsize].color[0] = _colormap[_current_color][0]; 
    _mycells[_popsize].color[1] = _colormap[_current_color][1]; 
    _mycells[_popsize].color[2] = _colormap[_current_color][2];
    _mycells[_popsize].colorindex = _current_color;

    _popsize++;
 
  }
  
  glutPostRedisplay();

}

void add_spheroid(int ind)
{
  char    *endptr = NULL;
  GLfloat  v[3], p[3];
  int      i = 0;
  
  /* read spheroid options 
   * x y z r n 
   */
  GLfloat c[3] = {  strtof(_spheroid_pars[ind],&endptr),
                    strtof(endptr,&endptr),
                    strtof(endptr,&endptr) };
  GLfloat spheriod_radius = strtof(endptr,&endptr);
  GLfloat cell_radius     = strtof(endptr,&endptr);
  GLint color = strtol(endptr,&endptr,10);
  GLint n     = strtol(endptr,&endptr,10);
  GLint dim   = strtol(endptr,&endptr,16);

  const int dimx = 0x04;
  const int dimy = 0x02;
  const int dimz = 0x01;

  while ( i < n ) 
  {
    if ( _mycellscap <= _popsize )
    {
      _mycells = (cell_sphere*)realloc(_mycells,2*_mycellscap*sizeof(cell_sphere));
      _mycellscap *= 2;    
    }

    if ( dim & dimx )
      v[0] = -1 + 2.0*(GLfloat)rand()/RAND_MAX; 
    else
      v[0] = 0.0;
    if ( dim & dimy )
      v[1] = -1 + 2.0*(GLfloat)rand()/RAND_MAX;
    else
      v[1] = 0.0;
    if ( dim & dimz )
      v[2] = -1 + 2.0*(GLfloat)rand()/RAND_MAX;
    else
      v[2] = 0.0;

    p[0] = c[0] + spheriod_radius*v[0];
    p[1] = c[1] + spheriod_radius*v[1];
    p[2] = c[2] + spheriod_radius*v[2];

    /* fail conditions */
    if ( p[0] < 0.0 || p[0] > _worldsize.x ||
         p[1] < 0.0 || p[1] > _worldsize.y ||
         p[2] < 0.0 || p[2] > _worldsize.z )
    {  continue; }

    if ( normsq3v(v) > 1 )
    { continue; }


    /* pass: add new cell */
    _mycells[_popsize].id = _popsize+1;
    _mycells[_popsize].position[0] = p[0]; 
    _mycells[_popsize].position[1] = p[1];
    _mycells[_popsize].position[2] = p[2];
    _mycells[_popsize].radius = cell_radius;
    _mycells[_popsize].color[0] = _colormap[color][0]; 
    _mycells[_popsize].color[1] = _colormap[color][1]; 
    _mycells[_popsize].color[2] = _colormap[color][2];
    _mycells[_popsize].colorindex = color;

    _popsize++;
    i++;
  }
}

void draw_bulk_area()
{

  glDepthMask(GL_FALSE);
  glBegin(GL_QUADS);
  glColor4fv (_bulk_area_color);
  if ( _pos_plane == 0 ) /* x-plane x = _x_plane */ 
  {
    glVertex3d(_x_plane + 0.1, _y_bulk, _z_bulk);
    glVertex3d(_x_plane + 0.1, _mouse_y, _z_bulk);
    glVertex3d(_x_plane + 0.1, _mouse_y, _mouse_z);
    glVertex3d(_x_plane + 0.1, _y_bulk, _mouse_z);
  }
  else if ( _pos_plane == 1 ) /* y-plane y = _y_plane */
  {
    glVertex3d(_x_bulk, _y_plane - 0.1, _z_bulk);
    glVertex3d(_mouse_x, _y_plane - 0.1, _z_bulk);
    glVertex3d(_mouse_x, _y_plane - 0.1, _mouse_z);
    glVertex3d(_x_bulk, _y_plane - 0.1, _mouse_z);
  }
  else if ( _pos_plane == 2 ) /* z-plane z = _z_plane */
  {
    glVertex3d(_x_bulk, _y_bulk, _z_plane + 0.1);
    glVertex3d(_mouse_x, _y_bulk, _z_plane + 0.1);
    glVertex3d(_mouse_x, _mouse_y, _z_plane + 0.1);
    glVertex3d(_x_bulk, _mouse_y, _z_plane + 0.1);
  }
  glEnd();
  glDepthMask(GL_TRUE);
}

void draw_mouse_shadow()
{
  glDepthMask(GL_FALSE);
  glBegin(GL_QUADS);
  if ( _pos_plane == 0 ) /* x-plane x = _x_plane */ 
  {
    glColor4fv (_mouse_shadow_color);
    glVertex3d(_x_plane + 0.1, _mouse_y - 0.8, _mouse_z - 0.8);
    glVertex3d(_x_plane + 0.1, _mouse_y - 0.8, _mouse_z + 0.8);
    glVertex3d(_x_plane + 0.1, _mouse_y + 0.8, _mouse_z + 0.8);
    glVertex3d(_x_plane + 0.1, _mouse_y + 0.8, _mouse_z - 0.8);
  }
  else if ( _pos_plane == 1 ) /* y-plane y = _y_plane */
  {
    glColor4fv (_mouse_shadow_color);
    glVertex3d(_mouse_x - 0.8, _y_plane - 0.1, _mouse_z - 0.8);
    glVertex3d(_mouse_x - 0.8, _y_plane - 0.1, _mouse_z + 0.8);
    glVertex3d(_mouse_x + 0.8, _y_plane - 0.1, _mouse_z + 0.8);
    glVertex3d(_mouse_x + 0.8, _y_plane - 0.1, _mouse_z - 0.8);
  }
  else if ( _pos_plane == 2 ) /* z-plane z = _z_plane */
  {
    glColor4fv (_mouse_shadow_color);
    glVertex3d(_mouse_x - 0.8, _mouse_y - 0.8, _z_plane + 0.1);
    glVertex3d(_mouse_x - 0.8, _mouse_y + 0.8, _z_plane + 0.1);
    glVertex3d(_mouse_x + 0.8, _mouse_y + 0.8, _z_plane + 0.1);
    glVertex3d(_mouse_x + 0.8, _mouse_y - 0.8, _z_plane + 0.1);
  }
  glEnd();
  glDepthMask(GL_TRUE);
}

int remove_cell(GLint x, GLint y)
{
  int window_height = glutGet(GLUT_WINDOW_HEIGHT);
  GLfloat depth;
  GLdouble model[16];
  GLdouble projection[16];
  GLint view[4];
  GLdouble x_coord,y_coord,z_coord;
  GLdouble min_dist = INFINITY;
  GLdouble dist;

  int to_remove;

  glReadPixels(x, window_height - y - 1, 1, 1, GL_DEPTH_COMPONENT, GL_FLOAT, &depth);

  glGetDoublev(GL_MODELVIEW_MATRIX,model);
  glGetDoublev(GL_PROJECTION_MATRIX,projection);
  glGetIntegerv(GL_VIEWPORT,view);

  /* get approximate nearest object coordinates */
  gluUnProject((GLdouble)x, (GLdouble)(window_height - y - 1),
      depth,model,projection,view,&x_coord,&y_coord,&z_coord);

  /* find cell with center closest to mouse coordinates */
  /* this will return a cell ID even if user clicked far away */
  for ( int i = 0; i < _popsize ; i++ )
  {
    dist = (_mycells[i].position[0] - x_coord)*(_mycells[i].position[0] - x_coord) +
         (_mycells[i].position[1] - y_coord)*(_mycells[i].position[1] - y_coord) +
         (_mycells[i].position[2] - z_coord)*(_mycells[i].position[2] - z_coord);
    if ( dist < min_dist )
    {
      min_dist = dist; 
      to_remove = i;
    }
  }

  remove_cell_index(to_remove);

  glutPostRedisplay();
  return 1;

}

void remove_cell_index(int to_remove)
{

  if ( to_remove < 0 )
  {
    return;
  }

  for ( int i = to_remove; i < _popsize -  1 ; i++ )
  {
    /* shift cell i+1 -> cell i */
    _mycells[i].position[0] = _mycells[i+1].position[0];
    _mycells[i].position[1] = _mycells[i+1].position[1]; 
    _mycells[i].position[2] = _mycells[i+1].position[2];
    _mycells[i].radius = _mycells[i+1].radius;
    _mycells[i].color[0] = _mycells[i+1].color[0]; 
    _mycells[i].color[1] = _mycells[i+1].color[1];
    _mycells[i].color[2] = _mycells[i+1].color[2];
    _mycells[i].colorindex = _mycells[i+1].colorindex;
  }
  _popsize--;

}

int move_plane(GLint x, GLint y)
{

  int window_height = glutGet(GLUT_WINDOW_HEIGHT);
  GLfloat depth;
  GLdouble model[16];
  GLdouble projection[16];
  GLint view[4];
  GLdouble new_x, new_y, new_z;
  static GLdouble old_x, old_y, old_z;

  glReadPixels(x, window_height - y - 1, 1, 1, GL_DEPTH_COMPONENT, GL_FLOAT, &depth);

  glGetDoublev(GL_MODELVIEW_MATRIX,model);
  glGetDoublev(GL_PROJECTION_MATRIX,projection);
  glGetIntegerv(GL_VIEWPORT,view);

  /* get approximate nearest object coordinates */
  gluUnProject((GLdouble)x, (GLdouble)(window_height - y - 1),
      0.0,model,projection,view,&new_x,&new_y,&new_z);

  switch(_pos_plane)
  {
    case 0:
      _x_plane += 0.02*_worldsize.x*sign(new_x - old_x)*(fabs(new_x - old_x) < 0.02*_worldsize.x);
      clamp(&_x_plane,0.0,_worldsize.x);
      break;
    case 1:
      _y_plane += 0.02*_worldsize.y*sign(new_y - old_y)*(fabs(new_y - old_y) < 0.02*_worldsize.y);
      clamp(&_y_plane,0.0,_worldsize.y);
      break;
    case 2:
      _z_plane += 0.02*_worldsize.z*sign(new_z - old_z)*(fabs(new_z - old_z) < 0.02*_worldsize.z);
      clamp(&_z_plane,0.0,_worldsize.z);
      break;
  }


  old_x = new_x;
  old_y = new_y;
  old_z = new_z;

  glutPostRedisplay();
  
  return 0;
}

GLvoid motion( GLint x, GLint y )
{
  switch (_action)
  {
    case MOVE_EYE:
      /* Adjust the eye position based on the mouse position */
      _azimAngle += (GLdouble) (x - _xStart);
      _incAngle -= (GLdouble) (y - _yStart);
      break;
    case TWIST_EYE:
      /* Adjust the eye twist based on the mouse position */
      _twistAngle = fmod(_twistAngle+(x - _xStart), 360.0);
      break;
    case ZOOM:
      /* Adjust the eye distance based on the mouse position */
      _distance -= (GLdouble) (y - _yStart)/10.0;
      break;
    case MOVE_PLANE:
      move_plane(x,y);
      break;
  }

  /* update position of mouse on plane */
  mouse_on_plane(x,y);

  /* Update the stored mouse position for later use */
  _xStart = x;
  _yStart = y;

  glutPostRedisplay();
}

GLvoid passive_motion( GLint x, GLint y)
{

  mouse_on_plane(x,y);

  glutPostRedisplay();
}

void mouse_on_plane( GLint x, GLint y)
{
  int window_height = glutGet(GLUT_WINDOW_HEIGHT);
  int window_width = glutGet(GLUT_WINDOW_WIDTH);
  GLfloat depth;
  GLdouble model[16];
  GLdouble projection[16];
  GLint view[4];

  if ( x < 1 | x > window_width | y < 1 | y > window_height )
  { return; }

  glReadPixels(x, window_height - y - 1, 1, 1, GL_DEPTH_COMPONENT, GL_FLOAT, &depth);

  glGetDoublev(GL_MODELVIEW_MATRIX,model);
  glGetDoublev(GL_PROJECTION_MATRIX,projection);
  glGetIntegerv(GL_VIEWPORT,view);

  /* get approximate nearest object coordinates */
  gluUnProject((GLdouble)x, (GLdouble)(window_height - y - 1),
      depth,model,projection,view,&_mouse_x,&_mouse_y,&_mouse_z);

  if ( _mouse_x >= 0 && _mouse_x <= _worldsize.x && 
       _mouse_y >= 0 && _mouse_y <= _worldsize.y && 
       _mouse_z >= 0 && _mouse_z <= _worldsize.z )
  {
    _within_world = 1;
  }
  else
  {
    _within_world = 0;
  }

}


void saveimage() 
{
    int window_height = glutGet(GLUT_WINDOW_HEIGHT);
    int window_width  = glutGet(GLUT_WINDOW_WIDTH);
    char *buffer;
    char filepath[64];
    FILE *fid;
    char mn[] = "BM";
    char creator[] = "SCPP";
    unsigned int offset = 0x36;  /* = 54: header=14 and DIB=40 */
    unsigned int dibsize = 0x28; /* = 40 */
    unsigned int zero = 0x0;     /* = 00 00 00 00 */
    unsigned int resolution = 2835; /* points per meter = 72dpi */
    GLsizei nr_channels = 3;         
    unsigned int image_size;
    unsigned int buffer_size = window_width * window_height * nr_channels; /* size in bytes */
    image_size = window_width * window_height * nr_channels + offset; /* total image size on disk */
    buffer = (char*)malloc(buffer_size * sizeof(char));
    glPixelStorei(GL_PACK_ALIGNMENT, 4);
    glReadBuffer(GL_FRONT);
    glReadPixels(0, 0, window_width, window_height, GL_BGR, GL_UNSIGNED_BYTE, buffer);

    /* filepath: sSTEPNUMBER_nPOPSIZE_SIGNAL.bmp */
    snprintf(filepath,64,"init.bmp");
    if ( ( fid = fopen(filepath,"w") ) == NULL )
    {
          fprintf(stderr, "ERROR : could not open file '%s'.\n", filepath );
    }
    printf("image: %s\n",filepath);
    printf("  window_size: %u x %u\n",window_width,window_height);
    printf("  buffer size: %u bytes\n",buffer_size);
    printf("  total size:  %u bytes\n",image_size);

    /* bitmap header: 14 bytes */
    fwrite(mn,1,2,fid);             /* bitmap image type: BM */
    fwrite(&image_size,4,1,fid);    /* total image size */
    fwrite(creator,1,4,fid);        /* user-defined bytes */
    fwrite(&offset,4,1,fid);        /* total offset, including DIB */

    /* DIB - BITMAPINFOHEADER */
    fwrite(&dibsize,4,1,fid);          /* DIB size: 40 bytes - 4 bytes */
    fwrite(&window_width,4,1,fid);     /* image width, in pixels  - 4 bytes */
    fwrite(&window_height,4,1,fid);    /* image height, in pixels - 4 bytes */
    fputc(1,fid);                      /* number of color planes: 1 */ 
    fputc(0,fid);                      /* pad with zeros to fill 2 bytes - 2 bytes */
    fputc(0x18,fid);                   /* number of bits/pixel n = 24 */
    fputc(0,fid);                      /* pad with zeros to fill 2 bytes - 2 bytes */
    fwrite(&zero,4,1,fid);             /* compression mode: 0 = no compression - 4 bytes */
    fwrite(&zero,4,1,fid);             /* image size: 0 when no compression - 4 bytes */
    fwrite(&resolution,4,1,fid);       /* res horiz - 4 bytes */
    fwrite(&resolution,4,1,fid);       /* res vert  - 4 bytes */
    fwrite(&zero,4,1,fid);             /* default number colors in palette: 0 to defaults to n^24 - 4 bytes */
    fwrite(&zero,4,1,fid);             /* default important colors: 0 - 4 bytes */

    /* image raw data, not memory aligned */
    fwrite(buffer,1,buffer_size,fid);   /* write image raw data */

    /* clean up */
    free(buffer);
    fclose(fid);
}

void on_exit() 
{
  FILE *out_fid;
  if ( ( out_fid = fopen(_out_filename, "w") ) == NULL)
  {
    /* _out_filename cannot be opened, perhaps because it was not set 
     * revert silently to stdout
     */
    out_fid = stdout;
  }
  for ( int i = 0; i < _popsize ; i++ )
  {
    fprintf(out_fid, "%d %f %f %f %f %d\n", _mycells[i].id,
                                           _mycells[i].position[0],
                                           _mycells[i].position[1],
                                           _mycells[i].position[2],
                                           _mycells[i].radius,
                                           _mycells[i].colorindex);
  }
  if ( out_fid != stdout )
    fclose(out_fid);
  fclose(_flog);
  free(_mycells);
}

void print_help( char* prog_name )
{
  printf( "\n\
Usage : " fmt(bb) "%s" fmt(plain) " -h\n\
   or : " fmt(bb) "%s" fmt(plain) " [options] " fmt(ul) "file" fmt(plain) "\n\n\
Create or edit an initial cell population dataset and print it to \
" fmt(ul) "file" fmt(plain) ". \n\
If " fmt(ul) "file" fmt(plain) " argument is missing, the content is printed to standard output.\n\
\n\
\tOptions\n\n\
\t" fmt(bb) "-h" fmt(plain) " or " fmt(bb) "--help" fmt(plain) "       Display this screen\n\
\t" fmt(bb) "-e" fmt(plain) " or " fmt(bb) "--edit" fmt(plain) " " fmt(ul) "file" fmt(plain) "  If file exists, first load existing\n\
\t                   data from file.\n\
\t" fmt(bb) "-s" fmt(plain) " or " fmt(bb) "--worldsize" fmt(plain) " \
" fmt(ul) "'x y z'" fmt(plain) " \n\
\t                   Specify the world size 'x y z' (default: 40 40 40)\n\
\t" fmt(bb) "-o" fmt(plain) " or " fmt(bb) "--spheroid" fmt(plain) " \
" fmt(ul) "'x y z s r c n d'" fmt(plain) " \n\
\t                   Generate a spheroid with n cells centered at x y z, \n\
\t                   of total radius s, cell radius r, color c (c=0 to 6),\n\
\t                   and along dimension d (octal value: 4:x, 2:y, 1:z). \n\
\t                   Up to 8 spheroids can be generated.\n\
\t                   \n\
\n\
\tRuntime commands\n\
\n\
\t" fmt(bb) "left arrow" fmt(plain) " \tpan left\n\
\t" fmt(bb) "right arrow" fmt(plain) "\tpan right\n\
\t" fmt(bb) "up arrow" fmt(plain) " \tpan up\n\
\t" fmt(bb) "down arrow" fmt(plain) "\tpan down\n\
\t" fmt(bb) "shift + up arrow" fmt(plain) "   \tzoom in (same as right mouse button)\n\
\t" fmt(bb) "shift + down arrow" fmt(plain) " \tzoom out (same as right mouse button)\n\
\t" fmt(bb) "-" fmt(plain) "," fmt(bb) "+" fmt(plain) "        \tDecrease, increase radius\n\
\t" fmt(bb) "a" fmt(plain) "          \tAdd a cell a mouse position            \n\
\t" fmt(bb) "d" fmt(plain) "          \tToggle cell delete mode\n\
\t" fmt(bb) "m" fmt(plain) "          \tToggle plane moving mode\n\
\t" fmt(bb) "r" fmt(plain) " or " fmt(bb) "R" fmt(plain) "     \tReset 3D view\n\
\t" fmt(bb) "s" fmt(plain) "          \tSave current frame as an image (bitmap)\n\
\t" fmt(bb) "w" fmt(plain) "          \tToggle worldborder\n\
\t" fmt(bb) "x" fmt(plain) ", " fmt(bb) "y" fmt(plain) ", or " fmt(bb) "z" fmt(plain) " \tActivate x, y, or z positioning plane\n\
\n\
\t" fmt(bb) "right mouse" fmt(plain) " \tzoom in or out\n\
\t" fmt(bb) "left mouse" fmt(plain) "  \tIn rotate mode: rotate view,\n\
\t                    \tIn plane mode: move plane along its normal\n\
\t" fmt(bb) "shift + left click" fmt(plain) "  \tIn add mode: add a new cell\n\
\t                In remove mode: delete cell\n\
\n\
\tPress ESC or q to quit.\n\n",\
   prog_name, prog_name );
}

