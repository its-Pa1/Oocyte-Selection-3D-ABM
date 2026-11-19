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

#include <stdlib.h>
#include <stdio.h>
#include <signal.h>
#include <string.h>
#include <math.h>
#include <getopt.h>
#include <png.h> // libpng should be installed [for macOS: brew install libpng, for Linux: sudo apt-get install libpng-dev]
#include <time.h> // For time-related functions
#include "view.h"
#include "make_colormap.h"
#include "read_display_parameters.h"
#include <opencv2/opencv.hpp>

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

/******* Helper function to render signal names with subscripts *******/
/**
 * @brief Render a signal name with subscript formatting
 * 
 * This function renders signal names like "S_nbd" in a LaTeX-style format
 * where the text after underscore is displayed as a subscript.
 * 
 * @param signal_name The signal name string (e.g., "S_nbd")
 * @param font The font to use for regular text (e.g., GLUT_BITMAP_TIMES_ROMAN_24)
 * @param subscript_font The font to use for subscript (e.g., GLUT_BITMAP_HELVETICA_18)
 * @return Total width in pixels consumed by rendering
 */
inline int renderSignalNameWithSubscript(const char* signal_name, void* font, void* subscript_font) {
  const char* underscore = strchr(signal_name, '_');
  int total_width = 0;
  
  if (underscore != NULL) {
    // Render text before underscore normally
    for (const char* c = signal_name; c < underscore; c++) {
      glutBitmapCharacter(font, *c);
      total_width += glutBitmapWidth(font, *c);
    }
    
    // Get current raster position for subscript offset
    GLfloat raster_pos[4];
    glGetFloatv(GL_CURRENT_RASTER_POSITION, raster_pos);
    
    // Move down for subscript (negative Y)
    glWindowPos2f(raster_pos[0], raster_pos[1] - 5);
    
    // Render subscript text
    for (const char* c = underscore + 1; *c != '\0'; c++) {
      glutBitmapCharacter(subscript_font, *c);
      total_width += glutBitmapWidth(subscript_font, *c);
    }
  } else {
    // No underscore, render normally
    for (const char* c = signal_name; *c != '\0'; c++) {
      glutBitmapCharacter(font, *c);
      total_width += glutBitmapWidth(font, *c);
    }
  }
  
  return total_width;
}

/******* global variables *******/

rgbacolor * _colormap;
GLfloat _colormap_len;
GLfloat _alpha = 1.0;
GLint   _color_model = COLOR_MODEL_NORMAL;
int _using_niche = 0;
int _long_option = OPT_NONE;
char *_viewpoint_str = NULL;
char _hide_str[64] = "";
char _length_unit[64] = "";
GLfloat _min_signal = 0.0, _max_signal = 100.0;
GLfloat _min_clip_signal = -INFINITY, _max_clip_signal = INFINITY;
GLfloat _font_color[3]         =  {1.0, 1.0, 1.0}; /* default font color is white */
GLfloat _border_color[3]       =  {1.0, 1.0, 1.0}; /* default font color is white */
// GLfloat _highlight_color[4][3] = {{1.0,0.9,0.0},   /* cell hoover */
//                                   {0.8,0.4,0.0},   /* cell ancestors */
//                                   {0.3,0.9,0.6},   /* cell progeny */
//                                   {0.2,0.2,0.2}};  /* background */
GLfloat _highlight_color[4][3] = {{0.0,0.8,1.0},   /* cell hoover - bright cyan */
                                  {0.0,0.5,0.7},   /* cell ancestors - darker cyan */
                                  {0.4,0.2,0.8},   /* cell progeny - purple */
                                  {0.2,0.2,0.2}};  /* background */
GLfloat _diffuse_light[4]      =  {0.3, 0.3, 0.3 ,1.0}; 
GLfloat _light1_pos[4];
int _s_res = 12; /* number of slices and stacks for spheres */
int _blend = 1;  /* transparency enabled by default  */
FILE * _traj_file = NULL;
FILE * _log_file  = NULL;
FILE * _colormap_file = NULL;
FILE * _display_parameter_file = NULL;
FILE * _lineage_file = NULL;
char * _norm_filename = NULL;
char * _traj_filename = NULL;
struct worldsize _worldsize = {0.0, 0.0, 0.0};
cell_sphere * _mycells = NULL;
int _popsize = 0;
float _time = 0.0;
int _line = 0;
int _pause = 0;
int _end_of_file = 0;
int _skip_frame = 1;
int _no_annotations = 0;  /* 0: show annotations (default), 1: hide all annotations */
int _colorbar_only = 0;   /* 0: show all annotations (default), 1: show only colorbar */
int _minimal_annotations = 0; /* 0: normal, 1: show only time and number of cells */
GLfloat _background_color[3] = {0.0f, 0.0f, 0.0f};
unsigned int _colormap_name = 0; /* default colormap: jet */
int _clipcell = 0;
int _draw_guides = GUIDE_BORDER | GUIDE_VIEWCUT; 
int _single_cell_tracking = 0; /* 0: none, 1: tracking */
int _single_cell_hoovering = HOOV_OFF; /* HOOV_DISABLED, HOOV_OFF, HOOV_ON */
int _nbr_cell_tracking = 0;
int _keyboard_input_mode = 0;
int _single_cell = 0;
int _single_cell_hoovered = 0;
timeseries _ts = {(GLfloat *)malloc(sizeof(GLfloat)),
                  (GLfloat *)malloc(sizeof(GLfloat)),
                  NULL,
                  (GLuint *)malloc(sizeof(GLuint)),
                  0.0,INITIAL_TS_RANGE,0.0,0.1,0,1,0};
timeseries _stats_ts = {(GLfloat *)malloc(sizeof(GLfloat)),
                  (GLfloat *)malloc(sizeof(GLfloat)),
                  (GLfloat *)malloc(sizeof(GLfloat)),
                  (GLuint *)malloc(sizeof(GLuint)),
                  0.0,INITIAL_TS_RANGE,0.0,0.1,0,1,0};
char _inputstr[128];
GLuint _cell_dlist;
float _speed = 0.02;  /* default is 0.02 frames per ms, i.e. wait 50 ms between two frames */
int _first_signal_col;
int _signal_index; 
int _orientation = 0;
int _force = 0;
int _velocity = 0;
int _draw_elements = DRAW_CELLS;
int _living_status = 0;
char* _signal_name;
char* _tags = NULL;
int _tagged_cell_highlighting = 0;
GLfloat _tag_color[4] = {0.0, 1.0, 0.0, 1.0};
char* _snapshot_file = NULL;
int _initial_step = -1;
char* _display_parameter_filename = NULL;
char* _lineage_filename = NULL;
char* _lineage_str = NULL;
int _main_window_size[2] = {500, 500};
int _main_menu;
int _menu_hoovering = 0;
int _shade_model = GL_SMOOTH; 
struct { 
  unsigned long *pos; 
  unsigned long i; 
  unsigned long size; 
  unsigned long maxsize; 
}  _fpos = {(unsigned long*)malloc(128*sizeof(unsigned long)), 0, 0, 128};
struct list_signals { unsigned int size; char name[MAX_NBR_SIGNAL][128]; };
struct list_signals _list_signals = { .size = 0, .name = {} };
char _signal_found = 0;
GLint _action;
GLint _main_window;
GLint _aux_window;
GLint _stats_window;
GLdouble _x_start = 0.0, _y_start = 0.0;
GLfloat _near_clip, _far_clip, _distance, _twist_angle, _inc_angle, _azim_angle;
GLfloat _pan_x, _pan_y;
GLfloat _fovy = 45.0;
GLfloat _x_plane = 0.0 , _y_plane = 0.0 , _z_plane = 0.0;
GLfloat _scale_pos[3] = {1.0f,1.0f,0.0f}; 
int _pos_plane = 2; /* 0: x, 1: y, 2: z */
int _move_plane = 0;
GLfloat _positioning_plane_color[3][4] = {{1.0, 0.9, 0.3, 0.2},
                                          {0.3, 1.0, 0.9, 0.2},
                                          {0.9, 0.3, 1.0, 0.2}};

GLubyte _raster_rect[16] = {0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff};
GLubyte _raster_upper_left[16] = {0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0x7f, 0x7f, 0x3f, 0x0f};
GLubyte _raster_upper_right[16] = {0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xfe, 0xfe, 0xfc, 0xf0};
GLubyte _raster_lower_left[16] = {0x0f, 0x3f, 0x7f, 0x7f, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff};
GLubyte _raster_lower_right[16] = {0xf0, 0xfc, 0xfe, 0xfe, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff};
GLubyte _raster_mu[9] = {0x40, 0x40, 0x74, 0x4c, 0x44, 0x44, 0x44, 0x44, 0x44};

char _filename[256] = "default_filename.png";
int _record_video = 0;
char* _video_filename = NULL;
int _video_frame_count = 0;
int _record_window = 0; // 0=main, 1=aux, 2=stats (default: main)
cv::VideoWriter _video_writer;
int _video_codec = cv::VideoWriter::fourcc('m','p','4','v'); // int _video_codec = -1; // it should choose automatically only on Windows
int _video_fps = 25;
int _show_neighborhood_sphere = 0; /* 0: off, 1: on */
/* MAIN */
int main(int argc, char** argv)
{

  /* Define allowed command-line options */
  const char * options_list = "hn:s:i:pwc:a:m:g:v:V:r";
  static struct option long_options_list[] = {
    { "help",     no_argument, NULL, 'h' },
    { "normalization",     required_argument, NULL, 'n' },
    { "speed",    required_argument, NULL, 's' },
    { "signal",   required_argument, NULL, 'i' },
    { "pause",    no_argument, NULL, 'p' },
    { "white",    no_argument, NULL, 'w' },
    { "colormap", required_argument, NULL, 'c' },
    { "alpha",    required_argument, NULL, 'a' },
    { "no-blend",    no_argument, &_blend, 0 },
    { "colormodel",    required_argument, NULL, 'm' },
    { "color-model",    required_argument, NULL, 'm' },
    { "tags",     required_argument, NULL, 'g'},
    { "tag-color", required_argument, &_long_option, OPT_TAG}, 
    { "viewangle", required_argument, NULL, 'v' },
    { "viewpoint", required_argument, NULL, 'V' },
    { "niche",  no_argument,  &_using_niche, 1 },
    { "hide",  required_argument,  &_long_option, OPT_HIDE },
    { "no-hoovering",  no_argument,  &_single_cell_hoovering, HOOV_DISABLED },
    { "hoovering",  no_argument,  &_single_cell_hoovering, HOOV_OFF },
    { "guides",  required_argument,  &_long_option, OPT_GUIDE },
    { "snapshot",  required_argument,  &_long_option, OPT_SNAPSHOT },
    { "step",  required_argument,  &_long_option, OPT_STEP },
    { "background-color", required_argument, &_long_option, OPT_BACKGROUND},
    { "window-size", required_argument, &_long_option, OPT_WINDOW_SIZE},
    { "animate", required_argument, &_long_option, OPT_ANIMATE},
    { "lineage", required_argument, &_long_option, OPT_LINEAGE},
    { "scale-bar-position", required_argument, &_long_option, OPT_SCALEBARPOS},
    { "shade-flat", no_argument, &_shade_model, GL_FLAT},
    { "shade-smooth", no_argument, &_shade_model, GL_SMOOTH},
    { "resolution-high", no_argument, &_s_res, 20},
    { "resolution-medium", no_argument, &_s_res, 12},
    { "resolution-low", no_argument, &_s_res, 5},
    { "no-annotations", no_argument, &_no_annotations, 1},
    { "colorbar-only", no_argument, &_colorbar_only, 1},
    { "minimal-annotations", no_argument, &_minimal_annotations, 1},
    {"record", optional_argument, NULL, 'r'},
    { 0, 0, 0, 0 }
  };

  /* Get actual values of the command-line options */
  int option;
  char signal_specified = 0;


  /* open log file before loading option */
  if ( ( _log_file = fopen("view.log","w") ) == NULL )
  {
    fprintf(stderr, "Error, could not open view.log.\n");
    exit(EXIT_FAILURE);
  }

  /* parse command line options */
  while ((option = getopt_long(argc, argv,
                               options_list, long_options_list, NULL)) != -1)
  {
    switch ( option )
    {
      case 'h' :
      {
        print_help( argv[0] );
        exit( EXIT_SUCCESS );
      }
      case 'n' :
      {
        if ( strcmp( optarg, "" ) == 0 )
        {
          fprintf(stderr, "ERROR : Option -n or --normalization : missing argument.\n" );
          exit( EXIT_FAILURE );
        }

        _norm_filename = (char *) malloc((strlen(optarg) + 1)*sizeof(char));
        snprintf( _norm_filename, strlen(optarg) + 1, "%s", optarg );
        printf("Normalization file: %s\n", optarg);
        break;
      }
      case 's' :
      {
        if ( strcmp( optarg, "" ) == 0 )
          {
            fprintf(stderr,
                "ERROR : Option -s or --speed : missing argument.\n" );
            exit( EXIT_FAILURE );
          }
        _speed = atof(optarg); /* in frames per millisecond */
        break;
      }
      case 'i' :
        if ( strcmp( optarg, "" ) == 0 )
        {
            fprintf(stderr,
                "ERROR : Option -i or --signal : missing argument.\n" );
            exit( EXIT_FAILURE );
        }
        _signal_name = (char *)malloc((strlen(optarg) + 1)*sizeof(char));
        snprintf( _signal_name, strlen(optarg) + 1, "%s", optarg );
        _signal_index = SIGNAL_NOT_FOUND; /* the index will need to be found in the header file */
        signal_specified = 1;
        break;
      case 'p':
        _pause = 1;
        break;
      case 'w':
        _background_color[0] = 1.0;
        _background_color[1] = 1.0;
        _background_color[2] = 1.0;
        break;
      case 'c':
        if ( strcmp( optarg, "" ) == 0 )
        {
            fprintf(stderr,
                "ERROR : Option -c or --colormap : missing argument.\n" );
            exit( EXIT_FAILURE );
        }
        _colormap_name = NB_COLORMAPS;
        for ( int i = 0; i < NB_COLORMAPS; i++ )
        {
          if ( strcmp( optarg, _list_colormaps.names[i] ) == 0 )
          {
            _colormap_name = i; 
            break;
          }
        }
        if ( _colormap_name == NB_COLORMAPS )
        {
          init_color_map_custom(optarg);
        }
        break;
      case 'a' :
        _alpha = atof(optarg);
        break;
      case 'm' :
        _color_model = atoi(optarg);
        break;
      case 'g': 
        if ( strcmp( optarg, "" ) == 0 )
        {
            fprintf(stderr,
                "ERROR : Option -g or --tags : missing argument.\n" );
            exit( EXIT_FAILURE );
        }
        _tags = (char *)malloc((strlen(optarg) + 1)*sizeof(char));
        snprintf( _tags, strlen(optarg) + 1, "%s", optarg );
        _tagged_cell_highlighting = 1;
        break;
      case 'v':
        _fovy = atof(optarg); /* view angle, in degrees */
        if ( _fovy < 0 )
          _fovy *= -1.0;
        if ( _fovy > 180 )
          _fovy = 180.0;
        break;
      case 'V':
        _viewpoint_str = optarg;        
        break;
      case 'r':
      {
          _record_video = 1;
          if (optarg != NULL && strcmp(optarg, "") != 0) {
              // Parse format: filename:window (e.g., "video:0" for main window)
              char *colon = strchr(optarg, ':');
              if (colon != NULL) {
                  *colon = '\0';
                  _video_filename = strdup(optarg);
                  char *window_str = colon + 1;
                  int window_num = atoi(window_str);
                  if (window_num >= 0 && window_num <= 2) {
                      _record_window = window_num;
                  } else {
                      _record_window = 0;
                  }
              } else {
                  _video_filename = strdup(optarg);
                  _record_window = 0;
              }
          } else {
                  time_t rawtime;
                  struct tm * timeinfo;
                  char timestamp[32];
                  char timestamped_filename[256];
                  
                  time(&rawtime);
                  timeinfo = localtime(&rawtime);
                  strftime(timestamp, sizeof(timestamp), "%Y%m%d_%H%M%S", timeinfo);
                  snprintf(timestamped_filename, 256, "simuscale_video_%s", timestamp);
              _video_filename = strdup(timestamped_filename);
              _record_window = 0;
          }
          
          // Add .mp4 extension if not present
          char full_filename[512];
          const char* ext = strrchr(_video_filename, '.');
          if (!ext || strcmp(ext, ".mp4") != 0) {
              snprintf(full_filename, 512, "%s.mp4", _video_filename);
          } else {
              strncpy(full_filename, _video_filename, 512);
          }
          
          // Initialize OpenCV video writer (will start when first frame is added)
          printf("Video recording enabled: %s\n", full_filename);
          
          const char* window_names[] = {"main", "aux", "stats"};
          printf("Recording from %s window\n", window_names[_record_window]);
          break;
      }
      case 0: /* long option only, no short option */
        char *endptr = optarg;
        switch(_long_option)
        {
          case OPT_HIDE:
            snprintf( _hide_str, strlen(optarg) + 1, "%s", optarg );
            break;
          case OPT_TAG:
            _tag_color[0]  = strtof(endptr,&endptr);
            _tag_color[1]  = strtof(endptr,&endptr);
            _tag_color[2]  = strtof(endptr,&endptr);
            _tag_color[3]  = strtof(endptr,&endptr);
            break;  
          case OPT_GUIDE:
            _draw_guides = strtol(optarg,NULL,10);
            break;
          case OPT_SNAPSHOT:
            _snapshot_file = optarg;
            break;
          case OPT_STEP:
            _skip_frame = strtol(optarg,NULL,10);
            _skip_frame++;
            _initial_step = _skip_frame;
            break;
          case OPT_BACKGROUND:
            _background_color[0] = strtof(endptr,&endptr);
            _background_color[1] = strtof(endptr,&endptr);
            _background_color[2] = strtof(endptr,&endptr);
            break;
          case OPT_WINDOW_SIZE:
            _main_window_size[0] = strtol(endptr,&endptr,10);
            _main_window_size[1] = strtol(endptr,&endptr,10);
            break;
          case OPT_ANIMATE:
            _display_parameter_filename = optarg;
            break;
          case OPT_LINEAGE:
            _lineage_filename = optarg;
            break;
          case OPT_SCALEBARPOS:
            _scale_pos[0] = strtof(endptr,&endptr);
            _scale_pos[1] = strtof(endptr,&endptr);
            _scale_pos[2] = strtof(endptr,&endptr);
            break;
        }
        _long_option = OPT_NONE;
    }
  }

  argc -= optind;
  argv += optind;

  if ( argc == 0 )
  {
    _traj_filename = strdup("trajectory.txt");
  }
  else
  {
    _traj_filename = strdup(argv[0]); 
  }
  
  if ( ( _traj_file = fopen(_traj_filename, "r") ) == NULL)
  {
    fprintf(stderr, "Error, file %s could not be opened.\n", _traj_filename);
    exit(EXIT_FAILURE);
  }
  fprintf(_log_file,"trajectory file: %s\n", _traj_filename);
  if ( _norm_filename == NULL ) 
  {
    _norm_filename = strdup("normalization.txt");
  }
  fprintf(_log_file,"normalization file: %s\n", _norm_filename);

  if ( _display_parameter_filename != NULL )
  {
    if ( ( _display_parameter_file = fopen(_display_parameter_filename, "r") ) == NULL)
    {
      fprintf(stderr, "Error in opening display parameter file, file %s could not be opened.\n", _display_parameter_filename);
      exit(EXIT_FAILURE);
    }
    fprintf(_log_file,"animate file: %s\n", _display_parameter_filename);
  }

  if ( _lineage_filename != NULL )
  {
    if ( ( _lineage_file = fopen(_lineage_filename, "r") ) == NULL)
    {
      fprintf(stderr, "Error, file %s could not be opened.\n", _lineage_filename);
      exit(EXIT_FAILURE);
    }
    fprintf(_log_file,"lineage file: %s\n", _lineage_filename);
  }

  if (!signal_specified)
  {
    _signal_name = (char *)malloc(sizeof(char));
    _signal_name[0] = 0;
    _signal_index = DEFAULT_SIGNAL; /* by default */
  }

  /* set font and border color to good contrast: black or white */
  _font_color[0] = _background_color[0] + _background_color[1] + _background_color[2] > 1.4 ? 0.0 : 1.0;
  _font_color[1] = _font_color[0];
  _font_color[2] = _font_color[0];
  _border_color[0] = _font_color[0];
  _border_color[1] = _font_color[0];
  _border_color[2] = _font_color[0];

  read_header();
  glutInit(&argc, argv);
  glutInitDisplayMode (GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH | GLUT_STENCIL | GLUT_MULTISAMPLE);

  /* create main view window */
  glutInitWindowSize (_main_window_size[0], _main_window_size[1]);
  glutInitWindowPosition (100, 100);
  _main_window = glutCreateWindow ( "Simuscale - View" );
  init ();

  // Load first frame before first display
  read_file_line();

  glutDisplayFunc(display);
  glutReshapeFunc(reshape);
  glutKeyboardFunc(keyboard);
  glutSpecialFunc(special_keyboard);
  glutMouseFunc(mouse);
  glutMotionFunc(motion);
  glutPassiveMotionFunc(passive_motion);
  create_menu();
  /* end create main view window */

  /* create auxiliary window _aux_window */
  glutInitWindowSize (500, 200);
  glutInitWindowPosition (130 + _main_window_size[0], 100);
  _aux_window = glutCreateWindow("Simuscale - View - Cell");
  glClearColor (_background_color[0],
                _background_color[1],
                _background_color[2],
                1.0f);
  glutDisplayFunc(auxdisplay);
  glutReshapeFunc(reshape);
  glutKeyboardFunc(keyboard);
  glutSpecialFunc(special_keyboard);
  glutMouseFunc(mouse);

  if ( _snapshot_file != NULL ) glutHideWindow();
  /* end - create auxiliary window _aux_window */ 

  /* create statistics window _stats_window */
  glutInitWindowSize (500, 200);
  glutInitWindowPosition (130 + _main_window_size[0], 350);
  _stats_window = glutCreateWindow("Simuscale - View - Stats");
  glClearColor (_background_color[0],
                _background_color[1],
                _background_color[2],
                1.0f);
  glutDisplayFunc(statsdisplay);
  glutReshapeFunc(reshape);
  glutKeyboardFunc(keyboard);
  glutSpecialFunc(special_keyboard);
  glutMouseFunc(mouse);

  if ( _snapshot_file != NULL ) glutHideWindow();
  /* end - create statistics window _stats_window */

  glutSetWindow(_main_window);

  atexit(on_exit);

  glutTimerFunc(0,animate,0); 
  glutMainLoop();

  return 0;
}


void read_min_max_signals()
{
  FILE * minmaxfile = fopen(_norm_filename, "r");
  if (minmaxfile == NULL)
  {
    fprintf(stderr,"No file called %s in current directory.\n", _norm_filename);
    fprintf(stderr,"Using default values for coloring: min=0.0 and max=100.0.\n");
    return;
  }

  int retval;
  char c;
  char str[128];
  retval = fscanf(minmaxfile, "%f", &_max_signal);
  while (retval != 1)
  {
    if (retval == EOF) {printf("File normalization.txt found but premature end of file.\n");return;}

    if (fscanf(minmaxfile, "%s", str) == 1)
    {
      if (str[0]=='#' || (strncmp(str, _signal_name, 127) != 0 && strlen(_signal_name)) ) 
      {
        /* ignore the whole line */
        /* printf("comment line :\n"); */
        do {retval = fscanf(minmaxfile, "%c", &c); /* printf("%c", c); */}  while ((c != '\n') && (retval != EOF));
        if (retval == EOF) {fprintf(_log_file,"File normalization.txt found but premature end of file.\n");return;}
      }
      else if ( strncmp(str, _signal_name, 127) == 0 || strlen(_signal_name) == 0 )
      {
        /* nothing to do */
      }
      else
      {
        /* line starting with a letter, for example */
        fprintf(stderr, "Unknown file format.\n");
        return;
      }
    }
    retval = fscanf(minmaxfile, "%f", &_max_signal);
  }
  _max_signal += 1e-6; /* to avoid clipping for the max value */

  retval = fscanf(minmaxfile, "%f", &_min_signal);
  while (retval != 1)
  {
    if (retval == EOF) {printf("File normalization.txt found but premature end of file.\n");return;}

    if (fscanf(minmaxfile, "%c", &c) == 1)
    {
      if (c=='#')
      {
        /* comment line, ignore the whole line */
        /* printf("comment line :\n"); */
        do {retval = fscanf(minmaxfile, "%c", &c); /* printf("%c", c); */}  while ((c != '\n') && (retval != EOF));
        if (retval == EOF) {printf("File normalization.txt found but premature end of file.\n");return;}
      }
      else
      {
        /* line strating with a letter, for example */
        fprintf(stderr, "Unknown file format.\n");
        return;
      }
    }
    retval = fscanf(minmaxfile, "%f", &_min_signal);
  }

  fclose(minmaxfile);
}

void read_header()
{
  int ind;
  int retval;
  char c;
  char str[128];


  while ( ( retval = fscanf(_traj_file,"%[!$#]", &c ) == 1 ) )
  {
    switch(c)
    {
      case '#' :
        if ( fscanf(_traj_file, " %*d : %*[a-z] %s", str) == 1 )
        {
          /* comment line, check if orientation is printed */
          if ( strncmp("orientation",str,11) == 0 )
          {
            _orientation = 1;
          }
          /* comment line, check if force is printed */
          if ( strncmp("force",str,11) == 0 )
          {
            _force = 1;
          }
          /* comment line, check if velocity is printed */
          if ( strncmp("velocity",str,11) == 0 )
          {
            _velocity = 1;
          }
          /* comment line, check if living status is printed */
          if ( strncmp("status",str,6) == 0 )
          {
            _living_status = 1;
          }
        }
        /* printf("comment line :\n"); */
        break; 
      case '!' :
        /* param line specifying _worldsize */
        fscanf(_traj_file, "%f %f %f \"%[^\"]\"", &(_worldsize.x), &(_worldsize.y), &(_worldsize.z), _length_unit );
        break;
      case '$' :
        /* signal name e.g. $ 8 = SIGNAL_1 */
        fscanf(_traj_file, "%d = %s", &ind, str); 
        if ( _list_signals.size == 0 ) 
        {
          _first_signal_col = ind;
          /* printf("_first_signal_col = %d\n", _first_signal_col); */
        }
        _list_signals.size++;
        if ( _list_signals.size < MAX_NBR_SIGNAL )
        {
          strncpy(_list_signals.name[_list_signals.size-1], str, 127);
        }
        /* first signal found is used as default signal */ 
        if (!_signal_found)
        {
          if ( strncmp(str,_signal_name,127) == 0 )
          {
            _signal_index = ind; 
            _signal_found = 1;
          }
          if ( _signal_index == DEFAULT_SIGNAL )
          {
            _signal_name = (char*)realloc(_signal_name,(strlen(str) + 1)*sizeof(char));
            _signal_index = ind;
            snprintf(_signal_name, strlen(str) + 1, "%s", str);
            _signal_found = 1;
          }
        }
        break;
      default :  
        fprintf(stderr, "Unknown file format.\n");
        exit(EXIT_FAILURE);
    }

    /* go to next line or exit if EOF is reached */
    do {retval = fscanf(_traj_file, "%c", &c); /* printf("%c", c); */ }  while ((c != '\n') && (retval != EOF));
    if (retval == EOF) {exit(0);}
    _line ++;

  }

  if ( _signal_index == SIGNAL_NOT_FOUND )
  {
    fprintf(stderr, "Signal %s not found.\n", _signal_name);
    exit(EXIT_FAILURE);
  }

  if ( _worldsize.x == 0.0 ) 
  {
    fprintf(stderr, "worldsize not found.\n");
    exit(EXIT_FAILURE);
  }

  /* printf("line %d time %f pop %d \n", _line, _time, _popsize);  */
  
}

void read_file_line()
{
  float signal;
  int i, col;
  int colorindex;
  int minclip = 0, maxclip = 0;
  int retval;
  char c;

  if ( _fpos.i >= _fpos.maxsize ) /* increase storage for file position markers */
  {
    _fpos.maxsize *= 2;
    _fpos.pos = (unsigned long *)realloc(_fpos.pos,_fpos.maxsize*sizeof(unsigned long));
  }

  /* save current position in _traj_file */
  _fpos.pos[_fpos.i] = ftell(_traj_file);
  /* printf("file: frame i = %lu, pos = %lu\n",_fpos.i,_fpos.pos[_fpos.i]);  */

  /* look forward to _popsize or EOF */
  retval = fscanf(_traj_file, "%f %d", &_time, &_popsize);

  if (retval == EOF) 
  {
    _end_of_file = 1;
    return;
  }
  else
  {
    _end_of_file = 0;
  }

  /* get back to current position */
  fseek(_traj_file, _fpos.pos[_fpos.i], SEEK_SET);      
  _fpos.i++;

  /* clear _mycells */
  if (_mycells != NULL) free(_mycells);
  _mycells = (cell_sphere *) malloc(_popsize * sizeof(cell_sphere));

  /* read _popsize lines in _traj_file   */
  for ( i = 0; i < _popsize; i++ )  
  {
    retval = fscanf(_traj_file, "%f %d", &_time, &_popsize);
    if ( retval == EOF )
    {
      fprintf(stderr, "Time %f contains only %d rows, but expects %d.\n", _time, i, _popsize);
      exit( EXIT_FAILURE );
    }

    fscanf(_traj_file, "%d %f %f %f", &_mycells[i].id, &_mycells[i].position[0], &_mycells[i].position[1], &_mycells[i].position[2]);

    if ( _orientation )
    {
      fscanf(_traj_file, "%f %f %f", &_mycells[i].orientation[0], &_mycells[i].orientation[1], &_mycells[i].orientation[2]); /* read and discard orientation vector */
    }

    if ( _force )
    {
      fscanf(_traj_file, "%f %f %f", &_mycells[i].force[0], &_mycells[i].force[1], &_mycells[i].force[2]); /* read and discard force vector */
    }

    if ( _velocity )
    {
      fscanf(_traj_file, "%f %f %f", &_mycells[i].velocity[0], &_mycells[i].velocity[1], &_mycells[i].velocity[2]); /* read and discard velocity vector */
    }

    fscanf(_traj_file, "%f", &_mycells[i].radius);
    fscanf(_traj_file, "%s",  _mycells[i].cell_tag);

    if ( _living_status ) /* column with living status present; read from it */
    {
      fscanf(_traj_file, "%*[ ] %c",  &_mycells[i].living_status);
    }
    else /* living status not present; all cell are presumed alive */
    {
      _mycells[i].living_status = 'A';
    } 

    col = 9 + 3*_orientation + 3*_force + 3*_velocity + _living_status;
    while ( col < _first_signal_col )
    {
      retval = fscanf(_traj_file, "%*s"); /* skip to next entry */
      if ( retval == EOF )
      {
        fprintf(stderr, "File contains only %d column, but needs %d.\n", col, _signal_index);
        exit( EXIT_FAILURE );
      }
      ++col;
    }

    while ( col < _signal_index )
    {
      retval = fscanf(_traj_file, "%f", &signal); /* skip to next entry */
      if ( retval == EOF )
      {
        fprintf(stderr, "File contains only %d column, but needs %d.\n", col, _signal_index);
        exit( EXIT_FAILURE );
      }
      ++col;
    }
    fscanf(_traj_file, "%f", &signal);

    /* get mother cell */
    _mycells[i].lineage = get_lineage(_mycells[i].id);

    /* compute color from signal */
    /* Check if signal is NaN (used for cells that don't have this signal type) */
    if (isnan(signal))
    {
      /* Make cell transparent for signals it doesn't participate in */
      _mycells[i].color[0] = 0.4;
      _mycells[i].color[1] = 0.4;
      _mycells[i].color[2] = 0.4;
      _mycells[i].color[3] = 0.4;  /* fully transparent */
    }
    else
    {
      /* Normal colormap assignment for valid signal values */
      colorindex = (int) _colormap_len * (signal - _min_signal) / (_max_signal - _min_signal);

      if (colorindex < 0) 
      {
        colorindex = 0; 
        minclip = 1;
      }
      else if (colorindex > _colormap_len - 1) 
      {
        colorindex = _colormap_len - 1; 
        maxclip = 1;
      }

      _mycells[i].color[0] = _colormap[colorindex].c[0];
      _mycells[i].color[1] = _colormap[colorindex].c[1];
      _mycells[i].color[2] = _colormap[colorindex].c[2];
      _mycells[i].color[3] = _alpha;
      
      if ( _color_model == COLOR_MODEL_FLAT )
      {
        darken(_mycells[i].color, 1.0); /* make color about a third 
                                         * darker to compensate for the
                                         * diffusion color in the lighting
                                         * model: see diffuse_light */
      }
    }

    _mycells[i].clip     = ( signal > _max_clip_signal ) || ( signal < _min_clip_signal );

    _mycells[i].signal = signal;

    do {retval = fscanf(_traj_file, "%c", &c); /* printf("%c", c); */}  while ((c != '\n') && (retval != EOF));
    if (retval == EOF) {exit(0);}

  } 
  /* end read _popsize lines in _traj_file */

  if (minclip) fprintf(_log_file,"Warning, at least one cell had min color clipping at t=%f.\n", _time);
  if (maxclip) fprintf(_log_file,"Warning, at least one cell had max color clipping at t=%f.\n", _time);

  _line++;
}

void read_display_parameters() 
{
  int retval;
  static int display_version = 0;
  struct display_pars disp = { 
    .fval = { _azim_angle, 
              _inc_angle, 
              _twist_angle, 
              _distance, 
              _pan_x, 
              _pan_y, 
              _x_plane, 
              _y_plane, 
              _z_plane }, 
    .sval = "", 
    .ival = (unsigned int)_draw_guides,
    .pause = 0 };
  strncpy(disp.sval,_signal_name,127);

  /* on first read, scan for animate file format
   * if first row is "version X", use version X,
   * otherwise use version 1
   */
  if ( display_version == 0 ) /* first time reading */
  {
    retval = fscanf(_display_parameter_file, "version %d\n", &display_version); 
    if ( retval == 0 )
    {
      display_version = 1; /* default value */
    }
  }

  if ( display_version == 2 )
  {
    display_par2(_display_parameter_file, _fpos.i - 1, &disp);
  }
  else if ( display_version == 1 )
  {
    display_par1(_display_parameter_file, _fpos.i - 1, &disp);
  }

  _azim_angle = disp.fval[0];
  _inc_angle  = disp.fval[1];
  _twist_angle = disp.fval[2];
  _distance   = disp.fval[3];
  _pan_x = disp.fval[4];
  _pan_y = disp.fval[5];
  _x_plane = disp.fval[6];
  _y_plane = disp.fval[7];
  _z_plane = disp.fval[8];
  _signal_name = strndup(disp.sval,127);
  _draw_guides = (int)disp.ival;
  _pause = disp.pause;
  update_signal(SIGNAL_NAME);
  read_min_max_signals();

  /* printf("--%s line %d  step=%lu, pause=%d\n", __func__, __LINE__, _fpos.i -1, _pause); */
  
}

void animate(int)
{
  
  /* printf("step %ld\n",_fpos.i); */
  if ( _fpos.i  == (unsigned long)_initial_step )
  {
    glutSetWindow(_main_window);
    glutPostRedisplay();
    glFinish();
    /* if _snapshot_file, save to bitmap and exit */
    if ( _snapshot_file != NULL )
    {
      printf("snapshot: step %ld\n",_fpos.i - 1);
      save_image();
      print_snapshot();
      exit( EXIT_SUCCESS );
    }
  }

  while ( _skip_frame-- )
  {
    read_file_line();
  }
  _skip_frame = 1;

  if ( _display_parameter_file != NULL )
  {
    /* read display parameters
     * from _display_parameter_filename 
     * and update display parameters */
    read_display_parameters(); 
  }

  if ( _pause == 0 && _end_of_file == 0 ) 
  {
    glutTimerFunc(1.0/_speed, animate, 0);
  }
  if ( _pause == 1 && _display_parameter_file != NULL )
  {
    standstill();
    glutTimerFunc(1.0/_speed, animate, 0);
  }


  /* glutPostRedisplay(); */
  glutSetWindow(_main_window);
  glutPostRedisplay();
  glutSetWindow(_aux_window);
  glutPostRedisplay();
  glutSetWindow(_stats_window);
  glutPostRedisplay();

  if (_record_video) {
      // Set the correct window before capturing
      switch(_record_window) {
          case 0: glutSetWindow(_main_window); break;
          case 1: glutSetWindow(_aux_window); break;
          case 2: glutSetWindow(_stats_window); break;
      }
      
      glFinish(); // Make sure rendering is complete
      add_frame_to_video();
  }

  // Check if we've reached end of file and stop recording
  if (_record_video && _end_of_file) {
      printf("Animation ended. Finalizing video...\n");
      cleanup_opencv_video();
      _record_video = 0;
  }


}
// 

void draw_neighborhood_sphere(int cell_index) {
    if (!_show_neighborhood_sphere) return;

    const double nbd_factor = 1.0;
    const double sphere_radius = _mycells[cell_index].radius * (1.0 + nbd_factor);

    const GLfloat cx = _mycells[cell_index].position[0];
    const GLfloat cy = _mycells[cell_index].position[1];
    const GLfloat cz = _mycells[cell_index].position[2];

    glPushMatrix();
    // Save/restore enables, line/poly modes, color buffer (blend), depth, stencil
    glPushAttrib(GL_ENABLE_BIT | GL_LINE_BIT | GL_COLOR_BUFFER_BIT | GL_POLYGON_BIT |
                 GL_DEPTH_BUFFER_BIT | GL_STENCIL_BUFFER_BIT);

    glDisable(GL_LIGHTING);
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

    // Don’t let stencil test from cells affect the wireframe
    glStencilFunc(GL_ALWAYS, STENCIL_LEGEND, -1);

    glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
    glLineWidth(2.0f);
    glColor4f(0.2f, 0.2f, 0.2f, 0.85f); // Dark gray for visibility against all colors
    // glColor4f(1.00f, 1.00f, 1.00f, 0.90f); // White for visibility against all colors (was magenta)

    glTranslatef(cx, cy, cz);
    glutWireSphere(sphere_radius, 20, 20);

    glPopAttrib();   // Restores GL_BLEND to previous (keeps plane transparency intact)
    glPopMatrix();
}


void init(void)
{

  if ( _colormap_name < NB_COLORMAPS )
  {
    unsigned long *c = (unsigned long *)malloc(5*sizeof(unsigned long)); 
    printf("colormap %s (%d)\n", _list_colormaps.names[_colormap_name], _colormap_name);
    for ( int i = 0; i < 5; i++ )
    { 
      c[i] = _list_colormaps.colors[_colormap_name][i];
    }
    _colormap_len = NB_COLORS;
    make_colormap_array(c, 5, &_colormap, _colormap_len );
    free(c);
  }
  else
  {
    /* custom */
  }

  read_min_max_signals();

  if ( _lineage_file != NULL )
  {
    read_lineage_file();
  }

  glPixelStorei(GL_UNPACK_ALIGNMENT, 1);

  glClearColor (_background_color[0],
                _background_color[1],
                _background_color[2],
                1.0f);




   /* Set up near_clip and far_clip so that ( far_clip - near_clip ) > maxObjectSize, */
   /* and determine the viewing distance (adjust for zooming) */
   _near_clip = 0.1; /* 0.1*_worldsize.y; */
   _far_clip = _near_clip + 4.0*_worldsize.y;
   reset_view();

   object_material();

   /* display list to store the geometry of a sphere */
   _cell_dlist = glGenLists(8);
   glNewList(_cell_dlist + GEOM_NORMAL, GL_COMPILE);
   glutSolidSphere(1.0, _s_res, _s_res);
   glEndList();

   glNewList(_cell_dlist + GEOM_SECTION, GL_COMPILE);
   sphere_section( .0f );
   glEndList();

   glNewList(_cell_dlist + GEOM_SECTION_X_14, GL_COMPILE);
   glPushMatrix(); 
   glRotatef( 90, 0.0f, -1.0f, 0.0f);
   sphere_section( 0.50f );
   glPopMatrix(); 
   glEndList();

   glNewList(_cell_dlist + GEOM_SECTION_X_HALF, GL_COMPILE);
   glPushMatrix(); 
   glRotatef( 90, 0.0f, -1.0f, 0.0f);
   sphere_section( .0f );
   glPopMatrix(); 
   glEndList();

   glNewList(_cell_dlist + GEOM_SECTION_X_34, GL_COMPILE);
   glPushMatrix(); 
   glRotatef( 90, 0.0f, -1.0f, 0.0f);
   sphere_section( -.50f );
   glPopMatrix(); 
   glEndList();

   glNewList(_cell_dlist + GEOM_SECTION_Y_14, GL_COMPILE);
   glPushMatrix(); 
   glRotatef( 90, -1.0f, 0.0f, 0.0f);
   sphere_section( 0.50f );
   glPopMatrix(); 
   glEndList();

   glNewList(_cell_dlist + GEOM_SECTION_Y_HALF, GL_COMPILE);
   glPushMatrix(); 
   glRotatef( 90, -1.0f, 0.0f, 0.0f);
   sphere_section( .0f );
   glPopMatrix(); 
   glEndList();

   glNewList(_cell_dlist + GEOM_SECTION_Y_34, GL_COMPILE);
   glPushMatrix(); 
   glRotatef( 90, -1.0f, 0.0f, 0.0f);
   sphere_section( -0.50f );
   glPopMatrix(); 
   glEndList();

   glNewList(_cell_dlist + GEOM_SECTION_Z_14, GL_COMPILE);
   glPushMatrix(); 
   glRotatef( 180, 1.0f, 0.0f, 0.0f);
   sphere_section( 0.50f );
   glPopMatrix(); 
   glEndList();

   glNewList(_cell_dlist + GEOM_SECTION_Z_HALF, GL_COMPILE);
   glPushMatrix(); 
   glRotatef( 180, 1.0f, 0.0f, 0.0f);
   sphere_section( .0f );
   glPopMatrix(); 
   glEndList();

   glNewList(_cell_dlist + GEOM_SECTION_Z_34, GL_COMPILE);
   glPushMatrix(); 
   glRotatef( 180, 1.0f, 0.0f, 0.0f);
   sphere_section( -0.50f );
   glPopMatrix(); 
   glEndList();

   arrow_geometry();
   scalebar_geometry();
   dead_cell_geometry();
}

void display(void)
{
  char legend_string[128];
  char legend_string2[128];
  char single_cell_legend_string[256];
  char plane_position_string[24];
  const char axis_name[3][7] = {"[x]", "[y]", "[z]"};
  int colorwidth = 520/_colormap_len - 1;
  static GLfloat menu_alpha = 0.0;

  /***** Clear color buffer, depth buffer, stencil buffer *****/
  /* depth buffer: used for single cell selection
   * stencil buffer: used to identify the background STENCIL_BACKGROUND 
   *                 non-selected cells STENCIL_CELL
   *                 selected cell STENCIL_TRACK
   *                 world margins STENCIL_BORDER
   *                 legends STENCIL_LEGEND
   */
  glClearStencil(STENCIL_BACKGROUND); /* background default stencil value: 0 */ 
  glClear (GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT | GL_STENCIL_BUFFER_BIT );


  /*************** Draw the legend texts and colorbar **********************/
  
  /* Draw legend texts and box (unless --no-annotations or --colorbar-only is set) */
  if ( !_no_annotations && !_colorbar_only )
  {
    glStencilFunc(GL_ALWAYS, STENCIL_LEGEND, -1); /* set stencil to STENCIL_LEGEND for legends */
    if (_minimal_annotations) {
      snprintf(legend_string, 128, "t=%.2f N=%d", _time, _popsize);
      glColor3fv (_font_color);
      glWindowPos2i(15, 30);
      for (size_t i = 0; i < strlen(legend_string); i++) glutBitmapCharacter(GLUT_BITMAP_HELVETICA_18, legend_string[i]);
    } else {
      snprintf(plane_position_string, 24, " (x>%4.1f,y<%4.1f,z>%4.1f)", _x_plane, _y_plane, _z_plane); 
      snprintf(legend_string, 128, "t=%.2f s=%ld N=%d%6s%7s", _time, _fpos.i - 1, _popsize, \
          _clipcell ? " clip" : "", _pause ? " paused" : "");
      snprintf(legend_string2, 128, "%s%s %s", _move_plane ? axis_name[_pos_plane] : "hide ",
          plane_position_string, _length_unit); 
      glColor3fv (_font_color);
      glWindowPos2i(15, 30); /* also sets current raster color to _font_color */
      for (size_t i = 0; i < strlen(legend_string); i++) glutBitmapCharacter(GLUT_BITMAP_HELVETICA_18, legend_string[i]);

      glWindowPos2i(15, 13); /* also sets current raster color to _font_color */
      for (size_t i = 0; i < strlen(legend_string2); i++) glutBitmapCharacter(GLUT_BITMAP_HELVETICA_18, legend_string2[i]);

      /*** if _length_unit is known, print it ***/
      if ( strlen(_length_unit) > 0 )
      {
        glBitmap(8,9,0.0,2.0,7.0,0.0,_raster_mu); /* bitmap mu */
        glutBitmapCharacter(GLUT_BITMAP_HELVETICA_18, 'm');
      }
      
      /****************** draw legend box *************************/
      glColor3f(0.8 - 0.7*_font_color[0],0.8 - 0.7*_font_color[1],0.8 - 0.7*_font_color[2]);
      glWindowPos2i(10, 9); 
      glBitmap(8, 16, 0.0, 0.0, 0, 16, _raster_lower_left);
      glBitmap(8, 16, 0.0, 0.0, 8, -16, _raster_upper_left);
      for(auto i = 0; i < 26; i++)
      {
        glBitmap(8, 16, 0.0, 0.0, 0, 16, _raster_rect);
        glBitmap(8, 16, 0.0, 0.0, 8, -16, _raster_rect);
      }
      glBitmap(8, 16, 0.0, 0.0, 0, 16, _raster_lower_right);
      glBitmap(8, 16, 0.0, 0.0, 8, 0, _raster_upper_right);
    }
  }

  /* Draw colorbar (unless --no-annotations or --minimal-annotations is set) */
  if ( !_no_annotations && !_minimal_annotations )
  {
    glStencilFunc(GL_ALWAYS, STENCIL_LEGEND, -1); /* set stencil to STENCIL_LEGEND for colorbar */

    /*************** Draw the signal colorbar *******************/
    for(auto i = 0; i < _colormap_len; i++)
    {
      glColor3fv(_colormap[i].c);
      glWindowPos2i(glutGet(GLUT_WINDOW_WIDTH) - 524 + (colorwidth + 1)*i, 10); 
      for ( int i = 0; i < colorwidth/4; i++ )
      {
        glBitmap(4, 16, 0.0, 0.0, 4, 0, _raster_rect);
      }
      glBitmap(colorwidth % 4, 16, 0.0, 0.0, 4, 0, _raster_rect);
    }

    glColor3fv (_font_color);
    snprintf(legend_string, 128, "%.2g", _min_signal);
    glWindowPos2i(glutGet(GLUT_WINDOW_WIDTH) - 520, 30); /* also sets current raster color to _font_color */
    for (size_t i = 0; i < strlen(legend_string); i++)
      glutBitmapCharacter(GLUT_BITMAP_TIMES_ROMAN_24, legend_string[i]);

    glColor3fv (_font_color);
    /* position signal name to the right of the min signal 
     * padding by 15 pixels for each bitmap character in min signal + extra 40px for subscripts
     * */
    glWindowPos2i(glutGet(GLUT_WINDOW_WIDTH) - 520 + 15*strlen(legend_string) +10, 32); /* also sets current raster color to _font_color */
    snprintf(legend_string, 128, "%s", _signal_name);
    if ( strlen(legend_string) > MAX_SIG_STRLEN ) /* print first three chars ... last (MAX_SIG_STRLEN - 6) */
    {
      for (auto i = 0; i < 3; i++)
        glutBitmapCharacter(GLUT_BITMAP_TIMES_ROMAN_24, legend_string[i]);
      for (auto i = 3; i < 6; i++)
        glutBitmapCharacter(GLUT_BITMAP_TIMES_ROMAN_24, '.');
      for (auto i = strlen(legend_string) - (MAX_SIG_STRLEN - 6); i < strlen(legend_string); i++)
        glutBitmapCharacter(GLUT_BITMAP_TIMES_ROMAN_24, legend_string[i]);
    }
    else
    {
      /* Render signal name with subscript formatting */
      renderSignalNameWithSubscript(legend_string, GLUT_BITMAP_TIMES_ROMAN_24, GLUT_BITMAP_HELVETICA_18);
    }

    glColor3fv (_font_color);
    snprintf(legend_string, 128, "%.2g", _max_signal);
    /* position max signal to the right of the colorbar
     * padding by 8 pixels for each bitmap character in max signal  
     * */
    glWindowPos2i(glutGet(GLUT_WINDOW_WIDTH) - 12*strlen(legend_string)-5, 30); /* also sets current raster color to _font_color */
    for (size_t i = 0; i < strlen(legend_string); i++)
      glutBitmapCharacter(GLUT_BITMAP_TIMES_ROMAN_24, legend_string[i]);
  } /* end of colorbar drawing */

  /*********************** Set the view **************************/
  glLoadIdentity ();             /* clear the matrix */
  /* viewing transformation  */
  /* gluLookAt (-1.5*_worldsize, -1.5*_worldsize, 1.5*_worldsize, _worldsize, _worldsize, 0.0, _worldsize, _worldsize, 1.5*_worldsize); */
  polar_view( _distance, _azim_angle, _inc_angle, _twist_angle );

  /******************** Draw the world borders ********************/
  glStencilFunc(GL_ALWAYS, STENCIL_BORDER, -1); /* set stencil to STENCIL_BORDER for world borders */
  glColor3fv (_border_color); // Black or White
  glLineWidth(1.0);
  glNormal3d(0,-1,0);
  glPolygonMode(GL_FRONT_AND_BACK, GL_LINE); // Mesh

  float world_size[3] = {_worldsize.x, _worldsize.y, _worldsize.z};
  float world_margins[3] = {0.05f*_worldsize.x, 0.05f*_worldsize.y, 0.05f*_worldsize.z}; 
  float total_world_size[3] = {world_size[0] + 2 * world_margins[0],
                               world_size[1] + 2 * world_margins[1],
                               world_size[2] + 2 * world_margins[2]};

  if ( _draw_guides & GUIDE_BORDER )
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

  /*********** draw scalebar ************/ 
  if ( _draw_guides & GUIDE_SCALEBAR )
  {
    glStencilFunc(GL_ALWAYS, STENCIL_SCALEBAR, -1); /* set stencil to STENCIL_SCALEBAR for world borders */
    glPushMatrix(); 
    glTranslatef( _scale_pos[0], _scale_pos[1], _scale_pos[2] ); 
    glCallList(_cell_dlist + GEOM_SCALEBAR ); 
    glPopMatrix(); 
  }

  /*********** set highlight hoovering properties ************/ 
  for ( size_t i = 0; i < (size_t)_popsize ; i++ ) 
  { 
    _mycells[i].highlight = 0;
  }
  if ( _single_cell_hoovering == HOOV_ON )  /* && _mycells[i].id == _single_cell_hoovered ) */
  {
    /************ reset highlight hoovering properties *******/
    for ( size_t i = 0; i < (size_t)_popsize ; i++ ) 
    { 
      if ( _mycells[i].id == _single_cell_hoovered )
      {
        _mycells[i].highlight = 1;
      
        /* ancestors: 2 */
        int index = i; ;
        while ( ( index = get_cell_index(_mycells[index].lineage) ) < _popsize )
        {
          _mycells[index].highlight = 2;
        }

        /* progeny: 3 */
        for ( size_t j = 0; j < (size_t)_popsize ; j++ ) 
        { 
          int index = get_cell_index(_mycells[j].lineage); /* mother cell index */
          if ( index < _popsize && (_mycells[index].highlight % 2) == 1 )
          {
            _mycells[j].highlight = 3; 
          }
        }

      }
    }

  }

  /******************* Draw the cells ***************************/
  for (int i = 0; i < _popsize; i++)
  {
    /* if _clipcell = true, do not draw color clipped cells */
    if ( _mycells[i].clip && _clipcell ) /* cell is clipped, do not draw */
    {
      continue;
    }

    /* hide (cut plane) cells that are above x-plane, below y-plane
     * and above z-plane */
    if ( ( _mycells[i].position[0] >= _x_plane + 0.5*_mycells[i].radius ) &&
         ( _mycells[i].position[1] <= _y_plane - 0.5*_mycells[i].radius ) &&
         ( _mycells[i].position[2] >= _z_plane + 0.5*_mycells[i].radius ) )
    {
      continue; /* do not draw the cell */
    }

    /* Set stencil and single cell legends */
    glStencilFunc(GL_ALWAYS, STENCIL_LEGEND, -1); /* set stencil to STENCIL_LEGEND for legends */

    /******** set stencil index for single cell tracking ********/
    if ( !_no_annotations && !_colorbar_only && _single_cell_tracking == 1 && _mycells[i].id == _single_cell )
    {

      /*********************** Draw tracked cell legend ***********************/
      /* Render prefix first */
      char prefix[128];
      snprintf(prefix, 128, "ID: %d%c, tag: %s, ", _single_cell, _mycells[i].living_status, _mycells[i].cell_tag);
      glColor3fv (_font_color);
      glWindowPos2i(15, glutGet(GLUT_WINDOW_HEIGHT) - 15);
      for (size_t j = 0; j < strlen(prefix); j++) 
        glutBitmapCharacter(GLUT_BITMAP_HELVETICA_12, prefix[j]);
      /* Render signal name with subscript */
      renderSignalNameWithSubscript(_signal_name, GLUT_BITMAP_HELVETICA_12, GLUT_BITMAP_HELVETICA_10);
      /* Render value */
      snprintf(prefix, 128, ": %.4g", _mycells[i].signal);
      GLfloat raster_pos[4];
      glGetFloatv(GL_CURRENT_RASTER_POSITION, raster_pos);
      glWindowPos2f(raster_pos[0], raster_pos[1]);
      for (size_t j = 0; j < strlen(prefix); j++) 
        glutBitmapCharacter(GLUT_BITMAP_HELVETICA_12, prefix[j]);
    }
    /*********************** Draw hoovered single cell legend ***********************/
    else if ( !_no_annotations && !_colorbar_only && _single_cell_hoovering == HOOV_ON && _mycells[i].id == _single_cell_hoovered )
    {

      /*********************** Draw hoovered cell legend ***********************/
      snprintf(single_cell_legend_string, 256, "ID: %d%c, tag: %s, mother: %d",
          _single_cell_hoovered, _mycells[i].living_status, _mycells[i].cell_tag, _mycells[i].lineage);
      glColor3fv (_font_color);
      glWindowPos2i(glutGet(GLUT_WINDOW_HEIGHT) - strlen(single_cell_legend_string)*6, glutGet(GLUT_WINDOW_HEIGHT) - 15); /* also sets current raster color to white */
      for (size_t i = 0; i < strlen(single_cell_legend_string); i++) glutBitmapCharacter(GLUT_BITMAP_HELVETICA_10, single_cell_legend_string[i]);
      /* Render signal name with subscript and value */
      glColor3fv (_font_color);
      glWindowPos2i(glutGet(GLUT_WINDOW_HEIGHT) - 100, glutGet(GLUT_WINDOW_HEIGHT) - 30);
      renderSignalNameWithSubscript(_signal_name, GLUT_BITMAP_HELVETICA_10, GLUT_BITMAP_8_BY_13);
      snprintf(single_cell_legend_string, 256, ": %.4g", _mycells[i].signal);
      GLfloat raster_pos2[4];
      glGetFloatv(GL_CURRENT_RASTER_POSITION, raster_pos2);
      glWindowPos2f(raster_pos2[0], raster_pos2[1]);
      for (size_t i = 0; i < strlen(single_cell_legend_string); i++) 
        glutBitmapCharacter(GLUT_BITMAP_HELVETICA_10, single_cell_legend_string[i]);
    }

    if ( _draw_elements & DRAW_CELLS ) {
        draw_cell(i);
        
        // ADD THE NEIGHBORHOOD SPHERE HERE - inside the main loop
        if (_single_cell_tracking == 1 && _mycells[i].id == _single_cell) {
            draw_neighborhood_sphere(i);
        }
    }

    /**** draw cell i *****/
    // if ( _draw_elements & DRAW_CELLS ) draw_cell(i);
    if ( _draw_elements & DRAW_FORCES ) draw_force(i);
    if ( _draw_elements & DRAW_VELOCITIES ) draw_velocity(i);

  }

  /******************* Draw positioning plane *******************/
  if ( _draw_guides & GUIDE_VIEWCUT )
  {
    glStencilFunc(GL_GEQUAL, STENCIL_CELL, 0xFF); /* not cell */
    draw_plane();
  }

  /******************** Draw the menu 'button' ***********************/
  if ( _menu_hoovering )
  {
    menu_alpha += 0.05;
  }
  else
  {
    menu_alpha *= 0.8;
  }
  glColor4f (_font_color[0],_font_color[1],_font_color[2],menu_alpha);
  glWindowPos2i(10, glutGet(GLUT_WINDOW_HEIGHT) - 20); /* also sets current raster color to _font_color */
  for (size_t i = 0; i < strlen(menu_txt); i++)
    glutBitmapCharacter(GLUT_BITMAP_HELVETICA_12, menu_txt[i]);
  glColor4f(0.8 - 0.7*_font_color[0],0.8 - 0.7*_font_color[0],0.8 - 0.7*_font_color[0],menu_alpha);
  glWindowPos2i(0, glutGet(GLUT_WINDOW_HEIGHT) - 32); 
  for(auto i = 0; i < 6; i++)
  {
    glBitmap(8, 16, 0.0, 0.0, 0, 16, _raster_rect);
    glBitmap(8, 16, 0.0, 0.0, 8, -16, _raster_rect);
  }
  glBitmap(8, 16, 0.0, 0.0, 0, 16, _raster_lower_right);
  glBitmap(8, 16, 0.0, 0.0, 8, -16, _raster_rect);

  /******************** end of display ***********************/


  glutSwapBuffers();

}

void draw_cell(GLint i) 
{
  static const GLfloat niche_color[] = {0.5,0.5,0.5,0.1};
  static const GLfloat orientation_color[] = {1.0,1.0,1.0,0.8};
  static const GLfloat dying_color[] = {0.4,0.2,0.3,0.1};
  static const GLfloat flat_contour_color[] = {0.0, 0.1, 0.2, 1.0};
  GLuint is_dying = 0;
  GLfloat radius;
  int cell_loc = locate_cell(i);
  int full_cell_display = 0;
  
  
  if (_mycells[i].position[2] <= 0.0 && _using_niche ) {
    glColor4fv(niche_color);
  }   
  /************** Highlight tagged cell with a different color *************************/
  else if ( _tagged_cell_highlighting && 
            _tags != NULL && 
            strncmp(_mycells[i].cell_tag,_tags,32) == 0 )
  {
    /* highlight color */
    glColor4fv(_tag_color);
  }
  else if ( _mycells[i].living_status == 'D' ) 
  {
    glColor3fv(dying_color);
    is_dying = 1;
  }
  else {
    glColor4fv(_mycells[i].color);
  }



  /******************* Draw the cell *******************/

  glPushMatrix(); /* save V */

  /* go to cell position */
  glTranslatef(_mycells[i].position[0], _mycells[i].position[1], _mycells[i].position[2]); /* current matrix is VT  */

  /* scale the cell with radius radius */
  radius = _mycells[i].radius; 
  if ( _color_model == COLOR_MODEL_FLAT ) 
  {
    radius = _mycells[i].radius/1.15; 
  }
  glScalef(radius, radius, radius );  

  /******** set stencil index for single cells ********/
  if ( _single_cell_tracking == 1 && _mycells[i].id == _single_cell )
  {
    glStencilFunc(GL_ALWAYS, STENCIL_TRACK, -1); 
  }
  else if ( _single_cell_hoovering == HOOV_ON && _mycells[i].id == _single_cell_hoovered )
  {
    glStencilFunc(GL_ALWAYS, STENCIL_HOOVER, -1); 
  }
  else if ( _tags != NULL && strncmp(_mycells[i].cell_tag,_tags,32) == 0 )
  {
    glStencilFunc(GL_ALWAYS, STENCIL_TAG, -1);
  }
  else
  {
    glStencilFunc(GL_ALWAYS, STENCIL_CELL, -1); 
  }

  /* draw the cell with its color */

  switch( cell_loc ) 
  {
    case   0: /* no draw */
      break;
    case   1: /* Z_14 */
    case  11: /* xYZ Z_14 */
    case 101: /* XyZ Z_14 */
    case 111: /* XYZ Z_14 */
      glCallList(_cell_dlist + GEOM_SECTION_Z_14 ); 
      break;
    case   2: /* Z_HALF */
    case  12:
    case 102:
    case 112:
    case  22: /* xYZ Z_HALF */
    case 122: /* xYZ Z_HALF */
    case 202: /* XyZ Z_HALF */
    case 212: /* XyZ Z_HALF */
    case 222: /* XYZ Z_HALF */
      glCallList(_cell_dlist + GEOM_SECTION_Z_HALF ); 
      break;
    case   3: /* Z_34 */
    case  13:
    case  23:
    case 103:
    case 113:
    case 123:
    case 203:
    case 213:
    case 223:
      glCallList(_cell_dlist + GEOM_SECTION_Z_34 ); 
      break;
    case  10: /* Y_14 */
    case 110: /* XYz Y_14 */
      glCallList(_cell_dlist + GEOM_SECTION_Y_14 ); 
      break;
    case  20: /* Y_HALF */
    case  21:
    case 120:
    case 121:
    case 220: /* XYz Y_HALF */
    case 221: /* XYz Y_HALF */
      glCallList(_cell_dlist + GEOM_SECTION_Y_HALF ); 
      break;
    case  30: /* Y_34 */
    case  31:
    case  32:
    case 130:
    case 131:
    case 132:
    case 230:
    case 231:
    case 232:
      break;
    case 100: /* X_14 */
      glCallList(_cell_dlist + GEOM_SECTION_X_14 ); 
      break;
    case 200: /* X_HALF */
    case 201:
    case 210:
    case 211:
      glCallList(_cell_dlist + GEOM_SECTION_X_HALF ); 
      break;
    case 300: /* X_34 */
    case 301:  
    case 302:  
    case 310:  
    case 311:  
    case 312:  
    case 320:  
    case 321:  
    case 322:  
      glCallList(_cell_dlist + GEOM_SECTION_X_34 ); 
      break;
    case  33: /* FULL */
    case 133: /* FULL */
    case 233: /* FULL */
    case 303: /* FULL */
    case 313: /* FULL */
    case 323: /* FULL */
    case 330: /* FULL */
    case 331: /* FULL */
    case 332: /* FULL */
    case 333: /* FULL */
    default: /* FULL any loc with 4: none of the above */
      glCallList(_cell_dlist + GEOM_NORMAL + is_dying*GEOM_DEAD); 
      full_cell_display = 1;
      break;
  }

  if ( _orientation == 1 && !is_dying && full_cell_display )
  {
    /* Add a small bump along cell orientation */
    glColor4fv(orientation_color);
    glTranslatef(0.7*_mycells[i].orientation[0], 0.7*_mycells[i].orientation[1], 0.7*_mycells[i].orientation[2]); /* current matrix is VT  */
    glScalef(0.5,0.5,0.5);
    glCallList(_cell_dlist + GEOM_NORMAL);
  }

  glPopMatrix(); /* current matrix is restored to V */

  /**************** contour if COLOR_MODEL_FLAT mode *************/
  if ( _color_model == COLOR_MODEL_FLAT ) 
  {
    /********** Draw a flat_contour_color-color circle around the cell **********/
    glPushMatrix(); /* save V */
    glDisable(GL_LIGHTING);
    glStencilFunc(GL_GREATER, STENCIL_CELL, 0xFF); /* avoid drawing on cells */
    glColor4fv(flat_contour_color);
    glTranslatef(_mycells[i].position[0], _mycells[i].position[1], _mycells[i].position[2]); /* current matrix is VT  */
    glScalef(_mycells[i].radius, _mycells[i].radius, _mycells[i].radius ); /* current matrix is VTS */
    glCallList(_cell_dlist + is_dying*GEOM_DEAD);
    glEnable(GL_LIGHTING);
    glPopMatrix(); /* current matrix is restored to V */
  }

  /**************** highlight cell if tracking enabled *************/
  if ( _single_cell_tracking == 1 && _mycells[i].id == _single_cell )
  {
    /************** Draw a yellow circle around tracked cell ****************/
    glStencilFunc(GL_NOTEQUAL, STENCIL_TRACK, 0xFF); /* tracked cell stencil */
    glPushMatrix(); /* save V */
    glDisable(GL_LIGHTING);
    glTranslatef(_mycells[i].position[0], _mycells[i].position[1], _mycells[i].position[2]); /* current matrix is VT  */
    glScalef(_mycells[i].radius, _mycells[i].radius, _mycells[i].radius ); /* current matrix is VTS */

    /* highlight color 1 */
    glColor3fv(_highlight_color[0]);
    glutSolidSphere(1.2, 2*_s_res, 2*_s_res);

    /* highlight color 4 */
    glColor3fv(_highlight_color[3]);
    glutSolidSphere(1.3, 2*_s_res, 2*_s_res);

    glPopMatrix(); /* current matrix is restored to V */
    glEnable(GL_LIGHTING);
  }
  else if ( _mycells[i].highlight > 0 )
  {
    switch(_mycells[i].highlight)
    {
      case 1:
      /* Draw a yellow circle around hoovered cell */
      hoover_cell(i,0);
      break;
      case 2: 
      /* Draw faint yellow circle around ancestors */ 
      hoover_cell(i,1);
      break;
      case 3:
      /* Draw faint green circle around progeny */
      hoover_cell(i,2);
      break;
    }
  }


}

int locate_cell(GLint i)
{
  /* cell position */
  GLfloat x = _mycells[i].position[0];
  GLfloat y = _mycells[i].position[1];
  GLfloat z = _mycells[i].position[2];
  
  /* x plane positions x4 < x3 < x2 < x1 */
  GLfloat x1 = _x_plane + 0.5*_mycells[i].radius;
  GLfloat x2 = _x_plane;
  GLfloat x3 = _x_plane - 0.5*_mycells[i].radius;
  GLfloat x4 = _x_plane -     _mycells[i].radius;

  /* y plane positions y4 > y3 > y2 > y1 */
  GLfloat y1 = _y_plane - 0.5*_mycells[i].radius;
  GLfloat y2 = _y_plane;
  GLfloat y3 = _y_plane + 0.5*_mycells[i].radius;
  GLfloat y4 = _y_plane +     _mycells[i].radius;
  
  /* z plane positions z4 < z3 < z2 < z1 */
  GLfloat z1 = _z_plane + 0.5*_mycells[i].radius;
  GLfloat z2 = _z_plane;
  GLfloat z3 = _z_plane - 0.5*_mycells[i].radius;
  GLfloat z4 = _z_plane -     _mycells[i].radius;

  /* return value: integer p 
   * p = XYZ, X,Y,Z in {0,1,2,3,4} */
  int p = 0;

  if ( x <= x1 && x > x2 ) p += 100;
  if ( x <= x2 && x > x3 ) p += 200;
  if ( x <= x3 && x > x4 ) p += 300;
  if ( x <= x4 )           p += 400;

  if ( y >= y1 && y < y2 ) p += 10;
  if ( y >= y2 && y < y3 ) p += 20;
  if ( y >= y3 && y < y4 ) p += 30;
  if ( y >= y4 )           p += 40;

  if ( z <= z1 && z > z2 ) p += 1;
  if ( z <= z2 && z > z3 ) p += 2;
  if ( z <= z3 && z > z4 ) p += 3;
  if ( z <= z4 )           p += 4;
  
  /* printf("--cell %d, p = %d\n",i,p); */
  return p;

}

void draw_force(GLint i)
{
  GLfloat force_color[3] = {1.0f,0.1f,0.1f};

  if ( _force == 1 && _mycells[i].living_status != 'D' )
  {
    GLfloat s = sqrt(normsq3v(_mycells[i].force));
    if ( s > 1e-5f )
    {
      /* Add a small arrow along cell force */
      glPushMatrix(); /* save V */
      GLfloat s = sqrt(normsq3v(_mycells[i].force));
      /* go to cell position */
      glTranslatef( _mycells[i].position[0] + _mycells[i].radius*_mycells[i].force[0]/s, 
                    _mycells[i].position[1] + _mycells[i].radius*_mycells[i].force[1]/s, 
                    _mycells[i].position[2] + _mycells[i].radius*_mycells[i].force[2]/s); 
      glColor3fv(force_color);
      GLfloat theta = 180.0/M_PI*acos( _mycells[i].force[0]/s );
      glRotatef(theta, 0.0f , - _mycells[i].force[2], _mycells[i].force[1]);
      glScalef(5.0*s,_mycells[i].radius,_mycells[i].radius);
      glCallList(_cell_dlist + GEOM_ARROW);
      glPopMatrix(); /* current matrix is restored to V */
    }
  }
}

void draw_velocity(GLint i)
{
  const GLfloat velocity_color[3] = {0.1f,0.5f,1.0f};

  if ( _velocity == 1 && _mycells[i].living_status != 'D' )
  {
    GLfloat s = sqrt(normsq3v(_mycells[i].velocity));
    if ( s > 1e-5f )
    {
      /* Add a small arrow along cell velocity */
      glPushMatrix(); /* save V */
      GLfloat s = sqrt(normsq3v(_mycells[i].velocity));
      /* go to cell position */
      glTranslatef( _mycells[i].position[0] + _mycells[i].radius*_mycells[i].velocity[0]/s, 
                    _mycells[i].position[1] + _mycells[i].radius*_mycells[i].velocity[1]/s, 
                    _mycells[i].position[2] + _mycells[i].radius*_mycells[i].velocity[2]/s); 
      glColor3fv(velocity_color);
      GLfloat theta = 180.0/M_PI*acos( _mycells[i].velocity[0]/s );
      glRotatef(theta, 0.0f , - _mycells[i].velocity[2], _mycells[i].velocity[1]);
      glScalef(5.0*s,_mycells[i].radius,_mycells[i].radius);
      glCallList(_cell_dlist + GEOM_ARROW);
      glPopMatrix(); /* current matrix is restored to V */
    }
  }
}

void hoover_cell(GLint id, int color)
{
  glStencilFunc(GL_NOTEQUAL, STENCIL_HOOVER, 0xFF); /* hoovered cell stencil */
  glPushMatrix(); /* save V */
  glDisable(GL_LIGHTING);
  glTranslatef(_mycells[id].position[0], _mycells[id].position[1], _mycells[id].position[2]); /* current matrix is VT  */

  /* transparent color */
  glColor4f(0.0,0.0,0.0,0.0);
  glutSolidSphere(_mycells[id].radius + 0.04, 2*_s_res, 2*_s_res);

  for ( int j = 0; j < 10; j++)
  {
    /* highlight color */
    glColor4f(_highlight_color[color][0],
              _highlight_color[color][1],
              _highlight_color[color][2],
              1 - j/10.0);
    glutSolidSphere(_mycells[id].radius + 0.2 + 0.02*j, 2*_s_res, 2*_s_res);
  }

  glPopMatrix(); /* current matrix is restored to V */
  glEnable(GL_LIGHTING);

  /* printf("--id: %d, cell loc: %d\n", _mycells[id].id, locate_cell(id)); */
}

size_t get_cell_index(GLint id)
{
  /* fetch the index of cell with ID = id */
  size_t i;
  for ( i = 0; i < (size_t)_popsize ; i++ )
  {
    if ( _mycells[i].id == id) 
    {
      return i;
    }
  }
  return i;
}

void auxdisplay()
{
  size_t index = -1;

  glClear (GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT | GL_STENCIL_BUFFER_BIT );

  glMatrixMode (GL_PROJECTION);
  glLoadIdentity ();
  gluOrtho2D(-0.1,1.1,-0.3,1.3);


  if ( _single_cell_tracking )
  {
    if ( (size_t)_single_cell != _ts.id ) /* reset the time series if
                                             cell ID has changed */
    {
      _ts.t = (GLfloat*)realloc(_ts.t,sizeof(GLfloat));
      _ts.y = (GLfloat*)realloc(_ts.y,sizeof(GLfloat));
      _ts.step = (GLuint*)realloc(_ts.step,sizeof(GLuint));
      _ts.size = 0;
      _ts.id = 0;
      _ts.maxsize = 1;
      _ts.ymax = 0.1;
      _ts.tmax = INITIAL_TS_RANGE;
      _ts.tmin = 0.0f;
    }
    
    /* fetch the index of cell with ID = _single_cell */
    for ( size_t i = 0; i < (size_t)_popsize ; i++ )
    {
      if ( _mycells[i].id == _single_cell) 
      {
        index = i;
        _ts.id = _single_cell;
        break;
      }
    }
    
    /* _fpos.i - 1 is the trajectory step */
    size_t current;
    if ( _ts.size >= _ts.maxsize ) /* increase storage for time series */
    {
      _ts.maxsize *= 2*_ts.maxsize; 
      _ts.t = (GLfloat*)realloc(_ts.t,_ts.maxsize*sizeof(GLfloat));
      _ts.y = (GLfloat*)realloc(_ts.y,_ts.maxsize*sizeof(GLfloat));
      _ts.step = (GLuint*)realloc(_ts.step,_ts.maxsize*sizeof(GLuint));
    }

    current = _ts.size;
    for ( size_t i = 0; i < _ts.size ; i++ )
    {
      if ( _ts.step[i] == _fpos.i - 1 ) /* step already exists */
      {  
        current = i;
        _ts.size--;
        break;
      }
      if ( _ts.step[i] > _fpos.i - 1 ) /* have to insert step */
      {
        for ( size_t j = _ts.size; j > i; j-- )
        {
          _ts.t[j] = _ts.t[j-1];
          _ts.y[j] = _ts.y[j-1];
          _ts.step[j] = _ts.step[j-1];
        }
        current = i;
        break;
      }
    }

    /* add or replace point to time series */
    _ts.t[current] = _time;
    _ts.step[current] = _fpos.i - 1;
    if ( _ts.id != 0 && index < (size_t)_popsize )
    {
      _ts.y[current] = _mycells[index].signal;
    }
    else /* tracked cell not found -- assign NAN to y */
    {
      _ts.y[current] = NAN;
    }
    if ( _ts.y[current] > _ts.ymax ) 
    { 
      _ts.ymax = _ts.y[current];
    }
    if ( _ts.t[current] > _ts.tmax )
    {
      _ts.tmax = _ts.t[current];
    }
    if ( _ts.t[current] < _ts.tmin )
    {
      _ts.tmin = _ts.t[current];
    }
    _ts.size++;

    /***** add custom legend *****/ 
    draw_ts(&_ts, current);

  }
  else if ( _ts.size > 0 ) 
  {
    _ts.maxsize = 1;
    _ts.t = (GLfloat*)realloc(_ts.t,_ts.maxsize*sizeof(GLfloat));
    _ts.y = (GLfloat*)realloc(_ts.y,_ts.maxsize*sizeof(GLfloat));
    _ts.size = 0;
  } 

  /***** add custom legend *****/
  glColor3fv(_font_color);
  /* Render "ID: %d, SIGNAL: " part first */
  char legend_prefix[64];
  snprintf(legend_prefix, 64, "ID: %d, SIGNAL: ", _single_cell);
  glRasterPos2f(-0.05, 1.2);
  for (size_t i = 0; i < strlen(legend_prefix); i++) 
    glutBitmapCharacter(GLUT_BITMAP_HELVETICA_12, legend_prefix[i]);
  /* Then render signal name with subscript formatting */
  renderSignalNameWithSubscript(_signal_name, GLUT_BITMAP_HELVETICA_12, GLUT_BITMAP_HELVETICA_10);

  glutSwapBuffers();
}

void statsdisplay()
{
  char legend_string[128];

  glClear (GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT | GL_STENCIL_BUFFER_BIT );

  glMatrixMode (GL_PROJECTION);
  glLoadIdentity ();
  gluOrtho2D(-0.1,1.1,-0.3,1.3);


  /* _fpos.i - 1 is the trajectory step */

  size_t current;
  if ( _stats_ts.size >= _stats_ts.maxsize ) /* increase storage for time series */
  {
    _stats_ts.maxsize *= 2*_stats_ts.maxsize; 
    _stats_ts.t = (GLfloat*)realloc(_stats_ts.t,_stats_ts.maxsize*sizeof(GLfloat));
    _stats_ts.y = (GLfloat*)realloc(_stats_ts.y,_stats_ts.maxsize*sizeof(GLfloat));
    _stats_ts.y2 = (GLfloat*)realloc(_stats_ts.y2,_stats_ts.maxsize*sizeof(GLfloat));
    _stats_ts.step = (GLuint*)realloc(_stats_ts.step,_stats_ts.maxsize*sizeof(GLuint));
  }

  /* shift the time series */
  current = _stats_ts.size;
  for ( size_t i = 0; i < _stats_ts.size ; i++ )
  {
    if ( _stats_ts.step[i] == _fpos.i - 1 ) /* step already exists */
    {  
      current = i;
      _stats_ts.size--;
      break;
    }
    if ( _stats_ts.step[i] > _fpos.i - 1 ) /* have to insert step */
    {
      for ( size_t j = _stats_ts.size; j > i; j-- )
      {
        _stats_ts.t[j] = _stats_ts.t[j-1];
        _stats_ts.y[j] = _stats_ts.y[j-1];
        _stats_ts.y2[j] = _stats_ts.y2[j-1];
        _stats_ts.step[j] = _stats_ts.step[j-1];
      }
      current = i;
      break;
    }
  }

  /* fetch the stats -- add or replace point to time series */
  _stats_ts.t[current]= _time;
  _stats_ts.step[current] = _fpos.i - 1;
  _stats_ts.y[current] = 0.0;
  _stats_ts.y2[current] = 0.0;

  /* compute signal average */
  for ( size_t i = 0; i < (size_t)_popsize ; i++ ) 
  {
    _stats_ts.y[current] += _mycells[i].signal;
  }
  _stats_ts.y[current] /= _popsize;
  /* end compute signal average */

  /* compute signal standard deviation */
  for ( size_t i = 0; i < (size_t)_popsize ; i++ ) 
  {
    _stats_ts.y2[current] += (_mycells[i].signal - _stats_ts.y[current])*(_mycells[i].signal - _stats_ts.y[current]);
  }
  _stats_ts.y2[current] /= (_popsize - 1);
  _stats_ts.y2[current] = sqrt(_stats_ts.y2[current]);
  /* end compute signal deviation */

  if ( _stats_ts.y[current] > _stats_ts.ymax ) 
  {
    _stats_ts.ymax = _stats_ts.y[current];
  }
  if ( _stats_ts.t[current] > _stats_ts.tmax )
  {
      _stats_ts.tmax = _stats_ts.t[current];
  }
  if ( _stats_ts.t[current] < _stats_ts.tmin )
  {
    _stats_ts.tmin = _stats_ts.t[current];
  }
  
  _stats_ts.size++;

  draw_ts(&_stats_ts, current);

  /***** add custom legend *****/
  glColor3fv(_font_color);
  /* Render "<signal_name> = value" with subscript formatting */
  glRasterPos2f(-0.05, 1.2);
  glutBitmapCharacter(GLUT_BITMAP_HELVETICA_12, '<');
  renderSignalNameWithSubscript(_signal_name, GLUT_BITMAP_HELVETICA_12, GLUT_BITMAP_HELVETICA_10);
  /* Get current raster position and continue with " = value" */
  GLfloat raster_pos[4];
  glGetFloatv(GL_CURRENT_RASTER_POSITION, raster_pos);
  glWindowPos2f(raster_pos[0], raster_pos[1]);
  snprintf(legend_string, 128, "> = %6.2g", _stats_ts.y[current]);
  for (size_t i = 0; i < strlen(legend_string); i++) 
    glutBitmapCharacter(GLUT_BITMAP_HELVETICA_12, legend_string[i]);

  glutSwapBuffers();
}

void create_menu()
{
  int load_colormap_menu;
  int display_signal_menu;
  char menu_entry[256];

  /* create colormap submenu */
  load_colormap_menu = glutCreateMenu(colormap_menu);
  glutAddMenuEntry("custom...", SUBMENU_COLORMAP_CUSTOM); 

  for ( int i = 0; i < NB_COLORMAPS; i++ )
  {
    glutAddMenuEntry(_list_colormaps.names[i], SUBMENU_COLORMAP_LIST_0 + i);
  }

  glutAddMenuEntry("old jet", SUBMENU_COLORMAP_JET);  
  glutAddMenuEntry("old spring", SUBMENU_COLORMAP_SPRING);
  glutAddMenuEntry("old neon", SUBMENU_COLORMAP_NEON);  

  /* create display signal submenu */
  display_signal_menu = glutCreateMenu(cell_signal_menu);
  for ( unsigned int i = 0; i < _list_signals.size; i++ )
  {
    snprintf(menu_entry, 255, "%s (%d)", _list_signals.name[i], i + 1);   
    glutAddMenuEntry(menu_entry, SUBMENU_DISPLAY_SIGNAL_0 + i);  
  }

  /* create main menu */
  _main_menu = glutCreateMenu(menu);
  glutAddSubMenu("Colormap", load_colormap_menu);  /* pos 1 */
  glutAddMenuEntry("Hide plane mode (m)", MENU_HIDE_PLANE); /* pos 2 */
  glutAddMenuEntry("Clip cells (c)", MENU_CLIP_CELLS); /* pos 3 */
  glutAddMenuEntry("Display next signal (])", MENU_NEXT_SIGNAL); /* pos 4 */
  glutAddMenuEntry("Display previous signal ([)", MENU_PREVIOUS_SIGNAL); /* pos 5 */
  glutAddSubMenu("Display signal", display_signal_menu); /* pos 6 */
  glutAddMenuEntry("Reset view (R)", MENU_RESET); /* pos 7 */
  glutAddMenuEntry("Help", MENU_HELP); /* pos 8 */
  glutAddMenuEntry("Quit (Q)", MENU_QUIT); /* last pos */
  /* end create main menu */

}


void menu(int value)
{

  switch(value) {
    case MENU_HIDE_PLANE: /* hide plane mode */
      _move_plane = 1 - _move_plane;
      glutPostRedisplay();
      break;
    case MENU_CLIP_CELLS: /* clip out-of-bound cells */
      if ( _clipcell == 0 )
      { 
        _clipcell = 1;  
        _min_clip_signal = _min_signal;
        _max_clip_signal = _max_signal;
      }
      else
      { 
        _clipcell = 0; 
      }
      standstill();
      glutTimerFunc(1.0/_speed,animate,0);
      glutPostRedisplay();
      break;
    case MENU_NEXT_SIGNAL:
      next_signal();
      read_min_max_signals();
      standstill();
      glutTimerFunc(1.0/_speed,animate,0);
      break;
    case MENU_PREVIOUS_SIGNAL:
      previous_signal();
      read_min_max_signals();
      standstill();
      glutTimerFunc(1.0/_speed,animate,0);
      break;
    case MENU_RESET:
      reset_view();
      glutPostRedisplay();
      break;
    case MENU_HELP:
      print_help(NULL);
      break;
    case MENU_QUIT:
      exit(0);
  }

  if ( _move_plane )
    glutChangeToMenuEntry(2,"Rotation mode (m)", MENU_HIDE_PLANE);
  else
    glutChangeToMenuEntry(2,"Hide plane mode (m)", MENU_HIDE_PLANE);
  
  if ( _clipcell )
    glutChangeToMenuEntry(3,"Show all cells (c)", MENU_CLIP_CELLS);
  else
    glutChangeToMenuEntry(3,"Clip cells (c)", MENU_CLIP_CELLS); /* pos 3 */

}

void colormap_menu(int value)
{
  
  switch(value) {
    case SUBMENU_COLORMAP_JET: /* load colormap */
      free(_colormap);
      init_color_map_jet();
      standstill();
      glutTimerFunc(1.0/_speed,animate,0);
      break;
    case SUBMENU_COLORMAP_SPRING: 
      free(_colormap);
      init_color_map_spring();
      standstill();
      glutTimerFunc(1.0/_speed,animate,0);
      break;
    case SUBMENU_COLORMAP_NEON: 
      free(_colormap);
      init_color_map_neon();
      standstill();
      glutTimerFunc(1.0/_speed,animate,0);
      break;
    case SUBMENU_COLORMAP_CUSTOM:
      free(_colormap);
      memset(_inputstr,0,128); 
      _keyboard_input_mode = KEY_COLORMAP;
      printf("enter colormap file name: ");
      fflush(stdout);
      standstill();
      glutTimerFunc(1.0/_speed,animate,0);
      break;
  }

  if ( value >= SUBMENU_COLORMAP_LIST_0 &&
       value <  SUBMENU_COLORMAP_LIST_0 + NB_COLORMAPS )
  {
    free(_colormap);
    unsigned long *c = (unsigned long *)malloc(5*sizeof(unsigned long)); 
    _colormap_name = value - SUBMENU_COLORMAP_LIST_0;
    printf("colormap %s\n", _list_colormaps.names[_colormap_name]);
    for ( int i = 0; i < 5; i++ )
    { 
      c[i] = _list_colormaps.colors[_colormap_name][i];
      /* printf("%lX ", c[i]); */
    }
    /* printf("\n"); */
    _colormap_len = NB_COLORS;
    make_colormap_array(c, 5, &_colormap, _colormap_len ); 
    /* printf("-- colormap: %f %f %f\n,",_colormap[0].c[0],_colormap[0].c[1],_colormap[0].c[2]); */
    standstill();
    glutTimerFunc(1.0/_speed,animate,0);
    free(c);
  }
}

void cell_signal_menu(int value)
{

  for ( unsigned int i = 0; i < _list_signals.size; i++ )
  {
    if ( (unsigned int)value == SUBMENU_DISPLAY_SIGNAL_0 + i )
    {
      set_signal(i);
      read_min_max_signals();
      standstill();
      glutTimerFunc(1.0/_speed,animate,0);
      break;
    }
  }

}

void draw_ts(timeseries *ts, size_t current) 
{
  GLfloat tt, tt1, s;
  char legend_string[128];
  GLint viewport[4];
  const int y_ticks = 6;
  const int x_ticks = 11;
  GLfloat ylims[2] = {_min_signal, _max_signal};
  GLfloat std_whisker_width = (float)glutGet(GLUT_WINDOW_WIDTH)/current;
  

  if ( (_max_signal - _min_signal) < 0.1*fabs(_min_signal) ) {
    ylims[0] -= 0.05;
    ylims[1] += 0.05;
  }

  clamp(&std_whisker_width,1.0,3.0);
  //printf("std_whisker_width: %f, current: %zu\n",std_whisker_width,current);

  /***** draw the standard deviation *****/
  if ( ts->y2 != NULL ) 
  {
    glEnable(GL_LINE_SMOOTH);
    glLineWidth(std_whisker_width);
    glBegin(GL_LINES);
      for ( size_t i = 0; i < current + 1 ; i++ )
      {
        if ( (tt =(ts->t[i] - ts->t[current])/(ts->tmax - ts->tmin) + 1.0) > 0.0 &&
              tt <= 1.0 )
        {
          tt1 = (ts->t[i-1] - ts->t[current])/(ts->tmax - ts->tmin) + 1.0;
          s = (GLfloat)tt;
          //s *= sqrt(s);
          glColor3f( (1.0 - 1.*s)*_background_color[0] + 0.5*s,
                     (1.0 - 1.*s)*_background_color[1] + 0.5*s,
                     (1.0 - 1.4*s)*_background_color[2] + 1.1*s);
          if ( tt - tt1 > 0.1 ) /* skip line if delta t is > 0.1 or 10% of axes */
          {
            glColor3fv(_background_color);
          }
          glVertex3d(tt, (ts->y[i] - ts->y2[i] - ylims[0])/(ylims[1] - ylims[0]), 0.0);
          glVertex3d(tt, (ts->y[i] + ts->y2[i] - ylims[0])/(ylims[1] - ylims[0]), 0.0);
        }
      }
    glEnd();
    glDisable(GL_LINE_SMOOTH);
  }

  /***** draw the time series *****/
  glEnable(GL_LINE_SMOOTH);
  glLineWidth(0.5);
  glBegin(GL_LINE_STRIP);
    for ( size_t i = 0; i < current + 1 ; i++ )
    {
      if ( (tt =(ts->t[i] - ts->t[current])/(ts->tmax - ts->tmin) + 1.0) > 0.0 &&
            tt <= 1.0 )
      {
        tt1 = (ts->t[i-1] - ts->t[current])/(ts->tmax - ts->tmin) + 1.0;
        s = (GLfloat)tt;
        //s *= sqrt(s);
        glColor3f( (1.0 - 1.6*s)*_background_color[0] + 0.9*s,
                   (1.0 - 1.6*s)*_background_color[1] + 0.9*s,
                   (1.0 - 2.0*s)*_background_color[2] + 1.5*s);
        if ( tt - tt1 > 0.1 ) /* skip line if delta t is > 0.1 or 10% of axes */
        {
          glColor3fv(_background_color);
        }
        glVertex3d(tt1, (ts->y[i-1] - ylims[0])/(ylims[1] - ylims[0]), 0.0);
        glVertex3d(tt, (ts->y[i] - ylims[0])/(ylims[1] - ylims[0]), 0.0);
      }
    }
  glEnd();
  glDisable(GL_LINE_SMOOTH);

  /***** draw the axes *****/
  glGetIntegerv(GL_VIEWPORT, viewport);
  glColor3fv(_border_color);
  glLineWidth( 1.0);
  glBegin(GL_LINE_STRIP);         /* xy-axes */ 
    glVertex3d(1.05, 0.0, 0.0);
    glVertex3d(0.0, 0.0, 0.0);
    glVertex3d(0.0, 1.05, 0.0);
  glEnd();
  glBegin(GL_LINES);              /* ticks */
    for ( size_t i = 0; i < x_ticks ; i++ )
    {
      glVertex3d((float)i/(x_ticks - 1), -5.0/viewport[3], 0.0);
      glVertex3d((float)i/(x_ticks - 1), 0.0, 0.0);                /* at (i/10.0) */
    }
    for ( size_t i = 0; i < y_ticks ; i++ )
    {
      glVertex3d(-5.0/viewport[2], (float)i/(y_ticks - 1), 0.0);
      glVertex3d(0.0, (float)i/(y_ticks - 1), 0.0);                /* at (0,i/5.0) */
    }
  glEnd();

  /***** print legend *****/
  glColor3fv(_font_color);

  for ( size_t i = 0; i < y_ticks ; i++ )
  {
    snprintf(legend_string, 128, "%6.2g",ylims[0] + (ylims[1] - ylims[0])*i/(y_ticks - 1.0));
    glRasterPos2f(-(50.0/viewport[2]), (float)i/(y_ticks - 1)); /* also sets current raster color */
    for (size_t i = 0; i < strlen(legend_string); i++) glutBitmapCharacter(GLUT_BITMAP_HELVETICA_12, legend_string[i]);
  }

  for ( size_t i = 0; i < x_ticks ; i++ )
  {
    snprintf(legend_string, 128, "%2.1f",_time - (ts->tmax - ts->tmin)*(i/(x_ticks - 1.0)));
    glRasterPos2f((float)(x_ticks - 1 - i)/(x_ticks - 1), -30.0/viewport[3]); /* also sets current raster color */
    for (size_t i = 0; i < strlen(legend_string); i++) glutBitmapCharacter(GLUT_BITMAP_HELVETICA_12, legend_string[i]);
  }
  
}

void reshape (int w, int h)
{
  glViewport (0, 0, (GLsizei) w, (GLsizei) h);
  glMatrixMode (GL_PROJECTION);
  glLoadIdentity ();
  /*gluPerspective(30.0, (1.0*w)/h, 0.5*_worldsize, 5*_worldsize); */
  /* gluPerspective( 45.0f, (1.0*w)/h, _near_clip, _far_clip ); */
  gluPerspective( _fovy, (1.0*w)/h, _near_clip, _far_clip );
  glMatrixMode (GL_MODELVIEW);
}


void reset_view()
{
  if ( _viewpoint_str != NULL )
  {
    char *endptr = _viewpoint_str;
    _azim_angle  = strtof(endptr,&endptr);
    _inc_angle   = strtof(endptr,&endptr);
    _twist_angle = strtof(endptr,&endptr);
    _distance    = strtof(endptr,&endptr);
    _pan_x       = strtof(endptr,&endptr);
    _pan_y       = strtof(endptr,&endptr);
  }
  else /* default view */
  {
    _distance = _near_clip + (_far_clip - _near_clip) / 1.5;
    _twist_angle = 0.0;	/* rotation of viewing volume (camera) */
    _inc_angle = 85.0;
    _azim_angle = -30.0;
    _pan_x = _worldsize.x/10.0; 
    _pan_y = 0.0;
  }
  if ( strlen(_hide_str) ) 
  {
    char *endptr = _hide_str;
    _x_plane = strtof(endptr,&endptr);
    _y_plane = strtof(endptr,&endptr);
    _z_plane = strtof(endptr,&endptr);
  }
}


void polar_view( GLfloat distance, GLfloat azimuth, GLfloat incidence, GLfloat twist)
{
  /* printf(" incidence %f azimuth %f\n", incidence, azimuth); */
  // glTranslatef( -_worldsize.x/2.0, -_worldsize.y/2.0, -distance);
  glTranslatef( _pan_x, _pan_y, 0.0f);
  glTranslatef( 0.0f, 0.0f, -distance + _worldsize.z/2.0 );
  glRotatef( -twist, 0.0f, 0.0f, 1.0f);
  glRotatef( -incidence, 1.0f, 0.0f, 0.0f);
  glRotatef( azimuth, 0.0f, 0.0f, 1.0f);
  glTranslatef( -_worldsize.x/2.0, -_worldsize.y/2.0, -_worldsize.z/2.0);
  glLightfv(GL_LIGHT1, GL_POSITION, _light1_pos);
}


void keyboard(unsigned char key, int, int)
{
    // Ensure we have proper window focus (Mac trackpad fix)
  // glutSetWindow(_main_window);
  if ( _keyboard_input_mode ) /* capture the key */
  {
    get_keyboard_input(key);
    return;
  }

  if ( key > '0' && key <= '9' ) 
  {
    set_signal(key - 49); /* '1' == 49 */
    read_min_max_signals();
    standstill();
    glutTimerFunc(1.0/_speed,animate,0);
  }

  switch (key)
  {
    case 'r': /* r or R : reset viewpoint */
      reset_view();
      glutPostRedisplay();
      break;
    case 'R':
      reset_view();
      glutPostRedisplay();
      break;
    case 32: /* space bar = pause the animation */
      if (_pause == 0)
      {  _pause = 1; }
      else
      {
        _pause = 0;
        glutTimerFunc(1.0/_speed,animate,0);
      }
      break;
    case 'j': /* decrease max signal */
      _max_signal -= 0.1*(_max_signal - _min_signal); 
      standstill();
      glutTimerFunc(1.0/_speed,animate,0);
      break;
    case 'k': /* increase max signal */
      _max_signal += 0.1*(_max_signal - _min_signal);
      standstill();
      glutTimerFunc(1.0/_speed,animate,0);
      break;
    case 'h': /* decrease min signal */
      _min_signal -= 0.1*(_max_signal - _min_signal); 
      standstill();
      glutTimerFunc(1.0/_speed,animate,0);
      break;
    case 'l': /* increase min signal */
      _min_signal += 0.1*(_max_signal - _min_signal); 
      standstill();
      glutTimerFunc(1.0/_speed,animate,0);
      break;
    case ']': /* display next signal */
      next_signal();
      read_min_max_signals();
      standstill();
      glutTimerFunc(1.0/_speed,animate,0);
      break;
    case '[':
      previous_signal();
      read_min_max_signals();
      standstill();
      glutTimerFunc(1.0/_speed,animate,0);
      break;
    case 'n': /* stop at next time step */
      _pause = 1;
      glutTimerFunc(0,animate,0);
      break;
    case 'b': /* back to previous step */
      _pause = 1;
      backtrack(1);
      glutTimerFunc(1.0/_speed,animate,0);
      break;
    case 'B': /* back to previous step */
      _pause = 1;
      backtrack(50);
      glutTimerFunc(1.0/_speed,animate,0);
      break;
    case 'N':
      _skip_frame = 50;
      printf("skipping 50 steps\n");
      glutTimerFunc(1.0/_speed,animate,0);
      break;
    case 'c': /* toggle clip: hide cells that have signal outside min/max signal */
      if ( _clipcell == 0 )
      { 
        _clipcell = 1;  
        _min_clip_signal = _min_signal;
        _max_clip_signal = _max_signal;
      }
      else
      { 
        _clipcell = 0; 
      }
      standstill();
      glutTimerFunc(1.0/_speed,animate,0);
      break;
    case 'w': /* toggle show/hide world border display */
      _draw_guides++;
      _draw_guides %= 8;
      standstill();
      glutTimerFunc(1.0/_speed,animate,0);
      break;
    case 't':
      if ( _single_cell_tracking == 0)
      {
        printf("enter cell ID>> ");
        fflush(stdout);
        memset(_inputstr,0,128); 
        _keyboard_input_mode = KEY_TRACKING;
      }
      else
      {
        printf("  tracking off\n");
        _single_cell_tracking = 0;
      }
      standstill();
      glutTimerFunc(1.0/_speed,animate,0);
      break;
    case 's':
      save_image();
      print_snapshot();
      break;
    case 'S': // Save PNG with timestamp
      {
      time_t rawtime;
      struct tm * timeinfo;
      char timestamp[32];
      char timestamped_filename[256];
      
      time(&rawtime);
      timeinfo = localtime(&rawtime);
      strftime(timestamp, sizeof(timestamp), "%Y%m%d_%H%M%S", timeinfo);
      snprintf(timestamped_filename, 256, "simuscale_%s_t%.2f.png", timestamp, _time);
      
      save_image(timestamped_filename, "png");
      printf("Saved PNG: %s\n", timestamped_filename);
      }
      break;
    case 'P':
      {
          // Create a simple input dialog using system dialogs
          #ifdef __APPLE__
          char command[512];
          snprintf(command, sizeof(command), 
                  "osascript -e 'display dialog \"Enter filename:\" default answer \"image.png\"' -e 'text returned of result'");
          
          FILE *fp = popen(command, "r");
          if (fp) {
              if (fgets(_filename, sizeof(_filename), fp) != NULL) {
                  // Remove newline
                  _filename[strcspn(_filename, "\n")] = 0;
                  
                  // Add .png extension if not present
                  const char* ext = strrchr(_filename, '.');
                  if (!ext || strcmp(ext, ".png") != 0) {
                      strcat(_filename, ".png");
                  }
                  save_image(_filename, "png");
                  printf("Saved PNG: %s\n", _filename);
              }
              pclose(fp);
          }
          #elif defined(__linux__)
          // Use zenity for Linux
          char command[512];
          snprintf(command, sizeof(command), 
                  "zenity --entry --title=\"Save Image\" --text=\"Enter filename:\" --entry-text=\"image.png\"");
          
          FILE *fp = popen(command, "r");
          if (fp) {
              if (fgets(_filename, sizeof(_filename), fp) != NULL) {
                  _filename[strcspn(_filename, "\n")] = 0;
                  const char* ext = strrchr(_filename, '.');
                  if (!ext || strcmp(ext, ".png") != 0) {
                      strcat(_filename, ".png");
                  }
                  save_image(_filename, "png");
                  printf("Saved PNG: %s\n", _filename);
              }
              pclose(fp);
          }
          #endif
      }
      break;
    case 'o': 
        _show_neighborhood_sphere = 1 - _show_neighborhood_sphere;
        printf("Neighborhood sphere: %s\n", _show_neighborhood_sphere ? "ON" : "OFF");
        glutPostRedisplay();
        break; 
    case 'F': /* full screen */
      glutFullScreen();
      break;
    case 'f': /* leave full screen */
      glutReshapeWindow(500,500);
      glutPositionWindow(100,100);
      break;
    case 'g': /* toggle tagged cell highlighting */
      _tagged_cell_highlighting = 1 - _tagged_cell_highlighting;
      glutPostRedisplay();
      break;
    case 'm':
      _move_plane = 1 - _move_plane;
      glutPostRedisplay();
      break;
    case 'v':
      _draw_elements++;
      _draw_elements %= 8;
      printf("--draw elements: %d\n",_draw_elements);
      standstill();
      glutPostRedisplay();
      break;
    case 'x': /* toggle x-plane positioning */
      _pos_plane = 0;
      glutPostRedisplay();
      break;
    case 'y': /* toggle y-plane positioning */
      _pos_plane = 1;
      glutPostRedisplay();
      break;
    case 'z': /* toggle z-plane positioning */
      _pos_plane = 2;
      glutPostRedisplay();
      break;
    case 27:  /* escape to quit */
    case 'q': /* q to quit */
    case 'Q': /* Q to quit */
      exit(0);
  }
}


void special_keyboard(int key, int, int)
{
  int modifiers = glutGetModifiers();
  switch (key)
  {
    case GLUT_KEY_UP:
      /* printf("Up key pressed\n"); */
      if ( modifiers & GLUT_ACTIVE_ALT )
      {
        _distance -= 0.025*_worldsize.x; /* zoom in */
      }
      else if ( modifiers & GLUT_ACTIVE_SHIFT )
      {
        _inc_angle -= 0.05*_worldsize.z; /* decrease inclination */
      }
      else
      {
        _pan_y -= 0.025*_worldsize.z;
      }
      glutPostRedisplay();
      break;
    case GLUT_KEY_DOWN:
      /* printf("Down key pressed\n"); */
      if ( modifiers & GLUT_ACTIVE_ALT )
      {
        _distance += 0.025*_worldsize.x; /* zoom out */
      }
      else if ( modifiers & GLUT_ACTIVE_SHIFT )
      {
        _inc_angle += 0.05*_worldsize.z; /* increase declination */
      }
      else
      {
        _pan_y += 0.025*_worldsize.z;
      }
      glutPostRedisplay();
      break;
    case GLUT_KEY_LEFT:
      /* printf("Left key pressed\n"); */
      if ( modifiers & GLUT_ACTIVE_SHIFT )
      {
        /* Adjust the azimuth eye position  */
        _azim_angle += 0.025*_worldsize.x;
      }
      else
      {
        _pan_x += 0.025*_worldsize.x; /* pan right */
      }
      glutPostRedisplay();
      break;
    case GLUT_KEY_RIGHT:
      /* printf("Right key pressed\n"); */
      if ( modifiers & GLUT_ACTIVE_SHIFT )
      {
        /* Adjust the azimuth eye position  */
        _azim_angle -= 0.025*_worldsize.x;
      }
      else
      {
        _pan_x -= 0.025*_worldsize.x;
      }
      glutPostRedisplay();
      break;
  }
}


void get_keyboard_input(unsigned char key)
{
  if ( _keyboard_input_mode == KEY_TRACKING && key == 13 ) /* carriage return */
  {
    _keyboard_input_mode = KEY_OFF;
    _single_cell_tracking = sscanf(_inputstr,"%d",&_single_cell); 
    if ( _single_cell_tracking ) printf("\n  tracking cell %d\n",_single_cell);
    else printf("\n  unknown ID\n");
    memset(_inputstr, 0, 128);
    standstill();
    glutTimerFunc(0,animate,0);
    return;
  }
  if ( _keyboard_input_mode == KEY_COLORMAP && key == 13 ) /* carriage return */
  {
    _keyboard_input_mode = KEY_OFF;
    /* input string -> _inputstr */
    printf("\n  loading colormap: %s\n", _inputstr);
    init_color_map_custom(_inputstr);
    memset(_inputstr, 0, 128);
    return;
  }
  if ( key == 8 || key == 16 || key == 127) /* backspace */ 
  {
    if ( strlen(_inputstr) > 0 )
    {
      _inputstr[strlen(_inputstr)-1] = 0;  
      PR_BS_CLEAR; fflush(stdout);
    }
    return;
  }

  /* get key */
  _inputstr[strlen(_inputstr)] = key;  
  putchar(key); fflush(stdout);
  return;
}

void update_signal(int from_what)
{
  if ( from_what == SIGNAL_INDEX )
  {
    _signal_name = (char*)realloc(_signal_name,(strlen(_list_signals.name[_signal_index - _first_signal_col]) + 1)*sizeof(char));
    snprintf(_signal_name, strlen(_list_signals.name[_signal_index - _first_signal_col]) + 1, "%s", _list_signals.name[_signal_index - _first_signal_col]); 
  }
  else if ( from_what == SIGNAL_NAME ) 
  {
    int index = _first_signal_col;
    while ( strncmp(_signal_name,_list_signals.name[index - _first_signal_col],127) != 0 )
    {
      index++;
    }
    if ( ( index - _first_signal_col ) < (int)_list_signals.size )
    {
      _signal_index = index;
    }
    else /* do not update _signal_index, revert _signal_index back to match _signal_index */
    {
      update_signal(SIGNAL_INDEX);
    }
  }
}

void next_signal()
{
      _signal_index++;
      _signal_index = ( ( _signal_index - _first_signal_col ) % _list_signals.size ) + _first_signal_col;
      update_signal(SIGNAL_INDEX);
}

void previous_signal()
{
      _signal_index--;
      _signal_index = ( ( _signal_index - _first_signal_col + _list_signals.size ) % _list_signals.size );
      _signal_index +=  _first_signal_col;
      update_signal(SIGNAL_INDEX);
}

void set_signal(int sig)
{
      sig %= _list_signals.size;
      _signal_index = sig + _first_signal_col;
      update_signal(SIGNAL_INDEX);
}

GLvoid mouse( GLint button, GLint state, GLint x, GLint y )
{
    // Force focus back to the current window (Mac trackpad fix)
  // glutSetWindow(glutGetWindow());
  if (state == GLUT_DOWN)
  {
    switch (button)
    {
      case GLUT_LEFT_BUTTON:
        if ( _move_plane )
        {
          if ( glutGetModifiers() & GLUT_ACTIVE_SHIFT )
          {
            _action = MOVE_PLANE_FINE;
          }
          else
          {
            _action = MOVE_PLANE;
          }
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
    _x_start = x;
    _y_start = y;

    if ( ( button == GLUT_LEFT_BUTTON ) &&  ( glutGetModifiers() & GLUT_ACTIVE_SHIFT ) )
    {
      _single_cell_tracking = select_cell(x,y);
    }

  }
  else
  {
    _action = MOVE_NONE;
  }

}

int select_cell(GLint x, GLint y)
{
  int window_height = glutGet(GLUT_WINDOW_HEIGHT);
  GLfloat depth;
  GLint stencil;
  GLdouble model[16];
  GLdouble projection[16];
  GLint view[4];
  GLdouble x_coord,y_coord,z_coord;
  GLdouble min_dist = INFINITY;
  GLdouble dist;

  glReadPixels(x, window_height - y - 1, 1, 1, GL_STENCIL_INDEX, GL_UNSIGNED_INT, &stencil);
  glReadPixels(x, window_height - y - 1, 1, 1, GL_DEPTH_COMPONENT, GL_FLOAT, &depth);

  /* printf("--stencil: %u\n\n",_single_cell); */

  glGetDoublev(GL_MODELVIEW_MATRIX,model);
  glGetDoublev(GL_PROJECTION_MATRIX,projection);
  glGetIntegerv(GL_VIEWPORT,view);

  /* get approximate nearest object coordinates */
  gluUnProject((GLdouble)x, (GLdouble)(window_height - y - 1),
      depth,model,projection,view,&x_coord,&y_coord,&z_coord);

  printf("object coordinates: (%f,%f,%f), ", x_coord, y_coord, z_coord);

  /* user did not click on untracked cell */
  if ( stencil != STENCIL_CELL && 
       stencil != STENCIL_HOOVER )  
  {
    printf(" no selection\n\n");
    glutSetWindow(_main_window);
    glutPostRedisplay();
    glutSetWindow(_aux_window);
    glutPostRedisplay();
    glutSetWindow(_stats_window);
    glutPostRedisplay();

    return 0;
  }

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
      _single_cell = _mycells[i].id;
    }
  }

  printf("ID: %u %c\n\n",_single_cell,_mycells[_single_cell].living_status);
  
  glutPostRedisplay();
  return 1;
}

void hoover_main(GLint x, GLint y)
{
  int window_height = glutGet(GLUT_WINDOW_HEIGHT);
  int window_width = glutGet(GLUT_WINDOW_WIDTH);
  GLfloat depth;
  GLint stencil;
  GLdouble model[16];
  GLdouble projection[16];
  GLint view[4];
  GLdouble x_coord,y_coord,z_coord;
  GLdouble min_dist = INFINITY;
  GLdouble dist;

  glutSetMenu(_main_menu);
  glutDetachMenu(GLUT_LEFT_BUTTON);

  /* if mouse is not in window, return immediately */
  if ( x < 0 || x > window_width || 
       y < 0 || y > window_height )
  {
    return ;
  }

  /* test - attach menu if mouse is in the upper left corner */
  if ( x < 56 && y < 32 )
  {
    glutSetMenu(_main_menu);
    glutAttachMenu(GLUT_LEFT_BUTTON);
    _menu_hoovering = 2;
    glutPostRedisplay();
    return ;
  }
  else if ( x < 84 && y < 48 )
  {
    _menu_hoovering = 1;
    glutPostRedisplay();
    return ;
  }
  else
  {
    _menu_hoovering = 0;
  }

  if ( _single_cell_hoovering == HOOV_DISABLED )
  {
    return; /* do not try to detect cells */
  }

  /**** execute only if hoovering is not disabled ****/
  _single_cell_hoovering = HOOV_OFF;

  /* get depth and stencil data */
  glReadPixels(x, window_height - y - 1, 1, 1, GL_STENCIL_INDEX, GL_UNSIGNED_INT, &stencil);
  glReadPixels(x, window_height - y - 1, 1, 1, GL_DEPTH_COMPONENT, GL_FLOAT, &depth);

  glGetDoublev(GL_MODELVIEW_MATRIX,model);
  glGetDoublev(GL_PROJECTION_MATRIX,projection);
  glGetIntegerv(GL_VIEWPORT,view);

  /* get approximate nearest object coordinates */
  gluUnProject((GLdouble)x, (GLdouble)(window_height - y - 1),
      depth,model,projection,view,&x_coord,&y_coord,&z_coord);


  /* user hoovered outside untracked cell */
  if ( stencil != STENCIL_CELL &&
       stencil != STENCIL_HOOVER )  
  {
    glutSetWindow(_main_window);
    glutPostRedisplay();
    glutSetWindow(_aux_window);
    glutPostRedisplay();
    glutSetWindow(_stats_window);
    glutPostRedisplay();

    return;
  }

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
      _single_cell_hoovered = _mycells[i].id;
    }
  }

  _single_cell_hoovering = HOOV_ON;
  glutPostRedisplay();
  return;
}

GLvoid motion( GLint x, GLint y )
{
  switch (_action)
  {
    case MOVE_EYE:
      /* Adjust the eye position based on the mouse position */
      _azim_angle += (GLdouble) (x - _x_start);
      _inc_angle -= (GLdouble) (y - _y_start);
      break;
    case TWIST_EYE:
      /* Adjust the eye twist based on the mouse position */
      _twist_angle = fmod(_twist_angle+(x - _x_start), 360.0);
      break;
    case ZOOM:
      /* Adjust the eye distance based on the mouse position */
      _distance -= (GLdouble) (y - _y_start)/10.0;
      break;
    case MOVE_PLANE:
    case MOVE_PLANE_FINE:
      move_plane(x,y);
      break;
  }

  fprintf(_log_file,"viewpoint: %g %g %g %g %g %g\n", _azim_angle, _inc_angle, _twist_angle, _distance, _pan_x, _pan_y);

  /* Update the stored mouse position for later use */
  _x_start = x;
  _y_start = y;

  glutPostRedisplay();
}

GLvoid passive_motion( GLint x, GLint y)
{
  hoover_main(x,y);
}

void backtrack( unsigned int nbrframes )
{
    int newframe = _fpos.i - nbrframes - 1;
    if ( newframe < 0 ) 
      _fpos.i = 0;
    else
      _fpos.i = newframe;
    /* printf("bactrack to %ld\n",_fpos.pos[_fpos.i]); */
    fseek(_traj_file, _fpos.pos[_fpos.i], SEEK_SET);      
}

void darken( float *color, float amount )
{
  color[0] *=  amount;
  color[1] *=  amount;
  color[2] *=  amount;
}

void lighten( float *color, float amount )
{
  color[0] += (1.0f - color[0])*amount;
  color[1] += (1.0f - color[0])*amount;
  color[2] += (1.0f - color[0])*amount;
}

GLfloat distsq3v( GLfloat *p1, GLfloat *p2)
{
  return       (p1[0] - p2[0])*(p1[0] - p2[0]) + \
               (p1[1] - p2[1])*(p1[1] - p2[1]) + \
               (p1[2] - p2[2])*(p1[2] - p2[2]);
}

GLfloat normsq3v( GLfloat *v)
{
  return v[0]*v[0] + v[1]*v[1] + v[2]*v[2];
}


void assign3v(const GLfloat *source, size_t size_source, GLfloat *target, size_t size_target)
{
  if ( size_target % size_source )
  {
    fprintf(stderr, "Error in assign3, size of source (%zu) does no divide size of target (%zu): %zu.\n", size_source, size_target, size_target % size_source); 
    exit( EXIT_FAILURE );
  }

  for (size_t i = 0; i < size_target/size_source; i++ )
  {
    for (size_t j = 0; j < 3*size_source; j++ )
    {
      target[3*size_source*i + j] = source[j];
    }
  }

}

void copy_triangle(const triangle source, triangle *target)
{
  target->v1[0] = source.v1[0];
  target->v1[1] = source.v1[1];
  target->v1[2] = source.v1[2];
  target->v2[0] = source.v2[0];
  target->v2[1] = source.v2[1];
  target->v2[2] = source.v2[2];
  target->v3[0] = source.v3[0];
  target->v3[1] = source.v3[1];
  target->v3[2] = source.v3[2];
  target->normal[0] = source.normal[0];
  target->normal[1] = source.normal[1];
  target->normal[2] = source.normal[2];
}

/* rotate3v rotates _in_ in 3D by angle _angle_ 
 * around axis _axis_
 * output is written in _out_
 * 
 * 3D rotation matrix around axis (x,y,z)
 * R = [ c+x2(1-c)  xy(1-c)-zs  xz(1-c)+ys ]
 *     [ xy(1-c)+zs c+y2(1-c)   yz(1-c)-xs ]
 *     [ xz(1-c)-ys yz(1-c)+xs  c+z2(1-c)  ]
 *
 */
void rotate3v(const GLfloat *axis, const GLfloat *angle, GLfloat *points, const size_t nbr_points)
{
  GLfloat x = axis[0];
  GLfloat y = axis[1];
  GLfloat z = axis[2];
  GLfloat px, py, pz, c, s, c1;

  /* printf("c = %f, s = %f, c1 = %f, angle = %f\n",c,s,c1,angle); */

  for ( size_t i = 0; i < nbr_points ; i++ )
  {
    px = points[3*i + 0];
    py = points[3*i + 1];
    pz = points[3*i + 2];
    c = cos(angle[i]);
    s = sin(angle[i]);
    c1 = 1 - c;

    points[3*i + 0] = (c + x*x*c1)*px + \
               (x*y*c1 - z*s)*py + \
               (x*z*c1 + y*s)*pz;
    points[3*i + 1] = (x*y*c1 + z*s)*px + \
               (c + y*y*c1)*py + \
               (y*z*c1 + x*s)*pz;
    points[3*i + 2] = (x*z*c1 - y*s)*px + \
               (y*z*c1 + x*s)*py + \
               (c + z*z*c1)*pz;  
  }
  /* printf("out: %f %f %f\n",out[0],out[1],out[2]); */

}

void translate3v(const GLfloat shift[3], GLfloat *points, size_t nbr_points)
{
  for ( size_t i = 0; i < nbr_points ; i++ )
  {
    points[3*i] += shift[0]; 
    points[3*i + 1] += shift[1]; 
    points[3*i + 2] += shift[2]; 
  }
}

/* normal to the plane defined by p1, p2, p3 */
void normal3v( GLfloat *p1, GLfloat *p2, GLfloat *p3, GLfloat *nv )
{
  GLfloat a1 = p2[0] - p1[0];
  GLfloat a2 = p2[1] - p1[1];
  GLfloat a3 = p2[2] - p1[2];
  GLfloat b1 = p3[0] - p1[0];
  GLfloat b2 = p3[1] - p1[1];
  GLfloat b3 = p3[2] - p1[2];
  GLfloat dyz  = a2*b3 - b2*a3;
  GLfloat dxz  = a1*b3 - b1*a3;
  GLfloat dxy  = a1*b2 - b1*a2;

  nv[0] =   dyz;
  nv[1] = - dxz;
  nv[2] =   dxy;

  GLfloat norm = sqrt(normsq3v(nv));

  nv[0] /= norm;
  nv[1] /= norm;
  nv[2] /= norm;

}

void object_material()
{
  GLfloat mat_ambient_refl[4];
  GLfloat mat_shininess[1];
  GLfloat mat_diffuse_color[4];
  GLfloat mat_emission[4];
  GLfloat *p;

  glShadeModel (_shade_model);

  /* color of the light */
  glLightfv(GL_LIGHT0, GL_AMBIENT, _white_light);
  glLightfv(GL_LIGHT0, GL_DIFFUSE, _diffuse_light);
  /* glLightfv(GL_LIGHT0, GL_SPECULAR, _white_light); */

  /* position of the light */
  GLfloat light_position[] = {0.0, 0.0, 0.0, 1.0}; /* will move with the viewpoint */
  /* w=0 thus parallel rays like the sun, called directional light source */
  /* The direction (along z axis here) is transformed by the current modelview matrix */
  /* There is no attenuation for a directional light source */
  glLightfv(GL_LIGHT0, GL_POSITION, light_position);

  /* lighting model */
  GLfloat global_ambient_light[] = {0.0, 0.0, 0.0, 1.0}; /* to see objects even if no light source */
  glLightModelfv(GL_LIGHT_MODEL_AMBIENT, global_ambient_light);
  glLightModeli(GL_LIGHT_MODEL_LOCAL_VIEWER, GL_FALSE); /* infinite viewpoint */
  glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, GL_FALSE); /* back faces are inside the spheres, never seen */

  /* material for the objects */
  switch(_color_model)
  {
   case COLOR_MODEL_FLAT :
     mat_shininess[0] = 0.0f;
     p = mat_diffuse_color; *p++ = +0.0; *p++ = +0.0; *p++ = +0.0; *p++ = _alpha; /* no diffusion == flat */
     p = mat_emission; *p++ = 0.0; *p++ = 0.0; *p++ = 0.0; *p++ = _alpha;
     glMaterialfv(GL_FRONT_AND_BACK /* GL_FRONT */, GL_DIFFUSE, mat_diffuse_color);
     glColorMaterial(GL_FRONT_AND_BACK /* GL_FRONT */, GL_EMISSION); /* now glColor changes emission color */
     break;
   case COLOR_MODEL_INVERSE :
     mat_shininess[0] = 4.0f;
     p = mat_diffuse_color; *p++ = -.9; *p++ = -.9; *p++ = -.9; *p++ = _alpha; /* negative diffusion increase 
                                                                                * light away from normal */
     p = mat_emission; *p++ = -0.0; *p++ = -0.0; *p++ = -0.0; *p++ = _alpha;
     glMaterialfv(GL_FRONT_AND_BACK /* GL_FRONT */, GL_DIFFUSE, mat_diffuse_color);
     glColorMaterial(GL_FRONT_AND_BACK /* GL_FRONT */, GL_EMISSION); /* now glColor changes emission color */
     break;
   case COLOR_MODEL_GLOW : 
     p = mat_ambient_refl; *p++ = 0.65; *p++ = 0.65; *p++ = 0.65; *p++ = _alpha;
     p = mat_diffuse_color; *p++ = -1.0; *p++ = -1.0; *p++ = -1.0; *p++ = _alpha; /* negative diffusion */
     mat_shininess[0] = 0.0f;
     p = _diffuse_light; *p++ = 0.65; *p++ = 0.65; *p++ = 0.65; *p++ = 1.0; 
     glMaterialfv(GL_FRONT_AND_BACK /* GL_FRONT */, GL_DIFFUSE, mat_diffuse_color);
     glMaterialfv(GL_FRONT_AND_BACK /* GL_FRONT */, GL_AMBIENT, mat_ambient_refl);
     glColorMaterial(GL_FRONT_AND_BACK /* GL_FRONT */, GL_EMISSION); /* now glColor changes emission color */
     glLightfv(GL_LIGHT0, GL_DIFFUSE, _diffuse_light);
     break;
   case COLOR_MODEL_NORMAL :
   default :
     p = mat_ambient_refl; *p++ = 0.0; *p++ = 0.0; *p++ = 0.0; *p++ = _alpha;
     p = mat_diffuse_color; *p++ = 0.5; *p++ = 0.5; *p++ = 0.5; *p++ = _alpha; 
     mat_shininess[0] = 0.0f;
     p = mat_emission; *p++ = 0.1; *p++ = 0.1; *p++ = 0.1; *p++ = _alpha;
     // glMaterialfv(GL_FRONT_AND_BACK /* GL_FRONT */, GL_DIFFUSE, mat_diffuse_color);
     // glMaterialfv(GL_FRONT_AND_BACK /* GL_FRONT */, GL_AMBIENT, mat_ambient_refl);
     glMaterialfv(GL_FRONT_AND_BACK, GL_SHININESS, mat_shininess);
     glColorMaterial(GL_FRONT_AND_BACK /* GL_FRONT */, GL_EMISSION); /* now glColor changes emission color */
     /**** directional light -- only with color model 0 ****/
     _light1_pos[0] = 0.67f*_worldsize.x;
     _light1_pos[1] = _worldsize.y;
     _light1_pos[2] = 0.67f*_worldsize.z;
     _light1_pos[3] = 1.0;
     GLfloat light1_attenuation = 0.1f;
     glLightfv(GL_LIGHT1, GL_DIFFUSE, _white_light);
     glLightfv(GL_LIGHT1, GL_POSITION, _light1_pos);
     glLightfv(GL_LIGHT1, GL_LINEAR_ATTENUATION, &light1_attenuation);
     glEnable(GL_LIGHT1);
     /***************************/
  }

  glEnable(GL_LIGHTING);
  glEnable(GL_LIGHT0);
  glEnable(GL_RESCALE_NORMAL);
  glEnable(GL_DEPTH_TEST);
  glEnable(GL_COLOR_MATERIAL);
  glEnable(GL_MULTISAMPLE);

  /* transparency */
  if ( _blend )
  { 
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
  }

  /* Enable stencil operations */
  glEnable(GL_STENCIL_TEST);
  glStencilOp(GL_KEEP, GL_KEEP, GL_REPLACE);

  glMaterialfv(GL_FRONT_AND_BACK /* GL_FRONT */, GL_SHININESS, mat_shininess);
  glMaterialfv(GL_FRONT_AND_BACK /* GL_FRONT */, GL_EMISSION, mat_emission);
}

void dead_cell_geometry(void)
{
   triangle *trgl = (triangle*)malloc(8*sizeof(triangle));

   GLfloat a = 1.0;

   /* triangle 0 */
   trgl[0].v1[0] =  a/2.0; trgl[0].v1[1] = -a/2.0; trgl[0].v1[2] = 0.0;
   trgl[0].v2[0] =  0.0;   trgl[0].v2[1] = 0.0;    trgl[0].v2[2] = a/sqrtf(2.0);
   trgl[0].v3[0] = -a/2.0; trgl[0].v3[1] = -a/2.0; trgl[0].v3[2] = 0.0;
   normal3v(trgl[0].v1, trgl[0].v2, trgl[0].v3, trgl[0].normal);

   /* triangle 1 */
   trgl[1].v1[0] =  a/2.0; trgl[1].v1[1] =  a/2.0; trgl[1].v1[2] = 0.0;
   trgl[1].v2[0] =  0.0;   trgl[1].v2[1] = 0.0;    trgl[1].v2[2] = a/sqrtf(2.0);
   trgl[1].v3[0] =  a/2.0; trgl[1].v3[1] = -a/2.0; trgl[1].v3[2] = 0.0;
   normal3v(trgl[1].v1, trgl[1].v2, trgl[1].v3, trgl[1].normal);

   /* triangle 2 */
   trgl[2].v1[0] =  -a/2.0; trgl[2].v1[1] =  a/2.0; trgl[2].v1[2] = 0.0;
   trgl[2].v2[0] =  0.0;    trgl[2].v2[1] = 0.0;    trgl[2].v2[2] = a/sqrtf(2.0);
   trgl[2].v3[0] =  a/2.0;  trgl[2].v3[1] =  a/2.0; trgl[2].v3[2] = 0.0;
   normal3v(trgl[2].v1, trgl[2].v2, trgl[2].v3, trgl[2].normal);

   /* triangle 3 */
   trgl[3].v1[0] =  -a/2.0; trgl[3].v1[1] = -a/2.0; trgl[3].v1[2] = 0.0;
   trgl[3].v2[0] =  0.0;    trgl[3].v2[1] = 0.0;    trgl[3].v2[2] = a/sqrtf(2.0);
   trgl[3].v3[0] = -a/2.0;  trgl[3].v3[1] =  a/2.0; trgl[3].v3[2] = 0.0;
   normal3v(trgl[3].v1, trgl[3].v2, trgl[3].v3, trgl[3].normal);

   /* triangle 4 */
   trgl[4].v1[0] = -a/2.0; trgl[4].v1[1] = -a/2.0; trgl[4].v1[2] = 0.0;
   trgl[4].v2[0] =  0.0;   trgl[4].v2[1] = 0.0;    trgl[4].v2[2] = -a/sqrtf(2.0);
   trgl[4].v3[0] =  a/2.0; trgl[4].v3[1] = -a/2.0; trgl[4].v3[2] = 0.0;
   normal3v(trgl[4].v1, trgl[4].v2, trgl[4].v3, trgl[4].normal);

   /* triangle 5 */
   trgl[5].v1[0] =  a/2.0; trgl[5].v1[1] = -a/2.0; trgl[5].v1[2] = 0.0;
   trgl[5].v2[0] =  0.0;   trgl[5].v2[1] = 0.0;    trgl[5].v2[2] = -a/sqrtf(2.0);
   trgl[5].v3[0] =  a/2.0; trgl[5].v3[1] =  a/2.0; trgl[5].v3[2] = 0.0;
   normal3v(trgl[5].v1, trgl[5].v2, trgl[5].v3, trgl[5].normal);

   /* triangle 6 */
   trgl[6].v1[0] =  a/2.0;  trgl[6].v1[1] =  a/2.0; trgl[6].v1[2] = 0.0;
   trgl[6].v2[0] =  0.0;    trgl[6].v2[1] = 0.0;    trgl[6].v2[2] = -a/sqrtf(2.0);
   trgl[6].v3[0] =  -a/2.0; trgl[6].v3[1] =  a/2.0; trgl[6].v3[2] = 0.0;
   normal3v(trgl[6].v1, trgl[6].v2, trgl[6].v3, trgl[6].normal);

   /* triangle 7 */
   trgl[7].v1[0] =  -a/2.0; trgl[7].v1[1] =  a/2.0; trgl[7].v1[2] = 0.0;
   trgl[7].v2[0] =  0.0;    trgl[7].v2[1] = 0.0;    trgl[7].v2[2] = -a/sqrtf(2.0);
   trgl[7].v3[0] = -a/2.0;  trgl[7].v3[1] = -a/2.0; trgl[7].v3[2] = 0.0;
   normal3v(trgl[7].v1, trgl[7].v2, trgl[7].v3, trgl[7].normal);

   triangularize(&trgl,8); 

   glNewList(_cell_dlist + GEOM_DEAD, GL_COMPILE);
   glBegin( GL_TRIANGLES );
   for ( int i = 0*6; i < 8*6; i++ )
   {
     glNormal3f(trgl[i].normal[0],trgl[i].normal[1],trgl[i].normal[2]);
     glVertex3f(trgl[i].v1[0],trgl[i].v1[1],trgl[i].v1[2]);
     glVertex3f(trgl[i].v2[0],trgl[i].v2[1],trgl[i].v2[2]);
     glVertex3f(trgl[i].v3[0],trgl[i].v3[1],trgl[i].v3[2]);
   }
   glEnd();
   glEndList();

#if 0 
   for ( int i = 0*6; i < 8*6; i++ )
   {
     printf("triangle %d\n normal: %f %f %f\n",i,trgl[i].normal[0],trgl[i].normal[1],trgl[i].normal[2]);
     printf(" coords: %f %f %f\n", trgl[i].v1[0],trgl[i].v1[1],trgl[i].v1[2]);
     printf("         %f %f %f\n", trgl[i].v2[0],trgl[i].v2[1],trgl[i].v2[2]);
     printf("         %f %f %f\n", trgl[i].v3[0],trgl[i].v3[1],trgl[i].v3[2]);
   }
#endif

   free(trgl);

}

void triangularize(triangle **trgl, size_t nbr_triangles)
{
   

   /*
       1_  
       |  \_
       |    \3_
       |     | \__
       4   6 |   __\0
       |     | _-
       |  _ /5/
       2/          
   */
  *trgl = (triangle *)realloc(*trgl,6*nbr_triangles*sizeof(triangle));
  triangle t1;

  triangle *tr = NULL;

  for ( size_t i = nbr_triangles; i-- ; )
  {
    copy_triangle((*trgl)[i],*trgl + 6*i);
  }
  

  for ( size_t i = 0; i < nbr_triangles ; i++ )
  {
    tr = *trgl + 6*i;
    copy_triangle(*tr, &t1);

    /* tr 0  right */
    tr[0].v1[0] = t1.v1[0]; 
    tr[0].v1[1] = t1.v1[1]; 
    tr[0].v1[2] = t1.v1[2]; /* x1, y1, z1 */
    tr[0].v2[0] = 0.5*(t1.v1[0] + t1.v2[0]); /* x2 */
    tr[0].v2[1] = 0.5*(t1.v1[1] + t1.v2[1]); /* y2 */
    tr[0].v2[2] = 0.5*(t1.v1[2] + t1.v2[2]); /* z2 */
    tr[0].v3[0] = 0.5*(t1.v1[0] + t1.v3[0]); /* x3 */
    tr[0].v3[1] = 0.5*(t1.v1[1] + t1.v3[1]); /* y3 */
    tr[0].v3[2] = 0.5*(t1.v1[2] + t1.v3[2]); /* z3 */
    assign3v(t1.normal, 1, tr[0].normal, 3);

    /* tr 1 upper left */
    tr[1].v1[0] = tr[0].v2[0]; /* x1 */ 
    tr[1].v1[1] = tr[0].v2[1]; /* y1 */
    tr[1].v1[2] = tr[0].v2[2]; /* z1 */
    tr[1].v2[0] = t1.v2[0]; /* x2 */
    tr[1].v2[1] = t1.v2[1]; /* y2 */
    tr[1].v2[2] = t1.v2[2]; /* z2 */
    tr[1].v3[0] = 0.5*(t1.v2[0] + t1.v3[0]); /* x3 */
    tr[1].v3[1] = 0.5*(t1.v2[1] + t1.v3[1]); /* y3 */
    tr[1].v3[2] = 0.5*(t1.v2[2] + t1.v3[2]); /* z3 */
    assign3v(t1.normal, 1, tr[1].normal, 3);

    /* tr 2 lower left */
    tr[2].v1[0] = tr[0].v3[0];  /* x1 */ 
    tr[2].v1[1] = tr[0].v3[1];  /* y1 */
    tr[2].v1[2] = tr[0].v3[2];  /* z1 */
    tr[2].v2[0] = tr[1].v3[0]; /* x2 */
    tr[2].v2[1] = tr[1].v3[1]; /* y2 */
    tr[2].v2[2] = tr[1].v3[2]; /* z2 */
    tr[2].v3[0] = t1.v3[0]; /* x3 */
    tr[2].v3[1] = t1.v3[1]; /* y3 */
    tr[2].v3[2] = t1.v3[2]; /* z3 */
    assign3v(t1.normal, 1, tr[2].normal, 3);

    /* tr 4 mid-right */
    GLfloat vz = sqrt(11*distsq3v(tr[0].v1, tr[0].v3)/16);
    GLfloat center[3] = {(t1.v1[0]+t1.v2[0]+t1.v3[0])/3.0f, \
                         (t1.v1[1]+t1.v2[1]+t1.v3[1])/3.0f, \
                         (t1.v1[2]+t1.v2[2]+t1.v3[2])/3.0f};
    /* printf(" center of triangle %zu: %f %f %f\n", i, center[0], center[1], center[2]); */

    tr[3].v1[0] = tr[0].v2[0];       /* x1 */ 
    tr[3].v1[1] = tr[0].v2[1];       /* y1 */
    tr[3].v1[2] = tr[0].v2[2];       /* z1 */
    tr[3].v2[0] = center[0] + vz*tr[0].normal[0]; /* x2 */
    tr[3].v2[1] = center[1] + vz*tr[0].normal[1]; /* y2 */
    tr[3].v2[2] = center[2] + vz*tr[0].normal[2]; /* z2 */
    tr[3].v3[0] = tr[0].v3[0];       /* x3 */
    tr[3].v3[1] = tr[0].v3[1];       /* y3 */
    tr[3].v3[2] = tr[0].v3[2];       /* z3 */
    normal3v(tr[3].v1, tr[3].v2, tr[3].v3, tr[3].normal);

    /* tr 4 mid-upper-left */
    tr[4].v1[0] = tr[1].v3[0]; /* x2 */
    tr[4].v1[1] = tr[1].v3[1]; /* y2 */
    tr[4].v1[2] = tr[1].v3[2]; /* z2 */
    tr[4].v2[0] = tr[3].v2[0]; /* x1 */
    tr[4].v2[1] = tr[3].v2[1]; /* y1 */
    tr[4].v2[2] = tr[3].v2[2]; /* z1 */
    tr[4].v3[0] = tr[0].v2[0];  /* x3 */ 
    tr[4].v3[1] = tr[0].v2[1];  /* y3 */
    tr[4].v3[2] = tr[0].v2[2];  /* z3 */
    normal3v(tr[4].v1, tr[4].v2, tr[4].v3, tr[4].normal);

    /* tr 5 mid-lower-left */
    tr[5].v1[0] = tr[0].v3[0];  /* x1 */ 
    tr[5].v1[1] = tr[0].v3[1];  /* y1 */
    tr[5].v1[2] = tr[0].v3[2];  /* z1 */
    tr[5].v2[0] = tr[3].v2[0]; /* x2 */
    tr[5].v2[1] = tr[3].v2[1]; /* y2 */
    tr[5].v2[2] = tr[3].v2[2]; /* z2 */
    tr[5].v3[0] = tr[1].v3[0]; /* x3 */
    tr[5].v3[1] = tr[1].v3[1]; /* y3 */
    tr[5].v3[2] = tr[1].v3[2]; /* z3 */
    normal3v(tr[5].v1, tr[5].v2, tr[5].v3, tr[5].normal);

  }

}

void arrow_geometry()
{
  /* display list to store the geometry of an arrow (3D unit vector = (1,0,0)) */
  glNewList(_cell_dlist + GEOM_ARROW, GL_COMPILE);
  glEnable(GL_LINE_SMOOTH);
  glLineWidth(1.0);
  glBegin(GL_LINES);
  glVertex3f(0.0,0.0,0.0);
  glVertex3f(0.98,0.0,0.0);
  glEnd(); /* endGL_LINES */
  glBegin(GL_TRIANGLE_FAN); /* arrow tip */
  glVertex3f(1.0,0.0,0.0);
  for ( int ii = 0; ii < 17; ii++)
  {
   glVertex3f(0.90,0.1*cos(2.0*M_PI*(GLfloat)ii/16.0),0.1*sin(2.0*M_PI*(GLfloat)ii/16.0));
  }
  glEnd(); /* end GL_TRIANGLE_FAN */
  glBegin(GL_TRIANGLE_FAN); /* arrow butt */
  glVertex3f(0.0,0.0,0.0);
  for ( int ii = 0; ii < 17; ii++)
  {
   glVertex3f(0.0,0.1*cos(2.0*M_PI*(GLfloat)ii/16.0),0.1*sin(2.0*M_PI*(GLfloat)ii/16.0));
  }
  glEnd(); /* end GL_TRIANGLE_FAN */
  glEndList();
}

void scalebar_geometry()
{
  /* typical animal cell is 10-20µm in diameter
   * cell radius is between 0.5 and 1.0 world unit
   * or, between 1.0 and 2.0 units in diameter
   * We set the default scale to 10µm per world unit
   */
  glNewList(_cell_dlist + GEOM_SCALEBAR, GL_COMPILE);
  glEnable(GL_LINE_SMOOTH);
  glLineWidth(4.0);
  glBegin(GL_LINES);
  glVertex3f(0.0,0.0,0.0);
  glVertex3f(5.0,0.0,0.0); /* 5.0 <=> 50µm */ 
  glVertex3f(0.0,0.0,0.0);
  glVertex3f(0.0,5.0,0.0); /* 5.0 <=> 50µm */ 
  glVertex3f(0.0,0.0,0.0);
  glVertex3f(0.0,0.0,5.0); /* 5.0 <=> 50µm */ 
  glEnd(); /* end GL_LINES */

#if 0
  glLineWidth(1.0);
  glTranslatef(0.0,1.0,0.0);
  glScalef(0.5/_worldsize.x,0.5/_worldsize.x,0.5/_worldsize.z);
  for (size_t i = 0; i < strlen(scale); i++) 
    glutStrokeCharacter(GLUT_STROKE_ROMAN, scale[i]);
#endif

  glEndList();

}

void sphere_section(GLfloat /* -1 < s < 1 */ s )
{
  const GLfloat s_angle = M_PI/_s_res;
  GLfloat b,b2,r,r2; 
  GLfloat n1; /* 1/norm of normal */
  int last = 0;
  b = cos(s_angle);
  if ( s > b )
    b = s;
  r = sqrt(1 /* hypot */ - \
                   b*b /* height */);
  glBegin(GL_TRIANGLE_FAN); /* North Pole */
  n1 = sqrt((1 - b)/2.0);
  glNormal3f( n1*cos(s_angle),
              n1*sin(s_angle),
              n1*r/(1 - b));
  glVertex3f(0.0f,0.0f,1.0f);
  glVertex3f(r, 0.0f, b);
  for ( int ii = 0; ii < _s_res ; ii++)
  {
    glVertex3f( r*cos(2.0*s_angle*(GLfloat)(ii+1)),
                r*sin(2.0*s_angle*(GLfloat)(ii+1)),
                b);
    glNormal3f( n1*cos(s_angle*(GLfloat)(2*ii+1)),
                n1*sin(s_angle*(GLfloat)(2*ii+1)),
                n1*r/(1 - b));
  }
  glEnd(); /* end North Pole */
  for ( int jj = 0; jj < _s_res - 2; jj++ )
  {
    b =  cos((jj+1)*s_angle); /* 1 - (jj+1)*stack_height; */
    b2 = cos((jj+2)*s_angle); /* 1 - (jj+2)*stack_height; */
    if ( s > b2 )
    {
      b2 = s;
      last = 1;
    }
    r  = sqrt( 1 - b*b );
    r2 = sqrt( 1 - b2*b2 );
    glBegin(GL_QUAD_STRIP); /* Belt */
    n1 = (b-b2)/sqrt( (b-b2)*(b-b2) + (r-r2)*(r-r2) );
    glNormal3f( n1*cos(s_angle),
                n1*sin(s_angle),
                n1*(r2 - r)/(b - b2));
    glVertex3f( r, 0.0f, b);
    glVertex3f( r2, 0.0f, b2);
    for ( int ii = 0; ii < _s_res; ii++)
    {
      glVertex3f( r*cos(2.0*s_angle*(GLfloat)(ii+1)),
                  r*sin(2.0*s_angle*(GLfloat)(ii+1)),
                  b);
      glVertex3f( r2*cos(2.0*s_angle*(GLfloat)(ii+1)),
                  r2*sin(2.0*s_angle*(GLfloat)(ii+1)),
                  b2);
      glNormal3f( n1*cos(s_angle*(GLfloat)(2*ii+1)),
                  n1*sin(s_angle*(GLfloat)(2*ii+1)),
                  n1*(r2 - r)/(b - b2));
    }
    glEnd(); /* Belt */
    if ( last ) break;
  }
  glBegin(GL_TRIANGLE_FAN); /* South section */
  glNormal3f(0.0f, 0.0f, -1.0f);
  glVertex3f(0.0f, 0.0f, b2);
  glVertex3f(r2, 0.0f, b2);
  for ( int ii = 0; ii < _s_res; ii++)
  {
    glVertex3f( r2*cos(2.0*M_PI*(GLfloat)(ii+1)/(_s_res)),
                r2*sin(2.0*M_PI*(GLfloat)(ii+1)/(_s_res)),
                b2);
    glNormal3f( 0.0f, 0.0f, -1.0f);
  }
  glEnd(); /* end South section */
  glBegin(GL_TRIANGLE_FAN); /* South section */
  glColor3f(0.f,0.f,0.f);
  glNormal3f(0.0f, 0.0f, -1.0f);
  glVertex3f(0.0f, 0.0f, b2-0.01f);
  glVertex3f(r2/3.0, 0.0f, b2-0.01f);
  for ( int ii = 0; ii < _s_res; ii++)
  {
    glVertex3f( r2/3.0*cos(2.0*M_PI*(GLfloat)(ii+1)/(_s_res)),
                r2/3.0*sin(2.0*M_PI*(GLfloat)(ii+1)/(_s_res)),
                b2-0.01f);
    glNormal3f( 0.0f, 0.0f, -1.0f);
  }
  glEnd(); /* end South section */
#if 0
  glBegin(GL_LINE_STRIP);
  glLineWidth(2.0);
  glColor3f(0.f,0.f,0.f);
  for ( int ii = 0; ii < _s_res + 1; ii++)
  {
    glVertex3f( r2*cos(2.0*M_PI*(GLfloat)(ii)/(_s_res)),
                r2*sin(2.0*M_PI*(GLfloat)(ii)/(_s_res)),
                b2);
  }
  glEnd();
#endif
}

void standstill()
{
  /* void standstill()
   * rewinds the trajectory file by one frame, 
   * so that the next time the display is updated,
   * it uses the same data.
   */ 
    if (_fpos.i > 0) 
    { 
      _fpos.i--;
    }
    else 
    { 
      _fpos.i = 0;
    }
    fseek(_traj_file, _fpos.pos[_fpos.i], SEEK_SET);      
}

void draw_plane()
{

  float world_size[3] = {_worldsize.x, _worldsize.y, _worldsize.z};

  /******************* Draw positioning plane *******************/
  glDisable(GL_LIGHTING);
  if ( _move_plane ) 
  {
    glPolygonMode(GL_FRONT_AND_BACK, GL_FILL); // Mesh
    switch(_pos_plane)
    {
      case 0:
        glBegin(GL_QUADS); /* x-plane */
        glColor4fv (_positioning_plane_color[0]);
        glVertex3d(_x_plane, 0.0, 0.0);
        glVertex3d(_x_plane, 0.0, world_size[2]);
        glVertex3d(_x_plane, world_size[1], world_size[2]);
        glVertex3d(_x_plane, world_size[1], 0.0); 
        glEnd(); /* end x-plane */
        break;
      case 1:
        glBegin(GL_QUADS); /* y-plane */
        glColor4fv (_positioning_plane_color[1]);
        glVertex3d(0.0, _y_plane, 0.0);
        glVertex3d(0.0, _y_plane, world_size[2]);
        glVertex3d(world_size[0], _y_plane, world_size[2]);
        glVertex3d(world_size[0], _y_plane, 0.0); 
        glEnd(); /* end y-plane */
        break;
      case 2:
        glBegin(GL_QUADS); /* z-plane */
        glColor4fv (_positioning_plane_color[2]);
        glVertex3d(0.0, 0.0, _z_plane);
        glVertex3d(world_size[0], 0.0, _z_plane);
        glVertex3d(world_size[0], world_size[1], _z_plane);
        glVertex3d(0.0, world_size[1], _z_plane);
        glEnd(); /* end z-plane */
        break;
    }
  }
  glColor3f(0.4,0.5,0.7);
  glLineWidth(2.0);
  glBegin(GL_LINES);

  /* draw x-axis */
  glVertex3d(0.0, _y_plane, _z_plane);
  glVertex3d(world_size[0], _y_plane, _z_plane);

  /* draw y-axis */
  glVertex3d(_x_plane, 0.0, _z_plane);
  glVertex3d(_x_plane, world_size[1], _z_plane);

  /* draw z-axis */
  glVertex3d(_x_plane, _y_plane, 0.0);
  glVertex3d(_x_plane, _y_plane, world_size[2]);
  glEnd();

  glLineWidth(2.0);
  glBegin(GL_LINES);
  /* x-ticks */
  for ( size_t i = 0; i < world_size[0] ; i++ )
  {
    glVertex3d((GLfloat)i, _y_plane, _z_plane);
    glVertex3d((GLfloat)i, _y_plane - 0.02*world_size[1], _z_plane);
  }

  /* y-ticks */
  for ( size_t i = 0; i < world_size[1] ; i++ )
  {
    glVertex3d(_x_plane, (GLfloat)i, _z_plane);
    glVertex3d(_x_plane + 0.02*world_size[0], (GLfloat)i, _z_plane);
  }

  /* z-ticks */
  for ( size_t i = 0; i < world_size[2] ; i++ )
  {
    glVertex3d(_x_plane, _y_plane, (GLfloat)i);
    glVertex3d(_x_plane, _y_plane - 0.02*world_size[1], (GLfloat)i);
  }

  glEnd();

  glEnable(GL_LIGHTING);

}


void move_plane(GLint x, GLint y)
{

  int window_height = glutGet(GLUT_WINDOW_HEIGHT);
  GLfloat depth;
  GLdouble model[16];
  GLdouble projection[16];
  GLint view[4];
  GLdouble new_x, new_y, new_z;
  GLfloat delta;
  static GLdouble old_x, old_y, old_z;

  glReadPixels(x, window_height - y - 1, 1, 1, GL_DEPTH_COMPONENT, GL_FLOAT, &depth);

  glGetDoublev(GL_MODELVIEW_MATRIX,model);
  glGetDoublev(GL_PROJECTION_MATRIX,projection);
  glGetIntegerv(GL_VIEWPORT,view);

  /* get approximate nearest object coordinates */
  gluUnProject((GLdouble)x, (GLdouble)(window_height - y - 1),
      0.0,model,projection,view,&new_x,&new_y,&new_z);

  if ( _action == MOVE_PLANE ) delta = 0.05;
  else delta = 0.005;

  switch(_pos_plane)
  {
    case 0:
      _x_plane += delta*_worldsize.x*sign(new_x - old_x)*(fabs(new_x - old_x) < delta*_worldsize.x);
      clamp(&_x_plane,0,_worldsize.x);
      break;
    case 1:
      _y_plane += delta*_worldsize.y*sign(new_y - old_y)*(fabs(new_y - old_y) < delta*_worldsize.y);
      clamp(&_y_plane,0,_worldsize.y);
      break;
    case 2:
      _z_plane += delta*_worldsize.z*sign(new_z - old_z)*(fabs(new_z - old_z) < delta*_worldsize.z);
      clamp(&_z_plane,0,_worldsize.z);
      break;
  }


  old_x = new_x;
  old_y = new_y;
  old_z = new_z;

  glutPostRedisplay();
  
}

void save_image(const char* filename, const char* file_type)
{
    int window_height = glutGet(GLUT_WINDOW_HEIGHT);
    int window_width  = glutGet(GLUT_WINDOW_WIDTH);
    char *buffer;
    FILE *fid;
    GLsizei _nr_channels = 3;         
    unsigned int _buffer_size = window_width * window_height * _nr_channels; /* size in bytes */
    buffer = (char*)malloc(_buffer_size * sizeof(char));
    glPixelStorei(GL_PACK_ALIGNMENT, 1);
    glReadBuffer(GL_FRONT);
    glReadPixels(0, 0, window_width, window_height, GL_RGB, GL_UNSIGNED_BYTE, buffer);

    if (strcmp(file_type, "png") == 0)
    {
        if ((fid = fopen(filename, "wb")) == NULL)
        {
            fprintf(stderr, "ERROR : could not open file '%s'.\n", filename);
            free(buffer);
            return;
        }

        png_structp png_ptr = png_create_write_struct(PNG_LIBPNG_VER_STRING, NULL, NULL, NULL);
        if (!png_ptr)
        {
            fprintf(stderr, "ERROR : could not create PNG write struct.\n");
            fclose(fid);
            free(buffer);
            return;
        }

        png_infop info_ptr = png_create_info_struct(png_ptr);
        if (!info_ptr)
        {
            fprintf(stderr, "ERROR : could not create PNG info struct.\n");
            png_destroy_write_struct(&png_ptr, (png_infopp)NULL);
            fclose(fid);
            free(buffer);
            return;
        }

        if (setjmp(png_jmpbuf(png_ptr)))
        {
            fprintf(stderr, "ERROR : error during PNG creation.\n");
            png_destroy_write_struct(&png_ptr, &info_ptr);
            fclose(fid);
            free(buffer);
            return;
        }

        png_init_io(png_ptr, fid);
        png_set_IHDR(png_ptr, info_ptr, window_width, window_height, 8, PNG_COLOR_TYPE_RGB, PNG_INTERLACE_NONE, PNG_COMPRESSION_TYPE_DEFAULT, PNG_FILTER_TYPE_DEFAULT);
        png_write_info(png_ptr, info_ptr);

        for (int y = 0; y < window_height; y++)
        {
            png_write_row(png_ptr, (png_bytep)(buffer + (window_height - y - 1) * window_width * _nr_channels));
        }

        png_write_end(png_ptr, NULL);
        png_destroy_write_struct(&png_ptr, &info_ptr);
        fclose(fid);
    }
    else
    {
        if ((fid = fopen(filename, "wb")) == NULL)
        {
            fprintf(stderr, "ERROR : could not open file '%s'.\n", filename);
            free(buffer);
            return;
        }

        char mn[] = "BM";
        char creator[] = "SCVW";
        unsigned int offset = 0x36;  /* = 54: header=14 and DIB=40 */
        unsigned int dibsize = 0x28; /* = 40 */
        unsigned int zero = 0x0;     /* = 00 00 00 00 */
        unsigned int resolution = 2835; /* points per meter = 72dpi */
        unsigned int image_size = window_width * window_height * _nr_channels + offset; /* total image size on disk */

        /* bitmap header: 14 bytes */
        fwrite(mn, 1, 2, fid);             /* bitmap image type: BM */
        fwrite(&image_size, 4, 1, fid);    /* total image size */
        fwrite(creator, 1, 4, fid);        /* user-defined bytes */
        fwrite(&offset, 4, 1, fid);        /* total offset, including DIB */

        /* DIB - BITMAPINFOHEADER */
        fwrite(&dibsize, 4, 1, fid);          /* DIB size: 40 bytes - 4 bytes */
        fwrite(&window_width, 4, 1, fid);     /* image width, in pixels  - 4 bytes */
        fwrite(&window_height, 4, 1, fid);    /* image height, in pixels - 4 bytes */
        fputc(1, fid);                      /* number of color planes: 1 */ 
        fputc(0, fid);                      /* pad with zeros to fill 2 bytes - 2 bytes */
        fputc(0x18, fid);                   /* number of bits/pixel n = 24 */
        fputc(0, fid);                      /* pad with zeros to fill 2 bytes - 2 bytes */
        fwrite(&zero, 4, 1, fid);             /* compression mode: 0 = no compression - 4 bytes */
        fwrite(&zero, 4, 1, fid);             /* image size: 0 when no compression - 4 bytes */
        fwrite(&resolution, 4, 1, fid);       /* res horiz - 4 bytes */
        fwrite(&resolution, 4, 1, fid);       /* res vert  - 4 bytes */
        fwrite(&zero, 4, 1, fid);             /* default number colors in palette: 0 to defaults to n^24 - 4 bytes */
        fwrite(&zero, 4, 1, fid);             /* default important colors: 0 - 4 bytes */

        /* image raw data, not memory aligned */
        fwrite(buffer, 1, _buffer_size, fid);   /* write image raw data */
        fclose(fid);
    }

    fprintf(_log_file, "saved image: %s\n", filename);
    printf("image: %s\n", filename);
    printf("  window_size: %u x %u\n", window_width, window_height);
    printf("  buffer size: %u bytes\n", _buffer_size);

    /* clean up */
    free(buffer);
}
void save_image() {
    char default_filename[256];
    snprintf(default_filename, 256, "s%06lu_n%d_%s.bmp", _fpos.i-1, _popsize, _signal_name);
    save_image(default_filename, "bmp");
}
void print_snapshot()
{
  char filepath[64];
  FILE *fid;

  /* filepath: sSTEPNUMBER_nPOPSIZE_SIGNAL.txt */
  if ( _snapshot_file != NULL )
  {
    snprintf(filepath,64,"%s_%06lu_n%d_%s.txt",_snapshot_file,_fpos.i-1,_popsize,_signal_name);
  }
  else 
  {
    snprintf(filepath,64,"s%06lu_n%d_%s.txt",_fpos.i-1,_popsize,_signal_name);
  }
  if ( ( fid = fopen(filepath,"w") ) == NULL )
  {
        fprintf(stderr, "ERROR : could not open file '%s'.\n", filepath );
        return;
  }
  fprintf(_log_file,"saved snapshot: %s\n",filepath);
  printf("data: %s\nto extract from %s:\n\n",filepath,_traj_filename);

  printf("  awk '($4 < %g || $5 > %g || $6 < %g) && $1 == %g {print}' %s\n",
      _x_plane,_y_plane,_z_plane,_time,_traj_filename);
  fprintf(_log_file,"  awk '($4 < %g || $5 > %g || $6 < %g) && $1 == %g {print}' %s\n",
      _x_plane,_y_plane,_z_plane,_time,_traj_filename);

  fprintf(fid,"# excluding x: %g:%g & y: %g:%g & z: %g:%g,\n",
      _x_plane,_worldsize.x, 0.0, _y_plane, _z_plane, _worldsize.z );
  fprintf(fid,"# excluding signal %s outside range %g:%g\n",
      _signal_name, _min_clip_signal, _max_clip_signal);

  for ( int i = 0; i < _popsize ; i++ )
  {
    /* if _clipcell = true, do not print clipped cells */
    if ( _mycells[i].clip && _clipcell ) 
    {
      continue; /* cell is clipped, do not print */
    }

    /* do not print cells that are above x-plane, below y-plane
     * and above z-plane */
    if ( ( _mycells[i].position[0] > _x_plane ) &&
         ( _mycells[i].position[1] < _y_plane ) &&
         ( _mycells[i].position[2] > _z_plane ) )
    {
      continue; /* do not print the cell */
    }

    fprintf(fid,"%g %d %d %g %g %g %g %g\n",_time, _popsize, _mycells[i].id, 
        _mycells[i].position[0], _mycells[i].position[1], _mycells[i].position[2],
        _mycells[i].radius, _mycells[i].signal);
  }

  /* clean up */
  fclose(fid);
}

int get_lineage(int daughter)
{
  if ( _lineage_str == NULL )
    return -1;

  char daughter_str[32];
  char *pos;
  int mother = -1;
  char mother_str[32];
  snprintf(daughter_str,31,")%d>",daughter);
  pos = strstr(_lineage_str,daughter_str);
  if ( sscanf(pos,")%*d>%d",&mother) )
  {
    return mother;
  }
  else
  {
    sscanf(pos,")%*d>%[a-zA-Z]",mother_str);
    if ( strncmp(mother_str,"root",4) == 0 )
    {
      /* mother is 'root', return id 0 */
      return 0; 
    }
  }
  
  /* no mother found, return -1 */
  return mother;
}

void read_lineage_file()
{
  size_t linecap;
  int retval;
  char c;
  do {retval = fscanf(_lineage_file, "%c", &c); }  while ((c != '\n') && (retval != EOF));
  getline(&_lineage_str, &linecap,_lineage_file);
}

void on_exit() 
{
  fclose(_traj_file);
  fclose(_log_file);
  if ( _display_parameter_file != (FILE *)NULL ) fclose(_display_parameter_file);
  if ( _lineage_file != (FILE *)NULL ) fclose(_lineage_file);
  free(_mycells);
  free(_colormap);
  free(_norm_filename);
  free(_traj_filename);
  free(_signal_name);
  free(_tags);
  free(_fpos.pos);
  free(_ts.t);
  free(_ts.y);
  free(_ts.step);
  free(_stats_ts.t);
  free(_stats_ts.y);
  free(_stats_ts.y2);
  free(_stats_ts.step);
  if (_record_video) {cleanup_opencv_video();}
  free(_video_filename);
  printf("bye...\n");
}

void print_help( char* prog_name )
{
  printf( "\
Usage : " fmt(bb) "%s" fmt(plain) " -h\n\
   or : " fmt(bb) "%s" fmt(plain) " [" options_opt "] [" file_opt "]\n\n\
\t" fmt(bb) "-a " num_opt ", " fmt(bb) "--alpha=" num_opt         "\t\tSet alpha (transparency) to " num_opt ", from 0.0 (transparent)\n\
"                                                                 "\t\t\t\t\tto 1.0 (opaque), default value: 1.0)\n\n\
\t" fmt(bb) "-c " name_opt "," fmt(bb) "--colormap=" name_opt     "\t\tSet colormap to " name_opt " ({jet}, spring or neon)\n\
"                                                                 "\t\t\t\t\tor to the colormap defined in file " name_opt "\n\n\
\t" fmt(bb) "-g " name_opt ", " fmt(bb) "--tags=" name_opt        "\t\tColor cells that have cell type " name_opt "\n\
"                                                                 "\t\t\t\t\t(default color: green, see option tag-color to change color)\n\n\
\t" fmt(bb) "-h" fmt(plain) ", " fmt(bb) "--help" fmt(plain)      "\t\t\tDisplay this screen\n\n\
\t" fmt(bb) "-i " name_opt ", " fmt(bb) "--signal=" name_opt      "\t\tSet display signal to " name_opt " (default: first signal in file)\n\n\
\t" fmt(bb) "-m " num_opt ", " fmt(bb) "--colormodel=" num_opt    "\tSet colormodel to " num_opt ":\n\
"                                                                 "\t\t\t\t\t0: default, cells appear as hard balls.\n\
"                                                                 "\t\t\t\t\t1: inverse, negative diffusive light, cells appear glowing.\n\
"                                                                 "\t\t\t\t\t2: flat, no shading, cells appear flat.\n\
"                                                                 "\t\t\t\t\t3: glow, same as 1, but stronger effect.\n\n\
\t" fmt(bb) "-n " name_opt ", " fmt(bb) "--normalization=" name_opt "\tUse signal normalization file " name_opt " (default: normalization.txt)\n\n\
\t" fmt(bb) "-p" fmt(plain) ", " fmt(bb) "--pause "    fmt(plain) "\t\t\tLaunch in pause mode on (resume with space bar, default: no pause)\n\n\
\t" fmt(bb) "-s " num_opt ", " fmt(bb) "--speed="  num_opt        "\t\tSet playing speed to " num_opt ", in frames per ms (default: 0.02)\n\n\
\t" fmt(bb) "-v " num_opt ", " fmt(bb) "--viewangle=" num_opt     "\t\tSet perspective viewing angle to " num_opt ",\n\
"                                                                 "\t\t\t\t\t(betweeen 0.0 and 180.0 degrees, default: 45.0)\n\n\
\t" fmt(bb) "-V " num_opt ", " fmt(bb) "--viewpoint=" num_opt     "\t\tSet default viewpoint to " num_opt "\n\
"                                                                 "\t\t\t\t\t" num_opt "='azimuth inclination twist distance pan_x pan_y'\n\n\
\t" fmt(bb) "-w" fmt(plain) ", " fmt(bb) "--white "    fmt(plain) "\t\t\tSet white background (default: black)\n\
"                                                                 "\t\t\t\t\tUse option --background-color to set custom color.\n\n\
\t" fmt(bb) "-r " name_opt ", " fmt(bb) "--record=" name_opt      "\t\tRecord animation as MP4 video with prefix " name_opt "\n\
"                                                                 "\t\t\t\t\tOptional format: prefix:window (0=main, 1=aux, 2=stats)\n\
"                                                                 "\t\t\t\t\tExample: --record=video:1 records aux window\n\
\t" fmt(bb) "--animate="                file_opt                  "\t\t\tSpecify viewpoint, signal, guide and hide view from " file_opt ".\n\
"                                                                 "\t\t\t\t\t" file_opt " must have 6 to 11 whitespace-separated columns\n\
"                                                                 "\t\t\t\t\t'azimuth inclination twist distance pan_x pan_y [signal [guide [hide]]]'\n\
"                                                                 "\t\t\t\t\t" file_opt " is read line by line and the view is updated at each\n\
"                                                                 "\t\t\t\t\tstep, or at regular intervals if the view is pause\n\
"                                                                 "\t\t\t\t\tSee the file view.log for how to set the viewpoint.\n\
"                                                                 "\t\t\t\t\tsignal must be a valid signal name, guide is a number between 0 and 7,\n\
"                                                                 "\t\t\t\t\tand hide is a triplet of the coordinates of the hiding planes\n\n\
\t" fmt(bb) "--background-color="       num_opt                   "\t\tSet background color to " num_opt " (RGB triplet)\n\n\
\t" fmt(bb) "--guides="                 num_opt                   "\t\t\tSet border/hide/scale bar visual guides to " num_opt "\n\
"                                                                 "\t\t\t\t\tnum border hide   scale bar\n\
"                                                                 "\t\t\t\t\t 0   ---    ---    --- \n\
"                                                                 "\t\t\t\t\t 1   -X-    ---    --- \n\
"                                                                 "\t\t\t\t\t 2   ---    -X-    --- \n\
"                                                                 "\t\t\t\t\t 3   -X-    -X-    --- \n\
"                                                                 "\t\t\t\t\t 4   ---    ---    -X- \n\
"                                                                 "\t\t\t\t\t 5   -X-    ---    -X- \n\
"                                                                 "\t\t\t\t\t 6   ---    -X-    -X- \n\
"                                                                 "\t\t\t\t\t 7   -X-    -X-    -X- \n\n\
\t" fmt(bb) "--hide="                   num_opt                   "\t\t\tSet hide plane positions to " num_opt " ='x y z'\n\n\
\t" fmt(bb) "--hoovering"    fmt(plain)                           "\t\t\tEnable single-cell hoovering (default)\n\n\
\t" fmt(bb) "--no-hoovering" fmt(plain)                           "\t\t\tDisable hoovering (useful when visualisation refreshes slowly)\n\n\
\t" fmt(bb) "--lineage="                file_opt                  "\t\t\tUse lineage " file_opt " to highlight cell lineage\n\n\
	" fmt(bb) "--no-annotations"       fmt(plain)                "\t\t\tHide all annotations (legend, colorbar, timing, cell numbers)\n\n\
	" fmt(bb) "--colorbar-only"        fmt(plain)                "\t\t\tShow only the colorbar (hide legend, timing, cell numbers)\n\n\
	" fmt(bb) "--minimal-annotations"  fmt(plain)                "\t\t\tShow only time and number of cells (hide legend, colorbar, cell numbers, etc.)\n\n\
\t" fmt(bb) "--resolution-low" fmt(plain) ", " fmt(bb) "--resolution-medium" fmt(plain) ", " fmt(bb) "--resolution-high" fmt(plain) "\n\
"                                                                 "\t\t\t\t\tSets the number of stack and slice cell sphere rendering\n\
"                                                                 "\t\t\t\t\t(low: 5, medium: 12, high: 20, default: medium)\n\n\
\t" fmt(bb) "--shade-flat" fmt(plain) ", " fmt(bb) "--shade-smooth" fmt(plain) "\tSet light shading to flat or smooth\n\n\
\t" fmt(bb) "--scale-bar-position="     num_opt                   "\tSet scale bar at position " num_opt " ='x y z')\n\n\
\t" fmt(bb) "--snapshot="               name_opt                  "\t\t\tSave image with prefix " name_opt ", and quit\n\n\
\t" fmt(bb) "--step="                   num_opt                   "\t\t\tForward to time step " num_opt " before displaying the view\n\n\
\t" fmt(bb) "--tag-color="              num_opt                   "\t\t\tSet tag color to " num_opt " (RGBA quadruplet).\n\
"                                                                 "\t\t\t\t\t(Use with option --tags)\n\n\
\t" fmt(bb) "--window-size="            num_opt                   "\t\tSet the main window size to " num_opt " ='width height'\n\n\
\nView the simulation stored in " file_opt ", or in trajectory.txt if " file_opt " is missing\n\n\
Runtime commands\n\n\
\tmouse\n\
\t" fmt(bb) "left button      " fmt(plain) "\tRotate viewpoint\n\
\t" fmt(bb) "right button     " fmt(plain) "\tZoom in or zoom out\n\
\t" fmt(bb) "shift+left click " fmt(plain) "\tSelect cell\n\
\t\t\t\tShift click on the background to\n\
\t\t\t\tunselect the cell\n\n\
\tkeyboard\n\
\t" fmt(bb) "r" fmt(plain) ", " fmt(bb) "R" fmt(plain) "\t\tReset perspective\n\
\t" fmt(bb) "Space    " fmt(plain) "\tPause, resume animation\n\
\t" fmt(bb) "◀ (left) " fmt(plain) "\tPan left\n\
\t" fmt(bb) "▶ (right)" fmt(plain) "\tPan right\n\
\t" fmt(bb) "▲ (up)   " fmt(plain) "\tPan up\n\
\t" fmt(bb) "▼ (down) " fmt(plain) "\tPan down\n\
\t" fmt(bb) "alt+▲    " fmt(plain) "\tZoom in (same as right button)\n\
\t" fmt(bb) "alt+▼    " fmt(plain) "\tZoom out (same as right button)\n\
\t" fmt(bb) "shift+◀  " fmt(plain) "\tRotate left\n\
\t" fmt(bb) "shift+▶  " fmt(plain) "\tRotate right\n\
\t" fmt(bb) "shift+▲  " fmt(plain) "\tIncrease inclination\n\
\t" fmt(bb) "shift+▼  " fmt(plain) "\tDecrease inclination\n\
\t" fmt(bb) "1" fmt(plain) "-" fmt(bb) "9" fmt(plain) "\t\tDisplay signal 1 to 9\n\
\t" fmt(bb) "[" fmt(plain) ", " fmt(bb) "]" fmt(plain) "\t\tCycle forward, backward through signals\n\
\t" fmt(bb) "b        " fmt(plain) "\tBacktrack to previous timestep\n\
\t" fmt(bb) "B        " fmt(plain) "\tBacktrack 50 timesteps\n\
\t" fmt(bb) "c        " fmt(plain) "\tToggle clipping of cells with signal outside range\n\
\t" fmt(bb) "f        " fmt(plain) "\tReset the initial window size, position and view\n\
\t" fmt(bb) "F        " fmt(plain) "\tGo full screen\n\
\t" fmt(bb) "h" fmt(plain) ", " fmt(bb) "l" fmt(plain) "\t\tDecrease (h) or increase (l) min signal\n\
\t" fmt(bb) "j" fmt(plain) ", " fmt(bb) "k" fmt(plain) "\t\tDecrease (j) or increase (k) max signal\n\
\t" fmt(bb) "m        " fmt(plain) "\tToggle the hide planes\n\
\t" fmt(bb) "n        " fmt(plain) "\tPause at the next timestep\n\
\t" fmt(bb) "N        " fmt(plain) "\tSkip over 50 timesteps\n\
\t" fmt(bb) "s        " fmt(plain) "\tSave the current frame as a BMP image\n\
\t" fmt(bb) "t        " fmt(plain) "\tToggle single cell tracking. Upon toggling,\n\
\t" fmt(bb) "         " fmt(plain) "\tenter cell id to enable single cell tracking.\n\
\t" fmt(bb) "         " fmt(plain) "\tPress t to disable tracking\n\
\t" fmt(bb) "w        " fmt(plain) "\tCycle through border/hide/scale bar visual guides\n\
\t" fmt(bb) "x" fmt(plain) ", " fmt(bb) "y" fmt(plain) ", " fmt(bb) "z" fmt(plain) "\t\tSelect the hide plane\n\
\t" fmt(bb) "ESC" fmt(plain) ", " fmt(bb) "q" fmt(plain) ", " fmt(bb) "Q" fmt(plain) "\tQuit\n\n",\
   prog_name, prog_name );
}

/* legacy colormaps
 * Pre-defined colormaps are now generated with make_colormap_array 
 * based on the list of available colormaps defined in view.h: _list_colormaps
 */
void init_color_map_neon()
{
  _colormap_len = 64;
  _colormap = (rgbacolor *)malloc(_colormap_len*sizeof(rgbacolor));
  _colormap[0].c[0] = 0.; _colormap[0].c[1] = 0.2; _colormap[0].c[2] = 0.333333; _colormap[0].c[3] = 0.0;
_colormap[1].c[0] = 0.; _colormap[1].c[1] = 0.233333; _colormap[1].c[2] = 0.3375; _colormap[1].c[3] = 0.0;
_colormap[2].c[0] = 0.; _colormap[2].c[1] = 0.266667; _colormap[2].c[2] = 0.341667; _colormap[2].c[3] = 0.0;
_colormap[3].c[0] = 0.; _colormap[3].c[1] = 0.3; _colormap[3].c[2] = 0.345833; _colormap[3].c[3] = 0.0;
_colormap[4].c[0] = 0.; _colormap[4].c[1] = 0.333333; _colormap[4].c[2] = 0.35; _colormap[4].c[3] = 0.0;
_colormap[5].c[0] = 0.; _colormap[5].c[1] = 0.366667; _colormap[5].c[2] = 0.354167; _colormap[5].c[3] = 0.0;
_colormap[6].c[0] = 0.; _colormap[6].c[1] = 0.4; _colormap[6].c[2] = 0.358333; _colormap[6].c[3] = 0.0;
_colormap[7].c[0] = 0.; _colormap[7].c[1] = 0.433333; _colormap[7].c[2] = 0.3625; _colormap[7].c[3] = 0.0;
_colormap[8].c[0] = 0.; _colormap[8].c[1] = 0.466667; _colormap[8].c[2] = 0.366667; _colormap[8].c[3] = 0.0;
_colormap[9].c[0] = 0.; _colormap[9].c[1] = 0.5; _colormap[9].c[2] = 0.370833; _colormap[9].c[3] = 0.0;
_colormap[10].c[0] = 0.; _colormap[10].c[1] = 0.533333; _colormap[10].c[2] = 0.375; _colormap[10].c[3] = 0.0;
_colormap[11].c[0] = 0.; _colormap[11].c[1] = 0.566667; _colormap[11].c[2] = 0.379167; _colormap[11].c[3] = 0.0;
_colormap[12].c[0] = 0.; _colormap[12].c[1] = 0.6; _colormap[12].c[2] = 0.383333; _colormap[12].c[3] = 0.0;
_colormap[13].c[0] = 0.; _colormap[13].c[1] = 0.633333; _colormap[13].c[2] = 0.3875; _colormap[13].c[3] = 0.0;
_colormap[14].c[0] = 0.; _colormap[14].c[1] = 0.666667; _colormap[14].c[2] = 0.391667; _colormap[14].c[3] = 0.0;
_colormap[15].c[0] = 0.; _colormap[15].c[1] = 0.7; _colormap[15].c[2] = 0.395833; _colormap[15].c[3] = 0.0;
_colormap[16].c[0] = 0.; _colormap[16].c[1] = 0.733333; _colormap[16].c[2] = 0.4; _colormap[16].c[3] = 0.0;
_colormap[17].c[0] = 0.0125; _colormap[17].c[1] = 0.7375; _colormap[17].c[2] = 0.425; _colormap[17].c[3] = 0.0;
_colormap[18].c[0] = 0.025; _colormap[18].c[1] = 0.741667; _colormap[18].c[2] = 0.45; _colormap[18].c[3] = 0.0;
_colormap[19].c[0] = 0.0375; _colormap[19].c[1] = 0.745833; _colormap[19].c[2] = 0.475; _colormap[19].c[3] = 0.0;
_colormap[20].c[0] = 0.05; _colormap[20].c[1] = 0.75; _colormap[20].c[2] = 0.5; _colormap[20].c[3] = 0.0;
_colormap[21].c[0] = 0.0625; _colormap[21].c[1] = 0.754167; _colormap[21].c[2] = 0.525; _colormap[21].c[3] = 0.0;
_colormap[22].c[0] = 0.075; _colormap[22].c[1] = 0.758333; _colormap[22].c[2] = 0.55; _colormap[22].c[3] = 0.0;
_colormap[23].c[0] = 0.0875; _colormap[23].c[1] = 0.7625; _colormap[23].c[2] = 0.575; _colormap[23].c[3] = 0.0;
_colormap[24].c[0] = 0.1; _colormap[24].c[1] = 0.766667; _colormap[24].c[2] = 0.6; _colormap[24].c[3] = 0.0;
_colormap[25].c[0] = 0.1125; _colormap[25].c[1] = 0.770833; _colormap[25].c[2] = 0.625; _colormap[25].c[3] = 0.0;
_colormap[26].c[0] = 0.125; _colormap[26].c[1] = 0.775; _colormap[26].c[2] = 0.65; _colormap[26].c[3] = 0.0;
_colormap[27].c[0] = 0.1375; _colormap[27].c[1] = 0.779167; _colormap[27].c[2] = 0.675; _colormap[27].c[3] = 0.0;
_colormap[28].c[0] = 0.15; _colormap[28].c[1] = 0.783333; _colormap[28].c[2] = 0.7; _colormap[28].c[3] = 0.0;
_colormap[29].c[0] = 0.1625; _colormap[29].c[1] = 0.7875; _colormap[29].c[2] = 0.725; _colormap[29].c[3] = 0.0;
_colormap[30].c[0] = 0.175; _colormap[30].c[1] = 0.791667; _colormap[30].c[2] = 0.75; _colormap[30].c[3] = 0.0;
_colormap[31].c[0] = 0.1875; _colormap[31].c[1] = 0.795833; _colormap[31].c[2] = 0.775; _colormap[31].c[3] = 0.0;
_colormap[32].c[0] = 0.2; _colormap[32].c[1] = 0.8; _colormap[32].c[2] = 0.8; _colormap[32].c[3] = 0.0;
_colormap[33].c[0] = 0.241667; _colormap[33].c[1] = 0.779167; _colormap[33].c[2] = 0.808333; _colormap[33].c[3] = 0.0;
_colormap[34].c[0] = 0.283333; _colormap[34].c[1] = 0.758333; _colormap[34].c[2] = 0.816667; _colormap[34].c[3] = 0.0;
_colormap[35].c[0] = 0.325; _colormap[35].c[1] = 0.7375; _colormap[35].c[2] = 0.825; _colormap[35].c[3] = 0.0;
_colormap[36].c[0] = 0.366667; _colormap[36].c[1] = 0.716667; _colormap[36].c[2] = 0.833333; _colormap[36].c[3] = 0.0;
_colormap[37].c[0] = 0.408333; _colormap[37].c[1] = 0.695833; _colormap[37].c[2] = 0.841667; _colormap[37].c[3] = 0.0;
_colormap[38].c[0] = 0.45; _colormap[38].c[1] = 0.675; _colormap[38].c[2] = 0.85; _colormap[38].c[3] = 0.0;
_colormap[39].c[0] = 0.491667; _colormap[39].c[1] = 0.654167; _colormap[39].c[2] = 0.858333; _colormap[39].c[3] = 0.0;
_colormap[40].c[0] = 0.533333; _colormap[40].c[1] = 0.633333; _colormap[40].c[2] = 0.866667; _colormap[40].c[3] = 0.0;
_colormap[41].c[0] = 0.575; _colormap[41].c[1] = 0.6125; _colormap[41].c[2] = 0.875; _colormap[41].c[3] = 0.0;
_colormap[42].c[0] = 0.616667; _colormap[42].c[1] = 0.591667; _colormap[42].c[2] = 0.883333; _colormap[42].c[3] = 0.0;
_colormap[43].c[0] = 0.658333; _colormap[43].c[1] = 0.570833; _colormap[43].c[2] = 0.891667; _colormap[43].c[3] = 0.0;
_colormap[44].c[0] = 0.7; _colormap[44].c[1] = 0.55; _colormap[44].c[2] = 0.9; _colormap[44].c[3] = 0.0;
_colormap[45].c[0] = 0.741667; _colormap[45].c[1] = 0.529167; _colormap[45].c[2] = 0.908333; _colormap[45].c[3] = 0.0;
_colormap[46].c[0] = 0.783333; _colormap[46].c[1] = 0.508333; _colormap[46].c[2] = 0.916667; _colormap[46].c[3] = 0.0;
_colormap[47].c[0] = 0.825; _colormap[47].c[1] = 0.4875; _colormap[47].c[2] = 0.925; _colormap[47].c[3] = 0.0;
_colormap[48].c[0] = 0.866667; _colormap[48].c[1] = 0.466667; _colormap[48].c[2] = 0.933333; _colormap[48].c[3] = 0.0;
_colormap[49].c[0] = 0.875; _colormap[49].c[1] = 0.499020; _colormap[49].c[2] = 0.8875; _colormap[49].c[3] = 0.0;
_colormap[50].c[0] = 0.883333; _colormap[50].c[1] = 0.531373; _colormap[50].c[2] = 0.841667; _colormap[50].c[3] = 0.0;
_colormap[51].c[0] = 0.891667; _colormap[51].c[1] = 0.563725; _colormap[51].c[2] = 0.795833; _colormap[51].c[3] = 0.0;
_colormap[52].c[0] = 0.9; _colormap[52].c[1] = 0.596078; _colormap[52].c[2] = 0.75; _colormap[52].c[3] = 0.0;
_colormap[53].c[0] = 0.908333; _colormap[53].c[1] = 0.628431; _colormap[53].c[2] = 0.704167; _colormap[53].c[3] = 0.0;
_colormap[54].c[0] = 0.916667; _colormap[54].c[1] = 0.660784; _colormap[54].c[2] = 0.658333; _colormap[54].c[3] = 0.0;
_colormap[55].c[0] = 0.925; _colormap[55].c[1] = 0.693137; _colormap[55].c[2] = 0.6125; _colormap[55].c[3] = 0.0;
_colormap[56].c[0] = 0.933333; _colormap[56].c[1] = 0.725490; _colormap[56].c[2] = 0.566667; _colormap[56].c[3] = 0.0;
_colormap[57].c[0] = 0.941667; _colormap[57].c[1] = 0.757843; _colormap[57].c[2] = 0.520833; _colormap[57].c[3] = 0.0;
_colormap[58].c[0] = 0.95; _colormap[58].c[1] = 0.790196; _colormap[58].c[2] = 0.475; _colormap[58].c[3] = 0.0;
_colormap[59].c[0] = 0.958333; _colormap[59].c[1] = 0.822549; _colormap[59].c[2] = 0.429167; _colormap[59].c[3] = 0.0;
_colormap[60].c[0] = 0.966667; _colormap[60].c[1] = 0.854902; _colormap[60].c[2] = 0.383333; _colormap[60].c[3] = 0.0;
_colormap[61].c[0] = 0.975; _colormap[61].c[1] = 0.887255; _colormap[61].c[2] = 0.3375; _colormap[61].c[3] = 0.0;
_colormap[62].c[0] = 0.983333; _colormap[62].c[1] = 0.919608; _colormap[62].c[2] = 0.291667; _colormap[62].c[3] = 0.0;
_colormap[63].c[0] = 0.991667; _colormap[63].c[1] = 0.951961; _colormap[63].c[2] = 0.245833; _colormap[63].c[3] = 0.0;
}

void init_color_map_spring()
{
  _colormap_len = 64;
  _colormap = (rgbacolor *)malloc(_colormap_len*sizeof(rgbacolor));
  _colormap[0].c[0] = 0.066667; _colormap[0].c[1] = 0.4; _colormap[0].c[2] = 0.533333; _colormap[0].c[3] = 0.0;
_colormap[1].c[0] = 0.096774; _colormap[1].c[1] = 0.415054; _colormap[1].c[2] = 0.522581; _colormap[1].c[3] = 0.0;
_colormap[2].c[0] = 0.126882; _colormap[2].c[1] = 0.430108; _colormap[2].c[2] = 0.511828; _colormap[2].c[3] = 0.0;
_colormap[3].c[0] = 0.156989; _colormap[3].c[1] = 0.445161; _colormap[3].c[2] = 0.501075; _colormap[3].c[3] = 0.0;
_colormap[4].c[0] = 0.187097; _colormap[4].c[1] = 0.460215; _colormap[4].c[2] = 0.490323; _colormap[4].c[3] = 0.0;
_colormap[5].c[0] = 0.217204; _colormap[5].c[1] = 0.475269; _colormap[5].c[2] = 0.479570; _colormap[5].c[3] = 0.0;
_colormap[6].c[0] = 0.247312; _colormap[6].c[1] = 0.490323; _colormap[6].c[2] = 0.468817; _colormap[6].c[3] = 0.0;
_colormap[7].c[0] = 0.277419; _colormap[7].c[1] = 0.505376; _colormap[7].c[2] = 0.458065; _colormap[7].c[3] = 0.0;
_colormap[8].c[0] = 0.307527; _colormap[8].c[1] = 0.520430; _colormap[8].c[2] = 0.447312; _colormap[8].c[3] = 0.0;
_colormap[9].c[0] = 0.337634; _colormap[9].c[1] = 0.535484; _colormap[9].c[2] = 0.436559; _colormap[9].c[3] = 0.0;
_colormap[10].c[0] = 0.367742; _colormap[10].c[1] = 0.550538; _colormap[10].c[2] = 0.425806; _colormap[10].c[3] = 0.0;
_colormap[11].c[0] = 0.397849; _colormap[11].c[1] = 0.565591; _colormap[11].c[2] = 0.415054; _colormap[11].c[3] = 0.0;
_colormap[12].c[0] = 0.427957; _colormap[12].c[1] = 0.580645; _colormap[12].c[2] = 0.404301; _colormap[12].c[3] = 0.0;
_colormap[13].c[0] = 0.458064; _colormap[13].c[1] = 0.595699; _colormap[13].c[2] = 0.393548; _colormap[13].c[3] = 0.0;
_colormap[14].c[0] = 0.488172; _colormap[14].c[1] = 0.610753; _colormap[14].c[2] = 0.382796; _colormap[14].c[3] = 0.0;
_colormap[15].c[0] = 0.518280; _colormap[15].c[1] = 0.625806; _colormap[15].c[2] = 0.372043; _colormap[15].c[3] = 0.0;
_colormap[16].c[0] = 0.548387; _colormap[16].c[1] = 0.640860; _colormap[16].c[2] = 0.361290; _colormap[16].c[3] = 0.0;
_colormap[17].c[0] = 0.578495; _colormap[17].c[1] = 0.655914; _colormap[17].c[2] = 0.350538; _colormap[17].c[3] = 0.0;
_colormap[18].c[0] = 0.608602; _colormap[18].c[1] = 0.670968; _colormap[18].c[2] = 0.339785; _colormap[18].c[3] = 0.0;
_colormap[19].c[0] = 0.638710; _colormap[19].c[1] = 0.686022; _colormap[19].c[2] = 0.329032; _colormap[19].c[3] = 0.0;
_colormap[20].c[0] = 0.668817; _colormap[20].c[1] = 0.701075; _colormap[20].c[2] = 0.318280; _colormap[20].c[3] = 0.0;
_colormap[21].c[0] = 0.698925; _colormap[21].c[1] = 0.716129; _colormap[21].c[2] = 0.307527; _colormap[21].c[3] = 0.0;
_colormap[22].c[0] = 0.729032; _colormap[22].c[1] = 0.731183; _colormap[22].c[2] = 0.296774; _colormap[22].c[3] = 0.0;
_colormap[23].c[0] = 0.759140; _colormap[23].c[1] = 0.746237; _colormap[23].c[2] = 0.286022; _colormap[23].c[3] = 0.0;
_colormap[24].c[0] = 0.789247; _colormap[24].c[1] = 0.761290; _colormap[24].c[2] = 0.275269; _colormap[24].c[3] = 0.0;
_colormap[25].c[0] = 0.819355; _colormap[25].c[1] = 0.776344; _colormap[25].c[2] = 0.264516; _colormap[25].c[3] = 0.0;
_colormap[26].c[0] = 0.849462; _colormap[26].c[1] = 0.791398; _colormap[26].c[2] = 0.253763; _colormap[26].c[3] = 0.0;
_colormap[27].c[0] = 0.879570; _colormap[27].c[1] = 0.806452; _colormap[27].c[2] = 0.243011; _colormap[27].c[3] = 0.0;
_colormap[28].c[0] = 0.909677; _colormap[28].c[1] = 0.821505; _colormap[28].c[2] = 0.232258; _colormap[28].c[3] = 0.0;
_colormap[29].c[0] = 0.939785; _colormap[29].c[1] = 0.836559; _colormap[29].c[2] = 0.221505; _colormap[29].c[3] = 0.0;
_colormap[30].c[0] = 0.969892; _colormap[30].c[1] = 0.851613; _colormap[30].c[2] = 0.210753; _colormap[30].c[3] = 0.0;
_colormap[31].c[0] = 1.; _colormap[31].c[1] = 0.866667; _colormap[31].c[2] = 0.2; _colormap[31].c[3] = 0.0;
_colormap[32].c[0] = 0.995833; _colormap[32].c[1] = 0.843750; _colormap[32].c[2] = 0.214583; _colormap[32].c[3] = 0.0;
_colormap[33].c[0] = 0.991667; _colormap[33].c[1] = 0.820833; _colormap[33].c[2] = 0.229167; _colormap[33].c[3] = 0.0;
_colormap[34].c[0] = 0.9875; _colormap[34].c[1] = 0.797917; _colormap[34].c[2] = 0.243750; _colormap[34].c[3] = 0.0;
_colormap[35].c[0] = 0.983333; _colormap[35].c[1] = 0.775; _colormap[35].c[2] = 0.258333; _colormap[35].c[3] = 0.0;
_colormap[36].c[0] = 0.979167; _colormap[36].c[1] = 0.752083; _colormap[36].c[2] = 0.272917; _colormap[36].c[3] = 0.0;
_colormap[37].c[0] = 0.975; _colormap[37].c[1] = 0.729167; _colormap[37].c[2] = 0.2875; _colormap[37].c[3] = 0.0;
_colormap[38].c[0] = 0.970833; _colormap[38].c[1] = 0.706250; _colormap[38].c[2] = 0.302083; _colormap[38].c[3] = 0.0;
_colormap[39].c[0] = 0.966667; _colormap[39].c[1] = 0.683333; _colormap[39].c[2] = 0.316667; _colormap[39].c[3] = 0.0;
_colormap[40].c[0] = 0.9625; _colormap[40].c[1] = 0.660417; _colormap[40].c[2] = 0.331250; _colormap[40].c[3] = 0.0;
_colormap[41].c[0] = 0.958333; _colormap[41].c[1] = 0.6375; _colormap[41].c[2] = 0.345833; _colormap[41].c[3] = 0.0;
_colormap[42].c[0] = 0.954167; _colormap[42].c[1] = 0.614583; _colormap[42].c[2] = 0.360417; _colormap[42].c[3] = 0.0;
_colormap[43].c[0] = 0.95; _colormap[43].c[1] = 0.591667; _colormap[43].c[2] = 0.375; _colormap[43].c[3] = 0.0;
_colormap[44].c[0] = 0.945833; _colormap[44].c[1] = 0.568750; _colormap[44].c[2] = 0.389583; _colormap[44].c[3] = 0.0;
_colormap[45].c[0] = 0.941667; _colormap[45].c[1] = 0.545833; _colormap[45].c[2] = 0.404167; _colormap[45].c[3] = 0.0;
_colormap[46].c[0] = 0.9375; _colormap[46].c[1] = 0.522917; _colormap[46].c[2] = 0.418750; _colormap[46].c[3] = 0.0;
_colormap[47].c[0] = 0.933333; _colormap[47].c[1] = 0.5; _colormap[47].c[2] = 0.433333; _colormap[47].c[3] = 0.0;
_colormap[48].c[0] = 0.929167; _colormap[48].c[1] = 0.477083; _colormap[48].c[2] = 0.447917; _colormap[48].c[3] = 0.0;
_colormap[49].c[0] = 0.925; _colormap[49].c[1] = 0.454167; _colormap[49].c[2] = 0.4625; _colormap[49].c[3] = 0.0;
_colormap[50].c[0] = 0.920833; _colormap[50].c[1] = 0.431250; _colormap[50].c[2] = 0.477083; _colormap[50].c[3] = 0.0;
_colormap[51].c[0] = 0.916667; _colormap[51].c[1] = 0.408333; _colormap[51].c[2] = 0.491667; _colormap[51].c[3] = 0.0;
_colormap[52].c[0] = 0.9125; _colormap[52].c[1] = 0.385417; _colormap[52].c[2] = 0.506250; _colormap[52].c[3] = 0.0;
_colormap[53].c[0] = 0.908333; _colormap[53].c[1] = 0.3625; _colormap[53].c[2] = 0.520833; _colormap[53].c[3] = 0.0;
_colormap[54].c[0] = 0.904167; _colormap[54].c[1] = 0.339583; _colormap[54].c[2] = 0.535417; _colormap[54].c[3] = 0.0;
_colormap[55].c[0] = 0.9; _colormap[55].c[1] = 0.316667; _colormap[55].c[2] = 0.55; _colormap[55].c[3] = 0.0;
_colormap[56].c[0] = 0.895833; _colormap[56].c[1] = 0.293750; _colormap[56].c[2] = 0.564583; _colormap[56].c[3] = 0.0;
_colormap[57].c[0] = 0.891667; _colormap[57].c[1] = 0.270833; _colormap[57].c[2] = 0.579167; _colormap[57].c[3] = 0.0;
_colormap[58].c[0] = 0.8875; _colormap[58].c[1] = 0.247917; _colormap[58].c[2] = 0.593750; _colormap[58].c[3] = 0.0;
_colormap[59].c[0] = 0.883333; _colormap[59].c[1] = 0.225; _colormap[59].c[2] = 0.608333; _colormap[59].c[3] = 0.0;
_colormap[60].c[0] = 0.879167; _colormap[60].c[1] = 0.202083; _colormap[60].c[2] = 0.622917; _colormap[60].c[3] = 0.0;
_colormap[61].c[0] = 0.875; _colormap[61].c[1] = 0.179167; _colormap[61].c[2] = 0.6375; _colormap[61].c[3] = 0.0;
_colormap[62].c[0] = 0.870833; _colormap[62].c[1] = 0.156250; _colormap[62].c[2] = 0.652083; _colormap[62].c[3] = 0.0;
_colormap[63].c[0] = 0.866667; _colormap[63].c[1] = 0.133333; _colormap[63].c[2] = 0.666667; _colormap[63].c[3] = 0.0;

}


void init_color_map_jet()
{
  _colormap_len = 64;
  _colormap = (rgbacolor *)malloc(_colormap_len*sizeof(rgbacolor));
  _colormap[0].c[0] = 0.; _colormap[0].c[1] = 0.; _colormap[0].c[2] = 0.5625; _colormap[0].c[3] = 0.0;
  _colormap[1].c[0] = 0.; _colormap[1].c[1] = 0.; _colormap[1].c[2] = 0.625; _colormap[1].c[3] = 0.0;
  _colormap[2].c[0] = 0.; _colormap[2].c[1] = 0.; _colormap[2].c[2] = 0.6875; _colormap[2].c[3] = 0.0;
  _colormap[3].c[0] = 0.; _colormap[3].c[1] = 0.; _colormap[3].c[2] = 0.75; _colormap[3].c[3] = 0.0;
  _colormap[4].c[0] = 0.; _colormap[4].c[1] = 0.; _colormap[4].c[2] = 0.8125; _colormap[4].c[3] = 0.0;
  _colormap[5].c[0] = 0.; _colormap[5].c[1] = 0.; _colormap[5].c[2] = 0.875; _colormap[5].c[3] = 0.0;
  _colormap[6].c[0] = 0.; _colormap[6].c[1] = 0.; _colormap[6].c[2] = 0.9375; _colormap[6].c[3] = 0.0;
  _colormap[7].c[0] = 0.; _colormap[7].c[1] = 0.; _colormap[7].c[2] = 1.; _colormap[7].c[3] = 0.0;
  _colormap[8].c[0] = 0.; _colormap[8].c[1] = 0.0625; _colormap[8].c[2] = 1.; _colormap[8].c[3] = 0.0;
  _colormap[9].c[0] = 0.; _colormap[9].c[1] = 0.125; _colormap[9].c[2] = 1.; _colormap[9].c[3] = 0.0;
  _colormap[10].c[0] = 0.; _colormap[10].c[1] = 0.1875; _colormap[10].c[2] = 1.; _colormap[10].c[3] = 0.0;
  _colormap[11].c[0] = 0.; _colormap[11].c[1] = 0.25; _colormap[11].c[2] = 1.; _colormap[11].c[3] = 0.0;
  _colormap[12].c[0] = 0.; _colormap[12].c[1] = 0.3125; _colormap[12].c[2] = 1.; _colormap[12].c[3] = 0.0;
  _colormap[13].c[0] = 0.; _colormap[13].c[1] = 0.375; _colormap[13].c[2] = 1.; _colormap[13].c[3] = 0.0;
  _colormap[14].c[0] = 0.; _colormap[14].c[1] = 0.4375; _colormap[14].c[2] = 1.; _colormap[14].c[3] = 0.0;
  _colormap[15].c[0] = 0.; _colormap[15].c[1] = 0.5; _colormap[15].c[2] = 1.; _colormap[15].c[3] = 0.0;
  _colormap[16].c[0] = 0.; _colormap[16].c[1] = 0.5625; _colormap[16].c[2] = 1.; _colormap[16].c[3] = 0.0;
  _colormap[17].c[0] = 0.; _colormap[17].c[1] = 0.625; _colormap[17].c[2] = 1.; _colormap[17].c[3] = 0.0;
  _colormap[18].c[0] = 0.; _colormap[18].c[1] = 0.6875; _colormap[18].c[2] = 1.; _colormap[18].c[3] = 0.0;
  _colormap[19].c[0] = 0.; _colormap[19].c[1] = 0.75; _colormap[19].c[2] = 1.; _colormap[19].c[3] = 0.0;
  _colormap[20].c[0] = 0.; _colormap[20].c[1] = 0.8125; _colormap[20].c[2] = 1.; _colormap[20].c[3] = 0.0;
  _colormap[21].c[0] = 0.; _colormap[21].c[1] = 0.875; _colormap[21].c[2] = 1.; _colormap[21].c[3] = 0.0;
  _colormap[22].c[0] = 0.; _colormap[22].c[1] = 0.9375; _colormap[22].c[2] = 1.; _colormap[22].c[3] = 0.0;
  _colormap[23].c[0] = 0.; _colormap[23].c[1] = 1.; _colormap[23].c[2] = 1.; _colormap[23].c[3] = 0.0;
  _colormap[24].c[0] = 0.0625; _colormap[24].c[1] = 1.; _colormap[24].c[2] = 0.9375; _colormap[24].c[3] = 0.0;
  _colormap[25].c[0] = 0.125; _colormap[25].c[1] = 1.; _colormap[25].c[2] = 0.875; _colormap[25].c[3] = 0.0;
  _colormap[26].c[0] = 0.1875; _colormap[26].c[1] = 1.; _colormap[26].c[2] = 0.8125; _colormap[26].c[3] = 0.0;
  _colormap[27].c[0] = 0.25; _colormap[27].c[1] = 1.; _colormap[27].c[2] = 0.75; _colormap[27].c[3] = 0.0;
  _colormap[28].c[0] = 0.3125; _colormap[28].c[1] = 1.; _colormap[28].c[2] = 0.6875; _colormap[28].c[3] = 0.0;
  _colormap[29].c[0] = 0.375; _colormap[29].c[1] = 1.; _colormap[29].c[2] = 0.625; _colormap[29].c[3] = 0.0;
  _colormap[30].c[0] = 0.4375; _colormap[30].c[1] = 1.; _colormap[30].c[2] = 0.5625; _colormap[30].c[3] = 0.0;
  _colormap[31].c[0] = 0.5; _colormap[31].c[1] = 1.; _colormap[31].c[2] = 0.5; _colormap[31].c[3] = 0.0;
  _colormap[32].c[0] = 0.5625; _colormap[32].c[1] = 1.; _colormap[32].c[2] = 0.4375; _colormap[32].c[3] = 0.0;
  _colormap[33].c[0] = 0.625; _colormap[33].c[1] = 1.; _colormap[33].c[2] = 0.375; _colormap[33].c[3] = 0.0;
  _colormap[34].c[0] = 0.6875; _colormap[34].c[1] = 1.; _colormap[34].c[2] = 0.3125; _colormap[34].c[3] = 0.0;
  _colormap[35].c[0] = 0.75; _colormap[35].c[1] = 1.; _colormap[35].c[2] = 0.25; _colormap[35].c[3] = 0.0;
  _colormap[36].c[0] = 0.8125; _colormap[36].c[1] = 1.; _colormap[36].c[2] = 0.1875; _colormap[36].c[3] = 0.0;
  _colormap[37].c[0] = 0.875; _colormap[37].c[1] = 1.; _colormap[37].c[2] = 0.125; _colormap[37].c[3] = 0.0;
  _colormap[38].c[0] = 0.9375; _colormap[38].c[1] = 1.; _colormap[38].c[2] = 0.0625; _colormap[38].c[3] = 0.0;
  _colormap[39].c[0] = 1.; _colormap[39].c[1] = 1.; _colormap[39].c[2] = 0.; _colormap[39].c[3] = 0.0;
  _colormap[40].c[0] = 1.; _colormap[40].c[1] = 0.9375; _colormap[40].c[2] = 0.; _colormap[40].c[3] = 0.0;
  _colormap[41].c[0] = 1.; _colormap[41].c[1] = 0.875; _colormap[41].c[2] = 0.; _colormap[41].c[3] = 0.0;
  _colormap[42].c[0] = 1.; _colormap[42].c[1] = 0.8125; _colormap[42].c[2] = 0.; _colormap[42].c[3] = 0.0;
  _colormap[43].c[0] = 1.; _colormap[43].c[1] = 0.75; _colormap[43].c[2] = 0.; _colormap[43].c[3] = 0.0;
  _colormap[44].c[0] = 1.; _colormap[44].c[1] = 0.6875; _colormap[44].c[2] = 0.; _colormap[44].c[3] = 0.0;
  _colormap[45].c[0] = 1.; _colormap[45].c[1] = 0.625; _colormap[45].c[2] = 0.; _colormap[45].c[3] = 0.0;
  _colormap[46].c[0] = 1.; _colormap[46].c[1] = 0.5625; _colormap[46].c[2] = 0.; _colormap[46].c[3] = 0.0;
  _colormap[47].c[0] = 1.; _colormap[47].c[1] = 0.5; _colormap[47].c[2] = 0.; _colormap[47].c[3] = 0.0;
  _colormap[48].c[0] = 1.; _colormap[48].c[1] = 0.4375; _colormap[48].c[2] = 0.; _colormap[48].c[3] = 0.0;
  _colormap[49].c[0] = 1.; _colormap[49].c[1] = 0.375; _colormap[49].c[2] = 0.; _colormap[49].c[3] = 0.0;
  _colormap[50].c[0] = 1.; _colormap[50].c[1] = 0.3125; _colormap[50].c[2] = 0.; _colormap[50].c[3] = 0.0;
  _colormap[51].c[0] = 1.; _colormap[51].c[1] = 0.25; _colormap[51].c[2] = 0.; _colormap[51].c[3] = 0.0;
  _colormap[52].c[0] = 1.; _colormap[52].c[1] = 0.1875; _colormap[52].c[2] = 0.; _colormap[52].c[3] = 0.0;
  _colormap[53].c[0] = 1.; _colormap[53].c[1] = 0.125; _colormap[53].c[2] = 0.; _colormap[53].c[3] = 0.0;
  _colormap[54].c[0] = 1.; _colormap[54].c[1] = 0.0625; _colormap[54].c[2] = 0.; _colormap[54].c[3] = 0.0;
  _colormap[55].c[0] = 1.; _colormap[55].c[1] = 0.; _colormap[55].c[2] = 0.; _colormap[55].c[3] = 0.0;
  _colormap[56].c[0] = 0.9375; _colormap[56].c[1] = 0.; _colormap[56].c[2] = 0.; _colormap[56].c[3] = 0.0;
  _colormap[57].c[0] = 0.875; _colormap[57].c[1] = 0.; _colormap[57].c[2] = 0.; _colormap[57].c[3] = 0.0;
  _colormap[58].c[0] = 0.8125; _colormap[58].c[1] = 0.; _colormap[58].c[2] = 0.; _colormap[58].c[3] = 0.0;
  _colormap[59].c[0] = 0.75; _colormap[59].c[1] = 0.; _colormap[59].c[2] = 0.; _colormap[59].c[3] = 0.0;
  _colormap[60].c[0] = 0.6875; _colormap[60].c[1] = 0.; _colormap[60].c[2] = 0.; _colormap[60].c[3] = 0.0;
  _colormap[61].c[0] = 0.625; _colormap[61].c[1] = 0.; _colormap[61].c[2] = 0.; _colormap[61].c[3] = 0.0;
  _colormap[62].c[0] = 0.5625; _colormap[62].c[1] = 0.; _colormap[62].c[2] = 0.; _colormap[62].c[3] = 0.0;
  _colormap[63].c[0] = 0.5; _colormap[63].c[1] = 0.; _colormap[63].c[2] = 0.; _colormap[63].c[3] = 0.0;
}

void init_color_map_custom(char *filename)
{
  int retval, i=0;
  int size = 4;
  char c;
  if ( ( _colormap_file = fopen(filename,"r") ) != NULL )
  {
    _colormap_name = NB_COLORMAPS; /* custom colormap */
    fprintf(_log_file,"loading colormap from %s\n", filename);
  }
  else
  {
    printf("Unknown colormap or file '%s', using colormap 'jet'.\n",filename);
    return;
  }
  
  _colormap = (rgbacolor *)malloc(size*sizeof(rgbacolor)); 
  if ( (retval = fscanf(_colormap_file,"#%c",&c)) == 1) 
  {
    fprintf(_log_file,"generated from colors ");
    do {retval = fscanf(_colormap_file, "%c", &c); fprintf(_log_file,"%c",c); }  while ((c != '\n') && (retval != EOF));
  }
  while ( (retval = fscanf(_colormap_file, "%f %f %f \n", 
          _colormap[i].c,_colormap[i].c+1,_colormap[i].c+2)) == 3 )
  {
    i++;
    if ( i >= size )
    {
      size *= 2;
      _colormap = (rgbacolor *)realloc(_colormap,size*sizeof(rgbacolor));
    }
  }
  fclose(_colormap_file);

  _colormap_len = i;
  _colormap = (rgbacolor *)realloc(_colormap,_colormap_len*sizeof(rgbacolor));
  fprintf(_log_file,"%d lines read.\n",i);


}

void init_opencv_video(const char* filename, int width, int height) {
    _video_writer.open(filename, _video_codec, _video_fps, cv::Size(width, height));
    if (!_video_writer.isOpened()) {
        fprintf(stderr, "Error: Could not open video writer for %s\n", filename);
        _record_video = 0;
    } else {
        printf("Video writer initialized: %dx%d @ %d fps -> %s\n", width, height, _video_fps, filename);
    }
}

void add_frame_to_video() {
    if (!_record_video) return;
    
    // Get current window dimensions
    int width = glutGet(GLUT_WINDOW_WIDTH);
    int height = glutGet(GLUT_WINDOW_HEIGHT);
    
    // Initialize video writer on first frame
    if (!_video_writer.isOpened()) {
        char full_filename[512];
        const char* ext = strrchr(_video_filename, '.');
        if (!ext || strcmp(ext, ".mp4") != 0) {
            snprintf(full_filename, 512, "%s.mp4", _video_filename);
        } else {
            strncpy(full_filename, _video_filename, 512);
        }
        init_opencv_video(full_filename, width, height);
        if (!_video_writer.isOpened()) return;
    }
    
    // Read pixels from OpenGL
    unsigned char* buffer = (unsigned char*)malloc(width * height * 3);
    glPixelStorei(GL_PACK_ALIGNMENT, 1);
    glReadBuffer(GL_FRONT);
    glReadPixels(0, 0, width, height, GL_RGB, GL_UNSIGNED_BYTE, buffer);
    
    // OpenCV Mat and flip vertically
    cv::Mat frame(height, width, CV_8UC3, buffer);
    cv::flip(frame, frame, 0); // Flip vertically
    cv::cvtColor(frame, frame, cv::COLOR_RGB2BGR); 
    
    // Write frame to video
    _video_writer.write(frame);
    
    free(buffer);
    _video_frame_count++;
    
    if (_video_frame_count % 30 == 0) {
        printf("Recorded frame %d\n", _video_frame_count);
    }
}

void cleanup_opencv_video() {
    if (_video_writer.isOpened()) {
        _video_writer.release();
        printf("Video recording completed: %d frames saved\n", _video_frame_count);
    }
}