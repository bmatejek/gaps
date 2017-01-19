// Source file for neuron visualizer



// include files

#include "RNDataStructures/RNDataStructures.h"
#include "R3Graphics/R3Graphics.h"
#include "fglut/fglut.h"
#include <vector>
#include <fstream>
#include <string>



// GLUT defines

#define ENTER 13
#define ESCAPE 27
#define SPACEBAR 32



// program arguments

static int print_debug = 0;
static int print_verbose = 0;
static const char *input_filename = NULL;
static const char *input_dataset = NULL;



// program variables

static RNScalar scaling[3] = {1, 1, 5};
static int resolution[3] = {-1, -1, -1};
static R3Affine transformation = R3null_affine;
static R3Viewer *viewer = NULL;
static R3Box world_box;



// voxel grids

static R3CharGrid *grid = NULL;



// projection variables

static R2Grid *selected_slice = NULL;
static R2Point selected_slice_position(RN_UNKNOWN, RN_UNKNOWN);
static RNInterval selected_slice_range(0, 0);
static RNScalar projection_scale = FLT_MAX;
static R3Mesh *mesh = NULL;


// GLUT variables

static int GLUTwindow = 0;
static int GLUTwindow_height = 1200;
static int GLUTwindow_width = 1200;
static int GLUTmouse[2] = { 0, 0 };
static int GLUTbutton[3] = { 0, 0, 0 };
static int GLUTmodifiers = 0;



// color arrays

static RNScalar background_color[] = { 0, 0, 0 };



////////////////////////////////////////////////////////////////////////
// Input/output functions
////////////////////////////////////////////////////////////////////////

static int ReadData(void)
{
   // start statistics
   RNTime start_time;
   start_time.Read();

   // read in voxel file
   grid = RNReadH5ImageFile(input_filename, input_dataset);
   if (!grid) { fprintf(stderr, "Failed to read %s from %s\n", input_dataset, input_filename); return 0; }
   
    resolution[RN_X] = grid->XResolution();
    resolution[RN_Y] = grid->YResolution();
    resolution[RN_Z] = grid->ZResolution();

   // print statistics
   if (print_verbose) {
      printf("Read voxel grids...\n");
      printf("  Time = %.2f seconds\n", start_time.Elapsed());
      printf("  Resolution = (%lu %lu %lu)\n", grid->XResolution(), grid->YResolution(), grid->ZResolution());
      printf("  Scaling = (%lf %lf %lf)\n", scaling[RN_X], scaling[RN_Y], scaling[RN_Z]);
   }

   // return success
   return 1;
}



////////////////////////////////////////////////////////////////////////
// Drawing utility functions
////////////////////////////////////////////////////////////////////////

static void DrawSlice(void)
{ 
   transformation.Push();
   printf("A\n");
   RNLoadRgb(0.0, 0.0, 0.0);
   grid->DrawSlice(RN_Z, 0);
   printf("D\n");
   transformation.Pop();
}




static void Draw3D(void)
{
   // set lights
   static GLfloat light0_position[] = { 3.0, 4.0, 5.0, 0.0 };
   glLightfv(GL_LIGHT0, GL_POSITION, light0_position);
   static GLfloat light1_position[] = { -3.0, -2.0, -3.0, 0.0 };
   glLightfv(GL_LIGHT1, GL_POSITION, light1_position);
    
    int index = GL_LIGHT1 + 1;
    
    R3PointLight point_light = R3PointLight(viewer->Camera().Origin(), RNwhite_rgb);
    mesh->Draw();
    point_light.Draw(index);
    

   // prologue
   //glDisable(GL_LIGHTING);

    //RNLoadRgb(RNred_rgb);
    //mesh->Draw();

   DrawSlice();
   
   //glEnable(GL_LIGHTING);
}



////////////////////////////////////////////////////////////////////////
// GLUT interface functions
////////////////////////////////////////////////////////////////////////

void GLUTStop(void)
{
   // destroy window
   glutDestroyWindow(GLUTwindow);

   // delete the neuron data
   RNTime start_time;
   start_time.Read();

   //delete[] machine_label;

   // print statistics
   if (print_verbose) {
      printf("Deleted data...\n");
      printf("  Time = %.2f seconds\n", start_time.Elapsed());
      fflush(stdout);
   }

   // exit
   exit(0);
}



void GLUTRedraw(void)
{
    // set viewing transformation
    viewer->Camera().Load();

   // clear window
   glClearColor(background_color[0], background_color[1], background_color[2], 1.0);
   glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

   // drawing varies on projection
   Draw3D();

   // swap buffers
   glutSwapBuffers();
}



void GLUTResize(int w, int h)
{
   // resize window
   glViewport(0, 0, w, h);

   // resize viewer viewport
   viewer->ResizeViewport(0, 0, w, h);

   // remember window size
   GLUTwindow_width = w;
   GLUTwindow_height = h;

   // redraw
   glutPostRedisplay();
}



void GLUTMotion3D(int x, int y)
{
   // compute mouse movement
   int dx = x - GLUTmouse[0];
   int dy = y - GLUTmouse[1];

   // world in hand navigation
   R3Point origin = world_box.Centroid();
   if (GLUTbutton[0]) viewer->RotateWorld(1.0, origin, x, y, dx, dy);
   else if (GLUTbutton[1]) viewer->ScaleWorld(1.0, origin, x, y, dx, dy);
}



void GLUTMotion(int x, int y)
{
   // invert y coordinate
   y = GLUTwindow_height - y;

   // different motion clicks for projection view
   GLUTMotion3D(x, y);

   // redisplay if a mouse was down
   if (GLUTbutton[0] || GLUTbutton[1] || GLUTbutton[2]) glutPostRedisplay();

   // remember mouse position
   GLUTmouse[0] = x;
   GLUTmouse[1] = y;
}



void GLUTMouse(int button, int state, int x, int y)
{
   // invert y coordinate
   y = GLUTwindow_height - y;

   // remember button state
   int b = (button == GLUT_LEFT_BUTTON) ? 0 : ((button == GLUT_MIDDLE_BUTTON) ? 1 : 2);
   GLUTbutton[b] = (state == GLUT_DOWN) ? 1 : 0;

   // remember modifiers
   GLUTmodifiers = glutGetModifiers();

   // remember mouse position
   GLUTmouse[0] = x;
   GLUTmouse[1] = y;

   // redraw
   glutPostRedisplay();
}



void GLUTInit(int *argc, char **argv)
{
   // open window
   glutInit(argc, argv);
   glutInitWindowPosition(10, 10);
   glutInitWindowSize(GLUTwindow_width, GLUTwindow_height);
   glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH);
   GLUTwindow = glutCreateWindow("Neuron Visualizer");

   // initialize background color
   glClearColor(background_color[0], background_color[1], background_color[2], 1.0);

   // initialize lights
   static GLfloat lmodel_ambient[] = { 0.2, 0.2, 0.2, 1.0 };
   glLightModelfv(GL_LIGHT_MODEL_AMBIENT, lmodel_ambient);
   glLightModeli(GL_LIGHT_MODEL_LOCAL_VIEWER, GL_TRUE);
   static GLfloat light0_diffuse[] = { 1.0, 1.0, 1.0, 1.0 };
   glLightfv(GL_LIGHT0, GL_DIFFUSE, light0_diffuse);
   glEnable(GL_LIGHT0);
   static GLfloat light1_diffuse[] = { 0.5, 0.5, 0.5, 1.0 };
   glLightfv(GL_LIGHT1, GL_DIFFUSE, light1_diffuse);
   glEnable(GL_LIGHT1);
   glEnable(GL_NORMALIZE);
   glEnable(GL_LIGHTING);

   // initialize graphics mode
   glEnable(GL_DEPTH_TEST);

   // initialize GLUT callback functions
   glutDisplayFunc(GLUTRedraw);
   glutReshapeFunc(GLUTResize);
   glutMouseFunc(GLUTMouse);
   glutMotionFunc(GLUTMotion);

   // initialize font
#if (RN_OS == RN_WINDOWSNT)
   int font = glGenLists(256);
   wglUseFontBitmaps(wglGetCurrentDC(), 0, 256, font);
   glListBase(font);
#endif
}



void GLUTMainLoop(void)
{
   // run main loop -- never returns
   glutMainLoop();
}



////////////////////////////////////////////////////////////////////////
// Viewer functions
////////////////////////////////////////////////////////////////////////

static R3Viewer *CreateViewer(void)
{
   // start statistics
   RNTime start_time;
   start_time.Read();


   if (world_box.IsEmpty()) RNAbort("Error in CreateViewer - box is empty");
   RNLength r = world_box.DiagonalRadius();
   if (r < 0 || RNIsInfinite(r)) RNAbort("Error in CreateViewer - r must be positive finite");

   // set up camera view looking down the z axis
   static R3Vector initial_camera_towards = R3Vector(-1.0, 1.0, -1.0);
   static R3Vector initial_camera_up = R3Vector(-0.3, 0.3, 0.85);
   R3Point initial_camera_origin = world_box.Centroid() - initial_camera_towards * 2.5 * r;
   R3Camera camera(initial_camera_origin, initial_camera_towards, initial_camera_up, 0.4, 0.4, 0.1 * r, 1000.0 * r);
   R2Viewport viewport(0, 0, GLUTwindow_width, GLUTwindow_height);
   R3Viewer *viewer = new R3Viewer(camera, viewport);
   

   // print statistics
   if (print_verbose) {
      printf("Created viewer ...\n");
      printf("  Time = %.2f seconds\n", start_time.Elapsed());
      printf("  Origin = %g %g %g\n", camera.Origin().X(), camera.Origin().Y(), camera.Origin().Z());
      printf("  Towards = %g %g %g\n", camera.Towards().X(), camera.Towards().Y(), camera.Towards().Z());
      printf("  Up = %g %g %g\n", camera.Up().X(), camera.Up().Y(), camera.Up().Z());
      printf("  Fov = %g %g\n", camera.XFOV(), camera.YFOV());
      printf("  Near = %g\n", camera.Near());
      printf("  Far = %g\n", camera.Far());
      fflush(stdout);
   }

   // return viewer
   return viewer;
}



////////////////////////////////////////////////////////////////////////
// Program argument parsing
////////////////////////////////////////////////////////////////////////

static int ParseArgs(int argc, char **argv)
{
   // parse arguments
   argc--; argv++;
   while (argc > 0) {
      if ((*argv)[0] == '-') {
         if (!strcmp(*argv, "-v")) print_verbose = 1;
         else if (!strcmp(*argv, "-debug")) print_debug = 1;
         else if (!strcmp(*argv, "-scaling")) {
             argv++; argc--;
             scaling[RN_X] = atof(*argv);
             argv++; argc--;
             scaling[RN_Y] = atof(*argv);
             argv++; argc--; 
             scaling[RN_Z] = atof(*argv); 
         }
         else { fprintf(stderr, "Invalid program argument: %s\n", *argv); return 0; }
      }
      else {
          if (!input_filename) { input_filename = *argv; }
          else if (!input_dataset) { input_dataset = *argv; }
          else { fprintf(stderr, "Invalid program argument: %s\n", *argv); return 0; }
      }
      argv++; argc--;
   }

   // error if there is no input name
   if (!input_filename || !input_dataset) { fprintf(stderr, "Need to supply a neuron data file\n"); return 0; }

   // return success
   return 1;
}



////////////////////////////////////////////////////////////////////////
// Main program
////////////////////////////////////////////////////////////////////////

int main(int argc, char **argv)
{
   // parse program arguments
   if (!ParseArgs(argc, argv)) exit(-1);

   /////////////////////////////////
   //// Read in the voxel files ////
   /////////////////////////////////

   if (!ReadData()) exit(-1);

    mesh = new R3Mesh();
    mesh->ReadFile("meshes/ascii.stl");

   // set world box
   world_box = R3Box(0, 0, 0, grid->XResolution() * scaling[RN_X], grid->YResolution() * scaling[RN_Y], grid->ZResolution() * scaling[RN_Z]);

   // get the transformation
   transformation = R3Affine(R4Matrix(scaling[RN_X], 0, 0, 0, 0, scaling[RN_Y], 0, 0, 0, 0, scaling[RN_Z], 0, 0, 0, 0, 1));

   // set projection scale variables
   for (int dim = 0; dim <= 2; ++dim) {
      RNScalar ratio = GLUTwindow_height / (scaling[dim] * grid->Resolution(dim));
      if (ratio < projection_scale) projection_scale = ratio;
   }
    
   // create viewer
   viewer = CreateViewer();
   if (!viewer) exit(-1);

   // initialize GLUT
   GLUTInit(&argc, argv);

   // run GLUT interface
   GLUTMainLoop();

   // return success
   return 1;
}