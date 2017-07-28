// Source file for ebro feature visualizer



// include files

#include "R3Graphics/R3Graphics.h"
#include "RNDataStructures/RNDataStructures.h"
#include "fglut/fglut.h"
#include <vector>



// GLUT defines

#define ENTER 13
#define ESCAPE 27
#define SPACEBAR 32



// program arguments

static int print_debug = 0;
static int print_verbose = 0;
static const char *prefix_one = NULL;
static const char *prefix_two = NULL;
static unsigned long threshold = 10000;
// window radius in nanometers
static RNScalar window_radius = 600;
static RNScalar resolution[3] = { 4, 4, 30 };
static R3Affine transformation;


// feature variables

static R3Grid *grid = NULL;
static int candidate_index = 0;



// transformation variables

static R3Viewer* viewer = NULL;
static R3Point selected_position;
static R3Box world_box;



// GLUT variables

static int GLUTwindow = 0;
static int GLUTwindow_height = 1200;
static int GLUTwindow_width = 1200;
static int GLUTmouse[2] = { 0, 0 };
static int GLUTbutton[3] = { 0, 0, 0 };
static int GLUTmodifiers = 0;



// style variables

static RNScalar background_color[3] = { 0, 0, 0 };



// save decisions

enum DECISION { NO, YES, UNDECIDED };
std::vector<enum DECISION> decisions = std::vector<enum DECISION>();




////////////////////////////////////////////////////////////////////////
// Input/output functions
////////////////////////////////////////////////////////////////////////

static int
ReadDecisionFile(void)
{
    char input_filename[4096];
    sprintf(input_filename, "gold/%s-%s-%lu-%dnm.gold", prefix_one, prefix_two, threshold, (int)(window_radius + 0.5));

    // open file
    FILE *fp = fopen(input_filename, "rb");
    if (!fp) return 0;

    int ndecisions;
    if (fread(&ndecisions, sizeof(int), 1, fp) != 1) return 0;

    for (int iv = 0; iv < ndecisions; ++iv) {
        int decision;
        if (fread(&decision, sizeof(int), 1, fp) != 1) return 0;
        decisions.push_back(DECISION(decision));
    }

    // close file
    fclose(fp);

    // return success
    return 1;
}



static int ReadFeature(int index)
{
    // free old grid memory
    if (grid) delete grid;

    // get the filename for this feature
    char feature_filename[4096];
    sprintf(feature_filename, "features/ebro/%s-%s/%lu-%dnm-%05d.h5", prefix_one, prefix_two, threshold, (int)(window_radius + 0.5), candidate_index);

    // read the feature into an R3Grid
    R3Grid **grids = RNReadH5File(feature_filename, "main");
    grid = grids[0];
    delete[] grids;

    // return success
    return 1;
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

    if (grid) delete grid;

    // print statistics
    if(print_verbose) {
        printf("Deleted data...\n");
        printf("  Time = %.2f seconds\n", start_time.Elapsed());
        fflush(stdout);
    }

    // exit
    exit(0);
}



void GLUTRedraw(void)
{
    // clear window
    glClearColor(background_color[0], background_color[1], background_color[2], 1.0);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    // set viewing transformation
    viewer->Camera().Load();

    // set lights
    static GLfloat light0_position[] = { 3.0, 4.0, 5.0, 0.0 };
    glLightfv(GL_LIGHT0, GL_POSITION, light0_position);
    static GLfloat light1_position[] = { -3.0, -2.0, -3.0, 0.0 };
    glLightfv(GL_LIGHT1, GL_POSITION, light1_position);

    // prologue
    glDisable(GL_LIGHTING);

    transformation.Push();

    // draw feature bounding box
    RNLoadRgb(RNwhite_rgb);
    world_box.Outline();

    // draw the actual points in 3D
    glBegin(GL_POINTS);
    for (int iz = 0; iz < grid->ZResolution(); ++iz) {
        for (int iy = 0; iy < grid->YResolution(); ++iy) {
            for (int ix = 0; ix < grid->XResolution(); ++ix) {
                int grid_value = (int)(grid->GridValue(ix, iy, iz) + 0.5);

                if (grid_value) {
                    if (RNRandomScalar() > 1.0 / 3.0) continue;
                }

                if (grid_value == 1) {
                    if (decisions[candidate_index] == YES) RNLoadRgb(RNgreen_rgb);
                    else RNLoadRgb(RNred_rgb);
                    glVertex3f(ix, iy, iz);
                }
                else if (grid_value == 2) {
                    if (decisions[candidate_index] == YES) RNLoadRgb(RNblue_rgb);
                    else RNLoadRgb(RNyellow_rgb);
                    glVertex3f(ix, iy, iz);
                }
                else if (grid_value == 3) {
                    if (decisions[candidate_index] == YES) RNLoadRgb(RNcyan_rgb);
                    else RNLoadRgb(RNmagenta_rgb);
                    glVertex3f(ix, iy,iz);
                }
            }
        }
    }
    glEnd();

    transformation.Pop();

    // write the feature
    char feature_label[4096];
    if (decisions[candidate_index] == YES) sprintf(feature_label, "Feature Visualizer - %d - YES - Predicted: \n", candidate_index);
    if (decisions[candidate_index] == NO) sprintf(feature_label, "Feature Visualizer - %d - NO - Predicted: \n", candidate_index);
    if (decisions[candidate_index] == UNDECIDED) sprintf(feature_label, "Feature Visualizer - %d - UNDECIDED - Predicted: \n", candidate_index);
    glutSetWindowTitle(feature_label);

    // epilogue
    glEnable(GL_LIGHTING);

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



void GLUTMotion(int x, int y)
{
    // invert y coordinate
    y = GLUTwindow_height - y;

    // remember mouse position
    int dx = x - GLUTmouse[0];
    int dy = y - GLUTmouse[1];

    // world in hand navigation
    R3Point origin = world_box.Centroid();
    if(GLUTbutton[0])
        viewer->RotateWorld(1.0, origin, x, y, dx, dy);
    else if(GLUTbutton[1])
        viewer->ScaleWorld(1.0, origin, x, y, dx, dy);
    else if(GLUTbutton[2]) {
        viewer->TranslateWorld(1.0, origin, x, y, dx, dy);
    }

    // redisplay if a mouse was down
    if(GLUTbutton[0] || GLUTbutton[1] || GLUTbutton[2]) glutPostRedisplay();

    // remember mouse position
    GLUTmouse[0] = x;
    GLUTmouse[1] = y;
}



void GLUTMouse(int button, int state, int x, int y)
{
    // invert y coordinate
    y = GLUTwindow_height - y;

    // remember mouse position
    GLUTmouse[0] = x;
    GLUTmouse[1] = y;

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



void GLUTSpecial(int key, int x, int y)
{
    // invert y coordinate
    y = GLUTwindow_height - y;

    // remember mouse position
    GLUTmouse[0] = x;
    GLUTmouse[1] = y;

    // remember modifiers
    GLUTmodifiers = glutGetModifiers();

    switch(key) {
        case GLUT_KEY_LEFT: {
            candidate_index--;
            if(candidate_index < 0) candidate_index = 0;
            ReadFeature(candidate_index);
            break;
        }

        case GLUT_KEY_RIGHT: {
            candidate_index++;
            if(candidate_index > (int)decisions.size() - 1)
                candidate_index = decisions.size() - 1;
            ReadFeature(candidate_index);
            break;
        }
    }

    // redraw
    glutPostRedisplay();
}



void GLUTKeyboard(unsigned char key, int x, int y)
{
    // invert y coordinate
    y = GLUTwindow_height - y;

    // remember mouse position
    GLUTmouse[0] = x;
    GLUTmouse[1] = y;

    // remember modifiers
    GLUTmodifiers = glutGetModifiers();

    // keys regardless of projection status
    switch(key) {
        case ESCAPE: {
            GLUTStop();
            break;
        }
    }

    // redraw
    glutPostRedisplay();
}



void GLUTInit(int* argc, char** argv)
{
    // open window
    glutInit(argc, argv);
    glutInitWindowPosition(10, 10);
    glutInitWindowSize(GLUTwindow_width, GLUTwindow_height);
    glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH);
    GLUTwindow = glutCreateWindow("Feature Visualizer");

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
    glutKeyboardFunc(GLUTKeyboard);
    glutSpecialFunc(GLUTSpecial);
    glutMouseFunc(GLUTMouse);
    glutMotionFunc(GLUTMotion);
}



void GLUTMainLoop(void)
{
    // run main loop -- never returns
    glutMainLoop();
}



////////////////////////////////////////////////////////////////////////
// Viewer functions
////////////////////////////////////////////////////////////////////////

static R3Viewer* CreateViewer(void)
{
    // start statistics
    RNTime start_time;
    start_time.Read();

    if(world_box.IsEmpty()) RNAbort("Error in CreateViewer - box is empty");
    RNLength r = world_box.DiagonalRadius();
    if(r < 0 || RNIsInfinite(r)) RNAbort("Error in CreateViewer - r must be positive finite");

    // set up camera view looking down the z axis
    static R3Vector initial_camera_towards = R3Vector(0.0, 0.0, -1.5);
    static R3Vector initial_camera_up = R3Vector(0.0, 1.0, 0.0);
    R3Point initial_camera_origin = world_box.Centroid() - initial_camera_towards * 2.5 * r;
    R3Camera camera(initial_camera_origin, initial_camera_towards, initial_camera_up, 0.4, 0.4, 0.1 * r, 1000.0 * r);
    R2Viewport viewport(0, 0, GLUTwindow_width, GLUTwindow_height);
    R3Viewer* viewer = new R3Viewer(camera, viewport);

    // print statistics
    if(print_verbose) {
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

static int ParseArgs(int argc, char** argv)
{
    // parse arguments
    argc--; argv++;
    while (argc > 0) {
        if ((*argv)[0] == '-') {
            if (!strcmp(*argv, "-v")) print_verbose = 1;
            else if (!strcmp(*argv, "-debug")) print_debug = 1;
            else if (!strcmp(*argv, "-threshold")) { argv++; argc--; threshold = atoi(*argv); }
            else if (!strcmp(*argv, "-window_radius")) { argv++; argc--; window_radius = atoi(*argv); }
            else { fprintf(stderr, "Invalid program argument: %s\n", *argv); return 0; }
        } else {
            if (!prefix_one) { prefix_one = *argv; } 
            else if (!prefix_two) { prefix_two = *argv; }
            else { fprintf(stderr, "Invalid program argument: %s\n", *argv); return 0; }
        }
        argv++; argc--;
    }

    // error if there is no input name
    if (!prefix_one || !prefix_two) { fprintf(stderr, "Need to supply a prefix file\n"); return 0; }

    // return success
    return 1;
}



////////////////////////////////////////////////////////////////////////
// Main program
////////////////////////////////////////////////////////////////////////

int main(int argc, char** argv)
{
    // parse program arguments
    if (!ParseArgs(argc, argv)) exit(-1);

    //////////////////////////////////////////////////
    //// Read in the voxel and ground truth files ////
    //////////////////////////////////////////////////

    // read in the ground truth file
    if (!ReadDecisionFile()) exit(-1);

    // read the first feature
    if (!ReadFeature(candidate_index)) exit(-1);

    if (print_verbose) {
        printf("Resolution = (%d %d %d)\n", grid->XResolution(), grid->YResolution(), grid->ZResolution());
    }

    // set world box
    world_box = R3Box(0, 0, 0, grid->XResolution(), grid->YResolution(), grid->ZResolution());
    
    // get the transformation
    transformation = R3Affine(R4Matrix(resolution[RN_X], 0, 0, 0, 0, resolution[RN_Y], 0, 0, 0, 0, resolution[RN_Z], 0, 0, 0, 0, 1));

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