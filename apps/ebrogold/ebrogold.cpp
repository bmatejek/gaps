// Source file for generating ground truth for EbroNet


// include files

#include "R3Graphics/R3Graphics.h"
#include "RNDataStructures/RNDataStructures.h"
#include "fglut/fglut.h"
#include <vector>



// constants

static const int ngrids = 2;



// GLUT defines

#define ENTER 13
#define ESCAPE 27
#define SPACEBAR 32



////////////////////////////////////////////////////////////////////////
// Helper struct defintions
////////////////////////////////////////////////////////////////////////

struct RNMeta {
    // instance variables
    int resolution[3];
    char prefix[4096];
    char boundary_filename[4096];
    char boundary_dataset[128];
    char image_filename[4096];
    char image_dataset[128];
    char rhoana_filename[4096];
    char rhoana_dataset[128];
    R3Box world_box;
};



// program arguments

static int print_verbose = 0;
static const char *prefixes[ngrids] = { NULL, NULL };



// program variables

static int resolution[3] = { 4, 4, 30 };
static RNMeta meta_data[ngrids];
static R3Affine transformation = R3null_affine;
static R3Viewer *viewer = NULL;
static R3Box world_box = R3null_box;



// voxel grids

static R3Grid *boundary_grids[ngrids] = { NULL, NULL };
static R3Grid *image_grids[ngrids] = { NULL, NULL };
static R3Grid *rhoana_grids[ngrids] = { NULL, NULL };



// GLUT variables

static int GLUTwindow = 0;
static int GLUTwindow_height = 1200;
static int GLUTwindow_width = 1200;
static int GLUTmouse[2] = { 0, 0 };
static int GLUTbutton[3] = { 0, 0, 0 };
static int GLUTmodifiers = 0;



// color arrays

static RNScalar background_color[3] = { 0, 0, 0 };





////////////////////////////////////////////////////////////////////////
// Input/output functions
////////////////////////////////////////////////////////////////////////

static int
ReadMetaData(const char *prefix, int index)
{
    // update the meta data prefix
    strncpy(meta_data[index].prefix, prefix, 4096);

    // get the meta data filename
    char meta_data_filename[4096];
    sprintf(meta_data_filename, "meta_data/%s.meta", prefix);

    // open the file
    FILE *fp = fopen(meta_data_filename, "r");
    if (!fp) { fprintf(stderr, "Failed to read %s\n", meta_data_filename); return 0; }

    // dummy variable
    char comment[4096];

    // read in requisite information
    if (!fgets(comment, 4096, fp)) return 0;
    if (fscanf(fp, "%dx%dx%d\n", &(meta_data[index].resolution[RN_X]), &(meta_data[index].resolution[RN_Y]), &(meta_data[index].resolution[RN_Z])) != 3) return 0;
    if (!fgets(comment, 4096, fp)) return 0;
    if (fscanf(fp, "%s %s\n", meta_data[index].boundary_filename, meta_data[index].boundary_dataset) != 2) return 0;
    if (!fgets(comment, 4096, fp)) return 0;
    if (fscanf(fp, "%s %s\n", meta_data[index].image_filename, meta_data[index].image_dataset) != 2) return 0;
    if (!fgets(comment, 4096, fp)) return 0;
    if (fscanf(fp, "%s %s\n", meta_data[index].rhoana_filename, meta_data[index].rhoana_dataset) != 2) return 0;

    // get the world box
    RNScalar xlow, ylow, zlow, xhigh, yhigh, zhigh;
    if (!fgets(comment, 4096, fp)) return 0;
    if (fscanf(fp, "(%lf,%lf,%lf)-(%lf,%lf,%lf)\n", &xlow, &ylow, &zlow, &xhigh, &yhigh, &zhigh) != 6) return 0;
    meta_data[index].world_box = R3Box(xlow, ylow, zlow, xhigh, yhigh, zhigh);

    // close the file
    fclose(fp);

    // return success
    return 1;
}



static int
ReadVoxelGrids(int index)
{
    // start statistics
    RNTime start_time;
    start_time.Read();

    // read in all of the voxel grids
    R3Grid **grids = RNReadH5File(meta_data[index].boundary_filename, meta_data[index].boundary_dataset);
    if (!grids) { fprintf(stderr, "Failed to read %s\n", meta_data[index].boundary_filename); return 0; }
    boundary_grids[index] = grids[0];
    delete[] grids;

    grids = RNReadH5File(meta_data[index].image_filename, meta_data[index].image_dataset);
    if (!grids) { fprintf(stderr, "Failed to read %s\n", meta_data[index].image_filename); return 0; }
    image_grids[index] = grids[0];
    delete[] grids;

    grids = RNReadH5File(meta_data[index].rhoana_filename, meta_data[index].rhoana_dataset);
    if (!grids) { fprintf(stderr, "Failed to read %s\n", meta_data[index].rhoana_filename); return 0; }
    rhoana_grids[index] = grids[0];
    delete[] grids;

    // print information
    if (print_verbose) {
        printf("Read h5 files in %0.2f seconds for %s:\n", start_time.Elapsed(), meta_data[index].prefix);
        printf("  Dimensions: (%d, %d, %d)\n", boundary_grids[index]->XResolution(), boundary_grids[index]->YResolution(), boundary_grids[index]->ZResolution());
        printf("  Resolution: (%d, %d, %d)\n", resolution[RN_X], resolution[RN_Y], resolution[RN_Z]);
    }


    // return success 
    return 1;
}




////////////////////////////////////////////////////////////////////////
// Program argument parsing
////////////////////////////////////////////////////////////////////////

static int 
ParseArgs(int argc, char** argv)
{
    // parse arguments
    argc--; argv++;
    while (argc > 0) {
        if ((*argv)[0] == '-') {
            if (!strcmp(*argv, "-v")) print_verbose = 1;
            else { fprintf(stderr, "Invalid program argument: %s\n", *argv); return 0; }
        } else {
            if (!prefixes[0]) { prefixes[0] = *argv; } 
            else if (!prefixes[1]) { prefixes[1] = *argv; }
            else { fprintf(stderr, "Invalid program argument: %s\n", *argv); return 0; }
        }
        argv++; argc--;
    }

    // error if there is no input name
    for (int iv = 0; iv < ngrids; ++iv) {
        if (!prefixes[iv]) { fprintf(stderr, "Need to supply a prefix\n"); return 0; }
    }

    // return success
    return 1;
}



////////////////////////////////////////////////////////////////////////
// Main program
////////////////////////////////////////////////////////////////////////

int 
main(int argc, char** argv) 
{
    // parse program arguments
    if (!ParseArgs(argc, argv)) exit(-1);

    /////////////////////////////////
    //// Read in the voxel files ////
    /////////////////////////////////

    for (int iv = 0; iv < ngrids; ++iv) {
        if (!ReadMetaData(prefixes[iv], iv)) exit(-1);
        if (!ReadVoxelGrids(iv)) exit(-1);
    }


    


    // return success
    return 1;
}
























/*// program arguments

static int print_debug = 0;
static int print_verbose = 0;
static const char* prefix = NULL;
static int window_width[3] = { -1, -1, -1 };



// feature variables

static std::vector<RNBoolean> ground_truth;
static R3Grid *grid = NULL;
static int feature_index = 0;



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



////////////////////////////////////////////////////////////////////////
// Input/output functions
////////////////////////////////////////////////////////////////////////

static int ReadGroundTruth(void)
{
    // get the ground truth file
    char ground_truth_filename[4096];
    sprintf(ground_truth_filename, "%s-ground-truth.txt", prefix);

    // open the file
    FILE *fp = fopen(ground_truth_filename, "r");
    if (!fp) { fprintf(stderr, "Failed to read %s\n", ground_truth_filename); return 0; }

    // create an empty vector for ground truth
    ground_truth = std::vector<RNBoolean>();
    
    // read in the file line by line
    char *line = NULL;
    size_t len = 0;
    while (getline(&line, &len, fp) != -1) {
        // add to the array of ground truth values
        ground_truth.push_back(atoi(line));
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
    sprintf(feature_filename, "%s/%05d-feature.h5", prefix, index);

    // read the feature into an R3Grid
    R3Grid **grids = RNReadH5File(feature_filename, "main");
    
    if (window_width[RN_X] > 0 && window_width[RN_Y] > 0 && window_width[RN_Z] > 0) {
        grid = new R3Grid(window_width[RN_X], window_width[RN_Y], window_width[RN_Z]);

        // iterate over all elements in the new grid
        for (int iz = 0; iz < window_width[RN_Z]; ++iz) {
            for (int iy = 0; iy < window_width[RN_Y]; ++iy) {
                for (int ix = 0; ix < window_width[RN_X]; ++ix) {
                    int iw = (int)(grids[0]->ZResolution() / (float)window_width[RN_Z] * iz);
                    int iv = (int)(grids[0]->YResolution() / (float)window_width[RN_Y] * iy);
                    int iu = (int)(grids[0]->XResolution() / (float)window_width[RN_X] * ix);

                    grid->SetGridValue(ix, iy, iz, grids[0]->GridValue(iu, iv, iw));
                }
            }
        }

        delete grids[0];
    }
    else {
        grid = grids[0];
    }

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

    // draw feature bounding box
    RNLoadRgb(RNwhite_rgb);
    world_box.Outline();

    // draw the actual points in 3D
    glBegin(GL_POINTS);
    for (int iz = 0; iz < grid->ZResolution(); ++iz) {
        for (int iy = 0; iy < grid->YResolution(); ++iy) {
            for (int ix = 0; ix < grid->XResolution(); ++ix) {
                int grid_value = (int)(grid->GridValue(ix, iy, iz) + 0.5);
                if (grid_value == 1) {
                    if (ground_truth[feature_index]) RNLoadRgb(RNgreen_rgb);
                    else RNLoadRgb(RNred_rgb);
                    glVertex3f(ix, iy, iz);
                }
                else if (grid_value == 2) {
                    if (ground_truth[feature_index]) RNLoadRgb(RNblue_rgb);
                    else RNLoadRgb(RNyellow_rgb);
                    glVertex3f(ix, iy, iz);
                }
            }
        }
    }
    glEnd();

    // write the feature
    char feature_label[4096];
    sprintf(feature_label, "Feature Visualizer - %d - Predicted: \n", feature_index);
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
            feature_index--;
            if(feature_index < 0) feature_index = 0;
            ReadFeature(feature_index);
            break;
        }

        case GLUT_KEY_RIGHT: {
            feature_index++;
            if(feature_index > (int)ground_truth.size() - 1)
                feature_index = ground_truth.size() - 1;
            ReadFeature(feature_index);
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
            else if (!strcmp(*argv, "-window_width")) {
                argv++; argc--; window_width[RN_X] = atoi(*argv);
                argv++; argc--; window_width[RN_Y] = atoi(*argv);
                argv++; argc--; window_width[RN_Z] = atoi(*argv);
            }
            else { fprintf(stderr, "Invalid program argument: %s\n", *argv); return 0; }
        } else {
            if (!prefix) { prefix = *argv; } 
            else { fprintf(stderr, "Invalid program argument: %s\n", *argv); return 0; }
        }
        argv++; argc--;
    }

    // error if there is no input name
    if(!prefix) { fprintf(stderr, "Need to supply a prefix file\n"); return 0; }

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
    if (!ReadGroundTruth()) exit(-1);

    // read the first feature
    if (!ReadFeature(feature_index)) exit(-1);

    if (print_verbose) {
        printf("Resolution = (%d %d %d)\n", grid->XResolution(), grid->YResolution(), grid->ZResolution());
    }

    // set world box
    world_box = R3Box(0, 0, 0, grid->XResolution(), grid->YResolution(), grid->ZResolution());
    
    // create viewer
    viewer = CreateViewer();
    if (!viewer) exit(-1);

    // initialize GLUT
    GLUTInit(&argc, argv);

    // run GLUT interface
    GLUTMainLoop();

    // return success
    return 1;
}*/