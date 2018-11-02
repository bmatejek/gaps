// Source file for skeleton visualizer



// include files

#include "R3Graphics/R3Graphics.h"
#include "RNDataStructures/RNDataStructures.h"
#include "fglut/fglut.h"
#include <fstream>
#include <string>
#include <vector>
#include <map>
#include <set>
#include <string>
#include <algorithm>
#include <dirent.h>


// GLUT defines

#define ENTER 13
#define ESCAPE 27
#define SPACEBAR 32




// program arguments

// I/O flags
static int print_debug = 0;
static int print_verbose = 0;
// dataset to examine
static const char *feature_set = NULL;
static const char *dataset = NULL;


// program variables

static RNScalar resolution[3] = { 3.6, 3.6, 40 };
static int grid_size[3] = { -1, -1, -1 };
static R3Affine transformation = R3null_affine;
static R3Viewer *viewer = NULL;
static R3Box world_box = R3null_box;



// voxel grids

static std::vector<R3Grid *> positive_examples = std::vector<R3Grid *>();
static std::vector<R3Grid *> negative_examples = std::vector<R3Grid *>();
static std::vector<R3Grid *> unknown_examples = std::vector<R3Grid *>();
static R3Grid *grid = NULL;


// display variables

static double background_color[3] = { 0.0, 0.0, 0.0 };



// GLUT variables

static int GLUTwindow = 0;
static int GLUTwindow_height = 1200;
static int GLUTwindow_width = 1200;
static int GLUTmouse[2] = { 0, 0 };
static int GLUTbutton[3] = { 0, 0, 0 };
static int GLUTmodifiers = 0;



// random access variables

static std::vector<long> *segmentations = NULL;
static long maximum_segmentation = -1;



// datasets for positives, negatives, and unknowns

static std::vector<std::string> positive_filenames;
static std::vector<std::string> negative_filenames;
static std::vector<std::string> unknown_filenames;



// display variables

static int show_bbox = 1;
static RNScalar downsample_rate = 2.0;
static int show_positives = 1;
static int show_negatives = 0;
static int show_unknowns = 0;

static int positive_index = 0;
static int negative_index = 0;
static int unknown_index = 0;



////////////////////////////////////////////////////////////////////////
// Helper functions
////////////////////////////////////////////////////////////////////////

static void IndexToIndices(int index, int& ix, int& iy, int& iz)
{
  // Set indices of grid value at index
  iz = index / (grid_size[RN_X] * grid_size[RN_Y]);
  iy = (index - iz * grid_size[RN_X] * grid_size[RN_Y]) / grid_size[RN_X];
  ix = index % grid_size[RN_X];
}



static long IndicesToIndex(long ix, long iy, long iz)
{
    return iz * grid_size[RN_X] * grid_size[RN_Y] + iy * grid_size[RN_X] + ix;
}



////////////////////////////////////////////////////////////////////////
// Processing functions
////////////////////////////////////////////////////////////////////////

static void ProcessFeature(void)
{
    RNTime start_time;
    start_time.Read();

    if (segmentations) delete[] segmentations;

    // get the maximum values for each grid
    maximum_segmentation = (long) (grid->Maximum() + 0.5) + 1;

    // create a vector for each valid ID
    segmentations = new std::vector<long>[maximum_segmentation];
    for (long iv = 0; iv < maximum_segmentation; ++iv)
        segmentations[iv] = std::vector<long>();

    // go through all voxels to see if it belongs to the boundary
    for (int iz = 1; iz < grid_size[RN_Z] - 1; ++iz) {
        for (int iy = 1; iy < grid_size[RN_Y] - 1; ++iy) {
            for (int ix = 1; ix < grid_size[RN_X] - 1; ++ix) {
                long iv = IndicesToIndex(ix, iy, iz);

                long segment = (long) (grid->GridValue(ix, iy, iz) + 0.5);

                // go through all six adjacent neighbors
                if ((long)(grid->GridValue(ix - 1, iy, iz) + 0.5) != segment ||
                    (long)(grid->GridValue(ix + 1, iy, iz) + 0.5) != segment ||
                    (long)(grid->GridValue(ix, iy - 1, iz) + 0.5) != segment ||
                    (long)(grid->GridValue(ix, iy + 1, iz) + 0.5) != segment ||
                    (long)(grid->GridValue(ix, iy, iz - 1) + 0.5) != segment ||
                    (long)(grid->GridValue(ix, iy, iz + 1) + 0.5) != segment)
                {
                    segmentations[segment].push_back(iv);
                }
            }
        }
    }
}



////////////////////////////////////////////////////////////////////////
// Input/output functions
////////////////////////////////////////////////////////////////////////


static int ReadData(void)
{
    // start statistics
    RNTime start_time;
    start_time.Read();

    long ngrids = 0;

    // read all of the filenames
    for (unsigned long is = 0; is < positive_filenames.size(); ++is) {
        std::string filename = positive_filenames[is];

        // read the example to see how many there are 
        FILE *fp = fopen(filename.c_str(), "rb");
        if (!fp) return 0;
        long nexamples;
        if (fread(&nexamples, sizeof(long), 1, fp) != 1) return 0;
        ngrids += nexamples;
        fclose(fp);

        // get the h5 filename
        char tmp_filename[4096];
        sprintf(tmp_filename, "%s", filename.c_str());
        char *extp = strrchr(tmp_filename, '.');
        *extp = '\0';
        char example_filename[4096];
        sprintf(example_filename, "%s-examples.h5", tmp_filename); 

        // read all of the grids
        R3Grid **grids = RNReadH5File(example_filename, "main");

        for (long ie = 0; ie < nexamples; ++ie) {
            if (ie % 10) { delete grids[ie]; continue; }
            positive_examples.push_back(grids[ie]);
        }

        // delete wrapper array
        delete[] grids;
    }

    // read all of the filenames
    for (unsigned long is = 0; is < negative_filenames.size(); ++is) {
        std::string filename = negative_filenames[is];

        // read the example to see how many there are 
        FILE *fp = fopen(filename.c_str(), "rb");
        if (!fp) return 0;
        long nexamples;
        if (fread(&nexamples, sizeof(long), 1, fp) != 1) return 0;
        ngrids += nexamples;
        fclose(fp);

        // get the h5 filename
        char tmp_filename[4096];
        sprintf(tmp_filename, "%s", filename.c_str());
        char *extp = strrchr(tmp_filename, '.');
        *extp = '\0';
        char example_filename[4096];
        sprintf(example_filename, "%s-examples.h5", tmp_filename); 

        // read all of the grids
        R3Grid **grids = RNReadH5File(example_filename, "main");

        for (long ie = 0; ie < nexamples; ++ie) {
            if (ie % 10) { delete grids[ie]; continue; }
            negative_examples.push_back(grids[ie]);
        }

        // delete wrapper array
        delete[] grids;
    }

    if (print_verbose) {
        printf("Read %ld grids in %0.2f seconds\n", ngrids, start_time.Elapsed());
    }

    // set the grid
    grid = positive_examples[positive_index];
    
    grid_size[RN_X] = grid->XResolution();
    grid_size[RN_Y] = grid->YResolution();
    grid_size[RN_Z] = grid->ZResolution();

    ProcessFeature();



    // return success
    return 1;
}





////////////////////////////////////////////////////////////////////////
// Drawing utility functions
////////////////////////////////////////////////////////////////////////



static void DrawSegment(int segment_index)
{
    glBegin(GL_POINTS);
    for (unsigned int iv = 0; iv < segmentations[segment_index].size(); ++iv) {
        // faster rendering with downsampling
        if (RNRandomScalar() > 1.0 / downsample_rate) continue;

        // get the coordinates from the linear index
        int ix, iy, iz;
        IndexToIndices(segmentations[segment_index][iv], ix, iy, iz);
        glVertex3f(ix, iy, iz);
    }
    glEnd();
}



static void DrawPointClouds(void)
{
    // push the transformation
    transformation.Push();

    if (show_positives) RNLoadRgb(RNgreen_rgb);
    if (show_negatives) RNLoadRgb(RNred_rgb);
    if (show_unknowns) RNLoadRgb(RNblue_rgb);

    DrawSegment(1);

    if (show_positives) RNLoadRgb(RNblue_rgb);
    if (show_negatives) RNLoadRgb(RNyellow_rgb);
    if (show_unknowns) RNLoadRgb(RNyellow_rgb);

    DrawSegment(2);

    // pop the transformation
    transformation.Pop();
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

    // drawing varies on projection
    viewer->Camera().Load();

    // set lights
    static GLfloat light0_position[] = { 3.0, 4.0, 5.0, 0.0 };
    glLightfv(GL_LIGHT0, GL_POSITION, light0_position);
    static GLfloat light1_position[] = { -3.0, -2.0, -3.0, 0.0 };
    glLightfv(GL_LIGHT1, GL_POSITION, light1_position);

    // prologue
    glDisable(GL_LIGHTING);


    // draw neuron data bounding box
    if(show_bbox) {
        RNLoadRgb(RNRgb(not background_color[0], not background_color[1], not background_color[2]));
        world_box.Outline();
    }

    // draw machine labels and skeletons
    DrawPointClouds();

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

    // compute mouse movement
    int dx = x - GLUTmouse[0];
    int dy = y - GLUTmouse[1];

    // world in hand navigation
    R3Point origin = world_box.Centroid();
    if(GLUTbutton[0])
        viewer->RotateWorld(1.0, origin, x, y, dx, dy);
    else if(GLUTbutton[1])
        viewer->ScaleWorld(1.0, origin, x, y, dx, dy);
    else if(GLUTbutton[2])
        viewer->TranslateWorld(1.0, origin, x, y, dx, dy);
    
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
            if (show_positives) {
                --positive_index;
                if (positive_index < 0) 
                    positive_index = 0;
                grid = positive_examples[positive_index];
            }
            if (show_negatives) {
                --negative_index;
                if (negative_index < 0)
                    negative_index = 0;
                grid = negative_examples[negative_index];
            }
            if (show_unknowns) {
                --unknown_index;
                if (unknown_index < 0)
                    unknown_index = 0;
                grid = unknown_examples[unknown_index];
            }

            ProcessFeature();

            break;
        }

        case GLUT_KEY_RIGHT: {
            if (show_positives) {
                ++positive_index;
                if (positive_index > (long)(positive_examples.size() - 1))
                    positive_index = positive_examples.size() - 1;
                grid = positive_examples[positive_index];
            }
            if (show_negatives) {
                ++negative_index;
                if (negative_index > (long)(negative_examples.size() - 1))
                    negative_index = negative_examples.size() - 1;
                grid = negative_examples[negative_index];
            }
            if (show_unknowns) {
                ++unknown_index;
                if (unknown_index > (long)(unknown_examples.size() - 1))
                    unknown_index = unknown_examples.size() - 1;
                grid = unknown_examples[unknown_index];
            }

            ProcessFeature();

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
        case 'B':
        case 'b': {
            show_bbox = 1 - show_bbox;
            break;
        }

        case 'N': 
        case 'n': {
            show_negatives = 1;
            show_positives = 0;
            show_unknowns = 0;

            grid = negative_examples[negative_index];

            ProcessFeature();

            break;
        }

        case 'P': 
        case 'p': {
            show_negatives = 0;
            show_positives = 1;
            show_unknowns = 0;

            grid = positive_examples[positive_index];
            
            ProcessFeature();

            break;
        }

        case 'U': 
        case 'u': {
            show_negatives = 0;
            show_positives = 0;
            show_unknowns = 1;
                
            grid = unknown_examples[unknown_index];

            ProcessFeature();

            break;
        }

        case ENTER: {
            background_color[0] = 1.0 - background_color[0];
            background_color[1] = 1.0 - background_color[1];
            background_color[2] = 1.0 - background_color[2];
            break;
        }


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
    GLUTwindow = glutCreateWindow("Skeleton Visualizer");

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

    if (world_box.IsEmpty()) RNAbort("Error in CreateViewer - box is empty");
    RNLength radius = world_box.DiagonalRadius();
    if (radius < 0 || RNIsInfinite(radius)) RNAbort("Error in CreateViewer - radius must be positive finite");

    // set up camera view looking down the z axis
    static R3Vector initial_camera_towards = R3Vector(0.0, 0.0, -1.5);
    static R3Vector initial_camera_up = R3Vector(0.0, 1.0, 0.0);
    R3Point initial_camera_origin = world_box.Centroid() - initial_camera_towards * 2.5 * radius;
    R3Camera camera(initial_camera_origin, initial_camera_towards, initial_camera_up, 0.4, 0.4, 0.1 * radius, 1000.0 * radius);
    R2Viewport viewport(0, 0, GLUTwindow_width, GLUTwindow_height);
    R3Viewer *viewer = new R3Viewer(camera, viewport);

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
    while(argc > 0) {
        if((*argv)[0] == '-') {
            if (!strcmp(*argv, "-v")) print_verbose = 1;
            else if (!strcmp(*argv, "-debug")) print_debug = 1;
            else if (!strcmp(*argv, "-resolution")) { 
                argv++; argc--; resolution[RN_X] = atof(*argv); 
                argv++; argc--; resolution[RN_Y] = atof(*argv);
                argv++; argc--; resolution[RN_Z] = atof(*argv);
            }
            else { fprintf(stderr, "Invalid program argument: %s\n", *argv); return 0; }
        } else {
            if (!feature_set) feature_set = *argv;
            else if (!dataset) dataset = *argv;
            else { fprintf(stderr, "Invalid program argument: %s\n", *argv); return 0; }
        }
        argv++; argc--;
    }

    // error if there is no input name
    if (!feature_set || !dataset) {
        fprintf(stderr, "Need to supply feature_set and dataset for data files\n");
        return 0;
    }

    // return success
    return 1;
}



////////////////////////////////////////////////////////////////////////
// Main program
////////////////////////////////////////////////////////////////////////

int main(int argc, char** argv)
{
    // parse program arguments
    if(!ParseArgs(argc, argv)) exit(-1);

    // get all of the subdirectories
    char positive_directory[4096];
    sprintf(positive_directory, "features/biological/%s/%s/positives", feature_set, dataset);
    char negative_directory[4096];
    sprintf(negative_directory, "features/biological/%s/%s/negatives", feature_set, dataset);
    char unknown_directory[4096];
    sprintf(unknown_directory, "features/biological/%s/%s/unknowns", feature_set, dataset);

    // create empty vectors
    positive_filenames = std::vector<std::string>();
    negative_filenames = std::vector<std::string>();
    unknown_filenames = std::vector<std::string>();

    // read all files in directory
    DIR *dir;
    struct dirent *ent;
    if ((dir = opendir(positive_directory)) != NULL) {
        while ((ent = readdir(dir)) != NULL) {
            if (ent->d_type != DT_REG) continue;
            // read the filename
            char filename[4096];
            sprintf(filename, "%s/%s", positive_directory, ent->d_name);
            
            // make sure that it is not h5
            std::string example_filename = std::string(filename);
            if (not (example_filename.substr(example_filename.find_last_of(".") + 1) == "examples")) continue;

            // add to the vector of files to read
            positive_filenames.push_back(example_filename);
        }
        closedir(dir);
    }
    if ((dir = opendir(negative_directory)) != NULL) {
        while ((ent = readdir(dir)) != NULL) {
            if (ent->d_type != DT_REG) continue;
            // read the filename
            char filename[4096];
            sprintf(filename, "%s/%s", negative_directory, ent->d_name);
            
            // make sure that it is not h5
            std::string example_filename = std::string(filename);
            if (not (example_filename.substr(example_filename.find_last_of(".") + 1) == "examples")) continue;

            // add to the vector of files to read
            negative_filenames.push_back(example_filename);
        }
        closedir(dir);
    }
    if ((dir = opendir(unknown_directory)) != NULL) {
        while ((ent = readdir(dir)) != NULL) {
            if (ent->d_type != DT_REG) continue;
            // read the filename
            char filename[4096];
            sprintf(filename, "%s/%s", unknown_directory, ent->d_name);
            
            // make sure that it is not h5
            std::string example_filename = std::string(filename);
            if (not (example_filename.substr(example_filename.find_last_of(".") + 1) == "examples")) continue;

            // add to the vector of files to read
            unknown_filenames.push_back(example_filename);
        }
        closedir(dir);
    }

    if (not positive_filenames.size() or not negative_filenames.size()) {
        printf("No featuers found...\n");
        exit(-1);
    }
    if (not unknown_filenames.size()) {
        printf("No unknowns. Continuing...\n");
    }


    ////////////////////////////////
    //// Read in the voxel files ////
    /////////////////////////////////

    if (!ReadData()) exit(-1);

    // set world box
    world_box = R3Box(0, 0, 0, resolution[RN_X] * grid_size[RN_X], resolution[RN_Y] * grid_size[RN_Y], resolution[RN_Z] * grid_size[RN_Z]);

    // get the transformation
    transformation = R3Affine(R4Matrix(resolution[RN_X], 0, 0, 0, 0, resolution[RN_Y], 0, 0, 0, 0, resolution[RN_Z], 0, 0, 0, 0, 1));

    // create viewer
    viewer = CreateViewer();
    if (!viewer) exit(-1);

    printf("Feature sizes: \n");
    printf("  %d %d %d\n", grid->XResolution(), grid->YResolution(), grid->ZResolution());

    // initialize GLUT
    GLUTInit(&argc, argv);

    // run GLUT interface
    GLUTMainLoop();

    // return success
    return 1;
}