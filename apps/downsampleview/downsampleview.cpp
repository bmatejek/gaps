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
#include <algorithm>


// GLUT defines

#define ENTER 13
#define ESCAPE 27
#define SPACEBAR 32
#define IB_Z 0
#define IB_Y 1
#define IB_X 2



// program arguments

// I/O flags
static int print_debug = 0;
static int print_verbose = 0;
// dataset to examine
static const char* prefix = NULL;



// program variables

static long resolution[3] = { 100, 100, 100 };
static long grid_size[3] = { -1, -1, -1 };
static R3Affine transformation = R3null_affine;
static R3Viewer *viewer = NULL;
static R3Box world_box = R3null_box;



// display variables

static int projection_dim = RN_Z;
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
static std::vector<long> *thinning_skeletons = NULL;
static std::vector<long> *medial_skeletons = NULL;
static std::vector<long> *teaser_skeletons = NULL;
static long maximum_segmentation = -1;



// display variables

static int show_bbox = 1;
static int segmentation_index = 1;
static int skeleton_type = 0;
static RNScalar downsample_rate = 2.0;




////////////////////////////////////////////////////////////////////////
// Input/output functions
////////////////////////////////////////////////////////////////////////


static int ReadData(void)
{
    // start statistics
    RNTime start_time;
    start_time.Read();

    char input_filename[4096];
    sprintf(input_filename, "skeletons/%s/topological-downsample-%ldx%ldx%ld.bytes", prefix, resolution[IB_X], resolution[IB_Y], resolution[IB_Z]);

    FILE *fp = fopen(input_filename, "rb"); 
    if (!fp) { fprintf(stderr, "Failed to read %s\n", input_filename); return 0; }

    if (fread(&(grid_size[RN_Z]), sizeof(long), 1, fp) != 1) return 0;
    if (fread(&(grid_size[RN_Y]), sizeof(long), 1, fp) != 1) return 0;
    if (fread(&(grid_size[RN_X]), sizeof(long), 1, fp) != 1) return 0;

    // read the number of segments
    if (fread(&maximum_segmentation, sizeof(long), 1, fp) != 1) return 0;
    segmentations = new std::vector<long>[maximum_segmentation];
    for (long iv = 0; iv < maximum_segmentation; ++iv) {
        segmentations[iv] = std::vector<long>();
        long nelements; 
        if (fread(&nelements, sizeof(long), 1, fp) != 1) return 0;
        for (long ie = 0; ie < nelements; ++ie) {
            long element;
            if (fread(&element, sizeof(long), 1, fp) != 1) return 0;
            segmentations[iv].push_back(element);
        }
    }

    fclose(fp);

    sprintf(input_filename, "skeletons/%s/topological-downsample-%ldx%ldx%ld-thinning-skeleton.pts", prefix, resolution[IB_X], resolution[IB_Y], resolution[IB_Z]);
    fp = fopen(input_filename, "rb");
    long skeleton_maximum_segmentation;
    if (fp) {
        if (fread(&skeleton_maximum_segmentation, sizeof(long), 1, fp) != 1) return 0;
        assert (skeleton_maximum_segmentation == maximum_segmentation);

        thinning_skeletons = new std::vector<long>[maximum_segmentation];
        for (long iv = 0; iv < maximum_segmentation; ++iv) {
            thinning_skeletons[iv] = std::vector<long>();

            long nelements; 
            if (fread(&nelements, sizeof(long), 1, fp) != 1) return 0;
            for (long ie = 0; ie < nelements; ++ie) {
                long element;
                if (fread(&element, sizeof(long), 1, fp) != 1) return 0;
                thinning_skeletons[iv].push_back(abs(element));
            }
        }  
        fclose(fp);
    }

    sprintf(input_filename, "skeletons/%s/topological-downsample-%ldx%ldx%ld-medial-axis-skeleton.pts", prefix, resolution[IB_X], resolution[IB_Y], resolution[IB_Z]);
    fp = fopen(input_filename, "rb");
    if (fp) {
        if (fread(&skeleton_maximum_segmentation, sizeof(long), 1, fp) != 1) return 0;
        assert (skeleton_maximum_segmentation == maximum_segmentation);

        medial_skeletons = new std::vector<long>[maximum_segmentation];
        for (long iv = 0; iv < maximum_segmentation; ++iv) {
            medial_skeletons[iv] = std::vector<long>();

            long nelements;
            if (fread(&nelements, sizeof(long), 1, fp) != 1) return 0;
            for (long ie = 0; ie < nelements; ++ie) {
                long element;
                if (fread(&element, sizeof(long), 1, fp) != 1) return 0;
                medial_skeletons[iv].push_back(abs(element));
            }
        }
        fclose(fp);
    }

    sprintf(input_filename, "skeletons/%s/topological-downsample-%ldx%ldx%ld-teaser-skeleton.pts", prefix, resolution[IB_X], resolution[IB_Y], resolution[IB_Z]);
    fp = fopen(input_filename, "rb");
    if (fp) {
        if (fread(&skeleton_maximum_segmentation, sizeof(long), 1, fp) != 1) return 0;
        assert (skeleton_maximum_segmentation == maximum_segmentation);

        teaser_skeletons = new std::vector<long>[maximum_segmentation];
        for (long iv = 0; iv < maximum_segmentation; ++iv) {
            teaser_skeletons[iv] = std::vector<long>();

            long nelements;
            if (fread(&nelements, sizeof(long), 1, fp) != 1) return 0;
            for (long ie = 0; ie < nelements; ++ie) {
                long element;
                if (fread(&element, sizeof(long), 1, fp) != 1) return 0;
                teaser_skeletons[iv].push_back(abs(element));
            }
        }
        fclose(fp);
    }

    // print statistics
    if(print_verbose) {
        printf("Read voxel grids...\n");
        printf("  Time = %.2f seconds\n", start_time.Elapsed());
        printf("  Grid Size = (%ld %ld %ld)\n", grid_size[RN_X], grid_size[RN_Y], grid_size[RN_Z]);
        printf("  Resolution = (%ld %ld %ld)\n", resolution[RN_X], resolution[RN_Y], resolution[RN_Z]);
    }

    // return success
    return 1;
}



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



static RNRgb Color(unsigned long value)
{
    RNScalar red = (RNScalar) (((107 * value) % 700) % 255) / 255.0;
    RNScalar green = (RNScalar) (((509 * value) % 900) % 255) / 255.0;
    RNScalar blue = (RNScalar) (((200 * value) % 777) % 255) / 255.0;

    return RNRgb(red, green, blue);
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


static void DrawSkeleton(int segment_index) {
    std::vector<long> *skeletons;
    if (skeleton_type == 0) skeletons = thinning_skeletons;
    else if (skeleton_type == 1) skeletons = medial_skeletons;
    else if (skeleton_type == 2) skeletons = teaser_skeletons;
    else return;

    glBegin(GL_POINTS);
    for (unsigned long iv = 0; iv < skeletons[segment_index].size(); ++iv) {
        // get the coordinates from the linear index
        int ix, iy, iz;
        IndexToIndices(skeletons[segment_index][iv], ix, iy, iz);
        glVertex3f(ix, iy, iz);
    }  
    glEnd();
}



static void DrawPointClouds(void)
{
    // push the transformation
    transformation.Push();

    glPointSize(2);
    RNLoadRgb(Color(segmentation_index));
    DrawSegment(segmentation_index);
    glPointSize(3);
    RNLoadRgb(RNRgb(1.0 - background_color[0], 1.0 - background_color[1], 1.0 - background_color[2]));
    DrawSkeleton(segmentation_index);

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

    delete[] segmentations;

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
            --segmentation_index;
            if (segmentation_index < 0)
                segmentation_index = 0;

            printf("Label: %d\n", segmentation_index);
            break;
        }

        case GLUT_KEY_RIGHT: {
            ++segmentation_index;
            if (segmentation_index >= maximum_segmentation)
                segmentation_index = maximum_segmentation - 1;

            printf("Label: %d\n", segmentation_index);
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

        case 'K':
        case 'k': {
            skeleton_type = (++skeleton_type) % 4;
        }

        case 'X':
        case 'x': {
            projection_dim = RN_X;
            break;
        }

        case 'Y':
        case 'y': {
            projection_dim = RN_Y;
            break;
        }

        case 'Z':
        case 'z': {
            projection_dim = RN_Z;
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
            if(!strcmp(*argv, "-v")) print_verbose = 1;
            else if(!strcmp(*argv, "-debug")) print_debug = 1;
            else { fprintf(stderr, "Invalid program argument: %s\n", *argv); return 0; }
        } else {
            if (!prefix) prefix = *argv;
            else { fprintf(stderr, "Invalid program argument: %s\n", *argv); return 0; }
        }
        argv++; argc--;
    }

    // error if there is no input name
    if (!prefix) {
        fprintf(stderr, "Need to supply prefix\n");
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

    /////////////////////////////////
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

    // initialize GLUT
    GLUTInit(&argc, argv);

    // run GLUT interface
    GLUTMainLoop();

    // return success
    return 1;
}