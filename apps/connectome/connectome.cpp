// Source file for connectome visualization



// include files

#include "R3Graphics/R3Graphics.h"
#include "RNDataStructures/RNDataStructures.h"
#include "fglut/fglut.h"
#include <fstream>
#include <vector>
#include <set>


// GLUT defines

#define ENTER 13
#define ESCAPE 27
#define SPACEBAR 32

static const int OR_Z = 0;
static const int OR_Y = 1;
static const int OR_X = 2;



// class declarations

struct RNMeta;



// program arguments

// I/O flags
static int print_debug = 0;
static int print_verbose = 0;
static const char* prefix = NULL;



// program variables

static double resolution[3] = { -1, -1, -1 };
static long grid_size[3] = { -1, -1, -1 };
static R3Affine transformation = R3null_affine;
static R3Viewer *viewer = NULL;
static R3Box world_box = R3null_box;



// voxel grids

static R3Grid *segmentation_grid = NULL;
static R3Grid *synapse_grid = NULL;



// connectome variables

static std::vector<long> *wirings = NULL;



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
static std::vector<long> *synapses = NULL;
static long maximum_synapse = -1;
static std::set<long> *synapses_per_segment = NULL;
static std::set<long> *segments_per_synapse = NULL;


// display variables

static const int ncolor_opts = 6;
static int show_bbox = 1;
static int segmentation_index = 1;
static int synapse_index = 1;
static int synapse_centric = 0;
static RNScalar downsample_rate = 2.0;
static int color_cycle = 0;
static const char *stages[2]  = { "thinning", "connectomes" };
static int stage_index = 1;



////////////////////////////////////////////////////////////////////////
// Input/output functions
////////////////////////////////////////////////////////////////////////

struct RNMeta {
    // instance variables
    double resolution[3];
    char prefix[4096];
    char segmentation_filename[4096];
    char segmentation_dataset[128];
    char synapse_filename[4096];
    char synapse_dataset[128];
};

// declare here
static RNMeta meta_data;



static int
ReadMetaData(const char *prefix)
{
    // update the meta data prefix
    strncpy(meta_data.prefix, prefix, 4096);

    // get the meta data filename
    char meta_data_filename[4096];
    sprintf(meta_data_filename, "meta/%s.meta", prefix);

    // open the file
    FILE *fp = fopen(meta_data_filename, "r");
    if (!fp) { fprintf(stderr, "Failed to read %s\n", meta_data_filename); return 0; }

    // dummy variable
    char comment[4096];

    // read in requisite information
    if (!fgets(comment, 4096, fp)) return 0;
    if (fscanf(fp, "%lfx%lfx%lf\n", &(meta_data.resolution[OR_X]), &(meta_data.resolution[OR_Y]), &(meta_data.resolution[OR_Z])) != 3) return 0;

    // read segmentation
    if (!fgets(comment, 4096, fp)) return 0;
    if (fscanf(fp, "%s %s\n", meta_data.segmentation_filename, meta_data.segmentation_dataset) != 2) return 0;

    // skip synapse
    if (!fgets(comment, 4096, fp)) return 0;
    if (fscanf(fp, "%s %s\n", meta_data.synapse_filename, meta_data.synapse_dataset) != 2) return 0;

    // update the global resolution
    for (int dim = 0; dim <= 2; ++dim)
        resolution[dim] = meta_data.resolution[dim];

    // close the file
    fclose(fp);

    // return success
    return 1;
}



static int ReadData(void)
{
    // start statistics
    RNTime start_time;
    start_time.Read();

    R3Grid **segmentation_grids = RNReadH5File(meta_data.segmentation_filename, meta_data.segmentation_dataset);
    segmentation_grid = segmentation_grids[0];
    delete[] segmentation_grids;

    R3Grid **synapse_grids = RNReadH5File(meta_data.synapse_filename, meta_data.synapse_dataset);
    synapse_grid = synapse_grids[0];
    delete[] synapse_grids;
    
    grid_size[OR_X] = segmentation_grid->XResolution();
    grid_size[OR_Y] = segmentation_grid->YResolution();
    grid_size[OR_Z] = segmentation_grid->ZResolution();

    // get the maximum values for each grid
    maximum_segmentation = (long) (segmentation_grid->Maximum() + 0.5) + 1;
    maximum_synapse = (long) (synapse_grid->Maximum() + 0.5) + 1;

    // print statistics
    if(print_verbose) {
        printf("Read voxel grids...\n");
        printf("  Time = %.2f seconds\n", start_time.Elapsed());
        printf("  Grid Size = (%ld %ld %ld)\n", grid_size[OR_X], grid_size[OR_Y], grid_size[OR_Z]);
        printf("  Resolution = (%0.2lf %0.2lf %0.2lf)\n", resolution[OR_X], resolution[OR_Y], resolution[OR_Z]);
    }

    // return success
    return 1;
}



static int
ReadConnectomeData(void)
{
    if (wirings) delete[] wirings;
    wirings = new std::vector<long>[maximum_segmentation];

    for (long iv = 0; iv < maximum_segmentation; ++iv) {
        wirings[iv] = std::vector<long>();
    
        char input_filename[4096];      
        sprintf(input_filename, "%s/%s/%06ld.pts", stages[stage_index], prefix, iv);

        FILE *fp = fopen(input_filename, "rb");
        if (!fp) { continue; }

        long nelements;
        long zres, yres, xres;
        if (fread(&zres, sizeof(long), 1, fp) != 1) return 0;
        if (fread(&yres, sizeof(long), 1, fp) != 1) return 0;
        if (fread(&xres, sizeof(long), 1, fp) != 1) return 0;
        if (fread(&nelements, sizeof(long), 1, fp) != 1) return 0;
        
        for (long ie = 0; ie < nelements; ++ie) {
            long element;
            if (fread(&element, sizeof(long), 1, fp) != 1) return 0;
            wirings[iv].push_back(element);
        }
        fclose(fp);
    }

    // return success
    return 1;
}


////////////////////////////////////////////////////////////////////////
// Helper functions
////////////////////////////////////////////////////////////////////////

static void IndexToIndices(long index, long& ix, long& iy, long& iz)
{
  // Set indices of grid value at index
  iz = index / (grid_size[OR_X] * grid_size[OR_Y]);
  iy = (index - iz * grid_size[OR_X] * grid_size[OR_Y]) / grid_size[OR_X];
  ix = index % grid_size[OR_X];
}



static long IndicesToIndex(long ix, long iy, long iz)
{
    return iz * grid_size[OR_X] * grid_size[OR_Y] + iy * grid_size[OR_X] + ix;
}



static RNRgb Color(unsigned long value)
{
    // allow alternating colors
    value += color_cycle;

    RNScalar red = (RNScalar) (((107 * value) % 700) % 255) / 255.0;
    RNScalar green = (RNScalar) (((509 * value) % 900) % 255) / 255.0;
    RNScalar blue = (RNScalar) (((200 * value) % 777) % 255) / 255.0;

    return RNRgb(red, green, blue);
}



////////////////////////////////////////////////////////////////////////
// Preprocessing functions
////////////////////////////////////////////////////////////////////////

static void Preprocessing(void)
{
    RNTime start_time;
    start_time.Read();

    // create a vector for each valid ID
    segmentations = new std::vector<long>[maximum_segmentation];
    synapses_per_segment = new std::set<long>[maximum_segmentation];
    for (long iv = 0; iv < maximum_segmentation; ++iv) {
        segmentations[iv] = std::vector<long>();
        synapses_per_segment[iv] = std::set<long>();
    }
    synapses = new std::vector<long>[maximum_synapse];
    segments_per_synapse = new std::set<long>[maximum_synapse];
    for (long iv = 0; iv < maximum_synapse; ++iv) {
        synapses[iv] = std::vector<long>();
        segments_per_synapse[iv] = std::set<long>();
    }


    // go through all voxels to see if it belongs to the boundary
    for (int iz = 1; iz < grid_size[OR_Z] - 1; ++iz) {
        for (int iy = 1; iy < grid_size[OR_Y] - 1; ++iy) {
            for (int ix = 1; ix < grid_size[OR_X] - 1; ++ix) {
                long iv = IndicesToIndex(ix, iy, iz);

                long segment = (long) (segmentation_grid->GridValue(ix, iy, iz) + 0.5);

                // skip background
                if (!segment) continue;

                // go through all six adjacent neighbors for the segmentations
                if ((long)(segmentation_grid->GridValue(ix - 1, iy, iz) + 0.5) != segment ||
                    (long)(segmentation_grid->GridValue(ix + 1, iy, iz) + 0.5) != segment ||
                    (long)(segmentation_grid->GridValue(ix, iy - 1, iz) + 0.5) != segment ||
                    (long)(segmentation_grid->GridValue(ix, iy + 1, iz) + 0.5) != segment ||
                    (long)(segmentation_grid->GridValue(ix, iy, iz - 1) + 0.5) != segment ||
                    (long)(segmentation_grid->GridValue(ix, iy, iz + 1) + 0.5) != segment)
                {
                    segmentations[segment].push_back(iv);
                }

                // add to this segment's list of synapses
                long synapse = (long) (synapse_grid->GridValue(ix, iy, iz) + 0.5);
                synapses_per_segment[segment].insert(synapse);
            }
        }
    }

    // go through all voxels to see if it belongs to the boundary
    for (int iz = 1; iz < grid_size[OR_Z] - 1; ++iz) {
        for (int iy = 1; iy < grid_size[OR_Y] - 1; ++iy) {
            for (int ix = 1; ix < grid_size[OR_X] - 1; ++ix) {
                long iv = IndicesToIndex(ix, iy, iz);

                long synapse = (long) (synapse_grid->GridValue(ix, iy, iz) + 0.5);
                
                // skip background
                if (!synapse) continue;

                // go through all six adjacent neighbors for the synapses
                if ((long)(synapse_grid->GridValue(ix - 1, iy, iz) + 0.5) != synapse ||
                    (long)(synapse_grid->GridValue(ix + 1, iy, iz) + 0.5) != synapse ||
                    (long)(synapse_grid->GridValue(ix, iy - 1, iz) + 0.5) != synapse ||
                    (long)(synapse_grid->GridValue(ix, iy + 1, iz) + 0.5) != synapse ||
                    (long)(synapse_grid->GridValue(ix, iy, iz - 1) + 0.5) != synapse ||
                    (long)(synapse_grid->GridValue(ix, iy, iz + 1) + 0.5) != synapse)
                {
                    synapses[synapse].push_back(iv);
                }

                // add to this synapses's list of segments
                long segment = (long) (segmentation_grid->GridValue(ix, iy, iz) + 0.5);
                segments_per_synapse[synapse].insert(segment);
            }
        }
    }

    if (print_verbose) {
        printf("Preprocessing...\n");
        printf("  Time = %0.2f seconds\n", start_time.Elapsed());
        printf("  Maximum Segment = %ld\n", maximum_segmentation);
    }

    delete segmentation_grid;
    delete synapse_grid;
}



////////////////////////////////////////////////////////////////////////
// Drawing utility functions
////////////////////////////////////////////////////////////////////////

static void DrawSegment(int segment_index)
{
    transformation.Push();
    glPointSize(1.0);
    glBegin(GL_POINTS);
    for (unsigned int iv = 0; iv < segmentations[segment_index].size(); ++iv) {
        // faster rendering with downsampling
        if (RNRandomScalar() > 1.0 / downsample_rate) continue;

        // get the coordinates from the linear index
        long ix, iy, iz;
        IndexToIndices(segmentations[segment_index][iv], ix, iy, iz);
        glVertex3f(ix, iy, iz);
    }
    glEnd();

    transformation.Pop();
}



static void DrawSynapse(int synapse_index)
{
    transformation.Push();
    glPointSize(1.0);
    glBegin(GL_POINTS);
    for (unsigned int iv = 0; iv < synapses[synapse_index].size(); ++iv) {
        // faster rendering with downsampling
        if (RNRandomScalar() > 1.0 / downsample_rate) continue;

        // get the coordinates from the linear index
        long ix, iy, iz;
        IndexToIndices(synapses[synapse_index][iv], ix, iy, iz);
        glVertex3f(ix, iy, iz);
    }
    glEnd();

    transformation.Pop();
}


static void DrawConnectome(int segment_index)
{
    // sizes for skeleton joints
    const double joint_size = 4;
    const double endpoint_size = 60;

    transformation.Push();

    glPointSize(joint_size);
    glBegin(GL_POINTS);
    for (unsigned long is = 0; is < wirings[segment_index].size(); ++is) {
        long iv = wirings[segment_index][is];
        bool endpoint = (iv < 0);
        if (endpoint) continue;
        long ix, iy, iz;
        IndexToIndices(iv, ix, iy, iz);
        glVertex3f(ix, iy, iz);
    }
    glEnd();
    glPointSize(1.0);

    transformation.Pop();

    for (unsigned long is = 0; is < wirings[segment_index].size(); ++is) {
        // get the coordinates from the linear index
        long iv = wirings[segment_index][is];
        bool endpoint = (iv < 0);
        if (not endpoint) continue;
        iv = -1 * iv;
        long ix, iy, iz;
        IndexToIndices(iv, ix, iy, iz);

        // convert to world coordinates since transformation is popped
        ix = resolution[OR_X] * ix;
        iy = resolution[OR_Y] * iy;
        iz = resolution[OR_Z] * iz;

        R3Sphere(R3Point(ix, iy, iz), endpoint_size).Draw();
    }
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

    if (synapse_centric) {
        // draw synapse
        RNLoadRgb(Color(synapse_index));
        DrawSynapse(synapse_index);
        for (std::set<long>::iterator it = segments_per_synapse[synapse_index].begin(); it != segments_per_synapse[synapse_index].end(); ++it) {
            long segment = *it;
            RNLoadRgb(Color(segment));
            DrawSegment(segment);
            RNLoadRgb(RNRgb(not background_color[0], not background_color[1], not background_color[2]));
            DrawConnectome(segment);
        }
    }
    else {
        // draw segmentation
        RNLoadRgb(Color(segmentation_index));
        DrawSegment(segmentation_index);
        RNLoadRgb(RNRgb(not background_color[0], not background_color[1], not background_color[2]));
        DrawConnectome(segmentation_index);
        for (std::set<long>::iterator it = synapses_per_segment[segmentation_index].begin(); it != synapses_per_segment[segmentation_index].end(); ++it) {
            long synapse = *it;
            RNLoadRgb(Color(synapse));
            DrawSynapse(synapse);
        }
    }

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
            if (synapse_centric) {
                --synapse_index;    
                if (synapse_index < 0)  
                    synapse_index = 0;

                printf("Synapse %d\n", synapse_index);
            }
            else {
                --segmentation_index;
                if (segmentation_index < 0)
                    segmentation_index = 0;

                printf("Segment %d\n", segmentation_index);                
            }

            break;
        }

        case GLUT_KEY_RIGHT: {
            if (synapse_centric) {
                ++synapse_index;
                if (synapse_index >= maximum_synapse)
                    synapse_index = maximum_synapse - 1;

                printf("Synapse %d\n", synapse_index);
            }
            else {
                ++segmentation_index;
                if (segmentation_index >= maximum_segmentation)
                    segmentation_index = maximum_segmentation - 1;

                printf("Segment %d\n", segmentation_index);
            }

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

        case 'C':
        case 'c': {
            color_cycle = (color_cycle + 1) % ncolor_opts;
            break;
        }

        case 'P':
        case 'p': {
            const R3Camera &camera = viewer->Camera();
            printf("Camera:\n");
            printf("  Origin:  (%lf, %lf, %lf)\n", camera.Origin().X(), camera.Origin().Y(), camera.Origin().Z());
            printf("  Towards: (%lf, %lf, %lf)\n", camera.Towards().X(), camera.Towards().Y(), camera.Towards().Z());
            printf("  Up:      (%lf, %lf, %lf)\n", camera.Up().X(), camera.Up().Y(), camera.Up().Z());
            printf("  XFOV:    %lf\n", camera.XFOV());
            printf("  YFOV:    %lf\n", camera.YFOV());
            printf("  Near:    %lf\n", camera.Near());
            printf("  Far:     %lf\n", camera.Far());
            const R2Viewport &viewport = viewer->Viewport();
            printf("Viewport:\n");
            printf("  XMin:    %d\n", viewport.XMin());
            printf("  YMin:    %d\n", viewport.YMin());
            printf("  Width:   %d\n", viewport.Width());
            printf("  Height:  %d\n", viewport.Height());
            printf("viewer = new R3Viewer(R3Camera(R3Point(%lf, %lf, %lf), R3Vector(%lf, %lf, %lf), R3Vector(%lf, %lf, %lf), %lf, %lf, %lf, %lf), R2Viewport(%d, %d, %d, %d));\n", camera.Origin().X(), camera.Origin().Y(), camera.Origin().Z(), 
                camera.Towards().X(), camera.Towards().Y(), camera.Towards().Z(), camera.Up().X(), camera.Up().Y(), camera.Up().Z(), camera.XFOV(), camera.YFOV(), camera.Near(), camera.Far(), viewport.XMin(), viewport.YMin(), viewport.Width(), viewport.Height());
            break;
        }

        case 'S': 
        case 's': {
            synapse_centric = 1 - synapse_centric;
            break;
        }

        case 'T':
        case 't': {
            stage_index = 1 - stage_index;
            ReadConnectomeData();
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
    GLUTwindow = glutCreateWindow("Wiring Diagram Visualizer");

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
            else if (!strcmp(*argv, "-start_index")) { argv++; argc--; segmentation_index = atoi(*argv); }
            else { fprintf(stderr, "Invalid program argument: %s\n", *argv); return 0; }
        } else {
            if (!prefix) prefix = *argv;
            else { fprintf(stderr, "Invalid program argument: %s\n", *argv); return 0; }
        }
        argv++; argc--;
    }

    // error if there is no input name
    if (!prefix) {
        fprintf(stderr, "Need to supply a prefix for data files\n");
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

    if (!ReadMetaData(prefix)) exit(-1);
    if (!ReadData()) exit(-1);
    if (!ReadConnectomeData()) exit(-1);

    // get all of the preprocessing time
    Preprocessing();
    // set world box
    world_box = R3Box(0, 0, 0, resolution[OR_X] * grid_size[OR_X], resolution[OR_Y] * grid_size[OR_Y], resolution[OR_Z] * grid_size[OR_Z]);

    // get the transformation
    transformation = R3Affine(R4Matrix(resolution[OR_X], 0, 0, 0, 0, resolution[OR_Y], 0, 0, 0, 0, resolution[OR_Z], 0, 0, 0, 0, 1));

    // create viewer
    viewer = CreateViewer();
    if (!viewer) exit(-1);

    // initialize GLUT
    GLUTInit(&argc, argv);

    // to start with different viewer add command here and recompile

    // run GLUT interface
    GLUTMainLoop();

    // return success
    return 1;
}