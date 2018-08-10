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
#define DELETE 127

// class declarations

struct RNMeta;
struct SWCEntry;
struct MergeCandidate;



// program arguments

// I/O flags
static int print_debug = 0;
static int print_verbose = 0;
// dataset to examine
static const char* prefix = NULL;
static char output_directory[4096];



// program variables

static RNScalar resolution[3] = { 6, 6, 30 };
static int grid_size[3] = { -1, -1, -1 };
static R3Affine transformation = R3null_affine;
static R3Viewer *viewer = NULL;
static R3Box world_box = R3null_box;



// voxel grids

static R3Grid *gold_grid = NULL;



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

static std::vector<long> *golds = NULL;
static long maximum_gold = -1;
static std::vector<R3Point> skeleton_endpoints = std::vector<R3Point>();



// display variables

static int show_bbox = 1;
static int example_index = 0;
static int gold_index = 0;
static const int cutoff = 500;
static int largest_segments[cutoff];
static RNScalar downsample_rate = 2.0;




////////////////////////////////////////////////////////////////////////
// Input/output functions
////////////////////////////////////////////////////////////////////////

struct RNMeta {
    // instance variables
    int resolution[3];
    char prefix[4096];
    char gold_filename[4096];
    char gold_dataset[128];
    char image_filename[4096];
    char image_dataset[128];
    char mask_filename[4096];
    char mask_dataset[128];
    char rhoana_filename[4096];
    char rhoana_dataset[128];
    R3Box world_box;
    R3Box scaled_box;
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
    if (fscanf(fp, "%dx%dx%d\n", &(meta_data.resolution[RN_X]), &(meta_data.resolution[RN_Y]), &(meta_data.resolution[RN_Z])) != 3) return 0;
    
    // skip affinities
    if (!fgets(comment, 4096, fp)) return 0;
    if (!fgets(comment, 4096, fp)) return 0;

    // skip boundaries
    if (!fgets(comment, 4096, fp)) return 0;
    if (!fgets(comment, 4096, fp)) return 0;

    // read in gold information
    if (!fgets(comment, 4096, fp)) return 0;
    if (fscanf(fp, "%s %s\n", meta_data.gold_filename, meta_data.gold_dataset) != 2) return 0;

    // read image
    if (!fgets(comment, 4096, fp)) return 0;
    if (!fgets(comment, 4096, fp)) return 0;

    // skip mask
    if (!fgets(comment, 4096, fp)) return 0;
    if (!fgets(comment, 4096, fp)) return 0;

    // read segmentation
    if (!fgets(comment, 4096, fp)) return 0;
    if (!fgets(comment, 4096, fp)) return 0;

    // skip synapse
    if (!fgets(comment, 4096, fp)) return 0;
    if (!fgets(comment, 4096, fp)) return 0;

    // update the global resolution
    for (int dim = 0; dim <= 2; ++dim)
        resolution[dim] = meta_data.resolution[dim];

    // skip world box
    if (!fgets(comment, 4096, fp)) return 0;
    if (!fgets(comment, 4096, fp)) return 0;
    meta_data.world_box = R3null_box;
    meta_data.scaled_box = R3null_box;

    // close the file
    fclose(fp);

    // return success
    return 1;
}



struct SWCEntry {
    // constructors
    SWCEntry(int sample_number, int structure_identifier, RNScalar x_position, RNScalar y_position, RNScalar z_position, RNScalar radius, int parent_sample) : 
    sample_number(sample_number), 
    structure_identifier(structure_identifier), 
    x_position(x_position), 
    y_position(y_position), 
    z_position(z_position), 
    radius(radius), 
    parent_sample(parent_sample)
    {
    }

    // access variables
    RNScalar X(void) { return x_position; }
    RNScalar Y(void) { return y_position; }
    RNScalar Z(void) { return z_position; }
    R3Point P(void) { return R3Point(x_position, y_position, z_position); }

    // instance variables
    int sample_number;
    int structure_identifier;
    RNScalar x_position;
    RNScalar y_position;
    RNScalar z_position;
    RNScalar radius;
    int parent_sample;
};



static int ReadData(void)
{
    // start statistics
    RNTime start_time;
    start_time.Read();

    R3Grid **gold_grids = RNReadH5File(meta_data.gold_filename, meta_data.gold_dataset);
    gold_grid = gold_grids[0];
    delete[] gold_grids;

    grid_size[RN_X] = gold_grid->XResolution();
    grid_size[RN_Y] = gold_grid->YResolution();
    grid_size[RN_Z] = gold_grid->ZResolution();

    // print statistics
    if(print_verbose) {
        printf("Read voxel grids...\n");
        printf("  Time = %.2f seconds\n", start_time.Elapsed());
        printf("  Grid Size = (%d %d %d)\n", gold_grid->XResolution(), gold_grid->YResolution(), gold_grid->ZResolution());
        printf("  Resolution = (%0.2lf %0.2lf %0.2lf)\n", resolution[RN_X], resolution[RN_Y], resolution[RN_Z]);
    }

    // return success
    return 1;
}



static int ReadLargestSegments(void)
{
    char input_filename[4096];
    sprintf(input_filename, "benchmarks/skeleton/%s-skeleton-benchmark-examples.bin", prefix);

    FILE *fp = fopen(input_filename, "rb");
    
    long input_cutoff;
    if (fread(&input_cutoff, sizeof(long), 1, fp) != 1) return 0;
    assert (cutoff == input_cutoff);
    for (long iv = 0; iv < cutoff; ++iv) {
        if (fread(&(largest_segments[iv]), sizeof(long), 1, fp) != 1) return 0;
    }

    fclose(fp);

    printf("Read %d examples...\n", cutoff);

    return 1;
}


////////////////////////////////////////////////////////////////////////
// Helper functions
////////////////////////////////////////////////////////////////////////

static void IndexToIndices(long index, long& ix, long& iy, long& iz)
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



static RNRgb Color(unsigned long value)
{
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

    // get the maximum values for each grid
    maximum_gold = (long) (gold_grid->Maximum() + 0.5) + 1;

    // create a vector for each valid ID
    golds = new std::vector<long>[maximum_gold];
    for (long iv = 0; iv < maximum_gold; ++iv) 
        golds[iv] = std::vector<long>();

    // go through all voxels to see if it belongs to the boundary
    for (int iz = 1; iz < grid_size[RN_Z] - 1; ++iz) {
        for (int iy = 1; iy < grid_size[RN_Y] - 1; ++iy) {
            for (int ix = 1; ix < grid_size[RN_X] - 1; ++ix) {
                long iv = IndicesToIndex(ix, iy, iz);

                long gold = (long) (gold_grid->GridValue(ix, iy, iz) + 0.5);

                if ((long)(gold_grid->GridValue(ix - 1, iy, iz) + 0.5) != gold ||
                    (long)(gold_grid->GridValue(ix + 1, iy, iz) + 0.5) != gold ||
                    (long)(gold_grid->GridValue(ix, iy - 1, iz) + 0.5) != gold ||
                    (long)(gold_grid->GridValue(ix, iy + 1, iz) + 0.5) != gold ||
                    (long)(gold_grid->GridValue(ix, iy, iz - 1) + 0.5) != gold ||
                    (long)(gold_grid->GridValue(ix, iy, iz + 1) + 0.5) != gold)
                {
                    golds[gold].push_back(iv);
                }
            }
        }
    }

    if (print_verbose) {
        printf("Preprocessing...\n");
        printf("  Time = %0.2f seconds\n", start_time.Elapsed());
        printf("  Maximum Gold = %ld\n", maximum_gold);
    }
}



////////////////////////////////////////////////////////////////////////
// Drawing utility functions
////////////////////////////////////////////////////////////////////////

static void DrawEndpoints()
{
    glPointSize(10.0);
    glBegin(GL_POINTS);
    RNLoadRgb(RNRgb(not background_color[0], not background_color[1], not background_color[2]));
    for (unsigned long ie = 0; ie < skeleton_endpoints.size(); ++ie) {
        R3Point endpoint = skeleton_endpoints[ie];
        glVertex3f(endpoint.X(), endpoint.Y(), endpoint.Z());
    }
    glEnd();
    glPointSize(1.0);

}



static void DrawGold(int segment_index)
{
    glBegin(GL_POINTS);
    for (unsigned int iv = 0; iv < golds[segment_index].size(); ++iv) {
        // faster rendering with downsampling
        if (RNRandomScalar() > 1.0 / downsample_rate) continue;

        // get the coordinates from the linear index
        long ix, iy, iz;
        IndexToIndices(golds[segment_index][iv], ix, iy, iz);
        glVertex3f(ix, iy, iz);
    }
    glEnd();
}



static void DrawPointClouds(void)
{
    // push the transformation
    transformation.Push();

    RNLoadRgb(Color(gold_index));
    DrawGold(gold_index);
    DrawEndpoints();

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

    if (gold_grid) delete gold_grid;

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



void ReadSkeletonEndpoints(void)
{
    char input_filename[4096];
    sprintf(input_filename, "%s/example-%d.pts", output_directory, example_index);

    FILE *fp = fopen(input_filename, "rb"); 
    if (fp) {
        long nendpoints;
        if (fread(&nendpoints, sizeof(long), 1, fp) != 1) fprintf(stderr, "Failed to read %s\n", input_filename);
        for (long ie = 0; ie < nendpoints; ++ie) {
            long zpoint, ypoint, xpoint;
            // read the three point locations
            if (fread(&zpoint, sizeof(long), 1, fp) != 1) fprintf(stderr, "Failed to read %s\n", input_filename);
            if (fread(&ypoint, sizeof(long), 1, fp) != 1) fprintf(stderr, "Failed to read %s\n", input_filename);
            if (fread(&xpoint, sizeof(long), 1, fp) != 1) fprintf(stderr, "Failed to read %s\n", input_filename);

            skeleton_endpoints.push_back(R3Point(xpoint, ypoint, zpoint));
        }
        fclose(fp);
    }   
}



void WriteSkeletonEndpoints(void)
{
    if (not skeleton_endpoints.size()) return;

    char output_filename[4096];
    sprintf(output_filename, "%s/example-%d.pts", output_directory, example_index);

    FILE *fp = fopen(output_filename, "wb");
    if (!fp) { fprintf(stderr, "Failed to write %s\n", output_filename); }

    long nendpoints = skeleton_endpoints.size();
    fwrite(&nendpoints, sizeof(long), 1, fp);
    for (long ie = 0; ie < nendpoints; ++ie) {
        R3Point endpoint = skeleton_endpoints[ie];
        long xpoint = (long) (endpoint.X() + 0.5);
        long ypoint = (long) (endpoint.Y() + 0.5);
        long zpoint = (long) (endpoint.Z() + 0.5);
        fwrite(&zpoint, sizeof(long), 1, fp);
        fwrite(&ypoint, sizeof(long), 1, fp);
        fwrite(&xpoint, sizeof(long), 1, fp);
    }

    fclose(fp);
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

    WriteSkeletonEndpoints();

    switch(key) {
        case GLUT_KEY_LEFT: {
            --example_index;
            if (example_index < 0)
                example_index = 0;

            break;
        }

        case GLUT_KEY_RIGHT: {
                ++example_index;
                if (example_index >= cutoff)
                    example_index = cutoff - 1;

            break;    
        }
    }

    // upload the new endpoints
    skeleton_endpoints.clear();
    ReadSkeletonEndpoints();

    gold_index = largest_segments[example_index];
    printf("Label: %d\n", gold_index);

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

        case 'E':
        case 'e': {
            R3Ray world_ray = viewer->WorldRay(x, y);

            RNLength closest_distance = DBL_MAX;
            R3Point closest_point;

            static const int buffer = 40;

            for (unsigned long ig = 0; ig < golds[gold_index].size(); ++ig) {
                long iv = golds[gold_index][ig];

                // convert to three coordinates
                long ix, iy, iz;
                IndexToIndices(iv, ix, iy, iz);

                R3Point location = R3Point(ix, iy, iz);
                location.Transform(transformation);

                RNLength distance = R3Distance(world_ray, location);

                if (distance < closest_distance) {
                    closest_distance = distance;
                    closest_point = R3Point(ix, iy, iz);
                }
            }

            // find the next closest point farther down the ray
            closest_distance = DBL_MAX;
            R3Point second_closest_point;

            for (unsigned long ig = 0 ; ig < golds[gold_index].size(); ++ig) {
                long iv = golds[gold_index][ig];

                long ix, iy, iz;
                IndexToIndices(iv, ix, iy, iz);

                R3Point location = R3Point(ix, iy, iz);
                location.Transform(transformation);

                RNLength distance = R3Distance(world_ray, location);
                
                long deltaz = (long) (closest_point.Z() - iz + 0.5);
                long deltay = (long) (closest_point.Y() - iy + 0.5);
                long deltax = (long) (closest_point.X() - ix + 0.5);

                RNScalar normalized_distance = sqrt(resolution[RN_Z] * resolution[RN_Z] * deltaz * deltaz + resolution[RN_Y] * resolution[RN_Y] * deltay * deltay + resolution[RN_X] * resolution[RN_X] * deltax * deltax);
                if (normalized_distance < buffer) continue;

                if (distance < closest_distance) {
                    closest_distance = distance;
                    second_closest_point = R3Point(ix, iy, iz);
                }
            }

            skeleton_endpoints.push_back((closest_point + second_closest_point) / 2);

            printf("  No. endpoints: %ld\n", skeleton_endpoints.size());

            break;
        }

        // CTRL + Z
        case 26: {
            if (not skeleton_endpoints.size()) break;
            skeleton_endpoints.pop_back();
            break;
        }

        case DELETE: {
            R3Ray world_ray = viewer->WorldRay(x, y);

            RNLength closest_distance = DBL_MAX;
            unsigned long delete_index = 0;

            for (unsigned long iv = 0; iv < skeleton_endpoints.size(); ++iv) {
                R3Point location = skeleton_endpoints[iv];
                location.Transform(transformation);

                RNLength distance = R3Distance(world_ray, location);
                if (distance < closest_distance) {
                    closest_distance = distance;
                    delete_index = iv;
                }
            }

            skeleton_endpoints.erase(skeleton_endpoints.begin() + delete_index);

            printf("  No. endpoints: %ld\n", skeleton_endpoints.size());

            break;

        }

        case ENTER: {
            background_color[0] = 1.0 - background_color[0];
            background_color[1] = 1.0 - background_color[1];
            background_color[2] = 1.0 - background_color[2];
            break;
        }

        case ESCAPE: {
            WriteSkeletonEndpoints();

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
            if(!prefix) prefix = *argv;
            else { fprintf(stderr, "Invalid program argument: %s\n", *argv); return 0; }
        }
        argv++; argc--;
    }

    // error if there is no input name
    if(!prefix) {
        fprintf(stderr, "Need to supply a prefix for data files\n");
        return 0;
    }
    //if (training and validation) { fprintf(stderr, "Need to choose either training or validation (or neither), not both\n"); return 0; }

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

    // get the output directory for the benchmark information
    sprintf(output_directory, "benchmarks/skeleton/%s/", prefix);

    /////////////////////////////////
    //// Read in the voxel files ////
    /////////////////////////////////

    if (!ReadMetaData(prefix)) exit(-1);
    if (!ReadData()) exit(-1);
    if (!ReadLargestSegments()) exit(-1);

    gold_index = largest_segments[example_index];

    // read the initial skeleton endpoints
    ReadSkeletonEndpoints();

    // get all of the preprocessing time
    Preprocessing();

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