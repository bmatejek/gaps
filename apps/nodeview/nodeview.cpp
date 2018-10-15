// Source file for visualizing node errors (small errors)



// include files

#include "R3Graphics/R3Graphics.h"
#include "RNDataStructures/RNDataStructures.h"
#include "fglut/fglut.h"
#include <fstream>
#include <vector>
#include <unordered_set>
#include <unordered_map>


// GLUT defines

#define ENTER 13
#define ESCAPE 27
#define SPACEBAR 32

static const int IB_Z = 0;
static const int IB_Y = 1;
static const int IB_X = 2;


// class declarations

struct RNMeta;



// program arguments

// I/O flags
static int print_debug = 0;
static int print_verbose = 0;
static const char* prefix = NULL;
static long threshold = 20000;



// program variables

static double resolution[3] = { -1, -1, -1 };
static long grid_size[3] = { -1, -1, -1 };
static R3Affine transformation = R3null_affine;
static R3Viewer *viewer = NULL;
static R3Box world_box = R3null_box;



// voxel grids

static R3Grid *rhoana_grid = NULL;
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

static std::vector<long> *segmentations = NULL;
static long maximum_segmentation = -1;
static long *nvoxels_per_segment = NULL;
static long *seg2gold_mapping = NULL;

static std::vector<long> segments = std::vector<long>();
static std::vector<long> *neighbors = NULL;

static long segment_index = 0;
static long neighbor_index = 0;





// display variables

static const int ncolor_opts = 6;
static int show_bbox = 1;
static RNScalar downsample_rate = 2.0;
static int color_cycle = 0;
static int show_ground_truth = 1;
static int show_neighbor = 1;


////////////////////////////////////////////////////////////////////////
// Input/output functions
////////////////////////////////////////////////////////////////////////

struct RNMeta {
    // instance variables
    double resolution[3];
    char prefix[4096];
    char gold_filename[4096];
    char gold_dataset[128];
    char image_filename[4096];
    char image_dataset[128];
    char mask_filename[4096];
    char mask_dataset[128];
    char rhoana_filename[4096];
    char rhoana_dataset[128];
    long grid_size[3];
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
    if (fscanf(fp, "%lfx%lfx%lf\n", &(meta_data.resolution[IB_X]), &(meta_data.resolution[IB_Y]), &(meta_data.resolution[IB_Z])) != 3) return 0;

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
    if (fscanf(fp, "%s %s\n", meta_data.rhoana_filename, meta_data.rhoana_dataset) != 2) return 0;

    // skip synapse
    if (!fgets(comment, 4096, fp)) return 0;
    if (!fgets(comment, 4096, fp)) return 0;

    // update the global resolution
    for (int dim = 0; dim <= 2; ++dim)
        resolution[dim] = meta_data.resolution[dim];

    // skip world box
    if (!fgets(comment, 4096, fp)) return 0;
    if (!fgets(comment, 4096, fp)) return 0;

    // read in requisite information
    if (!fgets(comment, 4096, fp)) return 0;
    if (fscanf(fp, "%ldx%ldx%ld\n", &(meta_data.grid_size[IB_X]), &(meta_data.grid_size[IB_Y]), &(meta_data.grid_size[IB_Z])) != 3) return 0;

    meta_data.world_box = R3null_box;
    meta_data.scaled_box = R3null_box;

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

    R3Grid **gold_grids = RNReadH5File(meta_data.gold_filename, meta_data.gold_dataset);
    gold_grid = gold_grids[0];
    delete[] gold_grids;

    R3Grid **rhoana_grids = RNReadH5File(meta_data.rhoana_filename, meta_data.rhoana_dataset);
    rhoana_grid = rhoana_grids[0];
    delete[] rhoana_grids;

    grid_size[IB_X] = rhoana_grid->XResolution();
    grid_size[IB_Y] = rhoana_grid->YResolution();
    grid_size[IB_Z] = rhoana_grid->ZResolution();

    // get the maximum values for each grid
    maximum_segmentation = (long) (rhoana_grid->Maximum() + 0.5) + 1;

    // print statistics
    if(print_verbose) {
        printf("Read voxel grids...\n");
        printf("  Time = %.2f seconds\n", start_time.Elapsed());
        printf("  Grid Size = (%ld %ld %ld)\n", grid_size[IB_X], grid_size[IB_Y], grid_size[IB_Z]);
        printf("  Resolution = (%0.2lf %0.2lf %0.2lf)\n", resolution[IB_X], resolution[IB_Y], resolution[IB_Z]);
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
  iz = index / (grid_size[IB_X] * grid_size[IB_Y]);
  iy = (index - iz * grid_size[IB_X] * grid_size[IB_Y]) / grid_size[IB_X];
  ix = index % grid_size[IB_X];
}



static long IndicesToIndex(long ix, long iy, long iz)
{
    return iz * grid_size[IB_Y] * grid_size[IB_X] + iy * grid_size[IB_X] + ix;
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

static long *Seg2GoldMapping(long *segmentation, long *gold, long nentries, double match_threshold=0.80, double nonzero_threshold=0.40)
{
    // find the maximum segmentation value
    long max_segmentation_value = 0;
    for (long iv = 0; iv < nentries; ++iv) {
        if (segmentation[iv] > max_segmentation_value)
            max_segmentation_value = segmentation[iv];
    }
    max_segmentation_value++;

    // find the maximum gold value
    long max_gold_value = 0;
    for (long iv = 0; iv < nentries; ++iv) {
        if (gold[iv] > max_gold_value) 
            max_gold_value = gold[iv];
    }
    max_gold_value++;

    // find the number of voxels per segment
    nvoxels_per_segment = new long[max_segmentation_value];
    for (long iv = 0; iv < max_segmentation_value; ++iv)
        nvoxels_per_segment[iv] = 0;
    for (long iv = 0; iv < nentries; ++iv)
        nvoxels_per_segment[segmentation[iv]]++;

    std::unordered_map<long, std::unordered_map<long, long> > seg2gold_overlap = std::unordered_map<long, std::unordered_map<long, long> >();   
    for (long is = 0; is < max_segmentation_value; ++is) {
        if (nvoxels_per_segment[is]) seg2gold_overlap.insert(std::pair<long, std::unordered_map<long, long> >(is, std::unordered_map<long, long>()));
    }
    
    for (long iv = 0; iv < nentries; ++iv) {
        seg2gold_overlap[segmentation[iv]][gold[iv]]++;
    }
    
    // create the mapping
    long *segmentation_to_gold = new long[max_segmentation_value];
    for (long is = 0; is < max_segmentation_value; ++is) {
        if (!nvoxels_per_segment[is]) { segmentation_to_gold[is] = 0; continue; }
        long gold_id = 0;
        long gold_max_value = 0;

        // only gets label of 0 if the number of non zero voxels is below threshold
        for (std::unordered_map<long, long>::iterator iter = seg2gold_overlap[is].begin(); iter != seg2gold_overlap[is].end(); ++iter) {
            if (not iter->first) continue;
            if (iter->second > gold_max_value) {
                gold_max_value = iter->second;
                gold_id = iter->first;
            }
        }
        
        // number of non zero pixels must be greater than the nonzero threshold
        if ((double)(nvoxels_per_segment[is] - seg2gold_overlap[is][0]) / nvoxels_per_segment[is] < nonzero_threshold) segmentation_to_gold[is] = 0;
        // the number of matching gold values divided by the number of non zero pixels must be greater than the match threshold or it is a merge error
        else if (gold_max_value / (double)(nvoxels_per_segment[is] - seg2gold_overlap[is][0]) < match_threshold) segmentation_to_gold[is] = -1;
        else segmentation_to_gold[is] = gold_id;
    }

    return segmentation_to_gold;
}



static void Preprocessing(void)
{
    RNTime start_time;
    start_time.Read();

    // create a vector for each valid ID
    segmentations = new std::vector<long>[maximum_segmentation];
    for (long iv = 0; iv < maximum_segmentation; ++iv)
        segmentations[iv] = std::vector<long>();

    // create input for the seg2gold mapping
    long nentries = grid_size[IB_X] * grid_size[IB_Y] * grid_size[IB_Z];
    long *gold = new long[nentries];
    long *segmentation = new long[nentries];
    for (long iv = 0; iv < nentries; ++iv) {
        long segment_index = (long)(rhoana_grid->GridValue(iv) + 0.5);
        long gold_index = (long)(gold_grid->GridValue(iv) + 0.5);

        gold[iv] = gold_index;
        segmentation[iv] = segment_index;
    }

    seg2gold_mapping = Seg2GoldMapping(segmentation, gold, nentries);

    std::unordered_set<long> *neighbor_set = new std::unordered_set<long>[maximum_segmentation];
    for (long is = 0; is < maximum_segmentation; ++is)
        neighbor_set[is] = std::unordered_set<long>();

    // go through all voxels to see if it belongs to the boundary
    for (int iz = 0; iz < grid_size[IB_Z]; ++iz) {
        for (int iy = 0; iy < grid_size[IB_Y]; ++iy) {
            for (int ix = 0; ix < grid_size[IB_X]; ++ix) {
                long segment = (long)(rhoana_grid->GridValue(ix, iy, iz) + 0.5);
                if (!segment) continue;

                bool boundary = false;

                if (ix == 0) boundary = true;
                else {
                    long neighbor = (long) (rhoana_grid->GridValue(ix - 1, iy, iz) + 0.5);
                    if (neighbor and neighbor != segment) { neighbor_set[segment].insert(neighbor); boundary = true; }
                }
                if (iy == 0) boundary = true;
                else {
                    long neighbor = (long) (rhoana_grid->GridValue(ix, iy - 1, iz) + 0.5);
                    if (neighbor and neighbor != segment) { neighbor_set[segment].insert(neighbor); boundary = true; }
                }
                if (iz == 0) boundary = true;
                else {
                    long neighbor = (long) (rhoana_grid->GridValue(ix, iy, iz - 1) + 0.5);
                    if (neighbor and neighbor != segment) { neighbor_set[segment].insert(neighbor); boundary = true; }
                }

                if (ix == grid_size[IB_X] - 1) boundary = true;
                else {
                    long neighbor = (long) (rhoana_grid->GridValue(ix + 1, iy, iz) + 0.5);
                    if (neighbor and neighbor != segment) { neighbor_set[segment].insert(neighbor); boundary = true; }
                }
                if (iy == grid_size[IB_Y] - 1) boundary = true;
                else {
                    long neighbor = (long) (rhoana_grid->GridValue(ix, iy + 1, iz) + 0.5);
                    if (neighbor and neighbor != segment) { neighbor_set[segment].insert(neighbor); boundary = true; }
                }
                if (iz == grid_size[IB_Z] - 1) boundary = true;
                else {
                    long neighbor = (long) (rhoana_grid->GridValue(ix, iy, iz + 1) + 0.5);
                    if (neighbor and neighbor != segment) { neighbor_set[segment].insert(neighbor); boundary = true; }
                }

                if (boundary) segmentations[segment].push_back(IndicesToIndex(ix, iy, iz));
            }
        }
    }

    neighbors = new std::vector<long>[maximum_segmentation];
    for (long is = 0; is < maximum_segmentation; ++is) {
        neighbors[is] = std::vector<long>();

        std::unordered_set<long>::iterator it;
        for (it = neighbor_set[is].begin(); it != neighbor_set[is].end(); ++it) 
            neighbors[is].push_back(*it);
    }

    delete[] gold;
    delete[] segmentation;
    delete[] neighbor_set;

    if (print_verbose) {
        printf("Preprocessing...\n");
        printf("  Time = %0.2f seconds\n", start_time.Elapsed());
        printf("  Maximum Segment = %ld\n", maximum_segmentation);
    }
}



static void AnalyzeSmallSegments(void)
{
    segments.clear();
    for (long is = 0 ; is < maximum_segmentation; ++is) {
        if (nvoxels_per_segment[is] < threshold) segments.push_back(is);
    }

    segment_index = 0;  
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
    if (rhoana_grid) delete rhoana_grid;

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
    
    RNLoadRgb(Color(segments[segment_index]));
    DrawSegment(segments[segment_index]);

    long neighbor = neighbors[segments[segment_index]][neighbor_index];

    if (show_neighbor) {
        if (seg2gold_mapping[neighbor] < 1 or seg2gold_mapping[segments[segment_index]] < 1) RNLoadRgb(RNblue_rgb);
        else if (seg2gold_mapping[segments[segment_index]] == seg2gold_mapping[neighbor]) RNLoadRgb(RNgreen_rgb);
        else RNLoadRgb(RNred_rgb);
        DrawSegment(neighbor);
    }

    // draw neuron data bounding box
    if(show_bbox) {
        RNLoadRgb(RNRgb(not background_color[0], not background_color[1], not background_color[2]));
        world_box.Outline();
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
            do {
                segment_index--;
                if (segment_index < 0) {
                    segment_index = 0;
                    break;
                }
            } while (show_ground_truth and seg2gold_mapping[segments[segment_index]] < 1);
            neighbor_index = 0;
            break;
        }

        case GLUT_KEY_RIGHT: {
            do {
                segment_index++;
                if (segment_index > (long)segments.size() - 1) {
                    segment_index = segments.size() - 1;
                    break;
                }
            } while (show_ground_truth and seg2gold_mapping[segments[segment_index]] < 1);
            neighbor_index = 0;
            break;    
        }

        case GLUT_KEY_UP: {
            neighbor_index++;
            if (neighbor_index > (long) neighbors[segments[segment_index]].size() - 1)
                neighbor_index = neighbors[segments[segment_index]].size() - 1;
            
            break;
        }

        case GLUT_KEY_DOWN: {
            neighbor_index--;
            if (neighbor_index < 0)
                neighbor_index = 0;

            break;
        }

    }

    printf("Showing segment %ld with neighbor %ld\n", segments[segment_index], neighbors[segments[segment_index]][neighbor_index]);

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

        case 'E': 
        case 'e': {
            threshold = threshold - 1000;
            if (threshold < 5000) threshold = 5000;
            AnalyzeSmallSegments();
            printf("Showing nodes smaller than %ld\n", threshold);
            break;
        }

        case 'G': 
        case 'g': {
            show_ground_truth = 1 - show_ground_truth;
            break;
        }

        case 'N':
        case 'n': {
            show_neighbor = 1 - show_neighbor;
            break;
        }

        case 'T': 
        case 't': {
            threshold = threshold + 1000;
            AnalyzeSmallSegments();
            printf("Showing nodes smaller than %ld\n", threshold);
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
            else if (!strcmp(*argv, "-threshold")) { argv++; argc--; threshold = atoi(*argv); }
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
    
    // get all of the preprocessing time
    Preprocessing();
    AnalyzeSmallSegments();


    // set world box
    world_box = R3Box(0, 0, 0, resolution[IB_X] * grid_size[IB_X], resolution[IB_Y] * grid_size[IB_Y], resolution[IB_Z] * grid_size[IB_Z]);

    // get the transformation
    transformation = R3Affine(R4Matrix(resolution[IB_X], 0, 0, 0, 0, resolution[IB_Y], 0, 0, 0, 0, resolution[IB_Z], 0, 0, 0, 0, 1));

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