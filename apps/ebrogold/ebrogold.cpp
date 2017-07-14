// Source file for generating ground truth for EbroNet


// include files

#include "R3Graphics/R3Graphics.h"
#include "RNDataStructures/RNDataStructures.h"
#include "fglut/fglut.h"
#include <vector>
#include <map>
#include <set>



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
    R3Box scaled_box;
};



struct MergeCandidate {
    MergeCandidate(unsigned long label_one, unsigned long label_two, unsigned long index_one, unsigned long index_two) :
    label_one(label_one),
    label_two(label_two),
    index_one(index_one),
    index_two(index_two)
    {}

    bool operator<(const MergeCandidate &other) const
    {
        if (this->label_one < other.label_one) return true;
        else if (this->label_one > other.label_one) return false;
        else if (this->label_two < other.label_two) return true;
        else if (this->label_two > other.label_two) return false;
        else return false;
    }

    unsigned long label_one;
    unsigned long label_two;
    unsigned long index_one;
    unsigned long index_two;
};




// program arguments

static int print_verbose = 0;
static const char *prefixes[ngrids] = { NULL, NULL };



// program variables

static int resolution[3] = { -1, -1, -1 };
static int grid_size[3] = { -1, -1, -1 };
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



// mapping and random access variables

static std::map<unsigned long, unsigned long> label_to_index[ngrids];
static unsigned long *index_to_label[ngrids] = { NULL, NULL };
static std::vector<unsigned long> *segmentations[ngrids] = { NULL, NULL };
static std::vector<MergeCandidate> candidates;



// display variables

static int show_bbox = 1;
static RNScalar downsample_rate = 6.0;
static int candidate_index = 0;



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

    // update the global resolution
    if (index == 0) {
        for (int dim = 0; dim <= 2; ++dim)
            resolution[dim] = meta_data[index].resolution[dim];
    }
    else {
        for (int dim = 0; dim <= 2; ++dim)
            rn_assertion(resolution[dim] == meta_data[index].resolution[dim]);
    }

    // get the world box
    RNScalar xmin, ymin, zmin, xmax, ymax, zmax;
    if (!fgets(comment, 4096, fp)) return 0;
    if (fscanf(fp, "(%lf,%lf,%lf)-(%lf,%lf,%lf)\n", &xmin, &ymin, &zmin, &xmax, &ymax, &zmax) != 6) return 0;
    meta_data[index].world_box = R3Box(xmin, ymin, zmin, xmax, ymax, zmax);
    meta_data[index].scaled_box = R3Box(xmin * resolution[RN_X], ymin * resolution[RN_Y], zmin * resolution[RN_Z], xmax * resolution[RN_X], ymax * resolution[RN_Y], zmax * resolution[RN_Z]);

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

    // set global grid size
    if (index == 0) {
        for (int dim = 0; dim <= 2; ++dim) 
            grid_size[dim] = rhoana_grids[index]->Resolution(dim);
    }
    else {
        for (int dim = 0; dim <= 2; ++dim)
            rn_assertion(grid_size[dim] == rhoana_grids[index]->Resolution(dim));
    }

    // print information
    if (print_verbose) {
        printf("Read h5 files in %0.2f seconds for %s:\n", start_time.Elapsed(), meta_data[index].prefix);
        printf("  Dimensions: (%d, %d, %d)\n", grid_size[RN_X], grid_size[RN_Y], grid_size[RN_Z]);
        printf("  Resolution: (%d, %d, %d)\n", resolution[RN_X], resolution[RN_Y], resolution[RN_Z]);
    }


    // return success 
    return 1;
}




////////////////////////////////////////////////////////////////////////
// Preprocessing functions
////////////////////////////////////////////////////////////////////////

static void LabelToIndexMapping(int grid_index)
{
    // start statistics
    RNTime start_time;
    start_time.Read();

    // temporary variable
    R3Grid *grid = rhoana_grids[grid_index];

    // find which labels are present
    unsigned long maximum_segmentation = (unsigned long)(grid->Maximum() + 0.5) + 1;
    RNBoolean *present_labels = new RNBoolean[maximum_segmentation];
    for (unsigned long iv = 0; iv < maximum_segmentation; ++iv) 
        present_labels[iv] = FALSE;

    // iterate over the entrie volume to find present labels
    for (long iv = 0; iv < grid->NEntries(); ++iv) {
        unsigned long label = (unsigned long)(grid->GridValue(iv) + 0.5);
        present_labels[label] = TRUE;
    }

    // crate the mapping from segment labels to indices
    label_to_index[grid_index] = std::map<unsigned long, unsigned long>();
    unsigned long nunique_labels = 1; /* 1 indexed for this to work */
    for (unsigned long iv = 1; iv < maximum_segmentation; ++iv) {
        if (present_labels[iv] && !label_to_index[grid_index][iv]) {
            label_to_index[grid_index][iv] = nunique_labels;
            nunique_labels++;
        }
    }

    // create the mapping from indices to labels
    index_to_label[grid_index] = new unsigned long[nunique_labels];
    nunique_labels = 1;
    for (unsigned long iv = 1; iv < maximum_segmentation; ++iv) {
        if (present_labels[iv]) {
            index_to_label[grid_index][nunique_labels] = iv;
            nunique_labels++;
        }
    }

    // free memory
    delete[] present_labels;

    // allocate memory for the segmentation 
    segmentations[grid_index] = new std::vector<unsigned long>[nunique_labels];
    for (unsigned long iv = 0; iv < nunique_labels; ++iv) {
        segmentations[grid_index][iv] = std::vector<unsigned long>();
    }

    // iterate over the entire volume
    int iv = 0;
    for (int iz = 0; iz < grid_size[RN_Z]; ++iz) {
        for (int iy = 0; iy < grid_size[RN_Y]; ++iy) {
            for (int ix = 0; ix < grid_size[RN_X]; ++ix) {
                unsigned long label = (unsigned long)(grid->GridValue(ix, iy, iz) + 0.5);
                if (!label) continue;

                // is this pixel boundary
                RNBoolean boundary = FALSE;
                if (ix > 0 && label != (unsigned long)(grid->GridValue(ix - 1, iy, iz) + 0.5)) boundary = TRUE;
                if (iy > 0 && label != (unsigned long)(grid->GridValue(ix, iy - 1, iz) + 0.5)) boundary = TRUE;
                if (iz > 0 && label != (unsigned long)(grid->GridValue(ix, iy, iz - 1) + 0.5)) boundary = TRUE;
                if (ix < grid_size[RN_X] - 1 && label != (unsigned long)(grid->GridValue(ix + 1, iy, iz) + 0.5)) boundary = TRUE;
                if (iy < grid_size[RN_Y] - 1 && label != (unsigned long)(grid->GridValue(ix, iy + 1, iz) + 0.5)) boundary = TRUE;
                if (iz < grid_size[RN_Z] - 1 && label != (unsigned long)(grid->GridValue(ix, iy, iz + 1) + 0.5)) boundary = TRUE;

                // get this index
                unsigned long mapped_index = label_to_index[grid_index][label];
                rn_assertion((0 < mapped_index) && (mapped_index < nunique_labels));

                // add to the vector
                if (boundary) segmentations[grid_index][mapped_index].push_back(iv);
            }
        }
    }

    // print statistics
    if (print_verbose) {
        printf("Created random access variables for %s in %0.2f seconds\n", meta_data[grid_index].prefix, start_time.Elapsed());
    }
}



static R3Point WorldToGrid(R3Point world_position, int index)
{
    // get the world box for this grid
    R3Box world_box = meta_data[index].world_box;

    // get the distance from the x, y, and z minima
    unsigned long xdiff = (unsigned long)(world_position.X() - world_box.XMin() + 0.5);
    unsigned long ydiff = (unsigned long)(world_position.Y() - world_box.YMin() + 0.5);
    unsigned long zdiff = (unsigned long)(world_position.Z() - world_box.ZMin() + 0.5);

    return R3Point(xdiff, ydiff, zdiff);
}



static R3Point GridToWorld(R3Point grid_position, int index)
{
    // get the world box for this grid
    R3Box world_box = meta_data[index].world_box;

    return R3Point(world_box.XMin() + grid_position.X(), world_box.YMin() + grid_position.Y(), world_box.ZMin() + grid_position.Z());
}



static void FindOverlapCandidates(void)
{
    // start statistics
    RNTime start_time;
    start_time.Read();

    // find all candidates that overlap between the two volumes
    R3Box intersection = meta_data[0].world_box;
    for (int iv = 0; iv < ngrids; ++iv) {
        intersection.Intersect(meta_data[iv].world_box);
    }

    unsigned long xmin = (unsigned long)(intersection.XMin() + 0.5);
    unsigned long ymin = (unsigned long)(intersection.YMin() + 0.5);
    unsigned long zmin = (unsigned long)(intersection.ZMin() + 0.5);
    unsigned long xmax = (unsigned long)(intersection.XMax() + 0.5);
    unsigned long ymax = (unsigned long)(intersection.YMax() + 0.5);
    unsigned long zmax = (unsigned long)(intersection.ZMax() + 0.5);

    // create a unique set of overlap candidates
    std::set<struct MergeCandidate> set = std::set<struct MergeCandidate>();
    candidates = std::vector<struct MergeCandidate>();
    for (unsigned long iz = zmin; iz <= zmax; ++iz) {
        for (unsigned long iy = ymin; iy <= ymax; ++iy) {
            for (unsigned long ix = xmin; ix <= xmax; ++ix) {
                // get the grid value for both candidates
                R3Point position_one = WorldToGrid(R3Point(ix, iy, iz), 0);
                R3Point position_two = WorldToGrid(R3Point(ix, iy, iz), 0);

                // get labels for both grids
                unsigned long label_one = (unsigned long)(rhoana_grids[0]->GridValue(position_one) + 0.5);
                unsigned long label_two = (unsigned long)(rhoana_grids[1]->GridValue(position_two) + 0.5);
                if (!label_one) continue;
                if (!label_two) continue;

                // get the array indices
                unsigned long index_one = label_to_index[0][label_one];
                unsigned long index_two = label_to_index[1][label_two];

                struct MergeCandidate candidate = MergeCandidate(label_one, label_two, index_one, index_two);

                if (set.find(candidate) == set.end()) {
                    set.insert(candidate);
                    candidates.push_back(candidate);
                }
            }
        }
    }

    // print statistics
    if (print_verbose) {
        printf("Calculated all examples in %0.2f seconds:\n", start_time.Elapsed());
        printf(  "Pairs found: %lu\n", candidates.size());
    }
}



////////////////////////////////////////////////////////////////////////
// Drawing utility functions
////////////////////////////////////////////////////////////////////////

static void 
IndexToIndices(int index, int& ix, int& iy, int& iz)
{
  // Set indices of grid value at index
  iz = index / (grid_size[RN_X] * grid_size[RN_Y]);
  iy = (index - iz * grid_size[RN_X] * grid_size[RN_Y]) / grid_size[RN_X];
  ix = index % grid_size[RN_X];
}



static void DrawIndividualSegment(unsigned long index, int grid_index)
{
    glBegin(GL_POINTS);
    for (unsigned int iv = 0; iv < segmentations[grid_index][index].size(); ++iv) {
        if (RNRandomScalar() > 1.0 / downsample_rate) continue;
        int ix, iy, iz;
        IndexToIndices(segmentations[grid_index][index][iv], ix, iy, iz);

        // convert to world coordinates
        R3Point world_position = GridToWorld(R3Point(ix, iy, iz), grid_index);
        printf("%d: %lf %lf %lf\n", iv, world_position.X(), world_position.Y(), world_position.Z());
        glVertex3f(world_position.X(), world_position.Y(), world_position.Z());
    }
    glEnd();
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

    for (int iv = 0; iv < ngrids; ++iv) {
        delete boundary_grids[iv];
        delete image_grids[iv];
        delete rhoana_grids[iv];
    }

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

    // push the transformation
    transformation.Push();

    // draw this candidate
    unsigned long index_one = candidates[candidate_index].index_one;
    RNLoadRgb(RNred_rgb);
    DrawIndividualSegment(index_one, 0);
    unsigned long index_two = candidates[candidate_index].index_two;
    RNLoadRgb(RNblue_rgb);
    DrawIndividualSegment(index_two, 1);

    // pop the transformation
    transformation.Pop();

    // draw large bounding box bounding box
    if(show_bbox) {
        RNLoadRgb(RNwhite_rgb);
        world_box.Outline();

        // draw all bounding boxes
        for (int iv = 0; iv < ngrids; ++iv) {
            meta_data[iv].scaled_box.Outline();
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
            --candidate_index;
            if (candidate_index < 0) candidate_index = 0;
            break;
        }

        case GLUT_KEY_RIGHT: {
            ++candidate_index;
            if (candidate_index >= candidates.size()) candidate_index = candidates.size();
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
    GLUTwindow = glutCreateWindow("Generate Gold Data for EbroNet");

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

static R3Viewer *
CreateViewer(void)
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
        LabelToIndexMapping(iv);
    }

    // find the potential merge candidates
    FindOverlapCandidates();

    // set world box
    world_box = meta_data[0].scaled_box;
    for (int iv = 1; iv < ngrids; ++iv) 
        world_box.Union(meta_data[iv].scaled_box);
    
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