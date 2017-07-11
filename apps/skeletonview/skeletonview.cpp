// Source file for skeleton visualizer



// include files

#include "R3Graphics/R3Graphics.h"
#include "RNDataStructures/RNDataStructures.h"
#include "fglut/fglut.h"
#include <fstream>
#include <string>
#include <vector>
#include <map>



// GLUT defines

#define ENTER 13
#define ESCAPE 27
#define SPACEBAR 32



// program arguments

static int print_debug = 0;
static int print_verbose = 0;
// maximum distance in nanometers
static int maximum_distance = 600;
static const char* prefix = NULL;



// program variables

static RNScalar resolution[3] = { 6, 6, 30 };
static int grid_size[3] = { -1, -1, -1 };
static R3Affine transformation = R3null_affine;
static R3Viewer* viewer = NULL;
static R3Box world_box;



// voxel grids

static R3Grid* grid = NULL;



// GLUT variables

static int GLUTwindow = 0;
static int GLUTwindow_height = 1200;
static int GLUTwindow_width = 1200;
static int GLUTmouse[2] = { 0, 0 };
static int GLUTbutton[3] = { 0, 0, 0 };
static int GLUTmodifiers = 0;



// color arrays

static RNScalar background_color[] = { 0, 0, 0 };



// struct name definitions

struct SWCEntry;
struct MergeCandidate;



// mapping and random access variables

static std::map<unsigned long, unsigned long> label_to_index;
static unsigned long *index_to_label = NULL;
static std::vector<unsigned long> *segmentations = NULL;
static std::vector<SWCEntry> *skeletons = NULL;
static std::vector<R3Point> *skeleton_endpoints = NULL;
static MergeCandidate *candidates = NULL;



// display variables

static int show_bbox = 1;
static int show_skeleton = 1;
static int show_merge_candidate = 1;
static unsigned int segmentation_index = 1;
static unsigned int candidate_index = 0;
static unsigned int ncandidates;
static RNScalar downsample_rate = 25.0;



////////////////////////////////////////////////////////////////////////
// Input/output functions
////////////////////////////////////////////////////////////////////////

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



static int ReadSWCFile(int index)
{
    char input_filename[4096];
    sprintf(input_filename, "skeletons/%s/tree_%lu.swc", prefix, index_to_label[index]);

    std::ifstream fd(input_filename);
    if(!fd.is_open()) {
        if (print_debug) fprintf(stderr, "Failed to read %s\n", input_filename);
        return 0;
    }

    std::string line;
    RNBoolean first_iteration = TRUE;
    while(std::getline(fd, line)) {
        if(first_iteration) {
            first_iteration = FALSE;
            continue;
        }

        int sample_number, structure_identifier, parent_sample;
        RNScalar x_position, y_position, z_position, radius;
        sscanf(line.c_str(), "%d %d %lf %lf %lf %lf %d", &sample_number, &structure_identifier, &x_position,
            &y_position, &z_position, &radius, &parent_sample);

        skeletons[index].push_back(SWCEntry(sample_number, structure_identifier, x_position, y_position, z_position, radius, parent_sample));
    }

    fd.close();

    // determine which points are endpoints
    RNBoolean* endpoints = new RNBoolean[skeletons[index].size()];
    for(unsigned int ie = 0; ie < skeletons[index].size(); ++ie) {
        endpoints[ie] = TRUE;
    }

    // not an endpoint if no children claim it as a parent
    for(unsigned int ie = 0; ie < skeletons[index].size(); ++ie) {
        int parent_sample = skeletons[index][ie].parent_sample;
        if(parent_sample == -1) continue;
        endpoints[parent_sample - 1] = FALSE;
    }

    for(unsigned int ie = 0; ie < skeletons[index].size(); ++ie) {
        // if it has no parent it is an endpoint
        int parent_sample = skeletons[index][ie].parent_sample;
        if(parent_sample == -1) endpoints[ie] = TRUE;
    }

    // add in all of the endpoints
    for(unsigned int ie = 0; ie < skeletons[index].size(); ++ie) {
        if(endpoints[ie]) {
            skeleton_endpoints[index].push_back(skeletons[index][ie].P());
        }
    }

    // free memory
    delete[] endpoints;

    // return success
    return 1;
}



struct MergeCandidate {
    // constructors
    MergeCandidate(void) :
    label_one(0),
    label_two(0),
    x(0),
    y(0),
    z(0),
    ground_truth(FALSE)
    {
    }

    MergeCandidate(unsigned long label_one, unsigned long label_two, unsigned long x, unsigned long y, unsigned long z, RNBoolean ground_truth) :
    label_one(label_one),
    label_two(label_two),
    x(x),
    y(y),
    z(z),
    ground_truth(ground_truth)
    {
    }

    // access variables
    unsigned long LabelOne(void) { return label_one; }
    unsigned long LabelTwo(void) { return label_two; }
    unsigned long X(void) { return x; }
    unsigned long Y(void) { return y; }
    unsigned long Z(void) { return z; }
    unsigned long GroundTruth(void) { return ground_truth; }

    // instance variables
    unsigned long label_one;   
    unsigned long label_two;
    unsigned long x;
    unsigned long y;
    unsigned long z;
    unsigned long ground_truth;
};



static int ReadData(void)
{
    // start statistics
    RNTime start_time;
    start_time.Read();

    // read in the grid
    char grid_filename[4096];
    sprintf(grid_filename, "rhoana/%s_rhoana.h5", prefix);

    R3Grid **grids = RNReadH5File(grid_filename, "main");
    grid = grids[0];
    delete[] grids;

    grid_size[RN_X] = grid->XResolution();
    grid_size[RN_Y] = grid->YResolution();
    grid_size[RN_Z] = grid->ZResolution();

    // print statistics
    if(print_verbose) {
        printf("Read voxel grid...\n");
        printf("  Time = %.2f seconds\n", start_time.Elapsed());
        printf("  Grid Size = (%d %d %d)\n", grid->XResolution(), grid->YResolution(), grid->ZResolution());
        printf("  Resolution = (%0.2lf %0.2lf %0.2lf)\n", resolution[RN_X], resolution[RN_Y], resolution[RN_Z]);
    }

    // return success
    return 1;
}



static int ReadMergeCandidates(void)
{
    // get the candidate filename
    char candidate_filename[4096];
    sprintf(candidate_filename, "skeletons/candidates/%s-%dnm_forward.candidates", prefix, maximum_distance);

    // open the file
    FILE *fp = fopen(candidate_filename, "rb");
    if (!fp) { fprintf(stderr, "Failed to read %s\n", candidate_filename); return 0; }

    if (fread(&ncandidates, sizeof(unsigned int), 1, fp) != 1) return 0;

    // read in all of the candidates
    candidates = new MergeCandidate[ncandidates];
    for (unsigned int iv = 0; iv < ncandidates; ++iv) {
        unsigned long label_one;
        unsigned long label_two;
        unsigned long x;
        unsigned long y;
        unsigned long z;
        unsigned long ground_truth;

        // read all of the data for this candidate
        if (fread(&label_one, sizeof(unsigned long), 1, fp) != 1) return 0;
        if (fread(&label_two, sizeof(unsigned long), 1, fp) != 1) return 0;
        if (fread(&z, sizeof(unsigned long), 1, fp) != 1) return 0;
        if (fread(&y, sizeof(unsigned long), 1, fp) != 1) return 0;
        if (fread(&x, sizeof(unsigned long), 1, fp) != 1) return 0;
        if (fread(&ground_truth, sizeof(unsigned long), 1, fp) != 1) return 0;

        candidates[iv] = MergeCandidate(label_one, label_two, x, y, z, ground_truth);
    }

    // close the file
    fclose(fp);

    // return success
    return 1;
}



////////////////////////////////////////////////////////////////////////
// Preprocessing functions
////////////////////////////////////////////////////////////////////////

static void LabelToIndexMapping(void)
{
    // find which labels are present
    unsigned long maximum_segmentation = (unsigned long)(grid->Maximum() + 0.5) + 1;
    RNBoolean *present_labels = new RNBoolean[maximum_segmentation];
    for (unsigned long iv = 0; iv < maximum_segmentation; ++iv)
        present_labels[iv] = FALSE;

    // iterate over the entire volume to find present labels
    for (long iv = 0; iv < grid->NEntries(); ++iv) {
        unsigned long label = (unsigned long)(grid->GridValue(iv) + 0.5);
        present_labels[label] = TRUE;
    }

    // create the mapping from segments labels to indices
    label_to_index = std::map<unsigned long, unsigned long>();
    unsigned long nunique_labels = 1; /* 1 indexed for this to work */
    for (unsigned long iv = 1; iv < maximum_segmentation; ++iv) {
        if (present_labels[iv] && !label_to_index[iv]) {
            label_to_index[iv] = nunique_labels;
            nunique_labels++;
        }
    }

    // create the mapping from indices to labels
    index_to_label = new unsigned long[nunique_labels];
    nunique_labels = 1;
    for (unsigned long iv = 1; iv < maximum_segmentation; ++iv) {
        if (present_labels[iv]) {
            index_to_label[nunique_labels] = iv;
            nunique_labels++;
        }
    }

    // free memory
    delete[] present_labels;

    // allocate memory for the segmentation vectors
    segmentations = new std::vector<unsigned long>[nunique_labels];
    skeletons = new std::vector<SWCEntry>[nunique_labels];
    skeleton_endpoints = new std::vector<R3Point>[nunique_labels];
    for (unsigned long iv = 0; iv < nunique_labels; ++iv) {
        segmentations[iv] = std::vector<unsigned long>();
        skeletons[iv] = std::vector<SWCEntry>();
        skeleton_endpoints[iv] = std::vector<R3Point>();

        ReadSWCFile(iv);
    }

    // iterate over the entire volume
    for (long iv = 0; iv < grid->NEntries(); ++iv) {
        unsigned long label = (unsigned long)(grid->GridValue(iv) + 0.5);
        // skip background
        if (label == 0) continue;

        // get this index
        unsigned long index = label_to_index[label];
        rn_assertion((0 < index) && (index < nunique_labels));

        // add to the vector
        segmentations[index].push_back(iv);
    }

    // remove the grid from memory 
    delete grid;
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



static void DrawSkeleton(int index)
{
    if(skeletons[index].size()) {
        RNLoadRgb(RNwhite_rgb);
        glLineWidth(5.0);

        glBegin(GL_LINES);
        // go through all of the entries
        for(unsigned int ie = 0; ie < skeletons[index].size(); ++ie) {
            SWCEntry entry = skeletons[index][ie];
            if(entry.parent_sample == -1) continue;
            SWCEntry parent = skeletons[index][entry.parent_sample - 1];
            rn_assertion(parent.sample_number == entry.parent_sample);

            // draw the line
            glVertex3f(entry.X(), entry.Y(), entry.Z());
            glVertex3f(parent.X(), parent.Y(), parent.Z());
        }
        glEnd();
        glLineWidth(1.0);
    }
}



static void DrawIndividualSegment(int index)
{
    glBegin(GL_POINTS);
    for (unsigned int iv = 0; iv < segmentations[index].size(); ++iv) {
        if (RNRandomScalar() > 1.0 / downsample_rate) continue;
        int ix, iy, iz;
        IndexToIndices(segmentations[index][iv], ix, iy, iz);
        glVertex3f(ix, iy, iz);
    }
    glEnd();
    if (show_skeleton) DrawSkeleton(index);
}



static void DrawSegmentations(void)
{
    // push the transformation
    transformation.Push();

    if (show_merge_candidate) {
        MergeCandidate candidate = candidates[candidate_index];
        if (candidate.GroundTruth()) {
            RNLoadRgb(RNblue_rgb);
            DrawIndividualSegment(label_to_index[candidate.LabelOne()]);
            RNLoadRgb(RNgreen_rgb);
            DrawIndividualSegment(label_to_index[candidate.LabelTwo()]);
        }
        else {
            RNLoadRgb(RNred_rgb);
            DrawIndividualSegment(label_to_index[candidate.LabelOne()]);
            RNLoadRgb(RNyellow_rgb);
            DrawIndividualSegment(label_to_index[candidate.LabelTwo()]);
        }

        // draw the central point
        RNLoadRgb(RNwhite_rgb);

        // draw the bounding box
        R3Box bounding_box = R3Box(
            candidate.X() - maximum_distance / resolution[RN_X], 
            candidate.Y() - maximum_distance / resolution[RN_Y], 
            candidate.Z() - maximum_distance / resolution[RN_Z], 
            candidate.X() + maximum_distance / resolution[RN_X], 
            candidate.Y() + maximum_distance / resolution[RN_Y], 
            candidate.Z() + maximum_distance / resolution[RN_Z]
        );
        bounding_box.Outline();

        // draw the central point
        RNScalar radius = 10.0;
        R3Box central_point = R3Box(
            candidate.X() - radius / resolution[RN_X],
            candidate.Y() - radius / resolution[RN_Y], 
            candidate.Z() - radius / resolution[RN_Z], 
            candidate.X() + radius / resolution[RN_X], 
            candidate.Y() + radius / resolution[RN_Y], 
            candidate.Z() + radius / resolution[RN_Z]
        );
        central_point.Draw();
    }
    else {
        RNLoadRgb(RNblue_rgb);
        DrawIndividualSegment(segmentation_index);    
    }

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
        RNLoadRgb(RNwhite_rgb);
        world_box.Outline();
    }

    // draw machine labels and skeletons
    DrawSegmentations();

    // epilogue
    glEnable(GL_LIGHTING);

    // write the title
    char title[4096];
    if (show_merge_candidate) {
        sprintf(title, "Skeleton Visualizer - %d\n", candidate_index);    
    }
    else {
        sprintf(title, "Skeleton Visualizer - %lu\n", index_to_label[segmentation_index]);
    }    
    glutSetWindowTitle(title);


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
            if (!show_merge_candidate) {
                segmentation_index--;
                if(segmentation_index < 1) segmentation_index = 1;
            }
            else {
                candidate_index--;
                if (candidate_index < 0) candidate_index = 0;
            }
            break;
        }

        case GLUT_KEY_RIGHT: {
            if (!show_merge_candidate) {
                segmentation_index++;
                if(segmentation_index >= label_to_index.size())
                    segmentation_index = label_to_index.size() - 1;
            }
            else {
                candidate_index++;
                if (candidate_index >= ncandidates) candidate_index = ncandidates - 1;
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
            show_merge_candidate = 1 - show_merge_candidate;
            break;
        }

        case 'S':
        case 's': {
            show_skeleton = 1 - show_skeleton;
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
    while(argc > 0) {
        if((*argv)[0] == '-') {
            if(!strcmp(*argv, "-v")) print_verbose = 1;
            else if(!strcmp(*argv, "-debug")) print_debug = 1;
            else if(!strcmp(*argv, "-max_distance")) { argv++; argc--; maximum_distance = atoi(*argv); } 
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

    if(!ReadData()) exit(-1);

    // map the labels to indices
    LabelToIndexMapping();

    // read all of the merge candidates
    ReadMergeCandidates();

    // set world box
    world_box = R3Box(0, 0, 0, resolution[RN_X] * grid_size[RN_X], resolution[RN_Y] * grid_size[RN_Y], resolution[RN_Z] * grid_size[RN_Z]);
    
    // get the transformation
    transformation = R3Affine(R4Matrix(resolution[RN_X], 0, 0, 0, 0, resolution[RN_Y], 0, 0, 0, 0, resolution[RN_Z], 0, 0, 0, 0, 1));

    // create viewer
    viewer = CreateViewer();
    if(!viewer) exit(-1);

    // initialize GLUT
    GLUTInit(&argc, argv);

    // run GLUT interface
    GLUTMainLoop();

    // return success
    return 1;
}