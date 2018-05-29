// Source file for nuclear visualizer



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


// class declarations
struct RNMeta;
struct SWCEntry;
struct SplitCandidate;



// program arguments

static int print_debug = 0;
static int print_verbose = 0;
// maximum distance in nanometers
static int network_distance = 600;
static int threshold = 20000;
static const char* prefix = NULL;



// program variables

static RNScalar resolution[3] = { 6, 6, 30 };
static int grid_size[3] = { -1, -1, -1 };
static R3Affine transformation = R3null_affine;
static R3Viewer *viewer = NULL;
static R3Box world_box = R3null_box;
static RNBoolean training = FALSE;
static RNBoolean validation = FALSE;



// voxel grids

static R3Grid *grid = NULL;
//static R3Grid *image_grid = NULL;
static R3Grid *gold_grid = NULL;
static int selected_slice_index = 0;
static int show_slice = 0;
static int projection_dim = RN_Z;



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

static std::map<unsigned long, unsigned long> label_to_index;
static unsigned long *index_to_label = NULL;
static std::vector<unsigned long> *segmentations = NULL;
static std::vector<SWCEntry> *skeletons = NULL;
static std::vector<R3Point> *skeleton_endpoints = NULL;
static SplitCandidate *candidates = NULL;



// display variables

static int show_bbox = 1;
static int show_skeleton = 1;
static int show_split_candidate = 1;
static int show_feature_box = 1;
static int show_only_positives = 0;
static int show_segmentation = 1;
static int show_legend = 1;
static int segmentation_index = 0;
static int candidate_index = 0;
static int ncandidates;
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
    sprintf(meta_data_filename, "meta_data/%s.meta", prefix);

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

    // skip gold
    if (!fgets(comment, 4096, fp)) return 0;
    if (fscanf(fp, "%s %s\n", meta_data.gold_filename, meta_data.gold_dataset) != 2) return 0;

    // read image
    if (!fgets(comment, 4096, fp)) return 0;
    if (fscanf(fp, "%s %s\n", meta_data.image_filename, meta_data.image_dataset) != 2) return 0;

    // skip mask
    if (!fgets(comment, 4096, fp)) return 0;
    if (fscanf(fp, "%s %s\n", meta_data.mask_filename, meta_data.mask_dataset) != 2) return 0;

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



struct SplitCandidate {
    // constructors
    SplitCandidate(void) :
    label(0),
    x(0),
    y(0),
    z(0),
    ground_truth(FALSE)
    {
    }

    SplitCandidate(unsigned long label, unsigned long x, unsigned long y, unsigned long z, RNBoolean ground_truth) :
    label(label),
    x(x),
    y(y),
    z(z),
    ground_truth(ground_truth)
    {
    }

    // access variables
    unsigned long Label(void) { return label; }
    unsigned long X(void) { return x; }
    unsigned long Y(void) { return y; }
    unsigned long Z(void) { return z; }
    unsigned long GroundTruth(void) { return ground_truth; }

    // instance variables
    unsigned long label;   
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

    R3Grid **grids = RNReadH5File(meta_data.rhoana_filename, meta_data.rhoana_dataset);
    grid = grids[0];
    delete[] grids;

    grid_size[RN_X] = grid->XResolution();
    grid_size[RN_Y] = grid->YResolution();
    grid_size[RN_Z] = grid->ZResolution();


    //R3Grid **image_grids = RNReadH5File(meta_data.image_filename, meta_data.image_dataset);
    //image_grid = image_grids[0];
    //delete[] image_grids;

    R3Grid **gold_grids = RNReadH5File(meta_data.gold_filename, meta_data.gold_dataset);
    gold_grid = gold_grids[0];
    delete[] gold_grids;

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



static int ReadSplitCandidates(void)
{
    // start statistics
    RNTime start_time;
    start_time.Read();

    // get the candidate filename
    char candidate_filename[4096];
    sprintf(candidate_filename, "features/nuclear/%s-%d-inference.candidates", prefix, threshold);

    // open the file
    FILE *fp = fopen(candidate_filename, "rb");
    if (!fp) { fprintf(stderr, "Failed to read %s\n", candidate_filename); return 0; }

    if (fread(&ncandidates, sizeof(int), 1, fp) != 1) return 0;

    int npositives = 0;
    int nnegatives = 0;
    // read in all of the candidates
    candidates = new SplitCandidate[ncandidates];
    for (int iv = 0; iv < ncandidates; ++iv) {
        unsigned long label;
        unsigned long x;
        unsigned long y;
        unsigned long z;
        unsigned long ground_truth;

        // read all of the data for this candidate
        if (fread(&label, sizeof(unsigned long), 1, fp) != 1) return 0;
        if (fread(&z, sizeof(unsigned long), 1, fp) != 1) return 0;
        if (fread(&y, sizeof(unsigned long), 1, fp) != 1) return 0;
        if (fread(&x, sizeof(unsigned long), 1, fp) != 1) return 0;
        if (fread(&ground_truth, sizeof(unsigned long), 1, fp) != 1) return 0;

        candidates[iv] = SplitCandidate(label, x, y, z, ground_truth);

        if (ground_truth) npositives++;
        else nnegatives++;
    }

    // close the file
    fclose(fp);

    // print statistics
    if (print_verbose) {
        printf("Read %s in %0.2f seconds\n", candidate_filename, start_time.Elapsed());
        printf(" Positive Candidates: %d\n", npositives);
        printf(" Negative Candidates: %d\n", nnegatives);
    }

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
    unsigned long nunique_labels = 0;
    for (unsigned long iv = 0; iv < maximum_segmentation; ++iv) {
        if (present_labels[iv] && !label_to_index[iv]) {
            label_to_index[iv] = nunique_labels;
            nunique_labels++;
        }
    }

    // create the mapping from indices to labels
    index_to_label = new unsigned long[nunique_labels];
    nunique_labels = 0;
    for (unsigned long iv = 0; iv < maximum_segmentation; ++iv) {
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
    int iv = 0;
    for (int iz = 0; iz < grid->ZResolution(); ++iz) {
        for (int iy = 0; iy < grid->YResolution(); ++iy) {
            for (int ix = 0; ix < grid->XResolution(); ++ix, ++iv) {
                unsigned long label = (unsigned long)(grid->GridValue(ix, iy, iz) + 0.5);

                // is this pixel boundary
                RNBoolean boundary = FALSE;
                if (ix > 0 && label != (unsigned long)(grid->GridValue(ix - 1, iy, iz) + 0.5)) boundary = TRUE;
                if (iy > 0 && label != (unsigned long)(grid->GridValue(ix, iy - 1, iz) + 0.5)) boundary = TRUE;
                if (iz > 0 && label != (unsigned long)(grid->GridValue(ix, iy, iz - 1) + 0.5)) boundary = TRUE;
                if (ix < grid_size[RN_X] - 1 && label != (unsigned long)(grid->GridValue(ix + 1, iy, iz) + 0.5)) boundary = TRUE;
                if (iy < grid_size[RN_Y] - 1 && label != (unsigned long)(grid->GridValue(ix, iy + 1, iz) + 0.5)) boundary = TRUE;
                if (iz < grid_size[RN_Z] - 1 && label != (unsigned long)(grid->GridValue(ix, iy, iz + 1) + 0.5)) boundary = TRUE;

                // get this index
                unsigned long index = label_to_index[label];
                rn_assertion((0 <= index) && (index < nunique_labels));

                // add to the vector
                segmentations[index].push_back(iv);
            }
        }
    }
}



////////////////////////////////////////////////////////////////////////
// Drawing utility functions
////////////////////////////////////////////////////////////////////////

static RNRgb
Color(unsigned long value)
{
    /* TODO remove hard coding */
    if (!value) value = 1;

    RNScalar red = (RNScalar) (((107 * value) % 700) % 255) / 255.0;
    RNScalar green = (RNScalar) (((509 * value) % 900) % 255) / 255.0;
    RNScalar blue = (RNScalar) (((200 * value) % 777) % 255) / 255.0;

    return RNRgb(red, green, blue);
}



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
    RNLoadRgb(RNwhite_rgb);

    if(skeletons[index].size()) {
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

    RNLoadRgb(RNwhite_rgb);
}



static void DrawIndividualSegment(int index, int gold_label) {
    glBegin(GL_POINTS);
    for (unsigned int iv = 0; iv < segmentations[index].size(); ++iv) {
        if (RNRandomScalar() > 1.0 / downsample_rate) continue;
        int gold_value = (int)(gold_grid->GridValue(segmentations[index][iv]) + 0.5);
        if (not gold_value) continue;

        if (gold_label != -1) {
            if (gold_value == gold_label) RNLoadRgb(RNgreen_rgb);
            else RNLoadRgb(RNblue_rgb);
        }

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

    if (show_split_candidate) {
        SplitCandidate candidate = candidates[candidate_index];
        if (candidate.GroundTruth()) {
            if (show_segmentation) DrawIndividualSegment(label_to_index[candidate.Label()], (int) (gold_grid->GridValue((float)candidate.X(), (float)candidate.Y(), (float)candidate.Z()) + 0.5));
        }
        else {
            RNLoadRgb(RNred_rgb);
            if (show_segmentation) DrawIndividualSegment(label_to_index[candidate.Label()], -1);
        }


        if (show_feature_box) {
            // draw the central point
            RNLoadRgb(RNwhite_rgb);

            // draw the bounding box
            R3Box bounding_box = R3Box(
                candidate.X() - network_distance / resolution[RN_X], 
                candidate.Y() - network_distance / resolution[RN_Y], 
                candidate.Z() - network_distance / resolution[RN_Z], 
                candidate.X() + network_distance / resolution[RN_X], 
                candidate.Y() + network_distance / resolution[RN_Y], 
                candidate.Z() + network_distance / resolution[RN_Z]
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
    }
    else {
        RNLoadRgb(Color(segmentation_index));
        if (show_segmentation) DrawIndividualSegment(segmentation_index, -1);    
    }

    // draw the slice if desired
    RNLoadRgb(RNwhite_rgb);
    //if (show_slice == 1) image_grid->DrawSlice(projection_dim, selected_slice_index);
    if (show_slice == 2) gold_grid->DrawColorSlice(projection_dim, selected_slice_index);
    else if (show_slice == 3) grid->DrawColorSlice(projection_dim, selected_slice_index);


    // pop the transformation
    transformation.Pop();
}


static void GLUTDrawText(const R2Point& position, const char *s)
{
   // draw text string s at position
   glRasterPos2d(position[0], position[1]);
   while (*s) glutBitmapCharacter(GLUT_BITMAP_HELVETICA_18, *(s++));
}

// Draw a legend in the lower-left corner
static void DrawLegend(void)
{
   // set projection matrix
   glMatrixMode(GL_PROJECTION);
   glPushMatrix();
   glLoadIdentity();
   gluOrtho2D(0, GLUTwindow_width, 0, GLUTwindow_height);

   // set model view matrix
   glMatrixMode(GL_MODELVIEW);
   glPushMatrix();
   glLoadIdentity();


   char legend[4096];

   if (show_split_candidate) {     
       // show split candidate information
       SplitCandidate candidate = candidates[candidate_index];
       sprintf(legend, "Candidate index: %u, Label: %lu", candidate_index, candidate.Label());
       GLUTDrawText(R2Point(10, 50), legend);
       sprintf(legend,  "Index: %lu", label_to_index[candidate.Label()]);
       GLUTDrawText(R2Point(10, 30), legend);

       // show control sequence 
       GLUTDrawText(R2Point(10, 10), "C - show single neurons");      
   } 
   else {
        sprintf(legend, "Segment %lu\n", index_to_label[segmentation_index]);
       GLUTDrawText(R2Point(10, 30), legend);
       GLUTDrawText(R2Point(10, 10), "C - show split candidates"); 
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

    if(show_legend){
        DrawLegend();
    }

    // epilogue
    glEnable(GL_LIGHTING);


    // write the title
    char title[4096];
    if (show_split_candidate) {
        sprintf(title, "Nuclear Visualizer (Split Candidates, %s) - %d\n", prefix, candidate_index);    
    }
    else {
        sprintf(title, "Nuclear Visualizer (Single Neurons, %s) - %lu\n", prefix, index_to_label[segmentation_index]);
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



static void DecrementIndex(void)
{
    if (candidate_index) --candidate_index;
}



static void IncrementIndex(void)
{
    if (candidate_index < ncandidates - 1) ++candidate_index;
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
        case GLUT_KEY_UP: {
            selected_slice_index += 1;
            if (selected_slice_index >= grid_size[projection_dim])
                selected_slice_index = grid_size[projection_dim] - 1;
            break;
        }

        case GLUT_KEY_DOWN: {
            selected_slice_index -= 1;
            if (selected_slice_index < 0)
                selected_slice_index = 0;
            break;
        }

        case GLUT_KEY_LEFT: {
            if (!show_split_candidate) {
                segmentation_index--;
                if(segmentation_index < 0) segmentation_index = 0;
            }
            else {
                DecrementIndex();
                while (show_only_positives && candidate_index && !candidates[candidate_index].ground_truth)
                    DecrementIndex();
            }
            break;
        }

        case GLUT_KEY_RIGHT: {
            if (!show_split_candidate) {
                segmentation_index++;
                if(segmentation_index >= (int)label_to_index.size())
                    segmentation_index = label_to_index.size() - 1;
            }
            else {
                IncrementIndex();
                while (show_only_positives && candidate_index < ncandidates - 1 && !candidates[candidate_index].ground_truth)
                    IncrementIndex();
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
            show_split_candidate = 1 - show_split_candidate;
            break;
        }

        case 'D':
        case 'd': {
            show_segmentation = 1 - show_segmentation;
            break;
        }

        case 'F':
        case 'f': {
            show_feature_box = 1 - show_feature_box;
            break;
        }

        case 'P':
        case 'p': {
            show_only_positives = 1 - show_only_positives;
            break;
        }

        case 'S':
        case 's': {
            show_skeleton = 1 - show_skeleton;
            break;
        }

        case 'W':
        case 'w': {
            show_slice = (++show_slice) % 4;
            break;
        }

        case 'L':
        case 'l': {
            show_legend = 1 - show_legend;
            break;
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
            else if (!strcmp(*argv, "-threshold")) { argv++; argc--; threshold = atoi(*argv); }
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

    if (!ReadMetaData(prefix)) exit(-1);
    if (!ReadData()) exit(-1);

    // map the labels to indices
    LabelToIndexMapping();

    // read all of the split candidates
    if (!ReadSplitCandidates()) return 0;

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