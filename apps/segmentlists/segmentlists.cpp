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


// class declarations

struct RNMeta;
struct SWCEntry;



// program arguments

// I/O flags
static int print_debug = 0;
static int print_verbose = 0;
// dataset to examine
static const char *prefix = NULL;
static const char *filename = NULL;



// program variables

static RNScalar resolution[3] = { 6, 6, 30 };
static int grid_size[3] = { -1, -1, -1 };
static R3Affine transformation = R3null_affine;
static R3Viewer *viewer = NULL;
static R3Box world_box = R3null_box;



// voxel grids

static R3Grid *rhoana_grid = NULL;
static R3Grid *gold_grid = NULL;



// skeleton variables

static std::vector<long> *topological_joints = NULL;
static std::vector<long> *topological_endpoints = NULL;
static std::vector<SWCEntry> *neutu_skeletons = NULL;



// display variables

//static int selected_slice_index = 0;
//static int show_slice = 0;
//static int projection_dim = RN_Z;
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
static std::vector<long> *golds = NULL;
static long maximum_segmentation = -1;
static long maximum_gold = -1;
static long *segmentation_to_gold = NULL;
static std::vector<long> *segmentations_per_gold = NULL;
static std::vector<std::vector<long> > segmentation_lists;





//static std::map<unsigned long, unsigned long> label_to_index;
//static unsigned long *index_to_label = NULL;
//static std::vector<unsigned long> *segmentations = NULL;
//static std::vector<R3Point> *skeleton_endpoints = NULL;



// display variables

static int cycle_colors = 0;
static int show_bbox = 1;
static int show_split_errors = 0;
static int show_skeleton = 1;
static int show_undetermined = 0;
//static int circular_index = 0;
//static int circle_size = 14;
//static int show_feature_box = 1;
static int show_only_positives = 0;
static int show_segmentation = 1;
static int segmentation_index = 1;
static int show_list_candidates = 1;
static int list_index = 0;
static int show_gold = 0;
static int gold_index = 1;
static int show_slice = 0;
static int selected_slice_index = 0;
//static int show_legend = 1;
//static int segmentation_index = 0;

//static int show_segments = 1;
//static int *final_labels = NULL;
//static int show_output = 0;
static RNScalar downsample_rate = 6.0;




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

    // read in gold information
    if (!fgets(comment, 4096, fp)) return 0;
    if (fscanf(fp, "%s %s\n", meta_data.gold_filename, meta_data.gold_dataset) != 2) return 0;

    // read image
    if (!fgets(comment, 4096, fp)) return 0;
    if (!fgets(comment, 4096, fp)) return 0;

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



static int ReadSWCFile(int label)
{
    char input_filename[4096];
    sprintf(input_filename, "skeletons/NeuTu/%s/tree_%d.swc", prefix, label);

    std::ifstream fd(input_filename);
    if(!fd.is_open()) {
        if (print_debug) fprintf(stderr, "Failed to read %s\n", input_filename);
        return 0;
    }

    std::string line;
    bool first_iteration = TRUE;
    while(std::getline(fd, line)) {
        if(first_iteration) {
            first_iteration = FALSE;
            continue;
        }

        int sample_number, structure_identifier, parent_sample;
        RNScalar x_position, y_position, z_position, radius;
        sscanf(line.c_str(), "%d %d %lf %lf %lf %lf %d", &sample_number, &structure_identifier, &x_position,
            &y_position, &z_position, &radius, &parent_sample);

        neutu_skeletons[label].push_back(SWCEntry(sample_number, structure_identifier, x_position, y_position, z_position, radius, parent_sample));
    }

    fd.close();

    // return success
    return 1;
}



static void ReadTopologicalSkeleton(int label)
{
    if (label > 5000) return;

    char input_filename[4096];
    sprintf(input_filename, "skeletons/topological-thinning/%s/skeleton-%d.pts", prefix, label);

    FILE *fp = fopen(input_filename, "rb");
    if (!fp) return;

    // read the number of labels
    long npts;
    if (fread(&npts, sizeof(long), 1, fp) != 1) return;

    for (long pt = 0; pt < npts; ++pt) {
        long index;
        if (fread(&index, sizeof(long), 1, fp) != 1) return;

        if (index < 0) {
            index = -1 * index;
            topological_joints[label].push_back(index);
            topological_endpoints[label].push_back(index);
        }
        else topological_joints[label].push_back(index);
    }

    fclose(fp);
}



static int ReadData(void)
{
    // start statistics
    RNTime start_time;
    start_time.Read();

    R3Grid **rhoana_grids = RNReadH5File(meta_data.rhoana_filename, meta_data.rhoana_dataset);
    rhoana_grid = rhoana_grids[0];
    delete[] rhoana_grids;

    grid_size[RN_X] = rhoana_grid->XResolution();
    grid_size[RN_Y] = rhoana_grid->YResolution();
    grid_size[RN_Z] = rhoana_grid->ZResolution();

    R3Grid **gold_grids = RNReadH5File(meta_data.gold_filename, meta_data.gold_dataset);
    gold_grid = gold_grids[0];
    delete[] gold_grids;

    // print statistics
    if(print_verbose) {
        printf("Read voxel grids...\n");
        printf("  Time = %.2f seconds\n", start_time.Elapsed());
        printf("  Grid Size = (%d %d %d)\n", rhoana_grid->XResolution(), rhoana_grid->YResolution(), rhoana_grid->ZResolution());
        printf("  Resolution = (%0.2lf %0.2lf %0.2lf)\n", resolution[RN_X], resolution[RN_Y], resolution[RN_Z]);
    }

    // return success
    return 1;
}



static int ReadSegmentCandidates(void)
{
    // open the file
    FILE *fp = fopen(filename, "rb");
    if (!fp) { fprintf(stderr, "Failed to read%s\n", filename); return 0; }

    long nlists;
    if (fread(&nlists, sizeof(long), 1, fp) != 1) return 0;

    // create a new list vector
    segmentation_lists = std::vector<std::vector<long> >(nlists);

    for (long ig = 0; ig < nlists; ++ig) {
        long nsegments;
        if (fread(&nsegments, sizeof(long), 1, fp) != 1) return 0;
        segmentation_lists.push_back(std::vector<long> (nsegments));

        for (long is = 0; is < nsegments; ++is) {
            long segment_id;
            if (fread(&segment_id, sizeof(long), 1, fp) != 1) return 0; 
            segmentation_lists[ig].push_back(segment_id);
        }
    }

    fclose(fp);

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



static long IndicesToIndex(long ix, long iy, long iz)
{
    return iz * grid_size[RN_X] * grid_size[RN_Y] + iy * grid_size[RN_X] + ix;
}



static RNRgb Color(unsigned long value)
{
    value = value;
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

    // constants for segmentation to gold mapping
    static const RNScalar low_threshold = 0.1;
    static const RNScalar high_threshold = 0.8;

    // get the maximum values for each grid
    maximum_segmentation = (long) (rhoana_grid->Maximum() + 0.5) + 1;
    maximum_gold = (long) (gold_grid->Maximum() + 0.5) + 1;

    // create a vector for each valid ID
    segmentations = new std::vector<long>[maximum_segmentation];
    for (long iv = 0; iv < maximum_segmentation; ++iv)
        segmentations[iv] = std::vector<long>();
    golds = new std::vector<long>[maximum_gold];
    for (long iv = 0; iv < maximum_gold; ++iv) 
        golds[iv] = std::vector<long>();

    // go through all voxels to see if it belongs to the boundary
    for (int iz = 1; iz < grid_size[RN_Z] - 1; ++iz) {
        for (int iy = 1; iy < grid_size[RN_Y] - 1; ++iy) {
            for (int ix = 1; ix < grid_size[RN_X] - 1; ++ix) {
                long iv = IndicesToIndex(ix, iy, iz);

                long segment = (long) (rhoana_grid->GridValue(ix, iy, iz) + 0.5);
                long gold = (long) (gold_grid->GridValue(ix, iy, iz) + 0.5);

                // go through all six adjacent neighbors
                if ((long)(rhoana_grid->GridValue(ix - 1, iy, iz) + 0.5) != segment ||
                    (long)(rhoana_grid->GridValue(ix + 1, iy, iz) + 0.5) != segment ||
                    (long)(rhoana_grid->GridValue(ix, iy - 1, iz) + 0.5) != segment ||
                    (long)(rhoana_grid->GridValue(ix, iy + 1, iz) + 0.5) != segment ||
                    (long)(rhoana_grid->GridValue(ix, iy, iz - 1) + 0.5) != segment ||
                    (long)(rhoana_grid->GridValue(ix, iy, iz + 1) + 0.5) != segment)
                {
                    segmentations[segment].push_back(iv);
                }
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

    // create a mapping from segmentation values to gold values
    long *nvoxels_per_segment = new long[maximum_segmentation];
    for (long iv = 0; iv < maximum_segmentation; ++iv)
        nvoxels_per_segment[iv] = 0;
    long **seg2gold_overlap = new long *[maximum_segmentation];
    for (long is = 0; is < maximum_segmentation; ++is) {
        seg2gold_overlap[is] = new long[maximum_gold];
        for (long ig = 0; ig < maximum_gold; ++ig) {
            seg2gold_overlap[is][ig] = 0;
        }
    }

    // iterate over every voxel
    for (long iv = 0; iv < rhoana_grid->NEntries(); ++iv) {
        long segment = (long) (rhoana_grid->GridValue(iv) + 0.5);
        long gold = (long) (gold_grid->GridValue(iv) + 0.5);
        nvoxels_per_segment[segment]++;
        seg2gold_overlap[segment][gold]++;
    }

    segmentation_to_gold = new long[maximum_segmentation];
    for (long is = 1; is < maximum_segmentation; ++is) {
        long gold_id = 0;
        long gold_max_value = 0;

        // only gets label of 0 if the number of non zero voxels is below threshold
        for (long ig = 1; ig < maximum_gold; ++ig) {
            if (seg2gold_overlap[is][ig] > gold_max_value) {
                gold_max_value = seg2gold_overlap[is][ig];
                gold_id = ig;
            }
        }

        // the number of non zero pixels must be greater than low threshold
        if (gold_max_value / (double)nvoxels_per_segment[is] < low_threshold) segmentation_to_gold[is] = 0;
        else if (gold_max_value / (double)(nvoxels_per_segment[is] - seg2gold_overlap[is][0]) > high_threshold) segmentation_to_gold[is] = gold_id;
        else segmentation_to_gold[is] = 0;
    }

    // create a cache of segmentations for every gold label
    segmentations_per_gold = new std::vector<long>[maximum_gold];
    for (long ig = 0; ig < maximum_gold; ++ig) 
        segmentations_per_gold[ig] = std::vector<long>();
    for (long is = 1; is < maximum_segmentation; ++is) {
        long gold_index = segmentation_to_gold[is];
        segmentations_per_gold[gold_index].push_back(is);
    }

    // create arrays for topological skeletons
    topological_joints = new std::vector<long>[maximum_segmentation];
    topological_endpoints = new std::vector<long>[maximum_segmentation];
    neutu_skeletons = new std::vector<SWCEntry>[maximum_segmentation];
    for (long iv = 0; iv < maximum_segmentation; ++iv) {
        topological_joints[iv] = std::vector<long>();
        topological_endpoints[iv] = std::vector<long>();
        neutu_skeletons[iv] = std::vector<SWCEntry>();
    }

    // read all of the topological skeletons
    for (long iv = 0; iv < maximum_segmentation; ++iv) {
        ReadSWCFile(iv);
        ReadTopologicalSkeleton(iv);
    }
    
    // free memory
    for (long is = 0; is < maximum_segmentation; ++is) 
        delete[] seg2gold_overlap[is];
    delete[] seg2gold_overlap;
    delete[] nvoxels_per_segment;

    if (print_verbose) {
        printf("Preprocessing...\n");
        printf("  Time = %0.2f seconds\n", start_time.Elapsed());
        printf("  Maximum Segment = %ld\n", maximum_segmentation);
        printf("  Maximum Gold = %ld\n", maximum_gold);
    }
}



////////////////////////////////////////////////////////////////////////
// Drawing utility functions
////////////////////////////////////////////////////////////////////////

static void DrawNeuTuSkeleton(int segment_index)
{
    RNRgb contrast_color = RNRgb(1.0 - background_color[0], 1.0 - background_color[1], 1.0 - background_color[2]);
    RNLoadRgb(contrast_color);

    glLineWidth(5.0);
    glBegin(GL_LINES);
    // go through all of the entries
    for(unsigned int ie = 0; ie < neutu_skeletons[segment_index].size(); ++ie) {
        SWCEntry entry = neutu_skeletons[segment_index][ie];
        if(entry.parent_sample == -1) continue;
        SWCEntry parent = neutu_skeletons[segment_index][entry.parent_sample - 1];
        rn_assertion(parent.sample_number == entry.parent_sample);

        // draw the line
        glVertex3f(entry.X(), entry.Y(), entry.Z());
        glVertex3f(parent.X(), parent.Y(), parent.Z());
    }
    glEnd();
    glLineWidth(1.0);
}



static void DrawToplogicalSkeleton(int segment_index)
{
    RNRgb contrast_color = RNRgb(1.0 - background_color[0], 1.0 - background_color[1], 1.0 - background_color[2]);
    RNLoadRgb(contrast_color);

    glPointSize(2.0);
    glBegin(GL_POINTS);
    // go through all points along the skeleton
    for (unsigned int iv = 0; iv < topological_joints[segment_index].size(); ++iv) {
        long index = topological_joints[segment_index][iv];

        // get the cartesian coordinates
        int ix, iy, iz;
        IndexToIndices(index, ix, iy, iz);

        // draw the veretx
        glVertex3f(ix, iy, iz);
    }
    glEnd();
}



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

    if (show_skeleton == 1) DrawToplogicalSkeleton(segment_index);
    if (show_skeleton == 2) DrawNeuTuSkeleton(segment_index);
}



static void DrawGold(int gold_index)
{
    glBegin(GL_POINTS);
    for (unsigned int iv = 0; iv < golds[gold_index].size(); ++iv) {
        // faster rendering with downsampling
        if (RNRandomScalar() > 1.0 / downsample_rate) continue;

        // get the coordinates from the linear index
        int ix, iy, iz;
        IndexToIndices(golds[gold_index][iv], ix, iy, iz);
        glVertex3f(ix, iy, iz);
    }
}



static void DrawPointClouds(void)
{
    // push the transformation
    transformation.Push();

    
    if (show_list_candidates and show_segmentation) {
        for (unsigned long is = 0; is < segmentation_lists[list_index].size(); ++is) {
            long segment_id = segmentation_lists[list_index][is];

            RNLoadRgb(Color(segment_id + cycle_colors));
            DrawSegment(segment_id);
        }
    }
    else if (show_list_candidates and not show_segmentation) {
        for (unsigned long is = 0; is < segmentation_lists[list_index].size(); ++is) {
            long segment_id = segmentation_lists[list_index][is];

            if (show_skeleton == 1) DrawToplogicalSkeleton(segment_id);
            else if (show_skeleton == 2) DrawNeuTuSkeleton(segment_id);
        }
    }
    else if (show_segmentation) {
        RNLoadRgb(Color(segmentation_index + cycle_colors));
        DrawSegment(segmentation_index);

        // show the neighbors with which it should merge
        long gold_index = segmentation_to_gold[segmentation_index];

        if (show_split_errors && gold_index) {
            for (unsigned int is = 0; is < segmentations_per_gold[gold_index].size(); ++is) {
                long neighbor_segment_index = segmentations_per_gold[gold_index][is];
                if (segmentation_index == neighbor_segment_index) continue;
                RNLoadRgb(Color(maximum_segmentation - segmentation_index + cycle_colors));
                DrawSegment(neighbor_segment_index);
            }
        }
    }

    if (show_gold) {
        RNLoadRgb(Color(gold_index));
        DrawGold(gold_index);
    }

    // pop the transformation
    transformation.Pop();
}



/*static void GLUTDrawText(const R2Point& position, const char *s)
{
   // draw text string s at position
   glRasterPos2d(position[0], position[1]);
   while (*s) glutBitmapCharacter(GLUT_BITMAP_HELVETICA_18, *(s++));
}*/



/*// Draw a legend in the lower-left corner
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

   if (show_merge_candidate) {     
       // show merge candidate information
       MergeCandidate candidate = candidates[candidate_index];
       sprintf(legend, "Candidate index: %u, Label One: %lu (blue/red), Label Two: %lu (green/yellow)", candidate_index, candidate.LabelOne(), candidate.LabelTwo());
       GLUTDrawText(R2Point(10, 50), legend);
       sprintf(legend,  "Index One: %lu, Index Two: %lu", label_to_index[candidate.LabelOne()], label_to_index[candidate.LabelTwo()]);
       GLUTDrawText(R2Point(10, 30), legend);

       // show control sequence 
       GLUTDrawText(R2Point(10, 10), "C - show single neurons");      
   } 
   else {
        sprintf(legend, "Segment %lu\n", index_to_label[segmentation_index]);
       GLUTDrawText(R2Point(10, 30), legend);
       GLUTDrawText(R2Point(10, 10), "C - show merge candidates"); 
   }
}*/


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

    if (rhoana_grid) delete rhoana_grid;
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

    /*if(show_legend){
        DrawLegend();
    }*/

    // epilogue
    glEnable(GL_LIGHTING);

    transformation.Push();
    if (show_slice == 1) rhoana_grid->DrawColorSlice(RN_Z, selected_slice_index);
    if (show_slice == 2) gold_grid->DrawColorSlice(RN_Z, selected_slice_index);
    transformation.Pop();


    // write the title
    /*char title[4096];
    if (show_merge_candidate) {
        sprintf(title, "Skeleton Visualizer (Merge Candidates, %s) - %d\n", prefix, candidate_index);    
    }
    else {
        sprintf(title, "Skeleton Visualizer (Single Neurons, %s) - %lu\n", prefix, index_to_label[segmentation_index]);
    }    
    glutSetWindowTitle(title);*/


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



/*static void DecrementIndex(void)
{
    if (candidate_index) --candidate_index;
}*/



/*static void IncrementIndex(void)
{
    if (candidate_index < ncandidates - 1) ++candidate_index;
}*/


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
            if (selected_slice_index >= grid_size[RN_Z])
                selected_slice_index = grid_size[RN_Z] - 1;
            break;
        }

        case GLUT_KEY_DOWN: {
            selected_slice_index -= 1;
            if (selected_slice_index < 0)
                selected_slice_index = 0;
            break;
        }

        case GLUT_KEY_LEFT: {
            if (GLUTmodifiers & GLUT_ACTIVE_SHIFT) {
                --gold_index;
                if (gold_index < 0) 
                    gold_index = 0;
            }
            else {
                if (show_list_candidates) {
                    --list_index;
                    if (list_index < 0)
                        list_index = 0;    
                    printf("Labels: ");
                    for (unsigned long is = 0; is < segmentation_lists[list_index].size(); ++is) {
                        printf("%ld ", segmentation_lists[list_index][is]);
                    }
                    printf("\n");
                }
                else {
                    --segmentation_index;
                    if (segmentation_index < 0)
                        segmentation_index = 0;
                    printf("Label: %d\n", segmentation_index);
                }
                break;
            }
        }

        case GLUT_KEY_RIGHT: {
            if (GLUTmodifiers & GLUT_ACTIVE_SHIFT) {
                ++gold_index;
                if (gold_index >= maximum_gold)
                    gold_index = maximum_gold - 1;
            }
            else {
                if (show_list_candidates) {
                    ++list_index;
                    if (list_index >= (long)segmentation_lists.size())
                        list_index = segmentation_lists.size() - 1;
                    printf("Labels: ");
                    for (unsigned long is = 0; is < segmentation_lists[list_index].size(); ++is) {
                        printf("%ld ", segmentation_lists[list_index][is]);
                    }
                    printf("\n");     
                }
                else {
                    ++segmentation_index;
                    if (segmentation_index >= maximum_segmentation)
                        segmentation_index = maximum_segmentation - 1;

                    printf("Label: %d\n", segmentation_index);
                }
                break;    
            }

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
            cycle_colors = (++cycle_colors) % 10;
            break;
        }

        case 'E':
        case 'e': {
            show_split_errors = 1 - show_split_errors;
            break;
        }

        case 'G':
        case 'g': {
            show_gold = 1 - show_gold;
            break;
        }

        case 'K': 
        case 'k': {
            show_skeleton = (show_skeleton + 1) % 3;
           break;
        }

        case 'L':
        case 'l': {
            show_list_candidates = 1 - show_list_candidates;
        }

        case 'P':
        case 'p': {
            show_only_positives = 1 - show_only_positives;
            break;
        }

        case 'S':
        case 's': {
            show_segmentation = 1 - show_segmentation;
            break;
        }

        case 'W':
        case 'w': {
            show_slice = (++show_slice) % 3;
            break;
        }

        case 'U':
        case 'u': {
            show_undetermined = 1 - show_undetermined;
            break;
        }

        case ENTER: {
            background_color[0] = 1.0 - background_color[0];
            background_color[1] = 1.0 - background_color[1];
            background_color[2] = 1.0 - background_color[2];
            break;
        }

        /*case 'C':
        case 'c': {
            show_merge_candidate = 1 - show_merge_candidate;
            break;
        }

        case 'D':
        case 'd': {
            show_segmentation = 1 - show_segmentation;
            break;
        }

        case 'E':
        case 'e': {
            show_segments = 1 - show_segments;
            break;
        }

        case 'F':
        case 'f': {
            show_feature_box = 1 - show_feature_box;
            break;
        }

        case 'G':
        case 'g': {
            circular_index = (circular_index + 1) % circle_size;
        }

        case 'M': 
        case 'm': {
            network_distance = network_distance + 200;
            printf("%d\n", network_distance);
            break;
        }

        case 'N':
        case 'n': {
            network_distance = network_distance - 200;
            printf("%d\n", network_distance);
            break;
        }

        case 'O': 
        case 'o': {
            show_output = 1 - show_output;
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
        }*/

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
            if(!prefix) prefix = *argv;
            else if (!filename) filename = *argv;
            else { fprintf(stderr, "Invalid program argument: %s\n", *argv); return 0; }
        }
        argv++; argc--;
    }

    // error if there is no input name
    if(!prefix && !filename) {
        fprintf(stderr, "Need to supply a prefix for data files and a filename with segmentation lists\n");
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

    /////////////////////////////////
    //// Read in the voxel files ////
    /////////////////////////////////
    
    if (!ReadSegmentCandidates()) exit(-1);
    if (!ReadMetaData(prefix)) exit(-1);
    if (!ReadData()) exit(-1);

    // get all of the preprocessing time
    Preprocessing();

    // read all of the merge candidates
    //if (!ReadMergeCandidates()) return 0;

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