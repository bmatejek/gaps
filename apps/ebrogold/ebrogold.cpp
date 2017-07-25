// Source file for generating ground truth for EbroNet


// include files

#include "R3Graphics/R3Graphics.h"
#include "RNDataStructures/RNDataStructures.h"
#include "fglut/fglut.h"
#include <vector>
#include <map>
#include <set>



// GLUT defines

#define ENTER 13
#define CTRLZ 26
#define ESCAPE 27
#define SPACEBAR 32



// grid enumeration defines

#define GRID_ONE 0
#define GRID_TWO 1
#define NGRIDS 2


////////////////////////////////////////////////////////////////////////
// Helper struct defintions
////////////////////////////////////////////////////////////////////////

struct RNMeta {
    // instance variables
    int resolution[3];
    char prefix[4096];
    char image_filename[4096];
    char image_dataset[128];
    char rhoana_filename[4096];
    char rhoana_dataset[128];
    R3Box world_box;
    R3Box scaled_box;
};



struct MergeCandidate {
    MergeCandidate(unsigned long label_one, unsigned long label_two, unsigned long index_one, unsigned long index_two, R3Point center) :
    label_one(label_one),
    label_two(label_two),
    index_one(index_one),
    index_two(index_two),
    center(center)
    {}

    unsigned long label_one;
    unsigned long label_two;
    unsigned long index_one;
    unsigned long index_two;
    R3Point center;
};




// program arguments

static int print_verbose = 0;
static const char *prefixes[NGRIDS] = { NULL, NULL };
static unsigned long threshold = 20000;
// window radius in nanometers
static RNScalar window_radius = 600;



// program variables

static int resolution[3] = { -1, -1, -1 };
static int grid_size[3] = { -1, -1, -1 };
static RNMeta meta_data[NGRIDS];
static R3Affine transformation = R3null_affine;
static R3Viewer *viewer = NULL;
static R3Box world_box = R3null_box;



// voxel grids

static R3Grid *image_grids[NGRIDS] = { NULL, NULL };
// combined image grid
static R3Grid *image_grid = NULL;
static R3Grid *rhoana_grids[NGRIDS] = { NULL, NULL };



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

static std::map<unsigned long, unsigned long> label_to_index[NGRIDS];
static unsigned long *index_to_label[NGRIDS] = { NULL, NULL };
static std::vector<unsigned long> *segmentations[NGRIDS] = { NULL, NULL };
static std::vector<MergeCandidate> candidates;
static std::vector<RNScalar> overlap_scores;



// display variables

static int show_bbox = 1;
static RNScalar downsample_rate = 3.0;
static int candidate_index = 0;



// slice display variables

static int show_slice = 1;
static int selected_slice_index[3] = { 0, 0, 0 };
static int projection_dim = RN_Z;



// save decisions

enum DECISION { YES, NO, UNDECIDED };
std::vector<enum DECISION> decisions = std::vector<enum DECISION>();




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
    
    // skip affinities
    if (!fgets(comment, 4096, fp)) return 0;
    if (!fgets(comment, 4096, fp)) return 0;

    // skip boundaries
    if (!fgets(comment, 4096, fp)) return 0;
    if (!fgets(comment, 4096, fp)) return 0;

    // skip gold
    if (!fgets(comment, 4096, fp)) return 0;
    if (!fgets(comment, 4096, fp)) return 0;

    // read image
    if (!fgets(comment, 4096, fp)) return 0;
    if (fscanf(fp, "%s %s\n", meta_data[index].image_filename, meta_data[index].image_dataset) != 2) return 0;

    // skip mask
    if (!fgets(comment, 4096, fp)) return 0;
    if (!fgets(comment, 4096, fp)) return 0;

    // read segmentation
    if (!fgets(comment, 4096, fp)) return 0;
    if (fscanf(fp, "%s %s\n", meta_data[index].rhoana_filename, meta_data[index].rhoana_dataset) != 2) return 0;

    // skip synapse
    if (!fgets(comment, 4096, fp)) return 0;
    if (!fgets(comment, 4096, fp)) return 0;

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
    R3Grid **grids = RNReadH5File(meta_data[index].image_filename, meta_data[index].image_dataset);
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

    image_grids[index]->SetWorldToGridTransformation(meta_data[index].world_box);

    // print information
    if (print_verbose) {
        printf("Read h5 files in %0.2f seconds for %s:\n", start_time.Elapsed(), meta_data[index].prefix);
        printf("  Dimensions: (%d, %d, %d)\n", grid_size[RN_X], grid_size[RN_Y], grid_size[RN_Z]);
        printf("  Resolution: (%d, %d, %d)\n", resolution[RN_X], resolution[RN_Y], resolution[RN_Z]);
    }


    // return success 
    return 1;
}



static int ReadOverlapCandidates(void)
{
    // get filename
    char overlap_filename[4096];
    sprintf(overlap_filename, "features/ebro/%s-%s-%lu-%dnm.candidates", prefixes[GRID_ONE], prefixes[GRID_TWO], threshold, (int)(window_radius + 0.5));

    // open file
    FILE *fp = fopen(overlap_filename, "rb"); 
    if (!fp) { fprintf(stderr, "Failed to open %s\n", overlap_filename); return 0; }

    // find the number of candidates
    unsigned long ncandidates;
    if (fread(&ncandidates, sizeof(unsigned long), 1, fp) != 1) return 0;

    // initialize empty candidate vector
    candidates = std::vector<struct MergeCandidate>();

    // read all candidates
    for (unsigned long iv = 0; iv < ncandidates; ++iv) {
        unsigned long label_one;
        unsigned long label_two;
        unsigned long xcenter;
        unsigned long ycenter;
        unsigned long zcenter;

        if (fread(&label_one, sizeof(unsigned long), 1, fp) != 1) return 0;
        if (fread(&label_two, sizeof(unsigned long), 1, fp) != 1) return 0;
        if (fread(&xcenter, sizeof(unsigned long), 1, fp) != 1) return 0;
        if (fread(&ycenter, sizeof(unsigned long), 1, fp) != 1) return 0;
        if (fread(&zcenter, sizeof(unsigned long), 1, fp) != 1) return 0;

        unsigned long index_one = label_to_index[GRID_ONE][label_one];
        unsigned long index_two = label_to_index[GRID_TWO][label_two];

        struct MergeCandidate candidate = MergeCandidate(label_one, label_two, index_one, index_two, R3Point(xcenter, ycenter, zcenter));
        candidates.push_back(candidate);
    }    

    // close file
    fclose(fp);

    if (print_verbose) printf("No. Candidates: %lu\n", ncandidates);

    // get count filename
    char count_filename[4096];
    sprintf(count_filename, "features/ebro/%s-%s-%lu-%dnm.counts", prefixes[GRID_ONE], prefixes[GRID_TWO], threshold, (int)(window_radius + 0.5));

    // open file
    fp = fopen(count_filename, "rb");
    if (!fp) { fprintf(stderr, "Failed to open %s\n", count_filename); return 0; }

    // find the number of candidates
    unsigned long ncount_candidates;
    if (fread(&ncount_candidates, sizeof(unsigned long), 1, fp) != 1) return 0;
    rn_assertion(ncount_candidates == ncandidates);

    overlap_scores = std::vector<RNScalar>();

    // read all scores
    for (unsigned long iv = 0; iv < ncount_candidates; ++iv) {
        unsigned long count_one;
        unsigned long count_two;
        unsigned long overlap_count;
        RNScalar score;

        if (fread(&count_one, sizeof(unsigned long), 1, fp) != 1) return 0;
        if (fread(&count_two, sizeof(unsigned long), 1, fp) != 1) return 0;
        if (fread(&overlap_count, sizeof(unsigned long), 1, fp) != 1) return 0;
        if (fread(&score, sizeof(RNScalar), 1, fp) != 1) return 0;

        overlap_scores.push_back(score);
    }

    // close file
    fclose(fp);

    // return success
    return 1;
}



static void
ReadDecisionFile(void)
{
    char input_filename[4096];
    sprintf(input_filename, "gold/%s-%s-%lu-%dnm.gold", prefixes[GRID_ONE], prefixes[GRID_TWO], threshold, (int)(window_radius + 0.5));

    // open file
    FILE *fp = fopen(input_filename, "rb");
    if (!fp) return;

    int ndecisions;
    if (fread(&ndecisions, sizeof(int), 1, fp) != 1) return;

    for (int iv = 0; iv < ndecisions; ++iv) {
        int decision;
        if (fread(&decision, sizeof(int), 1, fp) != 1) return;
        decisions.push_back(DECISION(decision));
        ++candidate_index;
    }

    // close file
    fclose(fp);
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
            for (int ix = 0; ix < grid_size[RN_X]; ++ix, ++iv) {
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
        printf("  Created random access variables in %0.2f seconds\n", start_time.Elapsed());
    }
}




static R3Point GridToWorld(R3Point grid_position, int index)
{
    // get the world box for this grid
    R3Box world_box = meta_data[index].world_box;

    return R3Point(world_box.XMin() + grid_position.X(), world_box.YMin() + grid_position.Y(), world_box.ZMin() + grid_position.Z());
}



static void GenerateCombinedImage(void)
{
    R3Box world_box = meta_data[GRID_ONE].world_box;
    for (int ig = 0; ig < NGRIDS; ++ig) 
        world_box.Union(meta_data[ig].world_box);


    // create a combined image grid
    int xwidth = (int)(world_box.XLength() + 0.5);
    int ywidth = (int)(world_box.YLength() + 0.5);
    int zwidth = (int)(world_box.ZLength() + 0.5);

    image_grid = new R3Grid(xwidth, ywidth, zwidth);
    image_grid->SetWorldToGridTransformation(world_box);

    int xmin = (int)(world_box.XMin() + 0.5);
    int ymin = (int)(world_box.YMin() + 0.5);
    int zmin = (int)(world_box.ZMin() + 0.5);
    int xmax = (int)(world_box.XMax() + 0.5);
    int ymax = (int)(world_box.YMax() + 0.5);
    int zmax = (int)(world_box.ZMax() + 0.5);

    for (int iz = zmin; iz < zmax; ++iz) {
        for (int iy = ymin; iy < ymax; ++iy) {
            for (int ix = xmin; ix < xmax; ++ix) {
                // get the position from this world coordinate
                R3Point position_one = image_grids[GRID_ONE]->GridPosition(ix, iy, iz);

                // if the position is not in image one, go to image two
                if (R3Intersects(image_grids[GRID_ONE]->GridBox(), position_one)) {
                    // get the image value
                    RNScalar image_value = image_grids[GRID_ONE]->GridValue(position_one);

                    // set the grid value
                    image_grid->SetGridValue(ix - xmin, iy - ymin, iz - zmin, image_value);
                }
                else {
                    // get the position from this world coordinate
                    R3Point position_two = image_grids[GRID_TWO]->GridPosition(ix, iy, iz);

                    // get the image value
                    RNScalar image_value = image_grids[GRID_TWO]->GridValue(position_two);

                    // set the grid value
                    image_grid->SetGridValue(ix - xmin, iy - ymin, iz - zmin, image_value);
                }
            }
        }
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
        glVertex3f(world_position.X(), world_position.Y(), world_position.Z());
    }
    glEnd();
}



static void DrawSlice()
{
    // reset color pallete to white
    RNLoadRgb(1.0, 1.0, 1.0);

    // draw the selected slice
    image_grid->DrawSlice(projection_dim, selected_slice_index[projection_dim]);
}



////////////////////////////////////////////////////////////////////////
// GLUT interface functions
////////////////////////////////////////////////////////////////////////

void GLUTStop(void)
{
    // save all of the decisions
    char output_filename[4096];
    sprintf(output_filename, "gold/%s-%s-%lu-%dnm.gold", prefixes[GRID_ONE], prefixes[GRID_TWO], threshold, (int)(window_radius + 0.5));

    // open file
    FILE *fp = fopen(output_filename, "wb");
    if (!fp) { fprintf(stderr, "Failed to write to %s\n", output_filename); }
    else {
        // write the number of decisions
        int ndecisions = decisions.size();
        fwrite(&ndecisions, sizeof(int), 1, fp);

        for (int iv = 0; iv < ndecisions; ++iv) {
            int decision = decisions[iv];
            fwrite(&decision, sizeof(int), 1, fp);
        }

        // close file
        fclose(fp);   
    }


    // destroy window
    glutDestroyWindow(GLUTwindow);

    // delete the neuron data
    RNTime start_time;
    start_time.Read();

    for (int iv = 0; iv < NGRIDS; ++iv) {
        delete image_grids[iv];
    }
    delete image_grid;

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
    
    // draw the bounding box for this location
    RNLoadRgb(RNwhite_rgb);
    RNScalar radius[3] = { window_radius / resolution[RN_X], window_radius / resolution[RN_Y], window_radius / resolution[RN_Z] };
    R3Point center = candidates[candidate_index].center;
    R3Box(center.X() - radius[RN_X], center.Y() - radius[RN_Y], center.Z() - radius[RN_Z], center.X() + radius[RN_X], center.Y() + radius[RN_Y], center.Z() + radius[RN_Z]).Outline();

    // draw the slice if needed
    if (show_slice) DrawSlice();

    // pop the transformation
    transformation.Pop();

    // draw large bounding box bounding box
    if(show_bbox) {
        RNLoadRgb(RNwhite_rgb);
        world_box.Outline();

        // draw all bounding boxes
        for (int iv = 0; iv < NGRIDS; ++iv) {
            meta_data[iv].scaled_box.Outline();
        }
    }

    // write the candidate index
    char candidate_label[4096];
    sprintf(candidate_label, "Ebro Gold - %d: %lf\n", candidate_index, overlap_scores[candidate_index]);
    glutSetWindowTitle(candidate_label);

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
            downsample_rate -= 1.0;
            if (downsample_rate < 1) downsample_rate = 1.0;
            break;
        }

        case GLUT_KEY_RIGHT: {
            downsample_rate += 1.0;
            break;
        }

        case GLUT_KEY_UP: {
	    selected_slice_index[projection_dim] += 5;
	    if (selected_slice_index[projection_dim] >= image_grid->Resolution(projection_dim)) 
                selected_slice_index[projection_dim] = image_grid->Resolution(projection_dim) - 1;
            break;
        }

        case GLUT_KEY_DOWN: {
	    selected_slice_index[projection_dim] -= 5;
	    if (selected_slice_index[projection_dim] < 0)
                selected_slice_index[projection_dim] = 0;
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
        /* decision keys */
        case 'A':
        case 'a': {
            decisions.push_back(YES);
            ++candidate_index;
            break;
        }

        case 'D':
        case 'd': {
            decisions.push_back(NO);
            ++candidate_index;
            break;
        }

        case SPACEBAR: {
            decisions.push_back(UNDECIDED);
            ++candidate_index;
            break;
        }


        case 'B':
        case 'b': {
            show_bbox = 1 - show_bbox;
            break;
        }

        case 'W':
        case 'w': {
            show_slice = 1 - show_slice;
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

        case CTRLZ: {
            if (decisions.size()) {
                decisions.pop_back();
                --candidate_index;
            }
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
            else if (!strcmp(*argv, "-threshold")) { argv++; argc--; threshold = atoi(*argv); }
            else if (!strcmp(*argv, "-window_radius")) { argv++; argc--; window_radius = atoi(*argv); }
            else { fprintf(stderr, "Invalid program argument: %s\n", *argv); return 0; }
        } else {
            if (!prefixes[GRID_ONE]) { prefixes[GRID_ONE] = *argv; } 
            else if (!prefixes[GRID_TWO]) { prefixes[GRID_TWO] = *argv; }
            else { fprintf(stderr, "Invalid program argument: %s\n", *argv); return 0; }
        }
        argv++; argc--;
    }

    // error if there is no input name
    for (int iv = 0; iv < NGRIDS; ++iv) {
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

    for (int iv = 0; iv < NGRIDS; ++iv) {
        if (!ReadMetaData(prefixes[iv], iv)) exit(-1);
        if (!ReadVoxelGrids(iv)) exit(-1);
        LabelToIndexMapping(iv);

        // delete rhoana grids to save space
        delete rhoana_grids[iv];
    }

    // find the potential merge candidates
    if (!ReadOverlapCandidates()) exit(-1);

    // set world box
    world_box = meta_data[GRID_ONE].scaled_box;
    for (int iv = 1; iv < NGRIDS; ++iv) 
        world_box.Union(meta_data[iv].scaled_box);
    
    // get the transformation
    transformation = R3Affine(R4Matrix(resolution[RN_X], 0, 0, 0, 0, resolution[RN_Y], 0, 0, 0, 0, resolution[RN_Z], 0, 0, 0, 0, 1));

    // merge the two images
    GenerateCombinedImage();

    // read the existing decision file
    ReadDecisionFile();

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
