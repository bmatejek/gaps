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
static int benchmark = 0;
static const char* prefix = NULL;



// program variables

static double resolution[3] = { -1, -1, -1 };
static long grid_size[3] = { -1, -1, -1 };
static R3Affine transformation = R3null_affine;
static R3Viewer *viewer = NULL;
static R3Box world_box = R3null_box;



// voxel grids

static R3Grid *grid = NULL;



// skeleton variables

static std::vector<long> *thinning_skeletons = NULL;
static std::vector<long> *medial_skeletons = NULL;
static std::vector<long> *teaser_skeletons = NULL;



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



// display variables

static const int ncolor_opts = 6;
static const int nskeleton_opts = 4;
static const int nteaser_scales = 6;
static const int nteaser_buffers = 5;
static const int nastar_expansions = 9;
static const int nresolutions = 18;

static int show_vectors = 1;
static int show_bbox = 1;
static int segmentation_index = 1;
static int skeleton_type = 0;
static RNScalar downsample_rate = 2.0;
static int color_cycle = 0;

// variables that enable cycling through skeletons
static long teaser_scales[nteaser_scales] = { 7, 9, 11, 13, 15, 17 };
static long teaser_buffers[nteaser_buffers] = { 1, 2, 3, 4, 5 };
static long resolutions[nresolutions][3] = {
    { 30, 30, 30 },
    { 40, 40, 40 },
    { 50, 50, 50 },
    { 60, 60, 60 },
    { 70, 70, 70 },
    { 80, 80, 80 },
    { 90, 90, 90 },
    { 100, 100, 100 },
    { 110, 110, 110 },
    { 120, 120, 120 },
    { 130, 130, 130 },
    { 140, 140, 140 },
    { 150, 150, 150 },
    { 160, 160, 160 },
    { 170, 170, 170 },
    { 180, 180, 180 },
    { 190, 190, 190 },
    { 200, 200, 200 },
};
static long astar_expansions[nastar_expansions] = { 0, 11, 13, 15, 17, 19, 21, 23, 25 };
static short teaser_scale_index = 2;
static short teaser_buffer_index = 1;
static short resolution_index = 5;
static short astar_expansion_index = 0;

static long astar_expansion = astar_expansions[astar_expansion_index];
static long teaser_scale = teaser_scales[teaser_scale_index];
static long teaser_buffer = teaser_buffers[teaser_buffer_index];
static long *downsample_resolution = resolutions[resolution_index];

static std::unordered_map<long, R3Vector> *thinning_endpoint_vectors = NULL;
static std::unordered_map<long, R3Vector> *medial_endpoint_vectors = NULL;



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

    // skip world box
    if (!fgets(comment, 4096, fp)) return 0;
    if (!fgets(comment, 4096, fp)) return 0;

    // read in requisite information
    if (!fgets(comment, 4096, fp)) return 0;
    if (fscanf(fp, "%ldx%ldx%ld\n", &(meta_data.grid_size[IB_X]), &(meta_data.grid_size[IB_Y]), &(meta_data.grid_size[IB_Z])) != 3) return 0;

    // update the global resolution
    for (int dim = 0; dim <= 2; ++dim)
        resolution[dim] = meta_data.resolution[dim];

    meta_data.world_box = R3null_box;
    meta_data.scaled_box = R3null_box;

    // close the file
    fclose(fp);

    // return success
    return 1;
}



static int
ReadSkeletonData(void)
{
    if (thinning_skeletons) { delete[] thinning_skeletons; thinning_skeletons = NULL; }
    if (medial_skeletons) { delete[] medial_skeletons; medial_skeletons = NULL; }
    if (teaser_skeletons) { delete[] teaser_skeletons; teaser_skeletons = NULL; }

    printf("Resolution: (%ld, %ld, %ld)\n", downsample_resolution[IB_X], downsample_resolution[IB_Y], downsample_resolution[IB_Z]);
    printf("  TEASER:\n");
    printf("    Scale: %ld\n", teaser_scale);
    printf("    Buffer: %ld\n", teaser_buffer);
    printf("  Thinning + Medial Axis:\n");
    printf("    A*: %ld\n", astar_expansion);


    char input_filename[4096];
    if (benchmark) sprintf(input_filename, "benchmarks/skeleton/%s-thinning-%03ldx%03ldx%03ld-upsample-%02ld-skeleton.pts", prefix, downsample_resolution[IB_X], downsample_resolution[IB_Y], downsample_resolution[IB_Z], astar_expansion);
    else sprintf(input_filename, "skeletons/%s/thinning-%03ldx%03ldx%03ld-upsample-%02ld-skeleton.pts", prefix, downsample_resolution[IB_X], downsample_resolution[IB_Y], downsample_resolution[IB_Z], astar_expansion);

    FILE *fp = fopen(input_filename, "rb");
    long skeleton_maximum_segmentation;
    long input_grid_size[3];
    if (fp) {
        if (fread(&input_grid_size[IB_Z], sizeof(long), 1, fp) != 1) return 0;
        if (fread(&input_grid_size[IB_Y], sizeof(long), 1, fp) != 1) return 0;
        if (fread(&input_grid_size[IB_X], sizeof(long), 1, fp) != 1) return 0;
        if (fread(&skeleton_maximum_segmentation, sizeof(long), 1, fp) != 1) return 0;
        assert (input_grid_size[IB_Z] == grid_size[IB_Z] and input_grid_size[IB_Y] == grid_size[IB_Y] and input_grid_size[IB_X] == grid_size[IB_X]);
        assert (skeleton_maximum_segmentation == maximum_segmentation);

        thinning_skeletons = new std::vector<long>[maximum_segmentation];
        for (long iv = 0; iv < maximum_segmentation; ++iv) {
            thinning_skeletons[iv] = std::vector<long>();

            long nelements; 
            if (fread(&nelements, sizeof(long), 1, fp) != 1) return 0;
            for (long ie = 0; ie < nelements; ++ie) {
                long element;
                if (fread(&element, sizeof(long), 1, fp) != 1) return 0;

                thinning_skeletons[iv].push_back(element);
            }
        }  
        fclose(fp);
    }

    if (benchmark) sprintf(input_filename, "benchmarks/skeleton/%s-medial-axis-%03ldx%03ldx%03ld-upsample-%02ld-skeleton.pts", prefix, downsample_resolution[IB_X], downsample_resolution[IB_Y], downsample_resolution[IB_Z], astar_expansion);
    else sprintf(input_filename, "skeletons/%s/medial-axis-%03ldx%03ldx%03ld-upsample-%02ld-skeleton.pts", prefix, downsample_resolution[IB_X], downsample_resolution[IB_Y], downsample_resolution[IB_Z], astar_expansion);

    fp = fopen(input_filename, "rb");
    if (fp) {
        if (fread(&input_grid_size[IB_Z], sizeof(long), 1, fp) != 1) return 0;
        if (fread(&input_grid_size[IB_Y], sizeof(long), 1, fp) != 1) return 0;
        if (fread(&input_grid_size[IB_X], sizeof(long), 1, fp) != 1) return 0;
        if (fread(&skeleton_maximum_segmentation, sizeof(long), 1, fp) != 1) return 0;
        assert (input_grid_size[IB_Z] == grid_size[IB_Z] and input_grid_size[IB_Y] == grid_size[IB_Y] and input_grid_size[IB_X] == grid_size[IB_X]);
        assert (skeleton_maximum_segmentation == maximum_segmentation);

        medial_skeletons = new std::vector<long>[maximum_segmentation];
        for (long iv = 0; iv < maximum_segmentation; ++iv) {
            medial_skeletons[iv] = std::vector<long>();

            long nelements;
            if (fread(&nelements, sizeof(long), 1, fp) != 1) return 0;
            for (long ie = 0; ie < nelements; ++ie) {
                long element;
                if (fread(&element, sizeof(long), 1, fp) != 1) return 0;
                medial_skeletons[iv].push_back(element);
            }
        }
        fclose(fp);
    }

    if (benchmark) sprintf(input_filename, "benchmarks/skeleton/%s-teaser-%03ldx%03ldx%03ld-upsample-%02ld-%02ld-00-skeleton.pts", prefix, downsample_resolution[IB_X], downsample_resolution[IB_Y], downsample_resolution[IB_Z], teaser_scale, teaser_buffer);
    else sprintf(input_filename, "skeletons/%s/teaser-%03ldx%03ldx%03ld-upsample-%02ld-%02ld-00-skeleton.pts", prefix, downsample_resolution[IB_X], downsample_resolution[IB_Y], downsample_resolution[IB_Z], teaser_scale, teaser_buffer);

    fp = fopen(input_filename, "rb");
    if (fp) {
        if (fread(&input_grid_size[IB_Z], sizeof(long), 1, fp) != 1) return 0;
        if (fread(&input_grid_size[IB_Y], sizeof(long), 1, fp) != 1) return 0;
        if (fread(&input_grid_size[IB_X], sizeof(long), 1, fp) != 1) return 0;
        if (fread(&skeleton_maximum_segmentation, sizeof(long), 1, fp) != 1) return 0;
        assert (input_grid_size[IB_Z] == grid_size[IB_Z] and input_grid_size[IB_Y] == grid_size[IB_Y] and input_grid_size[IB_X] == grid_size[IB_X]);
        assert (skeleton_maximum_segmentation == maximum_segmentation);

        teaser_skeletons = new std::vector<long>[maximum_segmentation];
        for (long iv = 0; iv < maximum_segmentation; ++iv) {
            teaser_skeletons[iv] = std::vector<long>();

            long nelements;
            if (fread(&nelements, sizeof(long), 1, fp) != 1) return 0;
            for (long ie = 0; ie < nelements; ++ie) {
                long element;
                if (fread(&element, sizeof(long), 1, fp) != 1) return 0;
                teaser_skeletons[iv].push_back(element);
            }
        }
        fclose(fp);
    }

    if (astar_expansion) return 1;

    if (benchmark) sprintf(input_filename, "benchmarks/skeleton/%s-thinning-%03ldx%03ldx%03ld-endpoint-vectors.vec", prefix, downsample_resolution[IB_X], downsample_resolution[IB_Y], downsample_resolution[IB_Z]);
    else sprintf(input_filename, "skeletons/%s/thinning-%03ldx%03ldx%03ld-endpoint-vectors.vec", prefix, downsample_resolution[IB_X], downsample_resolution[IB_Y], downsample_resolution[IB_Z]);

    fp = fopen(input_filename, "rb");
    if (fp) {
        if (fread(&input_grid_size[IB_Z], sizeof(long), 1, fp) != 1) return 0;
        if (fread(&input_grid_size[IB_Y], sizeof(long), 1, fp) != 1) return 0;
        if (fread(&input_grid_size[IB_X], sizeof(long), 1, fp) != 1) return 0;
        if (fread(&skeleton_maximum_segmentation, sizeof(long), 1, fp) != 1) return 0;
        assert (input_grid_size[IB_Z] == grid_size[IB_Z] and input_grid_size[IB_Y] == grid_size[IB_Y] and input_grid_size[IB_X] == grid_size[IB_X]);
        assert (skeleton_maximum_segmentation == maximum_segmentation);

        thinning_endpoint_vectors = new std::unordered_map<long, R3Vector>[maximum_segmentation];
        for (long is = 0; is < maximum_segmentation; ++is)
            thinning_endpoint_vectors[is] = std::unordered_map<long, R3Vector>();

        for (long is = 0; is < maximum_segmentation; ++is) {
            long nendpoints;
            if (fread(&nendpoints, sizeof(long), 1, fp) != 1) return 0;
            for (long ie = 0; ie < nendpoints; ++ie) {
                long endpoint;
                double vx, vy, vz;
                if (fread(&endpoint, sizeof(long), 1, fp) != 1) return 0;
                if (fread(&vz, sizeof(double), 1, fp) != 1) return 0;
                if (fread(&vy, sizeof(double), 1, fp) != 1) return 0;
                if (fread(&vx, sizeof(double), 1, fp) != 1) return 0;

                thinning_endpoint_vectors[is][endpoint] = R3Vector(vx, vy, vz);
            }
        }

        fclose(fp);
    }

    if (benchmark) sprintf(input_filename, "benchmarks/skeleton/%s-medial-axis-%03ldx%03ldx%03ld-endpoint-vectors.vec", prefix, downsample_resolution[IB_X], downsample_resolution[IB_Y], downsample_resolution[IB_Z]);
    else sprintf(input_filename, "skeletons/%s/medial-axis-%03ldx%03ldx%03ld-endpoint-vectors.vec", prefix, downsample_resolution[IB_X], downsample_resolution[IB_Y], downsample_resolution[IB_Z]);

    fp = fopen(input_filename, "rb");
    if (fp) {
        if (fread(&input_grid_size[IB_Z], sizeof(long), 1, fp) != 1) return 0;
        if (fread(&input_grid_size[IB_Y], sizeof(long), 1, fp) != 1) return 0;
        if (fread(&input_grid_size[IB_X], sizeof(long), 1, fp) != 1) return 0;
        if (fread(&skeleton_maximum_segmentation, sizeof(long), 1, fp) != 1) return 0;
        assert (input_grid_size[IB_Z] == grid_size[IB_Z] and input_grid_size[IB_Y] == grid_size[IB_Y] and input_grid_size[IB_X] == grid_size[IB_X]);
        assert (skeleton_maximum_segmentation == maximum_segmentation);

        medial_endpoint_vectors = new std::unordered_map<long, R3Vector>[maximum_segmentation];
        for (long is = 0; is < maximum_segmentation; ++is)
            medial_endpoint_vectors[is] = std::unordered_map<long, R3Vector>();
  
        for (long is = 0; is < maximum_segmentation; ++is) {
            long nendpoints;
            if (fread(&nendpoints, sizeof(long), 1, fp) != 1) return 0;
            for (long ie = 0; ie < nendpoints; ++ie) {
                long endpoint;
                double vx, vy, vz;
                if (fread(&endpoint, sizeof(long), 1, fp) != 1) return 0;
                if (fread(&vz, sizeof(double), 1, fp) != 1) return 0;
                if (fread(&vy, sizeof(double), 1, fp) != 1) return 0;
                if (fread(&vx, sizeof(double), 1, fp) != 1) return 0;

                medial_endpoint_vectors[is][endpoint] = R3Vector(vx, vy, vz);
            }
        }      

        fclose(fp);
    }

    // return success
    return 1;
}



static int ReadData(void)
{
    // start statistics
    RNTime start_time;
    start_time.Read();

    if (benchmark) {
        R3Grid **gold_grids = RNReadH5File(meta_data.gold_filename, meta_data.gold_dataset);
        grid = gold_grids[0];
        delete[] gold_grids;
    }
    else {
        R3Grid **rhoana_grids = RNReadH5File(meta_data.rhoana_filename, meta_data.rhoana_dataset);
        grid = rhoana_grids[0];
        delete[] rhoana_grids;
    }
    
    grid_size[IB_X] = grid->XResolution();
    grid_size[IB_Y] = grid->YResolution();
    grid_size[IB_Z] = grid->ZResolution();

    // get the maximum values for each grid
    maximum_segmentation = (long) (grid->Maximum() + 0.5) + 1;

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
    return iz * grid_size[IB_X] * grid_size[IB_Y] + iy * grid_size[IB_X] + ix;
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
    for (long iv = 0; iv < maximum_segmentation; ++iv)
        segmentations[iv] = std::vector<long>();

    // go through all voxels to see if it belongs to the boundary
    for (int iz = 1; iz < grid_size[IB_Z] - 1; ++iz) {
        for (int iy = 1; iy < grid_size[IB_Y] - 1; ++iy) {
            for (int ix = 1; ix < grid_size[IB_X] - 1; ++ix) {
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

    if (print_verbose) {
        printf("Preprocessing...\n");
        printf("  Time = %0.2f seconds\n", start_time.Elapsed());
        printf("  Maximum Segment = %ld\n", maximum_segmentation);
    }
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



static void DrawSkeleton(int segment_index)
{
    std::vector<long> *skeletons;
    if (skeleton_type == 0) skeletons = thinning_skeletons;
    else if (skeleton_type == 1) skeletons = medial_skeletons;
    else if (skeleton_type == 2) skeletons = teaser_skeletons;
    else return;
    if (!skeletons) return;

    std::unordered_map<long, R3Vector> *endpoint_vectors = NULL;
    if (skeleton_type == 0) endpoint_vectors = thinning_endpoint_vectors;
    else if (skeleton_type == 1) endpoint_vectors = medial_endpoint_vectors;

    // sizes for skeleton joints
    const double joint_size = 4;
    const double endpoint_size = 60;
    const double line_size = 2;
    const double line_length = 250;

    transformation.Push();

    glPointSize(joint_size);
    glBegin(GL_POINTS);
    for (unsigned long is = 0; is < skeletons[segment_index].size(); ++is) {
        // get the coordinates from the linear index
        long iv = skeletons[segment_index][is];
        bool endpoint = (iv < 0);
        if (endpoint) continue;
        long ix, iy, iz;
        IndexToIndices(iv, ix, iy, iz);
        glVertex3f(ix, iy, iz);
    }
    glEnd();
    glPointSize(1.0);

    transformation.Pop();

    for (unsigned long is = 0; is < skeletons[segment_index].size(); ++is) {
        // get the coordinates from the linear index
        long iv = skeletons[segment_index][is];
        bool endpoint = (iv < 0);
        if (!endpoint) continue;
        iv = -1 * iv;
        long ix, iy, iz;
        IndexToIndices(iv, ix, iy, iz);
        
        // convert to world coordinates since transformation is popped
        ix = resolution[IB_X] * ix;
        iy = resolution[IB_Y] * iy;
        iz = resolution[IB_Z] * iz;

        R3Sphere(R3Point(ix, iy, iz), endpoint_size).Draw();

        if (show_vectors) {
            // draw the vector if it exists
            glLineWidth(line_size);
            if (endpoint_vectors) {
                R3Vector vector = endpoint_vectors[segment_index][iv];

                // draw line
                glBegin(GL_LINES);
                glVertex3f(ix, iy, iz);
                glVertex3f(ix + line_length * vector.X(), iy + line_length * vector.Y(), iz + line_length * vector.Z());
                glEnd();

                R3Sphere(R3Point(ix, iy, iz) + line_length * vector, endpoint_size / 2.0).Draw();
            }
            glLineWidth(1.0);
        }
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
        RNLoadRgb(RNRgb(not background_color[0], not background_color[1], not background_color[2]));
        world_box.Outline();
    }

    // draw machine labels and skeletons
    glPointSize(1);
    RNLoadRgb(Color(segmentation_index));
    DrawSegment(segmentation_index);
    RNLoadRgb(RNRgb(1.0 - background_color[0], 1.0 - background_color[1], 1.0 - background_color[2]));
    DrawSkeleton(segmentation_index);

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
        case 'A':
        case 'a': {
            astar_expansion_index = (astar_expansion_index + 1) % nastar_expansions;

            astar_expansion = astar_expansions[astar_expansion_index];

            ReadSkeletonData();

            break;
        }

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

        case 'K':
        case 'k': {
            skeleton_type = (skeleton_type + 1) % nskeleton_opts;

            if (skeleton_type == 0) {
                printf("Thinning\n");
                printf("  A*: %0.2lf\n", astar_expansion / 10.0);
            }
            else if (skeleton_type == 1) {
                printf("Medial Axis\n");
                printf("  A*: %0.2lf\n", astar_expansion / 10.0);
            }
            else if (skeleton_type == 2) {
                printf("TEASER:\n");
                printf("  Scale: %0.2lf\n", teaser_scale / 10.0);
                printf("  Buffer: %ld\n", teaser_buffer);
            }

            break;
        }

        case 'R': 
        case 'r': {
            resolution_index = (resolution_index + 1) % nresolutions;
            
            downsample_resolution = resolutions[resolution_index];

            ReadSkeletonData();

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
            teaser_scale_index = (teaser_scale_index + 1) % nteaser_scales;

            teaser_scale = teaser_scales[teaser_scale_index];

            ReadSkeletonData();

            break;
        }

        case 'T': 
        case 't': {
            teaser_buffer_index = (teaser_buffer_index + 1) % nteaser_buffers;

            teaser_buffer = teaser_buffers[teaser_buffer_index];

            ReadSkeletonData();

            break;
        }

        case 'V': 
        case 'v': {
            show_vectors = 1 - show_vectors;
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
            else if (!strcmp(*argv, "-benchmark")) benchmark = 1;
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
    if (!ReadSkeletonData()) exit(-1);

    // get all of the preprocessing time
    Preprocessing();

    if (segmentation_index < 0) segmentation_index = 0;
    if (segmentation_index > maximum_segmentation - 1) segmentation_index = maximum_segmentation - 1;

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