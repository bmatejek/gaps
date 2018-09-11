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



// program variables

static double resolution[3] = { -1, -1, -1 };
static long grid_size[3] = { -1, -1, -1 };
static R3Affine transformation = R3null_affine;
static R3Viewer *viewer = NULL;
static R3Box world_box = R3null_box;



// voxel grids

static R3Grid *grid = NULL;



// skeleton variables

static std::vector<long> *thinning_endpoints = NULL;
static std::vector<long> *thinning_skeletons = NULL;

static std::vector<long> *medial_endpoints = NULL;
static std::vector<long> *medial_skeletons = NULL;

static std::vector<long> *teaser_endpoints = NULL;
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

static std::vector<long> *golds = NULL;
static long maximum_gold = -1;



// display variables

static const int ncolor_opts = 6;
static const int nskeleton_opts = 4;
static const int nteaser_scales = 6;
static const int nteaser_buffers = 5;
static const int nastar_expansions = 9;
static const int nresolutions = 18;

static int show_bbox = 1;
static int show_matches = 0;
static int gold_index = 1;
static int example_index = 0;
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


static std::set<long> *medial_auto_hits = NULL;
static std::set<long> *medial_gold_hits = NULL;
static std::set<long> *medial_auto_misses = NULL;
static std::set<long> *medial_gold_misses = NULL;

static std::set<long> *teaser_auto_hits = NULL;
static std::set<long> *teaser_gold_hits = NULL;
static std::set<long> *teaser_auto_misses = NULL;
static std::set<long> *teaser_gold_misses = NULL;

static std::set<long> *thinning_auto_hits = NULL;
static std::set<long> *thinning_gold_hits = NULL;
static std::set<long> *thinning_auto_misses = NULL;
static std::set<long> *thinning_gold_misses = NULL;


// variables for ground truth skeletons

static std::vector<R3Point> gold_endpoints = std::vector<R3Point>();
static const int cutoff = 500;
static long largest_segments[cutoff];



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
    meta_data.world_box = R3null_box;
    meta_data.scaled_box = R3null_box;

    // close the file
    fclose(fp);

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



void ReadSkeletonEndpoints(void)
{
    char input_filename[4096];
    sprintf(input_filename, "benchmarks/skeleton/%s/skeleton-endpoints-%05ld.pts", prefix, largest_segments[example_index]);

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

            gold_endpoints.push_back(R3Point(xpoint, ypoint, zpoint));
        }
        fclose(fp);
    }   
}



void WriteSkeletonEndpoints(void)
{
    if (not gold_endpoints.size()) return;

    char output_filename[4096];
    sprintf(output_filename, "benchmarks/skeleton/%s/skeleton-endpoints-%05ld.pts", prefix, largest_segments[example_index]);

    FILE *fp = fopen(output_filename, "wb");
    if (!fp) { fprintf(stderr, "Failed to write %s\n", output_filename); }

    long nendpoints = gold_endpoints.size();
    fwrite(&nendpoints, sizeof(long), 1, fp);
    for (long ie = 0; ie < nendpoints; ++ie) {
        R3Point endpoint = gold_endpoints[ie];
        long xpoint = (long) (endpoint.X() + 0.5);
        long ypoint = (long) (endpoint.Y() + 0.5);
        long zpoint = (long) (endpoint.Z() + 0.5);
        fwrite(&zpoint, sizeof(long), 1, fp);
        fwrite(&ypoint, sizeof(long), 1, fp);
        fwrite(&xpoint, sizeof(long), 1, fp);
    }

    fclose(fp);
}



static int
ReadSkeletonData(void)
{
    if (thinning_endpoints) { delete[] thinning_endpoints; thinning_endpoints = NULL; }
    if (thinning_skeletons) { delete[] thinning_skeletons; thinning_skeletons = NULL; }

    if (medial_endpoints) { delete[] medial_endpoints; medial_endpoints = NULL; }
    if (medial_skeletons) { delete[] medial_skeletons; medial_skeletons = NULL; }

    if (teaser_endpoints) { delete[] teaser_endpoints; teaser_endpoints = NULL; }
    if (teaser_skeletons) { delete[] teaser_skeletons; teaser_skeletons = NULL; }

    if (medial_auto_hits) { delete[] medial_auto_hits; medial_auto_hits = NULL; }
    if (medial_gold_hits) { delete[] medial_gold_hits; medial_gold_hits = NULL; }
    if (medial_auto_misses) { delete[] medial_auto_misses; medial_auto_misses = NULL; }
    if (medial_gold_misses) { delete[] medial_gold_misses; medial_gold_misses = NULL; }
    
    if (teaser_auto_hits) { delete[] teaser_auto_hits; teaser_auto_hits = NULL; }
    if (teaser_gold_hits) { delete[] teaser_gold_hits; teaser_gold_hits = NULL; }
    if (teaser_auto_misses) { delete[] teaser_auto_misses; teaser_auto_misses = NULL; }
    if (teaser_gold_misses) { delete[] teaser_auto_misses; teaser_auto_misses = NULL; }
    
    if (thinning_auto_hits) { delete[] thinning_auto_hits; thinning_auto_hits = NULL; }
    if (thinning_gold_hits) { delete[] thinning_gold_hits; thinning_gold_hits = NULL; }
    if (thinning_auto_misses) { delete[] thinning_auto_misses; thinning_auto_misses = NULL; }
    if (thinning_gold_misses) { delete[] thinning_gold_misses; thinning_gold_misses = NULL; }



    printf("Resolution: (%ld, %ld, %ld)\n", downsample_resolution[IB_X], downsample_resolution[IB_Y], downsample_resolution[IB_Z]);
    printf("  TEASER:\n");
    printf("    Scale: %ld\n", teaser_scale);
    printf("    Buffer: %ld\n", teaser_buffer);
    printf("  Thinning + Medial Axis:\n");
    printf("    A*: %ld\n", astar_expansion);


    char input_filename[4096];
    sprintf(input_filename, "benchmarks/skeleton/%s-thinning-%03ldx%03ldx%03ld-upsample-%02ld-skeleton.pts", prefix, downsample_resolution[IB_X], downsample_resolution[IB_Y], downsample_resolution[IB_Z], astar_expansion);

    FILE *fp = fopen(input_filename, "rb");
    long skeleton_maximum_gold;
    long input_grid_size[3];
    if (fp) {
        if (fread(&input_grid_size[IB_Z], sizeof(long), 1, fp) != 1) return 0;
        if (fread(&input_grid_size[IB_Y], sizeof(long), 1, fp) != 1) return 0;
        if (fread(&input_grid_size[IB_X], sizeof(long), 1, fp) != 1) return 0;
        if (fread(&skeleton_maximum_gold, sizeof(long), 1, fp) != 1) return 0;
        assert (input_grid_size[IB_Z] == grid_size[IB_Z] and input_grid_size[IB_Y] == grid_size[IB_Y] and input_grid_size[IB_X] == grid_size[IB_X]);
        assert (skeleton_maximum_gold == maximum_gold);

        thinning_endpoints = new std::vector<long>[maximum_gold];
        thinning_skeletons = new std::vector<long>[maximum_gold];
        for (long iv = 0; iv < maximum_gold; ++iv) {
            thinning_endpoints[iv] = std::vector<long>();
            thinning_skeletons[iv] = std::vector<long>();

            long nelements; 
            if (fread(&nelements, sizeof(long), 1, fp) != 1) return 0;
            for (long ie = 0; ie < nelements; ++ie) {
                long element;
                if (fread(&element, sizeof(long), 1, fp) != 1) return 0;

                if (element < 0) thinning_endpoints[iv].push_back(-1 * element);
                else thinning_skeletons[iv].push_back(element);
            }
        }  
        fclose(fp);
    }

    char matching_filename[4096];
    sprintf(matching_filename, "benchmarks/skeleton/matchings/%s-thinning-%03ldx%03ldx%03ld-%02ld-matching-pairs.pts", prefix, downsample_resolution[IB_X], downsample_resolution[IB_Y], downsample_resolution[IB_Z], astar_expansion);

    FILE *matching_fp = fopen(matching_filename, "rb");
    if (matching_fp) {
        long matching_maximum_gold;
        if (fread(&matching_maximum_gold, sizeof(long), 1, matching_fp) != 1) return 0;
        assert (matching_maximum_gold == maximum_gold);

        thinning_gold_hits = new std::set<long>[maximum_gold];
        thinning_auto_hits = new std::set<long>[maximum_gold];
        thinning_gold_misses = new std::set<long>[maximum_gold];
        thinning_auto_misses = new std::set<long>[maximum_gold];

        for (long iv = 0; iv < maximum_gold; ++iv) {
            long nmatches;
            if (fread(&nmatches, sizeof(long), 1, matching_fp) != 1) return 0;

            std::set<long> gold_matches = std::set<long>();
            std::set<long> auto_matches = std::set<long>();

            for (long ie = 0; ie < nmatches; ++ie) {
                long gold_index;
                long auto_index;
                if (fread(&gold_index, sizeof(long), 1, matching_fp) != 1) return 0;
                if (fread(&auto_index, sizeof(long), 1, matching_fp) != 1) return 0;

                gold_matches.insert(gold_index);
                auto_matches.insert(auto_index);
            }

            thinning_gold_hits[iv] = std::set<long>();
            thinning_auto_hits[iv] = std::set<long>();
            thinning_gold_misses[iv] = std::set<long>();
            thinning_auto_misses[iv] = std::set<long>();

            for (unsigned long ie = 0; ie < gold_endpoints.size(); ++ie) {
                if (gold_matches.find(ie) == gold_matches.end()) thinning_gold_misses[iv].insert(ie);
                else thinning_gold_hits[iv].insert(ie);
            }
            for (unsigned long ie = 0; ie < thinning_endpoints[iv].size(); ++ie) {
                if (auto_matches.find(ie) == auto_matches.end()) thinning_auto_misses[iv].insert(ie);
                else thinning_auto_hits[iv].insert(ie);
            }
        }

        fclose(matching_fp);
    }



    sprintf(input_filename, "benchmarks/skeleton/%s-medial-axis-%03ldx%03ldx%03ld-upsample-%02ld-skeleton.pts", prefix, downsample_resolution[IB_X], downsample_resolution[IB_Y], downsample_resolution[IB_Z], astar_expansion);

    fp = fopen(input_filename, "rb");
    if (fp) {
        if (fread(&input_grid_size[IB_Z], sizeof(long), 1, fp) != 1) return 0;
        if (fread(&input_grid_size[IB_Y], sizeof(long), 1, fp) != 1) return 0;
        if (fread(&input_grid_size[IB_X], sizeof(long), 1, fp) != 1) return 0;
        if (fread(&skeleton_maximum_gold, sizeof(long), 1, fp) != 1) return 0;
        assert (input_grid_size[IB_Z] == grid_size[IB_Z] and input_grid_size[IB_Y] == grid_size[IB_Y] and input_grid_size[IB_X] == grid_size[IB_X]);
        assert (skeleton_maximum_gold == maximum_gold);

        medial_endpoints = new std::vector<long>[maximum_gold];
        medial_skeletons = new std::vector<long>[maximum_gold];
        for (long iv = 0; iv < maximum_gold; ++iv) {
            medial_endpoints[iv] = std::vector<long>();
            medial_skeletons[iv] = std::vector<long>();

            long nelements;
            if (fread(&nelements, sizeof(long), 1, fp) != 1) return 0;
            for (long ie = 0; ie < nelements; ++ie) {
                long element;
                if (fread(&element, sizeof(long), 1, fp) != 1) return 0;
                if (element < 0) medial_endpoints[iv].push_back(-1 * element);
                else medial_skeletons[iv].push_back(element);
            }
        }
        fclose(fp);
    }


    sprintf(matching_filename, "benchmarks/skeleton/matchings/%s-medial-axis-%03ldx%03ldx%03ld-%02ld-matching-pairs.pts", prefix, downsample_resolution[IB_X], downsample_resolution[IB_Y], downsample_resolution[IB_Z], astar_expansion);

    matching_fp = fopen(matching_filename, "rb");
    if (matching_fp) {
        long matching_maximum_gold;
        if (fread(&matching_maximum_gold, sizeof(long), 1, matching_fp) != 1) return 0;
        assert (matching_maximum_gold == maximum_gold);

        medial_gold_hits = new std::set<long>[maximum_gold];
        medial_auto_hits = new std::set<long>[maximum_gold];
        medial_gold_misses = new std::set<long>[maximum_gold];
        medial_auto_misses = new std::set<long>[maximum_gold];

        for (long iv = 0; iv < maximum_gold; ++iv) {
            long nmatches;
            if (fread(&nmatches, sizeof(long), 1, matching_fp) != 1) return 0;

            std::set<long> gold_matches = std::set<long>();
            std::set<long> auto_matches = std::set<long>();

            for (long ie = 0; ie < nmatches; ++ie) {
                long gold_index;
                long auto_index;
                if (fread(&gold_index, sizeof(long), 1, matching_fp) != 1) return 0;
                if (fread(&auto_index, sizeof(long), 1, matching_fp) != 1) return 0;

                gold_matches.insert(gold_index);
                auto_matches.insert(auto_index);
            }

            medial_gold_hits[iv] = std::set<long>();
            medial_auto_hits[iv] = std::set<long>();
            medial_gold_misses[iv] = std::set<long>();
            medial_auto_misses[iv] = std::set<long>();

            for (unsigned long ie = 0; ie < gold_endpoints.size(); ++ie) {
                if (gold_matches.find(ie) == gold_matches.end()) medial_gold_misses[iv].insert(ie);
                else medial_gold_hits[iv].insert(ie);
            }
            for (unsigned long ie = 0; ie < medial_endpoints[iv].size(); ++ie) {
                if (auto_matches.find(ie) == auto_matches.end()) medial_auto_misses[iv].insert(ie);
                else medial_auto_hits[iv].insert(ie);
            }
        }

        fclose(matching_fp);
    }




    sprintf(input_filename, "benchmarks/skeleton/%s-teaser-%03ldx%03ldx%03ld-upsample-%02ld-%02ld-00-skeleton.pts", prefix, downsample_resolution[IB_X], downsample_resolution[IB_Y], downsample_resolution[IB_Z], teaser_scale, teaser_buffer);

    fp = fopen(input_filename, "rb");
    if (fp) {
        if (fread(&input_grid_size[IB_Z], sizeof(long), 1, fp) != 1) return 0;
        if (fread(&input_grid_size[IB_Y], sizeof(long), 1, fp) != 1) return 0;
        if (fread(&input_grid_size[IB_X], sizeof(long), 1, fp) != 1) return 0;
        if (fread(&skeleton_maximum_gold, sizeof(long), 1, fp) != 1) return 0;
        assert (input_grid_size[IB_Z] == grid_size[IB_Z] and input_grid_size[IB_Y] == grid_size[IB_Y] and input_grid_size[IB_X] == grid_size[IB_X]);
        assert (skeleton_maximum_gold == maximum_gold);

        teaser_endpoints = new std::vector<long>[maximum_gold];
        teaser_skeletons = new std::vector<long>[maximum_gold];
        for (long iv = 0; iv < maximum_gold; ++iv) {
            teaser_endpoints[iv] = std::vector<long>();
            teaser_skeletons[iv] = std::vector<long>();

            long nelements;
            if (fread(&nelements, sizeof(long), 1, fp) != 1) return 0;
            for (long ie = 0; ie < nelements; ++ie) {
                long element;
                if (fread(&element, sizeof(long), 1, fp) != 1) return 0;
                
                if (element < 0) teaser_endpoints[iv].push_back(-1 * element);
                else teaser_skeletons[iv].push_back(element);

            }
        }
        fclose(fp);
    }

    sprintf(matching_filename, "benchmarks/skeleton/matchings/%s-teaser-%03ldx%03ldx%03ld-%02ld-%02ld-00-matching-pairs.pts", prefix, downsample_resolution[IB_X], downsample_resolution[IB_Y], downsample_resolution[IB_Z], teaser_scale, teaser_buffer);

    matching_fp = fopen(matching_filename, "rb");
    if (matching_fp) {
        long matching_maximum_gold;
        if (fread(&matching_maximum_gold, sizeof(long), 1, matching_fp) != 1) return 0;
        assert (matching_maximum_gold == maximum_gold);

        teaser_gold_hits = new std::set<long>[maximum_gold];
        teaser_auto_hits = new std::set<long>[maximum_gold];
        teaser_gold_misses = new std::set<long>[maximum_gold];
        teaser_auto_misses = new std::set<long>[maximum_gold];

        for (long iv = 0; iv < maximum_gold; ++iv) {
            long nmatches;
            if (fread(&nmatches, sizeof(long), 1, matching_fp) != 1) return 0;

            std::set<long> gold_matches = std::set<long>();
            std::set<long> auto_matches = std::set<long>();

            for (long ie = 0; ie < nmatches; ++ie) {
                long gold_index;
                long auto_index;
                if (fread(&gold_index, sizeof(long), 1, matching_fp) != 1) return 0;
                if (fread(&auto_index, sizeof(long), 1, matching_fp) != 1) return 0;

                gold_matches.insert(gold_index);
                auto_matches.insert(auto_index);
            }

            teaser_gold_hits[iv] = std::set<long>();
            teaser_auto_hits[iv] = std::set<long>();
            teaser_gold_misses[iv] = std::set<long>();
            teaser_auto_misses[iv] = std::set<long>();

            for (unsigned long ie = 0; ie < gold_endpoints.size(); ++ie) {
                if (gold_matches.find(ie) == gold_matches.end()) teaser_gold_misses[iv].insert(ie);
                else teaser_gold_hits[iv].insert(ie);
            }
            for (unsigned long ie = 0; ie < teaser_endpoints[iv].size(); ++ie) {
                if (auto_matches.find(ie) == auto_matches.end()) teaser_auto_misses[iv].insert(ie);
                else teaser_auto_hits[iv].insert(ie);
            }
        }

        fclose(matching_fp);
    }



    // return success
    return 1;
}



static int ReadData(void)
{
    // start statistics
    RNTime start_time;
    start_time.Read();

    R3Grid **gold_grids = RNReadH5File(meta_data.gold_filename, meta_data.gold_dataset);
    grid = gold_grids[0];
    delete[] gold_grids;
    
    grid_size[IB_X] = grid->XResolution();
    grid_size[IB_Y] = grid->YResolution();
    grid_size[IB_Z] = grid->ZResolution();

    // get the maximum values for each grid
    maximum_gold = (long) (grid->Maximum() + 0.5) + 1;

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
    golds = new std::vector<long>[maximum_gold];
    for (long iv = 0; iv < maximum_gold; ++iv)
        golds[iv] = std::vector<long>();

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
                    golds[segment].push_back(iv);
                }
            }
        }
    }

    if (print_verbose) {
        printf("Preprocessing...\n");
        printf("  Time = %0.2f seconds\n", start_time.Elapsed());
        printf("  Maximum Segment = %ld\n", maximum_gold);
    }
}



////////////////////////////////////////////////////////////////////////
// Drawing utility functions
////////////////////////////////////////////////////////////////////////

static void DrawSegment(int segment_index)
{
    transformation.Push();

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

    transformation.Pop();
}



static void DrawSkeleton(int segment_index)
{
    std::vector<long> *endpoints, *skeletons;
    if (skeleton_type == 0) { endpoints = thinning_endpoints; skeletons = thinning_skeletons; }
    else if (skeleton_type == 1) { endpoints = medial_endpoints; skeletons = medial_skeletons; }
    else if (skeleton_type == 2) { endpoints = teaser_endpoints; skeletons = teaser_skeletons; }
    else return;

    if (!skeletons) return;
    if (!endpoints) return;

    // sizes for skeleton joints
    double joint_size = 3;
    double endpoint_size = 60;

    transformation.Push();

    glPointSize(joint_size);
    glBegin(GL_POINTS);
    for (unsigned long is = 0; is < skeletons[segment_index].size(); ++is) {
        // get the coordinates from the linear index
        long iv = skeletons[segment_index][is];
        long ix, iy, iz;
        IndexToIndices(iv, ix, iy, iz);
        glVertex3f(ix, iy, iz);
    }
    glEnd();

    transformation.Pop();

    std::set<long> hits, misses;
    if (skeleton_type == 0 and thinning_auto_hits and thinning_auto_misses) { hits = thinning_auto_hits[segment_index]; misses = thinning_auto_misses[segment_index]; }
    else if (skeleton_type == 1 and medial_auto_hits and medial_auto_misses) { hits = medial_auto_hits[segment_index]; misses = medial_auto_misses[segment_index]; }
    else if (skeleton_type == 2 and teaser_auto_hits and teaser_auto_misses) { hits = teaser_auto_hits[segment_index]; misses = teaser_auto_misses[segment_index]; }
    else return;

    std::set<long>::iterator it;
    if (show_matches) {
        RNLoadRgb(RNred_rgb);
        for (it = hits.begin(); it != hits.end(); ++it) {
            long index = *it;
            long iv = endpoints[segment_index][index];
            long ix, iy, iz;
            IndexToIndices(iv, ix, iy, iz);
            
            // convert to world coordinates since transformation is popped
            ix = resolution[IB_X] * ix;
            iy = resolution[IB_Y] * iy;
            iz = resolution[IB_Z] * iz;

            R3Sphere(R3Point(ix, iy, iz), endpoint_size).Draw();
        }  
    }

    RNLoadRgb(RNyellow_rgb);
    for (it = misses.begin(); it != misses.end(); ++it) {
        long index = *it;
        long iv = endpoints[segment_index][index];
        long ix, iy, iz;
        IndexToIndices(iv, ix, iy, iz);
        
        // convert to world coordinates since transformation is popped
        ix = resolution[IB_X] * ix;
        iy = resolution[IB_Y] * iy;
        iz = resolution[IB_Z] * iz;

        R3Sphere(R3Point(ix, iy, iz), endpoint_size).Draw();
    }
}



static void DrawGroundTruth(int segment_index)
{   
    double endpoint_size = 60;

    std::set<long> hits, misses;
    if (skeleton_type == 0 and thinning_gold_hits and thinning_gold_misses) { hits = thinning_gold_hits[segment_index]; misses = thinning_gold_misses[segment_index]; }
    else if (skeleton_type == 1 and medial_gold_hits and medial_gold_misses) { hits = medial_gold_hits[segment_index]; misses = medial_gold_misses[segment_index]; }
    else if (skeleton_type == 2 and teaser_gold_hits and teaser_gold_misses) { hits = teaser_gold_hits[segment_index]; misses = teaser_gold_misses[segment_index]; }
    else return;

    std::set<long>::iterator it;
    if (show_matches) {
        RNLoadRgb(RNgreen_rgb);
        for (it = hits.begin(); it != hits.end(); ++it) {
            long index = *it;
            R3Point endpoint = gold_endpoints[index];
            R3Sphere(R3Point(resolution[IB_X] * endpoint.X(), resolution[IB_Y] * endpoint.Y(), resolution[IB_Z] * endpoint.Z()), endpoint_size).Draw();
        }
    }
    RNLoadRgb(RNblue_rgb);
    for (it = misses.begin(); it != misses.end(); ++it) {
        long index = *it;
        R3Point endpoint = gold_endpoints[index];
        R3Sphere(R3Point(resolution[IB_X] * endpoint.X(), resolution[IB_Y] * endpoint.Y(), resolution[IB_Z] * endpoint.Z()), endpoint_size).Draw();
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
    RNLoadRgb(Color(gold_index));
    DrawSegment(gold_index);
    RNLoadRgb(RNRgb(1.0 - background_color[0], 1.0 - background_color[1], 1.0 - background_color[2]));
    DrawSkeleton(gold_index);
    DrawGroundTruth(gold_index);


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
            --example_index;
            if (example_index < 0)
                example_index = 0;

            gold_index = largest_segments[example_index];

            printf("Label: %d\n", gold_index);
            break;
        }

        case GLUT_KEY_RIGHT: {
            ++example_index;
            if (example_index >= cutoff)
                example_index = cutoff - 1;

            gold_index = largest_segments[example_index];

            printf("Label: %d\n", gold_index);
            break;    
        }
    }

    // upload the new endpoints
    gold_endpoints.clear();
    ReadSkeletonEndpoints();
    ReadSkeletonData();

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

            gold_endpoints.push_back((closest_point + second_closest_point) / 2);

            printf("  No. endpoints: %ld\n", gold_endpoints.size());

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

        case 'M':
        case 'm': {
            show_matches = 1 - show_matches;
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

        case 'R': 
        case 'r': {
            resolution_index = (resolution_index + 1) % nresolutions;
            
            downsample_resolution = resolutions[resolution_index];

            ReadSkeletonData();

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

        case DELETE: {
            R3Ray world_ray = viewer->WorldRay(x, y);

            RNLength closest_distance = DBL_MAX;
            unsigned long delete_index = 0;

            for (unsigned long iv = 0; iv < gold_endpoints.size(); ++iv) {
                R3Point location = gold_endpoints[iv];
                location.Transform(transformation);

                RNLength distance = R3Distance(world_ray, location);
                if (distance < closest_distance) {
                    closest_distance = distance;
                    delete_index = iv;
                }
            }

            gold_endpoints.erase(gold_endpoints.begin() + delete_index);

            printf("  No. endpoints: %ld\n", gold_endpoints.size());

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
    if (!ReadLargestSegments()) exit(-1);
    gold_index = largest_segments[example_index];
    ReadSkeletonEndpoints();
    if (!ReadSkeletonData()) exit(-1);

    // get all of the preprocessing time
    Preprocessing();

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