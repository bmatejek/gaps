// Source file for skeleton visualizer

// include files

#include "R3Graphics/R3Graphics.h"
#include "RNDataStructures/RNDataStructures.h"
#include <fstream>
#include <string>
#include <vector>
#include <algorithm>



// program arguments

static int print_debug = 0;
static int print_verbose = 0;
static const char* prefix = NULL;
static const char* gold_filename = NULL;
static const char* gold_dataset = NULL;
static const char* segmentation_filename = NULL;
static const char* segmentation_dataset = NULL;



// program variables

static RNScalar scaling[3] = { 1, 1, 5 };
static int resolution[3] = { -1, -1, -1 };



// voxel grids

static R3Grid* gold_grid = NULL;
static R3Grid* segmentation_grid = NULL;




// struct name definitions

struct SWCEntry;




// display machine label variables

static std::vector<int> *segmentations = NULL;
static std::vector<SWCEntry> *skeletons = NULL;
static std::vector<R3Point> *skeleton_endpoints = NULL;




// mapping functions

static RNBoolean** segmentation_neighbors = NULL;
static int* segmentation_to_gold = NULL;




// considered size restriction variables

static int normalized_distance = 100;
static RNScalar maximum_distance[3] = { -1, -1, -1 };
static int nclosest = 7;



////////////////////////////////////////////////////////////////////////
// Segment Functions
////////////////////////////////////////////////////////////////////////

struct SegmentPair {
    SegmentPair(void){ }
    SegmentPair(int segment_one, int segment_two, RNScalar distance, RNBoolean merge, R3Point center_point)
        : segment_one(segment_one)
        , segment_two(segment_two)
        , distance(distance)
        , merge(merge)
        , center_point(center_point)
        , bbox(R3null_box)
        , segment_one_points()
        , segment_two_points()
    {
    }

    int segment_one;
    int segment_two;
    RNScalar distance;
    RNBoolean merge;
    R3Point center_point;
    R3Box bbox;
    std::vector<int> segment_one_points;
    std::vector<int> segment_two_points;
};

std::vector<SegmentPair *> boundary_examples = std::vector<SegmentPair *>();
std::vector<RNScalar> *attributes = NULL;



////////////////////////////////////////////////////////////////////////
// Utility functions
////////////////////////////////////////////////////////////////////////

static RNScalar SquaredScaledDistance(R3Point p1, R3Point p2) 
{
    RNScalar xdiff = scaling[RN_X] * (p1.X() - p2.X());
    RNScalar ydiff = scaling[RN_Y] * (p1.Y() - p2.Y());
    RNScalar zdiff = scaling[RN_Z] * (p1.Z() - p2.Z());
    
    return xdiff * xdiff + ydiff * ydiff + zdiff * zdiff;
}



static RNScalar SquaredScaledDistance(int ix, int iy, int iz, int ii, int ij, int ik)
{
    return SquaredScaledDistance(R3Point(ix, iy, iz), R3Point(ii, ij, ik));
}



static R3Point IndexToPosition(int iv)
{
    int ix, iy, iz;
    segmentation_grid->IndexToIndices(iv, ix, iy, iz);
    return R3Point(ix, iy, iz);
}



static int PositionToIndex(R3Point p)
{
    int ix = (int)(p.X() + 0.5);
    int iy = (int)(p.Y() + 0.5);
    int iz = (int)(p.Z() + 0.5);
    
    int iv;
    segmentation_grid->IndicesToIndex(ix, iy, iz, iv);
    return iv;
}



////////////////////////////////////////////////////////////////////////
// Input/output functions
////////////////////////////////////////////////////////////////////////

struct SWCEntry {
    // constructors
    SWCEntry(int sample_number,
        int structure_identifier,
        RNScalar x_position,
        RNScalar y_position,
        RNScalar z_position,
        RNScalar radius,
        int parent_sample)
        : sample_number(sample_number)
        , structure_identifier(structure_identifier)
        , x_position(x_position)
        , y_position(y_position)
        , z_position(z_position)
        , radius(radius)
        , parent_sample(parent_sample)
    {
    }

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



static int ReadSWCFile(int segment)
{
    char input_filename[4096];
    sprintf(input_filename, "skeletons/%s/tree_%d.swc", prefix, segment);

    std::ifstream fd(input_filename);
    if(!fd.is_open()) {
        fprintf(stderr, "Failed to read %s\n", input_filename);
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

        skeletons[segment].push_back(
            SWCEntry(sample_number, structure_identifier, x_position, y_position, z_position, radius, parent_sample));
    }

    fd.close();

    // determine which points are endpoints
    RNBoolean* endpoints = new RNBoolean[skeletons[segment].size()];
    for(unsigned int ie = 0; ie < skeletons[segment].size(); ++ie) {
        endpoints[ie] = TRUE;
    }

    // not an endpoint if no children claim it as a parent
    for(unsigned int ie = 0; ie < skeletons[segment].size(); ++ie) {
        int parent_sample = skeletons[segment][ie].parent_sample;
        if(parent_sample == -1) continue;
        endpoints[parent_sample - 1] = FALSE;
    }

    for(unsigned int ie = 0; ie < skeletons[segment].size(); ++ie) {
        // if it has no parent it is an endpoint
        int parent_sample = skeletons[segment][ie].parent_sample;
        if(parent_sample == -1) endpoints[ie] = TRUE;
    }

    // add in all of the endpoints
    for(unsigned int ie = 0; ie < skeletons[segment].size(); ++ie) {
        if(endpoints[ie]) {
            skeleton_endpoints[segment].push_back(skeletons[segment][ie].P());
        }
    }

    // free memory
    delete[] endpoints;

    // return success
    return 1;
}



static int ReadData(void)
{
    // start statistics
    RNTime start_time;
    start_time.Read();

    // create default filenames
    if(!gold_filename) {
        char* filename = new char[4096];
        sprintf(filename, "gold/%s_gold.h5", prefix);
        gold_filename = filename;
        gold_dataset = "stack";
    }
    if(!segmentation_filename) {
        char* filename = new char[4096];
        sprintf(filename, "segmentations/%s_machine_labels_seg_28000.h5", prefix);
        segmentation_filename = filename;
        segmentation_dataset = "main";
    }

    R3Grid** golds = RNReadH5File(gold_filename, gold_dataset);
    if(!golds) { fprintf(stderr, "Failed to read %s from %s\n", gold_dataset, gold_filename); return 0; }
    gold_grid = golds[0];
    delete[] golds;

    R3Grid** segmentations = RNReadH5File(segmentation_filename, segmentation_dataset);
    if(!segmentations) { fprintf(stderr, "Failed to read %s from %s\n", segmentation_dataset, segmentation_filename); return 0; }
    segmentation_grid = segmentations[0];
    delete[] segmentations;

    // make sure the resolutions are the same
    rn_assertion(gold_grid->XResolution() == segmentation_grid->XResolution());
    rn_assertion(gold_grid->YResolution() == segmentation_grid->YResolution());
    rn_assertion(gold_grid->ZResolution() == segmentation_grid->ZResolution());

    resolution[RN_X] = gold_grid->XResolution();
    resolution[RN_Y] = gold_grid->YResolution();
    resolution[RN_Z] = gold_grid->ZResolution();

    // print statistics
    if(print_verbose) {
        printf("Read voxel grids...\n");
        printf("  Time = %.2f seconds\n", start_time.Elapsed());
        printf("  Resolution = (%d %d %d)\n", gold_grid->XResolution(), gold_grid->YResolution(), gold_grid->ZResolution());
        printf("  Scaling = (%0.2lf %0.2lf %0.2lf)\n", scaling[RN_X], scaling[RN_Y], scaling[RN_Z]);
    }

    // return success
    return 1;
}



////////////////////////////////////////////////////////////////////////
// Preprocessing functions
////////////////////////////////////////////////////////////////////////

static int MapSegmentationToGold(void)
{
    // get the number of segmentations and gold lables
    int max_segmentation = (int)(segmentation_grid->Maximum() + 0.5) + 1;
    int max_gold = (int)(gold_grid->Maximum() + 0.5) + 1;

    // create counters for segmentation and gold overlap
    int** segmentation_gold_overlap = new int*[max_segmentation];
    if (!segmentation_gold_overlap) return 0;
    for(int is = 0; is < max_segmentation; ++is) {
        segmentation_gold_overlap[is] = new int[max_gold];
        if (!segmentation_gold_overlap[is]) return 0;
        for(int ig = 0; ig < max_gold; ++ig) {
            segmentation_gold_overlap[is][ig] = 0;
        }
    }

    // iterate through all voxels
    for(int iv = 0; iv < segmentation_grid->NEntries(); ++iv) {
        int segmentation_id = (int)(segmentation_grid->GridValue(iv) + 0.5);
        int gold_id = (int)(gold_grid->GridValue(iv) + 0.5);

        segmentation_gold_overlap[segmentation_id][gold_id]++;
    }

    // create a mapping from segmentation ids to gold
    segmentation_to_gold = new int[max_segmentation];
    for(int is = 0; is < max_segmentation; ++is) {
        int max_gold_id = -1;
        int max_gold_value = -1;
        // skip the extracellular label
        for(int ig = 1; ig < max_gold; ++ig) {
            if(segmentation_gold_overlap[is][ig] > max_gold_value) {
                max_gold_value = segmentation_gold_overlap[is][ig];
                max_gold_id = ig;
            }
        }
        segmentation_to_gold[is] = max_gold_id;
    }

    // free memory
    for(int is = 0; is < max_segmentation; ++is) 
        delete[] segmentation_gold_overlap[is];
    delete[] segmentation_gold_overlap;
    
    // return success
    return 1;
}



static int ReadSWCFiles(void)
{
    // get the number of segmentations
    int max_segmentation = (int)(segmentation_grid->Maximum() + 0.5) + 1;
    
    // allocate memory for all segmentation vectors
    segmentations = new std::vector<int>[max_segmentation];
    if(!segmentations) return 0;
    skeletons = new std::vector<SWCEntry>[max_segmentation];
    if(!skeletons) return 0;
    skeleton_endpoints = new std::vector<R3Point>[max_segmentation];
    if(!skeleton_endpoints) return 0;

    for(int is = 0; is < max_segmentation; ++is) {
        segmentations[is] = std::vector<int>();
        skeletons[is] = std::vector<SWCEntry>();
        skeleton_endpoints[is] = std::vector<R3Point>();
    }

    // create boolean array of id existence
    RNBoolean *id_exists = new RNBoolean[max_segmentation];
    if (!id_exists) return 0;
    for(int is = 0; is < max_segmentation; ++is) 
        id_exists[is] = FALSE;

    for(int iv = 0; iv < segmentation_grid->NEntries(); ++iv) {
        int segmentation_id = (int)(segmentation_grid->GridValue(iv) + 0.5);
        id_exists[segmentation_id] = TRUE;
    }
    int nunique_segmentations = 0;
    for (int is = 0; is < max_segmentation; ++is) {
        if (id_exists[is]) nunique_segmentations++;
    }

    // read in all SWC files
    for(int is = 0; is < max_segmentation; ++is) {
        if(!id_exists[is]) continue;
        ReadSWCFile(is);
    }

    // print out the number of unique gold ids
    int max_gold = (int)(gold_grid->Maximum() + 0.5) + 1;
    RNBoolean* gold_id_exists = new RNBoolean[max_gold];
    for(int ig = 0; ig < max_gold; ++ig)
        gold_id_exists[ig] = FALSE;
    for(int iv = 0; iv < gold_grid->NEntries(); ++iv) {
        int gold_id = (int)(gold_grid->GridValue(iv) + 0.5);
        gold_id_exists[gold_id] = TRUE;
    }
    int nunique_gold = 0;
    for(int ig = 0; ig < max_gold; ++ig) {
        if(gold_id_exists[ig]) nunique_gold++;
    }

    // free memory
    delete[] id_exists;
    delete[] gold_id_exists;
    
    // no longer need gold grid
    delete gold_grid;

    if(print_verbose) {
        printf("  NSegments = %d\n", nunique_segmentations);
        printf("  NLabels = %d\n", nunique_gold);
    }

    // return success
    return 1;
}



static int PopulateSegmentationVectors(void)
{
    int max_segmentation = (int)(segmentation_grid->Maximum() + 0.5) + 1;

    // create neighbor array
    segmentation_neighbors = new RNBoolean *[max_segmentation];
    if (!segmentation_neighbors) return 0;
    for(int is1 = 0; is1 < max_segmentation; ++is1) {
        segmentation_neighbors[is1] = new RNBoolean[max_segmentation];
        if (!segmentation_neighbors[is1]) return 0;
        for(int is2 = 0; is2 < max_segmentation; ++is2) {
            segmentation_neighbors[is1][is2] = FALSE;
        }
    }

    // go through all voxels
    for(int ix = 0; ix < segmentation_grid->XResolution(); ++ix) {
        for(int iy = 0; iy < segmentation_grid->YResolution(); ++iy) {
            for(int iz = 0; iz < segmentation_grid->ZResolution(); ++iz) {
                // get the segmentation and grid value
                int segment_id = (int)(segmentation_grid->GridValue(ix, iy, iz) + 0.5);

                RNBoolean boundary = FALSE;
                if (ix != 0) {
                    int neighbor_id = (int)(segmentation_grid->GridValue(ix - 1, iy, iz) + 0.5);
                    segmentation_neighbors[segment_id][neighbor_id] = TRUE;
                    segmentation_neighbors[neighbor_id][segment_id] = TRUE;
                    if (neighbor_id != segment_id) boundary = TRUE;
                }
                if (iy != 0) {
                    int neighbor_id = (int)(segmentation_grid->GridValue(ix, iy - 1, iz) + 0.5);
                    segmentation_neighbors[segment_id][neighbor_id] = TRUE;
                    segmentation_neighbors[neighbor_id][segment_id] = TRUE;
                    if (neighbor_id != segment_id) boundary = TRUE;
                }
                if (iz != 0) {
                    int neighbor_id = (int)(segmentation_grid->GridValue(ix, iy, iz - 1) + 0.5);
                    segmentation_neighbors[segment_id][neighbor_id] = TRUE;
                    segmentation_neighbors[neighbor_id][segment_id] = TRUE;
                    if (neighbor_id != segment_id) boundary = TRUE;
                }
                if (ix != resolution[RN_X] - 1) {
                    int neighbor_id = (int)(segmentation_grid->GridValue(ix + 1, iy, iz) + 0.5);
                    if (neighbor_id != segment_id) boundary = TRUE;
                }
                if (iy != resolution[RN_Y] - 1) {
                    int neighbor_id = (int)(segmentation_grid->GridValue(ix, iy + 1, iz) + 0.5);
                    if (neighbor_id != segment_id) boundary = TRUE;
                }
                if (iz != resolution[RN_Z] - 1) {
                    int neighbor_id = (int)(segmentation_grid->GridValue(ix, iy, iz + 1) + 0.5);
                    if (neighbor_id != segment_id) boundary = TRUE;
                }

                // get the segment gold id
                if (!boundary) continue;
                int iv;
                segmentation_grid->IndicesToIndex(ix, iy, iz, iv);
                segmentations[segment_id].push_back(iv);
            }
        }
    }

    // return success
    return 1;
}



////////////////////////////////////////////////////////////////////////
// Calculate boundary examples
//////////////////////////////////////////////////////////////////////

static int CalculateExamples(void)
{
    // avoid pair duplication
    int max_segmentation = (int)(segmentation_grid->Maximum() + 0.5) + 1;
    RNBoolean** considered_merges = new RNBoolean*[max_segmentation];
    if (!considered_merges) return 0;
    RNBoolean** possible_misses = new RNBoolean*[max_segmentation];
    if (!possible_misses) return 0;
    for(int is1 = 0; is1 < max_segmentation; ++is1) {
        considered_merges[is1] = new RNBoolean[max_segmentation];
        if (!considered_merges[is1]) continue;
        possible_misses[is1] = new RNBoolean[max_segmentation];
        if (!possible_misses[is1]) continue;
        for(int is2 = 0; is2 < max_segmentation; ++is2) {
            considered_merges[is1][is2] = FALSE;
            possible_misses[is1][is2] = FALSE;
        }
    }

    // iterate over all segments
    RNTime start_time;
    start_time.Read();
    if(print_verbose) printf("Generating boundary examples...\n  ");
    for(int is = 0; is < max_segmentation; ++is) {
        if(print_verbose) RNProgressBar(is, max_segmentation);

        // only consider segments that exist
        if (!skeleton_endpoints[is].size()) continue;

        // initialize the priority queue
        SegmentPair tmp;
        RNMinBinaryHeap<SegmentPair *> heap = RNMinBinaryHeap<SegmentPair *>(&tmp, &(tmp.distance), max_segmentation);

        // allocate temporary memory
        SegmentPair *segment_data = new SegmentPair[max_segmentation];
        for(int in = 0; in < max_segmentation; ++in) {
            segment_data[in].distance = FLT_MAX;
            segment_data[in].segment_one = is;
            segment_data[in].segment_two = in;
            segment_data[in].merge = (segmentation_to_gold[is] == segmentation_to_gold[in]);
            heap.Insert(in, &(segment_data[in]));
        }

        // iterate over all endpoints
        for(unsigned int ie = 0; ie < skeleton_endpoints[is].size(); ++ie) {
            R3Point position = skeleton_endpoints[is][ie];

            int ix = (int)(position.X() + 0.5);
            int iy = (int)(position.Y() + 0.5);
            int iz = (int)(position.Z() + 0.5);

            int box_radius[3] = { (int)(maximum_distance[RN_X] + 0.5), (int)(maximum_distance[RN_Y] + 0.5), (int)(maximum_distance[RN_Z] + 0.5) };

            // consider elements in the bounding box
            for(int ik = iz - box_radius[RN_Z]; ik <= iz + box_radius[RN_Z]; ++ik) {
                if(ik < 0 || ik > resolution[RN_Z] - 1) continue;
                for(int ij = iy - box_radius[RN_Y]; ij <= iy + box_radius[RN_Y]; ++ij) {
                    if(ij < 0 || ij > resolution[RN_Y] - 1) continue;
                    for(int ii = ix - box_radius[RN_X]; ii <= ix + box_radius[RN_X]; ++ii) {
                        if(ii < 0 || ii > resolution[RN_X] - 1) continue;

                        // get this segment id
                        int neighbor_id = (int)(segmentation_grid->GridValue(ii, ij, ik) + 0.5);

                        // do not consider the same id
                        if(neighbor_id == is) continue;

                        // get the data for this neighbor
                        SegmentPair* neighbor_data = &(segment_data[neighbor_id]);

                        // should we update distance?
                        RNScalar distance = SquaredScaledDistance(ix, iy, iz, ii, ij, ik);
                        if(distance < neighbor_data->distance) {
                            neighbor_data->distance = distance;
                            // the new center point is this closest location
                            neighbor_data->center_point = R3Point(ii, ij, ik);
                            heap.DecreaseKey(neighbor_id, neighbor_data);
                        }
                    }
                }
            }
        }

        // add the 10 closest segments
        for(int ii = 0; ii < nclosest; ++ii) {
            if (heap.IsEmpty()) break;
            SegmentPair* pair = heap.DeleteMin();

            // make sure that the pair actually belonged to the bounding box (i.e. distance != FLT_MAX)
            if(pair->distance > 3 * normalized_distance * normalized_distance) break;

            // do not add if already considered
            if (considered_merges[pair->segment_one][pair->segment_two]) continue;

            // add this pair of segments
            boundary_examples.push_back(new SegmentPair(*pair));

            // mark as considered
            considered_merges[pair->segment_one][pair->segment_two] = TRUE;
            considered_merges[pair->segment_two][pair->segment_one] = TRUE;

            // no longer a possible miss
            possible_misses[pair->segment_one][pair->segment_two] = FALSE;
            possible_misses[pair->segment_two][pair->segment_one] = FALSE;
        }

        while(!heap.IsEmpty()) {
            SegmentPair* pair = heap.DeleteMin();

            // if already considered, no possible miss
            if(considered_merges[pair->segment_one][pair->segment_two]) continue;
            
            // make sure that the pair actually belonged to the bounding box (i.e. distance != FLT_MAX)
            if(pair->distance > 3 * normalized_distance * normalized_distance) continue;

            // update as a possible miss
            if(pair->merge) {
                possible_misses[pair->segment_one][pair->segment_two] = TRUE;
                possible_misses[pair->segment_two][pair->segment_one] = TRUE;
            }
        }

        // free memory
        delete[] segment_data;
    }
    if(print_verbose) printf("\nDone in %0.2f seconds.\n", start_time.Elapsed());


    // determine the number of non-considered false splits
    int nmissed_splits = 0;
    for(int is1 = 0; is1 < max_segmentation; ++is1) {
        for(int is2 = is1 + 1; is2 < max_segmentation; ++is2) {
            if(possible_misses[is1][is2]) nmissed_splits++;
        }
    }

    // determine the number of boundaries which should merge
    int nfalse_splits = 0;
    for (unsigned int ib = 0; ib < boundary_examples.size(); ++ib) {
        if (boundary_examples[ib]->merge) nfalse_splits++;
    }

    // free memory
    for(int is = 0; is < max_segmentation; ++is) {
        delete[] considered_merges[is];
        delete[] possible_misses[is];
    }
    delete[] considered_merges;
    delete[] possible_misses;

    // print out the number of entries and the number of false splits
    printf("Maximum Distance: %d\n", normalized_distance);
    printf("Closest: %d\n", nclosest);
    printf("  False Splits: %d\n", nfalse_splits);
    printf("  Total Pairs: %lu\n", boundary_examples.size());
    printf("  Missed splits: %d\n", nmissed_splits);

    return 1;
}



static int UpdateBoundingBox(void)
{
    // update the center points
    RNTime start_time;
    start_time.Read();
    if(print_verbose) printf("Updating bounding box centers...\n  ");
    for(unsigned int ie = 0; ie < boundary_examples.size(); ++ie) {
        if(print_verbose) RNProgressBar(ie, boundary_examples.size());
        SegmentPair *merge_candidate = boundary_examples[ie];

        // get the segment ids
        int segment_one = merge_candidate->segment_one;
        int segment_two = merge_candidate->segment_two;

        if(segmentation_neighbors[segment_one][segment_two]) {
            // get a list of vertices between these two segments
            std::vector<int>* points = (segmentations[segment_one].size() < segmentations[segment_two].size()) ? &(segmentations[segment_one]) : &(segmentations[segment_two]);
            std::vector<R3Point> boundary = std::vector<R3Point>();

            // go through all of the points
            for(unsigned int ip = 0; ip < points->size(); ++ip) {
                R3Point point = IndexToPosition((*points)[ip]);
                int ix = (int)(point.X() + 0.5);
                int iy = (int)(point.Y() + 0.5);
                int iz = (int)(point.Z() + 0.5);

                int segment_id = segmentation_grid->GridValue(ix, iy, iz);
                int other_segment_id = (segment_id == segment_one) ? segment_two : segment_one;

                // see if this has segment_two has a neighbor
                if(ix > 0) {
                    int neighbor_id = segmentation_grid->GridValue(ix - 1, iy, iz);
                    if(neighbor_id == other_segment_id) boundary.push_back(R3Point(ix - 0.5, iy, iz));
                }
                if(ix < resolution[RN_X] - 1) {
                    int neighbor_id = segmentation_grid->GridValue(ix + 1, iy, iz);
                    if(neighbor_id == other_segment_id) boundary.push_back(R3Point(ix + 0.5, iy, iz));
                }
                if(iy > 0) {
                    int neighbor_id = segmentation_grid->GridValue(ix, iy - 1, iz);
                    if(neighbor_id == other_segment_id) boundary.push_back(R3Point(ix, iy - 0.5, iz));
                }
                if(iy < resolution[RN_Y] - 1) {
                    int neighbor_id = segmentation_grid->GridValue(ix, iy + 1, iz);
                    if(neighbor_id == other_segment_id) boundary.push_back(R3Point(ix, iy + 0.5, iz));
                }
                if(iz > 0) {
                    int neighbor_id = segmentation_grid->GridValue(ix, iy, iz - 1);
                    if(neighbor_id == other_segment_id) boundary.push_back(R3Point(ix, iy, iz - 0.5));
                }
                if(iz < resolution[RN_Z] - 1) {
                    int neighbor_id = segmentation_grid->GridValue(ix, iy, iz + 1);
                    if(neighbor_id == other_segment_id) boundary.push_back(R3Point(ix, iy, iz + 0.5));
                }
            }

            // iterate through all points on the boundary
            R3Box bounding_box = R3null_box;
            for(unsigned int ip = 0; ip < boundary.size(); ++ip) {
                bounding_box.Union(boundary[ip]);
            }

            // set the center for this example
            R3Point center = bounding_box.Centroid();
            boundary_examples[ie]->center_point = center;

            RNScalar xradius = 1.5 * bounding_box.AxisRadius(RN_X);
            RNScalar yradius = 1.5 * bounding_box.AxisRadius(RN_Y);
            RNScalar zradius = 1.5 * bounding_box.AxisRadius(RN_Z);

            if(xradius < maximum_distance[RN_X]) xradius = maximum_distance[RN_X];
            if(yradius < maximum_distance[RN_Y]) yradius = maximum_distance[RN_Y];
            if(zradius < maximum_distance[RN_Z]) zradius = maximum_distance[RN_Z];

            boundary_examples[ie]->bbox = R3Box(center.X() - xradius, center.Y() - yradius, center.Z() - zradius, center.X() + xradius, center.Y() + yradius, center.Z() + zradius);
        } else {
            R3Point first_point = boundary_examples[ie]->center_point;

            // get the segment to which this index belongs
            int grid_value = (int)(segmentation_grid->GridValue(first_point) + 0.5);

            // get the grid value, if it equals segment one, iterate through segment two boundary points
            std::vector<int>* points = (grid_value == segment_one) ? &(segmentations[segment_two]) : &(segmentations[segment_one]);

            R3Point middle_point = R3null_point;
            RNScalar closest_distance = FLT_MAX;
            for(unsigned int ip = 0; ip < points->size(); ++ip) {
                R3Point second_point = IndexToPosition((*points)[ip]);
                RNScalar distance = SquaredScaledDistance(first_point, second_point);

                // update if the closest poitn
                if(distance < closest_distance) {
                    closest_distance = distance;
                    middle_point = (first_point + second_point) / 2;
                }
            }

            // create the bounding box
            RNScalar xradius = maximum_distance[RN_X];
            RNScalar yradius = maximum_distance[RN_Y];
            RNScalar zradius = maximum_distance[RN_Z];

            boundary_examples[ie]->center_point = middle_point;
            boundary_examples[ie]->bbox = R3Box(middle_point.X() - xradius, middle_point.Y() - yradius, middle_point.Z() - zradius, middle_point.X() + xradius, middle_point.Y() + yradius, middle_point.Z() + zradius);
        }
    }
    if(print_verbose) printf("\ndone in %0.2f seconds.\n", start_time.Elapsed());
    
    // return success
    return 1;
}



static int UpdatePointFeatures(void)
{
    // create a segment to pairs array
    int max_segmentation = (int)(segmentation_grid->Maximum() + 0.5) + 1;
    std::vector<int> *segments_to_pairs = new std::vector<int>[max_segmentation];
    if (!segments_to_pairs) return 0;
    for (int is = 0; is < max_segmentation; ++is)
        segments_to_pairs[is] = std::vector<int>();
    
    // go through every boundary example
    for (unsigned int ib = 0; ib < boundary_examples.size(); ++ib) {
        segments_to_pairs[boundary_examples[ib]->segment_one].push_back(ib);
        segments_to_pairs[boundary_examples[ib]->segment_two].push_back(ib);
    }
    
    // go through every coordinate
    RNTime start_time;
    start_time.Read();
    if (print_verbose) printf("Updated per bounding box segmentation points...\n  "); 
    for (int iz = 0; iz < resolution[RN_Z]; ++iz) {
        RNProgressBar(iz, resolution[RN_Z]);
        for (int iy = 0; iy < resolution[RN_Y]; ++iy) {
            for (int ix = 0; ix < resolution[RN_X]; ++ix) {
                int grid_value = (int)(segmentation_grid->GridValue(ix, iy, iz) + 0.5);
                
                R3Point position = R3Point(ix, iy, iz);
                for (unsigned int ip = 0; ip < segments_to_pairs[grid_value].size(); ++ip) {
                    SegmentPair *pair = boundary_examples[segments_to_pairs[grid_value][ip]];
                    
                    if (R3Intersects(pair->bbox, position)) {
                        if (grid_value == pair->segment_one) pair->segment_one_points.push_back(PositionToIndex(position));
                        else pair->segment_two_points.push_back(PositionToIndex(position));
                    }
                }
            }
        }
    }
    if (print_verbose) printf("\ndone in %0.2f seconds.\n", start_time.Elapsed());
    
    delete[] segments_to_pairs;
    
    // return success
    return 1;
}



////////////////////////////////////////////////////////////////////////
// Create features
////////////////////////////////////////////////////////////////////////

void AddBoxFeatures(void)
{
    RNTime start_time;
    start_time.Read();
    if (print_verbose) printf("Creating box features...\n  ");
    for (unsigned int ie = 0; ie < boundary_examples.size(); ++ie) {
        if (print_verbose) RNProgressBar(ie, boundary_examples.size());
        SegmentPair *pair = boundary_examples[ie];
        
        // simplify instance variable access
        int segment_one = pair->segment_one;
        int segment_two = pair->segment_two;
        R3Box bbox = pair->bbox;
        R3Point center_point = pair->center_point;

        
        //////////////////////////
        //// COUNTING METRICS ////
        //////////////////////////
     
        // add the bounding box volume
        attributes[ie].push_back((pair->segment_one_points.size() + pair->segment_two_points.size()) / bbox.Volume());
        if (pair->segment_one_points.size() > pair->segment_two_points.size())
            attributes[ie].push_back(pair->segment_one_points.size() / (RNScalar)pair->segment_two_points.size());
        else
            attributes[ie].push_back(pair->segment_two_points.size() / (RNScalar)pair->segment_one_points.size());
        if (segmentations[segment_one].size() > segmentations[segment_two].size()) {
            attributes[ie].push_back(segmentations[segment_one].size());
            attributes[ie].push_back(segmentations[segment_two].size());
        }
        else {
            attributes[ie].push_back(segmentations[segment_two].size());
            attributes[ie].push_back(segmentations[segment_one].size());
        }
            
        
        ///////////////////////////
        //// CURVATURE METRICS ////
        ///////////////////////////
        
        // find the furthest point from the central points for each segment        
        RNScalar segment_one_furthest_distance = -1;
        RNScalar segment_two_furthest_distance = -1;
        for (unsigned int ip = 0; ip < pair->segment_one_points.size(); ++ip) {
            R3Point segment_point = IndexToPosition(pair->segment_one_points[ip]);
            
            RNScalar distance = SquaredScaledDistance(center_point, segment_point);
            if (distance > segment_one_furthest_distance) {
                segment_one_furthest_distance = distance;
            }
        }
        for (unsigned int ip = 0; ip < pair->segment_two_points.size(); ++ip) {
            R3Point segment_point = IndexToPosition(pair->segment_two_points[ip]);
            
            RNScalar distance = SquaredScaledDistance(center_point, segment_point);
            if (distance > segment_two_furthest_distance) {
                segment_two_furthest_distance = distance;
            }
        }
        
        // create vectors of points with a small window of the maximum distance
        std::vector<R3Point> segment_one_points = std::vector<R3Point>();
        std::vector<R3Point> segment_two_points = std::vector<R3Point>();
        RNScalar window_size = 5.0;
        for (unsigned int ip = 0; ip < pair->segment_one_points.size(); ++ip) {
            R3Point segment_point = IndexToPosition(pair->segment_one_points[ip]);
            
            RNScalar distance = SquaredScaledDistance(center_point, segment_point);
            if (distance > segment_one_furthest_distance - window_size) {
                segment_one_points.push_back(segment_point);
            }
        }
        for (unsigned int ip = 0; ip < pair->segment_two_points.size(); ++ip) {
            R3Point segment_point = IndexToPosition(pair->segment_two_points[ip]);
            
            RNScalar distance = SquaredScaledDistance(center_point, segment_point);
            if (distance > segment_two_furthest_distance - window_size) {
                segment_two_points.push_back(segment_point);
            }
        }

        // run ntrial iterations
        int ntrials = 100;
        RNScalar average_angle = 0.0;
        std::vector<RNScalar> angles = std::vector<RNScalar>();
        for (int it = 0; it < ntrials; ++it) {
            int index_one = (int)(RNRandomScalar() * segment_one_points.size() + 0.5);
            int index_two = (int)(RNRandomScalar() * segment_two_points.size() + 0.5);
            
            R3Vector segment_one_vector = segment_one_points[index_one] - center_point;
            R3Vector segment_two_vector = segment_two_points[index_two] - center_point;
            
            segment_one_vector.Normalize();
            segment_two_vector.Normalize();
            
            // get the cosine (angle)
            RNScalar cosine_angle = segment_one_vector.Dot(segment_two_vector);
            
            average_angle += cosine_angle / ntrials;
            angles.push_back(cosine_angle);
        }
        RNScalar stddev_angle = 0.0;
        for (int it = 0; it < ntrials; ++it) {
            stddev_angle += (angles[it] - average_angle) * (angles[it] - average_angle);
        }
        stddev_angle /= (ntrials - 1);
        stddev_angle = sqrt(stddev_angle);
        
        attributes[ie].push_back(average_angle);
        attributes[ie].push_back(stddev_angle);
        
        
        //////////////////////////
        //// BOUNDARY METRICS ////
        //////////////////////////
        
        // calculate the number of voxels on the boundary between two segments
        std::vector<int> *points = (segmentations[segment_one].size() < segmentations[segment_two].size()) ? &(segmentations[segment_one]) : &(segmentations[segment_two]);
        
        // iterate through all points
        int nneighbors_locations = 0;
        for (unsigned int ip = 0; ip < points->size(); ++ip) {
            R3Point position = IndexToPosition((*points)[ip]);
            
            int ix = (int)(position.X() + 0.5);
            int iy = (int)(position.Y() + 0.5);
            int iz = (int)(position.Z() + 0.5);
            
            int segment_id = (int)(segmentation_grid->GridValue(ix, iy, iz) + 0.5);
            
            // go through all neighbors
            if (ix > 0) {
                int neighbor_id = (int)(segmentation_grid->GridValue(ix - 1, iy, iz) + 0.5);
                if (segment_id == neighbor_id) nneighbors_locations++;
            }
            if (ix < resolution[RN_X] - 1) {
                int neighbor_id = (int)(segmentation_grid->GridValue(ix + 1, iy, iz) + 0.5);
                if (segment_id == neighbor_id) nneighbors_locations++;
            }
            if (iy > 0) { 
                int neighbor_id = (int)(segmentation_grid->GridValue(ix, iy - 1, iz) + 0.5);
                if (segment_id == neighbor_id) nneighbors_locations++;
            }
            if (iy < resolution[RN_Z] - 1) { 
                int neighbor_id = (int)(segmentation_grid->GridValue(ix, iy + 1, iz) + 0.5);
                if (segment_id == neighbor_id) nneighbors_locations++;
            }
            
            if (iz > 0) {
                int neighbor_id = (int)(segmentation_grid->GridValue(ix, iy, iz - 1) + 0.5);
                if (segment_id == neighbor_id) nneighbors_locations++;
            }
            if (iz < resolution[RN_Z] - 1) {
                int neighbor_id = (int)(segmentation_grid->GridValue(ix, iy, iz + 1) + 0.5);
                if (segment_id == neighbor_id) nneighbors_locations++;
            }
        }
        
        attributes[ie].push_back(nneighbors_locations);
        if (segmentations[segment_one].size() < segmentations[segment_two].size()) {
            attributes[ie].push_back(segmentations[segment_one].size() / (RNScalar)nneighbors_locations);
            attributes[ie].push_back(segmentations[segment_two].size() / (RNScalar)nneighbors_locations);
        }
        else {
            attributes[ie].push_back(segmentations[segment_two].size() / (RNScalar)nneighbors_locations);
            attributes[ie].push_back(segmentations[segment_one].size() / (RNScalar)nneighbors_locations);            
        }
        
        
    }
    if (print_verbose) printf("\ndone in %0.2f seconds.\n", start_time.Elapsed());
}




int SaveFeatures(void)
{
    char attributes_filename[4096];
    sprintf(attributes_filename, "features/%s_%d_%d.features", prefix, normalized_distance, nclosest);
    
    // open file
    FILE *fp = fopen(attributes_filename, "wb");
    if (!fp) return 0;
    
    unsigned int nboundary_examples = boundary_examples.size();
    fwrite(&(nboundary_examples), sizeof(unsigned int), 1, fp);
    unsigned int nattributes = attributes[0].size();
    fwrite(&(nattributes), sizeof(unsigned int), 1, fp);
    
    // add all attributes to the file
    for (unsigned int ip = 0; ip < nboundary_examples; ++ip) {
        for (unsigned int ia = 0; ia < nattributes; ++ia) {
            fwrite(&(attributes[ip][ia]), sizeof(RNScalar), 1, fp);
        }
    }
    
    // close file
    fclose(fp);
    
    // return success
    return 1;
}



int CreateLabelFile(void)
{
    // create gold standard output
    char label_filename[4096];
    sprintf(label_filename, "features/%s_%d_%d.labels", prefix, normalized_distance, nclosest);

    // open file
    FILE* fp = fopen(label_filename, "wb");
    if(!fp) return 0;

    // create the gold standard lables
    unsigned int nboundary_examples = boundary_examples.size();
    fwrite(&(nboundary_examples), sizeof(unsigned int), 1, fp);

    for(unsigned int ie = 0; ie < boundary_examples.size(); ++ie) {
        SegmentPair *pair = boundary_examples[ie];

        fwrite(&(pair->segment_one), sizeof(int), 1, fp);
        fwrite(&(pair->segment_two), sizeof(int), 1, fp);
        fwrite(&(pair->merge), sizeof(int), 1, fp);
    }

    // close file
    fclose(fp);

    // return success
    return 1;
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
            else if(!strcmp(*argv, "-gold")) {
                argv++; argc--;
                gold_filename = *argv;
                argv++; argc--;
                gold_dataset = *argv;
            } else if(!strcmp(*argv, "-segmentation")) {
                argv++; argc--;
                segmentation_filename = *argv;
                argv++; argc--;
                segmentation_dataset = *argv;
            } else if(!strcmp(*argv, "-scaling")) {
                argv++; argc--;
                scaling[RN_X] = atof(*argv);
                argv++; argc--;
                scaling[RN_Y] = atof(*argv);
                argv++; argc--;
                scaling[RN_Z] = atof(*argv);
            } else if(!strcmp(*argv, "-closest")) { argv++; argc--; nclosest = atoi(*argv);
            } else if(!strcmp(*argv, "-max_distance")) { argv++; argc--; normalized_distance = atoi(*argv);
            } else { fprintf(stderr, "Invalid program argument: %s\n", *argv); return 0; }
        } else {
            if(!prefix) prefix = *argv;
            else { fprintf(stderr, "Invalid program argument: %s\n", *argv); return 0;
            }
        }
        argv++; argc--;
    }

    // error if there is no input name
    if(!prefix && (!gold_filename && !segmentation_filename)) { fprintf(stderr, "Need to supply a neuron data file\n"); return 0; }

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

    // set projection scale variables
    for(int dim = 0; dim <= 2; ++dim) {
        maximum_distance[dim] = normalized_distance / scaling[dim];
    }
    
    // create mapping from segmentation to gold
    if(!MapSegmentationToGold()) { fprintf(stderr, "Failed to allocate memory for segmentation-gold mapping.\n"); exit(-1); }
    
    // create mapping from index to ids
    if(!ReadSWCFiles()) { fprintf(stderr, "Failed to allocate memory for skeletons.\n"); exit(-1); }
    
    // create segmentation vectors
    if(!PopulateSegmentationVectors()) { fprintf(stderr, "Failed to allocate memory for segmentation vectors.\n"); exit(-1); }
    
    if(!CalculateExamples()) { fprintf(stderr, "Failed to calculate merge and split candidates.\n"); exit(-1); }
    
    if (!UpdateBoundingBox()) { fprintf(stderr, "Failed to update the bounding box.\n"); exit(-1); }
    
    if (!UpdatePointFeatures()) { fprintf(stderr, "Failed to update points features.\n"); exit(-1); }
    
    // make sure that no boundary pairs have the same indices
    for (unsigned int ip1 = 0; ip1 < boundary_examples.size(); ++ip1) {
        int pair_one_index_one = boundary_examples[ip1]->segment_one;
        int pair_one_index_two = boundary_examples[ip1]->segment_two;
        for (unsigned int ip2 = ip1 + 1; ip2 < boundary_examples.size(); ++ip2) {
            int pair_two_index_one = boundary_examples[ip2]->segment_one;
            int pair_two_index_two = boundary_examples[ip2]->segment_two;
            
            if (pair_one_index_one == pair_two_index_one) rn_assertion(pair_one_index_two != pair_two_index_two);
            if (pair_one_index_one == pair_two_index_two) rn_assertion(pair_one_index_two != pair_two_index_one);
        }
    }

    // allocate memory for attributes vector
    attributes = new std::vector<RNScalar>[boundary_examples.size()];
    for (unsigned int ip = 0; ip < boundary_examples.size(); ++ip) 
        attributes[ip] = std::vector<RNScalar>();
    
    //////////////////////////
    //// FEATURE CREATION ////
    //////////////////////////
    
    AddBoxFeatures();
    
    if (!SaveFeatures()) { fprintf(stderr, "Failed to write to feature attributes file\n"); exit(-1); }

    if(!CreateLabelFile()) { fprintf(stderr, "Failed to write to feature label file\n"); exit(-1); }

    // free memory
    delete[] segmentations;
    delete[] skeletons;
    delete[] skeleton_endpoints;

    // return success
    return 1;
}