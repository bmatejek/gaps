// Source file for the mesh viewer program

// include files
#include "RNDataStructures/RNDataStructures.h"
#include <vector>

// program variables

static const char* input_filename;
static const char* dataset_name;
static int print_verbose = 0;
static int print_debug = 0;
static R3Grid* machine_labels = NULL;
static std::vector<R3Point>* exterior_points = NULL;
static int resolution[3];
static int ntests = 100;

enum DIRECTION { UP, DOWN, NORTH, SOUTH, EAST, WEST, NDIRECTIONS };

enum RIDGE_STRENGTH { NO_RIDGE, WEAK_RIDGE, GOOD_RIDGE, STRONG_RIDGE };

void IndexToIndices(int index, int& ix, int& iy, int& iz) { machine_labels->IndexToIndices(index, ix, iy, iz); }

int IndicesToIndex(int ix, int iy, int iz)
{
    int index;
    machine_labels->IndicesToIndex(ix, iy, iz, index);
    return index;
}

int Neighbor(int iv, DIRECTION direction)
{
    int ix, iy, iz;
    IndexToIndices(iv, ix, iy, iz);

    if(direction == UP) {
        if(iz == resolution[RN_Z] - 1)
            return -1;
        else
            return IndicesToIndex(ix, iy, iz + 1);
    }
    if(direction == DOWN) {
        if(iz == 0)
            return -1;
        else
            return IndicesToIndex(ix, iy, iz - 1);
    }
    if(direction == NORTH) {
        if(ix == resolution[RN_X] - 1)
            return -1;
        else
            return IndicesToIndex(ix + 1, iy, iz);
    }
    if(direction == SOUTH) {
        if(ix == 0)
            return -1;
        else
            return IndicesToIndex(ix - 1, iy, iz);
    }
    if(direction == EAST) {
        if(iy == resolution[RN_Y] - 1)
            return -1;
        else
            return IndicesToIndex(ix, iy + 1, iz);
    }
    if(direction == WEST) {
        if(iy == 0)
            return -1;
        else
            return IndicesToIndex(ix, iy - 1, iz);
    }

    return -1;
}

////////////////////////////////////////////////////////////////////////
// Euclidean Distance Transform
////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////
// Program argument parsing
////////////////////////////////////////////////////////////////////////

static int ParseArgs(int argc, char** argv)
{
    // parse arguments
    argc--;
    argv++;
    while(argc > 0) {
        if((*argv)[0] == '-') {
            if(!strcmp(*argv, "-v"))
                print_verbose = 1;
            else if(!strcmp(*argv, "-debug"))
                print_debug = 1;
            else if(!strcmp(*argv, "-ntests")) {
                argv++;
                argc--;
                ntests = atoi(*argv);
            } else {
                fprintf(stderr, "Invalid program argument: %s", *argv);
                exit(1);
            }
        } else {
            if(!input_filename)
                input_filename = *argv;
            else if(!dataset_name)
                dataset_name = *argv;
            else {
                fprintf(stderr, "Invalid program argument: %s", *argv);
                exit(1);
            }
        }
        argv++;
        argc--;
    }

    // check filenames
    if(!input_filename || !dataset_name) {
        fprintf(stderr, "Usage: skeleton h5_file dataset_name [options]\n");
        return 0;
    }

    // return OK status
    return 1;
}

////////////////////////////////////////////////////////////////////////
// Main
////////////////////////////////////////////////////////////////////////

RNScalar f(int x, int i, int** g, int y)
{
    int g_i = g[i][y];

    return (x - i) * (x - i) + g_i * g_i;
}

int Sep(int i, int u, int** g, int y)
{
    int g_u = g[u][y];
    int g_i = g[i][y];

    return (u * u - i * i + g_u * g_u - g_i * g_i) / (2 * (u - i));
}

void TestDistanceTransform(R3Grid* distance_transform, R3Point*** boundary)
{
    int nentries = distance_transform->NEntries();
    for(int it = 0; it < ntests; ++it) {
        int iv = (int)(nentries * RNRandomScalar() + 0.5);

        RNScalar distance = distance_transform->GridValue(iv);

        int grid_value = (int)(machine_labels->GridValue(iv) + 0.5);

        int ix, iy, iz;
        IndexToIndices(iv, ix, iy, iz);
        R3Point position = R3Point(ix, iy, iz);

        R3Point nearest_point;

        // go through all points to see the correct distance
        RNScalar minimum_distance = FLT_MAX;
        for(unsigned int ib = 0; ib < exterior_points[grid_value].size(); ++ib) {
            R3Point boundary_point = exterior_points[grid_value][ib];

            RNScalar distance = R3Distance(position, boundary_point);
            if(distance < minimum_distance) {
                minimum_distance = distance;
                nearest_point = boundary_point;
            }
        }

        RNScalar epsilon = 10e-5;
        rn_assertion(fabs(distance - minimum_distance) < epsilon);

        RNScalar boundary_distance = R3Distance(position, boundary[ix][iy][iz]);
        rn_assertion(fabs(distance - boundary_distance) < epsilon);
    }
}

int sign(RNScalar arg1, RNScalar arg2)
{
    RNScalar diff = arg2 - arg1;
    if(diff == 0.0)
        return 0;
    else if(diff > 0.0)
        return 1;
    else
        return -1;
}

int main(int argc, char** argv)
{
    // measure the entire time
    RNTime total_time;
    total_time.Read();

    // parse program arguments
    if(!ParseArgs(argc, argv)) exit(-1);

    // read grid
    R3Grid** grids = RNReadH5File(input_filename, dataset_name);
    if(!grids) {
        fprintf(stderr, "Failed to read %s from %s\n", dataset_name, input_filename);
        exit(-1);
    }
    machine_labels = grids[0];
    delete[] grids;

    // update the resolutions
    for(int dim = 0; dim <= 2; ++dim) {
        resolution[dim] = machine_labels->Resolution(dim);
    }
    long long int nentries = resolution[RN_X] * resolution[RN_Y] * resolution[RN_Z];

    RNTime start_time;
    start_time.Read();

    int maximum_label = (int)(machine_labels->Maximum() + 0.5) + 1;
    exterior_points = new std::vector<R3Point>[maximum_label];
    for(int il = 0; il < maximum_label; ++il) {
        exterior_points[il] = std::vector<R3Point>();
    }

    // create boundary values
    long long int*** b = new long long int**[resolution[RN_X]];
    unsigned int*** zpass = new unsigned int**[resolution[RN_X]];
    unsigned int*** ypass = new unsigned int**[resolution[RN_X]];
    unsigned int*** xpass = new unsigned int**[resolution[RN_X]];
    for(int ix = 0; ix < resolution[RN_X]; ++ix) {
        b[ix] = new long long int*[resolution[RN_Y]];
        zpass[ix] = new unsigned int*[resolution[RN_Y]];
        ypass[ix] = new unsigned int*[resolution[RN_Y]];
        xpass[ix] = new unsigned int*[resolution[RN_Y]];
        for(int iy = 0; iy < resolution[RN_Y]; ++iy) {
            b[ix][iy] = new long long int[resolution[RN_Z]];
            zpass[ix][iy] = new unsigned int[resolution[RN_Z]];
            ypass[ix][iy] = new unsigned int[resolution[RN_Z]];
            xpass[ix][iy] = new unsigned int[resolution[RN_Z]];
        }
    }

    // use this value for infinity, no distance can exceed this value
    const RNScalar infinity = nentries * nentries * 3;

    // create boundary map
    for(int iv = 0; iv < nentries; ++iv) {
        // get the indices
        int ix, iy, iz;
        IndexToIndices(iv, ix, iy, iz);

        int grid_value = (int)(machine_labels->GridValue(iv) + 0.5);

        RNBoolean interior = TRUE;

        // go through all neighbors
        for(int in = 0; in < NDIRECTIONS; ++in) {
            int neighbor_index = Neighbor(iv, DIRECTION(in));
            if(neighbor_index == -1) continue;

            int neighbor_grid_value = (int)(machine_labels->GridValue(neighbor_index) + 0.5);

            if(grid_value != neighbor_grid_value) interior = FALSE;
        }

        if(interior)
            b[ix][iy][iz] = infinity;
        else {
            b[ix][iy][iz] = 0;
            exterior_points[grid_value].push_back(R3Point(ix, iy, iz));
        }
    }

    printf("Created boundary map in %0.2f seconds.\n", start_time.Elapsed());
    start_time.Read();

    // create the distance transform array
    long long int*** dt = new long long int**[resolution[RN_X]];
    for(int ix = 0; ix < resolution[RN_X]; ++ix) {
        dt[ix] = new long long int*[resolution[RN_Y]];
        for(int iy = 0; iy < resolution[RN_Y]; ++iy) {
            dt[ix][iy] = new long long int[resolution[RN_Z]];
            for(int iz = 0; iz < resolution[RN_Z]; ++iz) {
                dt[ix][iy][iz] = 0;
            }
        }
    }

    // get the distance transform along z
    for(int ix = 0; ix < resolution[RN_X]; ++ix) {
        for(int iy = 0; iy < resolution[RN_Y]; ++iy) {

            /* for every fixed row and column, compute DT */

            int k = 0;
            int* v = new int[resolution[RN_Z]];
            RNScalar* z = new RNScalar[resolution[RN_Z]];

            v[0] = 0;
            z[0] = -1 * infinity;
            z[1] = infinity;

            for(int q = 1; q < resolution[RN_Z]; ++q) {
            zlabel:
                RNScalar s = ((b[ix][iy][q] + q * q) - (b[ix][iy][v[k]] + v[k] * v[k])) / (RNScalar)(2 * q - 2 * v[k]);
                if(s <= z[k]) {
                    k = k - 1;
                    goto zlabel;
                } else {
                    k = k + 1;
                    v[k] = q;
                    z[k] = s;
                    z[k + 1] = infinity;
                }
            }

            k = 0;
            for(int q = 0; q < resolution[RN_Z]; ++q) {
                while(z[k + 1] < q) {
                    k = k + 1;
                }

                dt[ix][iy][q] = (q - v[k]) * (q - v[k]) + b[ix][iy][v[k]];
                zpass[ix][iy][q] = v[k];
            }

            delete[] v;
            delete[] z;
        }
    }

    // update b
    for(int ix = 0; ix < resolution[RN_X]; ++ix) {
        for(int iy = 0; iy < resolution[RN_Y]; ++iy) {
            for(int iz = 0; iz < resolution[RN_Z]; ++iz) {
                b[ix][iy][iz] = dt[ix][iy][iz];
            }
        }
    }

    // get the distance transform along y
    for(int ix = 0; ix < resolution[RN_X]; ++ix) {
        for(int iz = 0; iz < resolution[RN_Z]; ++iz) {

            /* for every fixed row and aisle, compute DT */

            int k = 0;
            int* v = new int[resolution[RN_Y]];
            RNScalar* z = new RNScalar[resolution[RN_Y]];
            v[0] = 0;
            z[0] = -1 * infinity;
            z[1] = infinity;

            for(int q = 1; q < resolution[RN_Y]; ++q) {
            ylabel:
                RNScalar s = ((b[ix][q][iz] + q * q) - (b[ix][v[k]][iz] + v[k] * v[k])) / (RNScalar)(2 * q - 2 * v[k]);
                if(s <= z[k]) {
                    k = k - 1;
                    goto ylabel;
                } else {
                    k = k + 1;
                    v[k] = q;
                    z[k] = s;
                    z[k + 1] = infinity;
                    ;
                }
            }

            k = 0;
            for(int q = 0; q < resolution[RN_Y]; ++q) {
                while(z[k + 1] < q) {
                    k = k + 1;
                }

                dt[ix][q][iz] = (q - v[k]) * (q - v[k]) + b[ix][v[k]][iz];
                ypass[ix][q][iz] = v[k];
            }

            delete[] v;
            delete[] z;
        }
    }

    // update b
    for(int ix = 0; ix < resolution[RN_X]; ++ix) {
        for(int iy = 0; iy < resolution[RN_Y]; ++iy) {
            for(int iz = 0; iz < resolution[RN_Z]; ++iz) {
                b[ix][iy][iz] = dt[ix][iy][iz];
            }
        }
    }

    // get the distance transform along x
    for(int iy = 0; iy < resolution[RN_Y]; ++iy) {
        for(int iz = 0; iz < resolution[RN_Z]; ++iz) {

            /* for every fixed column and aisle, compute DT */

            int k = 0;
            int* v = new int[resolution[RN_X]];
            RNScalar* z = new RNScalar[resolution[RN_X]];
            v[0] = 0;
            z[0] = -1 * infinity;
            z[1] = infinity;

            for(int q = 1; q < resolution[RN_X]; ++q) {
            xlabel:
                RNScalar s = ((b[q][iy][iz] + q * q) - (b[v[k]][iy][iz] + v[k] * v[k])) / (RNScalar)(2 * q - 2 * v[k]);
                if(s <= z[k]) {
                    k = k - 1;
                    goto xlabel;
                } else {
                    k = k + 1;
                    v[k] = q;
                    z[k] = s;
                    z[k + 1] = infinity;
                }
            }

            k = 0;
            for(int q = 0; q < resolution[RN_X]; ++q) {
                while(z[k + 1] < q) {
                    k = k + 1;
                }

                dt[q][iy][iz] = (q - v[k]) * (q - v[k]) + b[v[k]][iy][iz];
                xpass[q][iy][iz] = v[k];
            }

            delete[] v;
            delete[] z;
        }
    }

    printf("Created distance transform in %0.2f seconds.\n", start_time.Elapsed());

    start_time.Read();

    // create boundary points
    R3Point*** boundary = new R3Point**[resolution[RN_X]];
    R3Vector*** vectors = new R3Vector**[resolution[RN_X]];
    for(int ix = 0; ix < resolution[RN_X]; ++ix) {
        boundary[ix] = new R3Point*[resolution[RN_Y]];
        vectors[ix] = new R3Vector*[resolution[RN_Y]];
        for(int iy = 0; iy < resolution[RN_Y]; ++iy) {
            boundary[ix][iy] = new R3Point[resolution[RN_Z]];
            vectors[ix][iy] = new R3Vector[resolution[RN_Z]];
            for(int iz = 0; iz < resolution[RN_Z]; ++iz) {
                int xvalue = xpass[ix][iy][iz];
                int yvalue = ypass[xvalue][iy][iz];
                int zvalue = zpass[xvalue][yvalue][iz];

                boundary[ix][iy][iz] = R3Point(xvalue, yvalue, zvalue);
                vectors[ix][iy][iz] = boundary[ix][iy][iz] - R3Point(ix, iy, iz);
                vectors[ix][iy][iz].Normalize();
            }
        }
    }

    // free border and pass information
    for(int ix = 0; ix < resolution[RN_X]; ++ix) {
        for(int iy = 0; iy < resolution[RN_Y]; ++iy) {
            delete[] b[ix][iy];
            delete[] xpass[ix][iy];
            delete[] ypass[ix][iy];
            delete[] zpass[ix][iy];
        }
        delete[] b[ix];
        delete[] xpass[ix];
        delete[] ypass[ix];
        delete[] zpass[ix];
    }
    delete[] b;
    delete[] xpass;
    delete[] ypass;
    delete[] zpass;

    printf("Determined closest boundary point in %0.2f seconds.\n", start_time.Elapsed());

    start_time.Read();

    // create h5 output file
    R3Grid* distance_transform = new R3Grid(resolution[RN_X], resolution[RN_Y], resolution[RN_Z]);
    for(int iz = 0; iz < resolution[RN_Z]; ++iz) {
        for(int iy = 0; iy < resolution[RN_Y]; ++iy) {
            for(int ix = 0; ix < resolution[RN_X]; ++ix) {
                distance_transform->SetGridValue(ix, iy, iz, sqrt(dt[ix][iy][iz]));
            }
        }
    }

    for(int ix = 0; ix < resolution[RN_X]; ++ix) {
        for(int iy = 0; iy < resolution[RN_Y]; ++iy) {
            delete[] dt[ix][iy];
        }
        delete[] dt[ix];
    }
    delete[] dt;

    // test results
    if(ntests > 0) TestDistanceTransform(distance_transform, boundary);

    printf("Ran %d tests with no errors in %0.2f seconds.\n", ntests, start_time.Elapsed());

    start_time.Read();

    R3Grid* skeleton = new R3Grid(resolution[RN_X], resolution[RN_Y], resolution[RN_Z]);
    for(int iz = 0; iz < resolution[RN_Z]; ++iz) {
        for(int iy = 0; iy < resolution[RN_Y]; ++iy) {
            for(int ix = 0; ix < resolution[RN_X]; ++ix) {
                // traverse the vector one unit backwards
                R3Point neighbor = R3Point(ix, iy, iz) - vectors[ix][iy][iz];

                // get the vector at neighbor with interpolation
                int xhigh = (int)(ceil(neighbor.X()) + 0.5);
                int yhigh = (int)(ceil(neighbor.Y()) + 0.5);
                int zhigh = (int)(ceil(neighbor.Z()) + 0.5);

                int xlow = (int)(floor(neighbor.X()) + 0.5);
                int ylow = (int)(floor(neighbor.Y()) + 0.5);
                int zlow = (int)(floor(neighbor.Z()) + 0.5);

                if (xhigh > resolution[RN_X] - 1) continue;
                if (yhigh > resolution[RN_Y] - 1) continue;
                if (zhigh > resolution[RN_Z] - 1) continue;
                if (xhigh < 0) continue;
                if (yhigh < 0) continue;
                if (zhigh < 0) continue;

                RNScalar alphax = neighbor.X() - floor(neighbor.X());
                RNScalar alphay = neighbor.Y() - floor(neighbor.Y());
                RNScalar alphaz = neighbor.Z() - floor(neighbor.Z());

                R3Vector neighbor_vector = alphax * alphay * alphaz * vectors[xlow][ylow][zlow] +
                    alphax * alphay * (1 - alphaz) * vectors[xlow][ylow][zhigh] +
                    alphax * (1 - alphay) * alphaz * vectors[xlow][yhigh][zlow] +
                    alphax * (1 - alphay) * (1 - alphaz) * vectors[xlow][yhigh][zhigh] +
                    (1 - alphax) * alphay * alphaz * vectors[xhigh][ylow][zlow] +
                    (1 - alphax) * alphay * (1 - alphaz) * vectors[xhigh][ylow][zhigh] +
                    (1 - alphax) * (1 - alphay) * alphaz * vectors[xhigh][yhigh][zlow] +
                    (1 - alphax) * (1 - alphay) * (1 - alphaz) * vectors[xhigh][yhigh][zhigh];
                
                // normalize the neighbor vector
                neighbor_vector.Normalize();
                
                // get the dot product between the vectors
                RNScalar dot_product = neighbor_vector.Dot(vectors[ix][iy][iz]);
                
                // if the dot product is negative, add to skeleton
                if (dot_product < 0) skeleton->SetGridValue(ix, iy, iz, machine_labels->GridValue(ix, iy, iz));
                else skeleton->SetGridValue(ix, iy, iz, 0);
            }
        }
    }

    printf("Created skeletons in %0.2f seconds\n", start_time.Elapsed());

    start_time.Read();

    // get the root filename
    char root_filename[4096];
    strncpy(root_filename, input_filename, 4096);
    char *extp = strrchr(root_filename, '.');
    *extp = '\0';
    char *filename = strrchr(root_filename, '/');
    filename++;
    
    // output the distnace transform and the skeleton
    char dt_filename[4096];
    sprintf(dt_filename, "skeletons/%s_dt.h5", filename);
    char skel_filename[4096];
    sprintf(skel_filename, "skeletons/%s_skeleton.h5", filename);

    if(!RNWriteH5File(&distance_transform, 1, dt_filename, "main", FALSE)) {
        fprintf(stderr, "Failed to write distance transform to %s.\n", dt_filename);
        return 0;
    }
    if (!RNWriteH5File(&skeleton, 1, skel_filename, "main", TRUE)) {
        fprintf(stderr, "Failed to write skeleton to %s.\n", skel_filename);
        return 0;
    }

    printf("Output to h5 file in %0.2f seconds.\n", start_time.Elapsed());

    // free memory
    for(int ix = 0; ix < resolution[RN_X]; ++ix) {
        for(int iy = 0; iy < resolution[RN_Y]; ++iy) {
            delete[] vectors[ix][iy];
            delete[] boundary[ix][iy];
        }
        delete[] vectors[ix];
        delete[] boundary[ix];
    }
    delete[] vectors;
    delete[] boundary;
    delete machine_labels;
    delete distance_transform;
    delete skeleton;

    printf("Skeletonization completed in %0.2f seconds\n", total_time.Elapsed());

    // return success
    return 0;
}
