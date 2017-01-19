// Source file for neuron visualizer



// include files

#include "RNDataStructures/RNDataStructures.h"
#include <vector>
#include <fstream>
#include <string>



// GLUT defines

#define ENTER 13
#define ESCAPE 27
#define SPACEBAR 32



// program arguments

static const char *machine_labels_filename = NULL;
static const char *machine_labels_dataset = NULL;
static char *prefix = NULL;
static int print_verbose = 0;
static int print_debug = 0;



// voxel grids

static R3Grid *machine_labels_grid = NULL;



// struct name definitions

struct SWCEntry;



// display machine label variables 

static std::vector<R3Point> *swc_endpoints = NULL;
static std::vector<SWCEntry> *swc_skeleton = NULL;
static int *index_to_label = NULL;



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
    {}

    RNScalar X(void) { return x_position; }
    RNScalar Y(void) { return y_position; }
    RNScalar Z(void) { return z_position; }

    // instance variables
    int sample_number; 
    int structure_identifier;
    RNScalar x_position;
    RNScalar y_position;
    RNScalar z_position;
    RNScalar radius;
    int parent_sample;
};



static int ReadSWCFile(int index, int label)
{
    char input_filename[4096];
    sprintf(input_filename, "skeletons/%s/tree_%d.swc", prefix, label);

    std::ifstream fd(input_filename);
    if (!fd.is_open()) {
        fprintf(stderr, "Failed to read %s\n", input_filename);
        return 0;
    }

    std::string line;
    RNBoolean first_iteration = TRUE;
    while (std::getline(fd, line)) {
        if (first_iteration) { first_iteration = FALSE; continue; }

        int sample_number, structure_identifier, parent_sample;
        RNScalar x_position, y_position, z_position, radius;
        sscanf(line.c_str(), "%d %d %lf %lf %lf %lf %d", &sample_number, &structure_identifier, &x_position, &y_position, &z_position, &radius, &parent_sample);

        swc_skeleton[index].push_back(SWCEntry(sample_number, structure_identifier, x_position, y_position, z_position, radius, parent_sample));
    }
    fd.close();

    // determine which points are endpoints
    RNBoolean *endpoints = new RNBoolean[swc_skeleton[index].size()];
    for (unsigned int ie = 0; ie < swc_skeleton[index].size(); ++ie) {
        endpoints[ie] = TRUE;
    }

    // not an endpoint if no children claim it as a parent
    for (unsigned int ie = 0; ie < swc_skeleton[index].size(); ++ie) {
        int parent_sample = swc_skeleton[index][ie].parent_sample;
        if (parent_sample == -1) continue;
        endpoints[parent_sample - 1] = FALSE;
    }
    
    // if it has no parent it is an endpoint
    for (unsigned int ie = 0; ie < swc_skeleton[index].size(); ++ie) {
        int parent_sample = swc_skeleton[index][ie].parent_sample;
        if (parent_sample == -1) endpoints[ie] = TRUE;
    }
    
    swc_endpoints[index] = std::vector<R3Point>();
    for (unsigned int ie = 0; ie < swc_skeleton[index].size(); ++ie) {
        if (endpoints[ie]) {
            SWCEntry entry = swc_skeleton[index][ie];
            swc_endpoints[index].push_back(R3Point(entry.X(), entry.Y(), entry.Z()));
        }
    }
    
    // return success
    return 1;
}


static int ReadData(void)
{
   // start statistics
   RNTime start_time;
   start_time.Read();

   if (!machine_labels_filename) {
       char *filename = new char[4096];
       sprintf(filename, "machine_labels/%s_machine_labels_seg_28000.h5", prefix);
       machine_labels_filename = filename;
       machine_labels_dataset = "main";
   }

   // read in voxel files
   R3Grid **machine_labels = RNReadH5File(machine_labels_filename, machine_labels_dataset);
   if (!machine_labels) { fprintf(stderr, "Failed to read %s from %s\n", machine_labels_dataset, machine_labels_filename); return 0; }
   machine_labels_grid = machine_labels[0];
   delete[] machine_labels;

   // print statistics
   if (print_verbose) {
      printf("Read voxel grids...\n");
      printf("  Time = %.2f seconds\n", start_time.Elapsed());
   }

   // return success
   return 1;
}


static void ComputeDistances(void)
{
    int max_value = (int)(machine_labels_grid->Maximum() + 0.5) + 1;

    RNBoolean *element_found = new RNBoolean[max_value]; 
    for (int ie = 0; ie < max_value; ++ie) {
        element_found[ie] = FALSE;
    }
    for (int iv = 0; iv < machine_labels_grid->NEntries(); ++iv) {
        int grid_value = (int)(machine_labels_grid->GridValue(iv) + 0.5);
        element_found[grid_value] = TRUE;
    }

    // count the number of unique elements
    int nunique = 0;
    for (int ie = 0; ie < max_value; ++ie) {
        if (element_found[ie]) nunique++;
    }

    // create mapping variables
    index_to_label = new int[nunique];
    int seen = 0;
    for (int ie = 0; ie < max_value; ++ie) {
        if (element_found[ie]) {
            index_to_label[seen] = ie;
            ++seen;
        }
    }

    swc_skeleton = new std::vector<SWCEntry>[nunique];
    swc_endpoints = new std::vector<R3Point>[nunique];

    int *nboundary_elements = new int[nunique];
    for (int ie = 0; ie < nunique; ++ie) {
        swc_skeleton[ie] = std::vector<SWCEntry>();
        nboundary_elements[ie] = 0;
    }

    // try to create an SWC vector for each label
    for (int ie = 0; ie < nunique; ++ie) {
        ReadSWCFile(ie, index_to_label[ie]);
    }

    // compute the distances between endpoints
    int nendpoints = 0;
    std::vector<R3Point> all_endpoints = std::vector<R3Point>();
    for (int ii = 0; ii < nunique; ++ii) {
        nendpoints += swc_endpoints[ii].size();
        for (unsigned int ie = 0; ie < swc_endpoints[ii].size(); ++ie) {
            all_endpoints.push_back(swc_endpoints[ii][ie]);
        }
    }
    
    // compute all pairs of distances
    float **distances = new float *[nendpoints];
    for (int ie = 0; ie < nendpoints; ++ie) {
        distances[ie] = new float[nendpoints];
    }
    
    int *labels = new int[nendpoints];
    int index = 0;
    for (int ii = 0; ii < nunique; ++ii) {
        for (unsigned int ie = 0; ie < swc_endpoints[ii].size(); ++ie, ++index) {
            labels[index] = ii;
        }
    }
    
    for (unsigned int i = 0; i < all_endpoints.size(); ++i) {
        for (unsigned int j = 0; j < all_endpoints.size(); ++j) {
            distances[i][j] = R3Distance(all_endpoints[i], all_endpoints[j]);
        }
    }
    
    // write to file
    char output_filename[4096];
    sprintf(output_filename, "skeletons/%s_skeleton.dist", prefix);
    
    // open file
    FILE *fp = fopen(output_filename, "wb");
    if (!fp) { fprintf(stderr, "Failed to write to %s\n", output_filename); exit(-1); }
    
    // write the number of endpoints
    fwrite(&nendpoints, sizeof(unsigned int), 1, fp);
    fwrite(&(labels[0]), sizeof(int), nendpoints, fp);
    for (unsigned int i = 0; i < all_endpoints.size(); ++i) {
        fwrite(&(distances[i]), sizeof(float), all_endpoints.size(), fp);
    }
    
    // close file
    fclose(fp);
    
    printf("%d\n", nendpoints);

    // free memory
    delete[] nboundary_elements;
    delete[] element_found;
}



////////////////////////////////////////////////////////////////////////
// Program argument parsing
////////////////////////////////////////////////////////////////////////

static int ParseArgs(int argc, char **argv)
{
   // parse arguments
   argc--; argv++;
   while (argc > 0) {
      if ((*argv)[0] == '-') {
         if (!strcmp(*argv, "-v")) print_verbose = 1;
         else if (!strcmp(*argv, "-debug")) print_debug = 1;
         else if (!strcmp(*argv, "-machine_labels")) {
             argv++; argc--;
             machine_labels_filename = *argv;
             argv++; argc--;
             machine_labels_dataset = *argv;
         }
         else { fprintf(stderr, "Invalid program argument: %s\n", *argv); return 0; }
      }
      else {
         fprintf(stderr, "Invalid program argument: %s\n", *argv); return 0;
      }
      argv++; argc--;
   }

   // error if there is no input name
   if (!machine_labels_filename) { fprintf(stderr, "Need to supply a neuron data file\n"); return 0; }

   // return success
   return 1;
}



////////////////////////////////////////////////////////////////////////
// Main program
////////////////////////////////////////////////////////////////////////

int main(int argc, char **argv)
{
   // parse program arguments
   if (!ParseArgs(argc, argv)) exit(-1);

    // get the prefix
    char *filename;
    filename = strrchr((char *)machine_labels_filename, '/');
    filename++;
    char root_filename[4096];
    strncpy(root_filename, filename, 4096);
    char *extp = strrchr(root_filename, '.');
    *extp = '\0';
    extp = strchr(root_filename, '_');
    if (extp) *extp = '\0';
    prefix = root_filename;
    

   /////////////////////////////////
   //// Read in the voxel files ////
   /////////////////////////////////

   if (!ReadData()) exit(-1);

   ComputeDistances();

   // return success
   return 1;
}