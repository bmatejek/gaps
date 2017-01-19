// Source file for neuron visualizer



// include files

#include "RNDataStructures/RNDataStructures.h"
#include "R3Graphics/R3Graphics.h"
#include "fglut/fglut.h"
#include <vector>
#include <fstream>
#include <string>



// GLUT defines

#define ENTER 13
#define ESCAPE 27
#define SPACEBAR 32



// program arguments

static int print_debug = 0;
static int print_verbose = 0;
static const char *prefix = NULL;
static const char *affinities_filename = NULL;
static const char *affinities_dataset = NULL;
static const char *image_filename = NULL;
static const char *image_dataset = NULL;
static const char *gold_filename = NULL;
static const char *gold_dataset = NULL;
static const char *segmentation_filename = NULL;
static const char *segmentation_dataset = NULL;



// program variables

static RNScalar scaling[3] = {1, 1, 5};
static int resolution[3] = {-1, -1, -1};
static R3Affine transformation = R3null_affine;
static R3Viewer *viewer = NULL;
static R3Point selected_position;
static R3Box world_box;
static int selected_voxel[3];



// voxel grids

static R3Grid *affinity_grid[3] = { NULL, NULL, NULL };
static R3Grid *gold_grid = NULL;
static R3Grid *image_grid = NULL;
static R3Grid *segmentation_grid = NULL;



// projection variables

static R2Grid *affinity_selected_slice[3] = { NULL, NULL, NULL };
static R2Grid *gold_selected_slice = NULL;
static R2Grid *image_selected_slice = NULL;
static R2Grid *segmentation_selected_slice = NULL;
static RNInterval affinity_range[3];
static RNInterval gold_range;
static RNInterval image_range;
static RNInterval segmentation_range;
// always update the grid and the interval
static R2Grid *selected_slice = NULL;
static R2Point selected_slice_position(RN_UNKNOWN, RN_UNKNOWN);
static RNInterval selected_slice_range(0, 0);
static RNScalar projection_scale = FLT_MAX;



// GLUT variables

static int GLUTwindow = 0;
static int GLUTwindow_height = 1200;
static int GLUTwindow_width = 1200;
static int GLUTmouse[2] = { 0, 0 };
static int GLUTbutton[3] = { 0, 0, 0 };
static int GLUTmodifiers = 0;



// color arrays

static RNScalar background_color[] = { 0, 0, 0 };



// projection color variables

static int color_type = 0; // 0=gray, 1=red-green-blue



// display projection variables

static int show_projection_affinity = 0;
static int show_projection_image = 1;
static int show_projection_gold = 0;
static int show_projection_segmentation = 0;
static int projection_dim = RN_Z;
static int affinity_dim = RN_X;
static int selected_slice_index = 0;



// display dimenson variables

static int show_bbox = 1;
static int show_slice = 0;
static int projection = 0;
static int show_dimension_affinity = 0;
static int show_dimension_image = 1;
static int show_dimension_gold = 0;
static int show_dimension_segmentation = 0;



// struct name definitions

struct SWCEntry;



// display machine label variables 

static int segmentation_index = 0;
static int downsample_rate = 4;
static std::vector<R3Point> *segmentations = NULL;
static std::vector<int> *false_splits = NULL;



// mapping functions

static int *segmentation_to_gold = NULL;
static int *index_to_segmentation = NULL;
static int *segmentation_to_index = NULL;



////////////////////////////////////////////////////////////////////////
// Input/output functions
////////////////////////////////////////////////////////////////////////

static int ReadData(void)
{
   // start statistics
   RNTime start_time;
   start_time.Read();

   // create default filenames
   if (!affinities_filename) {
       char *filename = new char[4096];
       sprintf(filename, "affinities/%s_affinities.h5", prefix);
       affinities_filename = filename;
       affinities_dataset = "main";
   }
   if (!gold_filename) {
       char *filename = new char[4096];
       sprintf(filename, "gold/%s_gold.h5", prefix);
       gold_filename = filename;
       gold_dataset = "stack";
   }
   if (!image_filename) {
       char *filename = new char[4096];
       sprintf(filename, "images/%s_image.h5", prefix);
       image_filename = filename;
       image_dataset = "main";
   }
   if (!segmentation_filename) {
       char *filename = new char[4096];
       sprintf(filename, "segmentations/%s_machine_labels_seg_28000.h5", prefix);
       segmentation_filename = filename;
       segmentation_dataset = "main";
   }



   // read in voxel files
   R3Grid **affinities = RNReadH5File(affinities_filename, affinities_dataset);
   if (!affinities) { fprintf(stderr, "Failed to read %s from %s\n", affinities_dataset, affinities_filename); return 0; }
   for (int dim = 0; dim <= 2; ++dim) {
       affinity_grid[dim] = affinities[dim];
   }
   delete[] affinities;
   
   R3Grid **golds = RNReadH5File(gold_filename, gold_dataset);
   if (!golds) { fprintf(stderr, "Failed to read %s from %s\n", gold_dataset, gold_filename); return 0; }
   gold_grid = golds[0];
   delete[] golds;

   R3Grid **images = RNReadH5File(image_filename, image_dataset);
   if (!images) { fprintf(stderr, "Failed to read %s from %s\n", image_dataset, image_filename); return 0; }
   image_grid = images[0];
   delete[] images;

   R3Grid **segmentations = RNReadH5File(segmentation_filename, segmentation_dataset);
   if (!segmentations) { fprintf(stderr, "Failed to read %s from %s\n", segmentation_dataset, segmentation_filename); return 0; }
   segmentation_grid = segmentations[0];
   delete[] segmentations;
   
   

   // make sure the resolutions are the same
   rn_assertion(affinity_grid[RN_X]->XResolution() == affinity_grid[RN_Y]->XResolution());
   rn_assertion(affinity_grid[RN_Y]->XResolution() == affinity_grid[RN_Z]->XResolution());
   rn_assertion(affinity_grid[RN_Z]->XResolution() == gold_grid->XResolution());
   rn_assertion(gold_grid->XResolution() == image_grid->XResolution());
   rn_assertion(image_grid->XResolution() == segmentation_grid->XResolution());
   rn_assertion(affinity_grid[RN_X]->YResolution() == affinity_grid[RN_Y]->YResolution());
   rn_assertion(affinity_grid[RN_Y]->YResolution() == affinity_grid[RN_Z]->YResolution());
   rn_assertion(affinity_grid[RN_Z]->YResolution() == gold_grid->YResolution());
   rn_assertion(gold_grid->YResolution() == image_grid->YResolution());
   rn_assertion(image_grid->YResolution() == segmentation_grid->YResolution());
   rn_assertion(affinity_grid[RN_X]->ZResolution() == affinity_grid[RN_Y]->ZResolution());
   rn_assertion(affinity_grid[RN_Y]->ZResolution() == affinity_grid[RN_Z]->ZResolution());
   rn_assertion(affinity_grid[RN_Z]->ZResolution() == gold_grid->ZResolution());
   rn_assertion(gold_grid->ZResolution() == image_grid->ZResolution());
   rn_assertion(image_grid->ZResolution() == segmentation_grid->ZResolution());

    resolution[RN_X] = image_grid->XResolution();
    resolution[RN_Y] = image_grid->YResolution();
    resolution[RN_Z] = image_grid->ZResolution();

   // print statistics
   if (print_verbose) {
      printf("Read voxel grids...\n");
      printf("  Time = %.2f seconds\n", start_time.Elapsed());
      printf("  Resolution = (%d %d %d)\n", gold_grid->XResolution(), gold_grid->YResolution(), gold_grid->ZResolution());
      printf("  Scaling = (%lf %lf %lf)\n", scaling[RN_X], scaling[RN_Y], scaling[RN_Z]);
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
    int **segmentation_gold_overlap = new int *[max_segmentation];
    if (!segmentation_gold_overlap) return 0;
    for (int is = 0; is < max_segmentation; ++is) {
        segmentation_gold_overlap[is] = new int[max_gold];
        if (!segmentation_gold_overlap[is]) return 0;
        for (int ig = 0; ig < max_gold; ++ig) {
            segmentation_gold_overlap[is][ig] = 0;
        }
    }
    
    // iterate through all voxels
    for (int iv = 0; iv < segmentation_grid->NEntries(); ++iv) {
        int segmentation_id = (int)(segmentation_grid->GridValue(iv) + 0.5);
        int gold_id = (int)(gold_grid->GridValue(iv) + 0.5);
        
        segmentation_gold_overlap[segmentation_id][gold_id]++;
    }
    
    // create a mapping from segmentation ids to gold
    segmentation_to_gold = new int[max_segmentation];
    for (int is = 0; is < max_segmentation; ++is) {
        int max_gold_id = -1;
        int max_gold_value = -1;
        // skip the extracellular label
        for (int ig = 1; ig < max_gold; ++ig) {
            if (segmentation_gold_overlap[is][ig] > max_gold_value) {
                max_gold_value = segmentation_gold_overlap[is][ig];
                max_gold_id = ig;
            }
        }
        segmentation_to_gold[is] = max_gold_id;
    }
    
    // free memory
    for (int is = 0; is < max_segmentation; ++is)
        delete[] segmentation_gold_overlap[is];
    delete[] segmentation_gold_overlap;
    
    // return success
    return 1;
}



static int MapIndexToSegmentation(void)
{  
    // get the number of segmentations
    int max_segmentation = (int)(segmentation_grid->Maximum() + 0.5) + 1;
    
    // create boolean array of id existence
    RNBoolean *id_exists = new RNBoolean[max_segmentation];
    if (!id_exists) return 0;
    for (int is = 0; is < max_segmentation; ++is)
        id_exists[is] = FALSE;
    
    for (int iv = 0; iv < segmentation_grid->NEntries(); ++iv) {
        int segmentation_id = (int)(segmentation_grid->GridValue(iv) + 0.5);
        id_exists[segmentation_id] = TRUE;
    }
    
    // count the number of unique occurrences
    int nunique_segmentations = 0;
    for (int is = 0; is < max_segmentation; ++is) {
        if (id_exists[is]) nunique_segmentations++;
    }
    
    // create the mapping from indices to segmentations
    index_to_segmentation = new int[nunique_segmentations];
    if (!index_to_segmentation) return 0;
    segmentation_to_index = new int[max_segmentation];
    if (!segmentation_to_index) return 0;
    int number_visited = 0;
    for (int is = 0; is < max_segmentation; ++is) {
        if (id_exists[is]) {
            index_to_segmentation[number_visited] = is;
            segmentation_to_index[is] = number_visited;
            number_visited++;
        }
        else {
            segmentation_to_index[is] = -1;
        }
    }
    
    // allocate memory for all segmentation vectors
    segmentations = new std::vector<R3Point>[nunique_segmentations];
    if (!segmentations) return 0;
    false_splits = new std::vector<int>[nunique_segmentations];
    if (!false_splits) return 0;
    
    for (int is = 0; is < nunique_segmentations; ++is) {
        segmentations[is] = std::vector<R3Point>();
        false_splits[is] = std::vector<int>();
    }
    
    // print out the number of unique gold ids
    int max_gold = (int)(gold_grid->Maximum() + 0.5) + 1;
    RNBoolean *gold_id_exists = new RNBoolean[max_gold];
    for (int ig = 0; ig < max_gold; ++ig) {
        gold_id_exists[ig] = FALSE;
    }
    for (int iv = 0; iv < gold_grid->NEntries(); ++iv) {
        int gold_id = (int)(gold_grid->GridValue(iv) + 0.5);
        gold_id_exists[gold_id] = TRUE;
    }
    int nunique_gold = 0;
    for (int ig = 0; ig < max_gold; ++ig) {
        if (gold_id_exists[ig]) nunique_gold++;
    }
    
    if (print_verbose) {
        printf("  NSegments = %d\n", nunique_segmentations);
        printf("  NLabels = %d\n", nunique_gold);
    }
    
    // return success
    return 1;
}



static int PopulateSegmentationVectors(void)
{
    int max_segmentation = (int)(segmentation_grid->Maximum() + 0.5) + 1;
    int *nvoxels_per_segments = new int[max_segmentation];
    if (!nvoxels_per_segments) return 0;
    for (int is = 0; is < max_segmentation; ++is)
        nvoxels_per_segments[is] = 0;
    
    // go through all voxels
    for (int ix = 0; ix < segmentation_grid->XResolution(); ++ix) {
        for (int iy = 0; iy < segmentation_grid->YResolution(); ++iy) {
            for (int iz = 0; iz < segmentation_grid->ZResolution(); ++iz) {
                // get the segmentation and grid value
                int segmentation_id = (int)(segmentation_grid->GridValue(ix, iy, iz) + 0.5);
                
                // get the index - just checking...
                int index = segmentation_to_index[segmentation_id];
                rn_assertion(index != -1);
                
                
                // see if this element is on a boundary
                RNBoolean boundary = FALSE;
                if (ix == 0 || iy == 0 || iz == 0) boundary = TRUE;
                else if (ix == resolution[RN_X] - 1 || iy == resolution[RN_Y] - 1 || iz == resolution[RN_Z] - 1) boundary = TRUE;
                else {                        
                    int south_id = (int)(segmentation_grid->GridValue(ix - 1, iy, iz) + 0.5);
                    if (south_id != segmentation_id) boundary = TRUE;
                    int north_id = (int)(segmentation_grid->GridValue(ix + 1, iy, iz) + 0.5);
                    if (north_id != segmentation_id) boundary = TRUE;
                    int east_id = (int)(segmentation_grid->GridValue(ix, iy - 1, iz) + 0.5);
                    if (east_id != segmentation_id) boundary = TRUE;
                    int west_id = (int)(segmentation_grid->GridValue(ix, iy + 1, iz) + 0.5);  
                    if (west_id != segmentation_id) boundary = TRUE;
                    int down_id = (int)(segmentation_grid->GridValue(ix, iy, iz - 1) + 0.5);
                    if (down_id != segmentation_id) boundary = TRUE;
                    int up_id = (int)(segmentation_grid->GridValue(ix, iy, iz + 1) + 0.5);
                    if (up_id != segmentation_id) boundary = TRUE;
                }
                
                // update the number of voxels per segment
                nvoxels_per_segments[segmentation_id]++;
                if ((nvoxels_per_segments[segmentation_id] % downsample_rate) != 0) continue;
                
                // only draw boundaries
                if (!boundary) continue;
                
                // get the segment gold id               
                segmentations[index].push_back(R3Point(ix, iy, iz));
            }
        }
    }
    
    // update vector of false splits
    for (int is1 = 0; is1 < max_segmentation; ++is1) {
        if (segmentation_to_index[is1] == -1) continue;
        for (int is2 = is1 + 1; is2 < max_segmentation; ++is2) {
            if (segmentation_to_index[is2] == -1) continue;
            
            if (segmentation_to_gold[is1] == segmentation_to_gold[is2]) {
                int index_one = segmentation_to_index[is1];
                int index_two = segmentation_to_index[is2];
                false_splits[index_one].push_back(index_two);
                false_splits[index_two].push_back(index_one);
            }
        }
    }

    // free memory
    delete[] nvoxels_per_segments;
    
    // return success
    return 1;
}



////////////////////////////////////////////////////////////////////////
// Utility functions
////////////////////////////////////////////////////////////////////////

static R2Point ConvertGridToWorld(R2Point grid_point)
{
   // projection offsets
   int window_xdim = (projection_dim + 1) % 3;
   int window_ydim = (projection_dim + 2) % 3;
   int xoffset = (GLUTwindow_width - scaling[window_xdim] * projection_scale * image_grid->Resolution(window_xdim)) / 2;
   int yoffset = (GLUTwindow_height - scaling[window_ydim] * projection_scale * image_grid->Resolution(window_ydim)) / 2;

   // convert from grid to world coordinates
   return R2Point(grid_point.X() * projection_scale * scaling[window_xdim] + xoffset, grid_point.Y() * projection_scale * scaling[window_ydim] + yoffset);
}



static R2Point ConvertWorldToGrid(R2Point world_point)
{
   // projection offsets
   int window_xdim = (projection_dim + 1) % 3;
   int window_ydim = (projection_dim + 2) % 3;
   int xoffset = (GLUTwindow_width - scaling[window_xdim] * projection_scale * image_grid->Resolution(window_xdim)) / 2;
   int yoffset = (GLUTwindow_height - scaling[window_ydim] * projection_scale * image_grid->Resolution(window_ydim)) / 2;

   // convert from world to grid coordinates
   return R2Point((world_point.X() - xoffset) / (scaling[window_xdim] * projection_scale), (world_point.Y() - yoffset) / (scaling[window_ydim] * projection_scale));
}



static void UpdateSlices(void)
{
   // update selected affinities
   for (int dim = 0; dim <= 2; ++dim) {
      affinity_selected_slice[dim] = affinity_grid[dim]->Slice(projection_dim, selected_slice_index);
      affinity_range[dim] = affinity_selected_slice[dim]->Range();
   }

   // update selected images
   image_selected_slice = image_grid->Slice(projection_dim, selected_slice_index);
   image_range = image_selected_slice->Range();

   // update selected machine labels
   segmentation_selected_slice = segmentation_grid->Slice(projection_dim, selected_slice_index);
   segmentation_range = segmentation_selected_slice->Range();

   // update selected truths
   gold_selected_slice = gold_grid->Slice(projection_dim, selected_slice_index);
   gold_range = gold_selected_slice->Range();

   // update pointers to selected slice
   if (show_projection_affinity) {
      selected_slice = affinity_selected_slice[projection_dim];
      selected_slice_range = affinity_range[projection_dim];
   }
   else if (show_projection_image) {
      selected_slice = image_selected_slice;
      selected_slice_range = image_range;
   }
   else if (show_projection_segmentation) {
      selected_slice = segmentation_selected_slice;
      selected_slice_range = segmentation_range;
   }
   else if (show_projection_gold) {
      selected_slice = gold_selected_slice;
      selected_slice_range = gold_range;
   }
}



////////////////////////////////////////////////////////////////////////
// Drawing utility functions
////////////////////////////////////////////////////////////////////////

static RNRgb
Color(RNScalar value, RNInterval range)
{
   // Check for unknown value
   if (value == R2_GRID_UNKNOWN_VALUE) {
      if (color_type == 0) return RNRgb(1, 0.5, 0);
      else return RNblack_rgb;
   }


   if (show_projection_affinity || show_projection_image) {
        // Normalize color
        RNScalar value_min = range.Min();
        RNScalar value_width = range.Diameter();
        RNScalar value_scale = (value_width > 0) ? 1.0 / value_width : 1.0;
        RNScalar normalized_value = value_scale * (value - value_min);

        // Compute color
        RNRgb c(0, 0, 0);
        if (color_type == 0) {
            c[0] = normalized_value;
            c[1] = normalized_value;
            c[2] = normalized_value;
        }
        else {
            if (normalized_value < 0.5) {
                c[0] = 1 - 2 * normalized_value;
                c[1] = 2 * normalized_value;
            }
            else {
                c[1] = 1 - 2 * (normalized_value - 0.5);
                c[2] = 2 * (normalized_value - 0.5);
            }
        }

        // return color
        return c;
   }
   else {
        unsigned int index = (unsigned int) (value + 0.5);

        RNRgb color = RNRgb( ((107 * index) % 700) / 700.0, ((509 * index) % 900) / 900.0, ((200 * index) % 777) / 777.0 );

        // Return color
        return color;
   }
}



static void GLUTDrawText(const R2Point& position, const char *s)
{
   // draw text string s at position
   glRasterPos2d(position[0], position[1]);
   while (*s) glutBitmapCharacter(GLUT_BITMAP_HELVETICA_18, *(s++));
}



static void DrawIndividualSegment(int index)
{
    glBegin(GL_POINTS);
    for (unsigned int iv = 0; iv < segmentations[index].size(); ++iv) {
        glVertex3f(segmentations[index][iv].X(), segmentations[index][iv].Y(), segmentations[index][iv].Z());
    }        
    glEnd();
}


static void DrawSegmentation(void)
{
    // push the transformation
    transformation.Push();
   
    RNLoadRgb(RNblue_rgb);
    DrawIndividualSegment(segmentation_index);

    // add in false splits
    for (unsigned int iv = 0; iv < false_splits[segmentation_index].size(); ++iv) {
        int index = false_splits[segmentation_index][iv];
        RNLoadRgb(RNgreen_rgb);
        DrawIndividualSegment(index);
    }

    // pop the transformation
    transformation.Pop();
}



static void DrawSlice(void)
{
   transformation.Push();

   RNLoadRgb(1.0, 1.0, 1.0);
   // draw the selected slice
   if (show_dimension_affinity) affinity_grid[affinity_dim]->DrawSlice(projection_dim, selected_slice_index);
   if (show_dimension_gold) gold_grid->DrawSlice(projection_dim, selected_slice_index);
   if (show_dimension_image) image_grid->DrawSlice(projection_dim, selected_slice_index);
   if (show_dimension_segmentation) segmentation_grid->DrawSlice(projection_dim, selected_slice_index);

   transformation.Pop();
}



static void Draw2D(void)
{
   // prologue
   glDisable(GL_LIGHTING);

   // set projection matrix
   glMatrixMode(GL_PROJECTION);
   glPushMatrix();
   glLoadIdentity();
   gluOrtho2D(0, GLUTwindow_width, 0, GLUTwindow_height);

   // set model view matrix
   glMatrixMode(GL_MODELVIEW);
   glPushMatrix();
   glLoadIdentity();

   // draw value at selected position
   if ((selected_slice_position.X() != RN_UNKNOWN) && (selected_slice_position.Y() != RN_UNKNOWN)) {
      int ix = (int)(selected_slice_position.X() + 0.5);
      int iy = (int)(selected_slice_position.Y() + 0.5);
      RNScalar value = selected_slice->GridValue(ix, selected_slice->YResolution() - iy);

      // create hover text
      char buffer[1024];
      if (value != R2_GRID_UNKNOWN_VALUE) sprintf(buffer, "%d %d: %g", ix, iy, value);
      else sprintf(buffer, "%d %d: %s", ix, iy, "Unknown");
      RNLoadRgb(RNblue_rgb);

      // move point to world location
      R2Point world_position = ConvertGridToWorld(selected_slice_position);
      GLUTDrawText(world_position + 2 * R2ones_vector, buffer);
   }

   // draw the actual slice
   for (int iy = 1; iy < selected_slice->YResolution(); ++iy) {
      glBegin(GL_TRIANGLE_STRIP);
      for (int ix = 1; ix < selected_slice->XResolution(); ++ix) {
         for (int k = -1; k <= 0; ++k) {
            // convert from world to grid coordinates
            R2Point grid_position = R2Point(ix, selected_slice->YResolution() - (iy + k));
            R2Point world_position = ConvertGridToWorld(grid_position);

            // get the color for this area
            RNRgb color;
            if (show_projection_gold) {
               if (selected_slice->GridValue(ix, iy + k) < 0.5) color = RNred_rgb;
               else color = Color(selected_slice->GridValue(ix, iy + k), selected_slice_range);
            }
            else {
               color = Color(selected_slice->GridValue(ix, iy + k), selected_slice_range);
            }
            RNLoadRgb(color);

            // add the world position vertex
            glVertex2i(world_position.X(), world_position.Y());
         }
      }
      glEnd();
   }

   // reset projection matrix
   glMatrixMode(GL_PROJECTION);
   glPopMatrix();

   // reset model view matrix
   glMatrixMode(GL_MODELVIEW);
   glPopMatrix();

   // epilogue
   glEnable(GL_LIGHTING);
}



static void Draw3D(void)
{
   // set viewing transformation
   viewer->Camera().Load();

   // set lights
   static GLfloat light0_position[] = { 3.0, 4.0, 5.0, 0.0 };
   glLightfv(GL_LIGHT0, GL_POSITION, light0_position);
   static GLfloat light1_position[] = { -3.0, -2.0, -3.0, 0.0 };
   glLightfv(GL_LIGHT1, GL_POSITION, light1_position);

   // prologue
   glDisable(GL_LIGHTING);

   // draw neuron data bounding box
   if (show_bbox) {
      RNLoadRgb(RNwhite_rgb);
      world_box.Outline();
   }

   if (show_slice) DrawSlice();
   
   // draw machine labels and skeletons
   DrawSegmentation();

   // epilogue
   glEnable(GL_LIGHTING);
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

   //delete[] machine_label;

   // print statistics
   if (print_verbose) {
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
   if (projection) Draw2D();
   else Draw3D();

    // set window title
    char title[4096];
    sprintf(title, "Neuron Visualizer Slice %d with Machine Label %d", selected_slice_index, index_to_segmentation[segmentation_index]);
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



void GLUTMotion2D(int x, int y)
{
   if (GLUTbutton[0]) {
      // query
      R2Point world_point = R2Point(x, y);
      selected_slice_position = ConvertWorldToGrid(world_point);

      // get window dimensions
      int window_xdim = (projection_dim + 1) % 3;
      int window_ydim = (projection_dim + 2) % 3;

      // get the grid x and y coordinates
      int gridx = (int)(selected_slice_position.X() + 0.5);
      int gridy = (int)(selected_slice_position.Y() + 0.5);

      // cannot select a point outside of the image
      if (gridx < 0 || gridx >= image_grid->Resolution(window_xdim)) selected_slice_position.Reset(RN_UNKNOWN, RN_UNKNOWN);
      else if (gridy < 0 || gridy >= image_grid->Resolution(window_ydim)) selected_slice_position.Reset(RN_UNKNOWN, RN_UNKNOWN);

      glutPostRedisplay();
   }
   else if (GLUTbutton[1]) {
      // zoom
   }
   else if (GLUTbutton[2]) {
      // move voxel
      R2Point world_point = R2Point(x, y);
      R2Point grid_point = ConvertWorldToGrid(world_point);

      int window_xdim = (projection_dim + 1) % 3;
      int window_ydim = (projection_dim + 2) % 3;

      if (grid_point.X() < 0) grid_point.SetCoord(RN_X, 0);
      if (grid_point.Y() < 0) grid_point.SetCoord(RN_Y, 0);
      if (grid_point.X() > image_grid->Resolution(window_xdim) - 1) grid_point.SetCoord(RN_X, image_grid->Resolution(window_xdim) - 1);
      if (grid_point.Y() > image_grid->Resolution(window_ydim) - 1) grid_point.SetCoord(RN_Y, image_grid->Resolution(window_ydim) - 1);

      R3Point intersection;
      if (projection_dim == RN_X) intersection = R3Point(selected_slice_index, grid_point.X(), grid_point.Y());
      else if (projection_dim == RN_Y) intersection = R3Point(grid_point.Y(), selected_slice_index, grid_point.X());
      else intersection = R3Point(grid_point.X(), grid_point.Y(), selected_slice_index);

      // update the global selected position
      selected_position = intersection;
      selected_position.Transform(transformation);

      // update the selected voxel
      selected_voxel[RN_X] = (int)(intersection.X() + 0.5);
      selected_voxel[RN_Y] = (int)(intersection.Y() + 0.5);
      selected_voxel[RN_Z] = (int)(intersection.Z() + 0.5);
   }
}



void GLUTMotion3D(int x, int y)
{
   // compute mouse movement
   int dx = x - GLUTmouse[0];
   int dy = y - GLUTmouse[1];

   // world in hand navigation
   R3Point origin = world_box.Centroid();
   if (GLUTbutton[0]) viewer->RotateWorld(1.0, origin, x, y, dx, dy);
   else if (GLUTbutton[1]) viewer->ScaleWorld(1.0, origin, x, y, dx, dy);
   else if (GLUTbutton[2]) {
      if (glutGetModifiers() == GLUT_ACTIVE_SHIFT) {
         viewer->TranslateWorld(1.0, origin, x, y, dx, dy);
      }
      else {
         R3Vector towards = viewer->Camera().Towards();
         RNDimension dim = towards.MaxDimension();
         R3Vector normal = R3xyz_triad[dim];
         R3Plane plane(selected_position, normal);
         R3Ray viewer_ray = viewer->WorldRay(x, y);
         R3Point intersection;
         if (R3Intersects(viewer_ray, plane, &intersection)) selected_position = intersection;

         // confine point by x dimension
         if (selected_position.X() < world_box.XMin())
            selected_position.SetX(world_box.XMin());
         else if (selected_position.X() >= world_box.XMax())
            selected_position.SetX(world_box.XMax());
         // confine point by y dimension
         if (selected_position.Y() < world_box.YMin())
            selected_position.SetY(world_box.YMin());
         else if (selected_position.Y() >= world_box.YMax())
            selected_position.SetY(world_box.YMax());
         // confine point by z dimension
         if (selected_position.Z() < world_box.ZMin())
            selected_position.SetZ(world_box.ZMin());
         else if (selected_position.Z() >= world_box.ZMax())
            selected_position.SetZ(world_box.ZMax());

         // update selected voxel
         R3Point voxel = selected_position;
         voxel.InverseTransform(transformation);
         selected_voxel[RN_X] = (int)(voxel.X() + 0.5);
         selected_voxel[RN_Y] = (int)(voxel.Y() + 0.5);
         selected_voxel[RN_Z] = (int)(voxel.Z() + 0.5);
      }
   }
}



void GLUTMotion(int x, int y)
{
   // invert y coordinate
   y = GLUTwindow_height - y;

   // different motion clicks for projection view
   if (projection) GLUTMotion2D(x, y);
   else GLUTMotion3D(x, y);

   // redisplay if a mouse was down
   if (GLUTbutton[0] || GLUTbutton[1] || GLUTbutton[2]) glutPostRedisplay();

   // remember mouse position
   GLUTmouse[0] = x;
   GLUTmouse[1] = y;
}



void GLUTMouse2D(int button, int state, int x, int y)
{
   if (button == 0) {
      if (state == GLUT_DOWN) {
         // query 
         R2Point world_point = R2Point(x, y);
         selected_slice_position = ConvertGridToWorld(world_point);

         // get window dimensions
         int window_xdim = (projection_dim + 1) % 3;
         int window_ydim = (projection_dim + 2) % 3;

         // get the grid x and y coordinates
         int gridx = (int)(selected_slice_position.X() + 0.5);
         int gridy = (int)(selected_slice_position.Y() + 0.5);

         // cannot select a point outside of the image
         if (gridx < 0 || gridx >= image_grid->Resolution(window_xdim)) selected_slice_position.Reset(RN_UNKNOWN, RN_UNKNOWN);
         else if (gridy < 0 || gridy >= image_grid->Resolution(window_ydim)) selected_slice_position.Reset(RN_UNKNOWN, RN_UNKNOWN);

         glutPostRedisplay();
      }
      else {
         selected_slice_position.Reset(RN_UNKNOWN, RN_UNKNOWN);
         glutPostRedisplay();
      }
   }
}



void GLUTMouse3D(int button, int state, int x, int y)
{
   if (button == 2) {
      if (glutGetModifiers() != GLUT_ACTIVE_SHIFT) {
         if (state == GLUT_DOWN) {
            R3Vector towards = viewer->Camera().Towards();
            RNDimension dim = towards.MaxDimension();
            R3Vector normal = R3xyz_triad[dim];
            R3Plane plane(selected_position, normal);
            R3Ray viewer_ray = viewer->WorldRay(x, y);
            R3Point intersection;
            if (R3Intersects(viewer_ray, plane, &intersection)) selected_position = intersection;

            // confine point by x dimension
            if (selected_position.X() < world_box.XMin())
               selected_position.SetX(world_box.XMin());
            else if (selected_position.X() >= world_box.XMax())
               selected_position.SetX(world_box.XMax());
            // confine point by y dimension
            if (selected_position.Y() < world_box.YMin())
               selected_position.SetY(world_box.YMin());
            else if (selected_position.Y() >= world_box.YMax())
               selected_position.SetY(world_box.YMax());
            // confine point by z dimension
            if (selected_position.Z() < world_box.ZMin())
               selected_position.SetZ(world_box.ZMin());
            else if (selected_position.Z() >= world_box.ZMax())
               selected_position.SetZ(world_box.ZMax());

            // update selected voxel
            R3Point voxel = selected_position;
            voxel.InverseTransform(transformation);
            selected_voxel[RN_X] = (int)(voxel.X() + 0.5);
            selected_voxel[RN_Y] = (int)(voxel.Y() + 0.5);
            selected_voxel[RN_Z] = (int)(voxel.Z() + 0.5);
         }
      }
   }
}



void GLUTMouse(int button, int state, int x, int y)
{
   // invert y coordinate
   y = GLUTwindow_height - y;

   // process mouse button event
   if (projection) GLUTMouse2D(button, state, x, y);
   else GLUTMouse3D(button, state, x, y);

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

   switch (key) {
   case GLUT_KEY_UP: {
      selected_slice_index++;
      if (selected_slice_index > image_grid->Resolution(projection_dim) - 1) selected_slice_index = image_grid->Resolution(projection_dim) - 1;
      UpdateSlices();
      break; }

   case GLUT_KEY_DOWN: {
      selected_slice_index--;
      if (selected_slice_index < 0) selected_slice_index = 0;
      UpdateSlices();
      break; }
   
   case GLUT_KEY_LEFT: {
       segmentation_index--;
       if (segmentation_index < 0) segmentation_index = 0;    
       break; }

    case GLUT_KEY_RIGHT: {
        segmentation_index++;
        if (segmentation_index > (int)(segmentation_grid->Maximum() + 0.5)) segmentation_index = segmentation_grid->Maximum();         
        break; }
   }

   // redraw
   glutPostRedisplay();
}



void GLUTKeyboard2D(unsigned char key, int x, int y)
{
   switch (key) {
   case 'A':
   case 'a': {
      show_projection_affinity = 1;
      show_projection_image = 0;
      show_projection_segmentation = 0;
      show_projection_gold = 0;
      selected_slice = affinity_selected_slice[projection_dim];
      selected_slice_range = affinity_range[projection_dim];
      break; }

   case 'C':
   case 'c': {
       color_type = 1 - color_type;
       break; }

   case 'G':
   case 'g': {
      show_projection_affinity = 0;
      show_projection_image = 0;
      show_projection_segmentation = 0;
      show_projection_gold = 1;
      selected_slice = gold_selected_slice;
      selected_slice_range = gold_range;
      break; }

   case 'I':
   case 'i': {
      show_projection_affinity = 0;
      show_projection_image = 1;
      show_projection_segmentation = 0;
      show_projection_gold = 0;
      selected_slice = image_selected_slice;
      selected_slice_range = image_range;
      break; }

   case 'S':
   case 's': {
      show_projection_affinity = 0;
      show_projection_image = 0;
      show_projection_segmentation = 1;
      show_projection_gold = 0;
      selected_slice = segmentation_selected_slice;
      selected_slice_range = segmentation_range;
      break; }
   }
}



void GLUTKeyboard3D(unsigned char key, int x, int y)
{
   switch (key) {
   case 'A':
   case 'a': {
      show_dimension_affinity = 1;
      show_dimension_image = 0;
      show_dimension_gold = 0;
      show_dimension_segmentation = 0;
      selected_slice = affinity_selected_slice[projection_dim];
      selected_slice_range = affinity_range[projection_dim];
      break; }

   case 'B':
   case 'b': {
      show_bbox = 1 - show_bbox;
      break; }

   case 'G':
   case 'g': {
      show_dimension_affinity = 0;
      show_dimension_image = 0;
      show_dimension_segmentation = 0;
      show_dimension_gold = 1;
      selected_slice = gold_selected_slice;
      selected_slice_range = gold_range;
      break; }

   case 'I':
   case 'i': {
      show_dimension_affinity = 0;
      show_dimension_image = 1;
      show_dimension_segmentation = 0;
      show_dimension_gold = 0;
      selected_slice = image_selected_slice;
      selected_slice_range = image_range;
      break; }

   case 'S':
   case 's': {
      show_dimension_affinity = 0;
      show_dimension_image = 0;
      show_dimension_segmentation = 1;
      show_dimension_gold = 0;
      selected_slice = segmentation_selected_slice;
      selected_slice_range = segmentation_range;
      break; }


   case 'W':
   case 'w': {
      show_slice = 1 - show_slice;
      break; }
   }
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

   // different keys based on projection or normal
   if (projection) {
      GLUTKeyboard2D(key, x, y);
   }
   else {
      GLUTKeyboard3D(key, x, y);
   }


   // keys regardless of projection status
   switch (key) {
   case '1': {
      affinity_dim = RN_X;
      break; }

   case '2': {
      affinity_dim = RN_Y;
      break; }

   case '3': {
      affinity_dim = RN_Z;
      break; }

   case 'X':
   case 'x': {
      projection_dim = RN_X;
      if (selected_slice_index > image_grid->XResolution()) selected_slice_index = image_grid->XResolution() - 1;
      break; }

   case 'Y':
   case 'y': {
      projection_dim = RN_Y;
      if (selected_slice_index > image_grid->YResolution()) selected_slice_index = image_grid->YResolution() - 1;
      break; }

   case 'Z':
   case 'z': {
      projection_dim = RN_Z;
      if (selected_slice_index > image_grid->ZResolution()) selected_slice_index = image_grid->ZResolution() - 1;
      break; }

   case ENTER: {
      projection = 1 - projection;
      break; }

   case ESCAPE: {
      GLUTStop();
      break; }
   }

   // redraw
   glutPostRedisplay();
}



void GLUTInit(int *argc, char **argv)
{
   // open window
   glutInit(argc, argv);
   glutInitWindowPosition(10, 10);
   glutInitWindowSize(GLUTwindow_width, GLUTwindow_height);
   glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH);
   GLUTwindow = glutCreateWindow("Neuron Visualizer");

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

   // initialize font
#if (RN_OS == RN_WINDOWSNT)
   int font = glGenLists(256);
   wglUseFontBitmaps(wglGetCurrentDC(), 0, 256, font);
   glListBase(font);
#endif
}



void GLUTMainLoop(void)
{
   // run main loop -- never returns
   glutMainLoop();
}



////////////////////////////////////////////////////////////////////////
// Viewer functions
////////////////////////////////////////////////////////////////////////

static R3Viewer *CreateViewer(void)
{
   // start statistics
   RNTime start_time;
   start_time.Read();


   if (world_box.IsEmpty()) RNAbort("Error in CreateViewer - box is empty");
   RNLength r = world_box.DiagonalRadius();
   if (r < 0 || RNIsInfinite(r)) RNAbort("Error in CreateViewer - r must be positive finite");

   // set up camera view looking down the z axis
   static R3Vector initial_camera_towards = R3Vector(0.0, 0.0, -1.5);
   static R3Vector initial_camera_up = R3Vector(0.0, 1.0, 0.0);
   R3Point initial_camera_origin = world_box.Centroid() - initial_camera_towards * 2.5 * r;
   R3Camera camera(initial_camera_origin, initial_camera_towards, initial_camera_up, 0.4, 0.4, 0.1 * r, 1000.0 * r);
   R2Viewport viewport(0, 0, GLUTwindow_width, GLUTwindow_height);
   R3Viewer *viewer = new R3Viewer(camera, viewport);

   // print statistics
   if (print_verbose) {
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

static int ParseArgs(int argc, char **argv)
{
   // parse arguments
   argc--; argv++;
   while (argc > 0) {
      if ((*argv)[0] == '-') {
         if (!strcmp(*argv, "-v")) print_verbose = 1;
         else if (!strcmp(*argv, "-debug")) print_debug = 1;
         else if (!strcmp(*argv, "-affinities")) {
             argv++; argc--;
             affinities_filename = *argv;
             argv++; argc--;
             affinities_dataset = *argv;
         }
         else if (!strcmp(*argv, "-gold")) {
             argv++; argc--;
             gold_filename = *argv;
             argv++; argc--;
             gold_dataset = *argv;
         }
         else if (!strcmp(*argv, "-image")) { 
             argv++; argc--;
             image_filename = *argv;
             argv++; argc--;
             image_dataset = *argv;
         }
         else if (!strcmp(*argv, "-segmentation")) {
             argv++; argc--;
             segmentation_filename = *argv;
             argv++; argc--;
             segmentation_dataset = *argv;
         }
         else if (!strcmp(*argv, "-scaling")) {
             argv++; argc--;
             scaling[RN_X] = atof(*argv);
             argv++; argc--;
             scaling[RN_Y] = atof(*argv);
             argv++; argc--; 
             scaling[RN_Z] = atof(*argv); 
         }
         else { fprintf(stderr, "Invalid program argument: %s\n", *argv); return 0; }
      }
      else {
          if (!prefix) { prefix = *argv; }
         else { fprintf(stderr, "Invalid program argument: %s\n", *argv); return 0; }
      }
      argv++; argc--;
   }

   // error if there is no input name
   if (!prefix && (!affinities_filename && !image_filename && !gold_filename && !segmentation_filename)) { fprintf(stderr, "Need to supply a neuron data file\n"); return 0; }

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

   /////////////////////////////////
   //// Read in the voxel files ////
   /////////////////////////////////

   if (!ReadData()) exit(-1);

   // update pointers to selected slice
   UpdateSlices();

   // set world box
   world_box = R3Box(0, 0, 0, image_grid->XResolution() * scaling[RN_X], image_grid->YResolution() * scaling[RN_Y], image_grid->ZResolution() * scaling[RN_Z]);
   selected_position = world_box.Centroid();

   // get the transformation
   transformation = R3Affine(R4Matrix(scaling[RN_X], 0, 0, 0, 0, scaling[RN_Y], 0, 0, 0, 0, scaling[RN_Z], 0, 0, 0, 0, 1));

   // initialize selected voxel
   R3Point voxel = selected_position;
   voxel.InverseTransform(transformation);
   selected_voxel[RN_X] = (int)(voxel.X() + 0.5);
   selected_voxel[RN_Y] = (int)(voxel.Y() + 0.5);
   selected_voxel[RN_Z] = (int)(voxel.Z() + 0.5);
   selected_slice_index = selected_voxel[projection_dim];

   // set projection scale variables
   for (int dim = 0; dim <= 2; ++dim) {
      RNScalar ratio = GLUTwindow_height / (scaling[dim] * image_grid->Resolution(dim));
      if (ratio < projection_scale) projection_scale = ratio;
   }

    // create mapping from segmentation to gold
    if (!MapSegmentationToGold()) { fprintf(stderr, "Failed to allocate memory for segmentation-gold mapping.\n"); exit(-1); }

    // create mapping from index to ids
    if (!MapIndexToSegmentation()) { fprintf(stderr, "Failed to allocate memory for index-segmentation mapping.\n"); exit(-1); }

    // create segmentation vectors
    if (!PopulateSegmentationVectors()) { fprintf(stderr, "Failed to allocate memory for segmentation vectors.\n"); exit(-1); }

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
