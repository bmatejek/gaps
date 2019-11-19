// Source file for connectome visualization



// include files

#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <algorithm>
#include "R3Graphics/R3Graphics.h"
#include "RNDataStructures/RNDataStructures.h"
#include "fglut/fglut.h"



// GLUT defines

#define ENTER 13
#define ESCAPE 27
#define SPACEBAR 32

static const int OR_Z = 0;
static const int OR_Y = 1;
static const int OR_X = 2;



// program arguments

// I/O flags
static int print_debug = 0;
static int print_verbose = 0;
static const char* prefix = NULL;



// meta data variables

static double resolution[3] = { -1, -1, -1 };
static long volume_size[3] = { -1, -1, -1 };
static long block_size[3] = { -1, -1, -1 };
static char blocks_directory[4096];
static char synapses_directory[4096];
static char somae_directory[4096];
static char skeleton_directory[4096];



// current block

static int center_block_z_index = 0;
static int center_block_y_index = 0;
static int center_block_x_index = 0;



// program variables

static R3Affine transformation = R3null_affine;
static R3Viewer *viewer = NULL;
static R3Box world_box = R3null_box;



// display variables

static double background_color[3] = { 0.0, 0.0, 0.0 };
static int show_bbox = 1;
static const int ncolor_opts = 24;
static int color_cycle = 0;
static RNScalar downsample_rate = 0.5;
static double synapse_point_size = 100;
static double skeleton_point_size = 3;
static int show_axes = 1;



// GLUT variables

static int GLUTwindow = 0;
static int GLUTwindow_height = 720;
static int GLUTwindow_width = 1280;
static int GLUTmouse[2] = { 0, 0 };
static int GLUTbutton[3] = { 0, 0, 0 };
static int GLUTmodifiers = 0;



// variables for the neighboring cubes

static const short ncubes = 8;
static long offsets[ncubes][3] = {
    { 0, 0, 0 },
    { 0, 0, 1 },
    { 0, 1, 0 },
    { 0, 1, 1 },
    { 1, 0, 0 },
    { 1, 0, 1 },
    { 1, 1, 0 },
    { 1, 1, 1 }
};
static R3Box block_boxes[ncubes];


// random access variables

static long label = 0;
static long label_index = -1;
static std::vector<long> labels;
static std::unordered_map<long, std::vector<long> > *surfaces;
static std::unordered_map<long, std::vector<long> > *synapses;
static std::unordered_map<long, std::vector<long> > *skeletons;




////////////////////////////////////////////////////////////////////////
// Helper functions
////////////////////////////////////////////////////////////////////////

static void IndexToIndices(long index, long& ix, long& iy, long& iz)
{
  // Set indices of grid value at index
  iz = index / (block_size[OR_X] * block_size[OR_Y]);
  iy = (index - iz * block_size[OR_X] * block_size[OR_Y]) / block_size[OR_X];
  ix = index % block_size[OR_X];
}



static long IndicesToIndex(long ix, long iy, long iz)
{
    return iz * block_size[OR_X] * block_size[OR_Y] + iy * block_size[OR_X] + ix;
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
// Input/output functions
////////////////////////////////////////////////////////////////////////

static int ReadMetaData(void)
{
    // read the meta data
    char meta_filename[4096];
    snprintf(meta_filename, 4096, "meta/%s.meta", prefix);

    FILE *fp = fopen(meta_filename, "r");
    if (!fp) { fprintf(stderr, "Failed to read %s\n", meta_filename); return 0; }

    // dummy variable to ignore comments
    char comment[4096];

    // read in resolution information
    if (!fgets(comment, 4096, fp)) return 0;
    if (fscanf(fp, "%lfx%lfx%lf\n", &(resolution[OR_X]), &(resolution[OR_Y]), &(resolution[OR_Z])) != 3) return 0;

    // read in the volume size
    if (!fgets(comment, 4096, fp)) return 0;
    if (fscanf(fp, "%ldx%ldx%ld\n", &(volume_size[OR_X]), &(volume_size[OR_Y]), &(volume_size[OR_Z])) != 3) return 0;

    // read in the block size
    if (!fgets(comment, 4096, fp)) return 0;
    if (fscanf(fp, "%ldx%ldx%ld\n", &(block_size[OR_X]), &(block_size[OR_Y]), &(block_size[OR_Z])) != 3) return 0;

    // read in the blocks directory
    if (!fgets(comment, 4096, fp)) return 0;
    if (fscanf(fp, "%s\n", blocks_directory) != 1) return 0;

    // read in the synapses directory
    if (!fgets(comment, 4096, fp)) return 0;
    if (fscanf(fp, "%s\n", synapses_directory) != 1) return 0;

    // read in the somae directory
    if (!fgets(comment, 4096, fp)) return 0;
    if (fscanf(fp, "%s\n", somae_directory) != 1) return 0;

    // read in the skeleton directory
    if (!fgets(comment, 4096, fp)) return 0;
    if (fscanf(fp, "%s\n", skeleton_directory) != 1) return 0;

    // close the file
    fclose(fp);

    return 1;
}




static int ReadBlockData(long z_block_index, long y_block_index, long x_block_index, long linear_block_index)
{
    // keep statistics
    RNTime start_time;
    start_time.Read();

    // clear the previous segment information
    std::vector<long> block_labels = std::vector<long>();

    // get the segment file for this block
    char segment_filename[4096];
    snprintf(segment_filename, 4096, "%s/%s-input_labels-%04ldz-%04ldy-%04ldx.h5", blocks_directory, prefix, z_block_index, y_block_index, x_block_index);

    R3Grid **segmentations = RNReadH5File(segment_filename, "main");
    R3Grid *segmentation = segmentations[0];
    delete[] segmentations;

    int z_grid_size = segmentation->ZResolution();
    int y_grid_size = segmentation->YResolution();
    int x_grid_size = segmentation->XResolution();

    // go through all of the voxels and add surface points to surfaces mapping
    for (int iz = 1; iz < z_grid_size - 1; ++iz) {
        for (int iy = 1; iy < y_grid_size - 1; ++iy) {
            for (int ix = 1; ix < x_grid_size - 1; ++ix) {
                // get the linear index 
                long iv = IndicesToIndex(ix, iy, iz);

                // get the label for this location
                long label = (long) (segmentation->GridValue(ix, iy, iz) + 0.5);

                // skip extracellular material
                if (not label) continue;

                // has this label been seen so far
                if (surfaces[linear_block_index].find(label) == surfaces[linear_block_index].end()) {
                    surfaces[linear_block_index][label] = std::vector<long>();
                    block_labels.push_back(label);
                }

                // go through all six adjacent neighbors
                if ((long)(segmentation->GridValue(ix - 1, iy, iz) + 0.5) != label ||
                    (long)(segmentation->GridValue(ix + 1, iy, iz) + 0.5) != label ||
                    (long)(segmentation->GridValue(ix, iy - 1, iz) + 0.5) != label ||
                    (long)(segmentation->GridValue(ix, iy + 1, iz) + 0.5) != label ||
                    (long)(segmentation->GridValue(ix, iy, iz - 1) + 0.5) != label ||
                    (long)(segmentation->GridValue(ix, iy, iz + 1) + 0.5) != label)
                {
                    surfaces[linear_block_index][label].push_back(iv);
                }
            }
        }
    }

    delete segmentation;

    // sort the labels
    std::sort(block_labels.begin(), block_labels.end());

    printf("Processed file %s in %0.2f seconds.\n", segment_filename, start_time.Elapsed());

    // start statistics
    start_time.Read();

    // read the synapses
    char synapse_filename[4096];
    snprintf(synapse_filename, 4096, "%s/%s-synapses-%04ldz-%04ldy-%04ldx.pts", synapses_directory, prefix, z_block_index, y_block_index, x_block_index);

    FILE *fp = fopen(synapse_filename, "rb");
    if (!fp) { fprintf(stderr, "Failed to read %s.\n", synapse_filename); return 0; }

    long z_input_volume_size, y_input_volume_size, x_input_volume_size;
    long z_input_block_size, y_input_block_size, x_input_block_size;

    if (fread(&z_input_volume_size, sizeof(long), 1, fp) != 1) { fprintf(stderr, "Failed to read %s.\n", synapse_filename); return 0; }
    if (fread(&y_input_volume_size, sizeof(long), 1, fp) != 1) { fprintf(stderr, "Failed to read %s.\n", synapse_filename); return 0; }
    if (fread(&x_input_volume_size, sizeof(long), 1, fp) != 1) { fprintf(stderr, "Failed to read %s.\n", synapse_filename); return 0; }
    
    if (z_input_volume_size != volume_size[OR_Z]) { fprintf(stderr, "Failed to read %s.\n", synapse_filename); return 0; }
    if (y_input_volume_size != volume_size[OR_Y]) { fprintf(stderr, "Failed to read %s.\n", synapse_filename); return 0; }
    if (x_input_volume_size != volume_size[OR_X]) { fprintf(stderr, "Failed to read %s.\n", synapse_filename); return 0; }

    if (fread(&z_input_block_size, sizeof(long), 1, fp) != 1) { fprintf(stderr, "Failed to read %s.\n", synapse_filename); return 0; }
    if (fread(&y_input_block_size, sizeof(long), 1, fp) != 1) { fprintf(stderr, "Failed to read %s.\n", synapse_filename); return 0; }
    if (fread(&x_input_block_size, sizeof(long), 1, fp) != 1) { fprintf(stderr, "Failed to read %s.\n", synapse_filename); return 0; }

    if (z_input_block_size != block_size[OR_Z]) { fprintf(stderr, "Failed to read %s.\n", synapse_filename); return 0; }
    if (y_input_block_size != block_size[OR_Y]) { fprintf(stderr, "Failed to read %s.\n", synapse_filename); return 0; }
    if (x_input_block_size != block_size[OR_X]) { fprintf(stderr, "Failed to read %s.\n", synapse_filename); return 0; }

    long nneurons;
    if (fread(&nneurons, sizeof(long), 1, fp) != 1) { fprintf(stderr, "Failed to read %s.\n", synapse_filename); return 0; }

    for (long iv = 0; iv < nneurons; ++iv) {
        // get the label and number of synapses
        long label;
        long nsynapses;
        
        if (fread(&label, sizeof(long), 1, fp) != 1) { fprintf(stderr, "Failed to read %s.\n", synapse_filename); return 0; }
        if (fread(&nsynapses, sizeof(long), 1, fp) != 1) { fprintf(stderr, "Failed to read %s.\n", synapse_filename); return 0; }
    
        synapses[linear_block_index][label] = std::vector<long>();

        // ignore the global coordinates
        for (long is = 0; is < nsynapses; ++is) {
            long dummy_index;
            if (fread(&dummy_index, sizeof(long), 1, fp) != 1)  { fprintf(stderr, "Failed to read %s.\n", synapse_filename); return 0; }
        }

        // add the local coordinates
        for (long is = 0; is < nsynapses; ++is) {
            long linear_index;
            if (fread(&linear_index, sizeof(long), 1, fp) != 1)  { fprintf(stderr, "Failed to read %s.\n", synapse_filename); return 0; }
            synapses[linear_block_index][label].push_back(linear_index);
        }

    }

    // close file
    fclose(fp);

    printf("Processed file %s in %0.2f seconds.\n", synapse_filename, start_time.Elapsed());

    // start statistics
    start_time.Read();

    for (unsigned long iv = 0; iv < block_labels.size(); ++iv) {
        long neuron_id = block_labels[iv];

        // read the synapses
        char skeleton_filename[4096];
        snprintf(skeleton_filename, 4096, "%s/%s-skeleton-%04ldz-%04ldy-%04ldx-ID-%012ld.pts", skeleton_directory, prefix, z_block_index, y_block_index, x_block_index, neuron_id);

        FILE *fp = fopen(skeleton_filename, "rb");
        if (!fp) { continue; }

        long z_input_volume_size, y_input_volume_size, x_input_volume_size;
        long z_input_block_size, y_input_block_size, x_input_block_size;

        if (fread(&z_input_volume_size, sizeof(long), 1, fp) != 1) { fprintf(stderr, "Failed to read %s.\n", skeleton_filename); return 0; }
        if (fread(&y_input_volume_size, sizeof(long), 1, fp) != 1) { fprintf(stderr, "Failed to read %s.\n", skeleton_filename); return 0; }
        if (fread(&x_input_volume_size, sizeof(long), 1, fp) != 1) { fprintf(stderr, "Failed to read %s.\n", skeleton_filename); return 0; }
        
        if (z_input_volume_size != volume_size[OR_Z]) { fprintf(stderr, "Failed to read %s.\n", skeleton_filename); return 0; }
        if (y_input_volume_size != volume_size[OR_Y]) { fprintf(stderr, "Failed to read %s.\n", skeleton_filename); return 0; }
        if (x_input_volume_size != volume_size[OR_X]) { fprintf(stderr, "Failed to read %s.\n", skeleton_filename); return 0; }

        if (fread(&z_input_block_size, sizeof(long), 1, fp) != 1) { fprintf(stderr, "Failed to read %s.\n", skeleton_filename); return 0; }
        if (fread(&y_input_block_size, sizeof(long), 1, fp) != 1) { fprintf(stderr, "Failed to read %s.\n", skeleton_filename); return 0; }
        if (fread(&x_input_block_size, sizeof(long), 1, fp) != 1) { fprintf(stderr, "Failed to read %s.\n", skeleton_filename); return 0; }

        if (z_input_block_size != block_size[OR_Z]) { fprintf(stderr, "Failed to read %s.\n", skeleton_filename); return 0; }
        if (y_input_block_size != block_size[OR_Y]) { fprintf(stderr, "Failed to read %s.\n", skeleton_filename); return 0; }
        if (x_input_block_size != block_size[OR_X]) { fprintf(stderr, "Failed to read %s.\n", skeleton_filename); return 0; }

        // get the label and number of synapses
        long input_neuron_id;
        long nskeleton_points;
        
        if (fread(&input_neuron_id, sizeof(long), 1, fp) != 1) { fprintf(stderr, "Failed to read %s.\n", skeleton_filename); return 0; }
        if (fread(&nskeleton_points, sizeof(long), 1, fp) != 1) { fprintf(stderr, "Failed to read %s.\n", skeleton_filename); return 0; }
        if (input_neuron_id != neuron_id) { fprintf(stderr, "Failed to read %s.\n", skeleton_filename); return 0; }

        skeletons[linear_block_index][neuron_id] = std::vector<long>();
        // ignore the global coordinates
        for (long is = 0; is < nskeleton_points; ++is) {
            long dummy_index;
            if (fread(&dummy_index, sizeof(long), 1, fp) != 1)  { fprintf(stderr, "Failed to read %s.\n", skeleton_filename); return 0; }
        }

        // add the local coordinates
        for (long is = 0; is < nskeleton_points; ++is) {
            long linear_index;
            if (fread(&linear_index, sizeof(long), 1, fp) != 1)  { fprintf(stderr, "Failed to read %s.\n", skeleton_filename); return 0; }
            skeletons[linear_block_index][neuron_id].push_back(linear_index);
        }

        // close file
        fclose(fp);
    }

    printf("Processed skeleton files in %0.2f seconds.\n", start_time.Elapsed());

    // return success
    return 1;
}




////////////////////////////////////////////////////////////////////////
// Drawing utility functions
////////////////////////////////////////////////////////////////////////

static void DrawSurface(long linear_block_index)
{
    long zoffset = offsets[linear_block_index][OR_Z] * block_size[OR_Z];
    long yoffset = offsets[linear_block_index][OR_Y] * block_size[OR_Y];
    long xoffset = offsets[linear_block_index][OR_X] * block_size[OR_X];

    transformation.Push();
    glPointSize(1.0);
    glBegin(GL_POINTS);
    for (unsigned long iv = 0; iv < surfaces[linear_block_index][label].size(); ++iv) {
        // faster rendering with downsampling
        if (RNRandomScalar() > downsample_rate) continue;

        // get the coordinates from the linear index
        long ix, iy, iz;
        IndexToIndices(surfaces[linear_block_index][label][iv], ix, iy, iz);
        glVertex3f(ix + xoffset, iy + yoffset, iz + zoffset);
    }
    glEnd();
    transformation.Pop();
}



static void DrawSynapse(long linear_block_index)
{
    long zoffset = offsets[linear_block_index][OR_Z] * block_size[OR_Z];
    long yoffset = offsets[linear_block_index][OR_Y] * block_size[OR_Y];
    long xoffset = offsets[linear_block_index][OR_X] * block_size[OR_X];

    for (unsigned int iv = 0; iv < synapses[linear_block_index][label].size(); ++iv) {
        // get the coordinates from the linear index
        long ix, iy, iz;
        IndexToIndices(synapses[linear_block_index][label][iv], ix, iy, iz);

        // add the offsets
        ix += xoffset;
        iy += yoffset; 
        iz += zoffset;

        ix *= resolution[OR_X];
        iy *= resolution[OR_Y];
        iz *= resolution[OR_Z];

        R3Sphere(R3Point(ix, iy, iz), synapse_point_size).Draw();
    }
}



static void DrawSkeleton(long linear_block_index)
{
    long zoffset = offsets[linear_block_index][OR_Z] * block_size[OR_Z];
    long yoffset = offsets[linear_block_index][OR_Y] * block_size[OR_Y];
    long xoffset = offsets[linear_block_index][OR_X] * block_size[OR_X];

    transformation.Push();
    glPointSize(skeleton_point_size);
    glBegin(GL_POINTS);
    for (unsigned int iv = 0; iv < skeletons[linear_block_index][label].size(); ++iv) {
        // get the coordinates from the linear index
        long ix, iy, iz;
        IndexToIndices(skeletons[linear_block_index][label][iv], ix, iy, iz);
        glVertex3f(ix + xoffset, iy + yoffset, iz + zoffset);
    }
    glEnd();

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
    if (show_bbox) {
        for (long iu = 0; iu < ncubes; ++iu) {
            RNLoadRgb(RNRgb(not background_color[0], not background_color[1], not background_color[2]));
            block_boxes[iu].Outline();
        }
    }

    // draw segmentation
    for (long iu = 0; iu < ncubes; ++iu) {
        RNLoadRgb(Color(label));
        DrawSurface(iu);
        RNLoadRgb(RNRgb(not background_color[0], not background_color[1], not background_color[2]));
        DrawSynapse(iu);
        DrawSkeleton(iu);
    }

    if (show_axes) {
        R3Point origin = world_box.Centroid();
        R3Vector xaxis = world_box.XRadius() * R3Vector(1, 0, 0) / 2;
        R3Point xpoint = origin + xaxis;
        R3Vector yaxis = world_box.YRadius() * R3Vector(0, 1, 0) / 2;
        R3Point ypoint = origin + yaxis;
        R3Vector zaxis = world_box.ZRadius() * R3Vector(0, 0, 1) / 2;
        R3Point zpoint = origin + zaxis;

        // draw the x axis
        RNLoadRgb(RNred_rgb);
        glBegin(GL_LINES);
        glVertex3f(origin.X(), origin.Y(), origin.Z());
        glVertex3f(xpoint.X(), xpoint.Y(), xpoint.Z());
        glEnd();

        // draw the y axis
        RNLoadRgb(RNgreen_rgb);
        glBegin(GL_LINES);
        glVertex3f(origin.X(), origin.Y(), origin.Z());
        glVertex3f(ypoint.X(), ypoint.Y(), ypoint.Z());
        glEnd();


        // draw the
        RNLoadRgb(RNblue_rgb);
        glBegin(GL_LINES);
        glVertex3f(origin.X(), origin.Y(), origin.Z());
        glVertex3f(zpoint.X(), zpoint.Y(), zpoint.Z());
        glEnd();
    }

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
            --label_index;
            if (label_index < 0) label_index = 0;
            label = labels[label_index];

            printf("Neuron ID: %ld\n", label);                

            break;
        }

        case GLUT_KEY_RIGHT: {
            ++label_index;
            if (label_index > (long) labels.size() - 1) label_index = labels.size() - 1;
            label = labels[label_index];

            printf("Neuron ID %ld\n", label);

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
            show_axes = 1 - show_axes;
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

        case 'D': 
        case 'd': {
            if (not color_cycle) color_cycle = ncolor_opts - 1;
            else color_cycle = (color_cycle - 1) % ncolor_opts;
            break;
        }

        case 'E': 
        case 'e': {
            synapse_point_size -= 5;
            break;
        }

        case 'F':
        case 'f': {
            synapse_point_size += 5;
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

        case 'X':
        case 'x': {
            downsample_rate += 0.01;
            printf("Downsample Rate: %lf\n", downsample_rate);
            break;
        }

        case 'Z':
        case 'z': {
            downsample_rate -= 0.01;
            printf("Downsample Rate: %lf\n", downsample_rate);
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
    GLUTwindow = glutCreateWindow("Wiring Diagram Visualizer");

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
            else if (!strcmp(*argv, "-block")) {
                argv++; argc--;
                center_block_z_index = atoi(*argv);
                argv++; argc--;
                center_block_y_index = atoi(*argv);
                argv++; argc--;
                center_block_x_index = atoi(*argv);
            }
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

    // read the meta data for this prefix
    if (!ReadMetaData()) exit(-1);

    long nzblocks = (long) ceil(volume_size[OR_Z] / (double)block_size[OR_Z]);
    long nyblocks = (long) ceil(volume_size[OR_Y] / (double)block_size[OR_Y]);
    long nxblocks = (long) ceil(volume_size[OR_X] / (double)block_size[OR_X]);

    if (center_block_x_index > nxblocks - 2 || center_block_y_index > nyblocks - 2 || center_block_z_index > nzblocks - 2) {
        printf("Choose a block whose +1 neighbors are in the volume\n");
        return 0;
    }

    if (print_verbose) {
        printf("%s Dataset\n", prefix);
        printf("  Resolution: %0.1lf x %0.1lf x %0.1lf\n", resolution[OR_X], resolution[OR_Y], resolution[OR_Z]);
        printf("  Volume: %ld x %ld x %ld\n", volume_size[OR_X], volume_size[OR_Y], volume_size[OR_Z]);
        printf("  Block Sizes: %ld x %ld x %ld\n", block_size[OR_X], block_size[OR_Y], block_size[OR_Z]);
    }


    std::unordered_set<long> neuron_ids = std::unordered_set<long>();
   
    surfaces = new std::unordered_map<long, std::vector<long> >[ncubes];
    synapses = new std::unordered_map<long, std::vector<long> >[ncubes];
    skeletons = new std::unordered_map<long, std::vector<long> >[ncubes];

    // read the block data
    for (long iu = 0; iu < ncubes; ++iu) {
        // initialize the maps
        
        surfaces[iu] = std::unordered_map<long, std::vector<long> >();
        synapses[iu] = std::unordered_map<long, std::vector<long> >();
        skeletons[iu] = std::unordered_map<long, std::vector<long> >();
        
        long block_z_index = center_block_z_index + offsets[iu][OR_Z];
        long block_y_index = center_block_y_index + offsets[iu][OR_Y];
        long block_x_index = center_block_x_index + offsets[iu][OR_X];
        
        if (!ReadBlockData(block_z_index, block_y_index, block_x_index, iu)) exit(-1);
        
        for (std::unordered_map<long, std::vector<long> >::iterator it = surfaces[iu].begin(); it != surfaces[iu].end(); ++it) {
            neuron_ids.insert(it->first);
        }
        
        long zoffset = offsets[iu][OR_Z] * block_size[OR_Z];
        long yoffset = offsets[iu][OR_Y] * block_size[OR_Y];
        long xoffset = offsets[iu][OR_X] * block_size[OR_X];
        
        block_boxes[iu] = R3Box(resolution[OR_X] * xoffset, resolution[OR_Y] * yoffset, resolution[OR_Z] * zoffset, resolution[OR_X] * (xoffset + block_size[OR_X]), resolution[OR_Y] * (yoffset + block_size[OR_Y]), resolution[OR_Z] * (zoffset + block_size[OR_Z]));
    }

    // create a labels vector and add neurons to it
    labels = std::vector<long>();

    for (std::unordered_set<long>::iterator it = neuron_ids.begin(); it != neuron_ids.end(); ++it) {
        labels.push_back(*it);
    }
   
    sort(labels.begin(), labels.end());
    label_index = 0;
    label = labels[label_index];

    // set world box
    world_box = R3null_box;
    for (long iu = 0; iu < ncubes; ++iu) {
        world_box.Union(block_boxes[iu]);
    }

    // get the transformation
    transformation = R3Affine(R4Matrix(resolution[OR_X], 0, 0, 0, 0, resolution[OR_Y], 0, 0, 0, 0, resolution[OR_Z], 0, 0, 0, 0, 1));

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