// Source file for connectome visualization



// include files

#include <vector>
#include <unordered_map>
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



// current block

static int z_block_index = 0;
static int y_block_index = 0;
static int x_block_index = 0;



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



// GLUT variables

static int GLUTwindow = 0;
static int GLUTwindow_height = 720;
static int GLUTwindow_width = 1280;
static int GLUTmouse[2] = { 0, 0 };
static int GLUTbutton[3] = { 0, 0, 0 };
static int GLUTmodifiers = 0;



// random access variables

static long label = 0;
static long label_index = -1;
static std::vector<long> labels;
static std::unordered_map<long, std::vector<long> > surfaces;
static std::unordered_map<long, std::vector<long> > synapses;


/*static std::vector<long> synapses;
static std::vector<long> connectome;
static std::vector<long> somae;*/


// display variables

// static int segmentation_index = 1;
// static int skeleton_method = 0;
// static double skeleton_point_size = 3;



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

    // read in the block size
    if (!fgets(comment, 4096, fp)) return 0;
    if (fscanf(fp, "%s\n", blocks_directory) != 1) return 0;

    // read in the block size
    if (!fgets(comment, 4096, fp)) return 0;
    if (fscanf(fp, "%s\n", synapses_directory) != 1) return 0;

    // close the file
    fclose(fp);

    return 1;
}




static int ReadBlockData(void)
{
    // keep statistics
    RNTime start_time;
    start_time.Read();

    // clear the previous segment information
    labels = std::vector<long>();
    surfaces = std::unordered_map<long, std::vector<long> >();
    synapses = std::unordered_map<long, std::vector<long> >();

    // get the segment file for this block
    char segment_filename[4096];
    snprintf(segment_filename, 4096, "%s/%s-%04dz-%04dy-%04dx.h5", blocks_directory, prefix, z_block_index, y_block_index, x_block_index);

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
                if (surfaces.find(label) == surfaces.end()) {
                    surfaces[label] = std::vector<long>();
                    labels.push_back(label);
                }

                // go through all six adjacent neighbors
                if ((long)(segmentation->GridValue(ix - 1, iy, iz) + 0.5) != label ||
                    (long)(segmentation->GridValue(ix + 1, iy, iz) + 0.5) != label ||
                    (long)(segmentation->GridValue(ix, iy - 1, iz) + 0.5) != label ||
                    (long)(segmentation->GridValue(ix, iy + 1, iz) + 0.5) != label ||
                    (long)(segmentation->GridValue(ix, iy, iz - 1) + 0.5) != label ||
                    (long)(segmentation->GridValue(ix, iy, iz + 1) + 0.5) != label)
                {
                    surfaces[label].push_back(iv);
                }
            }
        }
    }

    delete segmentation;

    // sort the labels
    std::sort(labels.begin(), labels.end());

    printf("Processed file %s in %0.2f seconds.\n", segment_filename, start_time.Elapsed());

    // start statistics
    start_time.Read();

    // read the synapses
    char synapse_filename[4096];
    snprintf(synapse_filename, 4096, "%s/%s-%04dz-%04dy-%04dx.pts", synapses_directory, prefix, z_block_index, y_block_index, x_block_index);

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
    
        synapses[label] = std::vector<long>();

        // ignore the global coordinates
        for (long is = 0; is < nsynapses; ++is) {
            long dummy_index;
            if (fread(&dummy_index, sizeof(long), 1, fp) != 1)  { fprintf(stderr, "Failed to read %s.\n", synapse_filename); return 0; }
        }

        // add the local coordinates
        for (long is = 0; is < nsynapses; ++is) {
            long linear_index;
            if (fread(&linear_index, sizeof(long), 1, fp) != 1)  { fprintf(stderr, "Failed to read %s.\n", synapse_filename); return 0; }
            synapses[label].push_back(linear_index);
        }

    }

    // close file
    fclose(fp);

    printf("Processed file %s in %0.2f seconds.\n", synapse_filename, start_time.Elapsed());

    /* TODO would like to keep this when switching blocks */
    label_index = 0;
    label = labels[label_index];

    // return success
    return 1;
}




/*



static int ReadSynapseData(void)
{
    synapses = std::vector<long>();

    char synapse_filename[4096];
    sprintf(synapse_filename, "synapses/%s/%06d.pts", prefix, segmentation_index);

    FILE *fp = fopen(synapse_filename, "rb");
    if (!fp) { fprintf(stderr, "Failed to read %s.\n", synapse_filename); return 0; }

    long nsynapse_points;
    if (fread(&(grid_size[OR_Z]), sizeof(long), 1, fp) != 1) { fprintf(stderr, "Failed to read %s.\n", synapse_filename); return 0; }
    if (fread(&(grid_size[OR_Y]), sizeof(long), 1, fp) != 1) { fprintf(stderr, "Failed to read %s.\n", synapse_filename); return 0; }
    if (fread(&(grid_size[OR_X]), sizeof(long), 1, fp) != 1) { fprintf(stderr, "Failed to read %s.\n", synapse_filename); return 0; }
    if (fread(&nsynapse_points, sizeof(long), 1, fp) != 1) { fprintf(stderr, "Failed to read %s.\n", synapse_filename); return 0; }

    for (long iv = 0; iv < nsynapse_points; ++iv) {
        long voxel_index;
        if (voxel_index < 0) {
            voxel_index = -1 * voxel_index;
        }
        if (fread(&voxel_index, sizeof(long), 1, fp) != 1) { fprintf(stderr, "Failed to read %s.\n", synapse_filename); return 0; }
        synapses.push_back(voxel_index);
    }

    // return success 
    return 1;
}



static int ReadConnectomeData(void)
{
    connectome = std::vector<long>();

    char connectome_filename[4096];
    if (skeleton_method == 0) sprintf(connectome_filename, "connectomes/%s/%06d.pts", prefix, segmentation_index);
    else if (skeleton_method == 1) sprintf(connectome_filename, "skeletons/%s/%06d.pts", prefix, segmentation_index);
    else if (skeleton_method == 2) sprintf(connectome_filename, "baselines/topological-thinnings/%s/%06d.pts", prefix, segmentation_index);
    else if (skeleton_method == 3) sprintf(connectome_filename, "baselines/teasers/%s/%06d.pts", prefix, segmentation_index);
    else return 0;

    FILE *fp = fopen(connectome_filename, "rb"); 
    if (!fp) { fprintf(stderr, "Failed to read %s.\n", connectome_filename); return 0; }

    long nconnectome_points;
    if (fread(&(grid_size[OR_Z]), sizeof(long), 1, fp) != 1) { fprintf(stderr, "Failed to read %s.\n", connectome_filename); return 0; }
    if (fread(&(grid_size[OR_Y]), sizeof(long), 1, fp) != 1) { fprintf(stderr, "Failed to read %s.\n", connectome_filename); return 0; }
    if (fread(&(grid_size[OR_X]), sizeof(long), 1, fp) != 1) { fprintf(stderr, "Failed to read %s.\n", connectome_filename); return 0; }
    if (fread(&nconnectome_points, sizeof(long), 1, fp) != 1) { fprintf(stderr, "Failed to read %s.\n", connectome_filename); return 0; }
    printf("%ld\n", nconnectome_points);
    for (long iv = 0; iv < nconnectome_points; ++iv) {
        long voxel_index;
        if (fread(&voxel_index, sizeof(long), 1, fp) != 1)  { fprintf(stderr, "Failed to read %s.\n", connectome_filename); return 0; }
        if (voxel_index < 0) voxel_index = -1 * voxel_index;
        connectome.push_back(voxel_index);
    }

    // return success
    return 1;
}


static int ReadSomaeData(void)
{
    somae = std::vector<long>();

    char somae_filename[4096];
    sprintf(somae_filename, "volumetric_somae/surfaces/%s/%06d.pts", prefix, segmentation_index);
    
    FILE *fp = fopen(somae_filename, "rb"); 
    if (!fp) { fprintf(stderr, "Failed to read %s.\n", somae_filename); return 0; }

    long nsomae_points;
    if (fread(&(grid_size[OR_Z]), sizeof(long), 1, fp) != 1) { fprintf(stderr, "Failed to read %s.\n", somae_filename); return 0; }
    if (fread(&(grid_size[OR_Y]), sizeof(long), 1, fp) != 1) { fprintf(stderr, "Failed to read %s.\n", somae_filename); return 0; }
    if (fread(&(grid_size[OR_X]), sizeof(long), 1, fp) != 1) { fprintf(stderr, "Failed to read %s.\n", somae_filename); return 0; }
    if (fread(&nsomae_points, sizeof(long), 1, fp) != 1) { fprintf(stderr, "Failed to read %s.\n", somae_filename); return 0; }

    for (long iv = 0; iv < nsomae_points; ++iv) {
        long voxel_index;
        if (fread(&voxel_index, sizeof(long), 1, fp) != 1)  { fprintf(stderr, "Failed to read %s.\n", somae_filename); return 0; }
        if (voxel_index < 0) voxel_index = -1 * voxel_index;
        somae.push_back(voxel_index);
    }

    // return success
    return 1;
}*/





////////////////////////////////////////////////////////////////////////
// Drawing utility functions
////////////////////////////////////////////////////////////////////////

static void DrawSurface(void)
{
    transformation.Push();
    glPointSize(1.0);
    glBegin(GL_POINTS);
    for (unsigned long iv = 0; iv < surfaces[label].size(); ++iv) {
        // faster rendering with downsampling
        if (RNRandomScalar() > downsample_rate) continue;

        // get the coordinates from the linear index
        long ix, iy, iz;
        IndexToIndices(surfaces[label][iv], ix, iy, iz);
        glVertex3f(ix, iy, iz);
    }
    glEnd();
    transformation.Pop();
}



static void DrawSynapse(void)
{
    for (unsigned int iv = 0; iv < synapses[label].size(); ++iv) {
        // get the coordinates from the linear index
        long ix, iy, iz;
        IndexToIndices(synapses[label][iv], ix, iy, iz);

        ix *= resolution[OR_X];
        iy *= resolution[OR_Y];
        iz *= resolution[OR_Z];

        R3Sphere(R3Point(ix, iy, iz), synapse_point_size).Draw();
    }
}



/*static void DrawConnectome(void)
{
    transformation.Push();
    glPointSize(skeleton_point_size);
    glBegin(GL_POINTS);
    for (unsigned int iv = 0; iv < connectome.size(); ++iv) {
        // get the coordinates from the linear index
        long ix, iy, iz;
        IndexToIndices(connectome[iv], ix, iy, iz);
        glVertex3f(ix, iy, iz);
    }
    glEnd();

    transformation.Pop();
}


static void DrawSomae(void)
{
    transformation.Push();
    glPointSize(1.0);
    glBegin(GL_POINTS);
    for (unsigned int iv = 0; iv < somae.size(); ++iv) {
        // get the coordinates from the linear index
        long ix, iy, iz;
        IndexToIndices(somae[iv], ix, iy, iz);
        glVertex3f(ix, iy, iz);
    }
    glEnd();

    transformation.Pop();
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
        RNLoadRgb(RNRgb(not background_color[0], not background_color[1], not background_color[2]));
        world_box.Outline();
    }

    // draw segmentation
    RNLoadRgb(Color(label));
    DrawSurface();
    RNLoadRgb(RNRgb(not background_color[0], not background_color[1], not background_color[2]));
    DrawSynapse();
    //DrawConnectome();
    //DrawSomae();

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

        /*case 'E': 
        case 'e': {
            synapse_point_size -= 5;
            break;
        }

        case 'F':
        case 'f': {
            synapse_point_size += 5;
            break;
        }

        case 'O':
        case 'o': {
            skeleton_method = 0;
            ReadConnectomeData();
            break;
        }

        case 'I':
        case 'i': {
            skeleton_method = 2;
            ReadConnectomeData();
            break;
        }

        case 'N':
        case 'n': {
            skeleton_method = 4;
            ReadConnectomeData();
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
            skeleton_method = 1;
            ReadConnectomeData();
            break;
        }

        case 'T': 
        case 't': {
            skeleton_method = 3;
            ReadConnectomeData();
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
        }*/

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

    if (print_verbose) {
        printf("%s Dataset\n", prefix);
        printf("  Resolution: %0.1lf x %0.1lf x %0.1lf\n", resolution[OR_X], resolution[OR_Y], resolution[OR_Z]);
        printf("  Volume: %ld x %ld x %ld\n", volume_size[OR_X], volume_size[OR_Y], volume_size[OR_Z]);
        printf("  Block Sizes: %ld x %ld x %ld\n", block_size[OR_X], block_size[OR_Y], block_size[OR_Z]);
    }

    // read the block data
    if (!ReadBlockData()) exit(-1);

    /////////////////////////////////
    //// Read in the voxel files ////
    /////////////////////////////////
/*
    ReadSomaeData();
    ReadSurfaceData();
    ReadSynapseData();
    ReadConnectomeData();
*/

    // set world box
    world_box = R3Box(0, 0, 0, resolution[OR_X] * block_size[OR_X], resolution[OR_Y] * block_size[OR_Y], resolution[OR_Z] * block_size[OR_Z]);

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