// Source file for skeleton visualizer

// include files

#include "R3Graphics/R3Graphics.h"
#include "RNDataStructures/RNDataStructures.h"
#include "fglut/fglut.h"
#include <fstream>
#include <string>
#include <vector>
#include <algorithm>
#include <dirent.h>

// GLUT defines

#define ENTER 13
#define ESCAPE 27
#define SPACEBAR 32

// program arguments

static int print_debug = 0;
static int print_verbose = 0;
static const char* prefix = NULL;

// program variables

static int nnet_resolution[3] = { 51, 51, 51 };
static int resolution[3] = { -1, -1, -1 };
static R3Affine transformation = R3null_affine;
static R3Viewer* viewer = NULL;
static R3Point selected_position;
static R3Box world_box;
static int selected_voxel[3];

// voxel grids

static R3Grid *grid = NULL;
static R3Grid **grids = NULL;
static int ngrids = -1;



// always update the grid and the interval
static R2Grid* selected_slice = NULL;
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

static int projection_dim = RN_Z;
static int selected_slice_index = 0;

// display dimenson variables

static int show_bbox = 1;
static int show_slice = 0;
static int projection = 0;

// display machine label variables

static int feature_index = 0;
static RNBoolean *feature_labels = NULL;


////////////////////////////////////////////////////////////////////////
// Input/output functions
////////////////////////////////////////////////////////////////////////

static int ReadData(void)
{
    // start statistics
    RNTime start_time;
    start_time.Read();

    char directory_name[4096];
    sprintf(directory_name, "neural_input/%s", prefix);

    std::vector<std::string> filenames = std::vector<std::string>();

    DIR *dir;
    struct dirent *ent;
    if ((dir = opendir (directory_name)) != NULL) {
        // print all the files and directories within directory
        while ((ent = readdir (dir)) != NULL) {
            if (ent->d_name[0] == '.') continue;
            filenames.push_back(std::string(ent->d_name));
        }
        closedir (dir);
    } else {
        // could not open directory
        fprintf(stderr, "Failed to read %s\n", directory_name);
        return 0; 
    }
    
    // sort the filenames
    sort(filenames.begin(), filenames.end());
    
    ngrids = filenames.size() - 1;    
    grids = new R3Grid *[ngrids];
    
    for (int iv = 0; iv < ngrids; ++iv) {
        char filename[4096];
        sprintf(filename, "%s/%s", directory_name, filenames[iv].c_str());
        
        FILE *fp = fopen(filename, "rb");
        if (!fp) { fprintf(stderr, "Failed to read %s\n", filename); return 0; }
        
        grids[iv] = new R3Grid(nnet_resolution[RN_X], nnet_resolution[RN_Y], nnet_resolution[RN_Z]);
        if (!grids[iv]) { fprintf(stderr, "Failed to allocate memory for grid\n"); return 0; }
        
        for (int iz = 0; iz < nnet_resolution[RN_Z]; ++iz) {
            for (int iy = 0; iy < nnet_resolution[RN_Y]; ++iy) {
                for (int ix = 0; ix < nnet_resolution[RN_X]; ++ix) {
                    int grid_value;
                    if (fread(&grid_value, sizeof(int), 1, fp) != 1) { 
                        fprintf(stderr, "Failed to read %s\n", filename);
                        return 0;
                    }
                    
                    grids[iv]->SetGridValue(ix, iy, iz, grid_value);
                }
            }
        }
        
        fclose(fp);
    }
    
    resolution[RN_X] = grids[0]->XResolution();
    resolution[RN_Y] = grids[0]->YResolution();
    resolution[RN_Z] = grids[0]->ZResolution();
    
    grid = grids[0];
    
    // create feature labels
    feature_labels = new RNBoolean[ngrids];
    if (!feature_labels) { fprintf(stderr, "Failed to allocate memory for feature labels\n"); return 0; }
    
    // open feature filename
    char labels_filename[4096];
    sprintf(labels_filename, "%s/%s", directory_name, filenames[ngrids].c_str());
    
    FILE *labels_fp = fopen(labels_filename, "rb");
    if (!labels_fp) { fprintf(stderr, "Failed to read %s\n", labels_filename); return 0; }
    
    for (int iv = 0; iv < ngrids; ++iv) {
        if (fread(&(feature_labels[iv]), sizeof(int), 1, labels_fp) != 1) {
            fprintf(stderr, "Failed to read %s\n", labels_filename); return 0;
        }
    }
    
    // close file
    fclose(labels_fp);
    
    if (print_verbose) {
        printf("Read %d features grids...\n", ngrids);
        printf("  Time = %0.2f seconds\n", start_time.Elapsed());
        printf("  Resolution = (%d, %d, %d)\n", nnet_resolution[RN_X], nnet_resolution[RN_Y], nnet_resolution[RN_Z]);
    }

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
    int xoffset = (GLUTwindow_width - projection_scale * grid->Resolution(window_xdim)) / 2;
    int yoffset = (GLUTwindow_height - projection_scale * grid->Resolution(window_ydim)) / 2;

    // convert from grid to world coordinates
    return R2Point(grid_point.X() * projection_scale + xoffset,
        grid_point.Y() * projection_scale + yoffset);
}

static R2Point ConvertWorldToGrid(R2Point world_point)
{
    // projection offsets
    int window_xdim = (projection_dim + 1) % 3;
    int window_ydim = (projection_dim + 2) % 3;
    int xoffset = (GLUTwindow_width - projection_scale * grid->Resolution(window_xdim)) / 2;
    int yoffset = (GLUTwindow_height - projection_scale * grid->Resolution(window_ydim)) / 2;

    // convert from world to grid coordinates
    return R2Point((world_point.X() - xoffset) / (projection_scale),
        (world_point.Y() - yoffset) / (projection_scale));
}

static void UpdateSlices(void)
{
    selected_slice = grid->Slice(projection_dim, selected_slice_index);
    selected_slice_range = selected_slice->Range();
}

////////////////////////////////////////////////////////////////////////
// Drawing utility functions
////////////////////////////////////////////////////////////////////////

static RNRgb Color(RNScalar value, RNInterval range)
{
    // Check for unknown value
    if(value == R2_GRID_UNKNOWN_VALUE) {
        if(color_type == 0)
            return RNRgb(1, 0.5, 0);
        else
            return RNblack_rgb;
    }

    unsigned int index = (unsigned int)(value + 0.5);

    RNRgb color = RNRgb(((107 * index) % 700) / 700.0, ((509 * index) % 900) / 900.0, ((200 * index) % 777) / 777.0);

    // Return color
    return color;
}

static void GLUTDrawText(const R2Point& position, const char* s)
{
    // draw text string s at position
    glRasterPos2d(position[0], position[1]);
    while(*s) glutBitmapCharacter(GLUT_BITMAP_HELVETICA_18, *(s++));
}


static void DrawSlice(void)
{
    transformation.Push();

    RNLoadRgb(1.0, 1.0, 1.0);
    // draw the selected slice
    grid->DrawSlice(projection_dim, selected_slice_index);
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
    if((selected_slice_position.X() != RN_UNKNOWN) && (selected_slice_position.Y() != RN_UNKNOWN)) {
        int ix = (int)(selected_slice_position.X() + 0.5);
        int iy = (int)(selected_slice_position.Y() + 0.5);
        RNScalar value = selected_slice->GridValue(ix, selected_slice->YResolution() - iy);

        // create hover text
        char buffer[1024];
        if(value != R2_GRID_UNKNOWN_VALUE)
            sprintf(buffer, "%d %d: %g", ix, iy, value);
        else
            sprintf(buffer, "%d %d: %s", ix, iy, "Unknown");
        RNLoadRgb(RNblue_rgb);

        // move point to world location
        R2Point world_position = ConvertGridToWorld(selected_slice_position);
        GLUTDrawText(world_position + 2 * R2ones_vector, buffer);
    }

    // draw the actual slice
    for(int iy = 1; iy < selected_slice->YResolution(); ++iy) {
        glBegin(GL_TRIANGLE_STRIP);
        for(int ix = 1; ix < selected_slice->XResolution(); ++ix) {
            for(int k = -1; k <= 0; ++k) {
                // convert from world to grid coordinates
                R2Point grid_position = R2Point(ix, selected_slice->YResolution() - (iy + k));
                R2Point world_position = ConvertGridToWorld(grid_position);

                // get the color for this area
                RNRgb color = Color(selected_slice->GridValue(ix, iy + k), selected_slice_range);
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
    if(show_bbox) {
        RNLoadRgb(RNwhite_rgb);
        world_box.Outline();
    }

    // draw the actual points in 3D
    glBegin(GL_POINTS);
    for (int iz = 0; iz < grid->ZResolution(); ++iz) {
        for (int iy = 0; iy < grid->YResolution(); ++iy) {
            for (int ix = 0; ix < grid->XResolution(); ++ix) {
                int grid_value = (int)(grid->GridValue(ix, iy, iz) + 0.5);
                if (grid_value == 1) {
                    if (feature_labels[feature_index]) RNLoadRgb(RNgreen_rgb);
                    else RNLoadRgb(RNred_rgb);
                    glVertex3f(ix, iy, iz);
                }
                else if (grid_value == 2) {
                    if (feature_labels[feature_index]) RNLoadRgb(RNblue_rgb);
                    else RNLoadRgb(RNyellow_rgb);
                    glVertex3f(ix, iy, iz);
                }
            }
        }
    }
    glEnd();

    if(show_slice) DrawSlice();

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

    // delete[] machine_label;

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
    if(projection)
        Draw2D();
    else
        Draw3D();

    // set window title
    char title[4096];
    sprintf(title, "Feature visualizer");
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
    if(GLUTbutton[0]) {
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
        if(gridx < 0 || gridx >= grid->Resolution(window_xdim))
            selected_slice_position.Reset(RN_UNKNOWN, RN_UNKNOWN);
        else if(gridy < 0 || gridy >= grid->Resolution(window_ydim))
            selected_slice_position.Reset(RN_UNKNOWN, RN_UNKNOWN);

        glutPostRedisplay();
    } else if(GLUTbutton[1]) {
        // zoom
    } else if(GLUTbutton[2]) {
        // move voxel
        R2Point world_point = R2Point(x, y);
        R2Point grid_point = ConvertWorldToGrid(world_point);

        int window_xdim = (projection_dim + 1) % 3;
        int window_ydim = (projection_dim + 2) % 3;

        if(grid_point.X() < 0) grid_point.SetCoord(RN_X, 0);
        if(grid_point.Y() < 0) grid_point.SetCoord(RN_Y, 0);
        if(grid_point.X() > grid->Resolution(window_xdim) - 1)
            grid_point.SetCoord(RN_X, grid->Resolution(window_xdim) - 1);
        if(grid_point.Y() > grid->Resolution(window_ydim) - 1)
            grid_point.SetCoord(RN_Y, grid->Resolution(window_ydim) - 1);

        R3Point intersection;
        if(projection_dim == RN_X)
            intersection = R3Point(selected_slice_index, grid_point.X(), grid_point.Y());
        else if(projection_dim == RN_Y)
            intersection = R3Point(grid_point.Y(), selected_slice_index, grid_point.X());
        else
            intersection = R3Point(grid_point.X(), grid_point.Y(), selected_slice_index);

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
    if(GLUTbutton[0])
        viewer->RotateWorld(1.0, origin, x, y, dx, dy);
    else if(GLUTbutton[1])
        viewer->ScaleWorld(1.0, origin, x, y, dx, dy);
    else if(GLUTbutton[2]) {
        if(glutGetModifiers() == GLUT_ACTIVE_SHIFT) {
            viewer->TranslateWorld(1.0, origin, x, y, dx, dy);
        } else {
            R3Vector towards = viewer->Camera().Towards();
            RNDimension dim = towards.MaxDimension();
            R3Vector normal = R3xyz_triad[dim];
            R3Plane plane(selected_position, normal);
            R3Ray viewer_ray = viewer->WorldRay(x, y);
            R3Point intersection;
            if(R3Intersects(viewer_ray, plane, &intersection)) selected_position = intersection;

            // confine point by x dimension
            if(selected_position.X() < world_box.XMin())
                selected_position.SetX(world_box.XMin());
            else if(selected_position.X() >= world_box.XMax())
                selected_position.SetX(world_box.XMax());
            // confine point by y dimension
            if(selected_position.Y() < world_box.YMin())
                selected_position.SetY(world_box.YMin());
            else if(selected_position.Y() >= world_box.YMax())
                selected_position.SetY(world_box.YMax());
            // confine point by z dimension
            if(selected_position.Z() < world_box.ZMin())
                selected_position.SetZ(world_box.ZMin());
            else if(selected_position.Z() >= world_box.ZMax())
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
    if(projection)
        GLUTMotion2D(x, y);
    else
        GLUTMotion3D(x, y);

    // redisplay if a mouse was down
    if(GLUTbutton[0] || GLUTbutton[1] || GLUTbutton[2]) glutPostRedisplay();

    // remember mouse position
    GLUTmouse[0] = x;
    GLUTmouse[1] = y;
}

void GLUTMouse2D(int button, int state, int x, int y)
{
    if(button == 0) {
        if(state == GLUT_DOWN) {
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
            if(gridx < 0 || gridx >= grid->Resolution(window_xdim))
                selected_slice_position.Reset(RN_UNKNOWN, RN_UNKNOWN);
            else if(gridy < 0 || gridy >= grid->Resolution(window_ydim))
                selected_slice_position.Reset(RN_UNKNOWN, RN_UNKNOWN);

            glutPostRedisplay();
        } else {
            selected_slice_position.Reset(RN_UNKNOWN, RN_UNKNOWN);
            glutPostRedisplay();
        }
    }
}

void GLUTMouse3D(int button, int state, int x, int y)
{
    if(button == 2) {
        if(glutGetModifiers() != GLUT_ACTIVE_SHIFT) {
            if(state == GLUT_DOWN) {
                R3Vector towards = viewer->Camera().Towards();
                RNDimension dim = towards.MaxDimension();
                R3Vector normal = R3xyz_triad[dim];
                R3Plane plane(selected_position, normal);
                R3Ray viewer_ray = viewer->WorldRay(x, y);
                R3Point intersection;
                if(R3Intersects(viewer_ray, plane, &intersection)) selected_position = intersection;

                // confine point by x dimension
                if(selected_position.X() < world_box.XMin())
                    selected_position.SetX(world_box.XMin());
                else if(selected_position.X() >= world_box.XMax())
                    selected_position.SetX(world_box.XMax());
                // confine point by y dimension
                if(selected_position.Y() < world_box.YMin())
                    selected_position.SetY(world_box.YMin());
                else if(selected_position.Y() >= world_box.YMax())
                    selected_position.SetY(world_box.YMax());
                // confine point by z dimension
                if(selected_position.Z() < world_box.ZMin())
                    selected_position.SetZ(world_box.ZMin());
                else if(selected_position.Z() >= world_box.ZMax())
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

    // remember mouse position
    GLUTmouse[0] = x;
    GLUTmouse[1] = y;

    // process mouse button event
    if(projection)
        GLUTMouse2D(button, state, x, y);
    else
        GLUTMouse3D(button, state, x, y);

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
    case GLUT_KEY_UP: {
        selected_slice_index++;
        if(selected_slice_index > grid->Resolution(projection_dim) - 1)
            selected_slice_index = grid->Resolution(projection_dim) - 1;
        UpdateSlices();
        break;
    }

    case GLUT_KEY_DOWN: {
        selected_slice_index--;
        if(selected_slice_index < 0) selected_slice_index = 0;
        UpdateSlices();
        break;
    }

    case GLUT_KEY_LEFT: {
        feature_index--;
        if(feature_index < 0) feature_index = 0;
        grid = grids[feature_index];
        UpdateSlices();
        break;
    }

    case GLUT_KEY_RIGHT: {
        feature_index++;
        if(feature_index > ngrids - 1)
            feature_index = ngrids - 1;
        grid = grids[feature_index];
        UpdateSlices();
        break;
    }
    }

    // redraw
    glutPostRedisplay();
}

void GLUTKeyboard2D(unsigned char key, int x, int y)
{
    switch(key) {
 
    case 'C':
    case 'c': {
        color_type = 1 - color_type;
        break;
    }
    }
}

void GLUTKeyboard3D(unsigned char key, int x, int y)
{
    switch(key) {
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
    if(projection) {
        GLUTKeyboard2D(key, x, y);
    } else {
        GLUTKeyboard3D(key, x, y);
    }

    // keys regardless of projection status
    switch(key) {
    case ENTER: {
        projection = 1 - projection;
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
#if(RN_OS == RN_WINDOWSNT)
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

static R3Viewer* CreateViewer(void)
{
    // start statistics
    RNTime start_time;
    start_time.Read();

    if(world_box.IsEmpty()) RNAbort("Error in CreateViewer - box is empty");
    RNLength r = world_box.DiagonalRadius();
    if(r < 0 || RNIsInfinite(r)) RNAbort("Error in CreateViewer - r must be positive finite");

    // set up camera view looking down the z axis
    static R3Vector initial_camera_towards = R3Vector(0.0, 0.0, -1.5);
    static R3Vector initial_camera_up = R3Vector(0.0, 1.0, 0.0);
    R3Point initial_camera_origin = world_box.Centroid() - initial_camera_towards * 2.5 * r;
    R3Camera camera(initial_camera_origin, initial_camera_towards, initial_camera_up, 0.4, 0.4, 0.1 * r, 1000.0 * r);
    R2Viewport viewport(0, 0, GLUTwindow_width, GLUTwindow_height);
    R3Viewer* viewer = new R3Viewer(camera, viewport);

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
    argc--;
    argv++;
    while(argc > 0) {
        if((*argv)[0] == '-') {
            if(!strcmp(*argv, "-v"))
                print_verbose = 1;
            else if(!strcmp(*argv, "-debug"))
                print_debug = 1;
            else if(!strcmp(*argv, "-nnet_resolution")) {
                argv++;
                argc--;
                nnet_resolution[RN_X] = atoi(*argv);
                argv++;
                argc--;
                nnet_resolution[RN_Y] = atoi(*argv);
                argv++;
                argc--;
                nnet_resolution[RN_Z] = atoi(*argv);
            } else {
                fprintf(stderr, "Invalid program argument: %s\n", *argv);
                return 0;
            }
        } else {
            if(!prefix) {
                prefix = *argv;
            } else {
                fprintf(stderr, "Invalid program argument: %s\n", *argv);
                return 0;
            }
        }
        argv++;
        argc--;
    }

    // error if there is no input name
    if(!prefix) {
        fprintf(stderr, "Need to supply a neuron data file\n");
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

    if(!ReadData()) exit(-1);

    // update pointers to selected slice
    UpdateSlices();

    // set world box
    world_box = R3Box(0, 0, 0, grid->XResolution(), grid->YResolution(),
        grid->ZResolution());
    selected_position = world_box.Centroid();

    // get the transformation
    transformation =
        R3Affine(R4Matrix(1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1));

    // initialize selected voxel
    R3Point voxel = selected_position;
    voxel.InverseTransform(transformation);
    selected_voxel[RN_X] = (int)(voxel.X() + 0.5);
    selected_voxel[RN_Y] = (int)(voxel.Y() + 0.5);
    selected_voxel[RN_Z] = (int)(voxel.Z() + 0.5);
    selected_slice_index = selected_voxel[projection_dim];

    // set projection scale variables
    for(int dim = 0; dim <= 2; ++dim) {
        RNScalar ratio = GLUTwindow_height / (grid->Resolution(dim));
        if(ratio < projection_scale) projection_scale = ratio;
    }


    // create viewer
    viewer = CreateViewer();
    if(!viewer) exit(-1);

    // initialize GLUT
    GLUTInit(&argc, argv);

    // run GLUT interface
    GLUTMainLoop();

    // return success
    return 1;
}