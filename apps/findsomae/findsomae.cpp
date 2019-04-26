// Source file for connectome visualization



// include files

#include "R3Graphics/R3Graphics.h"
#include "RNDataStructures/RNDataStructures.h"
#include "fglut/fglut.h"
#include <fstream>
#include <vector>


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



// program variables

static double resolution[3] = { 30, 8, 8 };
static long grid_size[3] = { -1, -1, -1 };
static R3Affine transformation = R3null_affine;
static R3Viewer *viewer = NULL;
static R3Box world_box = R3null_box;



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

static std::vector<long> surfaces;
static long soma_index = -1;


// display variables

static const int ncolor_opts = 24;
static int show_bbox = 1;
static int segmentation_index = 1;
static RNScalar downsample_rate = 20.0;
static int color_cycle = 0;



////////////////////////////////////////////////////////////////////////
// Input/output functions
////////////////////////////////////////////////////////////////////////


static int ReadSurfaceData(void)
{
    surfaces = std::vector<long>();

    char surface_filename[4096];
    sprintf(surface_filename, "original_data/surfaces/%s/%06d.pts", prefix, segmentation_index);

    FILE *fp = fopen(surface_filename, "rb");
    if (!fp) { fprintf(stderr, "Failed to read %s.\n", surface_filename); return 0; }

    long nsurface_points;
    if (fread(&(grid_size[OR_Z]), sizeof(long), 1, fp) != 1) { fprintf(stderr, "Failed to read %s.\n", surface_filename); return 0; }
    if (fread(&(grid_size[OR_Y]), sizeof(long), 1, fp) != 1) { fprintf(stderr, "Failed to read %s.\n", surface_filename); return 0; }
    if (fread(&(grid_size[OR_X]), sizeof(long), 1, fp) != 1) { fprintf(stderr, "Failed to read %s.\n", surface_filename); return 0; }
    if (fread(&nsurface_points, sizeof(long), 1, fp) != 1) { fprintf(stderr, "Failed to read %s.\n", surface_filename); return 0; }
    
    for (long iv = 0; iv < nsurface_points; ++iv) {
        long voxel_index;
        if (fread(&voxel_index, sizeof(long), 1, fp) != 1) { fprintf(stderr, "Failed to read %s.\n", surface_filename); return 0; }
        surfaces.push_back(voxel_index);
    }

    // return success
    return 1;
}

static int ReadSomaData(void)
{
    soma_index = 1;

    char soma_filename[4096];
    sprintf(soma_filename, "somae/%s/%06d.pts", prefix, segmentation_index);

    FILE *fp = fopen(soma_filename, "rb");
    if (!fp) { soma_index = -1; return 0; }

    long nsoma_points;
    if (fread(&(grid_size[OR_Z]), sizeof(long), 1, fp) != 1) { fprintf(stderr, "Failed to read %s.\n", soma_filename); return 0; }
    if (fread(&(grid_size[OR_Y]), sizeof(long), 1, fp) != 1) { fprintf(stderr, "Failed to read %s.\n", soma_filename); return 0; }
    if (fread(&(grid_size[OR_X]), sizeof(long), 1, fp) != 1) { fprintf(stderr, "Failed to read %s.\n", soma_filename); return 0; }
    if (fread(&nsoma_points, sizeof(long), 1, fp) != 1) { fprintf(stderr, "Failed to read %s.\n", soma_filename); return 0; }
    if (fread(&soma_index, sizeof(long), 1, fp) != 1)  { fprintf(stderr, "Failed to read %s.\n", soma_filename); return 0; }

    fclose(fp);
    return 1;
}


static int WriteSomaData(void)
{
    if (soma_index == -1) return 0;

    char soma_filename[4096];
    sprintf(soma_filename, "somae/%s/%06d.pts", prefix, segmentation_index);

    FILE *fp = fopen(soma_filename, "wb");
    if (!fp) { soma_index = -1; return 0; }

    long nsoma_points = 1;
    if (fwrite(&(grid_size[OR_Z]), sizeof(long), 1, fp) != 1) { fprintf(stderr, "Failed to write %s.\n", soma_filename); return 0; }
    if (fwrite(&(grid_size[OR_Y]), sizeof(long), 1, fp) != 1) { fprintf(stderr, "Failed to read %s.\n", soma_filename); return 0; }
    if (fwrite(&(grid_size[OR_X]), sizeof(long), 1, fp) != 1) { fprintf(stderr, "Failed to read %s.\n", soma_filename); return 0; }
    if (fwrite(&nsoma_points, sizeof(long), 1, fp) != 1) { fprintf(stderr, "Failed to read %s.\n", soma_filename); return 0; }
    
    if (fwrite(&soma_index, sizeof(long), 1, fp) != 1)  { fprintf(stderr, "Failed to read %s.\n", soma_filename); return 0; }
    fclose(fp);
    return 1;
}



////////////////////////////////////////////////////////////////////////
// Helper functions
////////////////////////////////////////////////////////////////////////

static void IndexToIndices(long index, long& ix, long& iy, long& iz)
{
  // Set indices of grid value at index
  iz = index / (grid_size[OR_X] * grid_size[OR_Y]);
  iy = (index - iz * grid_size[OR_X] * grid_size[OR_Y]) / grid_size[OR_X];
  ix = index % grid_size[OR_X];
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
// Drawing utility functions
////////////////////////////////////////////////////////////////////////

static void DrawSurface(void)
{
    transformation.Push();
    glPointSize(1.0);
    glBegin(GL_POINTS);
    for (unsigned long iv = 0; iv < surfaces.size(); ++iv) {
        // faster rendering with downsampling
        if (RNRandomScalar() > 1.0 / downsample_rate) continue;

        // get the coordinates from the linear index
        long ix, iy, iz;
        IndexToIndices(surfaces[iv], ix, iy, iz);
        glVertex3f(ix, iy, iz);
    }
    glEnd();
    transformation.Pop();
}


static void DrawSoma(void)
{
    transformation.Push();
    glPointSize(10.0);
    glBegin(GL_POINTS);
    long ix, iy, iz;
    IndexToIndices(soma_index, ix, iy, iz);
    glVertex3f(ix, iy, iz);
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
    if(show_bbox) {
        RNLoadRgb(RNRgb(not background_color[0], not background_color[1], not background_color[2]));
        world_box.Outline();
    }

    // draw segmentation
    RNLoadRgb(Color(segmentation_index));
    DrawSurface();
    RNLoadRgb(RNRgb(not background_color[0], not background_color[1], not background_color[2]));

    if (soma_index != -1) DrawSoma();

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
            WriteSomaData();
            --segmentation_index;
            ReadSurfaceData();
            ReadSomaData();
            printf("Segment %d\n", segmentation_index);                

            break;
        }

        case GLUT_KEY_RIGHT: {
            WriteSomaData();
            ++segmentation_index;
            ReadSurfaceData();
            ReadSomaData();
            printf("Segment %d\n", segmentation_index);

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

        case 'E':
        case 'e': {
            R3Ray world_ray = viewer->WorldRay(x, y);

            RNLength closest_distance = DBL_MAX;
            R3Point closest_point;

            static const int buffer = 40;

            for (unsigned long is = 0; is < surfaces.size(); ++is) {
                long iv = surfaces[is];

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

            for (unsigned long is = 0 ; is < surfaces.size(); ++is) {
                long iv = surfaces[is];

                long ix, iy, iz;
                IndexToIndices(iv, ix, iy, iz);

                R3Point location = R3Point(ix, iy, iz);
                location.Transform(transformation);

                RNLength distance = R3Distance(world_ray, location);
                
                long deltaz = (long) (closest_point.Z() - iz + 0.5);
                long deltay = (long) (closest_point.Y() - iy + 0.5);
                long deltax = (long) (closest_point.X() - ix + 0.5);

                RNScalar normalized_distance = sqrt(resolution[OR_Z] * resolution[OR_Z] * deltaz * deltaz + resolution[OR_Y] * resolution[OR_Y] * deltay * deltay + resolution[OR_X] * resolution[OR_X] * deltax * deltax);
                if (normalized_distance < buffer) continue;

                if (distance < closest_distance) {
                    closest_distance = distance;
                    second_closest_point = R3Point(ix, iy, iz);
                }
            }
            R3Point soma_point = (closest_point + second_closest_point) / 2;
            long soma_ix = soma_point.X();
            long soma_iy = soma_point.Y();
            long soma_iz = soma_point.Z();

            soma_index = soma_iz * grid_size[OR_X] * grid_size[OR_Y] + soma_iy * grid_size[OR_X] + soma_ix;

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
            else if (!strcmp(*argv, "-start_index")) { argv++; argc--; segmentation_index = atoi(*argv); }
            else if (!strcmp(*argv, "-resolution")) { 
                argv++; argc--;
                resolution[OR_X] = atof(*argv);
                argv++; argc--;
                resolution[OR_Y] = atof(*argv);
                argv++; argc--;
                resolution[OR_Z] = atof(*argv);
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

    /////////////////////////////////
    //// Read in the voxel files ////
    /////////////////////////////////

    ReadSurfaceData();
    ReadSomaData();

    // set world box
    world_box = R3Box(0, 0, 0, resolution[OR_X] * grid_size[OR_X], resolution[OR_Y] * grid_size[OR_Y], resolution[OR_Z] * grid_size[OR_Z]);

    // get the transformation
    transformation = R3Affine(R4Matrix(resolution[OR_X], 0, 0, 0, 0, resolution[OR_Y], 0, 0, 0, 0, resolution[OR_Z], 0, 0, 0, 0, 1));

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