// Source file for the mesh viewer program

// include files
#include "RNDataStructures/RNDataStructures.h"
#include <fglut/fglut.h>

// program variables

static char *image_filename = NULL;
static char *image_dataset = NULL;
static char *segment_filename = NULL;
static char *segment_dataset = NULL;
static int print_verbose = 0;
static int print_debug = 0;
static int projection_dim = RN_Z;
static bool affinity = false;


// GLUT variables

static int GLUTwindow = 0;
static int GLUTwindow_height = 800;
static int GLUTwindow_width = 800;
static int GLUTmouse[2] = {0, 0};
static int GLUTbutton[3] = {0, 0, 0};
static int GLUTmodifiers = 0;



// grid viewing variables

static double alpha = 0.70;
static R3Grid *image_grid = NULL;
static R3Grid *segment_grid = NULL;

static R2Grid *selected_image_slice = NULL;
static R2Grid *selected_segment_slice = NULL;

static int selected_slice_index = 0;
static R2Box selected_slice_window = R2null_box;
static R2Point selected_slice_position(RN_UNKNOWN, RN_UNKNOWN);


////////////////////////////////////////////////////////////////////////
// Utility helper functions
////////////////////////////////////////////////////////////////////////

static RNRgb
Color(RNScalar value, bool image_type)
{
    // check for unknown value
    if (value == R2_GRID_UNKNOWN_VALUE || value == 0.0) { return RNblack_rgb; }

    // compute color
    RNRgb c(0, 0, 0);
    if (image_type) {
        c[0] = value / 255.0;
        c[1] = value / 255.0;
        c[2] = value / 255.0;
    }
    else {
        unsigned long integral_value = (unsigned long) (value + 0.5);

        c[0] = (((107 * integral_value) % 700) % 255) / 255.0;
        c[1] = (((509 * integral_value) % 900) % 255) / 255.0;
        c[2] = (((200 * integral_value) % 777) % 255) / 255.0;
    }

    // return color
    return c;
}



static void
SelectGrid(int index)
{
    RNTime start_time;
    start_time.Read();
    // check index
    if (index < 0)
        index = 0;
    if (index > image_grid->Resolution(projection_dim) - 1)
        index = image_grid->Resolution(projection_dim) - 1;

    R2Grid *image_slice = image_grid->Slice(projection_dim, index);
    R2Grid *segment_slice = segment_grid->Slice(projection_dim, index);

    // set window title
    char title[4096];
    sprintf(title, "Overlay Viewer %s %s Slice %d", image_filename, segment_filename, index);
    glutSetWindowTitle(title);

    // update display variables
    if (!selected_image_slice || (selected_image_slice->XResolution() != image_slice->XResolution()) || (selected_image_slice->YResolution() != image_slice->YResolution()))
    {
        RNScalar window_aspect = (double)GLUTwindow_width / (double)GLUTwindow_height;
        RNScalar grid_aspect = (double)image_grid->XResolution() / (double)image_grid->YResolution();
        R2Point origin = image_slice->GridBox().Centroid();
        R2Vector diagonal = image_slice->GridBox().Max() - origin;
        diagonal[0] *= window_aspect / grid_aspect;
        selected_slice_window = R2Box(origin - diagonal, origin + diagonal);
        selected_slice_position.Reset(RN_UNKNOWN, RN_UNKNOWN);
    }

    // remove the old slices
    if (selected_image_slice) delete selected_image_slice;
    if (selected_segment_slice) delete selected_segment_slice;

    // update selected grid
    selected_image_slice = image_slice;
    selected_segment_slice = segment_slice;
    selected_slice_index = index;
}



////////////////////////////////////////////////////////////////////////
// GLUT functions
////////////////////////////////////////////////////////////////////////

void GLUTDrawText(const R2Point &p, const char *s)
{
    // draw text string s and position p
    glRasterPos2d(p[0], p[1]);
    while (*s)
        glutBitmapCharacter(GLUT_BITMAP_HELVETICA_18, *(s++));
}



void GLUTStop(void)
{
    // destro window
    glutDestroyWindow(GLUTwindow);

    // exit
    exit(0);
}



void GLUTRedraw(void)
{
    // check grid
    if (!selected_image_slice || !selected_segment_slice) return;

    // clear window
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    // set projection matrix
    glMatrixMode(GL_PROJECTION);
    glPushMatrix();
    glLoadIdentity();
    gluOrtho2D(selected_slice_window.XMin(), selected_slice_window.XMax(), selected_slice_window.YMin(), selected_slice_window.YMax());

    // set model view matrix
    glMatrixMode(GL_MODELVIEW);
    glPushMatrix();
    glLoadIdentity();

    // draw grid values
    RNTime start_time;
    start_time.Read();
    int xmin = (selected_slice_window.XMin() > 1) ? selected_slice_window.XMin() : 1;
    int ymin = (selected_slice_window.YMin() > 1) ? selected_slice_window.YMin() : 1;
    int xmax = (selected_slice_window.XMax() + 1 < selected_image_slice->XResolution() - 1) ? selected_slice_window.XMax() + 1 : selected_image_slice->XResolution() - 1;
    int ymax = (selected_slice_window.YMax() + 1 < selected_image_slice->YResolution() - 1) ? selected_slice_window.YMax() + 1 : selected_image_slice->YResolution() - 1;
    
    for (int j = ymin; j <= ymax; j += 1)
    {
        glBegin(GL_TRIANGLE_STRIP);

        for (int i = xmin; i <= xmax; i += 1)
        {
            for (int k = -1; k <= 0; k++)
            {
                RNRgb image_color, seg_color;
                image_color = Color(selected_image_slice->GridValue(i, j + k), true);
                seg_color = Color(selected_segment_slice->GridValue(i, j + k), false);

                RNRgb color;
                double epsilon = 1e-2;
                if (seg_color[0] < epsilon and seg_color[1] < epsilon and seg_color[2] < epsilon) {
                    color = image_color;
                }
                else {
                    color = alpha * image_color + (1.0 - alpha) * seg_color;
                }
                

                RNLoadRgb(color);
                // this transforms the coordinate system to top
                glVertex2i(i, (ymax) - (j + k));
            }
        }
        glEnd();
    }

    // draw value at selected grid position
    if ((selected_slice_position.X() != RN_UNKNOWN) && (selected_slice_position.Y() != RN_UNKNOWN))
    {
        int ix = (int)(selected_slice_position.X() + 0.5);
        int iy = (int)(ymax - selected_slice_position.Y() + 0.5);

        char buffer[1024];
        sprintf(buffer, "%d %d: %d", ix, iy, (int)(selected_segment_slice->GridValue(ix, iy) + 0.5));

        RNLoadRgb(RNmagenta_rgb);
        R2Box(selected_slice_position - 0.5 * R2ones_vector, selected_slice_position + 0.5 * R2ones_vector);
        GLUTDrawText(selected_slice_position + 2 * R2ones_vector, buffer);
    }

    // reset projection matrix
    glMatrixMode(GL_PROJECTION);
    glPopMatrix();

    // reset model view matrix
    glMatrixMode(GL_MODELVIEW);
    glPopMatrix();

    // swap buffers
    glutSwapBuffers();
}



void GLUTResize(int w, int h)
{
    // resize window
    glViewport(0, 0, w, h);

    // remember window size
    GLUTwindow_width = w;
    GLUTwindow_height = h;

    // update selected grid window
    if (selected_image_slice)
    {
        RNScalar window_aspect = (double)GLUTwindow_width / (double)GLUTwindow_height;
        RNScalar grid_aspect = (double)selected_image_slice->XResolution() / (double)selected_image_slice->YResolution();
        R2Point origin = selected_image_slice->GridBox().Centroid();
        R2Vector diagonal = selected_image_slice->GridBox().Max() - origin;
        diagonal[0] *= window_aspect / grid_aspect;
        selected_slice_window = R2Box(origin - diagonal, origin + diagonal);
    }

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

    // view manipulation
    if (selected_image_slice)
    {
        if (GLUTbutton[0])
        {
            // query
            RNScalar px = x * selected_slice_window.XLength() / (double)GLUTwindow_width + selected_slice_window.XMin();
            RNScalar py = y * selected_slice_window.YLength() / (double)GLUTwindow_height + selected_slice_window.YMin();
            selected_slice_position.Reset(px, py);
            glutPostRedisplay();
        }
        else if (GLUTbutton[1])
        {
            // zoom
            RNScalar scale_factor = 1;
            scale_factor *= 1.0 - (double)dx / (double)GLUTwindow_width;
            scale_factor *= 1.0 - (double)dy / (double)GLUTwindow_height;
            scale_factor *= scale_factor;
            selected_slice_window.Inflate(scale_factor);
            glutPostRedisplay();
        }
        else if (GLUTbutton[2])
        {
            // pan
            RNScalar tx = -dx * selected_slice_window.XLength() / (double)GLUTwindow_width;
            RNScalar ty = -dy * selected_slice_window.YLength() / (double)GLUTwindow_height;
            selected_slice_window.Translate(R2Vector(tx, ty));
            glutPostRedisplay();
        }
    }

    // remember mouse position
    GLUTmouse[0] = x;
    GLUTmouse[1] = y;
}



void GLUTMouse(int button, int state, int x, int y)
{
    // invert y coordinate
    y = GLUTwindow_height - y;

    // process mouse button event
    if (button == 0)
    {
        if (state == GLUT_DOWN)
        {
            // query
            RNScalar px = x * selected_slice_window.XLength() / (double)GLUTwindow_width + selected_slice_window.XMin();
            RNScalar py = y * selected_slice_window.YLength() / (double)GLUTwindow_height + selected_slice_window.YMin();
            selected_slice_position.Reset(px, py);
            glutPostRedisplay();
        }
        else
        {
            selected_slice_position.Reset(RN_UNKNOWN, RN_UNKNOWN);
            glutPostRedisplay();
        }
    }
    else if ((button == 3) || (button == 4))
    {
        if (state == GLUT_DOWN)
        {
            // zoom with wheel
            RNScalar scale_factor = (button == 3) ? 0.9 : 1.1;
            selected_slice_window.Inflate(scale_factor);
            glutPostRedisplay();
        }
    }

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

    // process keyboard button event
    switch (key)
    {
    case GLUT_KEY_PAGE_UP:
    case GLUT_KEY_PAGE_DOWN:
    case GLUT_KEY_UP:
    case GLUT_KEY_DOWN:
    case GLUT_KEY_LEFT:
    case GLUT_KEY_RIGHT:
    {
        if (selected_image_slice)
        {
            int shift = 0;
            if (key == GLUT_KEY_PAGE_DOWN)
                shift = -10;
            else if (key == GLUT_KEY_PAGE_UP)
                shift = 10;
            else if (key == GLUT_KEY_UP)
                shift += 1;
            else if (key == GLUT_KEY_DOWN)
                shift -= 1;
            else if (key == GLUT_KEY_LEFT) {
                alpha -= 0.05;
                if (alpha < 0) alpha = 0.0;;
            }
            else if (key == GLUT_KEY_RIGHT) {
                alpha += 0.05;
                if (alpha > 1) alpha = 1.0;
            }
            SelectGrid(selected_slice_index + shift);
            glutPostRedisplay();
        }
        break;
    }
    }

    // remember mouse position
    GLUTmouse[0] = x;
    GLUTmouse[1] = y;

    // remember modifiers
    GLUTmodifiers = glutGetModifiers();

    // redraw
    glutPostRedisplay();
}



void GLUTKeyboard(unsigned char key, int x, int y)
{
    // invert y coordinate
    y = GLUTwindow_height - y;

    // process keyboard button event
    switch (key)
    {
        case 'X':
        case 'x':
        {
            projection_dim = RN_X;
            SelectGrid(0);
            break;
        }

        case 'Y':
        case 'y':
        {
            projection_dim = RN_Y;
            SelectGrid(0);
            break;
        }

        case 'Z':
        case 'z':
        {
            projection_dim = RN_Z;
            SelectGrid(0);
            break;
        }

        case 27:
        {
            // ESCAPE
            GLUTStop();
            break;
    }
    }

    // remember mouse position
    GLUTmouse[0] = x;
    GLUTmouse[1] = y;

    // remember modifiers
    GLUTmodifiers = glutGetModifiers();

    // redraw
    glutPostRedisplay();
}



void GLUTInit(int *argc, char **argv)
{
    // set window dimensions
    RNScalar aspect = (RNScalar)image_grid->YResolution() / (RNScalar)image_grid->XResolution();
    GLUTwindow_width = image_grid->XResolution();
    GLUTwindow_height = image_grid->YResolution();
    GLUTwindow_height = aspect * GLUTwindow_width;

    // open window
    glutInit(argc, argv);
    glutInitWindowPosition(100, 100);
    glutInitWindowSize(GLUTwindow_width, GLUTwindow_height);
    glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH); // | GLUT_STENCIL
    GLUTwindow = glutCreateWindow("Overlay Viewer");

    // initialize background color
    glClearColor(00, 0.0, 0.0, 1.0);

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
    // Select first grid
    SelectGrid(selected_slice_index);
    // Run main loop -- never returns
    glutMainLoop();
}



////////////////////////////////////////////////////////////////////////
// Program argument parsing
////////////////////////////////////////////////////////////////////////

static int
ParseArgs(int argc, char **argv)
{
    // parse arguments
    argc--; argv++;
    while (argc > 0)
    {
        if ((*argv)[0] == '-') {
            if (!strcmp(*argv, "-v")) print_verbose = 1;
            else if (!strcmp(*argv, "-debug")) print_debug = 1;
            else if (!strcmp(*argv, "-affinities")) affinity = true;
            else { fprintf(stderr, "Invalid program argument: %s\n", *argv); return 0; }
        }
        else {
            if (!image_filename) { image_filename = *argv; }
            else if (!image_dataset) { image_dataset = *argv; }
            else if (!segment_filename) { segment_filename = *argv; }
            else if (!segment_dataset) { segment_dataset = *argv; }
            else { fprintf(stderr, "Invalid program argument: %s\n", *argv); return 0; }
        } 
        argv++; argc--;
    }

    // check filenames
    if (!image_filename || !image_dataset || !segment_filename || !segment_dataset) { fprintf(stderr, "Usage: overlay image_filename image_dataset segment_filename segment_dataset [options]\n"); return 0; }

    // return OK status
    return 1;
}



////////////////////////////////////////////////////////////////////////
// Main
////////////////////////////////////////////////////////////////////////

int main(int argc, char **argv)
{
    // parse program arguments
    if (!ParseArgs(argc, argv))
        exit(-1);

    RNTime start_time;
    start_time.Read();

    R3Grid **image_grids = RNReadH5File(image_filename, image_dataset);
    R3Grid **segment_grids = RNReadH5File(segment_filename, segment_dataset);
    image_grid = image_grids[0];
    segment_grid = segment_grids[0];

    if ((image_grid->XResolution() != segment_grid->XResolution()) || (image_grid->YResolution() != segment_grid->YResolution()) || (image_grid->ZResolution() != segment_grid->ZResolution())) {
        fprintf(stderr, "Image and segmentation are variable sizes...\n");
        return -1;
    }

    if (affinity) {
        for (int iz = 0; iz < image_grid->ZResolution(); ++iz) {
            for (int iy = 0; iy < image_grid->YResolution(); ++iy) {
                for (int ix = 0; ix < image_grid->XResolution(); ++ix) {
                    image_grid->SetGridValue(ix, iy, iz, 255 * image_grid->GridValue(ix, iy, iz));
                }
            }
        }
    }

    // free memory
    delete[] image_grids;
    delete[] segment_grids;

    if (print_verbose) {
        printf("Read h5 files in %0.2f seconds\n", start_time.Elapsed());
        printf("  Resolution: %d %d %d\n", image_grid->XResolution(), image_grid->YResolution(), image_grid->ZResolution());
    }

    // initialize GLUT
    GLUTInit(&argc, argv);

    // run GLUT interface
    GLUTMainLoop();

    // return success
    return 0;
}