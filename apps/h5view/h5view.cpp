// Source file for the mesh viewer program

// include files
#include "RNDataStructures/RNDataStructures.h"
#include <fglut/fglut.h>

// program variables

static char *input_filename = NULL;
static char *dataset_name = NULL;
static int print_verbose = 0;
static int print_debug = 0;
static int projection_dim = RN_Z;
static int scatter = 0;
static int ngrids = 0;

// GLUT variables

static int GLUTwindow = 0;
static int GLUTwindow_height = 800;
static int GLUTwindow_width = 800;
static int GLUTmouse[2] = {0, 0};
static int GLUTbutton[3] = {0, 0, 0};
static int GLUTmodifiers = 0;

// grid viewing variables
static R3Grid **grids = NULL;
static R3Grid *grid = NULL;
static R2Grid *selected_slice = NULL;
static int selected_slice_index = -1;
static R2Box selected_slice_window = R2null_box;
static R2Point selected_slice_position(RN_UNKNOWN, RN_UNKNOWN);
static RNInterval selected_slice_range(0, 0);
static int color_type = 0; // 0=gray, 1=red-green-blue

////////////////////////////////////////////////////////////////////////
// Utility helper functions
////////////////////////////////////////////////////////////////////////

static RNRgb
Color(RNScalar value)
{
    // check for unknown value
    if (value == R2_GRID_UNKNOWN_VALUE)
    {
        if (color_type == 0)
            return RNRgb(1, 0.5, 0);
        else
            return RNblack_rgb;
    }
    else if (value == 0.0)
    {
        return RNblack_rgb;
    }

    // normalize color
    RNScalar value_min = selected_slice_range.Min();
    RNScalar value_width = selected_slice_range.Diameter();
    RNScalar value_scale = (value_width > 0) ? 1.0 / value_width : 1.0;
    RNScalar normalized_value = value_scale * (value - value_min);

    if (scatter)
        normalized_value = (RNScalar)((int)(normalized_value * 65536) % 359) / 359;

    // compute color
    RNRgb c(0, 0, 0);
    if (color_type == 0)
    {
        c[0] = normalized_value;
        c[1] = normalized_value;
        c[2] = normalized_value;
    }
    else
    {
        if (normalized_value < 0.5)
        {
            c[0] = 1 - 2 * normalized_value;
            c[1] = (2 * normalized_value);
        }
        else
        {
            c[1] = (1 - 2 * (normalized_value - 0.5));
            c[2] = (2 * (normalized_value - 0.5));
        }
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
    if (index > grid->Resolution(projection_dim) - 1)
        index = grid->Resolution(projection_dim) - 1;

    R2Grid *slice = grid->Slice(projection_dim, index);
    // set window title
    char title[4096];
    sprintf(title, "H5 Viewer %s Slice %d", input_filename, index);
    glutSetWindowTitle(title);

    // update display variables
    if (!selected_slice || (selected_slice->XResolution() != slice->XResolution()) || (selected_slice->YResolution() != slice->YResolution()) || (0 && (selected_slice->WorldToGridTransformation() != slice->WorldToGridTransformation())))
    {
        RNScalar window_aspect = (double)GLUTwindow_width / (double)GLUTwindow_height;
        RNScalar grid_aspect = (double)grid->XResolution() / (double)grid->YResolution();
        R2Point origin = slice->GridBox().Centroid();
        R2Vector diagonal = slice->GridBox().Max() - origin;
        diagonal[0] *= window_aspect / grid_aspect;
        selected_slice_window = R2Box(origin - diagonal, origin + diagonal);
        selected_slice_position.Reset(RN_UNKNOWN, RN_UNKNOWN);
    }

    // update min and max values
    selected_slice_range = grid->Range();

    // update selected grid
    selected_slice_index = index;
    selected_slice = slice;
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
    if (!selected_slice)
        return;

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
    int xmax = (selected_slice_window.XMax() + 1 < selected_slice->XResolution() - 1) ? selected_slice_window.XMax() + 1 : selected_slice->XResolution() - 1;
    int ymax = (selected_slice_window.YMax() + 1 < selected_slice->YResolution() - 1) ? selected_slice_window.YMax() + 1 : selected_slice->YResolution() - 1;
    
    for (int j = ymin; j <= ymax; j += 1)
    {
        glBegin(GL_TRIANGLE_STRIP);

        for (int i = xmin; i <= xmax; i += 1)
        {
            for (int k = -1; k <= 0; k++)
            {
                RNScalar value = selected_slice->GridValue(i, j + k);
                RNRgb color = Color(value);
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
        int iy = (int)(selected_slice_position.Y() + 0.5);
        RNScalar value = selected_slice->GridValue(ix, ymax - iy);
        char buffer[1024];
        if (value != R2_GRID_UNKNOWN_VALUE)
            sprintf(buffer, "%d %d : %g", ix, iy, value);
        else
            sprintf(buffer, "%d %d : %s", ix, iy, "Unknown");
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
    if (selected_slice)
    {
        RNScalar window_aspect = (double)GLUTwindow_width / (double)GLUTwindow_height;
        RNScalar grid_aspect = (double)selected_slice->XResolution() / (double)selected_slice->YResolution();
        R2Point origin = selected_slice->GridBox().Centroid();
        R2Vector diagonal = selected_slice->GridBox().Max() - origin;
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
    if (selected_slice)
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
    {
        if (selected_slice)
        {
            int shift = 0;
            if (key == GLUT_KEY_PAGE_DOWN)
                shift = -1;
            else if (key == GLUT_KEY_PAGE_UP)
                shift = 1;
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
    case '1':
    case '2':
    case '3':
    case '4':
    {
        int ia = key - '1';
        if (ia >= ngrids)
            break;
        grid = grids[ia];
        SelectGrid(selected_slice_index);
        break;
    }

    case 'C':
    case 'c':
    {
        color_type = ((color_type + 1) % 2);
        break;
    }

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
    RNScalar aspect = (RNScalar)grid->YResolution() / (RNScalar)grid->XResolution();
    GLUTwindow_width = grid->XResolution();
    GLUTwindow_height = grid->YResolution();
    GLUTwindow_height = aspect * GLUTwindow_width;

    // open window
    glutInit(argc, argv);
    glutInitWindowPosition(100, 100);
    glutInitWindowSize(GLUTwindow_width, GLUTwindow_height);
    glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH); // | GLUT_STENCIL
    GLUTwindow = glutCreateWindow("H5 Viewer");

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
    SelectGrid(0);
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
    argc--;
    argv++;
    while (argc > 0)
    {
        if ((*argv)[0] == '-')
        {
            if (!strcmp(*argv, "-v"))
                print_verbose = 1;
            else if (!strcmp(*argv, "-debug"))
                print_debug = 1;
            else if (!strcmp(*argv, "-scatter"))
                scatter = 1;
            else
            {
                fprintf(stderr, "Invalid program argument: %s", *argv);
                exit(1);
            }
        }
        else
        {
            if (!input_filename)
            {
                input_filename = *argv;
                argv++;
                argc--;
                dataset_name = *argv;
            }
            else
            {
                fprintf(stderr, "Invalid program argument: %s", *argv);
                exit(1);
            }
        }
        argv++;
        argc--;
    }

    // check filenames
    if (!input_filename || !dataset_name)
    {
        fprintf(stderr, "Usage: h5view h5file [dataset_name] [options]\n");
        return 0;
    }

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

    grids = RNReadH5File(input_filename, dataset_name);
    grid = grids[0];


    // initialize GLUT
    GLUTInit(&argc, argv);

    // run GLUT interface
    GLUTMainLoop();

    // return success
    return 0;
}