// Source file for GAPS scalar grid class



////////////////////////////////////////////////////////////////////////
// NOTE:
// Grid values are defined as samples at the grid positions ranging from
// (0, 0, 0) to (xres-1, yres-1, zres-1).  Grid values outside this range
// are undefined.
////////////////////////////////////////////////////////////////////////



// Include files

#include "R3Shapes/R3Shapes.h"



R3CharGrid::
R3CharGrid(unsigned long xresolution, unsigned long yresolution, unsigned long zresolution)
{
  // Set grid resolution
  grid_resolution[0] = xresolution;
  grid_resolution[1] = yresolution;
  grid_resolution[2] = zresolution;
  
  grid_row_size = xresolution;
  grid_sheet_size = grid_row_size * yresolution;
  grid_size = grid_sheet_size * zresolution;

  // Allocate grid values
  if (grid_size == 0) grid_values = NULL;
  else grid_values = new unsigned char [ grid_size ];
  assert(!grid_size || grid_values);

  // Set all values to zero
  for (unsigned long i = 0; i < grid_size; i++) grid_values[i] = 0;

  // Set transformations
  grid_to_world_transform = R3identity_affine;
  world_to_grid_transform = R3identity_affine;
  world_to_grid_scale_factor = 1.0;
  grid_to_world_scale_factor = 1.0;
}



R3CharGrid::
R3CharGrid(unsigned long xresolution, unsigned long yresolution, unsigned long zresolution, const R3Box& bbox)
{
  // Set grid resolution
  grid_resolution[0] = xresolution;
  grid_resolution[1] = yresolution;
  grid_resolution[2] = zresolution;
  grid_row_size = xresolution;
  grid_sheet_size = grid_row_size * yresolution;
  grid_size = grid_sheet_size * zresolution;

  // Allocate grid values
  if (grid_size == 0) grid_values = NULL;
  else grid_values = new unsigned char [ grid_size ];
  assert(!grid_size || grid_values);

  // Set all values to zero
  for (unsigned long i = 0; i < grid_size; i++) grid_values[i] = 0;

  // Set transformations
  SetWorldToGridTransformation(bbox);
}



R3CharGrid::
R3CharGrid(const R3CharGrid& voxels)
  : grid_values(NULL)
{
  // Copy everything
  *this = voxels;
}



R3CharGrid::
~R3CharGrid(void)
{
  // Deallocate memory for grid values
  if (grid_values) delete [] grid_values;
}



RNInterval R3CharGrid::
Range(void) const
{
  // Find smallest and largest values
  RNScalar minimum = FLT_MAX;
  RNScalar maximum = -FLT_MAX;
  unsigned char *grid_valuep = grid_values;
  for (unsigned long i = 0; i < grid_size; i++) {
    if (*grid_valuep < minimum) minimum = *grid_valuep;
    if (*grid_valuep > maximum) maximum = *grid_valuep;
    grid_valuep++;
  }
  return RNInterval(minimum, maximum);
}




R3Point R3CharGrid::
GridCentroid(void) const
{
  // Compute weighted sum
  RNScalar total_value = 0;
  R3Point centroid(0,0,0);
  unsigned char *grid_valuesp = grid_values;
  for (unsigned long k = 0; k < grid_resolution[2]; k++) {
    for (unsigned long j = 0; j < grid_resolution[1]; j++) {
      for (unsigned long i = 0; i < grid_resolution[0]; i++) {
        R3Vector position(i, j, k);
        RNScalar value = *(grid_valuesp++);
        centroid += value * position;
        total_value += value;
      }
    }
  }

  // Divide by total value
  if (total_value > 0) centroid /= total_value;

  // Return centroid
  return centroid;
}



R2Grid *R3CharGrid::
Slice(int dim, int grid_coordinate) const
{
  // Extract 2D grid along slice at given coordinate in given dimension

  // Compute parameters
  int dim1 = (dim+1)%3;
  int dim2 = (dim+2)%3;
  R3Box world_box = WorldBox();
  R2Box slice_box(world_box[0][dim1], world_box[0][dim2], world_box[1][dim1], world_box[1][dim2]);

  // Allocate 2D grid for slice
  R2Grid *slice = new R2Grid(grid_resolution[dim1], grid_resolution[dim2], slice_box);
  if (!slice) {
    fprintf(stderr, "Unable to allocate slice\n");
    return NULL;
  }

  // Fill values
  for (unsigned long i = 0; i < grid_resolution[dim1]; i++) {
    for (unsigned long j = 0; j < grid_resolution[dim2]; j++) {
      RNScalar value = 0.0;
      if (dim == RN_X) value = GridValue(grid_coordinate, i, j);
      else if (dim == RN_Y) value = GridValue(j, grid_coordinate, i);
      else value = GridValue(i, j, grid_coordinate);
      slice->SetGridValue(i, j, value);
    }
  }

  // Return slice
  return slice;
}



void R3CharGrid::
SetWorldToGridTransformation(const R3Affine& affine)
{
  // Set transformations
  world_to_grid_transform = affine;
  grid_to_world_transform = affine.Inverse();
  world_to_grid_scale_factor = affine.ScaleFactor();
  grid_to_world_scale_factor = (world_to_grid_scale_factor != 0) ? 1 / world_to_grid_scale_factor : 1.0;
}



void R3CharGrid::
SetWorldToGridTransformation(const R3Box& world_box)
{
  // Just checking
  if (grid_size == 0) return;
  if (world_box.NDimensions() < 3) return;

  // Compute grid origin
  R3Vector grid_diagonal(XResolution()-1, YResolution()-1, ZResolution()-1);
  R3Vector grid_origin = 0.5 * grid_diagonal;

  // Compute world origin
  R3Vector world_diagonal(world_box.XLength(), world_box.YLength(), world_box.ZLength());
  R3Vector world_origin = world_box.Centroid().Vector();

  // Compute scale
  RNScalar scale = FLT_MAX;
  RNScalar xscale = (world_diagonal[0] > 0) ? grid_diagonal[0] / world_diagonal[0] : FLT_MAX;
  if (xscale < scale) scale = xscale;
  RNScalar yscale = (world_diagonal[1] > 0) ? grid_diagonal[1] / world_diagonal[1] : FLT_MAX;
  if (yscale < scale) scale = yscale;
  RNScalar zscale = (world_diagonal[2] > 0) ? grid_diagonal[2] / world_diagonal[2] : FLT_MAX;
  if (zscale < scale) scale = zscale;
  if (scale == FLT_MAX) scale = 1;

  // Compute world-to-grid transformation
  R3Affine affine(R3identity_affine);
  affine.Translate(grid_origin);
  if (scale != 1) affine.Scale(scale);
  affine.Translate(-world_origin);

  // Set transformations
  SetWorldToGridTransformation(affine);
}



void R3CharGrid::
SetWorldToGridTransformation(const R3Point& world_origin, const R3Vector& world_axis1, const R3Vector& world_axis2, RNLength world_radius)
{
  // Just checking
  if (grid_size == 0) return;

  // Compute grid origin
  R3Vector grid_diagonal(XResolution()-1, YResolution()-1, ZResolution()-1);
  R3Vector grid_origin = 0.5 * grid_diagonal;
  RNScalar grid_radius = grid_origin[0];
  if (grid_origin[1] < grid_radius) grid_radius = grid_origin[1];
  if (grid_origin[2] < grid_radius) grid_radius = grid_origin[2];

  // Compute scale
  if (RNIsZero(world_radius)) return;
  if (RNIsZero(grid_radius)) return;
  RNScalar scale = grid_radius / world_radius;

  // Compute rotation
  R3Triad world_triad(world_axis1, world_axis2);
  R3Affine rotation(world_triad.InverseMatrix());

  // Compute world-to-grid transformation
  R3Affine affine(R3identity_affine);
  affine.Translate(grid_origin);
  affine.Transform(rotation);
  affine.Scale(scale);
  affine.Translate(-(world_origin.Vector()));

  // Set transformations
  SetWorldToGridTransformation(affine);
}



RNScalar R3CharGrid::
WorldSpacing(RNDimension dim) const
{
  // Return distance between grid samples in dimension
  assert((0 <= dim) && (dim <= 2));
  RNScalar m0 = grid_to_world_transform.Matrix()[0][dim];
  RNScalar m1 = grid_to_world_transform.Matrix()[1][dim];
  RNScalar m2 = grid_to_world_transform.Matrix()[2][dim];
  RNScalar spacing_squared = m0*m0 + m1*m1 + m2*m2;
  return sqrt(spacing_squared);
}



R3Point R3CharGrid::
WorldPosition(RNCoord x, RNCoord y, RNCoord z) const
{
  // Transform point from grid coordinates to world coordinates
  R3Point world_point(x, y, z);
  world_point.Transform(grid_to_world_transform);
  return world_point;

}



R3Point R3CharGrid::
GridPosition(RNCoord x, RNCoord y, RNCoord z) const
{
  // Transform point from world coordinates to grid coordinates
  R3Point grid_point(x, y, z);
  grid_point.Transform(world_to_grid_transform);
  return grid_point;
}



void R3CharGrid::
DrawSlice(RNDimension dim, unsigned long coord) const
{
    printf("A\n");
    printf("%d %d\n", dim, coord);
  // Check coordinates
  if ((dim < 0) || (dim > 2)) return;
  if ((coord < 0) || (coord >= Resolution(dim))) return;
    printf("B\n");
  // Get useful variables
  RNDimension dim1 = (dim+1) % 3;
 RNDimension dim2 = (dim+2) % 3;
  int width = Resolution(dim1);
  int height = Resolution(dim2);
  const int max_resolution = 1024;
  if (width > max_resolution) width = max_resolution;
  if (height > max_resolution) height = max_resolution;
    printf("C\n");
  // Define slice texture
  static const R3CharGrid *previous_grid[3] = { NULL, NULL, NULL };
  static unsigned long previous_coord[3] = { -1, -1, -1 };
  static GLuint texture_id[3] = { 0, 0, 0 };
  if ((this != previous_grid[dim]) || (coord != previous_coord[dim])) {
    if (texture_id[dim] != 0) glDeleteTextures(1, &texture_id[dim]);
    glGenTextures(1, &texture_id[dim]);
    glBindTexture(GL_TEXTURE_2D, texture_id[dim]);
    glTexEnvi(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_MODULATE);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    static GLfloat pixels[max_resolution * max_resolution];
    RNInterval range = Range();
    printf("%lf %lf\n", range.Min(), range.Max());
    if (range.Diameter() <= 0) range = RNunit_interval;
    GLfloat *pixelsp = pixels;
    printf("%d %d\n", height, width);
    for (int j = 0; j < height; j++) {
      for (int i = 0; i < width; i++) {
        RNCoord p[3];
        p[dim] = coord;
        p[dim1] = (i+0.5) * Resolution(dim1) / width;
        p[dim2] = (j+0.5) * Resolution(dim2) / height;
        RNScalar value = GridValue(p[0], p[1], p[2]);
        *(pixelsp++) = (GLfloat) ((value - range.Min()) / range.Diameter());
      }
    }
    glTexImage2D(GL_TEXTURE_2D, 0, 1, width, height, 0, GL_LUMINANCE, GL_FLOAT, pixels);
    previous_coord[dim] = coord;
    previous_grid[dim] = this;
  }

  // Set OpenGL modes
  assert(texture_id[dim]);
  glBindTexture(GL_TEXTURE_2D, texture_id[dim]);
  glEnable(GL_TEXTURE_2D);

  // Create quad
  R3Point p0, p1, p2, p3;
  p0[dim] = coord; p0[dim1] = 0;                  p0[dim2] = 0;
  p1[dim] = coord; p1[dim1] = Resolution(dim1)-1; p1[dim2] = 0;
  p2[dim] = coord; p2[dim1] = Resolution(dim1)-1; p2[dim2] = Resolution(dim2)-1;
  p3[dim] = coord; p3[dim1] = 0;                  p3[dim2] = Resolution(dim2)-1;

  // Draw quad 
  RNScalar one = 1.0;
  RNScalar zero = 0;
  R3BeginPolygon();
  R3LoadTextureCoords(zero, zero);
  R3LoadPoint(p0[0], p0[1], p0[2]);
  R3LoadTextureCoords(one, zero);
  R3LoadPoint(p1[0], p1[1], p1[2]);
  R3LoadTextureCoords(one, one);
  R3LoadPoint(p2[0], p2[1], p2[2]);
  R3LoadTextureCoords(zero, one);
  R3LoadPoint(p3[0], p3[1], p3[2]);
  R3EndPolygon();

  // Reset OpenGL modes
  glDisable(GL_TEXTURE_2D);
}