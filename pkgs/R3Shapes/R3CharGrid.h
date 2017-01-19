// Header file for GAPS scalar grid class



// Class definition

class R3CharGrid {
public:
  // Constructors
  R3CharGrid(unsigned long xresolution = 0, unsigned long yresolution = 0, unsigned long zresolution = 0);
  R3CharGrid(unsigned long xresolution, unsigned long yresolution, unsigned long zresolution, const R3Box& bbox);
  R3CharGrid(const R3CharGrid& grid);
  ~R3CharGrid(void);

  // Grid property functions
  unsigned long NEntries(void) const;
  unsigned long XResolution(void) const;
  unsigned long YResolution(void) const;
  unsigned long ZResolution(void) const;
  unsigned long Resolution(RNDimension dim) const;
  unsigned char Maximum(void) const;
  unsigned char Minimum(void) const;
  
  
  RNInterval Range(void) const;
  R3Box GridBox(void) const;
  R3Box WorldBox(void) const;
  R3Point GridCentroid(void) const;
  R3Point WorldCentroid(void) const;
  R2Grid *Slice(int dim, int grid_coordinate) const;


  // Transformation property functions
  const R3Affine& WorldToGridTransformation(void) const;
  const R3Affine& GridToWorldTransformation(void) const;
  RNScalar WorldToGridScaleFactor(void) const;
  RNScalar GridToWorldScaleFactor(void) const;
  RNScalar WorldSpacing(RNDimension dim) const;

  // Grid value access functions
  unsigned char GridValue(unsigned long index) const;
  unsigned char GridValue(unsigned long i, unsigned long j, unsigned long k) const;
  unsigned char WorldValue(RNCoord x, RNCoord y, RNCoord z) const;
  unsigned char WorldValue(const R3Point& world_point) const;
  unsigned char& operator()(unsigned long i, unsigned long j,unsigned long k);

  // Grid manipulation functions
  void SetGridValue(unsigned long index, unsigned char value);
  void SetGridValue(unsigned long i, unsigned long j, unsigned long k, unsigned char value);

  // Transformation manipulation functions
  void SetWorldToGridTransformation(const R3Affine& affine);
  void SetWorldToGridTransformation(const R3Box& world_box);
  void SetWorldToGridTransformation(const R3Point& world_origin, const R3Vector& world_axis1, const R3Vector& world_axis2, RNLength world_radius);

  // Transformation utility functions
  R3Point WorldPosition(const R3Point& grid_point) const;
  R3Point GridPosition(const R3Point& world_point) const;
  R3Point WorldPosition(RNCoord x, RNCoord y, RNCoord z) const;
  R3Point GridPosition(RNCoord x, RNCoord y, RNCoord z) const;

  // Visualization functions
  void DrawSlice(int dim, unsigned long coord) const;

  // Debugging functions
  const unsigned char *GridValues(void) const;
  void IndicesToIndex(int i, int j, int k, int& index);
  void IndexToIndices(int index, int& i, int& j, int& k);

private:
  R3Affine grid_to_world_transform;
  R3Affine world_to_grid_transform;
  RNScalar world_to_grid_scale_factor;
  RNScalar grid_to_world_scale_factor;
  unsigned char *grid_values;
  unsigned long grid_resolution[3];
  unsigned long grid_row_size;
  unsigned long grid_sheet_size;
  unsigned long grid_size;
};



// Useful constants

extern const float R3_GRID_KEEP_VALUE;



// Inline functions

inline unsigned long R3CharGrid::
NEntries(void) const
{
  // Return total number of entries
  return grid_size;
}



inline unsigned long R3CharGrid::
Resolution(RNDimension dim) const
{
  // Return resolution in dimension
  assert((0 <= dim) && (dim <= 2));
  return grid_resolution[dim];
}




inline unsigned long R3CharGrid::
XResolution(void) const
{
  // Return resolution in X dimension
  return grid_resolution[RN_X];
}




inline unsigned long R3CharGrid::
YResolution(void) const
{
  // Return resolution in Y dimension
  return grid_resolution[RN_Y];
}




inline unsigned long R3CharGrid::
ZResolution(void) const
{
  // Return resolution in Z dimension
  return grid_resolution[RN_Z];
}



inline unsigned char R3CharGrid::
Minimum(void) const
{
  // Return smallest value
  return (unsigned char) (Range().Min() + 0.5);
}



inline unsigned char R3CharGrid::
Maximum(void) const
{
  // Return largest value
  return (unsigned char) (Range().Max() + 0.5);
}



inline R3Box R3CharGrid::
GridBox(void) const
{
  // Return bounding box in grid coordinates
  return R3Box(0, 0, 0, grid_resolution[0]-1, grid_resolution[1]-1, grid_resolution[2]-1);
}



inline R3Box R3CharGrid::
WorldBox(void) const
{
  // Return bounding box in world coordinates
  R3Point p1(0, 0, 0);
  R3Point p2(grid_resolution[0]-1, grid_resolution[1]-1, grid_resolution[2]-1);
  return R3Box(WorldPosition(p1), WorldPosition(p2));
}



inline R3Point R3CharGrid::
WorldCentroid(void) const
{
  // Return centroid in world coordinates
  R3Point grid_centroid = GridCentroid();
  return WorldPosition(grid_centroid);
}



inline const R3Affine& R3CharGrid::
WorldToGridTransformation(void) const
{
  // Return transformation from world coordinates to grid coordinates
  return world_to_grid_transform;
}



inline const R3Affine& R3CharGrid::
GridToWorldTransformation(void) const
{
  // Return transformation from grid coordinates to world coordinates
  return grid_to_world_transform;
}



inline RNScalar R3CharGrid::
WorldToGridScaleFactor(void) const
{
  // Return transformation from world coordinates to grid coordinates
  return world_to_grid_scale_factor;
}



inline RNScalar R3CharGrid::
GridToWorldScaleFactor(void) const
{
  // Return transformation from world coordinates to grid coordinates
  return grid_to_world_scale_factor;
}



inline const unsigned char *R3CharGrid::
GridValues(void) const
{
  // Return pointer to grid values
  return grid_values;
}



inline unsigned char R3CharGrid::
GridValue(unsigned long index) const
{
  // Return value at grid point referenced by index
  assert((0 <= index) && (index < grid_size));
  return grid_values[index];
}




inline unsigned char R3CharGrid::
GridValue(unsigned long i, unsigned long j, unsigned long k) const
{
  // Return value at grid point
  assert((0 <= i) && (i < XResolution()));
  assert((0 <= j) && (j < YResolution()));
  assert((0 <= k) && (k < ZResolution()));
  return grid_values[k * grid_sheet_size + j * grid_row_size + i];
}
inline unsigned char& R3CharGrid::
operator()(unsigned long i, unsigned long j, unsigned long k) 
{
  // Return value at grid point
  assert((0 <= i) && (i < XResolution()));
  assert((0 <= j) && (j < YResolution()));
  assert((0 <= k) && (k < ZResolution()));
  return grid_values[k * grid_sheet_size + j * grid_row_size + i];
}

inline R3Point R3CharGrid::
WorldPosition(const R3Point& grid_point) const
{
  // Transform point from grid coordinates to world coordinates
  return WorldPosition(grid_point[0], grid_point[1], grid_point[2]);
}



inline R3Point R3CharGrid::
GridPosition(const R3Point& world_point) const
{
  // Transform point from world coordinates to grid coordinates
  return GridPosition(world_point[0], world_point[1], world_point[2]);
}



inline void R3CharGrid::
IndicesToIndex(int i, int j, int k, int& index)
{
  // Set index of grid value at (i, j, k) 
  index = k * grid_sheet_size + j * grid_row_size + i;
}


inline void R3CharGrid::
IndexToIndices(int index, int& i, int& j, int& k)
{
  // Set indices of grid value at index
  k = index / grid_sheet_size;
  j = (index - k * grid_sheet_size) / grid_row_size;
  i = index % grid_row_size;
}




inline void R3CharGrid::
SetGridValue(unsigned long index, unsigned char value)
{
  // Set value at grid point
  assert((0 <= index) && (index < grid_size));
  grid_values[index] = value;
}



inline void R3CharGrid::
SetGridValue(unsigned long i, unsigned long j, unsigned long k, unsigned char value)
{
  // Set value at grid point
  assert((0 <= i) && (i < XResolution()));
  assert((0 <= j) && (j < YResolution()));
  assert((0 <= k) && (k < ZResolution()));
  (*this)(i, j, k) = value;
}
