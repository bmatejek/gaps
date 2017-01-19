// Include file for abstract linear system class



// Class definition

class RNLinearSystem {
public:
  // Constructor/destructor
  RNLinearSystem(int m, int n);
  RNLinearSystem(const RNMatrix& A, const RNVector& B);
  RNLinearSystem(const RNLinearSystem& system);
  ~RNLinearSystem(void);

  // Property functions/operators
  virtual int NEquations(void) const = 0;
  virtual int NVariables(void) const = 0;
  virtual RNBoolean IsDense(void) const = 0;
  virtual RNBoolean IsSparse(void) const = 0;
  virtual RNLinearSystem *SparseCopy(void) const = 0;
  virtual RNLinearSystem *DenseCopy(void) const = 0;
  virtual RNBoolean operator==(const RNLinearSystem& matrix) const = 0;
  virtual RNBoolean operator!=(const RNLinearSystem& matrix) const = 0;

  // Entry access
  virtual RNScalar A(int i, int j) const = 0;
  virtual RNScalar B(int i) const = 0;

  // Manipulation
  virtual SetA(int i, int j, RNScalar a) = 0;
  virtual SetB(int i, RNScalar b) = 0;

  // Solution methods
  int Solve(RNVector& X) = 0;
};




