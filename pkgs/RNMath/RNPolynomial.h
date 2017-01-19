// Include file for polynomial classes



////////////////////////////////////////////////////////////////////////
// Just making sure
////////////////////////////////////////////////////////////////////////

#ifndef __RN_POLYNOMIAL__
#define __RN_POLYNOMIAL__



////////////////////////////////////////////////////////////////////////
// Select dependencies (can be set in compile flags by app)
////////////////////////////////////////////////////////////////////////

// #define RN_USE_SPLM
// #define RN_NO_SPLM
// #define RN_USE_MINPACK
// #define RN_NO_MINPACK
// #define RN_USE_CERES
// #define RN_NO_CERES
// #define RN_USE_CSPARSE
// #define RN_NO_CSPARSE



////////////////////////////////////////////////////////////////////////
// Predeclarations
////////////////////////////////////////////////////////////////////////

class RNPolynomial;
class RNPolynomialTerm;
class RNPolynomialSystemOfEquations;



////////////////////////////////////////////////////////////////////////
// System of polynomial equations
////////////////////////////////////////////////////////////////////////

class RNPolynomialSystemOfEquations {
public:
  // Constructor/destructor
  RNPolynomialSystemOfEquations(int nvariables = 0);
  RNPolynomialSystemOfEquations(const RNPolynomialSystemOfEquations& system);
  ~RNPolynomialSystemOfEquations(void);

  // Property functions
  RNScalar Degree(void) const;
  int NVariables(void) const;
  int NPolynomials(void) const;
  int NTerms(void) const;
  int NPartialDerivatives(void) const;

  // Access functions
  RNPolynomial *Polynomial(int k) const;

  // Manipulation functions
  void InsertPolynomial(RNPolynomial *polynomial);
  void RemovePolynomial(RNPolynomial *polynomial);

  // Evaluation functions
  void Evaluate(const RNScalar *x, RNScalar *y) const;

  // Optimization functions
  int Minimize(RNScalar *x, int solver = 0) const;

  // Print functions
  void Print(FILE *fp = stdout) const;

private:
  RNArray<RNPolynomial *> polynomials;
  int nvariables;
};



////////////////////////////////////////////////////////////////////////
// Polynomial
////////////////////////////////////////////////////////////////////////

class RNPolynomial {
public:
  // Constructor/destructor
  RNPolynomial(void);
  RNPolynomial(const RNPolynomial& polynomial);
  ~RNPolynomial(void);

  // Property functions
  RNScalar Degree(void) const;
  int NTerms(void) const;
  int NPartialDerivatives(void) const;
  RNBoolean IsConstant(void) const;
  RNBoolean IsLinear(void) const;
  RNBoolean IsQuadratic(void) const;

  // Access functions
  RNPolynomialSystemOfEquations *SystemOfEquations(void) const;
  RNPolynomialTerm *Term(int k) const;

  // Manipulation functions
  void Empty(void);
  void Negate(void);
  void Add(RNScalar constant);
  void Subtract(RNScalar constant);
  void Multiply(RNScalar factor);
  void Divide(RNScalar factor);
  void Add(const RNPolynomial& polynomial);
  void Subtract(const RNPolynomial& polynomial);
  void Multiply(const RNPolynomial& polynomial);

  // Assignment operators
  RNPolynomial& operator=(const RNPolynomial& polynomial);
  RNPolynomial& operator+=(const RNPolynomial& polynomial);
  RNPolynomial& operator-=(const RNPolynomial& polynomial);
  RNPolynomial& operator*=(const RNPolynomial& polynomial);
  RNPolynomial& operator+=(RNScalar a);
  RNPolynomial& operator-=(RNScalar a);
  RNPolynomial& operator*=(RNScalar a);
  RNPolynomial& operator/=(RNScalar a);
  
  // Arithmetic operators
  friend RNPolynomial operator-(const RNPolynomial& polynomial);
  friend RNPolynomial operator+(const RNPolynomial& polynomial1, const RNPolynomial& polynomial2);
  friend RNPolynomial operator+(const RNPolynomial& polynomial, RNScalar a);
  friend RNPolynomial operator+(RNScalar a, const RNPolynomial& polynomial);
  friend RNPolynomial operator-(const RNPolynomial& polynomial1, const RNPolynomial& polynomial2);
  friend RNPolynomial operator-(const RNPolynomial& polynomial, RNScalar a);
  friend RNPolynomial operator-(RNScalar a, const RNPolynomial& polynomial);
  friend RNPolynomial operator*(const RNPolynomial& polynomial1, const RNPolynomial& polynomial2);
  friend RNPolynomial operator*(const RNPolynomial& polynomial, RNScalar a);
  friend RNPolynomial operator*(RNScalar a, const RNPolynomial& polynomial);
  friend RNPolynomial operator/(const RNPolynomial& polynomial, RNScalar a);

  // Construction functions
  void AddTerm(RNScalar c, RNBoolean already_unique = FALSE);
  void AddTerm(RNScalar c, int v1, RNScalar e1, RNBoolean already_unique = FALSE);
  void AddTerm(RNScalar c, int n, const int *v, const RNScalar *e,
    RNBoolean already_sorted = FALSE, RNBoolean already_unique = FALSE);

  // Evaluation functions
  RNScalar Evaluate(const RNScalar *x) const;

  // Print functions
  void Print(FILE *fp = stdout) const;

public:
  // Internal functions
  RNPolynomialTerm *FindTermWithSameVariables(const RNPolynomialTerm *query) const;
  RNPolynomialTerm *FindTermWithVariables(int n, int *v, RNScalar *e) const;

private:
  friend class RNPolynomialSystemOfEquations;
  RNPolynomialSystemOfEquations *system;
  int system_index;
  RNArray<RNPolynomialTerm *> terms;
};



////////////////////////////////////////////////////////////////////////
// Polynomial term
////////////////////////////////////////////////////////////////////////

class RNPolynomialTerm {
public:
  // Constructor/destructor
  RNPolynomialTerm(RNScalar c = 0.0, int nv = 0, const int *v = NULL, const RNScalar *e = NULL, 
    RNBoolean already_sorted = FALSE, RNBoolean already_unique = FALSE);
  RNPolynomialTerm(const RNPolynomialTerm& term);
  ~RNPolynomialTerm(void);

  // Property functions
  RNScalar Degree(void) const;
  int NVariables(void) const;
  int NPartialDerivatives(void) const;
  RNBoolean IsConstant(void) const;
  RNBoolean IsLinear(void) const;
  RNBoolean IsQuadratic(void) const;

  // Access functions
  RNScalar Coefficient(void) const;
  int Variable(int k) const;
  RNScalar Exponent(int k) const;
  const int *Variables(void) const;
  const RNScalar *Exponents(void) const;
  RNPolynomial *Polynomial(void) const;

  // Manipulation functions
  void Empty(void);
  void Negate(void);
  void Multiply(RNScalar factor);
  void Divide(RNScalar factor);
  void SetCoefficient(RNScalar c);
  void SetVariable(int k, int v);
  void SetExponent(int k, RNScalar e);

  // Evaluation functions
  RNScalar PartialDerivative(const RNScalar *x, int variable) const;
  RNScalar Evaluate(const RNScalar *x) const;

  // Print functions
  void Print(FILE *fp = stdout) const;

public:
  // Internal functions
  RNBoolean HasSameVariables(const RNPolynomialTerm *query) const;
  RNBoolean HasVariables(int n, const int *v, const RNScalar *e) const;

private:
  friend class RNPolynomial;
  RNPolynomial *polynomial;
  int n;
  RNScalar c;
  int *v;
  RNScalar *e;
};



////////////////////////////////////////////////////////////////////////
// Polynomial system of equation solvers
////////////////////////////////////////////////////////////////////////

enum {
  RN_POLYNOMIAL_CERES_SOLVER,
  RN_POLYNOMIAL_MINPACK_SOLVER,
  RN_POLYNOMIAL_SPLM_SOLVER,
  RN_POLYNOMIAL_CSPARSE_SOLVER,
  RN_POLYNOMIAL_NUM_SOLVERS
};



////////////////////////////////////////////////////////////////////////
// Inline functions for polynomial system of equations
////////////////////////////////////////////////////////////////////////

inline int RNPolynomialSystemOfEquations::
NVariables(void) const
{
  // Return number of variables
  return nvariables;
}



inline int RNPolynomialSystemOfEquations::
NPolynomials(void) const
{
  // Return number of polynomials
  return polynomials.NEntries();
}



inline RNPolynomial *RNPolynomialSystemOfEquations::
Polynomial(int k) const
{
  // Return Kth polynomial
  return polynomials.Kth(k);
}



////////////////////////////////////////////////////////////////////////
// Inline functions for polynomial 
////////////////////////////////////////////////////////////////////////

inline int RNPolynomial::
NTerms(void) const
{
  // Return number of terms
  return terms.NEntries();
}



inline RNPolynomialTerm *RNPolynomial::
Term(int k) const
{
  // Return kth term
  return terms.Kth(k);
}



inline RNPolynomialSystemOfEquations *RNPolynomial::
SystemOfEquations(void) const
{
  // Return system of equations this polynomial belongs to
  return system;
}



inline void RNPolynomial::
AddTerm(RNScalar c, RNBoolean already_unique)
{
  // Add term
  AddTerm(c, 0, NULL, NULL, TRUE, already_unique);
}



inline void RNPolynomial::
AddTerm(RNScalar c, int v1, RNScalar e1, RNBoolean already_unique)
{
  // Add term
  AddTerm(c, 1, &v1, &e1, TRUE, already_unique);
}



inline void RNPolynomial::
Add(RNScalar constant)
{
  // Add constant 
  AddTerm(constant);
}



inline void RNPolynomial::
Subtract(RNScalar constant)
{
  // Subtract constant 
  Add(-constant);
}



inline RNPolynomial& RNPolynomial::
operator+=(const RNPolynomial& polynomial)
{
  // Add polynomial
  Add(polynomial);
  return *this;
}



inline RNPolynomial& RNPolynomial::
operator-=(const RNPolynomial& polynomial)
{
  // Subtract polynomial
  Subtract(polynomial);
  return *this;
}



inline RNPolynomial& RNPolynomial::
operator*=(const RNPolynomial& polynomial)
{
  // Multiply polynomial
  Multiply(polynomial);
  return *this;
}



inline RNPolynomial& RNPolynomial::
operator+=(RNScalar a)
{
  // Add constant
  Add(a);
  return *this;
}



inline RNPolynomial& RNPolynomial::
operator-=(RNScalar a)
{
  // Subtract constant
  Subtract(a);
  return *this;
}



inline RNPolynomial& RNPolynomial::
operator*=(RNScalar a)
{
  // Multiply by constant
  Multiply(a);
  return *this;
}



inline RNPolynomial& RNPolynomial::
operator/=(RNScalar a)
{
  // Divide by constant
  Divide(a);
  return *this;
}



inline RNPolynomial
operator-(const RNPolynomial& polynomial)
{
  // Return negated polynomial
  RNPolynomial result(polynomial);
  result.Multiply(-1.0);
  return result;
}



inline RNPolynomial 
operator+(const RNPolynomial& polynomial1, const RNPolynomial& polynomial2)
{
  // Return sum
  RNPolynomial result(polynomial1);
  result.Add(polynomial2);
  return result;
}



inline RNPolynomial 
operator+(const RNPolynomial& polynomial, RNScalar constant)
{
  // Return sum
  RNPolynomial result(polynomial);
  result.Add(constant);
  return result;
}



inline RNPolynomial 
operator+(RNScalar constant, const RNPolynomial& polynomial)
{
  // Return sum
  RNPolynomial result(polynomial);
  result.Add(constant);
  return result;
}



inline RNPolynomial 
operator-(const RNPolynomial& polynomial1, const RNPolynomial& polynomial2)
{
  // Return difference
  RNPolynomial result(polynomial1);
  result.Subtract(polynomial2);
  return result;
}



inline RNPolynomial 
operator-(const RNPolynomial& polynomial, RNScalar constant)
{
  // Return difference
  RNPolynomial result(polynomial);
  result.Subtract(constant);
  return result;
}



inline RNPolynomial 
operator-(RNScalar constant, const RNPolynomial& polynomial)
{
  // Return difference
  RNPolynomial result(polynomial);
  result.Negate();
  result.Add(constant);
  return result;
}



inline RNPolynomial 
operator*(const RNPolynomial& polynomial1, const RNPolynomial& polynomial2)
{
  // Return product
  RNPolynomial result(polynomial1);
  result.Multiply(polynomial2);
  return result;
}



inline RNPolynomial 
operator*(const RNPolynomial& polynomial, RNScalar a)
{
  // Return product
  RNPolynomial result(polynomial);
  result.Multiply(a);
  return result;
}



inline RNPolynomial 
operator*(RNScalar a, const RNPolynomial& polynomial)
{
  // Return product
  RNPolynomial result(polynomial);
  result.Multiply(a);
  return result;
}



inline RNPolynomial 
operator/(const RNPolynomial& polynomial, RNScalar a)
{
  // Return quotient
  RNPolynomial result(polynomial);
  result.Divide(a);
  return result;
}



////////////////////////////////////////////////////////////////////////
// Inline functions for polynomial term
////////////////////////////////////////////////////////////////////////

inline int RNPolynomialTerm::
NVariables(void) const
{
  // Return the number of variables
  return n;
}



inline int RNPolynomialTerm::
NPartialDerivatives(void) const
{
  // Return the number of partial derivatives
  return n;
}



inline int RNPolynomialTerm::
Variable(int k) const
{
  // Return the index of the kth variable
  return v[k];
}



inline RNScalar RNPolynomialTerm::
Exponent(int k) const
{
  // Return the exponent of the kth variable
  return e[k];
}



inline RNScalar RNPolynomialTerm::
Coefficient(void) const
{
  // Return the coefficient
  return c;
}



inline const int *RNPolynomialTerm::
Variables(void) const
{
  // Return the variables
  return v;
}



inline const RNScalar *RNPolynomialTerm::
Exponents(void) const
{
  // Return the exponents
  return e;
}



inline RNPolynomial *RNPolynomialTerm::
Polynomial(void) const
{
  // Return the polynomial this term is part of
  return polynomial;
}



inline void RNPolynomialTerm::
Negate(void)
{
  // Negate coefficient
  c = -c;
}



inline void RNPolynomialTerm::
Multiply(RNScalar factor) 
{
  // Multiply term by a constant factor
  c *= factor;
}



inline void RNPolynomialTerm::
Divide(RNScalar factor) 
{
  // Multiply term by a constant factor
  if (RNIsZero(factor)) return;
  c /= factor;
}



inline void RNPolynomialTerm::
SetCoefficient(RNScalar c)
{
  // Set coefficient
  this->c = c;
}



inline void RNPolynomialTerm::
SetVariable(int k, int v)
{
  // Set index of kth variable
  this->v[k] = v;
}



inline void RNPolynomialTerm::
SetExponent(int k, RNScalar e)
{
  // Set exponent of kth variable
  this->e[k] = e;
}



inline RNBoolean RNPolynomialTerm::
HasSameVariables(const RNPolynomialTerm *query) const
{
  // Check if query has same variables
  return HasVariables(query->NVariables(), query->Variables(), query->Exponents());
}


////////////////////////////////////////////////////////////////////////
// CSPARSE Stuff
////////////////////////////////////////////////////////////////////////

#ifdef RN_NO_CSPARSE
#undef RN_USE_CSPARSE
#endif
#ifdef RN_USE_CSPARSE

#include "CSparse/CSparse.h"

static int 
MinimizeCSPARSE(const RNPolynomialSystemOfEquations *system, RNScalar *io)
{
  // Get convenient variables
  const int n = system->NVariables();
  const int mm = system->NPolynomials();
  const int nz = system->NTerms();

  // Allocate matrix
  cs *a = cs_spalloc (0, n, nz, 1, 1);
  if (!a) {
    fprintf(stderr, "Unable to allocate cs matrix: %d %d\n", n, nz);
    return 0;
  }
    
  // Allocate B vector
  double *b = new double [ mm ];
  for (int i = 0; i < mm; i++) b[i] = 0;

  // Allocate X vector
  double *x = new double [ n ];
  for (int i = 0; i < n; i++) x[i] = 0;

  // Fill matrix
  int m =0;
  for (int i = 0; i < system->NPolynomials(); i++) {
    RNPolynomial *polynomial = system->Polynomial(i);

    // Count non-zero terms
    int nz = 0;
    for (int j = 0; j < polynomial->NTerms(); j++) {
      RNPolynomialTerm *term = polynomial->Term(j);
      RNScalar coefficient = term->Coefficient();
      for (int k = 1; k < term->NVariables(); k++) coefficient *= io[term->Variable(k)];
      if (coefficient != 0) nz++;
    }

    // Check if there are nonzero terms
    if (nz == 0) continue;

    // Add row to matrix if there are non-zero terms
    for (int j = 0; j < polynomial->NTerms(); j++) {
      RNPolynomialTerm *term = polynomial->Term(j);
      RNScalar coefficient = term->Coefficient();
      if (coefficient == 0.0) continue;
      if (term->NVariables() == 0) {
        b[m] -= coefficient;
      }
      else if (term->NVariables() == 1) {
        // This is approximate if exponent!=1
        int variable = term->Variable(0);
        assert((variable >= 0) && (variable < n));
        cs_entry(a, m, variable, coefficient);
      }
      else {
        // This is approximate 

        // Find variable with largest partial derivative
        int best_variable = -1;
        RNScalar best_partial_derivative = 0;
        for (int k = 0; k < term->NVariables(); k++) {
          assert((term->Variable(k) >= 0) && (term->Variable(k) < n));
          RNScalar partial_derivative = fabs(term->PartialDerivative(io, k));
          if (partial_derivative > best_partial_derivative) {
            best_partial_derivative = partial_derivative;
            best_variable = k;
          }
        }
        
        // Add term for variable with largest partial derivative
        if (best_partial_derivative != 0) {
          // Compute factor
          RNScalar factor = coefficient;
          for (int k = 0; k < term->NVariables(); k++) {
            if (k == best_variable) continue;
            int v = term->Variable(k);
            RNScalar e = term->Exponent(k);
            if (e == 1.0) factor *= io[v];
            else if (e == 2.0) factor *= io[v] * io[v];
            else if ((e < 0.0) && RNIsZero(io[v])) factor *= RN_INFINITY;
            else factor *= pow(io[v], e);
          }

          // Add term to matrix
          if (factor != 0) {
            cs_entry(a, m, best_variable, factor);
          }
        }
      }
    }

    // Increment row of matrix
    assert(m < mm);
    m++;
  }

  // Just checking
  assert(a->m == m);
  assert(a->n == n);
  assert(a->n == system->NVariables());
  assert(a->m <= system->NPolynomials());
  assert(a->nz <= system->NPartialDerivatives());

  // Setup aT * a * x = aT * b        
  cs *A = cs_compress(a);
  assert(A);
  cs *AT = cs_transpose (A, 1);
  assert(AT);
  cs *ATA = cs_multiply (AT, A);
  assert(ATA);
  cs_gaxpy(AT, b, x);

  // Solve linear system
  int status = cs_lusol (1, ATA, x, RN_EPSILON);
  // int status = cs_cholsol (1, ATA, x);
  if (status == 0) fprintf(stderr, "Error in CSPARSE solver\n");
  else { for (int i = 0; i < n; i++) io[i] = x[i]; }

  // Delete stuff
  cs_spfree(A);
  cs_spfree(AT);
  cs_spfree(ATA);
  cs_spfree(a);
  delete [] b;
  delete [] x;

  // Return status
  return status;
}

#else

static int 
MinimizeCSPARSE(const RNPolynomialSystemOfEquations *system, RNScalar *io)
{
  // Print error message
  fprintf(stderr, "Cannot minimize polynomial: CSparse solver disabled during compile.\n");
  fprintf(stderr, "Enable it by adding -DRN_USE_CSPARSE and -lCSparse to compilation and link commands.\n");
  return 0;
}

#endif




////////////////////////////////////////////////////////////////////////
// Minimization function
////////////////////////////////////////////////////////////////////////

inline int RNPolynomialSystemOfEquations::
Minimize(RNScalar *x, int solver) const
{
  // Check solver
  if (solver == RN_POLYNOMIAL_CSPARSE_SOLVER) return MinimizeCSPARSE(this, x);
  fprintf(stderr, "Polynomial solver not recognized: %d\n", solver);
  return 0;
}



#endif
