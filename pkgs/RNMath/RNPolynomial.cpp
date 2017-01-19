// Source file for polynomial classes



////////////////////////////////////////////////////////////////////////
// Include files
////////////////////////////////////////////////////////////////////////

#include "RNMath/RNMath.h"



////////////////////////////////////////////////////////////////////////
// System of polynomial equations
////////////////////////////////////////////////////////////////////////

RNPolynomialSystemOfEquations::
RNPolynomialSystemOfEquations(int nvariables)
  : polynomials(),
    nvariables(nvariables)
{
}



RNPolynomialSystemOfEquations::
RNPolynomialSystemOfEquations(const RNPolynomialSystemOfEquations& system)
  : polynomials(),
    nvariables(system.nvariables)
{
  // Copy polynomials
  for (int i = 0; i < system.NPolynomials(); i++) {
    RNPolynomial *polynomial = system.Polynomial(i);
    InsertPolynomial(new RNPolynomial(*polynomial));
  }
}



RNPolynomialSystemOfEquations::
~RNPolynomialSystemOfEquations(void)
{
  // Delete all polynomials
  while (NPolynomials() > 0) {
    RNPolynomial *polynomial = polynomials.Tail();
    RemovePolynomial(polynomial);
    delete polynomial;
  }
}



RNScalar RNPolynomialSystemOfEquations::
Degree(void) const
{
  // Return highest degree of any polynomial
  RNScalar max_degree = 0;
  for (int i = 0; i < NPolynomials(); i++) {
    RNPolynomial *polynomial = Polynomial(i);
    RNScalar degree = polynomial->Degree();
    if (degree > max_degree) max_degree = degree;
  }
  return max_degree;
}



int RNPolynomialSystemOfEquations::
NTerms(void) const
{
  // Return number of terms
  int count = 0;
  for (int i = 0; i < NPolynomials(); i++) {
    RNPolynomial *polynomial = Polynomial(i);
    count += polynomial->NTerms();
  }
  return count;
}



int RNPolynomialSystemOfEquations::
NPartialDerivatives(void) const
{
  // Return number of partial derivatives
  int count = 0;
  for (int i = 0; i < NPolynomials(); i++) {
    RNPolynomial *polynomial = Polynomial(i);
    count += polynomial->NPartialDerivatives();
  }
  return count;
}



void RNPolynomialSystemOfEquations::
InsertPolynomial(RNPolynomial *polynomial)
{
  // Just checking
  assert(!polynomial->system);
  assert(polynomial->system_index == -1);
  assert(!polynomials.FindEntry(polynomial));

  // Insert polynomial
  polynomial->system = this;
  polynomial->system_index = polynomials.NEntries();
  polynomials.Insert(polynomial);
}



void RNPolynomialSystemOfEquations::
RemovePolynomial(RNPolynomial *polynomial)
{
  // Just checking
  assert(polynomial->system == this);
  assert(polynomial->system_index >= 0);
  assert(polynomials.FindEntry(polynomial));

  // Remove polynomial
  RNArrayEntry *entry = polynomials.KthEntry(polynomial->system_index);
  assert(entry && (polynomials.EntryContents(entry) == polynomial));
  RNPolynomial *tail = polynomials.Tail();
  tail->system_index = polynomial->system_index;
  polynomials.EntryContents(entry) = tail;
  polynomials.RemoveTail();
  polynomial->system_index = -1;
  polynomial->system = NULL;
}



void RNPolynomialSystemOfEquations::
Evaluate(const RNScalar *x, RNScalar *y) const
{
  // Evaluate polynomials
  for (int i = 0; i < NPolynomials(); i++) {
    RNPolynomial *polynomial = Polynomial(i);
    y[i] = polynomial->Evaluate(x);
  }
}



void RNPolynomialSystemOfEquations::
Print(FILE *fp) const
{
  // Print polynomial
  printf("%d\n", NPolynomials());
  for (int i = 0; i < NPolynomials(); i++) {
    RNPolynomial *polynomial = Polynomial(i);
    polynomial->Print(fp);
  }
}



////////////////////////////////////////////////////////////////////////
// Polynomial
////////////////////////////////////////////////////////////////////////

RNPolynomial::
RNPolynomial(void)
  : system(NULL),
    system_index(-1),
    terms()
{
}



RNPolynomial::
RNPolynomial(const RNPolynomial& polynomial)
  : system(NULL),
    system_index(-1),
    terms()
{
  // Copy polynomial terms
  for (int i = 0; i < polynomial.NTerms(); i++) {
    RNPolynomialTerm *polynomial_term = polynomial.Term(i);
    RNPolynomialTerm *term = new RNPolynomialTerm(*polynomial_term);
    terms.Insert(term); 
    term->polynomial = this; 
  }
}



RNPolynomial::
~RNPolynomial(void)
{
  // Remove from system of equations
  if (system) system->RemovePolynomial(this);

  // Delete terms
  Empty();
}



RNScalar RNPolynomial::
Degree(void) const
{
  // Return highest degree of any term
  RNScalar max_degree = 0;
  for (int i = 0; i < NTerms(); i++) {
    RNPolynomialTerm *term = Term(i);
    RNScalar degree = term->Degree();
    if (degree > max_degree) max_degree = degree;
  }
  return max_degree;
}



int RNPolynomial::
NPartialDerivatives(void) const
{
  // Return number of partial derivatives
  int count = 0;
  for (int i = 0; i < NTerms(); i++) {
    RNPolynomialTerm *term = Term(i);
    count += term->NPartialDerivatives();
  }
  return count;
}



RNBoolean RNPolynomial::
IsConstant(void) const
{
  // Check each term
  for (int i = 0; i < NTerms(); i++) {
    RNPolynomialTerm *term = Term(i);
    if (!term->IsConstant()) return FALSE;
  }

  // Passed all tests
  return TRUE;
}



RNBoolean RNPolynomial::
IsLinear(void) const
{
  // Check each term
  for (int i = 0; i < NTerms(); i++) {
    RNPolynomialTerm *term = Term(i);
    if (!term->IsLinear()) return FALSE;
  }

  // Passed all tests
  return TRUE;
}



RNBoolean RNPolynomial::
IsQuadratic(void) const
{
  // Check each term
  for (int i = 0; i < NTerms(); i++) {
    RNPolynomialTerm *term = Term(i);
    if (!term->IsQuadratic()) return FALSE;
  }

  // Passed all tests
  return TRUE;
}



void RNPolynomial::
Empty(void)
{
  // Delete all terms
  for (int i = 0; i < NTerms(); i++) {
    RNPolynomialTerm *term = Term(i);
    delete term;
  }

  // Empty array of terms
  terms.Empty();
}



void RNPolynomial::
Negate(void)
{
  // Negate all terms
  for (int i = 0; i < NTerms(); i++) {
    RNPolynomialTerm *term = Term(i);
    term->Negate();
  }
}



void RNPolynomial::
Multiply(RNScalar factor)
{
  // Multiply all terms by factor
  for (int i = 0; i < NTerms(); i++) {
    RNPolynomialTerm *term = Term(i);
    term->Multiply(factor);
  }
}



void RNPolynomial::
Divide(RNScalar factor)
{
  // Multiply all terms by factor
  if (RNIsZero(factor)) return;
  Multiply(1.0/factor);
}



void RNPolynomial::
Add(const RNPolynomial& polynomial)
{
  // Add polynomial
  for (int i = 0; i < polynomial.NTerms(); i++) {
    RNPolynomialTerm *term = polynomial.Term(i);
    AddTerm(term->Coefficient(), term->NVariables(), term->Variables(), term->Exponents(), TRUE, TRUE);
  }
}



void RNPolynomial::
Subtract(const RNPolynomial& polynomial)
{
  // Subtract polynomial
  for (int i = 0; i < polynomial.NTerms(); i++) {
    RNPolynomialTerm *term = polynomial.Term(i);
    AddTerm(-(term->Coefficient()), term->NVariables(), term->Variables(), term->Exponents(), TRUE, TRUE);
  }
}



void RNPolynomial::
Multiply(const RNPolynomial& polynomial)
{
  // Convenient variables
  RNPolynomial polynomial1(*this);
  RNPolynomial polynomial2(polynomial);
  const int max_factors = 1024;
  int v [ max_factors ];
  RNScalar e [ max_factors ];
  int n = 0;

  // Start from scratch
  Empty();

  // For each term of this polynomial
  for (int i = 0; i < polynomial1.NTerms(); i++) {
    RNPolynomialTerm *term1 = polynomial1.Term(i);

    // Start with factors from this term
    n = 0;
    for (int k = 0; k < term1->NVariables(); k++) {
      if (n >= max_factors) break;
      v[n] = term1->Variable(k);
      e[n] = term1->Exponent(k);
      n++;
    }

    // For each term of other polynomial
    for (int j = 0; j < polynomial2.NTerms(); j++) {
      RNPolynomialTerm *term2 = polynomial2.Term(j);

      // Include factors from other term
      n = term1->NVariables();
      for (int k = 0; k < term2->NVariables(); k++) {
        if (n >= max_factors) break;
        v[n] = term2->Variable(k);
        e[n] = term2->Exponent(k);
        n++;
      }

      // Add term to this polynomial
      AddTerm(term1->Coefficient() * term2->Coefficient(), n, v, e);
    }
  }
}



RNPolynomial& RNPolynomial::
operator=(const RNPolynomial& polynomial)
{
  // Empty this polynomial
  Empty();

  // Copy polynomial terms
  for (int i = 0; i < polynomial.NTerms(); i++) {
    RNPolynomialTerm *polynomial_term = polynomial.Term(i);
    RNPolynomialTerm *term = new RNPolynomialTerm(*polynomial_term);
    terms.Insert(term); 
    term->polynomial = this; 
  }

  // Return this
  return *this;
}



void RNPolynomial::
AddTerm(RNScalar c, int n, const int *v, const RNScalar *e, 
  RNBoolean already_sorted, RNBoolean already_unique)
{
  // Check coefficient
  if (c == 0) return;

  // Check everything else
  assert((n == 0) || ((v != NULL) && (e != NULL)));

  // Add term
  RNPolynomialTerm *term = new RNPolynomialTerm(c, n, v, e, already_sorted, already_unique);
  RNPolynomialTerm *match = FindTermWithSameVariables(term);
  if (match) { 
    // Update existing term
    match->c += c;
    delete term; 
    if (RNIsZero(match->c)) {
      // Remove unneeded term
      terms.Remove(match);
      delete match;
    }
  }
  else { 
    // Insert new term
    terms.Insert(term); 
    term->polynomial = this; 
  }
}



RNScalar RNPolynomial::
Evaluate(const RNScalar *x) const
{
  // Return sum of terms
  RNScalar sum = 0;
  for (int i = 0; i < NTerms(); i++) {
    RNPolynomialTerm *term = Term(i);
    sum += term->Evaluate(x);
  }
  return sum;
}



RNPolynomialTerm *RNPolynomial::
FindTermWithSameVariables(const RNPolynomialTerm *query) const
{
  // Find matching term
  for (int i = 0; i < NTerms(); i++) {
    RNPolynomialTerm *term = Term(i);
    if (term->HasVariables(query->NVariables(), query->Variables(), query->Exponents())) return term;
  }

  // Did not find matching term
  return NULL;
}
  


RNPolynomialTerm *RNPolynomial::
FindTermWithVariables(int n, int *v, RNScalar *e) const
{
  // Find matching term
  for (int i = 0; i < NTerms(); i++) {
    RNPolynomialTerm *term = Term(i);
    if (term->HasVariables(n, v, e)) return term;
  }

  // Did not find matching term
  return NULL;
}
  


void RNPolynomial::
Print(FILE *fp) const
{
  // Print polynomial
  printf("%d\n", NTerms());
  for (int i = 0; i < NTerms(); i++) {
    RNPolynomialTerm *term = Term(i);
    fprintf(fp, "  ");
    term->Print(fp);
  }
}



////////////////////////////////////////////////////////////////////////
// Polynomial term
////////////////////////////////////////////////////////////////////////

RNPolynomialTerm::
RNPolynomialTerm(RNScalar _c, int _n, const int *_v, const RNScalar *_e, 
  RNBoolean already_sorted, RNBoolean already_unique)
  : n(0),
    c(_c),
    v(NULL),
    e(NULL)
{
  // Check number of variables
  if (_n > 0) {
    // Copy term info
    n = 0;
    v = new int [ _n ];
    e = new RNScalar [ _n ];
    for (int i = 0; i < _n; i++) {
      if (_e[i] != 0.0) {
        assert(n < _n);
        v[n] = _v[i];
        e[n] = _e[i];
        n++;
      }
    }

    // Sort variables
    if (!already_sorted) {
      for (int i = 1; i < n; i++) {
        for (int j = i; j > 0; j--) {
          if (v[j] < v[j-1]) {
            assert(j > 0);
            assert(j < _n);
            int swap = v[j-1];
            v[j-1] = v[j];
            v[j] = swap;
          }
        }
      }
    }

    // Check for duplicate variables
    if (!already_unique) {
      for (int i = 0; i < n-1; i++) {
        assert(i >= 0);
        assert(n >= 0);
        assert(i < n-1);
        if (v[i] == v[i+1]) {
          assert(i < n-1);
          e[i] += e[i+1];
          for (int j = i+1; j < n-1; j++) {
            assert(j >= 0);
            assert(j < n-1);
            v[j] = v[j+1];
            e[j] = e[j+1];
          }
          n -= 1;
          i--;
        }
      }
    }
  }
}



RNPolynomialTerm::
RNPolynomialTerm(const RNPolynomialTerm& term)
  : n(term.n),
    c(term.c),
    v(NULL),
    e(NULL)
{
  // Copy stuff from term
  if (n > 0) {
    v = new int [ n ];
    e = new RNScalar [ n ];
    for (int i = 0; i < n; i++) {
      v[i] = term.v[i];
      e[i] = term.e[i];
    }
  }
}



RNPolynomialTerm::
~RNPolynomialTerm(void)
{
  // Delete stuff
  if (v) delete [] v;
  if (e) delete [] e;
}




RNScalar RNPolynomialTerm::
Degree(void) const
{
  // Return degree
  RNScalar degree = 0;
  for (int i = 0; i < NVariables(); i++) degree += e[i];
  return degree;
}



void RNPolynomialTerm::
Empty(void)
{
  // Delete stuff
  if (v) delete [] v;
  if (e) delete [] e;
  c = 0;
  n = 0;
}



RNScalar RNPolynomialTerm::
PartialDerivative(const RNScalar *x, int variable) const
{
  // Evaluate the derivative with respect to variable
  RNScalar result = c;
  for (int i = 0; i < n; i++) {
    if (variable == v[i]) {
      if (e[i] == 1.0) result *= 1.0;
      else if (e[i] == 2.0) result *= 2.0 * x[v[i]];
      else if ((e[i] < 1.0) && RNIsZero(x[v[i]])) result *= RN_INFINITY;
      else result *= e[i] * pow(x[v[i]], e[i] - 1.0);
    }
    else {
      if (e[i] == 1.0) result *= x[v[i]];
      else if (e[i] == 2.0) result *= x[v[i]] * x[v[i]];
      else if ((e[i] < 0) && RNIsZero(x[v[i]])) result *= RN_INFINITY;
      else result *= pow(x[v[i]], e[i]);
    }
  }
  return result;
}



RNBoolean RNPolynomialTerm::
IsConstant(void) const
{
  // Return whether term is a constant
  if (n == 0) return TRUE;
  return FALSE;
}



RNBoolean RNPolynomialTerm::
IsLinear(void) const
{
  // Return whether term is linear
  if ((n == 1) && (e[0] == 1.0)) return TRUE;
  return FALSE;
}



RNBoolean RNPolynomialTerm::
IsQuadratic(void) const
{
  // Return whether term is quadratic
  if ((n == 1) && (e[0] == 2.0)) return TRUE;
  if ((n == 2) && (e[0] == 1.0) && (e[1] == 1.0)) return TRUE;
  return FALSE;
}



RNScalar RNPolynomialTerm::
Evaluate(const RNScalar *x) const
{
  // Evaluate the term
  RNScalar result = c;
  for (int i = 0; i < n; i++) {
    if (e[i] == 1.0) result *= x[v[i]];
    else if (e[i] == 2.0) result *= x[v[i]] *x[v[i]];
    else if ((e[i] < 0) && RNIsZero(x[v[i]])) result *= RN_INFINITY;
    else result *= pow(x[v[i]], e[i]);
  }
  return result;
}



RNBoolean RNPolynomialTerm::
HasVariables(int query_n, const int *query_v, const RNScalar *query_e) const
{
  // Check if query has same number of variables
  if (n != query_n) return FALSE;

  // Compare variables and exponents 
  for (int i = 0; i < n; i++) {
    // This assumes variables are sorted
    if (v[i] != query_v[i]) return FALSE;
    if (e[i] != query_e[i]) return FALSE;
  }

  // Passed all tests
  return TRUE;
}



void RNPolynomialTerm::
Print(FILE *fp) const
{
  // Print term
  fprintf(fp, "%g ", Coefficient());
  for (int j = 0; j < NVariables(); j++) {
    fprintf(fp, "(%d %g) ", Variable(j), Exponent(j));
  }
  fprintf(fp, "\n");
}



