#include "stdafx.h"

#ifdef _DEBUG
   #ifndef DBG_NEW
      #define DBG_NEW new ( _NORMAL_BLOCK , __FILE__ , __LINE__ )
      #define new DBG_NEW
   #endif
#endif  // _DEBUG

void CVector::Kill()
{
  if (v && n > 0)
    delete [] v;

  v = 0;
  n = 0;
}

void CVector::Create()
{
  if (n < 1)
  {
    n = 0;
    return;
  }

  v = new double [n];
  if (!v)
  {
    n = 0 ;
    return;
  }

  int i;
  for (i=0; i<n; i++)
    v[i] = 0.0;
}

CVector::CVector() : v(0), n(0)
{}

CVector::CVector(int nvals) : v(0), n(nvals)
{
  Create();
}

CVector::~CVector()
{
  Kill();
}

CVector::CVector(const CVector& x)
{
  n = x.n;
  v = new double [x.n];

  int i;
  for (i=0;i<x.n;i++)
    v[i] = x.v[i];
}

CVector& CVector::operator=(CVector& x)
{
  if (n != x.n)
  {
    Kill();
    n = x.n;
    Create();
  }

  int i;
  for (i=0;i<n;i++)
    v[i] = x.v[i];

  return *this;
}

CMatrix::CMatrix() : r(0), c(0), m(0), good(0)
{}

CMatrix::CMatrix(int rws, int cls) : r(rws), c(cls), m(0)
{
  Allocate();
}

CMatrix::CMatrix(const CMatrix &mat)
{
  r = mat.r;
  c = mat.c;

  Allocate();

  int i,j;
  for (i=0; i<r; i++)
    for (j=0; j<c; j++)
      m[i][j] = mat.m[i][j];
}

void CMatrix::Allocate()
{
  good = 0;
  if (r < 1 || c < 1)    // no good
  {
    r = c = 0;
    return;
  }

  m = new double *[r];
  if (!m)
    return;

  m[0] = new double[r*c];

  if (!m[0])
  {
    delete m;
    return;
  }

  int i;
  for (i=1; i<r; i++)
    m[i] = m[i-1] + c;

  for (i=0; i<r*c; i++)
    m[0][i] = 0.0;

  good = 1;
}

void CMatrix::Destroy()
{
  if (m)
  {
    delete [] m[0];
    delete [] m;
  }
  r = c = 0;
  m = 0;
  good = 0;
}

CMatrix::~CMatrix()
{
  Destroy();
}

void CMatrix::Zero()
{
  if (!good)
    return;

  int i,j;
  for (i=0; i<r; i++)
    for (j=0; j<c; j++)
      m[i][j] = 0.0;
}

CMatrix &CMatrix::operator=(const CMatrix &mat)
{
  if (r != mat.r && c != mat.c)
  {
    Destroy();

    r = mat.r;
    c = mat.c;

    Allocate();
  }

  if (good)
  {
    int i,j;
    for (i=0; i<r; i++)
      for (j=0; j<c; j++)
        m[i][j] = mat.m[i][j];
  }

  return *this;
}

CMatrix &CMatrix::operator+=(const CMatrix & mat)
{
  if (r == mat.R() && c == mat.C())
  {
    int i,j;
    for (i=0; i<r; i++)
      for (j=0; j<c; j++)
        m[i][j] += mat.m[i][j];
  }

  return *this;
}

CMatrix &CMatrix::operator-=(const CMatrix & mat)
{
  if (r == mat.R() && c == mat.C())
  {
    int i,j;
    for (i=0; i<r; i++)
      for (j=0; j<c; j++)
        m[i][j] -= mat.m[i][j];
  }

  return *this;
}

CMatrix &CMatrix::operator*=(const CMatrix & mat)
{
  if (c == mat.R())
  {
    int i,j;
    CMatrix res(r, mat.C());
    for (i=0; i<r; i++)
    {
      for (j=0; j<mat.C(); j++)
      {
        res.m[i][j] = 0.0;
        int k;
        for (k=0; k<c; k++)
          res.m[i][j] = m[i][j] * mat.m[k][j];
      }
    }

    Destroy();
    r = res.r;
    c = res.c;
    Allocate();

    for (i=0; i<r; i++)
      for (j=0; j<c; j++)
        m[i][j] = res.m[i][j];
  }

  return *this;
}

CMatrix operator+(const CMatrix &m1, const CMatrix &m2)
{
  return CMatrix(m1) += m2;
}

CMatrix operator-(const CMatrix &m1, const CMatrix &m2)
{
  return CMatrix(m1) -= m2;
}

CMatrix operator*(const CMatrix &m1, const CMatrix &m2)
{
  return CMatrix(m1) *= m2;
}
