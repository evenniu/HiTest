#pragma once

#ifdef BLADEMATH_EXPORTS
#define DLLEXPORT __declspec(dllexport)
#else // !defined (BLADEMATH_EXPORTS)
#define DLLEXPORT __declspec(dllimport)
#endif

class DLLEXPORT CVector
{
public:

  CVector();
  CVector(int nvals);
  CVector(const CVector&);       // copy constructor
  CVector& operator=(CVector&);  // vector assignment
  ~CVector();
  int N() const
  {
    return n;
  }

  double *v;
protected:
  void Kill();
  void Create();

  int n;
};

class DLLEXPORT CMatrix
{
public:
  CMatrix();
  ~CMatrix();
  CMatrix(int rws, int cls);
  CMatrix(const CMatrix& m);       // copy constructor

  CMatrix &operator=(const CMatrix& mat);
  CMatrix &operator+=(const CMatrix& mat);
  CMatrix &operator-=(const CMatrix& mat);
  CMatrix &operator*=(const CMatrix& mat);

  friend CMatrix operator+(const CMatrix &m1, const CMatrix &m2);
  friend CMatrix operator-(const CMatrix &m1, const CMatrix &m2);
  friend CMatrix operator*(const CMatrix &m1, const CMatrix &m2);

  double **m;
  int R() const
  {
    return r;
  }
  int C() const
  {
    return c;
  }
  int Valid()
  {
    return good;
  }
  void Zero();
private:
  int r,c,good;

  void Allocate();
  void Destroy();
};

