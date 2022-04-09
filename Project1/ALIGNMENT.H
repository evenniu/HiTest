#pragma once

#ifdef BLADEMATH_EXPORTS
#define DLLEXPORT __declspec(dllexport)
#else // !defined (BLADEMATH_EXPORTS)
#define DLLEXPORT __declspec(dllimport)
#endif

class DLLEXPORT CAlignment
{
public:

  CAlignment();
  CAlignment(int d);
  ~CAlignment();
  CAlignment(CAlignment &obj);                      // copy constructor
  CAlignment& operator= (const CAlignment & obj);  // assignment operator

  int n;  // dimension

  double m_mat[3][3];
  double m_morig[3];
  double m_borig[3];

  void Clear();
  bool Rotate(int z, double rot);
  bool MatrixFromAxis(double *zaxis);
  void MeasToBest(const double *meas, int point, double *best) const;
  void MeasToBestVector2Xd(double* meas, int point, double* best, ptrdiff_t numberOfElements) const;
  void BestToMeas(double* best, int point, double* meas) const;
  //void BugOut(CString &s);
};
