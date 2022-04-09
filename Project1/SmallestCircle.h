#pragma once

#ifdef BLADEMATH_EXPORTS
#define DLLEXPORT __declspec(dllexport)
#else // !defined (BLADEMATH_EXPORTS)
#define DLLEXPORT __declspec(dllimport)
#endif

class DLLEXPORT CSmallestCircle
{
public:

  int m_numPoints;
  double *m_X;
  double *m_Y;

  double m_center[2];
  double m_radius;
  double m_leTarget[2];
  double m_teTarget[2];
  double m_lePoint[2];
  double m_tePoint[2];
  double m_circumference;

  CSmallestCircle(int n, double *x, double *y, double *leTarget, double *teTarget);
  CSmallestCircle(int n, double *x, double *y);
  ~CSmallestCircle(void);

  void FindCircle();  
};
