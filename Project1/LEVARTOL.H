#pragma once

#include "nurb.h"

#ifdef BLADEMATH_EXPORTS
#define DLLEXPORT __declspec(dllexport)
#else // !defined (BLADEMATH_EXPORTS)
#define DLLEXPORT __declspec(dllimport)
#endif

class DLLEXPORT CLEVariableTol
{
public:
  CLEVariableTol();
  ~CLEVariableTol();

  // for the following, index 0 is lowertol, index 1 is uppertol

  double *m_x[2];       // x coordinates of points in polygon
  double *m_y[2];       // y coordinates of points in polygon
  double m_letol[2];    // tolerances for leading edge
  double m_stacktol[2]; // tolerances at stack point
  double m_minouttol;   // largest out of tolerance condition on minus side
  double m_maxouttol;   // largest out of tolerance condition on plus side
  CNurbCurve *m_tol[2]; // tolerance curves

  int m_numPoints;      // number of points in polygons
  double m_tnose;       // spline parameter for nominal nose point
  double m_t0;          // spline parameter for beginning of LE
  double m_t1;          // spline parameter for end of LE

  CCurve *m_nom;        // nominal must be a single closed curve

  void SetOffsets(double lt0, double lt1, double let)
  {
    m_t0 = lt0;
    m_t1 = lt1;
    m_tnose = let;
  }

  void SetTols(double lemtol, double leptol, double stmtol, double stptol)
  {
    m_letol[0]    = lemtol;
    m_letol[1]    = leptol;
    m_stacktol[0] = stmtol;
    m_stacktol[1] = stptol;
  }

  void SetNom(CCurve *c)
  {
    m_nom = c;
  }

  int CalcTol(int nle);

  int Ramp(int &n, int num, double t0, double *tol0, double t1, double *tol1);

  int OutTol(double *act, double *tol, double *out, int n);
};

