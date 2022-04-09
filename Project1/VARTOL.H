#pragma once

#ifdef BLADEMATH_EXPORTS
#define DLLEXPORT __declspec(dllexport)
#else // !defined (BLADEMATH_EXPORTS)
#define DLLEXPORT __declspec(dllimport)
#endif

#include "nurb.h"

class CSection;  // not sure why I need to do this

class DLLEXPORT CVariableTol
{
public:
  CVariableTol();
  ~CVariableTol();

  // for the following, index 0 is lowertol, index 1 is uppertol

  double *m_x[2];       // x coordinates of points in polygon
  double *m_y[2];       // y coordinates of points in polygon
  double m_letol[2];    // tolerances for leading edge
  double m_stacktol[2]; // tolerances at stack point
  double m_tetol[2];    // tolerances at trailing edge
  double m_chordtol[2]; // tolerances to allow for chord length differences
  double m_minouttol;   // largest out of tolerance condition on minus side
  double m_maxouttol;   // largest out of tolerance condition on plus side
  CNurbCurve *m_tol[2]; // tolerance curves

  int m_numPoints;      // number of points in polygons
  double m_tnose;       // spline parameter for nominal nose point
  double m_ttail;       // spline parameter for nominal tail point
  double m_leoff;       // le transition point distance
  double m_teoff;       // le transition point distance

  CCurve *m_nom;        // nominal must be a single closed curve
  CCurve *m_mcc;        // nominal mean camber curve

  bool m_teSquare;      // true if trailing edge is nominally square
  double m_teWidth;     // width of te when it is square

  void SetOffsets(double leo, double teo, double let, double tet, bool teSquare, double teWidth)
  {
    m_leoff = leo;
    m_teoff = teo;
    m_tnose = let;
    m_ttail = tet;
    m_teSquare = teSquare;
    m_teWidth = teWidth;
  }

  void SetTols(double lemtol, double leptol, double stmtol, double stptol,
               double temtol, double teptol, double cdmtol, double cdptol)
  {
    m_letol[0]    = lemtol;
    m_letol[1]    = leptol;
    m_stacktol[0] = stmtol;
    m_stacktol[1] = stptol;
    m_tetol[0]    = temtol;
    m_tetol[1]    = teptol;
    m_chordtol[0] = cdmtol;
    m_chordtol[1] = cdptol;
  }

  void SetNom(CCurve *c, CCurve *m)
  {
    m_nom = c;
    m_mcc = m;
  }

  int CalcTol(int nle, int nte, int nside);
 // int CalcTol(CSection *sect);
  int CalcTol(int nlo, double *xlo, double *ylo, int nhi, double *xhi, double *yhi);

  int Ramp(int &n, int num, double t0, double *tol0, double t1, double *tol1, int bias);
  int Stretch(int &n, int num, double t0, double t1, double *tol1,
              double *ctr, double *end, double *tolend);

  int OutTol(double *act, double *tol, double *out, int n);
};
