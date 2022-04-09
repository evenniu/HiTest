#pragma once

#ifdef BLADEMATH_EXPORTS
#define DLLEXPORT __declspec(dllexport)
#else // !defined (BLADEMATH_EXPORTS)
#define DLLEXPORT __declspec(dllimport)
#endif

class DLLEXPORT ODRSpline
{
public:

  ODRSpline(int np, double *px, double *py, double *pz, bool isEnglish, int cl, int tn=1, double th=0.0001);

  ~ODRSpline();
  int isGood()
  {
    return m_good;
  }
  int GetNctl()
  {
    return m_nctl;
  }
  double *GetKnots()
  {
    return m_pknots;
  }
  double *GetCx()
  {
    return m_pcx;
  }
  double *GetCy()
  {
    return m_pcy;
  }
  double *GetCk()
  {
    return m_pck;
  }

private:

  int transform_points(int n, double **x, double *X0);
  int thin(int n, double *x, double *y, double *k, int *loc_parms);
  int closed_curve(int m, double *x, double *y, double *k, double *t,
    double *cx, double *cy, double *ck, int kts);
  int init_spline(int m, double *x, double *y, double *k, double *t,
    double *cx, double *cy, double *ck);
  int closed_knots(int m, double *x, double *y, double *t);
  int setcmatrix(int m, double *a, double *d, double *b, double *t);
  int solve_ctrid(int m, double *a, double *d, double *b, double *r,
    double *x);
  int chord(int m, double *x, double *y, double *t);
  int setmatrix(int L, double *a, double *d, double *b, double *t);
  int setrhs(int m, double *t, double *P);

  bool m_english;
  int m_good;
  int m_npts;
  int m_ndim;
  int m_nctl;
  int m_closed;
  int m_tension;

  double m_scale_fac;
  double m_thin_tol;
  double m_X0[3];

  double *m_pknots;
  double *m_pcx;
  double *m_pcy;
  double *m_pck;
  double *m_pa;
  double *m_pb;
  double *m_pd;
  double *m_puu;
  double *m_pzz;
  double *m_pcoeffs;
  double *m_pbeta;

  CMatrix *m_pxm;
};
