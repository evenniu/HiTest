#pragma once

#include "Alignment.h"
#include "VarTol.h"
#include "LEVarTol.h"
#include "CONSTANTS.H"
#include "ToleranceSection.h"

#ifdef BLADEMATH_EXPORTS
#define DLLEXPORT __declspec(dllexport)
#else // !defined (BLADEMATH_EXPORTS)
#define DLLEXPORT __declspec(dllimport)
#endif

class DLLEXPORT CFitParams
{
public:
  int fitcurve[4];
  int weightcurve[4];
  BestFitAlgorithm algorithm;
  int rotfit;
  int tranfit;
  int fitToMiddleOfZone;
  int pivot;
  int usenominals;
  double lepercent;
  double tepercent;
  // these leading and trailing edge offsets are used for Leading Edge Arc Least Squares and Trailing Edge Arc Least
  // Squares
  double leoff1, leoff2;
  double teoff1, teoff2;
  double rotMTol;
  double rotPTol;
  double tranPTol[2];
  double tranMTol[2];
  double profileMTol;
  double profilePTol;
  double forgeAngle;

  // these zone region distances and tolerances are used for Complex Edge best-fits and Chord best-fits
  const ComplexZone* complexEdgeZone;
  const ComplexZone* chordZone;

  CFitParams();
  ~CFitParams()
  {
  }
};

class DLLEXPORT CBestFit
{
public:
  CBestFit(int np);
  ~CBestFit();

  int m_totalBad;
  int m_totalChecked;

  double m_lserr;
  double m_mindev[4]; // furthest distance on minus side
  double m_maxdev[4]; // furthest distance on plus side
  double m_meandev[4];
  double m_stddev[4];
  // following 3 arrays are for variable tolerancing.

  double* m_btxpt;   // closest x values to tolerance curves.
  double* m_btypt;   // closest y values to tolerance curves.
  double* m_outtol;  // amount out of tolerance can be + or -, 0=intol.
  int* m_bestPartOf; // after best fit, point originally part of

  bool Omega();                // fill omega array
  bool Omega(int i, double v); // fill omega array
  double Omega(int i)
  {
    if(m_omega && i >= 0 && i < m_numpts)
      return m_omega->v[i];
    return 0.0;
  }

  bool LeastSquaresFit(CFitParams& fp, double* nompt = 0, double* meapt = 0);
  // bool TwoPointFit(double *n1, double *n2, double *m1, double *m2);
  // bool SetMinMaxY(double twist);
  // bool GuillotineFit(CCurve *c, CFitParams &fp);
  // bool VectorFit(CFitParams &fp, bool init);

  void XYDelta(double* xy);
  void RotDelta(double angle);

  int NumPoints()
  {
    return m_numpts;
  }

  void PutVal(int i, double v1, double v2); // copy vector to vals array
  void PutNom(int i, double v1, double v2); // copy vector to noms array
  void PutT(int i, double t);
  void PutVec(int i, double v1, double v2); // copy vector to ijks array
  bool GetInf(int i, double* v);            // copy vector from infs array
  void PutInf(int i, double v1, double v2); // copy vector to vals array
  bool GetNom(int i, double* v);            // copy vector from noms array
  bool GetVec(int i, double* v);            // copy vector from ijks array
  bool GetVal(int i, double* v);            // copy vector from vals array
  bool GetT(int i, double* t);

  bool Align(CAlignment* v);       // copy alignment out
  bool ImportAlign(CAlignment* v); // copy alignment in
  const CAlignment* GetAlign() const
  {
    return &m_align;
  }
  CAlignment* GetAlign()
  {
    return &m_align;
  }
  void ClearAlign()
  {
    m_align.Clear();
  }

  CVariableTol* VarTol()
  {
    return m_vartol;
  }
  void VarTol(CVariableTol* vt)
  {
    m_vartol = vt;
  }
  int ApplyVarTol(double* minOut, double* maxOut);
  int ReturnVarTol(double* minOut, double* maxOut);

  CLEVariableTol* LEVarTol()
  {
    return m_LEvartol;
  }
  void LEVarTol(CLEVariableTol* vt)
  {
    m_LEvartol = vt;
  }
  int ApplyLEVarTol(double* minOut, double* maxOut);
  int ReturnLEVarTol(double* minOut, double* maxOut);

  void ReportFit(int rf)
  {
    m_reportFit = rf;
  }
  int ReportFit()
  {
    return m_reportFit;
  }
  bool ReturnFit(double* x, double* y, double* theta, int method = -1);

  // move the following back to protected after making access methods

  int m_numpts;
  int m_xlate;
  int m_reportFit;

  double m_MinMaxY, m_PTol, m_MTol;

  CMatrix* m_vals;
  CMatrix* m_infs;
  CMatrix* m_noms;
  CMatrix* m_ijks;
  CVector* m_t;

  CVector* m_omega;

  CAlignment m_align;

  CVariableTol* m_vartol;     // Variable Tolerance
  CLEVariableTol* m_LEvartol; // LE Variable Tolerance

  CFitParams m_fitParams;
  bool* m_valWasUsedInFit;

protected:
  //// vector fit members.
  // CMatrix *m_Q;
  // CMatrix *m_Q1;
  // CVector *m_b;
  // CVector *m_xbar;

  // CCurve *m_gcurve;
  // double (CBestFit::* zfunc)(double);

  // double GuillotineMetric(double theta);
  // double MinMaxMetric(double theta);
  // double zeroin(double a, double b, double tol, double *fval);

  //// vector internal fit functions
  // int Qfunc(double beta, double **Q);
  // int new_grads(int n, double *beta, double *resvec, double **gradbeta);
  // int qd_qrls(int first, int m, int n, double **A, double *b, double eps, double *sqe);
  // int resid(int n, double beta, double *resvec, double *sqe);
  // void compscale(int n, double **A, double *scale);
  // void doscale(int n, double **A, double scale);
  // int control(int n, double *beta, double *dist2, double *lnorm, double *lambda, double *bnorm);
  // double normgrad(int n, double **gradbeta, double *dist);
  // int levmq(int count, int n, double **A, double *b, double *dbeta, double delta, double *bnorm, double *lambda,
  //  double *linres);
  // int Householder(int m, int n, double **A, double *beta, double *r, int *piv, double eps);
  // void House_mult(int rank, int m, double **U, double *beta, double *y);
  // void House_backsub(double **R, double *r, int rank, double *x, double *b);
  // int update_XZ(double **X, double **Z, double **T, double *b);
  // int update_Y(double **X, double **Z, double **T, double **Q, double **Y, double *b);
  // void matset(int n, int np, double lambda, double **A, double *b, double **AL, double *bl, int first);
  // int QRls(int m, int n, double **A, double *b, double eps, double *sqe);
};