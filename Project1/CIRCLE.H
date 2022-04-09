#pragma once

#ifdef BLADEMATH_EXPORTS
#define DLLEXPORT __declspec(dllexport)
#else // !defined (BLADEMATH_EXPORTS)
#define DLLEXPORT __declspec(dllimport)
#endif

// WWW  warning C4266: 'void CCurve::Extreme(double)' : no override available for virtual member function from base 'CCurve'; function is hidden

class DLLEXPORT CCircle : public CCurve
{
public:
  CCircle(double *c, double *s, double *e, double r);
  CCircle(int n, double *x, double *y);
  ~CCircle();

  // during reporting, it can be useful to retrieve information about the curve;
  // I call these annotations
  virtual void addDoubleAnnotation(const char* key, double value) override final;
  virtual double getDoubleAnnotation(const char* key) const override final;
  virtual void addBoolAnnotation(const char* key, bool value) override final;
  virtual bool getBoolAnnotation(const char* key) const override final;

  int Type() const override
  {
    return rad <= 0.0 ? POINT_TYPE : CIRCLE_TYPE;
  }

  void Seqno(long s)
  {
    seqno = s;
  }
  long StoreSize()
  {
    return sizeof(short) + sizeof(double) + storesize;
  }
  int CalcPoint(double *xyz, double a, double *tan = 0, double *curv = 0) const override;
  void Write(FILE *);
  double ClosestPoint(double *tgt, double *bestxyz, double *bestt = 0, double *tangent = 0, double t0=0.0, double t1=0.0,
    int nseed = 100);
  double ClosestNominal(double * /*tgt*/, double * /*norm*/, double * /*bestxyz*/, double * /*bestt = 0*/, double * /*tangent = 0*/, double/* t0=0.0*/,
    double /*t1=0.0*/, int /*nseed = 10*/, int /*match = 0*/)
  {
    return 0.0;  // don't think we need this
  }
  int LineIntersect(double *xy, double *ij, double *sol, double t0=0.0,
    double t1=0.0, double *tsol=0, double *dist=0, bool isRay = false, bool useSol = false);
  void Extent(double *min, double *max);
  void Extent();
  int Extreme(double *ij, double *t, double *p = 0, double t0 = 0.0, double t1 = 0.0, int save = 1, bool endOK = false);
  double Rad()
  {
    return rad;
  }
  double X()
  {
    return center[0];
  }
  double Y()
  {
    return center[1];
  }
  double Z()
  {
    return center[2];
  }
  double *Center()
  {
    return center;
  }
  int Valid()
  {
    return good;
  }
  double T0() const override
  {
    return a0;
  }
  double T1() const override
  {
    return a1;
  }
  void SetT0(double t)
  {
    a0 = t;
  }
  void SetT1(double t)
  {
    a1 = t;
  }

  double ls_err;

protected:
  int FitCircle(int npts, double *xp, double *yp);
  int CircInit(int n, double *x, double *y);
  double norm(int n, double *dx);
  void multATx(int rows, int cols, double **A, double *x, double *y);
  int LSnormal(int m, int n, double **A, double *b, double *x);
  int LDLT(int n, double **A);
  int LDLTsolve(int n, double **A, double *b);
  void lower_backsub(int n, double **L, double *b);
  void upper_backsub(int n, double **L, double *y);

  int good;
  double center[3];
  double a0;
  double a1;
  double Min[2];
  double Max[2];
  double rad;
  long seqno;
  long storesize;  // bytes needed to write to binary file
};
