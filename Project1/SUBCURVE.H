#pragma once

#include <stdio.h>

#ifdef BLADEMATH_EXPORTS
#define DLLEXPORT __declspec(dllexport)
#else // !defined (BLADEMATH_EXPORTS)
#define DLLEXPORT __declspec(dllimport)
#endif

class DLLEXPORT CSubCurve : public CCurve
{
public:
  CSubCurve(FILE *fp, int n=0);
  CSubCurve(CCurve *c, double s, double e, double p);   // construct subcurve
  ~CSubCurve()
  {}

  int Type() const override
  {
    return SUBCURVE_TYPE;
  }
  bool isPeriodic() const override
  {
    return false;
  }

  int Valid()
  {
    return m_psb ? 1 : 0;
  }
  double T0() const override
  {
    return m_st0;
  }
  double T1() const override
  {
    return m_st1;
  }
  void SetT0(double t)
  {
    m_st0 = t;
  }
  void SetT1(double t)
  {
    m_st1 = t;
  }
  long StoreSize()
  {
    return 0;
  }  // this shouldn't get called
  void Write(FILE *fp);
  void Extent();
  void Extent(double *min, double *max);

  int CalcPoint(double *xyz, double a, double *tan = 0, double *curv = 0) const override;

  int Extreme(double *ij, double *t, double *p = 0,
    double t0 = 0.0, double t1 = 0.0, int save = 1, bool endOK = false);
  virtual void Extreme(double extr)
  {
    m_extreme = extr;
  }

  double ClosestPoint(double *tgt, double *bestxyz, double *bestt = 0, double *tangent = 0, double t0=0.0, double t1=0.0,
    int nseed = 100)
  {
    if (m_psb)
      return m_psb->ClosestPoint(tgt, bestxyz, bestt, tangent, t0, t1, nseed);
    return 0.0;
  }

  double ClosestNominal(double *tgt, double *norm, double *bestxyz,
    double *bestt = 0, double *tangent = 0, double t0=0.0,
    double t1=0.0, int nseed = 10, int match = 0)
  {
    if (m_psb)
      return m_psb->ClosestNominal(tgt, norm, bestxyz, bestt, tangent, t0, t1,
      nseed, match);
    return 0.0;
  }
  int LineIntersect(double *xy, double *ij, double *sol, double t0=0.0,
    double t1=0.0, double *tsol=0, double *dist=0, bool isRay = false, bool useSol = false)
  {
    if (m_psb)
      return m_psb->LineIntersect(xy, ij, sol, t0, t1, tsol, dist, isRay, useSol);
    return 0;
  }

  double TotalLength(int numPoints, double t0, double t1)
  {
    if (m_psb)
      return m_psb->TotalLength(numPoints, t0, t1);
    else
      return 0.0;
  }
  // protected:

  double m_min[2], m_max[2];
  double m_st0, m_st1;
  CCurve *m_psb;

  // during reporting, it can be useful to retrieve information about the curve;
  // I call these annotations
  virtual void addDoubleAnnotation(const char* key, double value) override final;
  virtual double getDoubleAnnotation(const char* key) const override final;
  virtual void addBoolAnnotation(const char* key, bool value) override final;
  virtual bool getBoolAnnotation(const char* key) const override final;
};
