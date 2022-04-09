#pragma once

#include <stdio.h>
#include "Alignment.h"
#include "AbstractCurve.h"

using namespace std;
struct DumpRecord
{
  double v[2];
  DumpRecord(double v0, double v1)
  {
    v[0] = v0;
    v[1] = v1;
  }
};
#include <vector>
typedef std::vector < DumpRecord * > CDumpRecordArray;

#ifdef BLADEMATH_EXPORTS
#define DLLEXPORT __declspec(dllexport)
#else // !defined (BLADEMATH_EXPORTS)
#define DLLEXPORT __declspec(dllimport)
#endif

#define NURB_TYPE        1
#define CIRCLE_TYPE      2
#define POINT_TYPE       3
#define SUBCURVE_TYPE    4
#define LINE_TYPE        5
#define CONCATCURVE_TYPE 6
#define HERMITE_TYPE     7

class DLLEXPORT CCurve : public Hexagon::Blade::Curve<2>
{

protected:

  bool m_isEnglish;       // english or metric?
  double m_period;        // only for periodic curves and subcurves
  double m_extreme;       // t value for extreme points, used for LE and TE
   CAlignment *m_palign;   // alignment to run calcs through (or null);
public:


  // required functions for derived classes

  CCurve();
  virtual ~CCurve(); // destructor
  virtual int Type() const = 0;
  virtual int CalcPoint(double* xyz, double t, double* tan = 0, double* curv = 0) const = 0;
  virtual int Valid() = 0;
  virtual long StoreSize() = 0;
  virtual void Write(FILE* fp) = 0;
  virtual void Extent(double* min, double* max) = 0;
  virtual double T0() const = 0;
  virtual double T1() const = 0;
  virtual void SetT0(double t) = 0;
  virtual void SetT1(double t) = 0;
  virtual int CalcPoints(double* outPoints, const double* inT, ptrdiff_t numberOfT, double* outTangents = nullptr,
                         double* outCurves = nullptr) const;
  virtual void FindClosestTValues(double* outT, const double* inPoints, ptrdiff_t numberOfPoints);
  virtual double compute_length(double t0, double t1);
  virtual int Extreme(double* ij, double* t, double* p = 0, double t0 = 0.0, double t1 = 0.0, int save = 1,
                      bool endOK = false);
  virtual void Extreme(double extr);
  virtual int MatchTan(double* tv, double* t, double* p, double t0 = 0.0, double t1 = 0.0, bool debug = false);
  virtual int HighPoint(double* ij, double* t, double* p = 0, double t0 = 0.0, double t1 = 0.0, bool endOK = false);
  virtual double TotalLength(int numPoints, double t0, double t1);
  virtual void Walk(double& t, double d, double* pt, int numPoints);
  virtual double arc_inverse(double length, double t0, int& index);
  virtual void CenterOfMass(double* CM, bool isGEI);
  virtual double ComputeArea(double t0 = 0.0, double t1 = 0.0);

  // generic functions

  int EndCamber(int n, double* x, double* y, int useextreme = 0);
  double ClosestPoint(double* tgt, double* bestxyz, double* bestt = 0, double* tangent = 0, double t0 = 0.0,
                      double t1 = 0.0, int nseed = 10);
  double NewClosestPoint(const double* tgt, double* bestxyz, double* bestt = 0, double* tangent = 0, double t0 = 0.0,
                         double t1 = 0.0, int nseed = 10);
  double ClosestNominal(double* tgt, double* norm, double* bestxyz, double* bestt = 0, double* tangent = 0,
                        double t0 = 0.0, double t1 = 0.0, int nseed = 10, int match = 0);

  int LineIntersect(double* xy, double* ij, double* sol, double t0 = 0.0, double t1 = 0.0, double* tsol = 0,
                    double* dist = 0, bool isRay = false, bool useSol = false);
  virtual int NewLineIntersect(double* xy, double* ij, double* sol, double t0 = 0.0, double t1 = 0.0, double* tsol = 0,
                               double* dist = 0, bool isRay = false, bool useSol = false,
                               const ptrdiff_t numberOfSeeds = 500) const final;

  virtual int CircIntersect(const double* xy, double r, double* sol, double t0 = 0.0, double t1 = 0.0, double* tsol = 0,
                            double* tanv = 0) const final;

  bool AdvanceByLength(double delta, double t0, double t1, double* pt, double* tsol);

  void Align(CAlignment* al);
  void IsEnglish(bool eng);
  bool IsEnglish();
  double Extreme() const;
  double Period();
  void BugOut(int numPoints = 100, double t0 = 0.0, double t1 = 0.0);

  // during reporting, it can be useful to retrieve information about the curve;
  // I call these annotations
  virtual void addDoubleAnnotation(const char* key, double value) = 0;
  virtual double getDoubleAnnotation(const char* key) const = 0;
  virtual void addBoolAnnotation(const char* key, bool value) = 0;
  virtual bool getBoolAnnotation(const char* key) const = 0;

  // Inherited via Curve
  virtual bool isPeriodic() const override;
  virtual void parametricBounds(double* const outBounds) const override;
  virtual void evaluate(const double* const t, const ptrdiff_t n, double* const outPoints,
                        double* const outDerivative, double* const outSecondDerivative) const override;
  virtual double t0() const override;
  virtual double t1() const override;
  virtual double period() const override;
};

extern int ACamberPoint(CCurve *a, CCurve *b, double *g, double *p1, double *p2,
                        double *r1, double ta, double tb0, double tb1,
                        double *tb = 0);

extern int CamberPoints(int n, CCurve *a, CCurve *b, double *x, double *y,
                        double ta0=0.0, double ta1=0.0,
                        double tb0=0.0, double tb1=0.0, int useextreme=0);
