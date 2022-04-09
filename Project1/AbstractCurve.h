#pragma once
#include <stdint.h>

#ifdef BLADEMATH_EXPORTS
#define DLLEXPORT __declspec(dllexport)
#else // !defined (BLADEMATH_EXPORTS)
#define DLLEXPORT __declspec(dllimport)
#endif

namespace Hexagon
{
namespace Blade
{

// curves are parametrized with a single parameter t and may be periodic
template <ptrdiff_t dimension>
struct Curve
{
  virtual ~Curve()
  {
  }

  virtual bool isPeriodic() const = 0;
  virtual double t0() const = 0;
  virtual double t1() const = 0;
  virtual double period() const = 0; // is allowed to throw an error if non-periodic
  virtual void parametricBounds(double* const outBounds) const = 0;
  virtual void evaluate(const double* const t, const ptrdiff_t n, double* const outPoints, double* const outDerivative,
                        double* const outSecondDerivative) const = 0;
};

// curves are parametrized with a single parameter t and may be periodic
template <ptrdiff_t dimension>
struct CurveWithPreferredDomain : public Curve<dimension>
{
  virtual ~CurveWithPreferredDomain()
  {
  }

  virtual bool isInPreferredDomain(double t) const = 0;
};

DLLEXPORT int extreme(const Hexagon::Blade::Curve<2>& curve, const double* ij, double* t, double* p = 0,
                      double t0 = 0.0, double t1 = 0.0, bool endOK = false);

DLLEXPORT double closestPoint(const Hexagon::Blade::Curve<2>& curve, const double* tgt, double* bestxyz,
                              double* bestt = 0, double* tangent = 0, double t0 = 0.0, double t1 = 0.0, int nseed = 10);

DLLEXPORT int circleIntersection(const Hexagon::Blade::Curve<2>& curve, const double* xy, double r, double* sol,
                                 double t0 = 0.0, double t1 = 0.0, double* tsol = 0, double* tanv = 0);

DLLEXPORT int lineIntersection(const Hexagon::Blade::Curve<2>& curve, const double* xy, const double* ij, double* sol,
                               double t0 = 0.0, double t1 = 0.0, double* tsol = 0, double* dist = 0, bool isRay = false,
                               bool useSol = false, const ptrdiff_t numberOfSeeds = 500);

DLLEXPORT void findClosestTValues(const Hexagon::Blade::Curve<2>& curve, double* outT, const double* inPoints,
                                  ptrdiff_t numberOfPoints);

DLLEXPORT double wrapToAbove(const double t, const double rangeStart, const double wholePeriod);

DLLEXPORT bool tIsInSubcurve(const double t, const Hexagon::Blade::Curve<2>& subCurve, const double wholePeriod);
} // namespace Blade
} // namespace Hexagon
