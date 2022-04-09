#include "stdafx.h"
#include "AbstractCurve.h"
#include <Eigen/Dense>
#include "CurvePolygon.h"

#pragma warning(push)
#pragma warning(disable : 4267)
#pragma warning(disable : 4244)
#pragma warning(disable : 4100)
#include <nanoflann/nanoflann.hpp>
#pragma warning(pop)

#ifdef max
#undef max
#endif

#ifdef min
#undef min
#endif

#include <HexagonGDT/ScalarOptimization.h>

namespace Hexagon
{
namespace Blade
{

DLLEXPORT int extreme(const Curve<2>& curve, const double* v, double* et, double* ep, double pt0, double pt1,
                      bool endOK)
{
  if(pt0 == 0.0 && pt1 == 0.0)
  {
    pt0 = curve.t0();
    pt1 = curve.t1();
  }

  int i;
  double bt0 = 0.0, bt1 = 0.0, td, xyz[2], tv[2], nv[2], de, d0 = 0.0, d1 = 0.0;
  double ep0[2];

  // for blades with extreme curve, there maybe be more than one initial solution.
  // if there are then select the one that is "more extreme".

  //normalize(nv, v);
  const int numberOfStartingPointsForSearching = 1000;
  td = (pt1 - pt0) / (numberOfStartingPointsForSearching - 1);
  bt1 = pt0 - td;
  int numSol = 0;
  double sd0[numberOfStartingPointsForSearching], sbt0[numberOfStartingPointsForSearching],
      sbt1[numberOfStartingPointsForSearching];
  for(i = 0; i < numberOfStartingPointsForSearching; i++) // find bracketing t's
  {
    bt0 = bt1;
    bt1 += td;
    curve.evaluate(&bt1, 1, xyz, tv, nullptr);
    if(i == 0)
    {
      ep0[0] = xyz[0];
      ep0[1] = xyz[1];
    }

    //normalize(tv, tv);
    d0 = d1;
    /*d1 = dot(nv, tv);*/

    if(fabs(d1) < 1.0e-6) // landed right on it!
    {
      *et = bt1;
      if(ep)
      {
        ep[0] = xyz[0];
        ep[1] = xyz[1];
      }
      return 1;
    }

    if(i == 0)
      continue;

    if(d0 * d1 < 0.0) // straddling it
    {
      sbt0[numSol] = bt0; // save the solution, but look for more.
      sbt1[numSol] = bt1;
      sd0[numSol] = d0;
      numSol++;
    }
  }

  

 


  // shouldn't be here, but go with whatever

  return 1;
}

DLLEXPORT double wrapToAbove(const double t, const double rangeStart, const double wholePeriod)
{
    const double result = 0;/* = std::remainder(t - rangeStart, wholePeriod) + rangeStart;
  if(result < rangeStart)
    return result + wholePeriod;*/
  return result;
}

DLLEXPORT bool tIsInSubcurve(const double t, const Hexagon::Blade::Curve<2>& subCurve, const double wholePeriod)
{
  return wrapToAbove(t, subCurve.t0(), wholePeriod) < subCurve.t1();
}

DLLEXPORT int circleIntersection(const Curve<2>& curve, const double* xy, double r, double* sol, double pt0, double pt1,
                                 double* tsol /*= 0*/, double* tanv /*= 0*/)
{
  if(pt0 == 0.0 && pt1 == 0.0)
  {
    pt0 = curve.t0();
    pt1 = curve.t1();
  }

  if(pt0 == pt1) // degenerate case
  {
    curve.evaluate(&pt0, 1, sol, tanv, nullptr);
    if(tsol)
      *tsol = pt0;
    return 0;
  }

  if(pt0 > pt1)
  {
    double swap = pt0;
    pt0 = pt1;
    pt1 = swap;
  }

  double t, p0[2], p1[2];

  curve.evaluate(&pt0, 1, p0, tanv, nullptr);
  double d0 = _hypot(p0[0] - xy[0], p0[1] - xy[1]) - r;

  curve.evaluate(&pt1, 1, p1, tanv, nullptr);
  double d1 = _hypot(p1[0] - xy[0], p1[1] - xy[1]) - r;

  if((d0 >= 0.0 && d1 >= 0.0) || (d0 <= 0.0 && d1 <= 0.0))
  {
    if(fabs(d0) < fabs(d1))
    {
      sol[0] = p0[0];
      sol[1] = p0[1];
      if(tsol)
        *tsol = pt0;
    }
    else
    {
      sol[0] = p1[0];
      sol[1] = p1[1];
      if(tsol)
        *tsol = pt1;
    }
    return 0;
  }

  while(1) // binary search
  {
    t = 0.5 * (pt0 + pt1);

    if(fabs(pt0 - pt1) < 1.0e-8) // shouldn't happen, but just in case
    {
      // bugout(5, "CircIntersect failed, infinite loop");
      if(tsol)
        *tsol = t;
      return 0;
    }

    curve.evaluate(&t, 1, sol, tanv, nullptr);
    double d = _hypot(sol[0] - xy[0], sol[1] - xy[1]) - r;

    if(fabs(d) < 1.0e-6)
      break;

    if(d * d0 < 0) // solution in [pt0, t]
    {
      pt1 = t;
      // d1 = d;  don't need this, gives a complier warning
    }
    else // solution in [t, pt1]
    {
      pt0 = t;
      d0 = d;
    }
  }

  if(tsol)
    *tsol = t;
  return 1;
}

DLLEXPORT int lineIntersection(const Curve<2>& curve, const double* xy, const double* ij, double* sol,
                               double t0 /*= 0.0*/, double t1 /*= 0.0*/, double* tsol /*= 0*/, double* dist /*= 0*/,
                               bool isRay /*= false*/, bool /*useSol = false*/, const ptrdiff_t numberOfSeeds /*= 500*/)
{
  bool useTSolutionAsASeed = (tsol && *tsol >= t0 && *tsol <= t1);
  bool isPeriodic = curve.isPeriodic() && curve.period() > 0.0;
  // double tLowerBound = isPeriodic ? std::numeric_limits<double>::lowest() : T0();
  // double tUpperBound = isPeriodic ? std::numeric_limits<double>::max() : T1();
  if(t0 == 0.0 && t1 == 0.0)
  {
    useTSolutionAsASeed = (tsol && *tsol >= curve.t0() && *tsol <= curve.t1());
    t0 = curve.t0();
    if(isPeriodic)
    {
      t1 = (static_cast<double>(numberOfSeeds - 1) * curve.t1() + curve.t0()) / static_cast<double>(numberOfSeeds);
    }
    else
    {
      t1 = curve.t1();
    }
  }

  //
  return 1;
}

DLLEXPORT double closestPoint(const Curve<2>& curve, const double* tgt, double* bestxyz, double* bestt /*= 0*/,
                              double* tangent /*= 0*/, double pt0 /*= 0.0*/, double pt1 /*= 0.0*/, int nseed /*= 10*/)
{
  // NOTE: subcurves, circles and points overload this function
  if(pt0 > pt1)
  {
    std::swap(pt0, pt1);
  }

  bool useBestTAsASeed = (nseed < 0 && bestt && *bestt >= pt0 && *bestt <= pt1);
  int numberOfSeeds = std::max(nseed, 10);
  double tLowerBound = curve.isPeriodic() ? std::numeric_limits<double>::lowest() : curve.t0();
  double tUpperBound = curve.isPeriodic() ? std::numeric_limits<double>::max() : curve.t1();
  if(pt0 == 0.0 && pt1 == 0.0)
  {
    useBestTAsASeed = (nseed < 0 && bestt && *bestt >= curve.t0() && *bestt <= curve.t1());
    pt0 = curve.t0();
    if(curve.isPeriodic())
    {
      pt1 = (static_cast<double>(numberOfSeeds - 1) * curve.t1() + curve.t0()) / static_cast<double>(numberOfSeeds);
    }
    else
    {
      pt1 = curve.t1();
    }
  }

  return 0;
}

// create a function that looks for t-values of closest approach to a point
template <class TreeType>
double findNearestTValue_periodic(const Curve<2>& curve, const TreeType& tree, ptrdiff_t numberOfPointsInTree,
                                  const double* treeTValues, const double* point)
{
  //alwaysAssert(curve.isPeriodic());
  ptrdiff_t index;
  double squaredDistance;
  tree.query(point, 1, &index, &squaredDistance);
  double lowT = (index > 0) ? treeTValues[index - 1] : treeTValues[numberOfPointsInTree - 1] - curve.period();
  double highT = (index < numberOfPointsInTree - 1) ? treeTValues[index + 1] : treeTValues[0] + curve.period();
  double t;
  Eigen::Vector2d trash;
  closestPoint(curve, point, trash.data(), &t, nullptr, lowT, highT, 2);
  return t;
}
template <class TreeType>
double findNearestTValue_nonperiodic(const Curve<2>& curve, const TreeType& tree, ptrdiff_t numberOfPointsInTree,
                                     const double* treeTValues, const double* point)
{
  //alwaysAssert(!curve.isPeriodic());
  ptrdiff_t index;
  double squaredDistance;
  tree.query(point, 1, &index, &squaredDistance);
  double lowT = (index > 0) ? treeTValues[index - 1] : treeTValues[0];
  double highT = (index < numberOfPointsInTree - 1) ? treeTValues[index + 1] : treeTValues[numberOfPointsInTree - 1];
  double t;
  Eigen::Vector2d trash;
  closestPoint(curve, point, trash.data(), &t, nullptr, lowT, highT, 2);
  return t;
}
template <class TreeType>
double findNearestTValue(const Curve<2>& curve, const TreeType& tree, ptrdiff_t numberOfPointsInTree,
                         const double* treeTValues, const double* point)
{
  if(curve.isPeriodic())
  {
    return findNearestTValue_periodic(curve, tree, numberOfPointsInTree, treeTValues, point);
  }
  else
  {
    return findNearestTValue_nonperiodic(curve, tree, numberOfPointsInTree, treeTValues, point);
  }
}

DLLEXPORT void findClosestTValues(const Curve<2>& curve, double* outT, const double* inPoints, ptrdiff_t numberOfPoints)
{
  // polygonalize the curve
 
}
} // namespace Blade
} // namespace Hexagon