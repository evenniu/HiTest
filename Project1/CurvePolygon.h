#pragma once
#include "AbstractCurve.h"
#include <Eigen/Dense>
#include <tuple>
#include <array>
#include <vector>

namespace Hexagon
{
namespace Blade
{

Eigen::Matrix2Xd polygonalize(const Curve<2>& curve, const ptrdiff_t initialN, double targetError);
std::tuple<Eigen::Matrix2Xd, Eigen::VectorXd> polygonalizeWithT(const Curve<2>& curve, const ptrdiff_t initialN,
                                                                double targetError);

double signedArea(const Eigen::Ref<const Eigen::Matrix2Xd>& polygon);
double sign(double val);
double length(const Eigen::Ref<const Eigen::Matrix2Xd>& points);

Eigen::Vector2d centroid(const Eigen::Ref<const Eigen::Matrix2Xd>& polygon);

Eigen::Matrix2Xd convexHull(const Eigen::Ref<const Eigen::Matrix2Xd>& points);

struct VoronoiDiagram
{
  Eigen::Matrix2Xd points;
  Eigen::Matrix2Xd voronoiVertices;
  std::vector<std::array<ptrdiff_t, 3>> mapFromVoronoiVerticesToInputPoints;
};
VoronoiDiagram computeVoronoiDiagram(const Eigen::Ref<const Eigen::Matrix2Xd>& inputPoints);
VoronoiDiagram computeFarthestSiteVoronoiDiagram(const Eigen::Ref<const Eigen::Matrix2Xd>& inputPoints);

// note: gravity points "down" and so the surface normal(s) of the found contact-point(s)
// is(are) similar to the input gravity guess and equal to the output gravity value
Eigen::Vector2d gravityDirection(const Eigen::Ref<const Eigen::Matrix2Xd>& polygon,
                                 const Eigen::Ref<const Eigen::Vector2d>& guessGravity);

struct ThicknessResult
{
  Eigen::Vector2d thicknessCenter;
  double thickness;
};
// this is the center and diameter of the max-inscribed circle within a blade section
ThicknessResult maxThickness(const Curve<2>& curve, const Curve<2>& convexCurve, const Curve<2>& concaveCurve,
                             const ptrdiff_t initialN, double targetError);
// this is the center and diameter of the two-point max thickness within a blade section
ThicknessResult maxThicknessTwoPoint(const Curve<2>& curve, const Curve<2>& convexCurve, const Curve<2>& concaveCurve,
                                     const ptrdiff_t initialN, double targetError);

struct Camber
{
  Eigen::Matrix2Xd camberPoints;
  Eigen::Matrix2Xd camberTangents;
  Eigen::VectorXd camberRadii;
  Eigen::Matrix3Xd closestTValues;
};
Camber camber(const Curve<2>& curve, const Curve<2>& convexCurve, const Curve<2>& concaveCurve,
              const ptrdiff_t initialN, double targetError, const Eigen::Ref<const Eigen::Vector2d>& halfSpace1Point,
              const Eigen::Ref<const Eigen::Vector2d>& halfSpace1Vector,
              const Eigen::Ref<const Eigen::Vector2d>& halfSpace2Point,
              const Eigen::Ref<const Eigen::Vector2d>& halfSpace2Vector);

template <class FunctionType>
std::tuple<double, double> findInterval(const FunctionType& f, double x0, double scale)
{
  if(!std::isfinite(x0) || !std::isfinite(scale) || scale <= 0.0)
  {
    throw std::logic_error("Invalid parameters.");
  }
  double a = x0 - scale;
  double b = x0 + scale;
  double fa = f(a);
  double fb = f(b);
  if(fa * fb < 0.0)
  {
    return std::make_tuple(a, b);
  }
  double f0 = f(x0);
  if(f0 * fa < 0.0)
  {
    return std::make_tuple(a, x0);
  }
  if(f0 * fb < 0.0)
  {
    return std::make_tuple(x0, b);
  }
  double a2 = x0 - 2.0 * scale;
  double fa2 = f(a2);
  if(fa2 * fa < 0.0)
  {
    return std::make_tuple(a2, a);
  }
  double b2 = x0 + 2.0 * scale;
  double fb2 = f(b2);
  if(fb2 * fb < 0.0)
  {
    return std::make_tuple(b, b2);
  }
  double a4 = x0 - 4.0 * scale;
  double fa4 = f(a4);
  if(fa4 * fa2 < 0.0)
  {
    return std::make_tuple(a4, a2);
  }
  double b4 = x0 + 4.0 * scale;
  double fb4 = f(b4);
  if(fb4 * fb2 < 0.0)
  {
    return std::make_tuple(b2, b4);
  }
  throw std::logic_error("Could not find interval.");
}

template <class FunctionType>
double findZero(const FunctionType& f, double a, double b, double tol = 1e-12)
{
  // this could be a fast function like Brent's method, but for now a binary search is adequate
  if(f(a) == 0.0)
    return a;
  if(f(b) == 0.0)
    return b;
  if(a >= b)
  {
    return std::numeric_limits<double>::quiet_NaN();
  }
  double fa = f(a);
  double fb = f(b);
  if(fa * fb >= 0.0)
  {
    return std::numeric_limits<double>::quiet_NaN();
  }
  while(b - a > tol * abs(abs(a) + abs(b)))
  {
    double c = 0.5 * (a + b);
    double fc = f(c);
    if(fc == 0.0)
      return c;
    if(fc * fa < 0.0)
    {
      b = c;
      fb = fc;
    }
    else
    {
      a = c;
      fa = fc;
    }
  };
  return 0.5 * (a + b);
}
}
}
