#pragma once
#include "AbstractCurve.h"
#include <Eigen/Dense>
#include "ArraySlicing.h"

namespace Hexagon
{
namespace Blade
{

template <ptrdiff_t dimension>
Eigen::Vector2d parametricBounds(const Curve<dimension>& curve)
{
  return Eigen::Vector2d(curve.t0(), curve.t1());
}

template <ptrdiff_t dimension>
Eigen::Matrix<double, dimension, 1> evaluate(const Curve<dimension>& curve, const double& t)
{
  Eigen::Matrix<double, dimension, 1> evaluatedPoint;
  curve.evaluate(&t, 1, evaluatedPoint.data(), nullptr, nullptr);
  return evaluatedPoint;
}

template <ptrdiff_t dimension>
Eigen::Matrix<double, dimension, Eigen::Dynamic> evaluate(const Curve<dimension>& curve,
                                                          const Eigen::Ref<const Eigen::VectorXd>& t)
{
  const Eigen::VectorXd contiguousT = t;
  Eigen::Matrix<double, dimension, Eigen::Dynamic> evaluatedPoints(dimension, t.size());
  curve.evaluate(contiguousT.data(), t.size(), evaluatedPoints.data(), nullptr, nullptr);
  return evaluatedPoints;
}

template <ptrdiff_t dimension>
std::pair<Eigen::Matrix<double, dimension, Eigen::Dynamic>, Eigen::Matrix<double, dimension, Eigen::Dynamic>>
evaluateWithDerivative(const Curve<dimension>& curve, const Eigen::Ref<const Eigen::VectorXd>& t)
{
  const Eigen::VectorXd contiguousT = t;
  Eigen::Matrix<double, dimension, Eigen::Dynamic> evaluatedPoints(dimension, t.size());
  Eigen::Matrix<double, dimension, Eigen::Dynamic> evaluatedDerivative(dimension, t.size());
  curve.evaluate(contiguousT.data(), t.size(), evaluatedPoints.data(), evaluatedDerivative.data(), nullptr);
  return std::make_pair(evaluatedPoints, evaluatedDerivative);
}

inline Eigen::VectorXd wrapToAbove(const Eigen::Ref<const Eigen::VectorXd>& t, const double rangeStart,
                                   const double wholePeriod)
{
  Eigen::ArrayXd result = t.unaryExpr(
      [=](const double value) -> double { return std::remainder(value - rangeStart, wholePeriod) + rangeStart; });
  result = (result < rangeStart).select(result + wholePeriod, result);
  return result;
}

template <ptrdiff_t dimension>
Eigen::ArrayXb tIsInSubcurve_eigen(const Eigen::Ref<const Eigen::VectorXd>& t,
                                   const Hexagon::Blade::Curve<dimension>& subCurve, const double wholePeriod)
{
  return wrapToAbove(t, subCurve.t0(), wholePeriod).array() < subCurve.t1();
}

Eigen::Isometry2d toIsometry2d();

CAlignment toCAlignment(const Eigen::Isometry2d& nominalToMeasuredTransformation);
} // namespace Blade
} // namespace Hexagon
