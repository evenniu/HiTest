#pragma once

#include "AbstractCurve.h"
#include <Eigen/Dense>
#include <memory>

namespace Hexagon
{
namespace Blade
{
struct FitOptions
{
  bool allowTranslation;
  bool allowRotation;
  std::shared_ptr<const Curve<1>> innerTolerance; // negatively-valued curve in most cases
  std::shared_ptr<const Curve<1>> outerTolerance; // positively-valued curve in most cases
  std::shared_ptr<const Curve<1>> weightTargetCurve;
  Eigen::VectorXd weightFittedPoints;

  // only used when allowRotation == true
  Eigen::VectorXd rotationLimits; // should be in radians

  // only used when allowTranslation == true
  Eigen::VectorXd translationXLimits;
  Eigen::VectorXd translationYLimits;

  // only used when allowTranslation == false
  Eigen::VectorXd targetCurvePivotPoint;
  Eigen::VectorXd pointsPivotPoint;

  FitOptions();
};

struct LinearDeviation
{
  Eigen::VectorXd nominalPoint;
  Eigen::VectorXd nominalTangentDirection;
  Eigen::VectorXd measuredPoint;

  LinearDeviation(const Eigen::Ref<const Eigen::Vector2d>& nominalPoint,
                  const Eigen::Ref<const Eigen::Vector2d>& nominalTangentDirection,
                  const Eigen::Ref<const Eigen::Vector2d>& measuredPoint)
    : nominalPoint(nominalPoint), nominalTangentDirection(nominalTangentDirection), measuredPoint(measuredPoint)
  {
  }
};

// returns a transformation, which when applied to the fitted points, aligns them to the target points
// specifically, fittedPoint1 becomes exactly aligned to targetPoint1,
// and the (fittedPoint2 - fittedPoint1) vector becomes aligned to (targetPoint2 - targetPoint1)
Eigen::Isometry2d twoPointBestFit(const Eigen::Ref<const Eigen::Vector2d>& targetPoint1,
                                  const Eigen::Ref<const Eigen::Vector2d>& targetPoint2,
                                  const Eigen::Ref<const Eigen::Vector2d>& fittedPoint1,
                                  const Eigen::Ref<const Eigen::Vector2d>& fittedPoint2);

// returns a transformation, which when applied to the input points, aligns them to the curve
Eigen::Isometry2d computeLeastSquaresBestFit(const Curve<2>& targetCurve,
                                             const Eigen::Ref<const Eigen::Matrix2Xd>& points,
                                             const Eigen::Isometry2d& guess, const FitOptions& options,
                                             const std::vector<LinearDeviation>& linearDeviations,
                                             const double inchSize);
Eigen::Isometry2d computeLeastSquaresBestFit(const std::vector<const Curve<2>*>& targetCurves,
                                             const std::vector<Eigen::Ref<const Eigen::Matrix2Xd>>& pointSets,
                                             const Eigen::Isometry2d& guess, const std::vector<FitOptions>& options,
                                             const std::vector<LinearDeviation>& linearDeviations,
                                             const double inchSize);

// returns a transformation, which when applied to the input points, aligns them to the curve
Eigen::Isometry2d computeMinMaxBestFit(const Curve<2>& targetCurve, const Eigen::Ref<const Eigen::Matrix2Xd>& points,
                                       const Eigen::Isometry2d& guess, const FitOptions& options,
                                       const std::vector<LinearDeviation>& linearDeviations, const double inchSize);
} // namespace Blade
} // namespace Hexagon
