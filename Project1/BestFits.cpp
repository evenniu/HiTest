#pragma warning(disable : 4714) // sometimes Eigen generates a warning saying something marked __forceinline is not
                                // inlined; this warning is fine to ignore
#include "stdafx.h"
#include "BestFits.h"
#include "CurvePolygon.h"
#include "ArraySlicing.h"
#include "EigenAbstractCurve.h"
#include "TemplateHermiteSpline.h"
#include "HermiteCurve.h"
#include "SUBCURVE.H"
#include <HexagonGDT/FitProfile.h>
#include <boost/uuid/uuid.hpp>
#include <boost/uuid/detail/sha1.hpp>
#include <boost/uuid/uuid_generators.hpp>
#include <numeric>

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

namespace Hexagon
{
namespace Blade
{

const double infinity = std::numeric_limits<double>::infinity();

FitOptions::FitOptions() : allowTranslation(true), allowRotation(true)
{
  rotationLimits = Eigen::Vector2d(-infinity, infinity);
  translationXLimits = Eigen::Vector2d(-infinity, infinity);
  translationYLimits = Eigen::Vector2d(-infinity, infinity);
  targetCurvePivotPoint = Eigen::Vector2d(0.0, 0.0);
  pointsPivotPoint = Eigen::Vector2d(0.0, 0.0);
}

Eigen::Isometry2d twoPointBestFit(const Eigen::Ref<const Eigen::Vector2d>& targetPoint1,
                                  const Eigen::Ref<const Eigen::Vector2d>& targetPoint2,
                                  const Eigen::Ref<const Eigen::Vector2d>& fittedPoint1,
                                  const Eigen::Ref<const Eigen::Vector2d>& fittedPoint2)
{
  const Eigen::Vector2d fittedVector = fittedPoint2 - fittedPoint1;
  const Eigen::Vector2d targetVector = targetPoint2 - targetPoint1;
  const double theta = signedAngleBetween2DVectors(fittedVector.data(), targetVector.data());
  return Eigen::Translation2d(targetPoint1) * Eigen::Rotation2Dd(theta) * Eigen::Translation2d(-fittedPoint1);
}

struct ConstructedMiddleOfToleranceZone : public Curve<2>
{
  const Curve<2>& nominal;
  const Curve<1>& innerTolerance;
  const Curve<1>& outerTolerance;

  ConstructedMiddleOfToleranceZone(const Curve<2>& nominal, const Curve<1>& innerTolerance,
                                   const Curve<1>& outerTolerance)
    : nominal(nominal), innerTolerance(innerTolerance), outerTolerance(outerTolerance)
  {
    assert(nominal.isPeriodic() == innerTolerance.isPeriodic());
    const double equalityThreshold = 1e-14 + 1e-14 * (nominal.t1() - nominal.t0());
    const auto essentiallyEqual = [=](double x, double y) { return std::abs(x - y) < equalityThreshold; };
    assert(essentiallyEqual(nominal.t0(), innerTolerance.t0()));
    assert(essentiallyEqual(nominal.t1(), innerTolerance.t1()));

    assert(nominal.isPeriodic() == outerTolerance.isPeriodic());
    assert(essentiallyEqual(nominal.t0(), outerTolerance.t0()));
    assert(essentiallyEqual(nominal.t1(), outerTolerance.t1()));

    // make sure the curve is counterclockwise
#ifndef NDEBUG
    auto polygon = polygonalize(nominal, 2048, 1e-3);
#endif
    assert(signedArea(polygon) > 0.0);
  }

  // Inherited via Curve
  virtual bool isPeriodic() const override
  {
    return nominal.isPeriodic();
  }
  virtual double t0() const override
  {
    return nominal.t0();
  }
  virtual double t1() const override
  {
    return nominal.t1();
  }
  virtual double period() const override
  {
    return nominal.period();
  }
  virtual void parametricBounds(double* const outBounds) const override
  {
    nominal.parametricBounds(outBounds);
  }
  virtual void evaluate(const double* const t, const ptrdiff_t n, double* const outPoints, double* const outDerivative,
                        double* const outSecondDerivative) const override
  {
    //alwaysAssert(outSecondDerivative == nullptr);

    // evaluate the nominal curve
    Eigen::Matrix2Xd nominalPoints(2, n);
    Eigen::Matrix2Xd nominalDerivatives(2, n);
    Eigen::Matrix2Xd nominalSecondDerivatives(2, n);
    nominal.evaluate(t, n, nominalPoints.data(), nominalDerivatives.data(), nominalSecondDerivatives.data());

    // construct the surface normals of the nominal curve
    Eigen::Matrix2Xd surfaceNormals(2, n);
    surfaceNormals.row(0) = nominalDerivatives.row(1);
    surfaceNormals.row(1) = -nominalDerivatives.row(0);
    const Eigen::VectorXd derivativeNorms = surfaceNormals.colwise().norm().transpose();
    surfaceNormals = surfaceNormals.array().rowwise() / derivativeNorms.array().transpose();

    // construct the derivatives of the surface normals
    const Eigen::VectorXd derivativeProductTerm =
        (nominalDerivatives.array() * nominalSecondDerivatives.array()).colwise().sum().transpose() /
        derivativeNorms.cwiseAbs2().array();
    Eigen::Matrix2Xd surfaceNormalDerivatives(2, n);
    surfaceNormalDerivatives.row(0) = nominalSecondDerivatives.row(1).cwiseQuotient(derivativeNorms.transpose()) -
                                      surfaceNormals.row(0).cwiseProduct(derivativeProductTerm.transpose());
    surfaceNormalDerivatives.row(1) = -nominalSecondDerivatives.row(0).cwiseQuotient(derivativeNorms.transpose()) -
                                      surfaceNormals.row(1).cwiseProduct(derivativeProductTerm.transpose());

    // evaluate the inner and outer tolerances
    Eigen::VectorXd innerTolerances(n);
    Eigen::VectorXd innerToleranceDerivatives(n);
    innerTolerance.evaluate(t, n, innerTolerances.data(), innerToleranceDerivatives.data(), nullptr);
    Eigen::VectorXd outerTolerances(n);
    Eigen::VectorXd outerToleranceDerivatives(n);
    outerTolerance.evaluate(t, n, outerTolerances.data(), outerToleranceDerivatives.data(), nullptr);

    // construct the results
    Eigen::Map<Eigen::Matrix2Xd> mappedOutPoints(outPoints, 2, n);
    mappedOutPoints = nominalPoints.array() + 0.5 * (surfaceNormals.array().rowwise() *
                                                     (outerTolerances + innerTolerances).array().transpose());
    if(outDerivative)
    {
      Eigen::Map<Eigen::Matrix2Xd> mappedOutDerivatives(outDerivative, 2, n);
      mappedOutDerivatives =
          nominalDerivatives.array() +
          0.5 * (surfaceNormals.array().rowwise() *
                     (outerToleranceDerivatives + innerToleranceDerivatives).array().transpose() +
                 surfaceNormalDerivatives.array().rowwise() * (outerTolerances + innerTolerances).array().transpose());
    }
  }
};

struct ConstructedToleranceZoneWeight : public Curve<1>
{
  const Curve<1>* originalWeight;
  const Curve<1>* innerTolerance;
  const Curve<1>* outerTolerance;
  const bool m_isPeriodic;
  const double m_t0, m_t1;

  ConstructedToleranceZoneWeight(const Curve<1>& originalWeight, const Curve<1>& innerTolerance,
                                 const Curve<1>& outerTolerance)
    : originalWeight(&originalWeight), innerTolerance(&innerTolerance), outerTolerance(&outerTolerance),
      m_isPeriodic(originalWeight.isPeriodic()), m_t0(originalWeight.t0()), m_t1(originalWeight.t1())
  {
    assert(originalWeight.isPeriodic() == innerTolerance.isPeriodic());
    assert(originalWeight.t0() == innerTolerance.t0());
    assert(originalWeight.t1() == innerTolerance.t1());

    assert(originalWeight.isPeriodic() == outerTolerance.isPeriodic());
    assert(originalWeight.t0() == outerTolerance.t0());
    assert(originalWeight.t1() == outerTolerance.t1());
  }

  ConstructedToleranceZoneWeight(const Curve<1>& innerTolerance, const Curve<1>& outerTolerance)
    : originalWeight(nullptr), innerTolerance(&innerTolerance), outerTolerance(&outerTolerance),
      m_isPeriodic(innerTolerance.isPeriodic()), m_t0(innerTolerance.t0()), m_t1(innerTolerance.t1())
  {
    assert(innerTolerance.isPeriodic() == outerTolerance.isPeriodic());
    assert(innerTolerance.t0() == outerTolerance.t0());
    assert(innerTolerance.t1() == outerTolerance.t1());
  }

  ConstructedToleranceZoneWeight(const Curve<1>& originalWeight)
    : originalWeight(&originalWeight), innerTolerance(nullptr), outerTolerance(nullptr),
      m_isPeriodic(originalWeight.isPeriodic()), m_t0(originalWeight.t0()), m_t1(originalWeight.t1())
  {
  }

  // Inherited via Curve
  virtual bool isPeriodic() const override
  {
    return m_isPeriodic;
  }
  virtual double t0() const override
  {
    return m_t0;
  }
  virtual double t1() const override
  {
    return m_t1;
  }
  virtual double period() const override
  {
    return m_t1 - m_t0;
  }
  virtual void parametricBounds(double* const outBounds) const override
  {
    outBounds[0] = m_t0;
    outBounds[1] = m_t1;
  }
  virtual void evaluate(const double* const t, const ptrdiff_t n, double* const outPoints, double* const outDerivative,
                        double* const outSecondDerivative) const override
  {
    Eigen::Map<Eigen::VectorXd> mappedPoints(outPoints, n);
    Eigen::Map<Eigen::VectorXd> mappedDerivative(outDerivative, n);
    Eigen::Map<Eigen::VectorXd> mappedSecondDerivative(outSecondDerivative, n);

    if(originalWeight && innerTolerance && outerTolerance)
    {
      // evaluate the inner and outer tolerances
      Eigen::VectorXd innerTolerances(n);
      Eigen::VectorXd innerToleranceDerivatives(n);
      Eigen::VectorXd innerToleranceSecondDerivatives(n);
      innerTolerance->evaluate(t, n, innerTolerances.data(), innerToleranceDerivatives.data(),
                               innerToleranceSecondDerivatives.data());
      Eigen::VectorXd outerTolerances(n);
      Eigen::VectorXd outerToleranceDerivatives(n);
      Eigen::VectorXd outerToleranceSecondDerivatives(n);
      outerTolerance->evaluate(t, n, outerTolerances.data(), outerToleranceDerivatives.data(),
                               outerToleranceSecondDerivatives.data());
      assert(((outerTolerances - innerTolerances).array() > 1e-6).all());

      Eigen::VectorXd originalWeights(n);
      Eigen::VectorXd originalWeightDerivatives(n);
      Eigen::VectorXd originalWeightSecondDerivatives(n);
      originalWeight->evaluate(t, n, originalWeights.data(), originalWeightDerivatives.data(),
                               originalWeightSecondDerivatives.data());

      const Eigen::VectorXd toleranceWeight = (outerTolerances - innerTolerances).cwiseInverse();
      const double sumWeights = originalWeights.sum();
      if(sumWeights == 0.0)
      {
        mappedPoints = toleranceWeight;
      }
      else
      {
        mappedPoints = (originalWeights / sumWeights).cwiseProduct(toleranceWeight);
      }
      if(!mappedPoints.allFinite())
      {
        THROW_IMPOSSIBLE_ERROR;
      }
      if(outDerivative)
      {
        const Eigen::VectorXd toleranceWeightDerivatives =
            -toleranceWeight.cwiseAbs2().cwiseProduct(outerToleranceDerivatives - innerToleranceDerivatives);
        if(sumWeights == 0.0)
        {
          mappedDerivative = toleranceWeightDerivatives;
        }
        else
        {
          const double sumWeights_dt = originalWeightDerivatives.sum();
          const Eigen::VectorXd normalizedOriginalWeights_dt =
              originalWeightDerivatives / sumWeights - originalWeights * (sumWeights_dt / (sumWeights * sumWeights));
          mappedDerivative = normalizedOriginalWeights_dt.cwiseProduct(toleranceWeight) +
                             (originalWeights / sumWeights).cwiseProduct(toleranceWeightDerivatives);
        }
      }
      if(outSecondDerivative)
      {
        mappedSecondDerivative.setConstant(std::numeric_limits<double>::quiet_NaN());
      }
    }
    else if(!originalWeight && innerTolerance && outerTolerance)
    {
      // evaluate the inner and outer tolerances
      Eigen::VectorXd innerTolerances(n);
      Eigen::VectorXd innerToleranceDerivatives(n);
      Eigen::VectorXd innerToleranceSecondDerivatives(n);
      innerTolerance->evaluate(t, n, innerTolerances.data(), innerToleranceDerivatives.data(),
                               innerToleranceSecondDerivatives.data());
      Eigen::VectorXd outerTolerances(n);
      Eigen::VectorXd outerToleranceDerivatives(n);
      Eigen::VectorXd outerToleranceSecondDerivatives(n);
      outerTolerance->evaluate(t, n, outerTolerances.data(), outerToleranceDerivatives.data(),
                               outerToleranceSecondDerivatives.data());
      assert(((outerTolerances - innerTolerances).array() > 1e-6).all());

      mappedPoints = (outerTolerances - innerTolerances).cwiseInverse();
      if(!mappedPoints.allFinite())
      {
        THROW_IMPOSSIBLE_ERROR;
      }
      if(outDerivative)
      {
        mappedDerivative =
            -mappedPoints.cwiseAbs2().cwiseProduct(outerToleranceDerivatives - innerToleranceDerivatives);
      }
      if(outSecondDerivative)
      {
        mappedSecondDerivative.setConstant(std::numeric_limits<double>::quiet_NaN());
      }
    }
    else if(originalWeight && !innerTolerance && !outerTolerance)
    {
      Eigen::VectorXd originalWeights(n);
      Eigen::VectorXd originalWeightDerivatives(n);
      Eigen::VectorXd originalWeightSecondDerivatives(n);
      originalWeight->evaluate(t, n, originalWeights.data(), originalWeightDerivatives.data(),
                               originalWeightSecondDerivatives.data());

      double sumWeights = originalWeights.sum();
      if(sumWeights == 0.0)
      {
        mappedPoints.setConstant(1.0);
      }
      else
      {
        mappedPoints = originalWeights / sumWeights;
      }
      if(!mappedPoints.allFinite())
      {
        THROW_IMPOSSIBLE_ERROR;
      }
      if(outDerivative)
      {
        if(sumWeights == 0.0)
        {
          mappedDerivative.setZero();
        }
        else
        {
          const double sumWeights_dt = originalWeightDerivatives.sum();
          mappedDerivative =
              originalWeightDerivatives / sumWeights - originalWeights * (sumWeights_dt / (sumWeights * sumWeights));
        }
      }
      if(outSecondDerivative)
      {
        mappedSecondDerivative.setConstant(std::numeric_limits<double>::quiet_NaN());
      }
    }
    else
    {
      THROW_IMPOSSIBLE_ERROR
    }
  }
};

struct Curve1DTranslator : public Hexagon::MetrologyBuildingBlocks::Curve1D
{
  const Curve<1>& input;

  Curve1DTranslator(const Curve<1>& input) : input(input)
  {
  }

  // Inherited via Curve1D
  virtual bool isPeriodic() const override
  {
    return input.isPeriodic();
  }
  virtual void parametricBounds(double* out) const override
  {
    input.parametricBounds(out);
  }
  virtual void evaluate(ptrdiff_t numberOfPoints, const double* t, double* out) const override
  {
    Eigen::VectorXd values(numberOfPoints);
    Eigen::VectorXd derivatives(numberOfPoints);
    Eigen::VectorXd secondDerivatives(numberOfPoints);
    input.evaluate(t, numberOfPoints, values.data(), derivatives.data(), secondDerivatives.data());
    Eigen::Map<Eigen::Matrix3Xd> mappedOut(out, 3, numberOfPoints);
    mappedOut.row(0) = values.transpose();
    mappedOut.row(1) = derivatives.transpose();
    mappedOut.row(2) = secondDerivatives.transpose();
  }
};

struct ConstantCurve1DTranslator : public Hexagon::MetrologyBuildingBlocks::Curve1D
{
  const double m_t0;
  const double m_t1;
  const bool m_isPeriodic;
  const double m_value;

  ConstantCurve1DTranslator(double t0, double t1, bool isPeriodic, double value)
    : m_t0(t0), m_t1(t1), m_isPeriodic(isPeriodic), m_value(value)
  {
  }

  // Inherited via Curve1D
  virtual bool isPeriodic() const override
  {
    return m_isPeriodic;
  }
  virtual void parametricBounds(double* out) const override
  {
    out[0] = m_t0;
    out[1] = m_t1;
  }
  virtual void evaluate(ptrdiff_t numberOfPoints, const double* /*t*/, double* out) const override
  {
    Eigen::Map<Eigen::Matrix3Xd> mappedOut(out, 3, numberOfPoints);
    mappedOut.row(0).setConstant(m_value);
    mappedOut.bottomRows<2>().setZero();
  }
};

class Curve2DTranslator : public Hexagon::FeatureCommunication::SearchableCurve2D
{
  typedef nanoflann::KDTreeEigenMatrixAdaptor<Eigen::MatrixX2d, 2> KDTree;

  const Curve<2>& input;
  Eigen::VectorXd curveT;
  Eigen::MatrixX2d curvePoints;
  std::shared_ptr<const KDTree> kdtree;
  std::shared_ptr<const Curve<2>> resampledCurve;
  const Curve<2>* splineForEvaluating;

public:
  Curve2DTranslator(const Curve<2>& input) : input(input)
  {
    const auto curvePolygon = polygonalizeWithT(input, 2048, 1e-4);
    curveT = std::get<1>(curvePolygon);
    curvePoints = std::get<0>(curvePolygon).transpose();
    kdtree = std::make_shared<KDTree>(2 /*dimensions*/, curvePoints, 10 /*max leaf*/);

    Eigen::Matrix<double, 5, Eigen::Dynamic> txyij(5, curveT.size());
    txyij.row(0) = curveT;
    txyij.middleRows<2>(1) = curvePoints.transpose();
    txyij.bottomRows<2>() = evaluateWithDerivative(input, curveT).second;
/*    const bool isNURBS = dynamic_cast<const CNurbCurve*>(&input) != nullptr ||
                         (dynamic_cast<const CSubCurve*>(&input) != nullptr &&
                          dynamic_cast<const CNurbCurve*>(dynamic_cast<const CSubCurve&>(input).m_psb))*/;
    const bool isNURBS = false;
    const bool isPeriodicHermite = dynamic_cast<const Hexagon::Blade::HermiteCurve<true, 2>*>(&input) != nullptr;
    const bool isNonperiodicHermite = dynamic_cast<const Hexagon::Blade::HermiteCurve<false, 2>*>(&input) != nullptr;
    const bool isHermiteOpenCCurve = dynamic_cast<const Hexagon::Blade::HermiteOpenCurve*>(&input) != nullptr;
    if(!isNURBS && !isPeriodicHermite && !isNonperiodicHermite && !isHermiteOpenCCurve &&
       std::isfinite(input.t1() - input.t0()))
    {
      if(input.isPeriodic())
      {
        txyij.conservativeResize(5, curveT.size() + 1);
        txyij.col(curveT.size()) = txyij.col(0);
        txyij(0, curveT.size()) = curveT[0] + input.period();
        resampledCurve = std::make_shared<const HermiteCurve<true, 2>>(txyij);
      }
      else
      {
        resampledCurve = std::make_shared<const HermiteCurve<false, 2>>(txyij);
      }
      //alwaysAssert(std::abs(resampledCurve->t0() - input.t0()) < 1e-10);
      //alwaysAssert(std::abs(resampledCurve->t1() - input.t1()) < 1e-10);
      splineForEvaluating = resampledCurve.get();
    }
    else
    {
      splineForEvaluating = &input;
    }
  }

  // Inherited via Curve2D
  virtual bool isPeriodic() const override
  {
    return input.isPeriodic();
  }
  virtual void parametricBounds(double* out) const override
  {
    input.parametricBounds(out);
  }
  virtual void evaluate(ptrdiff_t numberOfPoints, const double* t, double* out, double* outDerivative,
                        double* outSecondDerivative) const override
  {
    splineForEvaluating->evaluate(t, numberOfPoints, out, outDerivative, outSecondDerivative);
  }
  void findClosestTValues_periodic(ptrdiff_t numberOfPoints, const double* targetPoints, double* outTValues) const
  {
    Eigen::Map<const Eigen::Array2Xd> targets(targetPoints, 2, numberOfPoints);
    //alwaysAssert(input.isPeriodic());

    // search for each one
    for(int j = 0; j < numberOfPoints; j++)
    {
      const double* point = targetPoints + 2 * j;
      ptrdiff_t index;
      double squaredDistance;
      kdtree->query(point, 1, &index, &squaredDistance);
      const double lowT = (index > 0) ? curveT[index - 1] : curveT[curveT.size() - 1] - input.period();
      const double highT = (index < curveT.size() - 1) ? curveT[index + 1] : curveT[0] + input.period();
      Eigen::Vector2d trash;
      Hexagon::Blade::closestPoint(input, point, trash.data(), outTValues + j, nullptr, lowT, highT, 2);
    }
  }
  void findClosestTValues_nonperiodic(ptrdiff_t numberOfPoints, const double* targetPoints, double* outTValues) const
  {
    Eigen::Map<const Eigen::Array2Xd> targets(targetPoints, 2, numberOfPoints);
   // alwaysAssert(!input.isPeriodic());
   // alwaysAssert(input.t0() < input.t1());

    // search for each one
    for(int j = 0; j < numberOfPoints; j++)
    {
      const double* point = targetPoints + 2 * j;
      ptrdiff_t index;
      double squaredDistance;
      kdtree->query(point, 1, &index, &squaredDistance);
      const double lowT = (index > 0) ? curveT[index - 1] : input.t0();
      const double highT = (index < curveT.size() - 1) ? curveT[index + 1] : input.t1();
      Eigen::Vector2d trash;
      Hexagon::Blade::closestPoint(input, point, trash.data(), outTValues + j, nullptr, lowT, highT, 2);
    }
  }
  virtual void findClosestTValues(ptrdiff_t numberOfPoints, const double* targetPoints,
                                  double* outTValues) const override
  {
    if(input.isPeriodic())
    {
      findClosestTValues_periodic(numberOfPoints, targetPoints, outTValues);
    }
    else
    {
      findClosestTValues_nonperiodic(numberOfPoints, targetPoints, outTValues);
    }
  }
};

class LinearCurve : public Hexagon::FeatureCommunication::SearchableCurve2D
{
  Eigen::VectorXd nominalPoint;
  Eigen::VectorXd tangentDirection;

public:
  LinearCurve(const Eigen::Ref<const Eigen::Vector2d>& inNominalPoint,
              const Eigen::Ref<const Eigen::Vector2d>& inTangentDirection)
    : nominalPoint(inNominalPoint), tangentDirection(inTangentDirection.normalized())
  {
  }

  // Inherited via Curve2D
  virtual bool isPeriodic() const override
  {
    return false;
  }
  virtual void parametricBounds(double* out) const override
  {
    out[0] = -std::numeric_limits<double>::infinity();
    out[1] = std::numeric_limits<double>::infinity();
  }
  virtual void evaluate(ptrdiff_t numberOfPoints, const double* t, double* out, double* outDerivative,
                        double* outSecondDerivative) const override
  {
    Eigen::Map<const Eigen::VectorXd> mappedT(t, numberOfPoints);
    Eigen::Map<Eigen::Matrix2Xd> mappedOut(out, 2, numberOfPoints);
    Eigen::Map<Eigen::Matrix2Xd> mappedDerivative(outDerivative, 2, numberOfPoints);
    Eigen::Map<Eigen::Matrix2Xd> mappedSecondDerivative(outSecondDerivative, 2, numberOfPoints);

    mappedOut = nominalPoint.replicate(1, numberOfPoints) + tangentDirection * mappedT.transpose();
    mappedDerivative = tangentDirection.replicate(1, numberOfPoints);
    mappedSecondDerivative.setZero();
  }
  virtual void findClosestTValues(ptrdiff_t numberOfPoints, const double* targetPoints,
                                  double* outTValues) const override
  {
    Eigen::Map<const Eigen::Matrix2Xd> mappedTargetPoints(targetPoints, 2, numberOfPoints);
    Eigen::Map<Eigen::VectorXd> mappedTValues(outTValues, numberOfPoints);
    mappedTValues = (mappedTargetPoints.colwise() - nominalPoint).transpose() * tangentDirection;
  }
};

bool validateLicense(const void* input, void* output)
{
  // unfortunately we can't really check here to see if our license is valid...
  // that's done elsewhere and not available here

  const size_t uuidSize = boost::uuids::uuid::static_size();
  const static boost::uuids::uuid secretPassword = {
    82, 10, 235, 29, 69, 239, 203, 60, 45, 149, 216, 31, 7, 36, 196, 94
  };
  //{ 59, 219, 129, 158, 46, 28, 103, 251, 125, 64, 63, 64, 208, 12, 243, 122 }; // this is the secret password for
  // testing

  // ask the application to verify its license
  std::array<uint8_t, 2 * uuidSize> concatenatedData;
  std::copy_n(secretPassword.begin(), uuidSize, concatenatedData.begin());
  std::copy_n(static_cast<const uint8_t*>(input), uuidSize, concatenatedData.begin() + uuidSize);
  boost::uuids::detail::sha1 sha;
  sha.process_block(concatenatedData.data(), concatenatedData.data() + 2 * uuidSize);
  unsigned int response[5];
  sha.get_digest(response);
  for(unsigned int i = 0; i < 5; i++)
  {
    response[i] = _byteswap_ulong(response[i]);
  }
  auto typedOutput = stdext::make_checked_array_iterator<unsigned int*>(static_cast<unsigned int*>(output), 5);
  std::copy_n(response, 5, typedOutput);
  return true;
}

Eigen::Isometry2d computeBestFit(const std::vector<const Curve<2>*> targetCurves,
                                 const std::vector<Eigen::Ref<const Eigen::Matrix2Xd>>& unalignedPointsArrays,
                                 const Eigen::Isometry2d& guess, const std::vector<FitOptions>& optionsArray,
                                 const Hexagon::MetrologyBuildingBlocks::MathType mathType,
                                 const std::vector<LinearDeviation>& linearDeviations, double inchSize)
{
  const size_t numberOfCurves = targetCurves.size();
  if(numberOfCurves == 0)
  {
    throw std::logic_error("I'm not sure how to do a fit without any curves");
  }
  if(numberOfCurves != unalignedPointsArrays.size() || numberOfCurves != optionsArray.size())
  {
    throw std::logic_error("The number of target curves and the number of options structures must be the same");
  }
  const FitOptions& firstOptions = optionsArray.at(0);
  if(!firstOptions.allowTranslation && !firstOptions.allowRotation)
  {
    //throw std::logic_error("I'm not sure how to do a fit that doesn't allow translation and doesn't allow rotation");
      //don't throw exception, just do nothing
  }

  // what is the maximum number of points?
  const auto numberOfPointsComparison = [](const Eigen::Ref<const Eigen::Matrix2Xd>& a,
                                           const Eigen::Ref<const Eigen::Matrix2Xd>& b) { return a.cols() < b.cols(); };
  const ptrdiff_t maximumNumberOfPoints =
      std::max_element(unalignedPointsArrays.begin(), unalignedPointsArrays.end(), numberOfPointsComparison)->cols();
  const Eigen::Matrix2Xd zeros = Eigen::Matrix2Xd::Zero(2, std::max<ptrdiff_t>(10, maximumNumberOfPoints));

  // use the guess as an initial condition
  // also figure out which points have any weight at all
  // and the weight functions
  // and the middle of the tolerance zone
  // and the searchable curve
  // and the SurfaceData structures
  std::vector<Eigen::Matrix2Xd> pointsWorthUsingArrays;
  std::vector<Eigen::ArrayXd> weightsWorthUsingArrays;
  std::vector<std::shared_ptr<const ConstructedToleranceZoneWeight>> constructedWeightArray;
  std::vector<std::shared_ptr<const Hexagon::MetrologyBuildingBlocks::Curve1D>> weightTargetCurveArray;
  std::vector<std::shared_ptr<const Curve<2>>> constructedMiddleOfToleranceZoneArray;
  std::vector<std::shared_ptr<const FeatureCommunication::SearchableCurve2D>> middleOfToleranceZoneArray;
  std::vector<FeatureCommunication::SurfaceData> surfaceDataArray;
  ptrdiff_t totalNumberOfPoints = 0;
  for(size_t i = 0; i < numberOfCurves; i++)
  {
    if(optionsArray.at(i).weightFittedPoints.size() != unalignedPointsArrays.at(i).cols())
    {
      throw std::logic_error("I'm expecting the weights to be the same size as the number of points");
    }
    const Eigen::Matrix2Xd guessPoints = guess * unalignedPointsArrays.at(i);
    const Eigen::ArrayXb pointHasWeight = optionsArray.at(i).weightFittedPoints.array() != 0.0;
    pointsWorthUsingArrays.emplace_back(sliceColumns(guessPoints, pointHasWeight));
    weightsWorthUsingArrays.emplace_back(sliceVector(optionsArray.at(i).weightFittedPoints, pointHasWeight));

    // also figure out the weight functions and middle of the tolerance zone and the searchable curve
    std::shared_ptr<ConstructedToleranceZoneWeight> constructedWeight;
    std::shared_ptr<ConstructedMiddleOfToleranceZone> constructedMiddleOfToleranceZone;
    if(optionsArray.at(i).weightTargetCurve && optionsArray.at(i).innerTolerance && optionsArray.at(i).outerTolerance)
    {
      // construct the weight
      constructedWeight = std::make_shared<ConstructedToleranceZoneWeight>(*optionsArray.at(i).weightTargetCurve,
                                                                           *optionsArray.at(i).innerTolerance,
                                                                           *optionsArray.at(i).outerTolerance);
      constructedWeightArray.push_back(constructedWeight);
      weightTargetCurveArray.push_back(std::make_shared<Curve1DTranslator>(*constructedWeight));

      // construct the middle of the zone
      constructedMiddleOfToleranceZone = std::make_shared<ConstructedMiddleOfToleranceZone>(
          *targetCurves.at(i), *optionsArray.at(i).innerTolerance, *optionsArray.at(i).outerTolerance);
      constructedMiddleOfToleranceZoneArray.push_back(constructedMiddleOfToleranceZone);
      middleOfToleranceZoneArray.push_back(std::make_shared<Curve2DTranslator>(*constructedMiddleOfToleranceZone));
    }
    else if(optionsArray.at(i).innerTolerance && optionsArray.at(i).outerTolerance)
    {
      // construct the weight
      constructedWeight = std::make_shared<ConstructedToleranceZoneWeight>(*optionsArray.at(i).innerTolerance,
                                                                           *optionsArray.at(i).outerTolerance);
      constructedWeightArray.push_back(constructedWeight);
      weightTargetCurveArray.push_back(std::make_shared<Curve1DTranslator>(*constructedWeight));

      // construct the middle of the zone
      constructedMiddleOfToleranceZone = std::make_shared<ConstructedMiddleOfToleranceZone>(
          *targetCurves.at(i), *optionsArray.at(i).innerTolerance, *optionsArray.at(i).outerTolerance);
      constructedMiddleOfToleranceZoneArray.push_back(constructedMiddleOfToleranceZone);
      middleOfToleranceZoneArray.push_back(std::make_shared<Curve2DTranslator>(*constructedMiddleOfToleranceZone));
    }
    else if(optionsArray.at(i).weightTargetCurve)
    {
      constructedWeight = std::make_shared<ConstructedToleranceZoneWeight>(*optionsArray.at(i).weightTargetCurve);
      constructedWeightArray.push_back(constructedWeight);
      weightTargetCurveArray.push_back(std::make_shared<Curve1DTranslator>(*constructedWeight));
      middleOfToleranceZoneArray.push_back(std::make_shared<Curve2DTranslator>(*targetCurves.at(i)));
    }
    else
    {
      weightTargetCurveArray.push_back(nullptr);
      middleOfToleranceZoneArray.push_back(std::make_shared<Curve2DTranslator>(*targetCurves.at(i)));
    }

    // now the surface data structures
    Hexagon::FeatureCommunication::SurfaceData surfaceData;
    surfaceData.numberOfPoints = pointsWorthUsingArrays.back().cols();
    surfaceData.ballCenters = pointsWorthUsingArrays.back().data();
    surfaceData.hitVectors = zeros.data();
    surfaceData.tipRadii = zeros.data();
    surfaceDataArray.push_back(surfaceData);
    totalNumberOfPoints += surfaceData.numberOfPoints;
  }

  // if we don't have any points left, we're done
  if(totalNumberOfPoints == 0)
  {
    return guess;
  }

  // create the options structure for the fitting
  std::vector<const FeatureCommunication::SearchableCurve2D*> curvesForCalling;
  Hexagon::MetrologyBuildingBlocks::ProfileFittingOptions fittingOptions;
  fittingOptions.allowTranslation = firstOptions.allowTranslation;
  fittingOptions.allowRotation = firstOptions.allowRotation;
  fittingOptions.allowScaling = false;
  fittingOptions.dataPivotPoint = guess * firstOptions.pointsPivotPoint.head<2>();
  fittingOptions.nominalCurvePivotPoint = firstOptions.targetCurvePivotPoint;
  assert(fittingOptions.dataPivotPoint.allFinite());
  assert(fittingOptions.nominalCurvePivotPoint.allFinite());
  fittingOptions.mathType = mathType;
  fittingOptions.rotationLimits = firstOptions.rotationLimits;
  fittingOptions.translationXLimits = firstOptions.translationXLimits;
  fittingOptions.translationYLimits = firstOptions.translationYLimits;
  for(size_t i = 0; i < numberOfCurves; i++)
  {
    fittingOptions.weightFittedPointsArrays.emplace_back(weightsWorthUsingArrays.at(i));
    fittingOptions.weightNominalCurves.push_back(weightTargetCurveArray.at(i).get());
    curvesForCalling.push_back(middleOfToleranceZoneArray.at(i).get());
  }

  // put the linear deviations in
  std::vector<std::shared_ptr<const Eigen::VectorXd>> alignedDeviationPoints;
  ptrdiff_t numberOfLinearDeviations = 0;
  for(const auto& linearDeviation : linearDeviations)
  {
    // add to the weights and the curves
    numberOfLinearDeviations++;
    fittingOptions.weightFittedPointsArrays.emplace_back(Eigen::VectorXd::Ones(1));
    fittingOptions.weightNominalCurves.push_back(nullptr);
    middleOfToleranceZoneArray.push_back(
        std::make_shared<const LinearCurve>(linearDeviation.nominalPoint, linearDeviation.nominalTangentDirection));
    curvesForCalling.push_back(middleOfToleranceZoneArray.back().get());
    // add to the surface data array
    Hexagon::FeatureCommunication::SurfaceData surfaceData;
    surfaceData.numberOfPoints = 1;
    alignedDeviationPoints.emplace_back(new Eigen::VectorXd(guess * linearDeviation.measuredPoint.head<2>()));
    surfaceData.ballCenters = alignedDeviationPoints.back()->data();
    surfaceData.hitVectors = zeros.data();
    surfaceData.tipRadii = zeros.data();
    surfaceDataArray.push_back(surfaceData);
  }

  // figure out the initial deviations, so we can compute a reasonable guess for the linear deviation weights
  std::vector<Eigen::VectorXd> initialDeviations;
  ptrdiff_t totalNumberOfDeviations = 0;
  for(size_t i = 0; i < numberOfCurves; i++)
  {
    totalNumberOfDeviations += surfaceDataArray.at(i).numberOfPoints;
    initialDeviations.emplace_back(surfaceDataArray.at(i).numberOfPoints);
    Hexagon::MetrologyBuildingBlocks::computeDeviations(
        *curvesForCalling.at(i), Hexagon::MetrologyBuildingBlocks::Transform2D(Eigen::Isometry2d::Identity()),
        surfaceDataArray.at(i), initialDeviations.at(i));
  }
  Eigen::VectorXd allInitialDeviations(totalNumberOfDeviations), allWeights(totalNumberOfDeviations);
  ptrdiff_t k = 0;
  for(size_t i = 0; i < numberOfCurves; i++)
  {
    allInitialDeviations.segment(k, surfaceDataArray.at(i).numberOfPoints) = initialDeviations.at(i);
    allWeights.segment(k, surfaceDataArray.at(i).numberOfPoints) = fittingOptions.weightFittedPointsArrays.at(i);
    k += surfaceDataArray.at(i).numberOfPoints;
  }
  const double characteristicLinearDeviationSize = 0.0002 * inchSize;
  double linearDeviationWeight = 0.0;
  if(mathType == Hexagon::MetrologyBuildingBlocks::MathType::LEAST_SQUARES)
  {
    const double sumOfWeightedSquaredDeviations = allInitialDeviations.cwiseAbs2().cwiseProduct(allWeights).sum();
    const double sumOfSquaredLinearDeviations_lowerLimit =
        numberOfLinearDeviations * std::pow(characteristicLinearDeviationSize, 2.0);
    linearDeviationWeight = sumOfWeightedSquaredDeviations / sumOfSquaredLinearDeviations_lowerLimit;
  }
  else if(mathType == Hexagon::MetrologyBuildingBlocks::MathType::MIN_MAX)
  {
    const double maxOfWeightedDeviations = allInitialDeviations.cwiseAbs().cwiseProduct(allWeights).maxCoeff();
    const double maxOfLinearDeviations_lowerLimit = characteristicLinearDeviationSize;
    linearDeviationWeight = maxOfWeightedDeviations / maxOfLinearDeviations_lowerLimit;
  }
  for(size_t i = numberOfCurves; i < curvesForCalling.size(); i++)
  {
    fittingOptions.weightFittedPointsArrays.at(i).setConstant(linearDeviationWeight);
  }

  // restrict the rotation, translation, etc in certain cases
  if(firstOptions.allowTranslation)
  {
    // in the documentation, we say the pivot-points are ignored if translation is allowed
    // and so we can stomp on them all we like
    const auto appendMatrix = [](const Eigen::Matrix2Xd& a, const Eigen::Matrix2Xd& b) -> Eigen::Matrix2Xd {
      Eigen::Matrix2Xd result(2, a.cols() + b.cols());
      result << a, b;
      return result;
    };
    const Eigen::Matrix2Xd allPoints = std::accumulate(pointsWorthUsingArrays.begin(), pointsWorthUsingArrays.end(),
                                                       Eigen::Matrix2Xd(2, 0), appendMatrix);
    const Eigen::Vector2d meanPoint = allPoints.rowwise().mean();
    fittingOptions.dataPivotPoint = meanPoint;
    fittingOptions.nominalCurvePivotPoint = meanPoint;

    // merge in reasonable translation limits and rotation limits
    const double radiusFactor =
        (guess.matrix() - Eigen::Matrix3d::Identity()).cwiseAbs().maxCoeff() < 1e-8 ? 2.0 : 0.25;
    const double maxRadius = radiusFactor * (allPoints.colwise() - meanPoint).colwise().norm().maxCoeff();
    fittingOptions.translationXLimits[0] = std::max(firstOptions.translationXLimits[0], -maxRadius);
    fittingOptions.translationXLimits[1] = std::min(firstOptions.translationXLimits[1], maxRadius);
    fittingOptions.translationYLimits[0] = std::max(firstOptions.translationYLimits[0], -maxRadius);
    fittingOptions.translationYLimits[1] = std::min(firstOptions.translationYLimits[1], maxRadius);
    fittingOptions.rotationLimits[0] = std::max(firstOptions.rotationLimits[0], -M_PI_4);
    fittingOptions.rotationLimits[1] = std::min(firstOptions.rotationLimits[1], M_PI_4);
  }

  // make sure the weighting curves are either empty or filled
  const auto isSomething = [](const MetrologyBuildingBlocks::Curve1D* curve) -> bool { return curve != nullptr; };
  const bool allAreSomething =
      std::all_of(fittingOptions.weightNominalCurves.begin(), fittingOptions.weightNominalCurves.end(), isSomething);
  const bool noneAreSomething =
      std::none_of(fittingOptions.weightNominalCurves.begin(), fittingOptions.weightNominalCurves.end(), isSomething);
  if(noneAreSomething)
  {
    fittingOptions.weightNominalCurves.clear();
  }
  else if(!allAreSomething)
  {
    for(size_t i = 0; i < fittingOptions.weightNominalCurves.size(); i++)
    {
      if(fittingOptions.weightNominalCurves.at(i) == nullptr)
      {
        if(i < targetCurves.size())
        {
          weightTargetCurveArray.push_back(std::make_shared<const ConstantCurve1DTranslator>(
              targetCurves.at(i)->t0(), targetCurves.at(i)->t1(), targetCurves.at(i)->isPeriodic(), 1.0));
        }
        else
        {
          //alwaysAssert(std::dynamic_pointer_cast<const LinearCurve>(middleOfToleranceZoneArray.at(i)) != nullptr);
          weightTargetCurveArray.push_back(std::make_shared<const ConstantCurve1DTranslator>(
              -std::numeric_limits<double>::infinity(), std::numeric_limits<double>::infinity(), false, 1.0));
        }
        fittingOptions.weightNominalCurves[i] = weightTargetCurveArray.back().get();
      }
    }
  }

  // do the profile fit
  auto fittedTransformation = Hexagon::MetrologyBuildingBlocks::fit2DProfiles(curvesForCalling, surfaceDataArray,
                                                                              fittingOptions, validateLicense);
//  fittedTransformation.matrix
  const bool allFinite = fittedTransformation.matrix().allFinite();
  if(!allFinite)
  {
    throw std::logic_error("The fitting produced NaN in the result transformation");
  }

  // finish up, remembering to add the guess to the result because we started with the guess
  auto result = Eigen::Isometry2d::Identity();
  //result.matrix().
  result.matrix().topRows<2>() = fittedTransformation.matrix();
 //double num=  result(2,0);
 // num = result(0, 0);
 // num = result(0, 1);
 // num = result(0, 2);
 // num = result(1, 0);
 // num = result(1, 1);
 // num = result(1, 2);
  // thankfully, affine transforms are associative and so the final result is the fitted result times the guess
  return result * guess;
}

Eigen::Isometry2d computeLeastSquaresBestFit(const Curve<2>& targetCurve,
                                             const Eigen::Ref<const Eigen::Matrix2Xd>& points,
                                             const Eigen::Isometry2d& guess, const FitOptions& options,
                                             const std::vector<LinearDeviation>& linearDeviations,
                                             const double inchSize)
{
  return computeBestFit({ &targetCurve }, { points }, guess, { options },
                        MetrologyBuildingBlocks::MathType::LEAST_SQUARES, linearDeviations, inchSize);
}

Eigen::Isometry2d computeLeastSquaresBestFit(const std::vector<const Curve<2>*>& targetCurves,
                                             const std::vector<Eigen::Ref<const Eigen::Matrix2Xd>>& pointSets,
                                             const Eigen::Isometry2d& guess, const std::vector<FitOptions>& options,
                                             const std::vector<LinearDeviation>& linearDeviations,
                                             const double inchSize)
{
  return computeBestFit(targetCurves, pointSets, guess, options, MetrologyBuildingBlocks::MathType::LEAST_SQUARES,
                        linearDeviations, inchSize);
}

// returns a transformation, which when applied to the input points, aligns them to the curve
Eigen::Isometry2d computeMinMaxBestFit(const Curve<2>& targetCurve, const Eigen::Ref<const Eigen::Matrix2Xd>& points,
                                       const Eigen::Isometry2d& guess, const FitOptions& options,
                                       const std::vector<LinearDeviation>& linearDeviations,
                                       const double inchSize)
{
  return computeBestFit({ &targetCurve }, { points }, guess, { options }, MetrologyBuildingBlocks::MathType::MIN_MAX,
                        linearDeviations, inchSize);
}
} // namespace Blade
} // namespace Hexagon

#ifdef GOOGLETEST
#include <gtest/gtest.h>
#include "HermiteCurve.h"

namespace Hexagon
{
namespace Blade
{

template <ptrdiff_t dimension>
Eigen::Matrix<double, dimension, Eigen::Dynamic> computeFunction(const Curve<dimension>& function,
                                                                 const Eigen::Ref<const Eigen::VectorXd>& t)
{
  const Eigen::VectorXd tCopy = t;
  Eigen::Matrix<double, dimension, Eigen::Dynamic> f(dimension, t.size());
  function.evaluate(tCopy.data(), tCopy.size(), f.data(), nullptr, nullptr);
  return f;
}

template <ptrdiff_t dimension>
Eigen::Matrix<double, dimension, Eigen::Dynamic> computeDerivative(const Curve<dimension>& function,
                                                                   const Eigen::Ref<const Eigen::VectorXd>& t)
{
  const Eigen::VectorXd tCopy = t;
  Eigen::Matrix<double, dimension, Eigen::Dynamic> fTrash(dimension, t.size());
  Eigen::Matrix<double, dimension, Eigen::Dynamic> fDerivative(dimension, t.size());
  function.evaluate(tCopy.data(), tCopy.size(), fTrash.data(), fDerivative.data(), nullptr);
  return fDerivative;
}

template <ptrdiff_t dimension>
Eigen::Matrix<double, dimension, Eigen::Dynamic> computeSecondDerivative(const Curve<dimension>& function,
                                                                         const Eigen::Ref<const Eigen::VectorXd>& t)
{
  const Eigen::VectorXd tCopy = t;
  Eigen::Matrix<double, dimension, Eigen::Dynamic> fTrash(dimension, t.size());
  Eigen::Matrix<double, dimension, Eigen::Dynamic> fDerivativeTrash(dimension, t.size());
  Eigen::Matrix<double, dimension, Eigen::Dynamic> fSecondDerivative(dimension, t.size());
  function.evaluate(tCopy.data(), tCopy.size(), fTrash.data(), fDerivativeTrash.data(), fSecondDerivative.data());
  return fSecondDerivative;
}

template <ptrdiff_t dimension>
Eigen::Matrix<double, dimension, Eigen::Dynamic>
computeDerivativeNumerically(const Curve<dimension>& function, const Eigen::Ref<const Eigen::VectorXd>& t,
                             const double epsilon = 3e-6)
{
  auto y_p2 = computeFunction(function, t.array() + 2 * epsilon);
  auto y_p1 = computeFunction(function, t.array() + epsilon);
  auto y_m1 = computeFunction(function, t.array() - epsilon);
  auto y_m2 = computeFunction(function, t.array() - 2 * epsilon);
  return (-y_p2 + 8 * y_p1 - 8 * y_m1 + y_m2) / (12 * epsilon);
}

template <ptrdiff_t dimension>
Eigen::Matrix<double, dimension, Eigen::Dynamic>
computeSecondDerivativeNumerically(const Curve<dimension>& function, const Eigen::Ref<const Eigen::VectorXd>& t,
                                   const double epsilon = 3e-6)
{
  auto yPrime_p2 = computeDerivative(function, t.array() + 2 * epsilon);
  auto yPrime_p1 = computeDerivative(function, t.array() + epsilon);
  auto yPrime_m1 = computeDerivative(function, t.array() - epsilon);
  auto yPrime_m2 = computeDerivative(function, t.array() - 2 * epsilon);
  return (-yPrime_p2 + 8 * yPrime_p1 - 8 * yPrime_m1 + yPrime_m2) / (12 * epsilon);
}

template <ptrdiff_t dimension>
void testVectorFunctionFirstDerivative(const Curve<dimension>& function, const Eigen::Ref<const Eigen::VectorXd>& t,
                                       const double maximumAbsoluteDifference = 1e-6,
                                       const double maximumRelativeNorm = 1e-6)
{
  auto analytic = computeDerivative(function, t);
  auto numeric = computeDerivativeNumerically(function, t);

  EXPECT_GE(maximumAbsoluteDifference, (analytic - numeric).cwiseAbs().maxCoeff());
  EXPECT_GE(maximumRelativeNorm, (analytic - numeric).norm() / analytic.norm());
}

template <ptrdiff_t dimension>
void testVectorFunctionSecondDerivative(const Curve<dimension>& function, const Eigen::Ref<const Eigen::VectorXd>& t,
                                        const double maximumAbsoluteDifference = 1e-6,
                                        const double maximumRelativeNorm = 1e-6)
{
  auto analytic = computeSecondDerivative(function, t);
  auto numeric = computeSecondDerivativeNumerically(function, t);

  EXPECT_GE(maximumAbsoluteDifference, (analytic - numeric).cwiseAbs().maxCoeff());
  EXPECT_GE(maximumRelativeNorm, (analytic - numeric).norm() / analytic.norm());
}

template <ptrdiff_t dimension>
void testVectorFunctionDerivatives(const Curve<dimension>& function, const Eigen::Ref<const Eigen::VectorXd>& t)
{
  testVectorFunctionFirstDerivative(function, t);
  testVectorFunctionSecondDerivative(function, t);
}

TEST(BestFitDerivatives, HermiteCurve2D)
{
  // start with egg-shaped nominal curve
  Eigen::MatrixXd TXYIJ(5, 5);
  TXYIJ << 0.0, 1.0, 2.0, 4.0, 6.0, // t-values
      0.0, -1.0, 0.0, 2.0, 0.0,     // x-values
      1.0, 0.0, -1.0, 0.0, 1.0,     // y-values
      1.0, 0.0, -1.0, 0.0, 1.0,     // i-values
      0.0, 1.0, 0.0, -1.0, 0.0;     // j-values
  auto nominalCurve = createHermiteCurve2D(TXYIJ, true);

  // construct the list of test points
  // const std::vector<double> testT = { -2.0, -0.01, 0.0,  0.01, 0.5,  0.99, 1.0, 1.01, 1.99, 2.0,  2.01, 2.49, 2.5,
  //                                    2.51, 3.0,   3.99, 4.0,  4.01, 4.49, 4.5, 4.51, 5.0,  5.99, 6.0,  6.01, 9.0 };
  const std::vector<double> testT = { -1.99, -0.01, 0.01, 0.5,  0.99, 1.01, 1.99, 2.01, 2.49, 2.5, 2.51,
                                      3.0,   3.99,  4.01, 4.49, 4.5,  4.51, 5.0,  5.99, 6.01, 9.0 };

  // do the test
  testVectorFunctionDerivatives(*nominalCurve, Eigen::Map<const Eigen::VectorXd>(testT.data(), testT.size()));
}

TEST(BestFitDerivatives, ConstructedMiddleOfToleranceZoneDerivatives)
{
  // start with egg-shaped nominal curve
  Eigen::MatrixXd TXYIJ(5, 5);
  TXYIJ << 0.0, 1.0, 2.0, 4.0, 6.0, // t-values
      0.0, -1.0, 0.0, 2.0, 0.0,     // x-values
      1.0, 0.0, -1.0, 0.0, 1.0,     // y-values
      1.0, 0.0, -1.0, 0.0, 1.0,     // i-values
      0.0, 1.0, 0.0, -1.0, 0.0;     // j-values
  auto nominalCurve = createHermiteCurve2D(TXYIJ, true);

  // create inner and outer tolerance curves
  Eigen::Matrix3Xd innerTXI(3, 4);
  innerTXI << 0.0, 2.5, 4.5, 6.0, // t-values
      0.1, 0.2, 0.3, 0.1,         // tolerance values
      0.0, 0.1, 0.0, 0.0;         // tolerance derivative values
  auto innerToleranceCurve = createHermiteCurve1D(innerTXI, true);
  Eigen::Matrix3Xd outerTXI(3, 2);
  outerTXI << 0.0, 6.0, // t-values
      0.4, 0.4,         // tolerance values
      0.0, 0.0;         // tolerance derivative values
  auto outerToleranceCurve = createHermiteCurve1D(outerTXI, true);

  // construct the middle-of-tolerance-zone test curve
  ConstructedMiddleOfToleranceZone testCurve(*nominalCurve, *innerToleranceCurve, *outerToleranceCurve);

  // construct the list of test points
  // const std::vector<double> testT = { -2.0, -0.01, 0.0,  0.01, 0.5,  0.99, 1.0, 1.01, 1.99, 2.0,  2.01, 2.49, 2.5,
  //                                    2.51, 3.0,   3.99, 4.0,  4.01, 4.49, 4.5, 4.51, 5.0,  5.99, 6.0,  6.01, 9.0 };
  const std::vector<double> testT = { -1.99, -0.01, 0.01, 0.5,  0.99, 1.01, 1.99, 2.01, 2.49, 2.5, 2.51,
                                      3.0,   3.99,  4.01, 4.49, 4.5,  4.51, 5.0,  5.99, 6.01, 9.0 };

  // do the test
  testVectorFunctionFirstDerivative(testCurve, Eigen::Map<const Eigen::VectorXd>(testT.data(), testT.size()));
}

TEST(BestFitDerivatives, ConstructedToleranceZoneWeightDerivatives_noOriginalWeight)
{
  // create inner and outer tolerance curves
  Eigen::Matrix3Xd innerTXI(3, 4);
  innerTXI << 0.0, 2.5, 4.5, 6.0, // t-values
      0.1, 0.2, 0.3, 0.1,         // tolerance values
      0.0, 0.1, 0.0, 0.0;         // tolerance derivative values
  auto innerToleranceCurve = createHermiteCurve1D(innerTXI, true);
  Eigen::Matrix3Xd outerTXI(3, 2);
  outerTXI << 0.0, 6.0, // t-values
      0.4, 0.4,         // tolerance values
      0.0, 0.0;         // tolerance derivative values
  auto outerToleranceCurve = createHermiteCurve1D(outerTXI, true);

  // construct the middle-of-tolerance-zone test curve
  ConstructedToleranceZoneWeight testCurve(*innerToleranceCurve, *outerToleranceCurve);

  // construct the list of test points
  // const std::vector<double> testT = { -2.0, -0.01, 0.0,  0.01, 0.5,  2.49, 2.5,  2.51, 3.0,
  //                                    4.49, 4.5,   4.51, 5.0,  5.99, 6.0,  6.01, 9.0 };
  const std::vector<double> testT = { -2.0, -0.01, 0.01, 0.5, 2.49, 2.51, 3.0, 4.49, 4.51, 5.0, 5.99, 6.01, 9.0 };

  // do the test
  testVectorFunctionDerivatives(testCurve, Eigen::Map<const Eigen::VectorXd>(testT.data(), testT.size()));
}

TEST(BestFitDerivatives, ConstructedToleranceZoneWeightDerivatives_withOriginalWeight)
{
  // create inner and outer tolerance curves
  Eigen::Matrix3Xd innerTXI(3, 4);
  innerTXI << 0.0, 2.5, 4.5, 6.0, // t-values
      0.1, 0.2, 0.3, 0.1,         // tolerance values
      0.0, 0.1, 0.0, 0.0;         // tolerance derivative values
  auto innerToleranceCurve = createHermiteCurve1D(innerTXI, true);
  Eigen::Matrix3Xd outerTXI(3, 2);
  outerTXI << 0.0, 6.0, // t-values
      0.4, 0.4,         // tolerance values
      0.0, 0.0;         // tolerance derivative values
  auto outerToleranceCurve = createHermiteCurve1D(outerTXI, true);

  // create an original weight
  Eigen::Matrix3Xd weightTXI(3, 4);
  weightTXI << 0.0, 1.7, 4.2, 6.0, // t-values
      0.1, 0.2, 0.2, 0.1,          // tolerance values
      0.0, 0.0, 0.0, 0.0;          // tolerance derivative values
  auto originalWeightCurve = createHermiteCurve1D(weightTXI, true);

  // construct the middle-of-tolerance-zone test curve
  ConstructedToleranceZoneWeight testCurve(*originalWeightCurve, *innerToleranceCurve, *outerToleranceCurve);

  // construct the list of test points
  // const std::vector<double> testT = { -2.0, -0.01, 0.0,  0.01, 0.5, 1.69, 1.70, 1.71, 2.49, 2.5,  2.51, 3.0,
  //                                    4.19, 4.20,  4.21, 4.49, 4.5, 4.51, 5.0,  5.99, 6.0,  6.01, 9.0 };
  const std::vector<double> testT = { -2.0, -0.01, 0.01, 0.5,  1.69, 1.71, 2.49, 2.51, 3.0,
                                      4.19, 4.21,  4.49, 4.51, 5.0,  5.99, 6.01, 9.0 };

  // do the test
  testVectorFunctionDerivatives(testCurve, Eigen::Map<const Eigen::VectorXd>(testT.data(), testT.size()));
}
} // namespace Blade
} // namespace Hexagon

#endif
