#pragma once
#include "CURVE.H"
#include <memory>
#include <Eigen/Dense>


namespace Hexagon
{
namespace Blade
{

struct OpenHermiteSpline;
struct AnnotationMap;

class HermiteOpenCurve : public CCurve
{
public:
  static const wchar_t* fileDescriptor;
  static const long fileDescriptorLength = 16;

private:
  std::shared_ptr<OpenHermiteSpline> spline;
  std::shared_ptr<AnnotationMap> annotations;

public:
  HermiteOpenCurve(const Eigen::Ref<const Eigen::Matrix2Xd>& controlPoints,
                   const Eigen::Ref<const Eigen::Matrix2Xd>& tangents, const Eigen::Ref<const Eigen::VectorXd>& tValues,
                   bool isEnglish);
  HermiteOpenCurve(FILE* fp, bool isEnglish);
  virtual ~HermiteOpenCurve();

  virtual int Type() const override;
  virtual int CalcPoint(double* xyz, double t, double* tan = 0, double* curv = 0) const override;
  virtual int CalcPoints(double* outPoints, const double* inT, ptrdiff_t numberOfT, double* outTangents = nullptr,
                         double* outCurves = nullptr) const override;
  virtual int Valid() override;
  virtual long StoreSize() override;
  virtual void Write(FILE* fp) override;
  virtual void Extent(double* min, double* max) override;
  virtual double T0() const override;
  virtual double T1() const override;
  virtual void SetT0(double t) override;
  virtual void SetT1(double t) override;
  virtual double compute_length(double t0, double t1) override;
  virtual double TotalLength(int numPoints, double t0, double t1) override;
  virtual void Walk(double& t, double d, double* pt, int numPoints) override;
  virtual double arc_inverse(double length, double t0, int& index) override;
  virtual void CenterOfMass(double* CM, bool isGEI) override;
  virtual double ComputeArea(double t0 = 0.0, double t1 = 0.0) override;

  const OpenHermiteSpline& getSpline() const
  {
    return *spline;
  }

  // during reporting, it can be useful to retrieve information about the curve;
  // I call these annotations
  virtual void addDoubleAnnotation(const char* key, double value) override final;
  virtual double getDoubleAnnotation(const char* key) const override final;
  virtual void addBoolAnnotation(const char* key, bool value) override final;
  virtual bool getBoolAnnotation(const char* key) const override final;
};

Eigen::Matrix<double, 5, Eigen::Dynamic> createTXYIJ(const Eigen::Ref<const Eigen::Matrix2Xd>& points);
Eigen::Matrix<double, 7, Eigen::Dynamic> createTXYZIJK(const Eigen::Ref<const Eigen::Matrix3Xd>& points);

std::shared_ptr<const Curve<1>> createHermiteCurve1D(const Eigen::Ref<const Eigen::Matrix3Xd>& txi, bool isPeriodic);

std::shared_ptr<const Curve<2>>
createHermiteCurve2D(const Eigen::Ref<const Eigen::Matrix<double, 5, Eigen::Dynamic>>& txyij, bool isPeriodic);

std::shared_ptr<const Curve<3>>
createHermiteCurve3D(const Eigen::Ref<const Eigen::Matrix<double, 7, Eigen::Dynamic>>& txyzijk, bool isPeriodic);
}
}
