#include "stdafx.h"
#include "HermiteCurve.h"
#include "CurvePolygon.h"
#include "TemplateHermiteSpline.h"
#include <map>
#include <string>
#include <sstream>
#include <iomanip>


namespace Hexagon
{
namespace Blade
{

typedef HermiteSpline<2, PeriodicSplineType::Nonperiodic> SplineType;

const wchar_t* HermiteOpenCurve::fileDescriptor = L"hermiteOpenCurve";


struct OpenHermiteSpline : public SplineType
{
};

struct AnnotationMap
{
  std::map<std::string, double> doubleMap;
  std::map<std::string, bool> boolMap;

  AnnotationMap()
  {
  }
  ~AnnotationMap()
  {
  }
};

ptrdiff_t storeSize(const Eigen::Ref<const Eigen::MatrixXd>& matrix)
{
  return 2 * sizeof(int64_t) + matrix.size() * sizeof(double);
}

void writeMatrix(FILE* fp, const Eigen::Ref<const Eigen::MatrixXd>& matrix)
{
  int64_t numberOfRows = static_cast<int64_t>(matrix.rows());
  int64_t numberOfColumns = static_cast<int64_t>(matrix.cols());
  fwrite(&numberOfRows, sizeof(int64_t), 1, fp);
  fwrite(&numberOfColumns, sizeof(int64_t), 1, fp);
  fwrite(matrix.data(), sizeof(double), matrix.size(), fp);
}

Eigen::MatrixXd readMatrix(FILE* fp)
{
  int64_t numberOfRows = 0;
  fread(&numberOfRows, sizeof(int64_t), 1, fp);
  int64_t numberOfColumns = 0;
  fread(&numberOfColumns, sizeof(int64_t), 1, fp);
  Eigen::MatrixXd result(static_cast<ptrdiff_t>(numberOfRows), static_cast<ptrdiff_t>(numberOfColumns));
  fread(result.data(), sizeof(double), result.size(), fp);
  return result;
}


std::string serializeAnnotationMap(const AnnotationMap& map)
{
  std::ostringstream serializer(std::ostringstream::out | std::ostringstream::binary);

  // write the boolean map
  int64_t boolMapSize = map.boolMap.size();
  serializer.write(reinterpret_cast<const char*>(&boolMapSize), sizeof(int64_t));
  for(const auto& element : map.boolMap)
  {
    int64_t keySize = element.first.size();
    serializer.write(reinterpret_cast<const char*>(&keySize), sizeof(int64_t));
    serializer.write(element.first.data(), element.first.size());
    int64_t value = element.second;
    serializer.write(reinterpret_cast<const char*>(&value), sizeof(int64_t));
  }

  // write the double map
  int64_t doubleMapSize = map.doubleMap.size();
  serializer.write(reinterpret_cast<const char*>(&doubleMapSize), sizeof(int64_t));
  for(const auto& element : map.doubleMap)
  {
    int64_t keySize = element.first.size();
    serializer.write(reinterpret_cast<const char*>(&keySize), sizeof(int64_t));
    serializer.write(element.first.data(), element.first.size());
    double value = element.second;
    serializer.write(reinterpret_cast<const char*>(&value), sizeof(double));
  }

  // all done
  return serializer.str();
}

std::shared_ptr<AnnotationMap> deserializeAnnotationMap(const std::string& text)
{
  std::istringstream deserializer(text, std::istringstream::out | std::istringstream::binary);
  auto result = std::make_shared<AnnotationMap>();

  // read the boolean map
  int64_t boolMapSize;
  deserializer.read(reinterpret_cast<char*>(&boolMapSize), sizeof(int64_t));
  for(int64_t i = 0; i < boolMapSize; i++)
  {
    int64_t keySize;
    deserializer.read(reinterpret_cast<char*>(&keySize), sizeof(int64_t));
    std::string key(static_cast<size_t>(keySize), ' ');
    deserializer.read(&key[0], static_cast<size_t>(keySize));
    int64_t value;
    deserializer.read(reinterpret_cast<char*>(&value), sizeof(int64_t));
#pragma warning(suppress : 4800)
    result->boolMap[key] = static_cast<bool>(value);
  }

  // read the double map
  int64_t doubleMapSize;
  deserializer.read(reinterpret_cast<char*>(&doubleMapSize), sizeof(int64_t));
  for(int64_t i = 0; i < doubleMapSize; i++)
  {
    int64_t keySize;
    deserializer.read(reinterpret_cast<char*>(&keySize), sizeof(int64_t));
    std::string key(static_cast<size_t>(keySize), ' ');
    deserializer.read(&key[0], static_cast<size_t>(keySize));
    double value;
    deserializer.read(reinterpret_cast<char*>(&value), sizeof(double));
    result->doubleMap[key] = value;
  }

  // all done
  return result;
}

ptrdiff_t storeSize(const AnnotationMap& map)
{
  return sizeof(int64_t) + serializeAnnotationMap(map).size();
}

void writeAnnotationMap(FILE* fp, const AnnotationMap& map)
{
  auto serializedMap = serializeAnnotationMap(map);
  int64_t mapSize = serializedMap.size();
  fwrite(&mapSize, sizeof(int64_t), 1, fp);
  fwrite(serializedMap.data(), sizeof(char), static_cast<size_t>(mapSize), fp);
}

std::shared_ptr<AnnotationMap> readAnnotationMap(FILE* fp)
{
  int64_t mapSize;
  fread(&mapSize, sizeof(int64_t), 1, fp);
  std::string serializedMap(static_cast<size_t>(mapSize), ' ');
  fread(&serializedMap[0], sizeof(char), static_cast<size_t>(mapSize), fp);
  return deserializeAnnotationMap(serializedMap);
}


// during reporting, it can be useful to retrieve information about the curve;
// I call these annotations
void HermiteOpenCurve::addDoubleAnnotation(const char* key, double value)
{
  annotations->doubleMap[key] = value;
}
double HermiteOpenCurve::getDoubleAnnotation(const char* key) const
{
  auto it = annotations->doubleMap.find(key);
  if(it != annotations->doubleMap.end())
  {
    return it->second;
  }
  return std::numeric_limits<double>::quiet_NaN();
}
void HermiteOpenCurve::addBoolAnnotation(const char* key, bool value)
{
  annotations->boolMap[key] = value;
}
bool HermiteOpenCurve::getBoolAnnotation(const char* key) const
{
  auto it = annotations->boolMap.find(key);
  if(it != annotations->boolMap.end())
  {
    return it->second;
  }
  return false;
}

HermiteOpenCurve::HermiteOpenCurve(FILE* fp, bool isEnglish) : CCurve()
{
  // base-class data-member
  m_isEnglish = isEnglish;

  // check the descriptor
  wchar_t writtenDescriptor[fileDescriptorLength];
  fread(writtenDescriptor, sizeof(wchar_t), fileDescriptorLength, fp);
  //alwaysAssert(wcsncmp(writtenDescriptor, fileDescriptor, fileDescriptorLength) == 0);

  // read in the matrices
  spline = std::make_shared<OpenHermiteSpline>();
  spline->points = readMatrix(fp);
  spline->tangents = readMatrix(fp);
  spline->t = readMatrix(fp);

  // read in the annotations
  annotations = readAnnotationMap(fp);
}

HermiteOpenCurve::HermiteOpenCurve(const Eigen::Ref<const Eigen::Matrix2Xd>& controlPoints,
                                   const Eigen::Ref<const Eigen::Matrix2Xd>& tangents,
                                   const Eigen::Ref<const Eigen::VectorXd>& tValues, bool isEnglish)
  : CCurve()
{
  // base-class data-member
  m_isEnglish = isEnglish;

  spline = std::make_shared<OpenHermiteSpline>();
  spline->points = controlPoints;
  spline->tangents = tangents;
  spline->t = tValues;

  annotations = std::make_shared<AnnotationMap>();
}

HermiteOpenCurve::~HermiteOpenCurve()
{
}

int HermiteOpenCurve::Type() const
{
  return HERMITE_TYPE;
}

int HermiteOpenCurve::CalcPoint(double* xyz, double t, double* tan /*= 0*/, double* curv /*= 0*/) const
{
  return CalcPoints(xyz, &t, 1, tan, curv);
}

int HermiteOpenCurve::CalcPoints(double* outPoints, const double* inT, ptrdiff_t numberOfT,
                                 double* outTangents /*= nullptr*/, double* outCurves /*= nullptr*/) const
{
  Eigen::Map<const Eigen::VectorXd> eigenT(inT, numberOfT);

  // evaluate the point
  if(outPoints)
  {
    Eigen::Map<Eigen::Matrix2Xd> eigenPoints(outPoints, 2, numberOfT);
    evaluateVector<0, SplineType>(*spline, eigenT, eigenPoints);
    if(m_palign)
      m_palign->MeasToBestVector2Xd(outPoints, 1, outPoints, numberOfT);
  }

  // evaluate the tangent (I believe this is really a derivative, not a normalized tangent)
  if(outTangents)
  {
    Eigen::Map<Eigen::Matrix2Xd> eigenTangents(outTangents, 2, numberOfT);
    evaluateVector<1, SplineType>(*spline, eigenT, eigenTangents);
    if(m_palign)
      m_palign->MeasToBestVector2Xd(outTangents, 0, outTangents, numberOfT);
  }

  // evaluate the curvature (I believe this really a second derivative, not a true curvature)
  if(outCurves)
  {
    Eigen::Map<Eigen::Matrix2Xd> eigenCurvatures(outCurves, 2, numberOfT);
    evaluateVector<2, SplineType>(*spline, eigenT, eigenCurvatures);
    if(m_palign)
      m_palign->MeasToBestVector2Xd(outCurves, 0, outCurves, numberOfT);
  }

  // all done
  return 1;
}


int HermiteOpenCurve::Valid()
{
  return static_cast<bool>(spline);
}

long HermiteOpenCurve::StoreSize()
{
  if(spline)
  {
    return static_cast<long>(fileDescriptorLength * sizeof(wchar_t) + storeSize(spline->points) +
                             storeSize(spline->tangents) + storeSize(spline->t) + storeSize(*annotations));
  }
  else
  {
    throw std::logic_error("It does not make sense to calculate the store size of a HermiteCurve without a spline.");
  }
}

void HermiteOpenCurve::Write(FILE* fp)
{
  if(spline)
  {
    fwrite(fileDescriptor, sizeof(wchar_t), fileDescriptorLength, fp);
    writeMatrix(fp, spline->points);
    writeMatrix(fp, spline->tangents);
    writeMatrix(fp, spline->t);
    writeAnnotationMap(fp, *annotations);
  }
  else
  {
    throw std::logic_error("It does not make sense to write to a file a HermiteCurve without a spline.");
  }
}

void HermiteOpenCurve::Extent(double* min, double* max)
{
  Eigen::Map<Eigen::Vector2d> eigenMin(min);
  Eigen::Map<Eigen::Vector2d> eigenMax(max);

  Eigen::Matrix2Xd polygon = polygonalize(*this, 2048, 1e-4);
  eigenMin = polygon.rowwise().minCoeff();
  eigenMax = polygon.rowwise().maxCoeff();
}

double HermiteOpenCurve::T0() const
{
  return computeLeftDomain<SplineType>(*spline);
}

double HermiteOpenCurve::T1() const
{
  return computeRightDomain<SplineType>(*spline);
}

void HermiteOpenCurve::SetT0(double /*t*/)
{
  throw std::logic_error("It does not make sense to set T0 on a HermiteCurve.");
}

void HermiteOpenCurve::SetT1(double /*t*/)
{
  throw std::logic_error("It does not make sense to set T1 on a HermiteCurve.");
}

double HermiteOpenCurve::compute_length(double /*t0*/, double /*t1*/)
{
  throw std::logic_error("The method or operation is not implemented.");
}

double HermiteOpenCurve::TotalLength(int /*numPoints*/, double /*t0*/, double /*t1*/)
{
  throw std::logic_error("The method or operation is not implemented.");
}

void HermiteOpenCurve::Walk(double& /*t*/, double /*d*/, double* /*pt*/, int /*numPoints*/)
{
  throw std::logic_error("The method or operation is not implemented.");
}

double HermiteOpenCurve::arc_inverse(double /*length*/, double /*t0*/, int& /*index*/)
{
  throw std::logic_error("The method or operation is not implemented.");
}

void HermiteOpenCurve::CenterOfMass(double* /*CM*/, bool /*isGEI*/)
{
  throw std::logic_error("The method or operation is not implemented.");
}

double HermiteOpenCurve::ComputeArea(double /*t0 = 0.0*/, double /*t1 = 0.0*/)
{
  throw std::logic_error("The method or operation is not implemented.");
}

std::shared_ptr<const Curve<1>>
createHermiteCurve1D(const Eigen::Ref<const Eigen::Matrix<double, 3, Eigen::Dynamic>>& txi, bool isPeriodic)
{
  if(isPeriodic)
  {
    return std::make_shared<HermiteCurve<true, 1>>(txi);
  }
  else
  {
    return std::make_shared<HermiteCurve<false, 1>>(txi);
  }
}

std::shared_ptr<const Curve<2>>
createHermiteCurve2D(const Eigen::Ref<const Eigen::Matrix<double, 5, Eigen::Dynamic>>& txyij, bool isPeriodic)
{
  if(isPeriodic)
  {
    return std::make_shared<HermiteCurve<true, 2>>(txyij);
  }
  else
  {
    return std::make_shared<HermiteCurve<false, 2>>(txyij);
  }
}

std::shared_ptr<const Curve<3>>
createHermiteCurve3D(const Eigen::Ref<const Eigen::Matrix<double, 7, Eigen::Dynamic>>& txyzijk, bool isPeriodic)
{
  if(isPeriodic)
  {
    return std::make_shared<HermiteCurve<true, 3>>(txyzijk);
  }
  else
  {
    return std::make_shared<HermiteCurve<false, 3>>(txyzijk);
  }
}
}
}