#include "stdafx.h"
#include "MeanCamberCurve.h"
#include "HermiteCurve.h"
#include "SUBCURVE.H"
#include "CurvePolygon.h"
#include "SECTION.H"
#include "NOMINALSECTION.H"
#include "ToleranceSection.h"
#include "Analysis.h"
#include "CIRCLE.H"
#include "TemplateHermiteSpline.h"


#ifdef min
#undef min
#endif

#pragma warning(push)
#pragma warning(disable : 4018) // nanoflann generates warning C4018
#pragma warning(disable : 4458) // nanoflann generates warning C4458
#pragma warning(disable : 4267) // nanoflann generates warning C4267
#pragma warning(disable : 4244) // nanoflann generates warning C4244
#pragma warning(disable : 4127) // nanoflann generates warning C4127
#pragma warning(disable : 4100) // nanoflann generates warning C4100
#include <nanoflann/nanoflann.hpp>
#include "NURB.H"
#pragma warning(pop)

namespace Hexagon
{
namespace Blade
{
// forward-declare a matrix-saving function
void savetxt(const Eigen::Ref<const Eigen::MatrixXd>& input, const std::string& fileName);
}
}


namespace Hexagon
{
namespace Blade
{

void writeMeanCamberCurveParameters2016(const MeanCamberCurveParameters2016* params, FILE* fp)
{
  // write a zero character if there is no parameters object
  if(!params)
  {
    fwrite("\0", sizeof(char), 1, fp);
    return;
  }

  // write a 'p' character otherwise (for "params")
  fwrite("p", sizeof(char), 1, fp);

  // now write the nose data
  if(!params->noseBackoff)
  {
    // write a zero character for no nose backoff
    fwrite("\0", sizeof(char), 1, fp);
  }
  else
  {
    // write an 'n' character for having a "nose" backoff
    fwrite("n", sizeof(char), 1, fp);
    assert(params->noseBackoff->point);
    assert(params->noseBackoff->normal);
    fwrite(params->noseBackoff->point->data(), sizeof(double), 2, fp);
    fwrite(params->noseBackoff->normal->data(), sizeof(double), 2, fp);
    fwrite(&params->noseBackoff->backoffDistance, sizeof(double), 1, fp);
  }
  fwrite(&params->noseStart, sizeof(double), 1, fp);
  fwrite(&params->noseEnd, sizeof(double), 1, fp);

  // and the tail data
  if(!params->tailBackoff)
  {
    // write a zero character for no tail backoff
    fwrite("\0", sizeof(char), 1, fp);
  }
  else
  {
    // write an 't' character for having a "tail" backoff
    fwrite("t", sizeof(char), 1, fp);
    assert(params->tailBackoff->point);
    assert(params->tailBackoff->normal);
    fwrite(params->tailBackoff->point->data(), sizeof(double), 2, fp);
    fwrite(params->tailBackoff->normal->data(), sizeof(double), 2, fp);
    fwrite(&params->tailBackoff->backoffDistance, sizeof(double), 1, fp);
  }
  fwrite(&params->tailStart, sizeof(double), 1, fp);
  fwrite(&params->tailEnd, sizeof(double), 1, fp);

  // all done
}

std::shared_ptr<CamberBackoff> readCamberBackoff(FILE* fp)
{
  auto result = std::make_shared<CamberBackoff>();
  result->point.reset(new Eigen::Array2d());
  result->normal.reset(new Eigen::Array2d());
  fread(result->point->data(), sizeof(double), 2, fp);
  fread(result->normal->data(), sizeof(double), 2, fp);
  fread(&result->backoffDistance, sizeof(double), 1, fp);
  return result;
}

std::shared_ptr<const MeanCamberCurveParameters2016> readMeanCamberCurveParameters2016(FILE* fp)
{
  // read one character (which determines whether there is an object or not)
  char objectIndicator;
  fread(&objectIndicator, sizeof(char), 1, fp);

  // is there an object?
  if(objectIndicator != 'p')
  {
    return nullptr;
  }

  // there is an object; allocate it
  auto result = std::make_shared<MeanCamberCurveParameters2016>();

  // now look to see if there is a nose backoff
  fread(&objectIndicator, sizeof(char), 1, fp);
  if(objectIndicator == 'n')
  {
    result->noseBackoff = readCamberBackoff(fp);
  }
  fread(&result->noseStart, sizeof(double), 1, fp);
  fread(&result->noseEnd, sizeof(double), 1, fp);

  // now look to see if there is a tail backoff
  fread(&objectIndicator, sizeof(char), 1, fp);
  if(objectIndicator == 't')
  {
    result->tailBackoff = readCamberBackoff(fp);
  }
  fread(&result->tailStart, sizeof(double), 1, fp);
  fread(&result->tailEnd, sizeof(double), 1, fp);

  // all done
  return result;
}

typedef HermiteSpline<2, PeriodicSplineType::Nonperiodic> SplineType;

struct OpenHermiteSpline : public SplineType
{
};

std::vector<ptrdiff_t> findSignTransitions(const Eigen::Ref<const Eigen::ArrayXd>& x)
{
  auto boundaryCrossed = [](double a, double b) { return (a < 0.0) != (b < 0.0); };

  std::vector<ptrdiff_t> transitions;
  for(ptrdiff_t i = 0; i < x.size(); i++)
  {
    ptrdiff_t iPlus1 = (i + 1) % x.size();
    if(boundaryCrossed(x[i], x[iPlus1]))
    {
      transitions.push_back(i);
    }
  }

  return transitions;
}

double findTransition(CCurve* wholeCurve, const Eigen::Ref<const Eigen::Array2d>& halfSpacePoint,
                      const Eigen::Ref<const Eigen::Array2d>& halfSpaceNormal, double tLow, double tHigh)
{
  // define a function that has a zero at the transition
  auto signedDistance = [&](double t) -> double 
  {
    Eigen::Array2d point;
    wholeCurve->CalcPoint(point.data(), t);
    return ((point - halfSpacePoint) * halfSpaceNormal).sum();
  };

  // make sure the low and high values make sense
  if(tLow >= tHigh)
  {
    tHigh += wholeCurve->Period();
  }

  // find the zero of the function between the low and high values
  return findZero(signedDistance, tLow, tHigh);
}

Eigen::ArrayXd findTransitions(CCurve* wholeCurve, const Eigen::Ref<const Eigen::Array2d>& halfSpacePoint,
                               const Eigen::Ref<const Eigen::Array2d>& halfSpaceNormal)
{
  // evaluate the total curve finely
  const ptrdiff_t samples = 10000;
  Eigen::ArrayXd t = Eigen::ArrayXd::LinSpaced(samples, wholeCurve->T0(), wholeCurve->T1());
  Eigen::Array2Xd points(2, samples);
  wholeCurve->CalcPoints(points.data(), t.data(), samples);

#ifndef NDEBUG
  // save some data
  savetxt(points.transpose(), "D:/temp/sampledCurve.txt");
#endif

  // evaluate the signed distances to the half-space boundaries
  Eigen::ArrayXd boundaryDistances =
      ((points.colwise() - halfSpacePoint).colwise() * halfSpaceNormal).colwise().sum().transpose();

  // figure out the transitions
  auto coarseTransitions = findSignTransitions(boundaryDistances);
  Eigen::ArrayXd transitions(coarseTransitions.size());
  for(ptrdiff_t i = 0; i < transitions.size(); i++)
  {
    transitions[i] = findTransition(wholeCurve, halfSpacePoint, halfSpaceNormal, t[coarseTransitions[i]],
                                    t[(coarseTransitions[i] + 1) % samples]);
  }

  // all done
  return transitions;
}

Eigen::Array2d findTwoTransitions(CCurve* wholeCurve, const Eigen::Ref<const Eigen::Array2d>& curvePoint,
                                  const Eigen::Ref<const Eigen::Array2d>& curveNormal, double& inout_backoffDistance)
{
  Eigen::ArrayXd transitions =
      findTransitions(wholeCurve, curvePoint - inout_backoffDistance * curveNormal, curveNormal);

  while(transitions.size() != 2)
  {
    if(transitions.size() < 2)
    {
      throw std::logic_error("I lost the curve somehow");
    }
    inout_backoffDistance *= 1.3;
    transitions = findTransitions(wholeCurve, curvePoint - inout_backoffDistance * curveNormal, curveNormal);
  }

  return transitions.head<2>();
}

std::tuple<double, double, double, double, double, double>
splitCurveBasedOnHalfSpaces(const MeanCamberCurveParameters2016& params)
{
  // figure out the transitions
  double noseBackoffDistance = std::numeric_limits<double>::quiet_NaN();
  Eigen::Array2d noseTransitions;
  if(params.noseBackoff)
  {
    noseBackoffDistance = params.noseBackoff->backoffDistance;
    noseTransitions = findTwoTransitions(params.wholeCurve, *params.noseBackoff->point, *params.noseBackoff->normal,
                                         noseBackoffDistance);
  }
  else
  {
    noseTransitions[0] = std::min(params.noseStart, params.noseEnd);
    noseTransitions[1] = std::max(params.noseStart, params.noseEnd);
  }
  double tailBackoffDistance = std::numeric_limits<double>::quiet_NaN();
  Eigen::Array2d tailTransitions;
  if(params.tailBackoff)
  {
    tailBackoffDistance = params.tailBackoff->backoffDistance;
    tailTransitions = findTwoTransitions(params.wholeCurve, *params.tailBackoff->point, *params.tailBackoff->normal,
                                         tailBackoffDistance);
  }
  else
  {
    tailTransitions[0] = std::min(params.tailStart, params.tailEnd);
    tailTransitions[1] = std::max(params.tailStart, params.tailEnd);
  }

  // figure out the segments
  const double period = params.wholeCurve->Period();
  auto nextPast = [&](double x, double pastThis) { return x - std::floor((x - pastThis) / period) * period; };
  auto constructSegments = [&](double segment1rawStart, double segment1rawEnd, double segment2rawStart,
                               double segment2rawEnd) 
  {
    double segment1start = segment1rawStart;
    double segment1end = nextPast(segment1rawEnd, segment1start);
    double segment2start = nextPast(segment2rawStart, segment1end);
    double segment2end = nextPast(segment2rawEnd, segment2start);
    return std::make_tuple(segment1start, segment1end, segment2start, segment2end);
  };
  double segment1start, segment1end, segment2start, segment2end;
  std::tie(segment1start, segment1end, segment2start, segment2end) =
      constructSegments(noseTransitions[0], tailTransitions[0], tailTransitions[1], noseTransitions[1]);
  auto doSegmentsMakeSense = [&]() { return segment2end - period < segment1start; };
  bool segmentsMakeSense = doSegmentsMakeSense();
  if(!segmentsMakeSense)
  {
    std::tie(segment1start, segment1end, segment2start, segment2end) =
        constructSegments(noseTransitions[1], tailTransitions[0], tailTransitions[1], noseTransitions[0]);
    segmentsMakeSense = doSegmentsMakeSense();
  }
  //alwaysAssert(segmentsMakeSense);

  // all done
  return std::make_tuple(segment1start, segment1end, segment2start, segment2end, noseBackoffDistance,
                         tailBackoffDistance);
}



std::unique_ptr<HermiteOpenCurve> computeMeanCamber(const MeanCamberCurveParameters2016& params,
                                                    const Eigen::Ref<const Eigen::Array2d>& limitsA,
                                                    const Eigen::Ref<const Eigen::Array2d>& limitsB, bool isEnglish)
{
  // wrap the inputs so I can call the Voronoi-camber routine
  CSubCurve curveA(params.wholeCurve, limitsA[0], limitsA[1], params.wholeCurve->Period());
  CSubCurve curveB(params.wholeCurve, limitsB[0], limitsB[1], params.wholeCurve->Period());

  // figure out the half-space defined by the limits
  Eigen::Vector2d pointA0, pointA1, pointB0, pointB1;
  params.wholeCurve->CalcPoint(pointA0.data(), limitsA[0]);
  params.wholeCurve->CalcPoint(pointA1.data(), limitsA[1]);
  params.wholeCurve->CalcPoint(pointB0.data(), limitsB[0]);
  params.wholeCurve->CalcPoint(pointB1.data(), limitsB[1]);
  Eigen::Vector2d point1, vector1, point2, vector2;
  point1 = 0.5 * (pointA0 + pointB1);
  point2 = 0.5 * (pointA1 + pointB0);
  vector1[0] = pointA0[1] - pointB1[1];
  vector1[1] = -pointA0[0] + pointB1[0];
  vector1.normalize();
  vector2[0] = pointA1[1] - pointB0[1];
  vector2[1] = -pointA1[0] + pointB0[0];
  vector2.normalize();
  if((point1-point2).dot(vector1) < 0.0)
  {
    vector1 *= -1.0;
  }
  if((point2-point1).dot(vector2) < 0.0)
  {
    vector2 *= -1.0;
  }

  // check correctness
  if(params.noseBackoff)
  {
    assert((params.noseBackoff->normal->matrix() - vector1).norm() < 1e-8);
    Eigen::Vector2d point =
        *params.noseBackoff->point - *params.noseBackoff->normal * params.noseBackoff->backoffDistance;
    assert(std::abs((point - point1).dot(vector1)) < 1e-8);
  }
  if(params.tailBackoff)
  {
    assert((params.tailBackoff->normal->matrix() - vector2).norm() < 1e-8);
    Eigen::Vector2d point =
        *params.tailBackoff->point - *params.tailBackoff->normal * params.tailBackoff->backoffDistance;
    assert(std::abs((point - point2).dot(vector2)) < 1e-8);
  }

  // call the Voronoi-camber routine
  Camber camberResult;
  try
  {
    camberResult = camber(*params.wholeCurve, curveA, curveB, 2048, 1e-5, point1, vector1, point2, vector2);
  }
  catch(std::exception& e)
  {
    throw std::logic_error("Camber line calculation failed, with error '" + std::string(e.what()) + "'.");
  }
  catch(...)
  {
    throw std::logic_error("Camber line calculation failed, with unknown error.");
  }

  // convert the camber result into an OpenSpline
  const ptrdiff_t N = camberResult.camberPoints.cols();
  Eigen::Matrix2Xd camberNormals(2, N);
  camberNormals.row(0) = camberResult.camberTangents.row(1);
  camberNormals.row(1) = -camberResult.camberTangents.row(0);
  auto meanCamberSpline = initialSpline<PeriodicSplineType::Nonperiodic>(camberResult.camberPoints, camberNormals);
  auto result = std::make_unique<HermiteOpenCurve>(meanCamberSpline.points, meanCamberSpline.tangents,
                                                   meanCamberSpline.t, isEnglish);
  result->addDoubleAnnotation("LeadingCamberRadius", camberResult.camberRadii[0]);
  result->addDoubleAnnotation("TrailingCamberRadius", camberResult.camberRadii[N - 1]);
  return result;
}

std::unique_ptr<HermiteOpenCurve> extrapolateMeanCamber(CCurve* wholeCurve, HermiteOpenCurve* meanCamber,
                                                        bool extrapolateLeading, bool extrapolateTrailing)
{
  // access core polygonalization
  Eigen::Matrix2Xd curvePolygon;
  Eigen::VectorXd curveT;
  std::tie(curvePolygon, curveT) = polygonalizeWithT(*wholeCurve, 2048, 1e-5);

  // make a KD tree for the core polygon
  typedef nanoflann::KDTreeEigenMatrixAdaptor<Eigen::MatrixX2d, 2> KDTree;
  Eigen::MatrixX2d curvePolygonTranspose = curvePolygon.transpose();
  KDTree kdtree(2 /*dimensions*/, curvePolygonTranspose, 10 /*max leaf*/);
  auto closestIndex = [&kdtree](const Eigen::Ref<const Eigen::Array2d>& inputPoint) -> ptrdiff_t
  {
    ptrdiff_t index;
    double distanceSquared;
    kdtree.query(inputPoint.data(), 1, &index, &distanceSquared);
    return static_cast<ptrdiff_t>(index);
  };
  auto closestT = [&](const Eigen::Ref<const Eigen::Array2d>& inputPoint) -> double
  {
    return curveT[closestIndex(inputPoint)];
  };
  auto closestPoint = [&](const Eigen::Ref<const Eigen::Array2d>& inputPoint) -> Eigen::Vector2d
  {
    return curvePolygon.col(closestIndex(inputPoint));
  };

  // gather information about the ends
  double t0 = meanCamber->T0();
  double t1 = meanCamber->T1();
  Eigen::Vector2d point0;
  Eigen::Vector2d tangent0;
  Eigen::Vector2d point1;
  Eigen::Vector2d tangent1;
  meanCamber->CalcPoint(point0.data(), t0, tangent0.data());
  meanCamber->CalcPoint(point1.data(), t1, tangent1.data());
  tangent0.normalize();
  tangent0 *= -1.0;

  // construct the extrapolation for the 0th end, adding an initial guess and a large number of seeds
  Eigen::Vector2d intersectionPoint0;
  double intersectionPointTGuess0 = closestT(point0);
  auto intersectionSuccess0 = wholeCurve->NewLineIntersect(point0.data(), tangent0.data(), intersectionPoint0.data(),
                                                           0.0, 0.0, &intersectionPointTGuess0, nullptr, true);
  double extrapolationT0 = -(intersectionPoint0 - point0).norm();
  if(!extrapolateLeading)
  {
    intersectionSuccess0 = 0;
  }

  // construct the extrapolation for the 1th end, adding an initial guess and a large number of seeds
  Eigen::Vector2d intersectionPoint1;
  double intersectionPointTGuess1 = closestT(point1);
  auto intersectionSuccess1 = wholeCurve->NewLineIntersect(point1.data(), tangent1.data(), intersectionPoint1.data(),
                                                           0.0, 0.0, &intersectionPointTGuess1, nullptr, true);
  double extrapolationT1 = (intersectionPoint1 - point1).norm();
  if(!extrapolateTrailing)
  {
    intersectionSuccess1 = 0;
  }

  // construct the final spline
  const auto& formerSpline = dynamic_cast<const SplineType&>(meanCamber->getSpline());
  ptrdiff_t N = formerSpline.t.size();
  ptrdiff_t N_new = N + intersectionSuccess0 + intersectionSuccess1;
  Eigen::VectorXd t(N_new);
  Eigen::Matrix2Xd points(2, N_new);
  Eigen::Matrix2Xd tangents(2, N_new);
  if(intersectionSuccess0)
  {
    t[0] = 0.0;
    t.segment(1, N) = -extrapolationT0 + formerSpline.t.array();
    points.col(0) = intersectionPoint0;
    points.block(0, 1, 2, N) = formerSpline.points;
    tangents.col(0) = -tangent0;
    tangents.block(0, 1, 2, N) = formerSpline.tangents;
  }
  else
  {
    t.head(N) = formerSpline.t.array();
    points.leftCols(N) = formerSpline.points;
    tangents.leftCols(N) = formerSpline.tangents;
  }
  if(intersectionSuccess1)
  {
    t[N_new - 1] = t[N_new - 2] + extrapolationT1;
    points.col(N_new - 1) = intersectionPoint1;
    tangents.col(N_new - 1) = tangent1;
  }

  // create the extrapolated spline
  auto result = std::make_unique<HermiteOpenCurve>(points, tangents, t, meanCamber->IsEnglish());

  // add some annotations
  result->addDoubleAnnotation("leadingInscribedX", formerSpline.points(0, 0));
  result->addDoubleAnnotation("leadingInscribedY", formerSpline.points(1, 0));
  result->addDoubleAnnotation("leadingInscribedR", meanCamber->getDoubleAnnotation("LeadingCamberRadius"));
  result->addDoubleAnnotation("trailingInscribedX", formerSpline.points(0, N - 1));
  result->addDoubleAnnotation("trailingInscribedY", formerSpline.points(1, N - 1));
  result->addDoubleAnnotation("trailingInscribedR", meanCamber->getDoubleAnnotation("TrailingCamberRadius"));

  // all done
  return result;
}

CCurve* createMeanCamberCurve2016(const MeanCamberCurveParameters2016& initialParams, bool isEnglish)
{
  Eigen::Array2d limitsA;
  Eigen::Array2d limitsB;
  double finalNoseBackoffDistance;
  double finalTailBackoffDistance;
  std::tie(limitsA[0], limitsA[1], limitsB[0], limitsB[1], finalNoseBackoffDistance, finalTailBackoffDistance) =
      splitCurveBasedOnHalfSpaces(initialParams);
  MeanCamberCurveParameters2016 finalParams = initialParams;
  if(initialParams.noseBackoff)
  {
    finalParams.noseBackoff = std::make_shared<CamberBackoff>(*initialParams.noseBackoff);
    finalParams.noseBackoff->backoffDistance = finalNoseBackoffDistance;
  }
  if(initialParams.tailBackoff)
  {
    finalParams.tailBackoff = std::make_shared<CamberBackoff>(*initialParams.tailBackoff);
    finalParams.tailBackoff->backoffDistance = finalTailBackoffDistance;
  }
  auto nonExtrapolatedMeanCamber = computeMeanCamber(finalParams, limitsA, limitsB, isEnglish);
  auto extrapolatedMeanCamber =
      extrapolateMeanCamber(finalParams.wholeCurve, nonExtrapolatedMeanCamber.get(),
                            std::isfinite(finalNoseBackoffDistance), std::isfinite(finalTailBackoffDistance));

  // add information that might be reported later
  extrapolatedMeanCamber->addBoolAnnotation("hasHalfSpaces", true);
  extrapolatedMeanCamber->addDoubleAnnotation("leadingConvexHalfSpaceIntersection", limitsA[0]);
  extrapolatedMeanCamber->addDoubleAnnotation("leadingConcaveHalfSpaceIntersection", limitsB[1]);
  extrapolatedMeanCamber->addDoubleAnnotation("trailingConvexHalfSpaceIntersection", limitsA[1]);
  extrapolatedMeanCamber->addDoubleAnnotation("trailingConcaveHalfSpaceIntersection", limitsB[0]);

  // all done
  return extrapolatedMeanCamber.release();
}

CCurve* createMeasuredMeanCamberCurve2016(CCurve* measuredWholeCurve, CCurve* measuredLeadingCurve,
                                          CCurve* measuredTrailingCurve,
                                          const MeanCamberCurveParameters2016* nominalMCLParams, bool isEnglish)
{
  assert(nominalMCLParams);

  MeanCamberCurveParameters2016 params = *nominalMCLParams;
  params.wholeCurve = measuredWholeCurve;
  params.noseStart = measuredLeadingCurve->T0();
  params.noseEnd = measuredLeadingCurve->T1();
  params.tailStart = measuredTrailingCurve->T0();
  params.tailEnd = measuredTrailingCurve->T1();

  // fix the nose point for measured curve
  Eigen::Array2Xd measuredPolygon = polygonalize(*params.wholeCurve, 2048, 1e-5);
  if(nominalMCLParams->noseBackoff)
  {
    Eigen::ArrayXd noseDistances = ((measuredPolygon.colwise() - *nominalMCLParams->noseBackoff->point).colwise() *
                                    (*nominalMCLParams->noseBackoff->normal))
                                       .colwise()
                                       .sum();
    ptrdiff_t measuredNoseIndex;
    noseDistances.maxCoeff(&measuredNoseIndex);

    params.noseBackoff->point.reset(new Eigen::Array2d(measuredPolygon.col(measuredNoseIndex)));
  }

  // fix the tail point for measured curve
  if(nominalMCLParams->tailBackoff)
  {
    Eigen::ArrayXd tailDistances = ((measuredPolygon.colwise() - *nominalMCLParams->tailBackoff->point).colwise() *
                                    (*nominalMCLParams->tailBackoff->normal))
                                       .colwise()
                                       .sum();
    ptrdiff_t measuredTailIndex;
    tailDistances.maxCoeff(&measuredTailIndex);

    params.tailBackoff->point.reset(new Eigen::Array2d(measuredPolygon.col(measuredTailIndex)));
  }

  // create the mean camber curve
  return createMeanCamberCurve2016(params, isEnglish);
}

MeanCamberResult createMeanCamberCurveOldMeasured(const MeanCamberCurveParameters& params)
{
  double ler = params.ler;
  double ter = params.ter;
  double voff = params.voff;
  double uoff = params.uoff;
  double mtle = params.mtle;
  double mtte = params.mtte;
  auto section = params.section;
  CCurve* whole = NULL;
  //CCurve *whole = section->MeaCurve();
  //CCurve *lec = section->MeaPart(LEC);
  //CCurve *tec = section->MeaPart(TEC);
  //CCurve *cvc = section->MeaPart(CVC);
  //CCurve *ccc = section->MeaPart(CCC);
  CCurve *mcc = nullptr;
  bool skipPitch = false;

  //// build mean camber
  //int letype = section->LEType();
  //int tetype = section->TEType();
  int letype = 0;
  int tetype = 0;


  if (letype != EDGE_NORMAL || tetype != EDGE_NORMAL)
  {
  //  double xyz[3], rim[3];
  //  CCurve *uselec = section->MeaPart(LEC);
  //  if (letype == EDGE_SQUARE)
  //  {
  //    xyz[0] = params.analysisSection->m_nose[0];  //changed how square ends are handled, trying to use the old way
  //    xyz[1] = params.analysisSection->m_nose[1];
  //    xyz[2] = 0.0;
  //    rim[0] = xyz[0] + 1.0;
  //    rim[1] = xyz[1];
  //    uselec = new CCircle(xyz, rim, rim, 1.0);
  //    //bugout(0, _T("10 %f %f LE"), xyz[0], xyz[1]);
  //  }
  //  else if (letype == EDGE_PARTIAL)
  //    uselec = 0;

  //  CCurve *usetec = tec;
  //  if (tetype == EDGE_SQUARE)
  //  {
  //    xyz[0] = params.analysisSection->m_tail[0];  //changed how square ends are handled, trying to use the old way
  //    xyz[1] = params.analysisSection->m_tail[1];
  //    xyz[2] = 0.0;
  //    rim[0] = xyz[0] + 1.0;
  //    rim[1] = xyz[1];
  //    usetec = new CCircle(xyz, rim, rim, 1.0);
  //    //bugout(0, _T("10 %f %f TE"), xyz[0], xyz[1]);
  //  }
  //  else if (tetype == EDGE_PARTIAL)
  //    usetec = 0;

  //  int numMeaMCLPoints = myGetProfileInt(L"MeasuredMCLPoints", 75);

  //  // If both ends are partial, CNurbCurve will need help deciding which end is LE.  The start of the nominal MCL is passed in to be used when needed
  //  double nomStart[2];
  //  params.section->NomPart(MCC)->CalcPoint(nomStart, section->NomPart(MCC)->T0());

   // mcc = new CNurbCurve(numMeaMCLPoints, ccc, cvc, uselec, usetec, 3, nomStart);

  /*  if (usetec && usetec != tec)
      delete usetec;

    if (uselec && uselec != lec)
      delete uselec;*/
  }
  else
  {
    if (ler < 0.0)
      ler = 0.75*voff;
    if (ter < 0.0)
      ter = 0.75*uoff;

    if (mcc)
      delete mcc;

    // ADJUSTENDS
    int numMeaMCLPoints = 75;//myGetProfileInt(L"MeasuredMCLPoints", 75);
    mcc=  new CNurbCurve(numMeaMCLPoints, mtle, mtte, ler, ter, whole, 0);
  }

  MeanCamberResult result;
  result.meanCamberCurve = mcc;
  result.skipPitch = skipPitch;
  return result;
}

MeanCamberResult createMeasuredMeanCamberCurve(const MeanCamberCurveParameters& params, bool isEnglish)
{
    MeanCamberResult result;
    result.skipPitch = false;
    return result;
  // check sanity; the sections are supposed to all be the same section
  //std::wstring sectionName = params.section->m_name;
  //std::wstring analysisName = params.analysisSection->m_sectName;
  //alwaysAssert(sectionName == analysisName);

  //if(params.section->nominalMCLParams)
  //{
  //  MeanCamberResult result;
  //  result.meanCamberCurve =
  //      createMeasuredMeanCamberCurve2016(params.section->MeaCurve(), params.section->MeaPart(LEC),
  //                                        params.section->MeaPart(TEC), params.section->nominalMCLParams, isEnglish);
  //  result.skipPitch = false;
  //  return result;
  //}
  //else
  //{
  //  return createMeanCamberCurveOldMeasured(params);
  //}
}
}
}