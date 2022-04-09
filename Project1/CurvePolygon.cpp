#define _SCL_SECURE_NO_WARNINGS
#pragma warning(disable : 4714)

#include "CurvePolygon.h"
//#include "AlwaysAssert.h"
#include "ArraySlicing.h"
#include "TemplateHermiteSpline.h"
#include "HermiteCurve.h"
#include "EigenAbstractCurve.h"
#include <exception>
#include <array>
#include <mutex>
#include <iostream>
#include <fstream>
#include <map>
#include <set>
#include <iomanip>
#include <string>

extern "C" {
#include <QHull/libqhull.h>
#include <QHull/qset.h>
#include <QHull/mem.h>
#include <QHull/geom.h>
//#include <QHull/ioQhull.h>
};

#pragma warning(push)
#pragma warning(disable : 4018) // nanoflann generates warning C4018
#pragma warning(disable : 4458) // nanoflann generates warning C4458
#pragma warning(disable : 4267) // nanoflann generates warning C4267
#pragma warning(disable : 4244) // nanoflann generates warning C4244
#pragma warning(disable : 4127) // nanoflann generates warning C4127
#pragma warning(disable : 4100) // nanoflann generates warning C4100
#include <nanoflann/nanoflann.hpp>
#pragma warning(pop)

#pragma warning(disable : 4503) // boost generated warning

#pragma warning(push)
#pragma warning(disable : 4100) // boost generated warning
#include <boost/config/warning_disable.hpp>
#include <boost/geometry.hpp>
#include <boost/geometry/geometries/point.hpp>
#include <boost/geometry/geometries/segment.hpp>
#pragma warning(pop)



#ifdef max
#undef max
#endif

#ifdef min
#undef min
#endif


#include <HexagonGDT/ScalarOptimization.h>


namespace QHullInternals
{

// QHull is not thread-safe; we can use this lock to ensure multiple threads do not access it
std::mutex QHullLock;

struct QHullMemoryDeallocator
{
  ~QHullMemoryDeallocator()
  {
    qh_freeqhull(qh_ALL);
    int curlong, totlong;
    qh_memfreeshort(&curlong, &totlong);
  }
};

Eigen::Matrix2Xd computeConvexHull(const Eigen::Ref<const Eigen::Matrix2Xd>& inputPoints)
{
  Eigen::Matrix2Xd inputCopy = inputPoints;
  Eigen::Matrix2Xd answer;

  std::lock_guard<std::mutex> lock(QHullLock);
  // call the QHull library
  {
    const int numberOfDimensions = 2;
    const int numberOfPoints = static_cast<int>(inputPoints.cols());
    const boolT isMalloc = 0;
    const char* QHullCommand = "qhull Qt";
    QHullMemoryDeallocator deallocator; // this object will handle QHull deallocation
    int returnValue = qh_new_qhull(numberOfDimensions, numberOfPoints, inputCopy.data(), isMalloc,
                                   (char *)QHullCommand, NULL, stderr);
    if(returnValue)
    {
      throw std::runtime_error("Failed to compute convex hull.");
    }

    // qh_triangulate turns rectangular facets into two triangular facets, for example
    qh_triangulate();

    // traverse the facets in counterclockwise order, yielding a final result set
    // this traversal code is translated from https://github.com/scipy/scipy/blob/master/scipy/spatial/qhull.pyx
    // in particular, the _get_extremes_2d() function
#pragma region Translated code from https://github.com/scipy/scipy/blob/master/scipy/spatial/qhull.pyx
    facetT *facet, *startfacet, *nextfacet;
    vertexT *vertexA, *vertexB;
    std::vector<ptrdiff_t> hullPointIndices;

    qh_qh->visit_id += 1;
    qh_qh->vertex_visit += 1;

    facet = qh_qh->facet_list;
    startfacet = facet;
    while(facet)
    {
      if(facet->visitid == qh_qh->visit_id)
      {
        throw std::logic_error("Internal error in QHull: loop in facet list");
      }

      if(facet->toporient)
      {
        vertexA = static_cast<vertexT*>(facet->vertices->e[0].p);
        vertexB = static_cast<vertexT*>(facet->vertices->e[1].p);
        nextfacet = static_cast<facetT*>(facet->neighbors->e[0].p);
      }
      else
      {
        vertexB = static_cast<vertexT*>(facet->vertices->e[0].p);
        vertexA = static_cast<vertexT*>(facet->vertices->e[1].p);
        nextfacet = static_cast<facetT*>(facet->neighbors->e[1].p);
      }

      if(vertexA->visitid != qh_qh->vertex_visit)
      {
        vertexA->visitid = qh_qh->vertex_visit;
        hullPointIndices.push_back(static_cast<ptrdiff_t>(qh_pointid(vertexA->point)));
      }

      if(vertexB->visitid != qh_qh->vertex_visit)
      {
        vertexB->visitid = qh_qh->vertex_visit;
        hullPointIndices.push_back(static_cast<ptrdiff_t>(qh_pointid(vertexB->point)));
      }

      facet->visitid = qh_qh->visit_id;
      facet = nextfacet;

      if(facet == startfacet)
        break;
    }
#pragma endregion

    // construct the answer
    answer.resize(2, static_cast<ptrdiff_t>(hullPointIndices.size()));
    for(size_t i = 0; i < hullPointIndices.size(); i++)
    {
      answer.col(static_cast<ptrdiff_t>(i)) = inputPoints.col(hullPointIndices.at(i));
    }

    // QHull memory freed by QHullMemoryDeallocator at the end of this scoping block
  }

  // all done
  return answer;
}

Hexagon::Blade::VoronoiDiagram computeVoronoiDiagram(const Eigen::Ref<const Eigen::Matrix2Xd>& inputPoints)
{
  Eigen::Matrix2Xd inputCopy = inputPoints;
  Hexagon::Blade::VoronoiDiagram result;

  std::lock_guard<std::mutex> lock(QHullLock);
  // call the QHull library; this Voronoi code is vaguely based on
  // https://github.com/scipy/scipy/blob/master/scipy/spatial/qhull.pyx
  // in particular, the Voronoi class and the _get_voronoi_diagram() function
  {
    const int numberOfDimensions = 2;
    const int numberOfPoints = static_cast<int>(inputPoints.cols());
    const boolT isMalloc = 0;
    const  char*  QHullCommand = "qhull v Qbb Qc Qz Qt";
    QHullMemoryDeallocator deallocator; // this object will handle QHull deallocation
    int returnValue = qh_new_qhull(numberOfDimensions, numberOfPoints, inputCopy.data(), isMalloc,
                                   (char *)QHullCommand, NULL, stderr);
    if(returnValue)
    {
      throw std::runtime_error("Failed to compute Voronoi diagram.");
    }

    // start to compute the answer
    result.points = inputPoints;

    // iterate over the facet list to see which input points are nearest each Voronoi vertex;
    // the facets are Delaunay triangles, I believe
    // apparently, there are at most 2n-5 Voronoi vertices
    result.voronoiVertices = Eigen::Matrix2Xd::Constant(2, 2 * inputPoints.cols(),
                                                        std::numeric_limits<double>::quiet_NaN());
    result.mapFromVoronoiVerticesToInputPoints.resize(static_cast<size_t>(2 * inputPoints.cols()));
    facetT *facet = qh_qh->facet_list;
    ptrdiff_t vertexIndex = 0;
    while(facet && facet->next)
    {
      assert(vertexIndex < result.voronoiVertices.cols());

      // the Voronoi vertex is the facet center of the Delaunay triangle
      auto pointDeleter = [](pointT* point)
      {
        qh_memfree(point, qh_qh->center_size);
      };
      std::shared_ptr<pointT> voronoiVertex(qh_facetcenter(facet->vertices), pointDeleter);
      result.voronoiVertices.col(vertexIndex) = Eigen::Map<const Eigen::Vector2d>(voronoiVertex.get());

      // find the points at the corners of the Delaunay triangle
      assert(3 == qh_setsize(facet->vertices));
      for(ptrdiff_t i = 0; facet->vertices->e[i].p != nullptr; i++)
      {
        vertexT* point = static_cast<vertexT*>(facet->vertices->e[i].p);
        ptrdiff_t pointIndex = qh_pointid(point->point);
        assert(pointIndex >= 0 && pointIndex <= inputPoints.cols());
        result.mapFromVoronoiVerticesToInputPoints[vertexIndex][i] = pointIndex;
      }
      facet = facet->next;
      vertexIndex++;
    }
    result.voronoiVertices.conservativeResize(2, vertexIndex);
    result.mapFromVoronoiVerticesToInputPoints.resize(static_cast<size_t>(vertexIndex));

    // QHull memory freed by QHullMemoryDeallocator at the end of this scoping block
  }

  // all done
  return result;
}

Hexagon::Blade::VoronoiDiagram computeFarthestSiteVoronoiDiagram(const Eigen::Ref<const Eigen::Matrix2Xd>& inputPoints)
{
  Eigen::Matrix2Xd inputCopy = inputPoints;
  Hexagon::Blade::VoronoiDiagram result;

  std::lock_guard<std::mutex> lock(QHullLock);
  // call the QHull library; this Voronoi code is vaguely based on
  // https://github.com/scipy/scipy/blob/master/scipy/spatial/qhull.pyx
  // in particular, the Voronoi class and the _get_voronoi_diagram() function
  {
    const int numberOfDimensions = 2;
    const int numberOfPoints = static_cast<int>(inputPoints.cols());
    const boolT isMalloc = 0;
    const char* QHullCommand = "qhull v Qbb Qt Qu"; // the Qu makes the Voronoi diagram be farthest-site
    QHullMemoryDeallocator deallocator; // this object will handle QHull deallocation
    int returnValue = qh_new_qhull(numberOfDimensions, numberOfPoints, inputCopy.data(), isMalloc,
        (char *)QHullCommand, NULL, stderr);
    if(returnValue)
    {
      throw std::runtime_error("Failed to compute Voronoi diagram.");
    }

    // start to compute the answer
    result.points = inputPoints;

    // iterate over the facet list to see which input points are nearest each Voronoi vertex;
    // the facets are Delaunay triangles, I believe
    // count the facets first
    facetT *facet = qh_qh->facet_list;
    ptrdiff_t vertexIndex = 0;
    while(facet && facet->next)
    {
      facet = facet->next;
      vertexIndex++;
    }
    const ptrdiff_t numberOfVoronoiVertices = vertexIndex;

    // allocate memory
    result.voronoiVertices = Eigen::Matrix2Xd::Constant(2, numberOfVoronoiVertices,
                                                        std::numeric_limits<double>::quiet_NaN());
    result.mapFromVoronoiVerticesToInputPoints.resize(numberOfVoronoiVertices);
    facet = qh_qh->facet_list;
    vertexIndex = 0;
    while(facet && facet->next)
    {
      assert(vertexIndex < result.voronoiVertices.cols());

      // the Voronoi vertex is the facet center of the Delaunay triangle
      auto pointDeleter = [](pointT* point)
      {
        qh_memfree(point, qh_qh->center_size);
      };
      std::shared_ptr<pointT> voronoiVertex(qh_facetcenter(facet->vertices), pointDeleter);
      result.voronoiVertices.col(vertexIndex) = Eigen::Map<const Eigen::Vector2d>(voronoiVertex.get());

      // find the points at the corners of the Delaunay triangle
      assert(3 == qh_setsize(facet->vertices));
      for(ptrdiff_t i = 0; facet->vertices->e[i].p != nullptr; i++)
      {
        vertexT* point = static_cast<vertexT*>(facet->vertices->e[i].p);
        ptrdiff_t pointIndex = qh_pointid(point->point);
        assert(pointIndex >= 0 && pointIndex <= inputPoints.cols());
        //double pointDistance = (result.voronoiVertices.col(vertexIndex) - result.points.col(pointIndex)).norm();
        result.mapFromVoronoiVerticesToInputPoints[vertexIndex][i] = pointIndex;
      }
      facet = facet->next;
      vertexIndex++;
    }

    // QHull memory freed by QHullMemoryDeallocator at the end of this scoping block
  }

  // all done
  return result;
}
}


namespace Hexagon
{
namespace Blade
{

  
class SampledCurve2D
{
  std::vector<double> txyijSamples;

public:
  typedef Eigen::Matrix<double, 5, Eigen::Dynamic> TXYIJ;

  // modifiers to add or change data
  void add(const Eigen::Ref<const TXYIJ>& newData)
  {
    txyijSamples.insert(txyijSamples.end(), newData.data(), newData.data() + newData.size());
  }
  void add(const Eigen::Ref<const Eigen::VectorXd>& newT, const Eigen::Ref<const Eigen::Matrix2Xd>& newXY,
           const Eigen::Ref<const Eigen::Matrix2Xd>& newIJ)
  {
    assert(newT.size() == newXY.cols());
    assert(newT.size() == newIJ.cols());
    TXYIJ newTXYIJ(5, newT.size());
    newTXYIJ << newT.transpose(), newXY, newIJ;
    add(newTXYIJ);
  }
  void overwrite(const Eigen::Ref<const TXYIJ>& newData)
  {
    txyijSamples.assign(newData.data(), newData.data() + newData.size());
  }
  void mergeTwoSortedRanges(ptrdiff_t firstIndexOfSecondSequence)
  {
    ptrdiff_t N = txyij().cols();
    Eigen::VectorXi unsortedIndices = Eigen::VectorXi::LinSpaced(N, 0, static_cast<int>(N) - 1);
    Eigen::VectorXi sortedIndices = Eigen::VectorXi::LinSpaced(N, 0, static_cast<int>(N) - 1);
    Eigen::VectorXd tCopy = t();
    auto comparison = [&](int a, int b) { return tCopy[a] < tCopy[b]; };
    std::merge(unsortedIndices.data(), unsortedIndices.data() + firstIndexOfSecondSequence,
               unsortedIndices.data() + firstIndexOfSecondSequence, unsortedIndices.data() + N, sortedIndices.data(),
               comparison);
    TXYIJ sortedResult(5, sortedIndices.size());
    for(ptrdiff_t i = 0; i < sortedIndices.size(); i++)
    {
      sortedResult.col(i) = txyij().col(sortedIndices[i]);
    }
    overwrite(sortedResult);
  }

  // provide views into the data
  Eigen::Ref<const TXYIJ> txyij() const
  {
    return Eigen::Map<const TXYIJ>(txyijSamples.data(), 5, txyijSamples.size() / 5);
  }
  Eigen::Ref<const Eigen::VectorXd> t() const
  {
    return txyij().row(0).transpose();
  }
  Eigen::Ref<const Eigen::Matrix2Xd> xy() const
  {
    return txyij().middleRows<2>(1);
  }
  Eigen::Ref<const Eigen::Matrix2Xd> ij() const
  {
    return txyij().bottomRows<2>();
  }
};

void evaluateAndResort(const Curve<2>& curve, const Eigen::Ref<const Eigen::VectorXd>& inputT,
                       SampledCurve2D& sampledCurve)
{
  const Eigen::VectorXd t = inputT;
  Eigen::Matrix2Xd newPoints(2, inputT.size());
  Eigen::Matrix2Xd newTangents(2, inputT.size());

  // evaluate
  curve.evaluate(t.data(), t.size(), newPoints.data(), newTangents.data(), nullptr);

  // add the new values
  ptrdiff_t oldSize = sampledCurve.t().size();
  sampledCurve.add(t, newPoints, newTangents);

  // resort
  sampledCurve.mergeTwoSortedRanges(oldSize);
}

// Returns an approximation to tan(x) that doesn't diverge for any angles.
// It is a good approximation from about -1 to 1 radians,
// and a fantastic approximation from about -0.5 to 0.5 radians
Eigen::ArrayXd nondivergentTangent(const Eigen::Ref<const Eigen::ArrayXd>& x)
{
  return x + x.pow(3.0) / 3.0 + 2.0 * x.pow(5.0) / 15.0;
}

Eigen::VectorXd signedAngleBetween(const Eigen::Ref<const Eigen::Array2Xd>& vector0,
                                   const Eigen::Ref<const Eigen::Array2Xd>& vector1)
{
  Eigen::ArrayXd dot = (vector0 * vector1).colwise().sum().transpose();
  Eigen::ArrayXd cross = (vector0.row(0) * vector1.row(1) - vector0.row(1) * vector1.row(0)).transpose();
  Eigen::ArrayXd result = cross.binaryExpr(dot, std::ptr_fun(std::atan2<double, double>));
  return result;
}

Eigen::VectorXd estimatedError(const Eigen::Ref<const Eigen::Matrix2Xd>& xy,
                               const Eigen::Ref<const Eigen::Matrix2Xd>& ij, double scale)
{
  const ptrdiff_t N = xy.cols();
  assert(N == ij.cols());
  Eigen::Matrix2Xd delta = xy.rightCols(N - 1) - xy.leftCols(N - 1);
  Eigen::ArrayXd theta0 = signedAngleBetween(delta, ij.leftCols(N - 1));
  Eigen::ArrayXd theta1 = signedAngleBetween(delta, ij.rightCols(N - 1));
  return 0.5 * delta.colwise().norm().array().transpose() *
             (nondivergentTangent(theta0).cwiseAbs() + nondivergentTangent(theta1).cwiseAbs()) +
         delta.colwise().squaredNorm().array().transpose() / scale;
}

Eigen::VectorXd refinedTValues(const Eigen::Ref<const Eigen::VectorXd>& t,
                               const Eigen::Ref<const Eigen::VectorXd>& errorEstimates, double targetError,
                               ptrdiff_t refinementRatio)
{
  Eigen::VectorXd answer;
  for(ptrdiff_t i = 0; i < errorEstimates.size(); i++)
  {
    if(errorEstimates[i] > targetError)
    {
      answer.conservativeResize(answer.size() + refinementRatio);
      answer.tail(refinementRatio) =
          Eigen::VectorXd::LinSpaced(refinementRatio + 2, t[i], t[i + 1]).segment(1, refinementRatio);
    }
  }
  return answer;
}

struct MemoizedPolygon
{
  SampledCurve2D sampledCurve;
  double targetError;
};

double length(const Eigen::Ref<const Eigen::Matrix2Xd>& points)
{
  const ptrdiff_t N = points.cols();
  const Eigen::Matrix2Xd delta = points.rightCols(N - 1) - points.leftCols(N - 1);
  return delta.colwise().norm().sum();
}

std::shared_ptr<MemoizedPolygon> corePolygonalize(const Curve<2>& curve, const ptrdiff_t initialN, double targetError)
{
  // construct the initial sample
  Eigen::Vector2d tBounds;
  curve.parametricBounds(tBounds.data());
  Eigen::VectorXd t0 = Eigen::VectorXd::LinSpaced(initialN, tBounds[0], tBounds[1]);
  SampledCurve2D initialResult;
  evaluateAndResort(curve, t0, initialResult);

  // estimate the errors
  const double scale = length(initialResult.xy());
  Eigen::ArrayXd errors = estimatedError(initialResult.xy(), initialResult.ij(), scale);

  // loop until the errors are sufficiently small
  while((errors > targetError).any())
  {
    const ptrdiff_t N = initialResult.t().size();
    Eigen::ArrayXd possibleNewT = 0.5 * (initialResult.t().segment(0, N - 1) + initialResult.t().segment(1, N - 1));
    Eigen::Array<bool, Eigen::Dynamic, 1> useNewT = errors > targetError;
    Eigen::ArrayXd newT = sliceVector(possibleNewT, useNewT);
    evaluateAndResort(curve, newT, initialResult);
    errors = estimatedError(initialResult.xy(), initialResult.ij(), scale);
  }

  // for periodic curves, remove the last sample
  const ptrdiff_t N = initialResult.t().size();
  auto finalResult = std::make_shared<MemoizedPolygon>();
  finalResult->targetError = targetError;
  if(curve.isPeriodic())
  {
    finalResult->sampledCurve.add(initialResult.txyij().leftCols(N - 1));
  }
  else
  {
    finalResult->sampledCurve.add(initialResult.txyij());
  }

  // all done
  return finalResult;
}

std::shared_ptr<const MemoizedPolygon> getMemoizedPolygon(const Curve<2>& curve, const ptrdiff_t initialN,
                                                          double targetError)
{
  // recompute
  auto memoizedPolygon = corePolygonalize(curve, initialN, targetError);

  //// save the recomputed thing
  //memoizedPolygons.emplace(curve, memoizedPolygon);

  // all done
  return memoizedPolygon;
}

std::tuple<Eigen::Matrix2Xd, Eigen::VectorXd> polygonalizeWithT(const Curve<2>& curve, const ptrdiff_t initialN,
                                                                double targetError)
{
  auto polygon = getMemoizedPolygon(curve, initialN, targetError);
  return std::make_tuple(Eigen::Matrix2Xd(polygon->sampledCurve.xy()), Eigen::VectorXd(polygon->sampledCurve.t()));
}

Eigen::Matrix2Xd polygonalize(const Curve<2>& curve, const ptrdiff_t initialN, double targetError)
{
  auto polygon = getMemoizedPolygon(curve, initialN, targetError);
  return polygon->sampledCurve.xy();
}



// this is the center and diameter of the max-inscribed circle within a blade section
ThicknessResult maxThickness(const Curve<2>& curve, const Curve<2>& convexCurve, const Curve<2>& concaveCurve,
                             const ptrdiff_t initialN, double targetError)
{
  auto polygon = getMemoizedPolygon(curve, initialN, targetError);
  auto voronoiDiagram = QHullInternals::computeVoronoiDiagram(polygon->sampledCurve.xy());
  Eigen::VectorXd t = polygon->sampledCurve.t();
  Eigen::Matrix2Xd points = polygon->sampledCurve.xy();

  // label the points according to their subcurves
  // -1 for concave side, +1 for convex side, and 0 for neither
  Eigen::VectorXi sideLabel(t.size());
  for(ptrdiff_t i = 0; i < t.size(); i++)
  {
    if(tIsInSubcurve(t[i], convexCurve, curve.period()))
    {
      sideLabel[i] = 1;
      continue;
    }
    if(tIsInSubcurve(t[i], concaveCurve, curve.period()))
    {
      sideLabel[i] = -1;
      continue;
    }
    sideLabel[i] = 0;
  }

  // construct a KD tree of the points under consideration
  typedef nanoflann::KDTreeEigenMatrixAdaptor<Eigen::MatrixX2d, 2> KDTree;
  Eigen::MatrixX2d pointsTranspose = points.transpose();
  KDTree kdtree(2 /*dimensions*/, pointsTranspose, 10 /*max leaf*/);

  // look through the Voronoi vertices, finding the one that is maximal size, within the subset
  // that is nearest to both sides

  ptrdiff_t indexOfMaxDistance = -1;
  double maxDistance = -1.0;
  for(ptrdiff_t i = 0; i < voronoiDiagram.voronoiVertices.cols(); i++)
  {
    // look at the neighbors
    const auto& neighborSet = voronoiDiagram.mapFromVoronoiVerticesToInputPoints[i];
    std::vector<int> vertexSides;
    for(const auto& pointIndex : neighborSet)
    {
      if(pointIndex < sideLabel.size())
      {
        vertexSides.push_back(sideLabel[pointIndex]);
      }
    }
    bool neighborsConvex = std::any_of(vertexSides.begin(), vertexSides.end(), [](const int& label){ return label == 1; });
    bool neighborsConcave = std::any_of(vertexSides.begin(), vertexSides.end(), [](const int& label){ return label == -1; });
    if(!neighborsConcave || !neighborsConvex)
    {
      continue;
    }
    // see if this is a current-max-thickness point
    double curDistance = (voronoiDiagram.voronoiVertices.col(i) -
                          voronoiDiagram.points.col(*neighborSet.begin())).norm();
    if(curDistance <= maxDistance)
    {
      continue;
    }
    // double-check to make sure the neighbors are closest
    // but only check if all the other checks pass, because the KD tree query is relatively slow
    ptrdiff_t closestIndex;
    double closestDistanceSquared;
    kdtree.query(voronoiDiagram.voronoiVertices.col(i).data(), 1, &closestIndex, &closestDistanceSquared);
    bool neighborsAreClosest = sqrt(closestDistanceSquared) >= (1.0 - 1e-6) * curDistance;
    if(!neighborsAreClosest)
    {
      continue;
    }
    // this is a best-so-far max-thickness point
    maxDistance = curDistance;
    indexOfMaxDistance = i;
  }

  // all done
  ThicknessResult result;
  result.thicknessCenter = voronoiDiagram.voronoiVertices.col(indexOfMaxDistance);
  result.thickness = 2.0 * maxDistance;
  return result;
}

// this is the center and diameter of the two-point max thickness within a blade section
ThicknessResult maxThicknessTwoPoint(const Curve<2>& curve, const Curve<2>& convexCurve, const Curve<2>& concaveCurve,
                                     const ptrdiff_t initialN, double targetError)
{
  auto polygon = getMemoizedPolygon(curve, initialN, targetError);
  auto voronoiDiagram = QHullInternals::computeVoronoiDiagram(polygon->sampledCurve.xy());
  Eigen::VectorXd t = polygon->sampledCurve.t();
  Eigen::Matrix2Xd points = polygon->sampledCurve.xy();

  // label the points according to their subcurves
  // -1 for concave side, +1 for convex side, and 0 for neither
  Eigen::VectorXi sideLabel(t.size());
  for(ptrdiff_t i = 0; i < t.size(); i++)
  {
    if(tIsInSubcurve(t[i], convexCurve, curve.period()))
    {
      sideLabel[i] = 1;
      continue;
    }
    if(tIsInSubcurve(t[i], concaveCurve, curve.period()))
    {
      sideLabel[i] = -1;
      continue;
    }
    sideLabel[i] = 0;
  }

  // construct a KD tree of the points under consideration
  typedef nanoflann::KDTreeEigenMatrixAdaptor<Eigen::MatrixX2d, 2> KDTree;
  Eigen::MatrixX2d pointsTranspose = points.transpose();
  KDTree kdtree(2 /*dimensions*/, pointsTranspose, 10 /*max leaf*/);

  // look through the Voronoi vertices, finding the one that is maximal size, within the subset
  // that is nearest to both sides

  ptrdiff_t indexOfMaxDistance = -1;
  double maxDistance = -1.0;
  for(ptrdiff_t i = 0; i < voronoiDiagram.voronoiVertices.cols(); i++)
  {
    // look at the neighbors
    const auto& neighborSet = voronoiDiagram.mapFromVoronoiVerticesToInputPoints[i];
    std::vector<int> vertexSides;
    for(const auto& pointIndex : neighborSet)
    {
      if(pointIndex < sideLabel.size())
      {
        vertexSides.push_back(sideLabel[pointIndex]);
      }
    }
    bool neighborsConvex = std::any_of(vertexSides.begin(), vertexSides.end(), [](const int& label){ return label == 1; });
    bool neighborsConcave = std::any_of(vertexSides.begin(), vertexSides.end(), [](const int& label){ return label == -1; });
    if(!neighborsConcave || !neighborsConvex)
    {
      continue;
    }
    // see if this is a current-max-thickness point
    double curDistance = (voronoiDiagram.voronoiVertices.col(i) -
                          voronoiDiagram.points.col(*neighborSet.begin())).norm();
    if(curDistance <= maxDistance)
    {
      continue;
    }
    // double-check to make sure the neighbors are closest
    // but only check if all the other checks pass, because the KD tree query is relatively slow
    ptrdiff_t closestIndex;
    double closestDistanceSquared;
    kdtree.query(voronoiDiagram.voronoiVertices.col(i).data(), 1, &closestIndex, &closestDistanceSquared);
    bool neighborsAreClosest = sqrt(closestDistanceSquared) >= (1.0 - 1e-6) * curDistance;
    if(!neighborsAreClosest)
    {
      continue;
    }
    // this is a best-so-far max-thickness point
    maxDistance = curDistance;
    indexOfMaxDistance = i;
  }

  // where are the curve touches?
  const auto& neighborSet = voronoiDiagram.mapFromVoronoiVerticesToInputPoints[indexOfMaxDistance];
  std::vector<int> vertexSides;
  for(const auto& pointIndex : neighborSet)
  {
    if(pointIndex < sideLabel.size())
    {
      vertexSides.push_back(sideLabel[pointIndex]);
    }
  }
  const Curve<2> *curve1, *curve2; // the curve where there is a 1-point touch and a 2-point touch respectively
  std::array<ptrdiff_t, 2> touch2index;
  ptrdiff_t touch1index;
  if(std::count(vertexSides.begin(), vertexSides.end(), 1) == 2)
  {
    // the convex side touches twice
    curve2 = &convexCurve;
    curve1 = &concaveCurve;

    // sort out the touches
    ptrdiff_t touch1indexIndex = std::find(vertexSides.begin(), vertexSides.end(), -1) - vertexSides.begin();
    ptrdiff_t touch2firstIndexIndex = std::find(vertexSides.begin(), vertexSides.end(), 1) - vertexSides.begin();
    ptrdiff_t touch2secondIndexIndex =
        std::find(vertexSides.begin() + touch2firstIndexIndex + 1, vertexSides.end(), 1) - vertexSides.begin();
    touch1index = neighborSet[touch1indexIndex];
    touch2index[0] = neighborSet[touch2firstIndexIndex];
    touch2index[1] = neighborSet[touch2secondIndexIndex];
  }
  else if(std::count(vertexSides.begin(), vertexSides.end(), -1) == 2)
  {
    // the concave side touches twice
    curve2 = &concaveCurve;
    curve1 = &convexCurve;

    // sort out the touches
    ptrdiff_t touch1indexIndex = std::find(vertexSides.begin(), vertexSides.end(), 1) - vertexSides.begin();
    ptrdiff_t touch2firstIndexIndex = std::find(vertexSides.begin(), vertexSides.end(), -1) - vertexSides.begin();
    ptrdiff_t touch2secondIndexIndex =
        std::find(vertexSides.begin() + touch2firstIndexIndex + 1, vertexSides.end(), -1) - vertexSides.begin();
    touch1index = neighborSet[touch1indexIndex];
    touch2index[0] = neighborSet[touch2firstIndexIndex];
    touch2index[1] = neighborSet[touch2secondIndexIndex];
  }
  else
  {
    throw std::logic_error("This should be impossible.");
  }

  // on the curve where the ball touches in 2 places, try to expand the ball to touch in only 1 place
  const Eigen::Vector2d pointA = points.col(touch1index);
  const Eigen::Vector2d pointB = 0.5 * (points.col(touch2index[0]) + points.col(touch2index[1]));
  Eigen::Vector2d expansionDirection = (pointB - pointA).normalized();
  Eigen::Vector2d tBounds(t[touch2index[0]], t[touch2index[1]]);
  std::sort(tBounds.data(), tBounds.data() + tBounds.size());
  double newTouchT;
  Eigen::Vector2d newTouchXY;
  if(extreme(curve, expansionDirection.data(), &newTouchT, newTouchXY.data(), tBounds[0], tBounds[1]))
  {
    // we found an extreme point within the two-point window, so use it
    ThicknessResult twoPointResult;
    twoPointResult.thicknessCenter = 0.5 * (pointA + newTouchXY);
    twoPointResult.thickness = (pointA - newTouchXY).norm();
    return twoPointResult;
  }
  else
  {
    // no extreme point found; use the max-inscribed thickness
    ThicknessResult twoPointResult;
    twoPointResult.thicknessCenter = 0.5 * (pointA + pointB);
    twoPointResult.thickness = (pointA - pointB).norm();
    return twoPointResult;
  }
}


double mean(const std::vector<double>& values)
{
  Eigen::Map<const Eigen::VectorXd> eigenValues(values.data(), static_cast<ptrdiff_t>(values.size()));
  return eigenValues.mean();
}

// the following function is inspired by
// http://stackoverflow.com/questions/1903954/is-there-a-standard-sign-function-signum-sgn-in-c-c
double sign(double val)
{
  return static_cast<int>(0.0 < val) - static_cast<int>(val < 0.0);
}

struct PreciseCamberPoint
{
  bool preciseCamberPointFound;
  Eigen::VectorXd camberPoint;
  Eigen::VectorXd camberTangent;
  double camberRadius;
  Eigen::VectorXd neighboringTValues;
};
PreciseCamberPoint computePreciseCamberPoint(
    const Curve<2>& curve, const Curve<2>& convexCurve, const Curve<2>& concaveCurve,
    const Eigen::Ref<const Eigen::Vector2d>& guessPoint, const std::array<ptrdiff_t, 3>& neighboringPointIndices,
    const Eigen::Ref<const Eigen::Matrix2Xd>& /*points*/, const Eigen::Ref<const Eigen::VectorXd>& pointsTValues,
    const nanoflann::KDTreeEigenMatrixAdaptor<Eigen::MatrixX2d, 2>& kdtree, const double polygonSign)
{
  // create a non-result that can be used to return early
  PreciseCamberPoint noResult;
  noResult.preciseCamberPointFound = false;

  // calculate the best-so-far points and tangents
  Eigen::Matrix<ptrdiff_t, 3, 1> neighboringIndices(neighboringPointIndices[0], neighboringPointIndices[1],
                                                    neighboringPointIndices[2]);
  Eigen::Vector3d neighboringTValues(pointsTValues[neighboringPointIndices[0]],
                                     pointsTValues[neighboringPointIndices[1]],
                                     pointsTValues[neighboringPointIndices[2]]);
  Eigen::Matrix2Xd neighboringPoints(2, 3);
  Eigen::Matrix2Xd neighboringTangents(2, 3);
  curve.evaluate(neighboringTValues.data(), 3, neighboringPoints.data(), neighboringTangents.data(), nullptr);

  // figure out which points are concave vs convex
  auto isWithinConvex = [&](double t) -> bool 
  {
    return tIsInSubcurve(t, convexCurve, curve.period());
  };
  auto isWithinConcave = [&](double t) -> bool 
  {
    return tIsInSubcurve(t, concaveCurve, curve.period());
  };
  auto tWithinConvex = [&](double t) -> double
  {
    return wrapToAbove(t, convexCurve.t0(), curve.period());
  };
  auto tWithinConcave = [&](double t) -> double
  {
    return wrapToAbove(t, concaveCurve.t0(), curve.period());
  };
  Eigen::Matrix<bool, 3, 1> isConvex(isWithinConvex(neighboringTValues[0]), isWithinConvex(neighboringTValues[1]),
                                     isWithinConvex(neighboringTValues[2]));
  Eigen::Matrix<bool, 3, 1> isConcave(isWithinConcave(neighboringTValues[0]), isWithinConcave(neighboringTValues[1]),
                                      isWithinConcave(neighboringTValues[2]));
  Eigen::Matrix<ptrdiff_t, Eigen::Dynamic, 1> convexNeighborIndices = sliceVector(neighboringIndices, isConvex);
  Eigen::VectorXd convexNeighborTValues = sliceVector(neighboringTValues, isConvex);
  Eigen::Matrix2Xd convexNeighborPoints = sliceColumns(neighboringPoints, isConvex);
  Eigen::Matrix2Xd convexNeighborTangents = sliceColumns(neighboringTangents, isConvex);
  Eigen::Matrix<ptrdiff_t, Eigen::Dynamic, 1> concaveNeighborIndices = sliceVector(neighboringIndices, isConcave);
  Eigen::VectorXd concaveNeighborTValues = sliceVector(neighboringTValues, isConcave);
  Eigen::Matrix2Xd concaveNeighborPoints = sliceColumns(neighboringPoints, isConcave);
  Eigen::Matrix2Xd concaveNeighborTangents = sliceColumns(neighboringTangents, isConcave);

  // figure out which way to allow translation of the guess point
  Eigen::Vector2d tangentDirection =
      (convexNeighborTangents.rowwise().mean().normalized() - concaveNeighborTangents.rowwise().mean().normalized())
          .normalized();
  Eigen::Vector2d translationDirection(tangentDirection[1], -tangentDirection[0]);
  auto translateCenter = [&](double t) -> Eigen::Vector2d { return guessPoint + t * translationDirection; };
  Eigen::Vector2d centerPoint = translateCenter(0.0);

  // precisely determine closest points
  Eigen::VectorXd convexPreciseNeighborTValues = convexNeighborTValues;
  Eigen::Matrix2Xd convexPreciseNeighborPoints(2, convexNeighborPoints.cols());
  Eigen::Matrix2Xd convexPreciseNeighborTangents(2, convexNeighborTangents.cols());
  for(ptrdiff_t i = 0; i < convexNeighborTValues.size(); i++)
  {
    convexPreciseNeighborTValues[i] = tWithinConvex(convexPreciseNeighborTValues[i]);
    closestPoint(convexCurve, centerPoint.data(), convexPreciseNeighborPoints.col(i).data(),
                 &convexPreciseNeighborTValues[i], convexPreciseNeighborTangents.col(i).data(),
                 tWithinConvex(pointsTValues[std::max<ptrdiff_t>(0, convexNeighborIndices[i] - 1)]),
                 tWithinConvex(pointsTValues[std::min(pointsTValues.size() - 1, convexNeighborIndices[i] + 1)]));
    if(!isWithinConvex(convexPreciseNeighborTValues[i]))
    {
      return noResult;
    }
  }
  Eigen::VectorXd concavePreciseNeighborTValues = concaveNeighborTValues;
  Eigen::Matrix2Xd concavePreciseNeighborPoints(2, concaveNeighborPoints.cols());
  Eigen::Matrix2Xd concavePreciseNeighborTangents(2, concaveNeighborTangents.cols());
  for(ptrdiff_t i = 0; i < concaveNeighborTValues.size(); i++)
  {
    concavePreciseNeighborTValues[i] = tWithinConcave(concavePreciseNeighborTValues[i]);
    closestPoint(concaveCurve, centerPoint.data(), concavePreciseNeighborPoints.col(i).data(),
                 &concavePreciseNeighborTValues[i], concavePreciseNeighborTangents.col(i).data(),
                 tWithinConcave(pointsTValues[std::max<ptrdiff_t>(0, concaveNeighborIndices[i] - 1)]),
                 tWithinConcave(pointsTValues[std::min(pointsTValues.size() - 1, concaveNeighborIndices[i] + 1)]));
    if(!isWithinConcave(concavePreciseNeighborTValues[i]))
    {
      return noResult;
    }
  }
  Eigen::VectorXd convexDistances = (convexPreciseNeighborPoints.colwise() - centerPoint).colwise().norm().transpose();
  Eigen::VectorXd concaveDistances =
      (concavePreciseNeighborPoints.colwise() - centerPoint).colwise().norm().transpose();
  double sphereRadius = std::min(concaveDistances.minCoeff(), convexDistances.minCoeff());

  // make a function that has a zero when the center point is equidistant from the two sides
  auto centerPointZeroingFunction = [&](double t) -> double
  {
    // this function has a zero when the center point is equidistant from the two sides
    Eigen::Vector2d point = translateCenter(t);
    Eigen::VectorXd convexDistances = (convexPreciseNeighborPoints.colwise() - point).colwise().norm().transpose();
    Eigen::VectorXd concaveDistances = (concavePreciseNeighborPoints.colwise() - point).colwise().norm().transpose();
    return concaveDistances.minCoeff() - convexDistances.minCoeff();
  };

  // move the center point to match the minimum concave and convex distances
  centerPoint = translateCenter(Hexagon::MetrologyBuildingBlocks::scalarZero(
      centerPointZeroingFunction, -0.5 * sphereRadius, 0.0, 0.5 * sphereRadius,
      Hexagon::MetrologyBuildingBlocks::rationalInterpolationStep));

  // readjust sphereRadius
  convexDistances = (convexPreciseNeighborPoints.colwise() - centerPoint).colwise().norm().transpose();
  concaveDistances = (concavePreciseNeighborPoints.colwise() - centerPoint).colwise().norm().transpose();
  sphereRadius = std::min(concaveDistances.minCoeff(), convexDistances.minCoeff());

  // OK, I think we're done optimizing here
  // double-check the KD-tree to make sure the sphere radius and the KD-tree's closest point are the same thing
  ptrdiff_t closestPointIndex;
  double distanceSquared;
  kdtree.query(centerPoint.data(), 1, &closestPointIndex, &distanceSquared);
  const double sphereRadiusCutoff = (1.0 - 1e-10) * std::sqrt(distanceSquared);
  if(sphereRadius > sphereRadiusCutoff)
  {
    return noResult;
  }

  // look at the closest tangents and make sure the convex sides are really orthogonal to the tangent
  Eigen::Matrix2Xd convexPreciseDelta = convexPreciseNeighborPoints.colwise() - centerPoint;
  convexPreciseDelta.colwise().normalize();
  convexPreciseNeighborTangents.colwise().normalize();
  Eigen::Matrix2Xd convexPreciseNeighborNormals(2, convexNeighborTangents.cols());
  convexPreciseNeighborNormals.row(0) = polygonSign * convexPreciseNeighborTangents.row(1);
  convexPreciseNeighborNormals.row(1) = -polygonSign * convexPreciseNeighborTangents.row(0);
  if(((convexPreciseDelta.array() * convexPreciseNeighborNormals.array()).colwise().sum() < 0.3).any())
  {
    // the center-point is not "inside" the curve
    return noResult;
  }
  Eigen::VectorXd convexOrthogonalities =
      (convexPreciseDelta.array() * convexPreciseNeighborTangents.array()).colwise().sum().cwiseAbs().transpose();
  ptrdiff_t bestConvexPoint;
  const double bestConvexOrthogonality = convexOrthogonalities.minCoeff(&bestConvexPoint);
  if(bestConvexOrthogonality > 1e-4)
  {
    // the inscribed circle does not precisely tangent-contact the curve
    return noResult;
  }

  // look at the closest tangents and make sure both the concave sides are really orthogonal to the tangent
  Eigen::Matrix2Xd concavePreciseDelta = concavePreciseNeighborPoints.colwise() - centerPoint;
  concavePreciseDelta.colwise().normalize();
  concavePreciseNeighborTangents.colwise().normalize();
  Eigen::Matrix2Xd concavePreciseNeighborNormals(2, concaveNeighborTangents.cols());
  concavePreciseNeighborNormals.row(0) = polygonSign * concavePreciseNeighborTangents.row(1);
  concavePreciseNeighborNormals.row(1) = -polygonSign * concavePreciseNeighborTangents.row(0);
  if(((concavePreciseDelta.array() * concavePreciseNeighborNormals.array()).colwise().sum() < 0.3).any())
  {
    // the center-point is not "inside" the curve
    return noResult;
  }
  Eigen::VectorXd concaveOrthogonalities =
      (concavePreciseDelta.array() * concavePreciseNeighborTangents.array()).colwise().sum().cwiseAbs().transpose();
  ptrdiff_t bestConcavePoint;
  const double bestConcaveOrthogonality = concaveOrthogonalities.minCoeff(&bestConcavePoint);
  if(bestConcaveOrthogonality > 1e-4)
  {
    // the inscribed circle does not precisely tangent-contact the curve
    return noResult;
  }

  // we're happy with the optimized result; report it
  PreciseCamberPoint result;
  result.camberPoint = centerPoint;
  result.camberRadius = sphereRadius;
  result.camberTangent =
      (convexPreciseNeighborTangents.col(bestConvexPoint) - concavePreciseNeighborTangents.col(bestConcavePoint))
          .normalized();
  result.neighboringTValues.resize(convexPreciseNeighborTValues.size() + concavePreciseNeighborTValues.size());
  result.neighboringTValues.head(convexPreciseNeighborTValues.size()) = convexPreciseNeighborTValues;
  result.neighboringTValues.tail(concavePreciseNeighborTValues.size()) = concavePreciseNeighborTValues;
  result.preciseCamberPointFound = true;
  return result;
}

Camber camber(const Curve<2>& curve, const Curve<2>& convexCurve, const Curve<2>& concaveCurve,
              const ptrdiff_t initialN, double targetError, const Eigen::Ref<const Eigen::Vector2d>& halfSpace1Point,
              const Eigen::Ref<const Eigen::Vector2d>& halfSpace1Vector,
              const Eigen::Ref<const Eigen::Vector2d>& halfSpace2Point,
              const Eigen::Ref<const Eigen::Vector2d>& halfSpace2Vector)
{
  // access core polygonalization
  auto polygon = getMemoizedPolygon(curve, initialN, targetError);
  Eigen::VectorXd polygon_t = polygon->sampledCurve.t();
  Eigen::Matrix2Xd polygon_points = polygon->sampledCurve.xy();
  Eigen::Matrix2Xd polygon_tangents = polygon->sampledCurve.ij();

  // figure out whether the polygon is clockwise or counterclockwise
  double polygonSign = sign(signedArea(polygon_points));

  // delete leading / trailing edges
  auto isWithinConvex = [&](double t) -> bool 
  {
    return tIsInSubcurve(t, convexCurve, curve.period());
  };
  auto isWithinConcave = [&](double t) -> bool 
  {
    return tIsInSubcurve(t, concaveCurve, curve.period());
  };
  typedef Eigen::Matrix<bool, Eigen::Dynamic, 1> BoolVector;
  BoolVector isPartOfVoronoi(polygon_t.size());
  for(ptrdiff_t i = 0; i < polygon_t.size(); i++)
  {
    isPartOfVoronoi[i] = isWithinConvex(polygon_t[i]) || isWithinConcave(polygon_t[i]);
  }
  ptrdiff_t pointsN = isPartOfVoronoi.count();
  Eigen::VectorXd pointsTValues(pointsN);
  Eigen::Matrix2Xd points(2, pointsN);
  Eigen::Matrix2Xd tangents(2, pointsN);
  ptrdiff_t k = 0;
  for(ptrdiff_t i = 0; i < polygon_t.size(); i++)
  {
    if(isPartOfVoronoi[i])
    {
      pointsTValues[k] = polygon_t[i];
      points.col(k) = polygon_points.col(i);
      tangents.col(k) = polygon_tangents.col(i);
      k++;
    }
  }
  //alwaysAssert(k == pointsN);

  // label the points according to their subcurves
  // -1 for concave side, +1 for convex side, and 0 for neither
  Eigen::VectorXi sideLabel(pointsN);
  for(ptrdiff_t i = 0; i < pointsN; i++)
  {
    if(isWithinConvex(pointsTValues[i]))
    {
      sideLabel[i] = 1;
      continue;
    }
    if(isWithinConcave(pointsTValues[i]))
    {
      sideLabel[i] = -1;
      continue;
    }
    sideLabel[i] = 0;
  }


  // do the Voronoi diagram
  auto voronoiDiagram = QHullInternals::computeVoronoiDiagram(points);


  // look through the Voronoi vertices, finding the ones that touch both sides
  Eigen::Matrix<bool, Eigen::Dynamic, 1> isCamber(voronoiDiagram.voronoiVertices.cols());
  Eigen::VectorXd closestDistances(voronoiDiagram.voronoiVertices.cols());
  typedef nanoflann::KDTreeEigenMatrixAdaptor<Eigen::MatrixX2d, 2> KDTree;
  Eigen::MatrixX2d pointsTranspose = points.transpose();
  KDTree kdtree(2 /*dimensions*/, pointsTranspose, 10 /*max leaf*/);
  for(ptrdiff_t i = 0; i < voronoiDiagram.voronoiVertices.cols(); i++)
  {
    // look at the neighbors
    Eigen::Vector2d candidateVertex = voronoiDiagram.voronoiVertices.col(i);
    const auto& neighborSet = voronoiDiagram.mapFromVoronoiVerticesToInputPoints[i];
    std::vector<double> distancesToPoints;
    std::vector<double> normalParallelisms;
    std::vector<int> closestPointSides;
    bool neighborsAreNormalParallel = true;
    for(const auto& pointIndex : neighborSet)
    {
      if(pointIndex < sideLabel.size())
      {
        closestPointSides.push_back(sideLabel[pointIndex]);
        Eigen::Vector2d delta = candidateVertex - points.col(pointIndex);
        distancesToPoints.push_back(delta.norm());
        Eigen::Vector2d surfaceNormal(-polygonSign * tangents(1, pointIndex), polygonSign * tangents(0, pointIndex));
        double normalParallelism = delta.normalized().dot(surfaceNormal.normalized());
        neighborsAreNormalParallel = neighborsAreNormalParallel && normalParallelism > 0.9;
        normalParallelisms.push_back(normalParallelism);
      }
    }

    // make sure that three actual points are neighbors
    if(distancesToPoints.size() < 3)
    {
      isCamber[i] = false;
      continue;
    }

    // make sure the distances are all similar
    bool distancesAreSimilar1 = std::abs(distancesToPoints[0] - distancesToPoints[1]) < 1e-6 * distancesToPoints[0];
    bool distancesAreSimilar2 = std::abs(distancesToPoints[0] - distancesToPoints[2]) < 1e-6 * distancesToPoints[0];

    // make sure both curves are neighbors
    bool neighborsConvex =
      std::any_of(closestPointSides.begin(), closestPointSides.end(), [](const int& label) { return label == 1; });
    bool neighborsConcave =
      std::any_of(closestPointSides.begin(), closestPointSides.end(), [](const int& label) { return label == -1; });

    // make sure it's within the half-spaces
    bool isWithinHalfSpace1 = (candidateVertex - halfSpace1Point).dot(halfSpace1Vector) < 0.0;
    bool isWithinHalfSpace2 = (candidateVertex - halfSpace2Point).dot(halfSpace2Vector) < 0.0;

    // now we can say whether the candidate is a camber-point or not
    isCamber[i] = distancesAreSimilar1 && distancesAreSimilar2 && neighborsAreNormalParallel && neighborsConcave &&
                  neighborsConvex && isWithinHalfSpace1 && isWithinHalfSpace2;

    // double-check to make sure the neighbors are closest
    // but only check if all the other checks pass, because the KD tree query is relatively slow
    if(isCamber[i])
    {
      ptrdiff_t closestIndex;
      double closestDistanceSquared;
      kdtree.query(candidateVertex.data(), 1, &closestIndex, &closestDistanceSquared);
      const double closestDistance = sqrt(closestDistanceSquared);
      bool neighborsAreClosest = closestDistance >= (1.0 - 1e-6) * distancesToPoints[0];
      isCamber[i] = isCamber[i] && neighborsAreClosest;
      closestDistances[i] = closestDistance;
    }
  }

  // compute the more-precise camber points and tangents and neighboring T-values
  std::vector<PreciseCamberPoint> preciseCamberPoints;
  for(ptrdiff_t i = 0; i < voronoiDiagram.voronoiVertices.cols(); i++)
  {
    if(isCamber[i])
    {
      // compute the precise value of the camber point
      auto preciseResult = computePreciseCamberPoint(
          curve, convexCurve, concaveCurve, voronoiDiagram.voronoiVertices.col(i),
          voronoiDiagram.mapFromVoronoiVerticesToInputPoints[i], points, pointsTValues, kdtree, polygonSign);
      if(preciseResult.preciseCamberPointFound)
      {
        preciseCamberPoints.push_back(preciseResult);
      }
    }
  }

  // assemble the (unsorted) results structure
  Camber unsortedResult;
  ptrdiff_t N = preciseCamberPoints.size();
  //alwaysAssert(N > 0);
  unsortedResult.camberPoints.resize(2, N);
  unsortedResult.camberRadii.resize(N);
  unsortedResult.camberTangents.resize(2, N);
  unsortedResult.closestTValues.resize(3, N);
  for(ptrdiff_t i = 0; i < N; i++)
  {
    unsortedResult.camberPoints.col(i) = preciseCamberPoints[i].camberPoint;
    unsortedResult.camberRadii[i] = preciseCamberPoints[i].camberRadius;
    unsortedResult.camberTangents.col(i) = preciseCamberPoints[i].camberTangent;
    unsortedResult.closestTValues.col(i) = preciseCamberPoints[i].neighboringTValues;
  }

  // sort the results by the (convex - concave) t-value average
  Eigen::VectorXd tConvexMinusConcave(N);
  for(ptrdiff_t i = 0; i < N; i++)
  {
    Eigen::Vector3d closestTValues = unsortedResult.closestTValues.col(i);
    std::vector<double> convexT;
    std::vector<double> concaveT;
    for(ptrdiff_t j = 0; j < closestTValues.size(); j++)
    {
      double closestTValue = closestTValues[j];
      if(isWithinConvex(closestTValue))
      {
        double tWithinConvex = wrapToAbove(closestTValue, convexCurve.t0(), curve.period());
        convexT.push_back(tWithinConvex);
      }
      if(isWithinConcave(closestTValue))
      {
        double tWithinConcave = wrapToAbove(closestTValue, concaveCurve.t0(), curve.period());
        concaveT.push_back(tWithinConcave);
      }
    }

    // check sanity
    //alwaysAssert(convexT.size() > 0);
    //alwaysAssert(convexT.size() < 3);
    //alwaysAssert(concaveT.size() > 0);
    //alwaysAssert(concaveT.size() < 3);

    tConvexMinusConcave[i] = mean(convexT) - mean(concaveT);
  }
  typedef Eigen::Matrix<ptrdiff_t, Eigen::Dynamic, 1> IndexVector;
  IndexVector argSortTConvexMinusConcave = IndexVector::LinSpaced(N, 0, N - 1);
  auto comparisonTConvexMinusConcave = [&](const ptrdiff_t& index1, const ptrdiff_t& index2)
  {
    return tConvexMinusConcave[index1] < tConvexMinusConcave[index2];
  };
  std::sort(argSortTConvexMinusConcave.data(), argSortTConvexMinusConcave.data() + N, comparisonTConvexMinusConcave);

  // the convex-minus-concave t-value average can cause some 'flipping' from point to point in certain cases
  // which yields very-bad behavior when at the ends (where it is most likely)
  IndexVector argSort = IndexVector::LinSpaced(N, 0, N - 1);
  if(N < 10)
  {
    argSort = argSortTConvexMinusConcave;
  }
  else
  {
    IndexVector coarsePointIndices = IndexVector::LinSpaced(10, 0, N - 1);
    const Eigen::Matrix2Xd coarsePoints = sliceColumnsWithIndices(
        sliceColumnsWithIndices(unsortedResult.camberPoints, argSortTConvexMinusConcave), coarsePointIndices);
    const Eigen::Matrix2Xd coarseTangents = sliceColumnsWithIndices(
        sliceColumnsWithIndices(unsortedResult.camberTangents, argSortTConvexMinusConcave), coarsePointIndices);
    Eigen::Matrix2Xd coarseNormals(2, 10);
    coarseNormals.row(0) = coarseTangents.row(1);
    coarseNormals.row(1) = -coarseTangents.row(0);
    auto coarseSpline = initialSpline<PeriodicSplineType::Nonperiodic>(coarsePoints, coarseNormals);
    // in the below, it doesnt matter if the curve is actually English or not, but I say it is English
    auto coarseCamberCurve =
       std::make_unique<HermiteOpenCurve>(coarseSpline.points, coarseSpline.tangents, coarseSpline.t, true);
    Eigen::VectorXd tValuesAssociatedWithCoarseSpline(N);
    findClosestTValues(*coarseCamberCurve, tValuesAssociatedWithCoarseSpline.data(), unsortedResult.camberPoints.data(),
                       N);
    // fix any past-the-end t-values using a simple linear approximation
    Eigen::Matrix2Xd closestPoints, closestDerivatives;
    std::tie(closestPoints, closestDerivatives) =
        evaluateWithDerivative(*coarseCamberCurve, tValuesAssociatedWithCoarseSpline);
    const Eigen::VectorXd horizontalMisalignments =
        ((unsortedResult.camberPoints - closestPoints).array() * closestDerivatives.colwise().normalized().array())
            .colwise()
            .sum()
            .transpose();
    const Eigen::VectorXd adjustmentsToTValues =
        horizontalMisalignments.array() / closestDerivatives.colwise().norm().transpose().array();
    tValuesAssociatedWithCoarseSpline += adjustmentsToTValues;
    auto finalComparison = [&](const ptrdiff_t& index1, const ptrdiff_t& index2) 
    {
      return tValuesAssociatedWithCoarseSpline[index1] < tValuesAssociatedWithCoarseSpline[index2];
    };
    std::sort(argSort.data(), argSort.data() + N, finalComparison);
  }

  // assemble the sorted result
  Camber sortedResult = unsortedResult;
  sortedResult.camberPoints = sliceColumnsWithIndices(unsortedResult.camberPoints, argSort);
  sortedResult.camberRadii = sliceVectorWithIndices(unsortedResult.camberRadii, argSort);
  sortedResult.camberTangents = sliceColumnsWithIndices(unsortedResult.camberTangents, argSort);
  sortedResult.closestTValues = sliceColumnsWithIndices(unsortedResult.closestTValues, argSort);

  // all done
  return sortedResult;
}

namespace bg = boost::geometry;
typedef bg::model::point<double, 2, bg::cs::cartesian> Point;
typedef bg::model::segment<Point> Segment;
typedef std::pair<Segment, ptrdiff_t> TreeValue; // the ptrdiff_t is the index of the first point in the segment
typedef bg::index::rtree<TreeValue, bg::index::rstar<16>> RTree;


double projectPointOntoSegmentBarycentricCoordinate(const Eigen::Ref<const Eigen::Vector2d>& point,
                                                    const Eigen::Ref<const Eigen::Vector2d>& pointOfSegment1,
                                                    const Eigen::Ref<const Eigen::Vector2d>& pointOfSegment2)
{
  // this algorithm is based on Boost; see <boost/geometry/strategies/cartesian/distance_projected_point.hpp>
  // this algorithm projects the point onto the infinite line, but returns one of the endpoints if the projection
  // is outside the limited segment
  Eigen::Vector2d v = pointOfSegment2 - pointOfSegment1;
  Eigen::Vector2d w = point - pointOfSegment1;
  double c1 = w.dot(v);
  if(c1 <= 0.0)
  {
    return 0.0; // projection to pointOfSegment1
  }
  double c2 = v.dot(v);
  if(c2 <= c1)
  {
    return 1.0; // projection to pointOfSegment2
  }
  double b = c1 / c2;
  assert(std::isfinite(b));
  assert(b <= 1.0);
  assert(b >= 0.0);
  return b; // projection onto somewhere in the middle of the segment
}

Eigen::Vector2d projectPointOntoSegment(const Eigen::Ref<const Eigen::Vector2d>& point,
                                        const Eigen::Ref<const Eigen::Vector2d>& pointOfSegment1,
                                        const Eigen::Ref<const Eigen::Vector2d>& pointOfSegment2)
{
  double alpha = projectPointOntoSegmentBarycentricCoordinate(point, pointOfSegment1, pointOfSegment2);
  return (1.0 - alpha) * pointOfSegment1 + alpha * pointOfSegment2;
}

double signedArea(const Eigen::Ref<const Eigen::Matrix2Xd>& polygon)
{
  Eigen::Ref<const Eigen::ArrayXd> x = polygon.row(0).transpose();
  Eigen::Ref<const Eigen::ArrayXd> y = polygon.row(1).transpose();

  // construct (for convenience) circularly-shifted vectors
  const Eigen::ArrayXd xShiftOne = circularlyShiftLeft(x, 1);
  const Eigen::ArrayXd yShiftOne = circularlyShiftLeft(y, 1);
  
  // compute the signed area
  return 0.5 * (x * yShiftOne - xShiftOne * y).sum();
}

Eigen::Vector2d centroid(const Eigen::Ref<const Eigen::Matrix2Xd>& polygon)
{
  ptrdiff_t numberOfPoints = polygon.cols();
  assert(numberOfPoints > 2);

  Eigen::Ref<const Eigen::ArrayXd> x = polygon.row(0).transpose();
  Eigen::Ref<const Eigen::ArrayXd> y = polygon.row(1).transpose();

  // construct (for convenience) circularly-shifted vectors
  Eigen::ArrayXd xShiftOne(numberOfPoints);
  xShiftOne.topRows(numberOfPoints - 1) = x.bottomRows(numberOfPoints - 1);
  xShiftOne[numberOfPoints - 1] = x[0];
  Eigen::ArrayXd yShiftOne(numberOfPoints);
  yShiftOne.topRows(numberOfPoints - 1) = y.bottomRows(numberOfPoints - 1);
  yShiftOne[numberOfPoints - 1] = y[0];

  Eigen::ArrayXd shiftedDifferenceProduct = x * yShiftOne - xShiftOne * y;

  // compute the signed area and the centroid
  Eigen::Vector2d result;
  double A = 0.5 * shiftedDifferenceProduct.sum();
  result[0] = ((x + xShiftOne) * shiftedDifferenceProduct).sum() / (6.0 * A);
  result[1] = ((y + yShiftOne) * shiftedDifferenceProduct).sum() / (6.0 * A);

  // all done
  return result;
}

Eigen::Matrix2Xd convexHull(const Eigen::Ref<const Eigen::Matrix2Xd>& points)
{
  return QHullInternals::computeConvexHull(points);
}
VoronoiDiagram computeVoronoiDiagram(const Eigen::Ref<const Eigen::Matrix2Xd>& inputPoints)
{
  return QHullInternals::computeVoronoiDiagram(inputPoints);
}
VoronoiDiagram computeFarthestSiteVoronoiDiagram(const Eigen::Ref<const Eigen::Matrix2Xd>& inputPoints)
{
  return QHullInternals::computeFarthestSiteVoronoiDiagram(inputPoints);
}

void savetxt(const Eigen::Ref<const Eigen::MatrixXd>& input, const std::string& fileName)
{
  std::ofstream file(fileName, std::ios::out);
  file << std::setprecision(15) << input << std::flush;
  file.close();
}

Eigen::Vector2d gravityDirection(const Eigen::Ref<const Eigen::Matrix2Xd>& polygon,
                                 const Eigen::Ref<const Eigen::Vector2d>& guessGravity)
{
  Eigen::Matrix2Xd hull = convexHull(polygon);
  Eigen::Vector2d centerOfGravity = centroid(polygon);

  //savetxt(polygon, "C:\\Users\\daniel.wilcox\\Documents\\BladeSamples\\BRK-472\\polygon.txt");
  //savetxt(hull, "C:\\Users\\daniel.wilcox\\Documents\\BladeSamples\\BRK-472\\hull.txt");
  //savetxt(centerOfGravity, "C:\\Users\\daniel.wilcox\\Documents\\BladeSamples\\BRK-472\\centerOfGravity.txt");
  //savetxt(guessGravity, "C:\\Users\\daniel.wilcox\\Documents\\BladeSamples\\BRK-472\\guessGravity.txt");

  // evaluate the distances from facets to the center of gravity
  // but only those facets that agree somewhat with the guessGravity
  Eigen::Vector2d gravitySurfaceNormal;
  double minimumDistance = std::numeric_limits<double>::max();

  for(ptrdiff_t i = 0; i < hull.cols(); i++)
  {
    // hull points are guaranteed to be in counterclockwise order
    Eigen::Vector2d tangent = hull.col((i + 1) % hull.cols()) - hull.col(i);
    //double tangentLength = tangent.norm();
    tangent.normalize();
    Eigen::Vector2d surfaceNormal; // points out from the material
    surfaceNormal[0] = tangent[1];
    surfaceNormal[1] = -tangent[0];

    if(surfaceNormal.dot(guessGravity) > 0.0)
    {
      double distance = -(centerOfGravity - hull.col(i)).dot(surfaceNormal);
      if(distance < minimumDistance)
      {
        minimumDistance = distance;
        gravitySurfaceNormal = surfaceNormal;
      }
    }
  }

  // all done
  return gravitySurfaceNormal;
}
}
}
