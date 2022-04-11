#include "stdafx.h"
#include "SectionCurve.h"
#include "CURVE.H"
#include "CurvePolygon.h"
#include "SUBCURVE.H"

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
//SectionCurve nominalSectionCurve(const CSection *section)
//{
//  SectionCurve answer;
// //answer.whole = section->NomCurve();
//  //answer.leading = section->NomPart(LEC);
//  //answer.trailing = section->NomPart(TEC);
//  //answer.concave = section->NomPart(CCC);
//  //answer.convex = section->NomPart(CVC);
//  //answer.meanCamber = section->NomPart(MCC);
//  //answer.leType = section->LEType();
//  //answer.teType = section->TEType();
//  return answer;
//}
//
SectionCurve measuredSectionCurve(const CSection* section)
{
  SectionCurve answer;
  answer.whole = section->MeaCurve();
  answer.leading = section->MeaPart(LEC);
  answer.trailing = section->MeaPart(TEC);
  answer.concave = section->MeaPart(CCC);
  answer.convex = section->MeaPart(CVC);
  answer.meanCamber = section->MeaPart(MCC);
  answer.leType = section->LEType();
  answer.teType = section->TEType();
  return answer;
}

Eigen::Vector2d gravityGuess(const SectionCurve& curve)
{
  Eigen::Vector2d convexMidpoint;
  curve.convex->CalcPoint(convexMidpoint.data(), 0.5 * (curve.convex->T0() + curve.convex->T1()));

  Eigen::Vector2d closestConcavePoint;
  curve.concave->ClosestPoint(convexMidpoint.data(), closestConcavePoint.data());

  // rough gravity guess
  return (closestConcavePoint - convexMidpoint).normalized();
}

std::unique_ptr<CSubCurve> makeHalfCurve(const SectionCurve& curve, unsigned int which_piece)
{
  // check validity
  if(which_piece != CVC && which_piece != CCC)
  {
    return nullptr;
  }
  const auto* curvePiece = which_piece == CVC ? curve.convex : curve.concave;

  // make sure at most one edge is partial
  if(curve.leType == EDGE_PARTIAL && curve.teType == EDGE_PARTIAL)
  {
    return nullptr;
  }

  // identify the t-points of the nominal endpoints
  double nom_t_leading = curve.leading->Extreme();
  double nom_t_trailing = curve.trailing->Extreme();

  // get the boundaries of the piece; identify a t-point on the desired half
  const double t0 = curvePiece->T0();
  const double t1 = curvePiece->T1();
  const double period = curve.whole->Period();
  const double nom_t_included = 0.5 * (t0 + t1);
  const auto distanceToRounded = [](double x) -> double { return std::abs(x - round(x)); };

  // change things if edges are partial
  if(curve.leType == EDGE_PARTIAL)
  {
    // pick the one farther away from the trailing point
    const double distance0 = distanceToRounded((t0 - nom_t_trailing) / period);
    const double distance1 = distanceToRounded((t1 - nom_t_trailing) / period);
    nom_t_leading = (distance0 > distance1) ? t0 : t1;
  }
  else if(curve.teType == EDGE_PARTIAL)
  {
    // pick the one farther away from the leading point
    const double distance0 = distanceToRounded((t0 - nom_t_leading) / period);
    const double distance1 = distanceToRounded((t1 - nom_t_leading) / period);
    nom_t_trailing = (distance0 > distance1) ? t0 : t1;
  }

  // identify a t-point on the desired half
  // create the desired nominal t-range
  std::pair<double, double> nom_t_range = FitRangeToPoint(nom_t_leading, nom_t_trailing, nom_t_included, period);

  // build the subcurves
  return std::make_unique<CSubCurve>(curve.whole, nom_t_range.first, nom_t_range.second, period);
}
} // namespace Blade
} // namespace Hexagon
