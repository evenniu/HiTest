#pragma once

#include <memory>
#include <Eigen/Dense>

class CCurve;
class CSection;
class CNominalSection;
class CAnalysisSect;
class CFlavor;

namespace Hexagon
{
namespace Blade
{

struct CamberBackoff
{
  std::shared_ptr<Eigen::Array2d> point;
  std::shared_ptr<Eigen::Array2d> normal;
  double backoffDistance;
};

struct MeanCamberCurveParameters2016
{
  CCurve* wholeCurve;
  std::shared_ptr<CamberBackoff> noseBackoff;
  std::shared_ptr<CamberBackoff> tailBackoff;

  // if the backoff doesn't exist, then the LEC or TEC type is EDGE_PARTIAL
  // and the nose and tail t-values should be put here
  double noseStart, noseEnd;
  double tailStart, tailEnd;
};

void writeMeanCamberCurveParameters2016(const MeanCamberCurveParameters2016* params, FILE* fp);
std::shared_ptr<const MeanCamberCurveParameters2016> readMeanCamberCurveParameters2016(FILE* fp);

CCurve* createMeanCamberCurve2016(const MeanCamberCurveParameters2016& params, bool isEnglish);

struct MeanCamberCurveParameters
{
  CSection* section;
  CAnalysisSect* analysisSection;
  CFlavor* flavor;
  double ler, ter;
  double voff, uoff;
  double mtle, mtte;
};

struct MeanCamberResult
{
  CCurve* meanCamberCurve;
  bool skipPitch;
};

MeanCamberResult createMeasuredMeanCamberCurve(const MeanCamberCurveParameters& params, bool isEnglish);
CCurve* createMeasuredMeanCamberCurve2016(CCurve* measuredWholeCurve, CCurve* measuredLeadingCurve,
                                          CCurve* measuredTrailingCurve,
                                          const MeanCamberCurveParameters2016* nominalMCLParams, bool isEnglish);
}
}