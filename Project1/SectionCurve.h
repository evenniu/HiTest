#pragma once
#include <Eigen/Dense>

class CCurve;
class CSubCurve;
class CSection;

namespace Hexagon
{
namespace Blade
{
struct SectionCurve
{
  CCurve* whole;
  CCurve* leading;
  CCurve* trailing;
  CCurve* concave;
  CCurve* convex;
  CCurve* meanCamber;
  int leType, teType;
};

SectionCurve nominalSectionCurve(const CSection* section);
//
SectionCurve measuredSectionCurve(const CSection* section);

Eigen::Vector2d gravityGuess(const SectionCurve& curve);

std::unique_ptr<CSubCurve> makeHalfCurve(const SectionCurve& curve, unsigned int which_piece);

} // namespace Blade
} // namespace Hexagon
