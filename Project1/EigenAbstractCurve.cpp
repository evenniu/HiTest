#include "stdafx.h"
#include "EigenAbstractCurve.h"

namespace Hexagon
{
namespace Blade
{

Eigen::Isometry2d toIsometry2d()
{
  //if(input.n != 2)
  //{
  //  throw std::logic_error("It does not make sense to convert a 3D CAlignment into an Eigen::Isometry2d");
  //}
  //const Eigen::Vector2d b(input.m_borig[0], input.m_borig[1]);
  //const Eigen::Vector2d m(input.m_morig[0], input.m_morig[1]);
  //Eigen::Matrix2d R;
  //R << input.m_mat[0][0], input.m_mat[0][1], input.m_mat[1][0], input.m_mat[1][1];
  Eigen::Isometry2d result = Eigen::Isometry2d::Identity();
  //result.translation() = b - R * m;
  //result.linear() = R;
  return result;
}


CAlignment toCAlignment(const Eigen::Isometry2d& nominalToMeasuredTransformation)
{
  CAlignment result;

  result.m_morig[0] = result.m_morig[1] = result.m_morig[2] = 0.0;
  result.m_borig[0] = nominalToMeasuredTransformation.translation()[0];
  result.m_borig[1] = nominalToMeasuredTransformation.translation()[1];
  result.m_borig[2] = 0.0;
  result.m_mat[0][0] = nominalToMeasuredTransformation.linear()(0, 0);
  result.m_mat[0][1] = nominalToMeasuredTransformation.linear()(0, 1);
  result.m_mat[1][0] = nominalToMeasuredTransformation.linear()(1, 0);
  result.m_mat[1][1] = nominalToMeasuredTransformation.linear()(1, 1);
  result.m_mat[0][2] = result.m_mat[1][2] = result.m_mat[2][0] = result.m_mat[2][1] = 0.0;
  result.m_mat[2][2] = 1.0;

  return result;
}
} // namespace Blade
} // namespace Hexagon