#pragma once
#include <Eigen/Dense>
#include "AbstractCurve.h"


namespace Hexagon
{
namespace Blade
{

enum class PeriodicSplineType
{
  Nonperiodic = 0,
  Periodic = 1
};

template <ptrdiff_t templateRangeDimension, PeriodicSplineType templatePeriodicType>
struct HermiteSpline
{
  static const PeriodicSplineType periodicType = templatePeriodicType;
  static const ptrdiff_t rangeDimension = templateRangeDimension;
  static_assert(rangeDimension > 0, "The range dimension must always be positive.");

  typedef Eigen::VectorXd VectorType;
  typedef Eigen::Ref<VectorType> VectorRefType;
  typedef const Eigen::Ref<const VectorType>& ConstVectorRefType;

  typedef Eigen::Matrix<double, rangeDimension, Eigen::Dynamic> MatrixType;
  typedef Eigen::Ref<MatrixType> MatrixRefType;
  typedef const Eigen::Ref<const MatrixType>& ConstMatrixRefType;

  VectorType t;
  MatrixType points;
  MatrixType tangents;
};

template <class ArrayType>
bool isIncreasing(const ArrayType& t)
{
  auto i = t.size();
  for(i = 1; i < t.size(); i++)
  {
    if(t[i] <= t[i - 1])
      return false;
  }
  return true;
}

template <class SplineType>
double computeLeftDomain(const SplineType& spline)
{
  assert(spline.t.size() > 0);
  return spline.t[0];
}

template <class SplineType>
double computeRightDomain(const SplineType& spline)
{
  assert(spline.t.size() > 0);
  return spline.t[spline.t.size() - 1];
}

template <class SplineType>
double computePeriod(const SplineType& spline)
{
  static_assert(SplineType::periodicType == PeriodicSplineType::Periodic, "SplineType must be periodic");
  return computeRightDomain(spline) - computeLeftDomain(spline);
}

template <PeriodicSplineType periodicType>
struct CorrectTRelativeSize
{
};
template <>
struct CorrectTRelativeSize<PeriodicSplineType::Periodic>
{
  static const ptrdiff_t value = 1;
};
template <>
struct CorrectTRelativeSize<PeriodicSplineType::Nonperiodic>
{
  static const ptrdiff_t value = 0;
};

// this implementation can benefit from a scalar implementation of Hermite Spline (derivative) evaluation

// this is inspired by http://eigenjoy.com/2011/09/09/binary-search-revisited/
// Return 1 << log_2(list_size-1), or 0 if list_size == 1.
// This sets the initial value of b in fbsearch().
// Meanwhile, the actual solution is based on
// http://stackoverflow.com/questions/466204/rounding-up-to-nearest-power-of-2
inline ptrdiff_t init_bit(ptrdiff_t list_size)
{
  if(list_size == 1)
    return 0;
  list_size -= 1;
  ptrdiff_t power = 1;
  while(list_size >>= 1)
    power <<= 1;
  return power;
}

// this is inspired by http://eigenjoy.com/2011/09/09/binary-search-revisited/
// which doesn't need to keep track of much "extra" information and so it scales
// nicely to a bulk binary search
inline void binarySearch(const Eigen::Ref<const Eigen::VectorXd>& haystack, const double needle, ptrdiff_t& index)
{
  assert(isIncreasing(haystack));
  if(haystack.size() == 0)
  {
    index = -1;
    return;
  }
  index = 0;

  for(ptrdiff_t b = init_bit(haystack.size()); b; b >>= 1)
  {
    ptrdiff_t j = index | b;
    if(haystack.size() <= j)
      continue;
    if(haystack[j] <= needle)
      index = j;
  }
  index = (haystack[0] <= needle) ? index : -1;
}

// this is inspired by http://eigenjoy.com/2011/09/09/binary-search-revisited/
// which doesn't need to keep track of much "extra" information and so it scales
// nicely to a bulk binary search
inline void binarySearch(const Eigen::Ref<const Eigen::VectorXd>& haystack,
                         const Eigen::Ref<const Eigen::VectorXd>& needles,
                         Eigen::Ref<Eigen::Matrix<ptrdiff_t, Eigen::Dynamic, 1>> indices)
{
  assert(isIncreasing(haystack));
  assert(needles.size() == indices.size());
  if(haystack.size() == 0)
  {
    indices.setConstant(-1);
    return;
  }
  indices.setZero();

  for(ptrdiff_t b = init_bit(haystack.size()); b; b >>= 1)
  {
    for (ptrdiff_t element = 0; element < needles.size(); element++)
    {
      ptrdiff_t j = indices[element] | b;
      if (haystack.size() <= j) continue;
      if (haystack[j] <= needles[element]) indices[element] = j;
    }
  }
  indices = (haystack[0] <= needles.array()).select(indices, -1);
}

template <ptrdiff_t derivativeOrder>
double compute_h00(const double /*alpha*/)
{
  throw std::logic_error("Never call the unspecialized compute_h00");
}
template <>
inline double compute_h00<0>(const double alpha)
{
  return (1.0 + 2.0 * alpha) * (1.0 - alpha) * (1.0 - alpha);
}
template <>
inline double compute_h00<1>(const double alpha)
{
  // compute the first derivative (this is continuous across knot boundaries)
  return 6.0 * alpha * (alpha - 1.0); // = 6*alpha^2 - 6*alpha
}
template <>
inline double compute_h00<2>(const double alpha)
{
  // compute the second derivative (this is discontinuous across knot boundaries, but at least always finite)
  return 12.0 * alpha - 6.0;
}

template <ptrdiff_t derivativeOrder>
double compute_h10(const double /*alpha*/)
{
  throw std::logic_error("Never call the unspecialized compute_h10");
}
template <>
inline double compute_h10<0>(const double alpha)
{
  return alpha * (1.0 - alpha) * (1.0 - alpha);
}
template <>
inline double compute_h10<1>(const double alpha)
{
  // compute the first derivative (this is continuous across knot boundaries)
  return (3.0 * alpha - 1.0) * (alpha - 1.0); // = 3*alpha^2 - 4*alpha + 1
}
template <>
inline double compute_h10<2>(const double alpha)
{
  // compute the second derivative (this is discontinuous across knot boundaries, but at least always finite)
  return 6.0 * alpha - 4.0;
}

template <ptrdiff_t derivativeOrder>
double compute_h01(const double /*alpha*/)
{
  throw std::logic_error("Never call the unspecialized compute_h01");
}
template <>
inline double compute_h01<0>(const double alpha)
{
  return alpha * alpha * (3.0 - 2.0 * alpha);
}
template <>
inline double compute_h01<1>(const double alpha)
{
  // compute the first derivative (this is continuous across knot boundaries)
  return 6.0 * alpha * (1.0 - alpha); // = -6*alpha^2 + 6*alpha
}
template <>
inline double compute_h01<2>(const double alpha)
{
  // compute the second derivative (this is discontinuous across knot boundaries, but at least always finite)
  return -12.0 * alpha + 6.0;
}

template <ptrdiff_t derivativeOrder>
double compute_h11(const double /*alpha*/)
{
  throw std::logic_error("Never call the unspecialized compute_h11");
}
template <>
inline double compute_h11<0>(const double alpha)
{
  return alpha * alpha * (alpha - 1.0);
}
template <>
inline double compute_h11<1>(const double alpha)
{
  // compute the first derivative (this is continuous across knot boundaries)
  return alpha * (3.0 * alpha - 2.0); // = 3*alpha^2 - 2*alpha
}
template <>
inline double compute_h11<2>(const double alpha)
{
  // compute the second derivative (this is discontinuous across knot boundaries, but at least always finite)
  return 6.0 * alpha - 2.0;
}

template <ptrdiff_t derivativeOrder>
void compute_h00(const Eigen::Ref<const Eigen::ArrayXd>& /*alpha*/, Eigen::Ref<Eigen::ArrayXd> /*output*/)
{
  throw std::logic_error("Never call the unspecialized compute_h00");
}
template <>
inline void compute_h00<0>(const Eigen::Ref<const Eigen::ArrayXd>& alpha, Eigen::Ref<Eigen::ArrayXd> output)
{
  output = (1.0 + 2.0 * alpha) * (1.0 - alpha).abs2();
}
template <>
inline void compute_h00<1>(const Eigen::Ref<const Eigen::ArrayXd>& alpha, Eigen::Ref<Eigen::ArrayXd> output)
{
  // compute the first derivative (this is continuous across knot boundaries)
  output = 6.0 * alpha * (alpha - 1.0); // = 6*alpha^2 - 6*alpha
}
template <>
inline void compute_h00<2>(const Eigen::Ref<const Eigen::ArrayXd>& alpha, Eigen::Ref<Eigen::ArrayXd> output)
{
  // compute the second derivative (this is discontinuous across knot boundaries, but at least always finite)
  output = 12.0 * alpha - 6.0;
}

template <ptrdiff_t derivativeOrder>
void compute_h10(const Eigen::Ref<const Eigen::ArrayXd>& /*alpha*/, Eigen::Ref<Eigen::ArrayXd> /*output*/)
{
  throw std::logic_error("Never call the unspecialized compute_h10");
}
template <>
inline void compute_h10<0>(const Eigen::Ref<const Eigen::ArrayXd>& alpha, Eigen::Ref<Eigen::ArrayXd> output)
{
  output = alpha * (1.0 - alpha).abs2();
}
template <>
inline void compute_h10<1>(const Eigen::Ref<const Eigen::ArrayXd>& alpha, Eigen::Ref<Eigen::ArrayXd> output)
{
  // compute the first derivative (this is continuous across knot boundaries)
  output = (3.0 * alpha - 1.0) * (alpha - 1.0); // = 3*alpha^2 - 4*alpha + 1
}
template <>
inline void compute_h10<2>(const Eigen::Ref<const Eigen::ArrayXd>& alpha, Eigen::Ref<Eigen::ArrayXd> output)
{
  // compute the second derivative (this is discontinuous across knot boundaries, but at least always finite)
  output = 6.0 * alpha - 4.0;
}

template <ptrdiff_t derivativeOrder>
void compute_h01(const Eigen::Ref<const Eigen::ArrayXd>& /*alpha*/, Eigen::Ref<Eigen::ArrayXd> /*output*/)
{
  throw std::logic_error("Never call the unspecialized compute_h01");
}
template <>
inline void compute_h01<0>(const Eigen::Ref<const Eigen::ArrayXd>& alpha, Eigen::Ref<Eigen::ArrayXd> output)
{
  output = alpha.abs2() * (3.0 - 2.0 * alpha);
}
template <>
inline void compute_h01<1>(const Eigen::Ref<const Eigen::ArrayXd>& alpha, Eigen::Ref<Eigen::ArrayXd> output)
{
  // compute the first derivative (this is continuous across knot boundaries)
  output = 6.0 * alpha * (1.0 - alpha); // = -6*alpha^2 + 6*alpha
}
template <>
inline void compute_h01<2>(const Eigen::Ref<const Eigen::ArrayXd>& alpha, Eigen::Ref<Eigen::ArrayXd> output)
{
  // compute the second derivative (this is discontinuous across knot boundaries, but at least always finite)
  output = -12.0 * alpha + 6.0;
}

template <ptrdiff_t derivativeOrder>
void compute_h11(const Eigen::Ref<const Eigen::ArrayXd>& /*alpha*/, Eigen::Ref<Eigen::ArrayXd> /*output*/)
{
  throw std::logic_error("Never call the unspecialized compute_h11");
}
template <>
inline void compute_h11<0>(const Eigen::Ref<const Eigen::ArrayXd>& alpha, Eigen::Ref<Eigen::ArrayXd> output)
{
  output = alpha.abs2() * (alpha - 1.0);
}
template <>
inline void compute_h11<1>(const Eigen::Ref<const Eigen::ArrayXd>& alpha, Eigen::Ref<Eigen::ArrayXd> output)
{
  // compute the first derivative (this is continuous across knot boundaries)
  output = alpha * (3.0 * alpha - 2.0); // = 3*alpha^2 - 2*alpha
}
template <>
inline void compute_h11<2>(const Eigen::Ref<const Eigen::ArrayXd>& alpha, Eigen::Ref<Eigen::ArrayXd> output)
{
  // compute the second derivative (this is discontinuous across knot boundaries, but at least always finite)
  output = 6.0 * alpha - 2.0;
}

template <ptrdiff_t derivativeOrder = 0, class SplineType>
void evaluateVector(const SplineType& spline, typename SplineType::ConstVectorRefType inputT,
                    typename SplineType::MatrixRefType output)
{
  static_assert(derivativeOrder == 0 || derivativeOrder == 1 || derivativeOrder == 2,
                "Evaluation of Hermite splines is only defined for derivative orders 0, 1, and 2");

  typedef Eigen::Matrix<ptrdiff_t, Eigen::Dynamic, 1> IndexVector;

  const double leftDomain = computeLeftDomain(spline);
  const double rightDomain = computeRightDomain(spline);
  const double period = rightDomain - leftDomain;

  // check input sanity
  assert(isIncreasing(spline.t));
  assert(inputT.size() == output.cols());
  assert(spline.points.cols() == spline.tangents.cols());
  assert(spline.t.size() > 1);
  assert(spline.t.size() == spline.points.cols() + CorrectTRelativeSize<SplineType::periodicType>::value);

  // for periodic splines, wrap to the main period
#pragma warning(push)
#pragma warning(disable:4127) // ignore conditional-expression-is-constant warning
  Eigen::ArrayXd t = inputT;
  if(SplineType::periodicType == PeriodicSplineType::Periodic)
  {
    Eigen::ArrayXd ratio = (t - leftDomain) / period;
    Eigen::ArrayXd alpha = ratio - ratio.floor();
    t = leftDomain + alpha * period;
    //auto wrapToMainPeriod = [leftDomain, period](double x) 
    //{
    //  double ratio = (x - leftDomain) / period;
    //  double alpha = ratio - std::floor(ratio);
    //  return leftDomain + alpha * period;
    //};
    //t = t.unaryExpr(wrapToMainPeriod);
  }
#pragma warning(pop)

  // find the spline sub-sections
  assert((t >= leftDomain).all());
  assert((t <= rightDomain).all());
  IndexVector indices(t.size());
  binarySearch(spline.t, t, indices);

  // verify we can access both index i and and index i+1
  assert((indices.array() >= 0).all());
  assert(((indices.array() < spline.t.size() - 1) || (t == rightDomain)).all());

  // verify the bounds we think we have
#ifndef NDEBUG
  for(ptrdiff_t i = 0; i < t.size(); i++)
  {
    assert(spline.t[indices[i]] <= t[i]);
    if(t[i] < rightDomain)
    {
      assert(spline.t[indices[i] + 1] > t[i]);
    }
    else
    {
      assert(spline.t[indices[i]] == t[i]);
    }
  }
#endif

  // at the end of the domain, shift left
  indices = (t == rightDomain).select(indices.array() - 1, indices);

  // assemble the spline coefficients
  // see https://en.wikipedia.org/wiki/Cubic_Hermite_spline for much of the symbol meanings
  Eigen::ArrayXd valueScale(t.size());
  Eigen::ArrayXd derivativeScale(t.size());
  Eigen::ArrayXd alpha(t.size()); // affine transformation of each interval to [0, 1]
  Eigen::Array<double, SplineType::rangeDimension, Eigen::Dynamic> pk(SplineType::rangeDimension, t.size());
  Eigen::Array<double, SplineType::rangeDimension, Eigen::Dynamic> mk(SplineType::rangeDimension, t.size());
  Eigen::Array<double, SplineType::rangeDimension, Eigen::Dynamic> pk1(SplineType::rangeDimension, t.size());
  Eigen::Array<double, SplineType::rangeDimension, Eigen::Dynamic> mk1(SplineType::rangeDimension, t.size());
  for(ptrdiff_t i = 0; i < t.size(); i++)
  {
    ptrdiff_t k = indices[i];
    double dt = spline.t[k + 1] - spline.t[k];
    alpha[i] = (t[i] - spline.t[k]) / dt; // affine transformation to [0, 1) interval
    valueScale[i] = pow(dt, -static_cast<double>(derivativeOrder));
    derivativeScale[i] = dt * valueScale[i];
    pk.col(i) = spline.points.col(k);
    mk.col(i) = spline.tangents.col(k);
    pk1.col(i) = spline.points.col((k + 1) % spline.points.cols());
    mk1.col(i) = spline.tangents.col((k + 1) % spline.points.cols());
  }
  assert((derivativeScale > 0.0).all());
  assert((alpha >= 0.0).all());
  assert((alpha <= 1.0).all());
  assert(((alpha < 1.0) || (t == rightDomain)).all());

  // compute the spline basis functions
  Eigen::ArrayXd h00(t.size());
  Eigen::ArrayXd h10(t.size());
  Eigen::ArrayXd h01(t.size());
  Eigen::ArrayXd h11(t.size());
  compute_h00<derivativeOrder>(alpha, h00);
  compute_h10<derivativeOrder>(alpha, h10);
  compute_h01<derivativeOrder>(alpha, h01);
  compute_h11<derivativeOrder>(alpha, h11);

  // compute the results
  output = pk.rowwise() * (h00 * valueScale).transpose() + mk.rowwise() * (h10 * derivativeScale).transpose() +
           pk1.rowwise() * (h01 * valueScale).transpose() + mk1.rowwise() * (h11 * derivativeScale).transpose();
}

template <class SplineType>
void evaluateScalarAndDerivatives(const SplineType& spline, double inputT,
                                  Eigen::Ref<Eigen::Matrix<double, SplineType::rangeDimension, 1>> output,
                                  Eigen::Ref<Eigen::Matrix<double, SplineType::rangeDimension, 1>> outDerivative,
                                  Eigen::Ref<Eigen::Matrix<double, SplineType::rangeDimension, 1>> outSecondDerivative)
{
  typedef Eigen::Matrix<double, SplineType::rangeDimension, 1> Vector;

  const double leftDomain = computeLeftDomain(spline);
  const double rightDomain = computeRightDomain(spline);
  const double period = rightDomain - leftDomain;

  // check input sanity
  assert(isIncreasing(spline.t));
  assert(spline.points.cols() == spline.tangents.cols());
  assert(spline.t.size() > 1);
  assert(spline.t.size() == spline.points.cols() + CorrectTRelativeSize<SplineType::periodicType>::value);

// for periodic splines, wrap to the main period
#pragma warning(push)
#pragma warning(disable : 4127) // ignore conditional-expression-is-constant warning
  double t = inputT;
  if(SplineType::periodicType == PeriodicSplineType::Periodic)
  {
    double ratio = (t - leftDomain) / period;
    double alpha = ratio - std::floor(ratio);
    t = leftDomain + alpha * period;
  }
#pragma warning(pop)

  // find the spline sub-sections
  assert(t >= leftDomain);
  assert(t <= rightDomain);
  ptrdiff_t index = -1;
  binarySearch(spline.t, t, index);

  // verify we can access both index i and and index i+1
  assert(index >= 0);
  assert((index < spline.t.size() - 1) || (t == rightDomain));

  // verify the bounds we think we have
  assert(spline.t[index] <= t);
  if(t < rightDomain)
  {
    assert(spline.t[index + 1] > t);
  }
  else
  {
    assert(spline.t[index] == t);
  }

  // at the end of the domain, shift left
  if(t == rightDomain)
  {
    index--;
  }

  // assemble the spline coefficients
  // see https://en.wikipedia.org/wiki/Cubic_Hermite_spline for much of the symbol meanings
  const double dt = spline.t[index + 1] - spline.t[index];
  const double alpha = (t - spline.t[index]) / dt; // affine transformation to [0, 1) interval
  const Vector pk = spline.points.col(index);
  const Vector mk = spline.tangents.col(index);
  const Vector pk1 = spline.points.col((index + 1) % spline.points.cols());
  const Vector mk1 = spline.tangents.col((index + 1) % spline.points.cols());
  assert(dt > 0.0);
  assert(alpha >= 0.0);
  assert(alpha <= 1.0);
  assert((alpha < 1.0) || (t == rightDomain));

  // compute the spline basis functions
  const double h00_0 = compute_h00<0>(alpha);
  const double h10_0 = compute_h10<0>(alpha);
  const double h01_0 = compute_h01<0>(alpha);
  const double h11_0 = compute_h11<0>(alpha);
  const double h00_1 = compute_h00<1>(alpha);
  const double h10_1 = compute_h10<1>(alpha);
  const double h01_1 = compute_h01<1>(alpha);
  const double h11_1 = compute_h11<1>(alpha);
  const double h00_2 = compute_h00<2>(alpha);
  const double h10_2 = compute_h10<2>(alpha);
  const double h01_2 = compute_h01<2>(alpha);
  const double h11_2 = compute_h11<2>(alpha);

  // compute the results
  const double inv_dt = 1.0 / dt;
  const double inv_dt2 = inv_dt * inv_dt;
  output = pk * h00_0 + mk * h10_0 * dt + pk1 * h01_0 + mk1 * h11_0 * dt;
  outDerivative = pk * h00_1 * inv_dt + mk * h10_1 + pk1 * h01_1 * inv_dt + mk1 * h11_1;
  outSecondDerivative = pk * h00_2 * inv_dt2 + mk * h10_2 * inv_dt + pk1 * h01_2 * inv_dt2 + mk1 * h11_2 * inv_dt;
}

template <class SplineType>
void evaluateScalarAndDerivative(const SplineType& spline, double inputT,
                                 Eigen::Ref<Eigen::Matrix<double, SplineType::rangeDimension, 1>> output,
                                 Eigen::Ref<Eigen::Matrix<double, SplineType::rangeDimension, 1>> outDerivative)
{
  typedef Eigen::Matrix<double, SplineType::rangeDimension, 1> Vector;

  const double leftDomain = computeLeftDomain(spline);
  const double rightDomain = computeRightDomain(spline);
  const double period = rightDomain - leftDomain;

  // check input sanity
  assert(isIncreasing(spline.t));
  assert(spline.points.cols() == spline.tangents.cols());
  assert(spline.t.size() > 1);
  assert(spline.t.size() == spline.points.cols() + CorrectTRelativeSize<SplineType::periodicType>::value);

// for periodic splines, wrap to the main period
#pragma warning(push)
#pragma warning(disable : 4127) // ignore conditional-expression-is-constant warning
  double t = inputT;
  if(SplineType::periodicType == PeriodicSplineType::Periodic)
  {
    double ratio = (t - leftDomain) / period;
    double alpha = ratio - std::floor(ratio);
    t = leftDomain + alpha * period;
  }
#pragma warning(pop)

  // find the spline sub-sections
  assert(t >= leftDomain);
  assert(t <= rightDomain);
  ptrdiff_t index = -1;
  binarySearch(spline.t, t, index);

  // verify we can access both index i and and index i+1
  assert(index >= 0);
  assert((index < spline.t.size() - 1) || (t == rightDomain));

  // verify the bounds we think we have
  assert(spline.t[index] <= t);
  if(t < rightDomain)
  {
    assert(spline.t[index + 1] > t);
  }
  else
  {
    assert(spline.t[index] == t);
  }

  // at the end of the domain, shift left
  if(t == rightDomain)
  {
    index--;
  }

  // assemble the spline coefficients
  // see https://en.wikipedia.org/wiki/Cubic_Hermite_spline for much of the symbol meanings
  const double dt = spline.t[index + 1] - spline.t[index];
  const double alpha = (t - spline.t[index]) / dt; // affine transformation to [0, 1) interval
  const Vector pk = spline.points.col(index);
  const Vector mk = spline.tangents.col(index);
  const Vector pk1 = spline.points.col((index + 1) % spline.points.cols());
  const Vector mk1 = spline.tangents.col((index + 1) % spline.points.cols());
  assert(dt > 0.0);
  assert(alpha >= 0.0);
  assert(alpha <= 1.0);
  assert((alpha < 1.0) || (t == rightDomain));

  // compute the spline basis functions
  const double h00_0 = compute_h00<0>(alpha);
  const double h10_0 = compute_h10<0>(alpha);
  const double h01_0 = compute_h01<0>(alpha);
  const double h11_0 = compute_h11<0>(alpha);
  const double h00_1 = compute_h00<1>(alpha);
  const double h10_1 = compute_h10<1>(alpha);
  const double h01_1 = compute_h01<1>(alpha);
  const double h11_1 = compute_h11<1>(alpha);

  // compute the results
  const double inv_dt = 1.0 / dt;
  output = pk * h00_0 + mk * h10_0 * dt + pk1 * h01_0 + mk1 * h11_0 * dt;
  outDerivative = pk * h00_1 * inv_dt + mk * h10_1 + pk1 * h01_1 * inv_dt + mk1 * h11_1;
}

template <class SplineType>
void evaluateScalar(const SplineType& spline, double inputT,
                    Eigen::Ref<Eigen::Matrix<double, SplineType::rangeDimension, 1>> output)
{
  typedef Eigen::Matrix<double, SplineType::rangeDimension, 1> Vector;

  const double leftDomain = computeLeftDomain(spline);
  const double rightDomain = computeRightDomain(spline);
  const double period = rightDomain - leftDomain;

  // check input sanity
  assert(isIncreasing(spline.t));
  assert(spline.points.cols() == spline.tangents.cols());
  assert(spline.t.size() > 1);
  assert(spline.t.size() == spline.points.cols() + CorrectTRelativeSize<SplineType::periodicType>::value);

// for periodic splines, wrap to the main period
#pragma warning(push)
#pragma warning(disable : 4127) // ignore conditional-expression-is-constant warning
  double t = inputT;
  if(SplineType::periodicType == PeriodicSplineType::Periodic)
  {
    double ratio = (t - leftDomain) / period;
    double alpha = ratio - std::floor(ratio);
    t = leftDomain + alpha * period;
  }
#pragma warning(pop)

  // find the spline sub-sections
  assert(t >= leftDomain);
  assert(t <= rightDomain);
  ptrdiff_t index = -1;
  binarySearch(spline.t, t, index);

  // verify we can access both index i and and index i+1
  assert(index >= 0);
  assert((index < spline.t.size() - 1) || (t == rightDomain));

  // verify the bounds we think we have
  assert(spline.t[index] <= t);
  if(t < rightDomain)
  {
    assert(spline.t[index + 1] > t);
  }
  else
  {
    assert(spline.t[index] == t);
  }

  // at the end of the domain, shift left
  if(t == rightDomain)
  {
    index--;
  }

  // assemble the spline coefficients
  // see https://en.wikipedia.org/wiki/Cubic_Hermite_spline for much of the symbol meanings
  const double dt = spline.t[index + 1] - spline.t[index];
  const double alpha = (t - spline.t[index]) / dt; // affine transformation to [0, 1) interval
  const Vector pk = spline.points.col(index);
  const Vector mk = spline.tangents.col(index);
  const Vector pk1 = spline.points.col((index + 1) % spline.points.cols());
  const Vector mk1 = spline.tangents.col((index + 1) % spline.points.cols());
  assert(dt > 0.0);
  assert(alpha >= 0.0);
  assert(alpha <= 1.0);
  assert((alpha < 1.0) || (t == rightDomain));

  // compute the spline basis functions
  const double h00_0 = compute_h00<0>(alpha);
  const double h10_0 = compute_h10<0>(alpha);
  const double h01_0 = compute_h01<0>(alpha);
  const double h11_0 = compute_h11<0>(alpha);

  // compute the results
  output = pk * h00_0 + mk * h10_0 * dt + pk1 * h01_0 + mk1 * h11_0 * dt;
}

template <bool templateIsPeriodic, ptrdiff_t dimension>
class HermiteCurve : public Curve<dimension>
{
public:
  static const PeriodicSplineType periodicType =
      templateIsPeriodic ? PeriodicSplineType::Periodic : PeriodicSplineType::Nonperiodic;
  typedef HermiteSpline<dimension, periodicType> SplineType;
  typedef Eigen::Matrix<double, 1 + 2 * dimension, Eigen::Dynamic> TPointsTangents;

private:
  SplineType spline;

public:
  HermiteCurve(const Eigen::Ref<const TPointsTangents>& tpointstangents)
  {
    // construct the spline object
    const ptrdiff_t N = tpointstangents.cols();
    spline.t = tpointstangents.row(0).transpose();
    if(isPeriodic())
    {
      spline.points = tpointstangents.middleRows<dimension>(1).leftCols(N - 1);
      spline.tangents = tpointstangents.bottomRows<dimension>().leftCols(N - 1);
    }
    else
    {
      spline.points = tpointstangents.middleRows<dimension>(1);
      spline.tangents = tpointstangents.bottomRows<dimension>();
    }
  }

  // Inherited via Curve<dimension>
  virtual bool isPeriodic() const override final
  {
    return templateIsPeriodic;
  }
  virtual double t0() const override final
  {
    return computeLeftDomain(spline);
  }
  virtual double t1() const override final
  {
    return computeRightDomain(spline);
  }
  virtual double period() const override final
  {
    return t1() - t0();
  }
  virtual void parametricBounds(double* const outBounds) const override final
  {
    outBounds[0] = computeLeftDomain(spline);
    outBounds[1] = computeRightDomain(spline);
  }
  virtual void evaluate(const double* const t, const ptrdiff_t n, double* const outPoints, double* const outDerivative,
                        double* const outSecondDerivative) const override final
  {
    if(n == 1 && outPoints && outDerivative)
    {
      typedef Eigen::Matrix<double, dimension, 1> Vector;
      if(!outSecondDerivative)
      {
        evaluateScalarAndDerivative(spline, t[0], Eigen::Map<Vector>(outPoints), Eigen::Map<Vector>(outDerivative));
      }
      else
      {
        evaluateScalarAndDerivatives(spline, t[0], Eigen::Map<Vector>(outPoints), Eigen::Map<Vector>(outDerivative),
                                     Eigen::Map<Vector>(outSecondDerivative));
      }
    }
    else
    {
      typedef Eigen::Matrix<double, dimension, Eigen::Dynamic> Matrix;
      if(outPoints)
      {
        evaluateVector(spline, Eigen::Map<const Eigen::VectorXd>(t, n), Eigen::Map<Matrix>(outPoints, dimension, n));
      }
      if(outDerivative)
      {
        evaluateVector<1>(spline, Eigen::Map<const Eigen::VectorXd>(t, n),
                          Eigen::Map<Matrix>(outDerivative, dimension, n));
      }
      if(outSecondDerivative)
      {
        evaluateVector<2>(spline, Eigen::Map<const Eigen::VectorXd>(t, n),
                          Eigen::Map<Matrix>(outSecondDerivative, dimension, n));
      }
    }
  }
};

template <class Derived>
void cumulativeSum(const Eigen::DenseBase<Derived>& input, Eigen::Ref<Eigen::VectorXd> output)
{
  assert(input.size() == output.size());
  if(input.size() == 0)
    return;
  output[0] = input[0];
  for(ptrdiff_t i = 1; i < input.size(); i++)
  {
    output[i] = input[i] + output[i - 1];
  }
}

inline Eigen::Matrix2d makeRotate90()
{
  Eigen::Matrix2d R;
  R(0, 0) = 0.0;
  R(1, 0) = 1.0;
  R(0, 1) = -1.0;
  R(1, 1) = 0.0;
  return R;
}

inline Eigen::Vector2d makeTangent(const Eigen::Ref<const Eigen::Vector2d>& normal,
                                   const Eigen::Ref<const Eigen::Vector2d>& approximateTangent)
{
  Eigen::Vector2d result = makeRotate90() * normal;
  if(result.dot(approximateTangent) < 0.0)
  {
    result = -result;
  }
  return result;
}

template <PeriodicSplineType templatePeriodicType>
HermiteSpline<2, templatePeriodicType> initialSpline(const Eigen::Ref<const Eigen::Array2Xd>& points,
                                                     const Eigen::Ref<const Eigen::Array2Xd>& normals)
{
  // require at least a few points in the spline
  assert(points.cols() > 2);

  HermiteSpline<2, templatePeriodicType> spline;
  spline.points = points;
  assert(spline.points.allFinite());

#pragma warning(push)
#pragma warning(disable : 4127) // ignore conditional-expression-is-constant warning

  // initialize the deltaTs
  Eigen::VectorXd deltaTs(points.cols() + CorrectTRelativeSize<templatePeriodicType>::value - 1);
  deltaTs.head(points.cols() - 1) =
      (points.leftCols(points.cols() - 1) - points.rightCols(points.cols() - 1)).matrix().colwise().norm().transpose();
  if(templatePeriodicType == PeriodicSplineType::Periodic)
  {
    deltaTs[points.cols() - 1] = (points.col(0) - points.col(points.cols() - 1)).matrix().norm();
  }
  assert(deltaTs.allFinite());

  // initialize spline.t
  spline.t.resize(deltaTs.size() + 1);
  spline.t[0] = 0.0;
  cumulativeSum(deltaTs, spline.t.tail(deltaTs.size()));
  assert(spline.t.allFinite());

  // initialize spline.tangents
  spline.tangents.resize(2, points.cols());
  if(templatePeriodicType == PeriodicSplineType::Periodic)
  {
    for(ptrdiff_t i = 0; i < points.cols(); i++)
    {
      // use a periodic central difference for the approximated tangent
      Eigen::Vector2d delta = 0.5 * (spline.points.col((i + 1) % points.cols()) -
                                     spline.points.col((i - 1 + points.cols()) % points.cols()));
      spline.tangents.col(i) = makeTangent(normals.col(i), delta);
    }
  }
  else
  {
    {
      // use a forward difference for the first point
      Eigen::Vector2d delta = spline.points.col(1) - spline.points.col(0);
      spline.tangents.col(0) = makeTangent(normals.col(0), delta);
    }
    for(ptrdiff_t i = 1; i < points.cols() - 1; i++)
    {
      // use a central difference for the middle points
      Eigen::Vector2d delta = 0.5 * (spline.points.col(i + 1) - spline.points.col(i - 1));
      spline.tangents.col(i) = makeTangent(normals.col(i), delta);
    }
    {
      // use a backward difference for the last point
      Eigen::Vector2d delta = spline.points.col(points.cols() - 1) - spline.points.col(points.cols() - 2);
      spline.tangents.col(points.cols() - 1) = makeTangent(normals.col(points.cols() - 1), delta);
    }
  }
  assert(spline.tangents.allFinite());

#pragma warning(pop)

  // all done
  return spline;
}
}
}