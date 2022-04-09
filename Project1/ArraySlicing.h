#pragma once
#include <Eigen/Dense>

namespace Eigen
{
typedef Array<bool, Dynamic, 1> ArrayXb;
typedef Array<ptrdiff_t, Dynamic, 1> ArrayXp;
} // namespace Eigen

namespace Hexagon
{
namespace Blade
{

template <class ArrayType>
Eigen::Array<typename ArrayType::Scalar, Eigen::Dynamic, 1> sliceVector(const ArrayType& inputArray,
                                                                        const Eigen::Ref<const Eigen::ArrayXb>& slice)
{
  static_assert(ArrayType::ColsAtCompileTime == 1, "we require vectors for this function");
  assert(inputArray.size() == slice.size());
  ptrdiff_t N = slice.count();
  Eigen::Array<typename ArrayType::Scalar, Eigen::Dynamic, 1> result(N);
  ptrdiff_t k = 0;
  for(ptrdiff_t i = 0; i < slice.size(); i++)
  {
    if(slice[i])
    {
      result[k] = inputArray[i];
      k++;
    }
  }
  assert(k == N);
  return result;
}

template <class ArrayType>
Eigen::Array<typename ArrayType::Scalar, Eigen::Dynamic, 1>
sliceVectorWithIndices(const ArrayType& inputArray, const Eigen::Ref<const Eigen::ArrayXp>& slice)
{
  static_assert(ArrayType::ColsAtCompileTime == 1, "we require vectors for this function");
  assert((slice >= 0).all());
  assert((slice < inputArray.size()).all());
  ptrdiff_t N = slice.size();
  Eigen::Array<typename ArrayType::Scalar, Eigen::Dynamic, 1> result(N);
  for(ptrdiff_t i = 0; i < N; i++)
  {
    result[i] = inputArray[slice[i]];
  }
  return result;
}

inline void assignToSliceVector(Eigen::Ref<Eigen::VectorXd> outputArray, const Eigen::Ref<const Eigen::ArrayXb>& slice,
                                const Eigen::Ref<const Eigen::VectorXd>& toWrite)
{
  assert(outputArray.size() == slice.size());
  assert(slice.count() == toWrite.size());
  ptrdiff_t k = 0;
  for(ptrdiff_t i = 0; i < slice.size(); i++)
  {
    if(slice[i])
    {
      outputArray[i] = toWrite[k];
      k++;
    }
  }
  assert(k == slice.count());
}

template <class ArrayType>
Eigen::Array<typename ArrayType::Scalar, Eigen::Dynamic, Eigen::Dynamic>
sliceColumns(const ArrayType& inputArray, const Eigen::Ref<const Eigen::ArrayXb>& slice)
{
  assert(inputArray.cols() == slice.size());
  ptrdiff_t N = slice.count();
  Eigen::Array<typename ArrayType::Scalar, Eigen::Dynamic, Eigen::Dynamic> result(inputArray.rows(), N);
  ptrdiff_t k = 0;
  for(ptrdiff_t i = 0; i < slice.size(); i++)
  {
    if(slice[i])
    {
      result.col(k) = inputArray.col(i);
      k++;
    }
  }
  assert(k == N);
  return result;
}

template <class ArrayType>
Eigen::Array<typename ArrayType::Scalar, Eigen::Dynamic, Eigen::Dynamic>
sliceColumnsWithIndices(const ArrayType& inputArray, const Eigen::Ref<const Eigen::ArrayXp>& slice)
{
  assert((slice >= 0).all());
  assert((slice < inputArray.cols()).all());
  ptrdiff_t N = slice.size();
  Eigen::Array<typename ArrayType::Scalar, Eigen::Dynamic, Eigen::Dynamic> result(inputArray.rows(), N);
  for(ptrdiff_t i = 0; i < slice.size(); i++)
  {
    result.col(i) = inputArray.col(slice[i]);
  }
  return result;
}

inline void assignToSliceColumns(Eigen::Ref<Eigen::MatrixXd> outputArray, const Eigen::Ref<const Eigen::ArrayXb>& slice,
                                 const Eigen::Ref<const Eigen::MatrixXd>& toWrite)
{
  assert(outputArray.cols() == slice.size());
  assert(outputArray.rows() == toWrite.rows());
  assert(slice.count() == toWrite.cols());
  ptrdiff_t k = 0;
  for(ptrdiff_t i = 0; i < slice.size(); i++)
  {
    if(slice[i])
    {
      outputArray.col(i) = toWrite.col(k);
      k++;
    }
  }
  assert(k == slice.count());
}

template <class ArrayType>
Eigen::Array<typename ArrayType::Scalar, Eigen::Dynamic, 1> circularlyShiftLeft(const ArrayType& input,
                                                                                ptrdiff_t shiftAmount)
{
  static_assert(ArrayType::ColsAtCompileTime == 1, "we require vectors for this function");
  const ptrdiff_t N = input.size();
  const ptrdiff_t shiftBranch = shiftAmount / N;
  shiftAmount -= shiftBranch * N;
  if(shiftAmount < 0)
  {
    shiftAmount += N;
  }
  assert(shiftAmount > 0);
  assert(shiftAmount < N);
  Eigen::Array<typename ArrayType::Scalar, Eigen::Dynamic, 1> result(N);
  result.head(N - shiftAmount) = input.tail(N - shiftAmount);
  result.tail(shiftAmount) = input.head(shiftAmount);
  return result;
}
} // namespace Blade
} // namespace Hexagon