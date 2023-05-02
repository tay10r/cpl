/**
 * @file cpl_math.h
 *
 * @brief A library for various math functions.
 * */

#ifndef CPL_MATH_H
#define CPL_MATH_H

#include <stddef.h>

/* Come up with a value for CPL_MATH_FUNC if user did not specify. */

/* If requested, consider all functions to have internal linkage. */

#if !defined(CPL_MATH_FUNC) && defined(CPA_MATH_STATIC)
#define CPL_MATH_FUNC static
#endif

/* If visible to C++ compiler, use extern "C" so that linking doesn't get messed up. */

#if !defined(CPL_MATH_FUNC) && defined(__cplusplus)
#define CPL_MATH_FUNC extern "C"
#endif

/* No value guessed, so just leave it blank. */

#ifndef CPL_MATH_FUNC
#define CPL_MATH_FUNC
#endif

/**
 * @defgroup cpl_math_api Math API
 *
 * @brief An API for various math functions.
 * */

/**
 * @brief Computes the reciprocal square root of a value.
 *
 * @param x The value to compute the reciprocal square root of.
 *
 * @return The reciprocal square root of @p x.
 *
 * @ingroup cpl_math_api
 * */
CPL_MATH_FUNC float
cpl_math_rsqrt(float x);

/**
 * @brief Computes the dot product between two vectors.
 *
 * @param num_dims The number of dimensions in each vector.
 *
 * @param a The left operand.
 *
 * @param b The right operand.
 *
 * @return The dot product between @p a and @p b.
 *
 * @ingroup cpl_math_api
 * */
CPL_MATH_FUNC float
cpl_math_dot(int num_dims, const float* a, const float* b);

/**
 * @brief Normalizes a column vector.
 *
 * @param num_dims The number of dimensions in the column vector.
 *
 * @param in The column vector to normalize.
 *
 * @param out The normalized column vector.
 *
 * @ingroup cpl_math_api
 * */
CPL_MATH_FUNC
void
cpl_math_normalize(const int num_dims, const float* in, float* out);

/**
 * @brief Computes the trace of a square matrix.
 *
 * @param size The number of rows and columns of the square matrix.
 *
 * @param m The matrix to calculate the trace of.
 *
 * @return The trace of @p m.
 *
 * @ingroup cpl_math_api
 * */
CPL_MATH_FUNC float
cpl_math_trace(int size, const float* m);

/**
 * @brief Computes the determinant of a 2x2 matrix.
 *
 * @param m The matrix to get the determinant of.
 *
 * @return The determinant of @p m.
 *
 * @ingroup cpl_math_api
 * */
CPL_MATH_FUNC float
cpl_math_det2x2(const float* m);

/**
 * @brief Computes the determinant of a 3x3 matrix.
 *
 * @param m The matrix to get the determinant of.
 *
 * @return The determinant of @p m.
 *
 * @ingroup cpl_math_api
 * */
CPL_MATH_FUNC float
cpl_math_det3x3(const float* m);

/**
 * @brief Performs a QR decomposition using the Gram-Schmidt process.
 *
 * @note This method is not numerically stable.
 *
 * @param num_dims The number of dimensions in the dataset.t
 *
 * @param a The matrix being decomposed.
 *          Must be a square matrix with a row and column count of @p num_dims.
 *
 * @param q The orthogonal output matrix.
 *          Must be able to fit a square matrix with a row and column count of @p num_dims.
 *
 * @param r The upper triangular output matrix.
 *          Must be able to fit a square matrix with a row and column count of @p num_dims.
 *
 * @param tmp_u A temporary vector for the process.
 *              Must be able to fit @p num_dims single-precision floats.
 *
 * @ingroup cpl_math_internals
 * */
CPL_MATH_FUNC void
cpl_math_gram_schmidt(const int num_dims, const float* a, float* q, float* r, float* tmp_u);

#if defined(CPL_MATH_INTERNALS) || defined(CPL_MATH_IMPLEMENTATION)

#ifndef CPL_MATH_RSQRT_ITERATIONS
#define CPL_MATH_RSQRT_ITERATIONS 5
#endif

/* These are function prototypes that are not part of the public API.
 * They are exposed here only for testing.
 */

/**
 * @defgroup cpl_math_internals MATH Internals
 *
 * @brief The internal API for the MATH estimation library.
 * */

/**
 * @brief This function evaluations the projection operator defined by the Gram-Schmidt process and subtracts it from
 *        the output operand.
 *
 * @note The matrices passed to this function should be in column-major storage order.
 *       This is the same storage order used by libraries such as Eigen, Armadillo, and GLM.
 *
 * @param num_dims The number of dimensions in the MATH estimation.
 *
 * @param u The first operand in the projection operator.
 *
 * @param a The second operand in the projection operator.
 *          This operand appears only in the numerator.
 *
 * @param out The result of the projection is subtracted from this output parameter.
 *
 * @ingroup cpl_math_internals
 * */
CPL_MATH_FUNC void
cpl_math__gram_schmidt_proj_sub(const int num_dims, const float* u, const float* a, float* out);

#endif /* defined(CPL_MATH_INTERNALS) || defined(CPL_MATH_IMPLEMENTATION) */

#ifdef CPL_MATH_IMPLEMENTATION

/* Below this point is the implementation of the library.
 * Only look here if you're curious on how it works.
 */

/**
 * @brief Computes the trace of a square matrix.
 *
 * @param size The number of rows and columns of the square matrix.
 *
 * @param m The matrix to calculate the trace of.
 *
 * @return The trace of @p m.
 *
 * @ingroup cpl_math_api
 * */
CPL_MATH_FUNC float
cpl_math_trace(const int size, const float* m)
{
  float result;

  int i;

  result = 0;

  for (i = 0; i < size; i++)
    result += m[i * size + i];

  return result;
}

CPL_MATH_FUNC float
cpl_math_det2x2(const float* m)
{
  return m[0] * m[3] - m[1] * m[2];
}

CPL_MATH_FUNC float
cpl_math_det3x3(const float* m)
{
  float result;

  float tmp2x2[4];

  tmp2x2[0] = m[4];
  tmp2x2[1] = m[5];
  tmp2x2[2] = m[7];
  tmp2x2[3] = m[8];

  result = m[0] * cpl_math_det2x2(tmp2x2);

  tmp2x2[0] = m[1];
  tmp2x2[1] = m[2];
  tmp2x2[2] = m[7];
  tmp2x2[3] = m[8];

  result -= m[3] * cpl_math_det2x2(tmp2x2);

  tmp2x2[0] = m[1];
  tmp2x2[1] = m[2];
  tmp2x2[2] = m[4];
  tmp2x2[3] = m[5];

  return result + m[6] * cpl_math_det2x2(tmp2x2);
}

CPL_MATH_FUNC void
cpl_math_gram_schmidt(const int num_dims, const float* a, float* q, float* r, float* tmp_u)
{
  int row_idx;

  int column_idx;

  /* First compute Q, which consists of the column vectors resulting from the normalized projection of A. */

  /* Compute u_1, which is just the first column of A */
  for (row_idx = 0; row_idx < num_dims; row_idx++)
    tmp_u[row_idx] = a[row_idx];

  /* Compute e_1 */
  cpl_math_normalize(num_dims, tmp_u, q);

  /* Compute u_n */
  for (column_idx = 1; column_idx < num_dims; column_idx++) {

    float* u1 = tmp_u + column_idx * num_dims;

    const float* a1 = a + column_idx * num_dims;

    for (row_idx = 0; row_idx < num_dims; row_idx++)
      u1[row_idx] = a1[row_idx];

    for (row_idx = 0; row_idx < column_idx; row_idx++)
      cpl_math__gram_schmidt_proj_sub(num_dims, tmp_u + row_idx * num_dims, a1, u1);

    cpl_math_normalize(num_dims, u1, q + column_idx * num_dims);
  }

  /* compute r */
  for (column_idx = 0; column_idx < num_dims; column_idx++) {

    for (row_idx = 0; row_idx <= column_idx; row_idx++)
      r[column_idx * num_dims + row_idx] = cpl_math_dot(num_dims, a + column_idx * num_dims, q + row_idx * num_dims);

    for (row_idx = column_idx + 1; row_idx < num_dims; row_idx++)
      r[column_idx * num_dims + row_idx] = 0;
  }
}

CPL_MATH_FUNC float
cpl_math_rsqrt(const float x)
{
  int iteration;

  float result;

  result = 1;

  for (iteration = 0; iteration < CPL_MATH_RSQRT_ITERATIONS; iteration++)
    result = result * (1.5f - 0.5f * x * result * result);

  return result;
}

CPL_MATH_FUNC float
cpl_math_dot(int num_dims, const float* a, const float* b)
{
  float sum;

  int dim_idx;

  sum = 0;

  for (dim_idx = 0; dim_idx < num_dims; dim_idx++)
    sum += a[dim_idx] * b[dim_idx];

  return sum;
}

CPL_MATH_FUNC
void
cpl_math_normalize(const int num_dims, const float* in, float* out)
{
  float sum_of_squares;

  float rcp_magnitude;

  int dim_idx;

  sum_of_squares = 0;

  for (dim_idx = 0; dim_idx < num_dims; dim_idx++)
    sum_of_squares += in[dim_idx] * in[dim_idx];

  rcp_magnitude = cpl_math_rsqrt(sum_of_squares);

  for (dim_idx = 0; dim_idx < num_dims; dim_idx++)
    out[dim_idx] = in[dim_idx] * rcp_magnitude;
}

CPL_MATH_FUNC void
cpl_math__gram_schmidt_proj_sub(const int num_dims, const float* u, const float* a, float* out)
{
  float inner_dot_numerator;

  float inner_dot_denominator;

  int dim_idx;

  /* The quotient between the two inner products. */
  float quotient;

  inner_dot_numerator = 0;

  inner_dot_denominator = 0;

  for (dim_idx = 0; dim_idx < num_dims; dim_idx++) {

    inner_dot_numerator += u[dim_idx] * a[dim_idx];

    inner_dot_denominator += u[dim_idx] * u[dim_idx];
  }

  quotient = inner_dot_numerator / inner_dot_denominator;

  for (dim_idx = 0; dim_idx < num_dims; dim_idx++)
    out[dim_idx] -= u[dim_idx] * quotient;
}

#endif /* CPL_MATH_IMPLEMENTATION */

#endif /* CPL_MATH_H */
