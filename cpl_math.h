/**
 * @file cpl_math.h
 *
 * @brief A library for various math functions.
 * */

#ifndef CPL_MATH_H
#define CPL_MATH_H

/**
 * @defgroup cpl_math_api Math API
 *
 * @brief An API for various math functions.
 * */

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
 * @brief Indicates if a value is equal to zero or negative zero.
 *
 * @ingroup cpl_math_api
 * */
#define CPL_MATH_IS_ZERO(x) (((x) == 0.0f) || ((x) == -0.0f))

/**
 * @brief The value of pi in radians.
 *
 * @ingroup cpl_math_api
 * */
#define CPL_MATH_PI 3.14159265359f

/**
 * @brief The value of tau in radians.
 *
 * @ingroup cpl_math_api
 * */
#define CPL_MATH_TAU 6.28318530717f

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
 * @brief Computes the square root of a value.
 *
 * @param x The value to compute the square root of.
 *          The absolute value is considered instead of the value is negative.
 *
 * @return The square root of the given value @p x.
 *         The square root of zero and negative zero is zero.
 *
 * @ingroup cpl_math_api
 * */
CPL_MATH_FUNC float
cpl_math_sqrt(float x);

/**
 * @brief Computes the cube root of a value.
 *
 * @param x The value to compute the cube root of.
 *
 * @returns The approximated cube root of @p x.
 *
 * @ingroup cpl_math_api
 * */
CPL_MATH_FUNC float
cpl_math_cbrt(float x);

/**
 * @brief Converts degrees to radians.
 *
 * @param x The value in degrees to convert into radians.
 *
 * @return The value of @p x in radians.
 *
 * @ingroup cpl_math_api
 * */
CPL_MATH_FUNC float
cpl_math_deg2rad(float x);

/**
 * @brief Converts radians to degrees.
 *
 * @param x The value in radians to convert into degrees.
 *
 * @return The value of @p x in degrees.
 *
 * @ingroup cpl_math_api
 * */
CPL_MATH_FUNC float
cpl_math_rad2deg(float x);

/**
 * @brief Converts a set of xy coordinates to polar coordinates.
 *
 * @param xy_coords The XY coordinates to convert to polar coordinates.
 *
 * @param polar_coords The array to place the radius and angle into.
 *
 * @ingroup cpl_math_api
 * */
CPL_MATH_FUNC void
cpl_math_to_polar(const float* xy_coords, float* polar_coords);

/**
 * @brief Solves a quadratic equation for roots.
 *
 * @param in The input coefficients, which should consist of all three coefficients of the polynomial.
 *
 * @param roots The real components of each root of the equation.
 *
 * @return The number of real roots.
 *
 * @ingroup cpl_math_api
 * */
CPL_MATH_FUNC int
cpl_math_solve_quadratic(const float* in, float* roots);

/**
 * @brief Solves a cubic equation for roots.
 *
 * @param in The coefficients of the cubic equation, starting from the highest degree to the lowest degree.
 *
 * @param roots The roots of the polynomial.
 *
 * @return The number of real roots.
 *
 * @ingroup cpl_math_api
 * */
CPL_MATH_FUNC int
cpl_math_solve_cubic(const float* in, float* roots);

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

#define CPL_MATH_ABS(x) (((x) < 0) ? -(x) : (x))

#ifndef CPL_MATH_CBRT_THRESHOLD
#define CPL_MATH_CBRT_THRESHOLD 0.00001f
#endif

#ifndef CPL_MATH_SQRT_THRESHOLD
#define CPL_MATH_SQRT_THRESHOLD 0.00001f
#endif

#ifndef CPL_MATH_RSQRT_THRESHOLD
#define CPL_MATH_RSQRT_THRESHOLD 0.00001f
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

CPL_MATH_FUNC float
cpl_math_deg2rad(const float x)
{
  return x * (CPL_MATH_PI / 180.0f);
}

CPL_MATH_FUNC float
cpl_math_rad2deg(const float x)
{
  return x * (180.0f / CPL_MATH_PI);
}

CPL_MATH_FUNC void
cpl_math_to_polar(const float* xy_coords, float* polar_coords)
{
  float x;
  float y;
  float r;
  float angle;

  x = xy_coords[0];
  y = xy_coords[1];

  r = cpl_math_sqrt(x * x + y * y);

  angle = 0;

  polar_coords[0] = r;
  polar_coords[1] = angle;
}

CPL_MATH_FUNC int
cpl_math_solve_quadratic(const float* in, float* roots)
{
  /* The value of -b */
  float neg_b;

  /* The value of 1 / 2a */

  float rcp_2a;

  /* The value of the descriminant -> b^2 - 4ac */
  float discriminant;

  /* The square root of the descriminant */
  float discriminant_sqrt;

  neg_b = -in[1];

  rcp_2a = 1 / (2 * in[0]);

  discriminant = in[1] * in[1] - 4 * in[0] * in[2];

  if (CPL_MATH_IS_ZERO(discriminant)) {
    roots[0] = neg_b * rcp_2a;
    roots[1] = roots[0];
    return 1;
  }

  discriminant_sqrt = cpl_math_sqrt(discriminant);

  roots[0] = (neg_b + discriminant_sqrt) * rcp_2a;
  roots[1] = (neg_b - discriminant_sqrt) * rcp_2a;

  return (discriminant < -0.0f) ? 0 : 2;
}

CPL_MATH_FUNC int
cpl_math_solve_cubic(const float* in, float* roots)
{
  /* An implementation of the Cardano method */

  /* The coefficients to the polynomial. */
  float a;
  float b;
  float c;
  float d;

  /* Discriminant coefficients */
  float q;
  float r;
  float discriminant;
  float discriminant_sqrt_real;
  float discriminant_sqrt_imag;

  /* Temporary variables */
  float s_real;
  float s_imag;
  float t_real;
  float t_imag;
  float b_div_3a;

  a = in[0];
  b = in[1];
  c = in[2];
  d = in[3];

  b_div_3a = b / (3 * a);

  q = (3 * a * c - b * b) / (9 * a * a);

  r = (9 * a * b * c - 27 * a * a * d - 2 * b * b * b) / (54 * a * a * a);

  discriminant = q * q * q + r * r;

  discriminant_sqrt_real = cpl_math_sqrt(discriminant);
  discriminant_sqrt_imag = (discriminant < -0.0f) ? 1.0f : 0.0f;

  /* TODO : solve cube root of complex numbers. */

  s_real = cpl_math_cbrt(r + discriminant_sqrt_real);
  t_real = cpl_math_cbrt(r - discriminant_sqrt_real);

  s_imag = discriminant_sqrt_imag;
  t_imag = discriminant_sqrt_imag;

  /* z = -0.5f * (s + t) - b_div_3a; */

  /* h = 0.86602540378f * (s - t); */

  if (CPL_MATH_IS_ZERO(discriminant)) {
  }

  (void)roots;
  (void)b_div_3a;
  (void)t_imag;
  (void)s_imag;
  (void)s_real;
  (void)t_real;

  return 0;
}

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
  float current_result;
  float previous_result;
  float abs_x;

  abs_x = (x < 0) ? -x : x;

  current_result = previous_result = 1;

  while (1) {

    current_result = previous_result * (1.5f - 0.5f * abs_x * previous_result * previous_result);

    if (CPL_MATH_ABS(current_result - previous_result) < CPL_MATH_RSQRT_THRESHOLD)
      break;

    previous_result = current_result;
  }

  return current_result;
}

CPL_MATH_FUNC float
cpl_math_sqrt(const float x)
{
  float current_result;
  float previous_result;
  float abs_x;

  if (CPL_MATH_IS_ZERO(x))
    return 0.0f;

  abs_x = (x < 0) ? -x : x;

  current_result = previous_result = abs_x * 0.5f;

  while (1) {

    current_result = 0.5f * (previous_result + abs_x / previous_result);

    if (CPL_MATH_ABS(current_result - previous_result) < CPL_MATH_SQRT_THRESHOLD)
      break;

    previous_result = current_result;
  }

  return current_result;
}

CPL_MATH_FUNC float
cpl_math_cbrt(const float x)
{
  float previous_result;
  float current_result;

  current_result = previous_result = x * (1.0f / 3.0f);

  while (1) {
    current_result = previous_result - (previous_result * previous_result * previous_result - x) /
                                         (3 * previous_result * previous_result);

    if (CPL_MATH_ABS(current_result - previous_result) < CPL_MATH_CBRT_THRESHOLD)
      break;

    previous_result = current_result;
  }

  return current_result;
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
