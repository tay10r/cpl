#include <gtest/gtest.h>

#define CPL_MATH_INTERNALS

#include "../cpl_math.h"

#include <random>

#include <cmath>

TEST(cpl_math, cbrt)
{
  std::mt19937 rng(1234);

  std::uniform_real_distribution<float> dist(-100'000, 100'000);

  for (int i = 0; i < 1000; i++) {

    const float x = dist(rng);

    const float y_expected = std::cbrt(x);

    const float y_actual = cpl_math_cbrt(x);

    EXPECT_NEAR(y_expected, y_actual, 0.001f);
  }
}

TEST(cpl_math, sqrt)
{
  std::mt19937 rng(1234);

  std::uniform_real_distribution<float> dist(-1'000'000, 1'000'000);

  for (int i = 0; i < 1000; i++) {

    const float x = dist(rng);

    const float y_expected = std::sqrt(std::abs(x));

    const float y_actual = cpl_math_sqrt(x);

    EXPECT_NEAR(y_expected, y_actual, 0.001f);
  }

  EXPECT_EQ(cpl_math_sqrt(0.0f), 0.0f);

  EXPECT_EQ(cpl_math_sqrt(-0.0f), 0.0f);
}

TEST(cpl_math, solve_quadratic)
{
  { // two real solutions
    std::vector<float> in{ 1, -3, -4 };

    std::vector<float> out(2);

    EXPECT_EQ(cpl_math_solve_quadratic(in.data(), out.data()), 2);

    EXPECT_NEAR(out[0], 4.0f, 0.001f);
    EXPECT_NEAR(out[1], -1.0f, 0.001f);
  }

  { // imaginary roots
    std::vector<float> in{ 1, 2, 10 };

    std::vector<float> out(2);

    EXPECT_EQ(cpl_math_solve_quadratic(in.data(), out.data()), 0);
  }
}

TEST(cpl_math, solve_cubic)
{
  const std::vector<float> in{ 8, 4, 5, 2 };

  std::vector<float> roots(3);

  cpl_math_solve_cubic(in.data(), roots.data());

  EXPECT_NEAR(roots[0], -0.4222f, 0.001f);
  EXPECT_NEAR(roots[1], -0.0389f, 0.001f);
  EXPECT_NEAR(roots[2], -0.0389f, 0.001f);

  /*
  EXPECT_EQ(roots[0].i, 0.0f);
  EXPECT_NEAR(roots[1].i, 0.76853f, 0.001f);
  EXPECT_NEAR(roots[2].i, -0.76853f, 0.001f);
  */
}

TEST(cpl_math, trace)
{
  std::vector<float> m{
    // clang-format off
    1, 11,  6,
    0,  5, 12,
    3,  2, -5
    // clang-format on
  };

  const float result = cpl_math_trace(3, m.data());

  EXPECT_NEAR(result, 1, 0.0001);
}

TEST(cpl_math, det3x3)
{
  const std::vector<float> m{
    // clang-format off
    4, 1, 7,
    2, 5, 8,
    3, 6, 4
    // clang-format off
  };

  const float result = cpl_math_det3x3(m.data());

  EXPECT_NEAR(result, -117.0f, 0.001f);
}

TEST(cpl_math, gram_schmidt_qr_decomposition)
{
  const std::vector<float> a{ 1, 1, 0, 1, 0, 1, 0, 1, 1 };

  std::vector<float> q(9);

  std::vector<float> r(9);

  std::vector<float> tmp_u(9);

  cpl_math_gram_schmidt(3, a.data(), q.data(), r.data(), tmp_u.data());

  const std::vector<float> expected_q{
    // clang-format off
     1 / std::sqrt(2),  1 / std::sqrt(2), 0,
     1 / std::sqrt(6), -1 / std::sqrt(6), 2 / std::sqrt(6),
    -1 / std::sqrt(3),  1 / std::sqrt(3), 1 / std::sqrt(3)
    // clang-format on
  };

  const std::vector<float> expected_r{
    // clang-format off
    2 / std::sqrt(2),                0,                0,
    1 / std::sqrt(2), 3 / std::sqrt(6),                0,
    1 / std::sqrt(2), 1 / std::sqrt(6), 2 / std::sqrt(3)
    // clang-format on
  };

  for (size_t i = 0; i < 9; i++)
    EXPECT_NEAR(q[i], expected_q[i], 0.01f);

  for (size_t i = 0; i < 9; i++)
    EXPECT_NEAR(r[i], expected_r[i], 0.01f);
}
