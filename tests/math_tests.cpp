#include <gtest/gtest.h>

#define CPL_MATH_INTERNALS

#include "../cpl_math.h"

#include <cmath>

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
