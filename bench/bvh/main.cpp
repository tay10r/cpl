#include <benchmark/benchmark.h>

#include "../../cpl_bvh.h"

#include <iostream>

#include <random>
#include <set>

#include <cstddef>

#define PRIMITIVE_COUNT (256 * 1024)

// Discretize Centroids to Bins

static void
bench_bin_indices(benchmark::State& state)
{
  const std::size_t count = PRIMITIVE_COUNT;

  std::mt19937 rng(1234);

  // Randomly shuffle primitive indices.
  // Ensure all indices are unique.

  std::vector<int> primitive_indices(count);

  {
    std::uniform_int_distribution<int> index_dist(0, count - 1);

    std::set<int> existing;

    for (int i = 0; i < count; i++) {

      while (true) {

        const auto index = index_dist(rng);

        if (existing.find(index) != existing.end())
          continue;

        primitive_indices[i] = index;

        existing.emplace(index);

        break;
      }
    }
  }

  // Randomly create centroids

  std::vector<float> xyz(count * 3);

  std::uniform_real_distribution<float> dist(-10, 10);

  for (std::size_t i = 0; i < count * 3; i++)
    xyz[i] = dist(rng);

  // Based on distribution parameters. Probably not exact but close enough.
  const float centroid_bounds[6]{ -10, 10, -10, 10, -10, 10 };

  std::vector<int> bin_indices(count);

  const int bin_axis = 0;

  for (auto _ : state) {

    float bounds[6]{};

    const int begin = 0;
    const int end = static_cast<int>(count);

    cpl_bvh_bin_indices(xyz.data(), centroid_bounds, primitive_indices.data(), bin_axis, begin, end, &bin_indices[0]);
  }
}

BENCHMARK(bench_bin_indices);

// Centroid Global Bounds

static void
bench_centroid_bounds(benchmark::State& state)
{
  std::mt19937 rng(1234);

  const std::size_t count = PRIMITIVE_COUNT;

  std::vector<float> xyz(count * 3);

  std::uniform_real_distribution<float> dist(-10, 10);

  for (std::size_t i = 0; i < count * 3; i++)
    xyz[i] = dist(rng);

  for (auto _ : state) {

    float bounds[6]{};

    cpl_bvh_centroid_bounds(xyz.data(), static_cast<int>(count), bounds);
  }
}

BENCHMARK(bench_centroid_bounds);

// Bounding Box Global Bounds

static void
bench_bbox_bounds(benchmark::State& state)
{
  std::mt19937 rng(1234);

  const std::size_t count = PRIMITIVE_COUNT;

  std::vector<float> xyz(count * 6);

  std::uniform_real_distribution<float> dist(-10, 10);

  for (std::size_t i = 0; i < count * 6; i++)
    xyz[i] = dist(rng);

  for (auto _ : state) {

    float bounds[6]{};

    cpl_bvh_bbox_bounds(xyz.data(), static_cast<int>(count), bounds);
  }
}

BENCHMARK(bench_bbox_bounds);
