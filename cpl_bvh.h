/**
 * @file cpl_bvh.h
 *
 * @brief A self contained header for building and traversing BVHs.
 *
 * @details
 *  The algorithm used to construct the BVH is based on Ingo Wald's 2007 paper:
 *
 *    "On fast Construction of SAH-based Bounding Volume Hierarchies"
 *
 *  The BVH can be used for accelerating the process of finding ray-triangle intersections or ray-sphere intersections.
 * */

#pragma once

/**
 * @defgroup cpl_bvh_api BVH API
 *
 * @brief The API for constructing and traversing BVHs.
 * */

#ifdef CPL_BVH_STATIC
#define CPL_BVH_FUNC static
#endif

#ifndef CPL_BVH_FUNC
#ifdef __cplusplus
#define CPL_BVH_FUNC extern "C"
#else
#define CPL_BVH_FUNC
#endif
#endif

/**
 * @brief This is the type of the callback that is used to test primitives for intersection with rays.
 *
 * @param user_data An optional pointer to store additional intersection info.
 *
 * @param ray_org The origin point of the ray.
 *
 * @param ray_dir The direction of the ray.
 *
 * @param tmin The minimum distance of the ray.
 *
 * @param tmax The maximum distance of the ray.
 *
 * @param primitive_index The index of the primitive that is being tested for intersection.
 *
 * @return The distance between the ray and the primitive.
 *         A value of infinity can be used to indicate there was no intersection.
 *
 * @ingroup cpl_bvh_api
 * */
typedef float (*cpl_bvh_intersector)(void* user_data,
                                     const float* ray_org,
                                     const float* ray_dir,
                                     const float tmin,
                                     const float tmax,
                                     int primitive_index);

/**
 * @brief This structure describes a single node of the BVH.
 *
 * @ingroup cpl_bvh_api
 * */
struct cpl_bvh_node
{
  /**
   * @brief The axis-aligned bounds of the node.
   *        Each axis has a min and max value in this array.
   *        Min and max values are side by side, so:
   *          - The X-axis takes up the first two elements.
   *          - The Y-axis takes up the third and fourth elements.
   *          - The Z-axis takes up the fifth and sixth elements.
   * */
  float bounds[6];

  /**
   * @brief This field can either indicate the location of the first
   *        primitive or the location of the first child node.
   *        The meaning depends on whether or not this node points to primitives.
   *        See the @ref cpl_bvh_node::size field on how to detect this.
   * */
  int offset;

  /**
   * @brief This field indicates the number of primitives that are direct children the node.
   *        If this value is set to zero, then the children of this node are also nodes.
   *        Since child nodes always come in pairs, a value of zero also implies there are two child nodes.
   * */
  int size;
};

/**
 * @brief Calling this function will construct a BVH from a series of centroids and primitive bounding boxes.
 *
 * @details
 *  The process of this algorithm is to recursively split the scene into two bounding boxes, which leads to the least
 *  amount of cost based on the surface area heuristic.
 *
 * @param primitive_count The number of primitives that the BVH is being constructed for.
 *                        This can mean, for example, the number of triangles or spheres that are in the scene that the
 *                        BVH is being made for.
 *
 * @param bboxes The bounding boxes of each of the primitives.
 *               Each bounding box must contain the min and max values of each axis.
 *               The min and max values for each axis must be side-by-side in memory, and the min value comes first.
 *               Just to be clear, these means the layout in memory is:
 *                 - bbox[0] -> x_min
 *                 - bbox[1] -> x_max
 *                 - bbox[2] -> y_min
 *                 - bbox[3] -> y_max
 *                 - bbox[4] -> z_min
 *                 - bbox[5] -> z_max
 *
 * @param primitive_indices This is used to indicate the order of primitives in the BVH.
 *                          It should be pre-allocated to fit an index for each primitive.
 *                          After calling this function, the primitives and their attributes
 *                          should be sorted according to this order.
 *
 * @param node_count The number of nodes that were initialized during the BVH construction process.
 *                   This can be useful for serialization or debugging, but is not needed otherwise.
 *                   If this value is not needed, the parameter can be set to a null pointer.
 *
 * @return The array of BVH nodes, where the first node is the root node.
 *         This can be used to accelerate ray-primitive intersections and nearest-neighbor searches.
 *         To release the memory allocated by this structure, just pass it to `void free(void*)`.
 *
 * @ingroup cpl_bvh_api
 * */
CPL_BVH_FUNC struct cpl_bvh_node*
cpl_bvh_build(const int primitive_count, const float* bboxes, int* primitive_indices, int* node_count);

/**
 * @brief This function can be used to test for intersection between a ray and a series of geometric primitives.
 *
 * @param nodes The BVH nodes made from the primitives, in a call to @ref cpl_bvh_build.
 *
 * @param primitive_indices The indices of each of the primitives.
 *
 * @param ray_org The origin point of the ray, as a 3D vector.
 *
 * @param ray_dir The direction of the ray, as a 3D vector.
 *
 * @param tmin The minimum distance that a ray is allowed to intersect a bounding box or primitive.
 *
 * @param tmax The maximum distance that a ray is allowed to intersect a bounding box or primitive.
 *
 * @param user_data Optional data to pass to the callback function.
 *                  This may be used to store additional information on the intersections.
 *
 * @param intersector The function used to test for intersection between a ray and a geometric primitive.
 *
 * @ingroup cpl_bvh_api
 * */
CPL_BVH_FUNC void
cpl_bvh_intersect(const struct cpl_bvh_node* nodes,
                  const int* primitive_indices,
                  const float* ray_org,
                  const float* ray_dir,
                  float tmin,
                  float tmax,
                  void* user_data,
                  cpl_bvh_intersector intersector);

#if defined(CPL_BVH_IMPLEMENTATION) || defined(CPL_BVH_BENCHMARKING)

/* All functions below this point are only exposed for profiling.
 * They are not considered part of the public API.
 */

CPL_BVH_FUNC void
cpl_bvh_itoa(const int primitive_count, int* primitive_indices);

CPL_BVH_FUNC void
cpl_bvh_centroids(const float* bounds, int count, float* xyz);

CPL_BVH_FUNC void
cpl_bvh_centroid_bounds(const float* xyz, const int* primitive_indices, int begin, int end, float* bounds);

CPL_BVH_FUNC void
cpl_bvh_bbox_bounds(const float* bboxes, const int* primitive_indices, int begin, int end, float* bounds);

CPL_BVH_FUNC void
cpl_bvh_bin_indices(const float* centroids,
                    const float* centroid_bounds,
                    const int* primitive_indices,
                    const int binning_axis,
                    int begin,
                    int end,
                    int* bin_indices);

CPL_BVH_FUNC void
cpl_bvh_bin_update(const float* bboxes,
                   const int* primitive_indices,
                   const int* bin_indices,
                   int begin,
                   int end,
                   int* bin_counters,
                   float* bin_bboxes);

CPL_BVH_FUNC void
cpl_bvh_find_split(const int* bin_counters,
                   const float* bin_bboxes,
                   int* split_index,
                   float* left_bounds,
                   float* right_bounds);

CPL_BVH_FUNC int
cpl_bvh_partition(const int pivot_bin_index, int begin, int end, int* primitive_indices, int* bin_indices);

CPL_BVH_FUNC void
cpl_bvh_build_recursive(const float* centroids,
                        const float* bboxes,
                        int begin,
                        int end,
                        int* bin_indices,
                        int* bin_counters,
                        float* bin_bboxes,
                        int* primitive_indices,
                        struct cpl_bvh_node* nodes,
                        int* node_count);

#endif /* defined(CPL_BVH_IMPLEMENTATION) || defined(CPL_BVH_BENCHMARKING) */

#ifdef CPL_BVH_IMPLEMENTATION

/**
 * @defgroup cpl_bvh_internals Internals
 *
 * @brief The documentation on the implementation of the bin SAH BVH builder.
 * */

#ifdef _MSC_VER
#define _USE_MATH_DEFINES
#endif

#include <assert.h>
#include <limits.h>
#include <math.h>
#include <stdlib.h>

#ifdef CPL_BVH_LOGGING
#include <stdio.h>
#endif

#ifdef CPL_BVH_LOGGING
#define CPL_BVH_LOG(...)                                                                                               \
  do {                                                                                                                 \
    printf("%s:%d:\n", "cpl_bvh.h", __LINE__);                                                                         \
    printf("  ");                                                                                                      \
    printf(__VA_ARGS__);                                                                                               \
    printf("\n");                                                                                                      \
  } while (0)
#else
#define CPL_BVH_LOG(...)
#endif

#ifndef CPL_BVH_BIN_COUNT
#define CPL_BVH_BIN_COUNT 16
#endif

#ifndef CPL_BVH_MAX_LEAF_SIZE
#define CPL_BVH_MAX_LEAF_SIZE 4
#endif

#ifndef CPL_BVH_STACK_SIZE
#define CPL_BVH_STACK_SIZE 32
#endif

#define CPL_BVH_MIN(a, b) ((a) < (b) ? (a) : (b))

#define CPL_BVH_MAX(a, b) ((a) > (b) ? (a) : (b))

#define CPL_BVH_ABS(x) (((x) < 0) ? -(x) : (x))

CPL_BVH_FUNC struct cpl_bvh_node*
cpl_bvh_build(const int primitive_count, const float* bboxes, int* primitive_indices, int* node_count_ptr)
{
  if (primitive_count == 0)
    return NULL;

  struct cpl_bvh_node* nodes = malloc(sizeof(struct cpl_bvh_node) * ((size_t)primitive_count) * 2u);
  if (nodes == NULL)
    return NULL;

  float* centroids = malloc(sizeof(float) * 3u * ((size_t)primitive_count));
  if (!centroids) {
    free(nodes);
    return NULL;
  }

  int* bin_indices = malloc(sizeof(int) * ((size_t)primitive_count));
  if (!bin_indices) {
    free(centroids);
    free(nodes);
    return NULL;
  }

  cpl_bvh_itoa(primitive_count, primitive_indices);

  cpl_bvh_centroids(bboxes, primitive_count, centroids);

  int bin_counters[CPL_BVH_BIN_COUNT];

  float bin_bboxes[CPL_BVH_BIN_COUNT * 6];

  int node_count = 0;

  cpl_bvh_build_recursive(centroids,
                          bboxes,
                          0,
                          primitive_count,
                          bin_indices,
                          bin_counters,
                          bin_bboxes,
                          primitive_indices,
                          nodes,
                          &node_count);

  free(bin_indices);

  free(centroids);

  if (node_count_ptr)
    *node_count_ptr = node_count;

  return nodes;
}

CPL_BVH_FUNC void
cpl_bvh_itoa(const int primitive_count, int* primitive_indices)
{
  for (int i = 0; i < primitive_count; i++)
    primitive_indices[i] = i;
}

CPL_BVH_FUNC void
cpl_bvh_centroids(const float* bboxes, const int primitive_count, float* xyz)
{
  for (int i = 0; i < primitive_count; i++) {
    const float* bbox = bboxes + i * 6;
    float* centroid = xyz + i * 3;
    centroid[0] = (bbox[0] + bbox[1]) * 0.5f;
    centroid[1] = (bbox[2] + bbox[3]) * 0.5f;
    centroid[2] = (bbox[4] + bbox[5]) * 0.5f;
  }
}

CPL_BVH_FUNC void
cpl_bvh_bbox_bounds(const float* bboxes, const int* primitive_indices, const int begin, const int end, float* bounds)
{
  bounds[0] = INFINITY;
  bounds[1] = -INFINITY;
  bounds[2] = INFINITY;
  bounds[3] = -INFINITY;
  bounds[4] = INFINITY;
  bounds[5] = -INFINITY;

  for (int i = begin; i < end; i++) {

    const int j = primitive_indices[i] * 6;

    bounds[0] = CPL_BVH_MIN(bounds[0], bboxes[j + 0]);
    bounds[1] = CPL_BVH_MAX(bounds[1], bboxes[j + 1]);
    bounds[2] = CPL_BVH_MIN(bounds[2], bboxes[j + 2]);
    bounds[3] = CPL_BVH_MAX(bounds[3], bboxes[j + 3]);
    bounds[4] = CPL_BVH_MIN(bounds[4], bboxes[j + 4]);
    bounds[5] = CPL_BVH_MAX(bounds[5], bboxes[j + 5]);
  }
}

CPL_BVH_FUNC void
cpl_bvh_centroid_bounds(const float* centroids,
                        const int* primitive_indices,
                        const int begin,
                        const int end,
                        float* bounds)
{
  bounds[0] = INFINITY;
  bounds[1] = -INFINITY;
  bounds[2] = INFINITY;
  bounds[3] = -INFINITY;
  bounds[4] = INFINITY;
  bounds[5] = -INFINITY;

  for (int i = begin; i < end; i++) {

    const float* centroid = centroids + primitive_indices[i] * 3;

    const float x = centroid[0];
    const float y = centroid[1];
    const float z = centroid[2];

    bounds[0] = CPL_BVH_MIN(bounds[0], x);
    bounds[1] = CPL_BVH_MAX(bounds[1], x);
    bounds[2] = CPL_BVH_MIN(bounds[2], y);
    bounds[3] = CPL_BVH_MAX(bounds[3], y);
    bounds[4] = CPL_BVH_MIN(bounds[4], z);
    bounds[5] = CPL_BVH_MAX(bounds[5], z);
  }
}

CPL_BVH_FUNC void
cpl_bvh_bin_indices(const float* xyz,
                    const float* centroid_bounds,
                    const int* primitive_indices,
                    const int binning_axis,
                    const int begin,
                    const int end,
                    int* bin_indices)
{
  const float* axis_bounds = centroid_bounds + binning_axis * 2;

  const float axis_min = axis_bounds[0];
  const float axis_max = axis_bounds[1];

  assert(CPL_BVH_ABS(axis_min) != INFINITY);
  assert(CPL_BVH_ABS(axis_max) != INFINITY);

  assert((axis_max - axis_min) != 0.0f);

  const float axis_scale = ((float)CPL_BVH_BIN_COUNT) / (axis_max - axis_min);

  for (int i = begin; i < end; i++) {

    const int j = primitive_indices[i];

    const float value = (xyz[(j * 3) + binning_axis] - axis_min) * axis_scale;

    assert(value >= 0.0f);

    bin_indices[i] = CPL_BVH_MIN((int)value, CPL_BVH_BIN_COUNT - 1);
  }
}

CPL_BVH_FUNC void
cpl_bvh_bin_update(const float* bboxes,
                   const int* primitive_indices,
                   const int* bin_indices,
                   int begin,
                   int end,
                   int* bin_counters,
                   float* bin_bboxes)
{
  for (int i = 0; i < CPL_BVH_BIN_COUNT; i++) {

    bin_counters[i] = 0;

    float* bbox = bin_bboxes + i * 6;

    bbox[0] = INFINITY;
    bbox[1] = -INFINITY;
    bbox[2] = INFINITY;
    bbox[3] = -INFINITY;
    bbox[4] = INFINITY;
    bbox[5] = -INFINITY;
  }

  for (int i = begin; i < end; i++) {

    const int bin_index = bin_indices[i];

    assert((bin_index >= 0) && (bin_index < CPL_BVH_BIN_COUNT));

    bin_counters[bin_index]++;

    const float* bbox = bboxes + primitive_indices[i] * 6;

    float* bin_bbox = bin_bboxes + bin_index * 6;

    bin_bbox[0] = CPL_BVH_MIN(bin_bbox[0], bbox[0]);
    bin_bbox[1] = CPL_BVH_MAX(bin_bbox[1], bbox[1]);
    bin_bbox[2] = CPL_BVH_MIN(bin_bbox[2], bbox[2]);
    bin_bbox[3] = CPL_BVH_MAX(bin_bbox[3], bbox[3]);
    bin_bbox[4] = CPL_BVH_MIN(bin_bbox[4], bbox[4]);
    bin_bbox[5] = CPL_BVH_MAX(bin_bbox[5], bbox[5]);
  }
}

CPL_BVH_FUNC void
cpl_bvh_find_split(const int* bin_counters, const float* bin_bboxes, int* split_index, float* l_bounds, float* r_bounds)
{
  float l_accum_bounds[6 * CPL_BVH_BIN_COUNT];

  float r_accum_bounds[6 * CPL_BVH_BIN_COUNT];

  int l_accum_counter[CPL_BVH_BIN_COUNT];

  int r_accum_counter[CPL_BVH_BIN_COUNT];

  l_accum_counter[0] = bin_counters[0];

  l_accum_bounds[0] = bin_bboxes[0];
  l_accum_bounds[1] = bin_bboxes[1];
  l_accum_bounds[2] = bin_bboxes[2];
  l_accum_bounds[3] = bin_bboxes[3];
  l_accum_bounds[4] = bin_bboxes[4];
  l_accum_bounds[5] = bin_bboxes[5];

  r_accum_counter[CPL_BVH_BIN_COUNT - 1] = bin_counters[CPL_BVH_BIN_COUNT - 1];

  r_accum_bounds[(CPL_BVH_BIN_COUNT - 1) * 6 + 0] = bin_bboxes[(CPL_BVH_BIN_COUNT - 1) * 6 + 0];
  r_accum_bounds[(CPL_BVH_BIN_COUNT - 1) * 6 + 1] = bin_bboxes[(CPL_BVH_BIN_COUNT - 1) * 6 + 1];
  r_accum_bounds[(CPL_BVH_BIN_COUNT - 1) * 6 + 2] = bin_bboxes[(CPL_BVH_BIN_COUNT - 1) * 6 + 2];
  r_accum_bounds[(CPL_BVH_BIN_COUNT - 1) * 6 + 3] = bin_bboxes[(CPL_BVH_BIN_COUNT - 1) * 6 + 3];
  r_accum_bounds[(CPL_BVH_BIN_COUNT - 1) * 6 + 4] = bin_bboxes[(CPL_BVH_BIN_COUNT - 1) * 6 + 4];
  r_accum_bounds[(CPL_BVH_BIN_COUNT - 1) * 6 + 5] = bin_bboxes[(CPL_BVH_BIN_COUNT - 1) * 6 + 5];

  for (int i = 1; i < CPL_BVH_BIN_COUNT - 1; i++) {

    l_accum_counter[i] = l_accum_counter[i - 1] + bin_counters[i];

    float* l_bounds = &l_accum_bounds[i * 6];

    float* l_bounds_prev = &l_accum_bounds[(i - 1) * 6];

    const float* l_bin_bbox = bin_bboxes + i * 6;

    l_bounds[0] = CPL_BVH_MIN(l_bounds_prev[0], l_bin_bbox[0]);
    l_bounds[1] = CPL_BVH_MAX(l_bounds_prev[1], l_bin_bbox[1]);
    l_bounds[2] = CPL_BVH_MIN(l_bounds_prev[2], l_bin_bbox[2]);
    l_bounds[3] = CPL_BVH_MAX(l_bounds_prev[3], l_bin_bbox[3]);
    l_bounds[4] = CPL_BVH_MIN(l_bounds_prev[4], l_bin_bbox[4]);
    l_bounds[5] = CPL_BVH_MAX(l_bounds_prev[5], l_bin_bbox[5]);

    const int j = (CPL_BVH_BIN_COUNT - 1) - i;

    r_accum_counter[j] = r_accum_counter[j + 1] + bin_counters[j];

    float* r_bounds = &r_accum_bounds[j * 6];

    float* r_bounds_prev = &r_accum_bounds[(j + 1) * 6];

    const float* r_bin_bbox = bin_bboxes + j * 6;

    r_bounds[0] = CPL_BVH_MIN(r_bounds_prev[0], r_bin_bbox[0]);
    r_bounds[1] = CPL_BVH_MAX(r_bounds_prev[1], r_bin_bbox[1]);
    r_bounds[2] = CPL_BVH_MIN(r_bounds_prev[2], r_bin_bbox[2]);
    r_bounds[3] = CPL_BVH_MAX(r_bounds_prev[3], r_bin_bbox[3]);
    r_bounds[4] = CPL_BVH_MIN(r_bounds_prev[4], r_bin_bbox[4]);
    r_bounds[5] = CPL_BVH_MAX(r_bounds_prev[5], r_bin_bbox[5]);
  }

  float best_cost = INFINITY;

  int best_split = -1;

  for (int i = 0; i < CPL_BVH_BIN_COUNT - 1; i++) {

    const float* l_bounds = &l_accum_bounds[i * 6];

    const float* r_bounds = &r_accum_bounds[(i + 1) * 6];

    const float l_scale[3] = { CPL_BVH_ABS(l_bounds[1] - l_bounds[0]),
                               CPL_BVH_ABS(l_bounds[3] - l_bounds[2]),
                               CPL_BVH_ABS(l_bounds[5] - l_bounds[4]) };

    const float r_scale[3] = { CPL_BVH_ABS(r_bounds[1] - r_bounds[0]),
                               CPL_BVH_ABS(r_bounds[3] - r_bounds[2]),
                               CPL_BVH_ABS(r_bounds[5] - r_bounds[4]) };

    const float l_cost = l_scale[0] * l_scale[1] * l_scale[2] * ((float)l_accum_counter[i]);

    const float r_cost = r_scale[0] * r_scale[1] * r_scale[2] * ((float)r_accum_counter[i + 1]);

    const float cost = l_cost + r_cost;

    if (cost < best_cost) {
      best_split = i;
      best_cost = cost;
    }
  }

  assert(best_split >= 0);

  *split_index = best_split + 1;

  for (int i = 0; i < 6; i++) {

    l_bounds[i] = l_accum_bounds[best_split * 6 + i];

    r_bounds[i] = r_accum_bounds[(best_split + 1) * 6 + i];
  }
}

/**
 * @brief This function sorts the primitive indices and bin indices so that each element with a respective bin index
 *        less than that of the pivot bin index is put into the beginning of the array.
 *
 * @param pivot_bin_index The bin index that is considered the threshold to partition the arrays with.
 *                        All elements with a bin index greater than or equal to this value are put into the second
 *                        half of the array.
 *
 * @param begin The index of the first element to sort.
 *
 * @param end The index that marks the end of the array. The element belonging to this index is not included or used.
 *
 * @param bin_indices The bin indices to sort.
 *
 * @return The index of the first element in the second half of the array.
 *
 * @ingroup cpl_bvh_internals
 * */
CPL_BVH_FUNC int
cpl_bvh_partition(const int pivot_bin_index, int begin, int end, int* primitive_indices, int* bin_indices)
{
  int i = begin;

  int split_offset = end;

  while (1) {

    while ((i < split_offset) && (bin_indices[split_offset - 1] >= pivot_bin_index))
      split_offset--;

    while ((i < split_offset) && (bin_indices[i] < pivot_bin_index))
      i++;

    if (i == split_offset)
      break;

    /* bin_indices[split_offset - 1] is part of lower group */

    /* bin_indices[i] is part of upper group */

    split_offset--;

    int tmp = bin_indices[i];
    bin_indices[i] = bin_indices[split_offset];
    bin_indices[split_offset] = tmp;

    tmp = primitive_indices[i];
    primitive_indices[i] = primitive_indices[split_offset];
    primitive_indices[split_offset] = tmp;
  }

  return split_offset;
}

CPL_BVH_FUNC void
cpl_bvh_build_recursive(const float* centroids,
                        const float* bboxes,
                        int begin,
                        int end,
                        int* bin_indices,
                        int* bin_counters,
                        float* bin_bboxes,
                        int* primitive_indices,
                        struct cpl_bvh_node* nodes,
                        int* node_count)
{
  CPL_BVH_LOG("recursing to range [%d, %d]", begin, end);

  if ((end - begin) <= CPL_BVH_MAX_LEAF_SIZE) {

    const int node_index = *node_count;

    CPL_BVH_LOG("reached max leaf size at %d, creating leaf node at index %d", end - begin, node_index);

    struct cpl_bvh_node* node = nodes + node_index;

    cpl_bvh_bbox_bounds(bboxes, primitive_indices, begin, end, node->bounds);

    node->offset = begin;

    node->size = end - begin;

    (*node_count)++;

    return;
  }

  float centroid_bounds[6];

  cpl_bvh_centroid_bounds(centroids, primitive_indices, begin, end, centroid_bounds);

  float scale[3];

  scale[0] = centroid_bounds[1] - centroid_bounds[0];
  scale[1] = centroid_bounds[3] - centroid_bounds[2];
  scale[2] = centroid_bounds[5] - centroid_bounds[4];

  int binning_axis = 0;

  if (CPL_BVH_ABS(scale[1]) > CPL_BVH_ABS(scale[0]))
    binning_axis = 1;

  if (CPL_BVH_ABS(scale[2]) > CPL_BVH_ABS(scale[binning_axis]))
    binning_axis = 2;

  CPL_BVH_LOG("binning axis '%c' = [%f, %f]",
              'x' + binning_axis,
              centroid_bounds[binning_axis * 2],
              centroid_bounds[binning_axis * 2 + 1]);

  CPL_BVH_LOG("centroid bounds = [%f, %f], [%f, %f], [%f, %f]",
              centroid_bounds[0],
              centroid_bounds[1],
              centroid_bounds[2],
              centroid_bounds[3],
              centroid_bounds[4],
              centroid_bounds[5]);

  assert(scale[binning_axis] > 0.0f);

  cpl_bvh_bin_indices(centroids, centroid_bounds, primitive_indices, binning_axis, begin, end, bin_indices);

  cpl_bvh_bin_update(bboxes, primitive_indices, bin_indices, begin, end, bin_counters, bin_bboxes);

  int split_index = -1;

  float l_bounds[6];

  float r_bounds[6];

  cpl_bvh_find_split(bin_counters, bin_bboxes, &split_index, l_bounds, r_bounds);

  const int split_offset = cpl_bvh_partition(split_index, begin, end, primitive_indices, bin_indices);

  CPL_BVH_LOG("splitting at bin %d, offset %d", split_index, split_offset);

  struct cpl_bvh_node* l_node = nodes + (*node_count);

  struct cpl_bvh_node* r_node = l_node + 1;

  l_node->size = 0;
  r_node->size = 0;

  for (int i = 0; i < 6; i++) {
    l_node->bounds[i] = l_bounds[i];
    r_node->bounds[i] = r_bounds[i];
  }

  CPL_BVH_LOG("created node %d with bounds = [%f, %f], [%f, %f], [%f, %f]",
              (*node_count),
              l_node->bounds[0],
              l_node->bounds[1],
              l_node->bounds[2],
              l_node->bounds[3],
              l_node->bounds[4],
              l_node->bounds[5]);

  CPL_BVH_LOG("created node %d with bounds = [%f, %f], [%f, %f], [%f, %f]",
              (*node_count) + 1,
              r_node->bounds[0],
              r_node->bounds[1],
              r_node->bounds[2],
              r_node->bounds[3],
              r_node->bounds[4],
              r_node->bounds[5]);

  *node_count += 2;

  l_node->offset = *node_count;

  cpl_bvh_build_recursive(centroids,
                          bboxes,
                          begin,
                          split_offset,
                          bin_indices,
                          bin_counters,
                          bin_bboxes,
                          primitive_indices,
                          nodes,
                          node_count);

  r_node->offset = *node_count;

  cpl_bvh_build_recursive(
    centroids, bboxes, split_offset, end, bin_indices, bin_counters, bin_bboxes, primitive_indices, nodes, node_count);
}

CPL_BVH_FUNC void
cpl_bvh_intersect(const struct cpl_bvh_node* nodes,
                  const int* primitive_indices,
                  const float* ray_org,
                  const float* ray_dir,
                  float tmin,
                  float tmax,
                  void* user_data,
                  const cpl_bvh_intersector intersector)
{
  const float inv_dir[3] = { 1.0f / ray_dir[0], 1.0f / ray_dir[1], 1.0f / ray_dir[2] };

  int node_stack[CPL_BVH_STACK_SIZE];

  int top = 0;

  if (nodes[0].size == 0) {
    top = 2;
    node_stack[0] = 0;
    node_stack[1] = 1;
  } else {
    /* Just a single node containing primitives. */
    top = 1;
    node_stack[0] = 0;
  }

  while (top > 0) {

    const int node_index = node_stack[--top];

    const struct cpl_bvh_node* node = nodes + node_index;

    float tmp_tmin = tmin;
    float tmp_tmax = tmax;

    for (int i = 0; i < 3; i++) {

      const float t1 = (node->bounds[i * 2 + 0] - ray_org[i]) * inv_dir[i];
      const float t2 = (node->bounds[i * 2 + 1] - ray_org[i]) * inv_dir[i];

      tmp_tmin = CPL_BVH_MAX(tmp_tmin, CPL_BVH_MIN(t1, t2));
      tmp_tmax = CPL_BVH_MIN(tmp_tmax, CPL_BVH_MAX(t1, t2));
    }

    if (tmp_tmin > tmp_tmax)
      continue;

    if (node->size == 0) {

      const int child_offset = node->offset;

      node_stack[top++] = child_offset;

      if (nodes[child_offset].size == 0)
        node_stack[top++] = child_offset + 1;

      continue;
    }

    for (int i = 0; i < node->size; i++) {

      const int primitive_index = primitive_indices[node->offset + i];

      const float hit_distance = intersector(user_data, ray_org, ray_dir, tmin, tmax, primitive_index);

      if (hit_distance < tmax)
        tmax = hit_distance;
    }
  }
}

#endif /* CPL_BVH_IMPLEMENTATION */
