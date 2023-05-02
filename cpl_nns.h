/**
 * @file cpl_nns.h
 *
 * @brief A library for performing nearest neighbor searches on point clouds, such as those produced by LiDAR scanners
 *        or N-body simulations.
 * */
#pragma once

#ifndef CPL_NNS_H
#define CPL_NNS_H

#if !defined(CPL_NNS_FUNC) && defined(CPL_NNS_STATIC)
#define CPL_NNS_FUNC static
#endif

/**
 * @defgroup cpl_nns_api Nearest Neighbor Search API
 *
 * @brief Provides a means of performing nearest neighbor searches on K-dimensional data.
 * */

/**
 * @brief This is the callback that is invoked when a neighboring point is found.
 *
 * @param user_data An optional pointer to pass from the API function call.
 *
 * @param node The node within the tree that is considered a neighboring point.
 *
 * @returns The flag to indicate whether or not the search should terminate.
 *          A value of zero indicates that the search should not terminate.
 *          A value of non-zero will terminate the tree traversal.
 *
 * @ingroup cpl_nns_api
 * */
typedef int (*cpl_nns_search_cb)(void* user_data, const int* node);

/**
 * @brief Builds the morton codes for each point and sorts them from smallest to greatest.
 *
 * @param num_points The number of points in the point cloud.
 *
 * @param num_dims The number of values in each point. Typical values are 2 or 3.
 *
 * @param search_radius The maximum radius in which neighbors will be searched.
 *                      This cannot change after the tree has been built.
 *
 * @param points The array of points to build the search tree for.
 *
 * @param tree The tree to write the sorted morton codes to.
 *             This must contain pre-allocated memory that can fit at least:
 *               `(num_dims + 1) * num_points`
 *             integers.
 *
 * @ingroup cpl_nns_api
 * */
CPL_NNS_FUNC void
cpl_nns_prepare(int num_points, int num_dims, float search_radius, const float* points, int* tree);

CPL_NNS_FUNC void
cpl_nns_search(int num_points,
               int num_dims,
               float search_radius,
               const float* query_point,
               const int* tree,
               void* user_data,
               cpl_nns_search_cb cb);

#endif /* CPL_NNS_H */
