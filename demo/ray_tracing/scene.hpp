#pragma once

#include <cpl_bvh.h>

#include <vector>

class scene final
{
public:
  scene() = default;

  scene(const scene&) = delete;

  scene(scene&&) = delete;

  scene& operator=(const scene&) = delete;

  scene& operator=(scene&&) = delete;

  ~scene();

  bool open(const char* obj_path);

  const cpl_bvh_node* get_nodes() const { return m_nodes; }

  int get_primitive_count() const { return m_primitive_count; }

  const float* get_position_buffer() const { return m_position_buffer.data(); }

  const int* get_primitive_indices() const { return m_primitive_indices.data(); }

protected:
  /// Constructs a BVH from previously loaded primitives.
  void commit_bvh();

private:
  /// X, Y, and Z position coordinates.
  std::vector<float> m_position_buffer;

  /// UV texture coordinates per vertex.
  std::vector<float> m_uv_buffer;

  /// The normal coordinates per vertex.
  std::vector<float> m_normal_buffer;

  /// The number of primitives in the scene.
  int m_primitive_count{ 0 };

  /// The BVH nodes for the scene.
  cpl_bvh_node* m_nodes{ nullptr };

  /// The number of BVH nodes there are.
  int m_node_count = 0;

  /// The primitive indices for each leaf node in the BVH.
  std::vector<int> m_primitive_indices;
};
