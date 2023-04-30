#include "scene.hpp"

#include <tiny_obj_loader.h>

#include <fstream>
#include <sstream>

scene::~scene()
{
  free(m_nodes);
}

bool
scene::open(const char* obj_path)
{
  tinyobj::ObjReader obj_reader;

  if (!obj_reader.ParseFromFile(obj_path))
    return false;

  const auto& attrib = obj_reader.GetAttrib();

  const auto& shapes = obj_reader.GetShapes();

  for (size_t i = 0; i < shapes.size(); i++) {

    const auto index_count = shapes[i].mesh.indices.size();

    const auto shape_tri_count = index_count / 3;

    const auto* indices_ptr = shapes[i].mesh.indices.data();

    const auto prev_tri_count = m_position_buffer.size() / 9;

    m_position_buffer.resize((prev_tri_count + shape_tri_count) * 9);

    for (size_t tri_index = 0; tri_index < shape_tri_count; tri_index++) {

      const size_t index_offset = tri_index * 3;

      const auto i0 = indices_ptr[index_offset + 0];
      const auto i1 = indices_ptr[index_offset + 1];
      const auto i2 = indices_ptr[index_offset + 2];

      const auto* v0 = attrib.vertices.data() + i0.vertex_index * 3;
      const auto* v1 = attrib.vertices.data() + i1.vertex_index * 3;
      const auto* v2 = attrib.vertices.data() + i2.vertex_index * 3;

      auto* dst = &m_position_buffer[(prev_tri_count + tri_index) * 9];

      dst[0] = v0[0];
      dst[1] = v0[1];
      dst[2] = v0[2];

      dst[3] = v1[0];
      dst[4] = v1[1];
      dst[5] = v1[2];

      dst[6] = v2[0];
      dst[7] = v2[1];
      dst[8] = v2[2];
    }

    m_primitive_count += static_cast<int>(shape_tri_count);
  }

  std::vector<float> bboxes(m_primitive_count * 6);

  for (int i = 0; i < m_primitive_count; i++) {

    const float* p = m_position_buffer.data() + i * 9;

    float* bbox = &bboxes[i * 6];

    bbox[0] = p[0];
    bbox[1] = p[0];
    bbox[2] = p[1];
    bbox[3] = p[1];
    bbox[4] = p[2];
    bbox[5] = p[2];

    for (int j = 1; j < 3; j++) {
      bbox[0] = std::min(bbox[0], p[j * 3 + 0]);
      bbox[1] = std::max(bbox[1], p[j * 3 + 0]);
      bbox[2] = std::min(bbox[2], p[j * 3 + 1]);
      bbox[3] = std::max(bbox[3], p[j * 3 + 1]);
      bbox[4] = std::min(bbox[4], p[j * 3 + 2]);
      bbox[5] = std::max(bbox[5], p[j * 3 + 2]);
    }
  }

  m_primitive_indices.resize(m_primitive_count);

  m_nodes = cpl_bvh_build(m_primitive_count, bboxes.data(), &m_primitive_indices[0], &m_node_count);

  return m_nodes != nullptr;
}
