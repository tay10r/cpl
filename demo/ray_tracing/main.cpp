#include <algorithm>
#include <iostream>
#include <string>
#include <vector>

#include <cstdlib>

#include <glm/glm.hpp>

#include <glm/gtc/type_ptr.hpp>

#include "stb_image_write.h"

#include "scene.hpp"
#include "scene_descriptor.hpp"
#include "texture.hpp"

const std::vector<scene_descriptor> g_scene_descriptors{
  { "cornell_box", "Cornell Box", CORNELL_BOX_PATH, { 0, 1, 2.1 } },
  { "teapot", "Teapot", TEAPOT_PATH, { 0, 40, 90 } },
  { "sponza", "Sponza", SPONZA_PATH, { 0, 50, 0 }, { 1, 0, 0 } },
  { "suzzane", "Suzzane", "/mnt/e/Models//McGuire Test Scenes/suzzane/untitled.obj", { 0, 0, 3 } }
};

struct intersection_info
{
  float distance;

  int cost;

  const float* position_buffer;
};

float
intersect_tri(void* info_ptr,
              const float* ray_org_ptr,
              const float* ray_dir_ptr,
              const float tmin,
              const float tmax,
              const int primitive_index)
{
  auto* info = static_cast<intersection_info*>(info_ptr);

  info->cost++;

  const float* tri_ptr = info->position_buffer + primitive_index * 9;

  const auto p0 = glm::make_vec3(tri_ptr);
  const auto p1 = glm::make_vec3(tri_ptr + 3);
  const auto p2 = glm::make_vec3(tri_ptr + 6);

  const auto ray_org = glm::make_vec3(ray_org_ptr);
  const auto ray_dir = glm::make_vec3(ray_dir_ptr);

  const auto e1 = p0 - p1;
  const auto e2 = p2 - p0;
  const auto n = cross(e1, e2);

  const auto c = p0 - ray_org;

  const auto r = cross(ray_dir, c);

  const auto inv_det = 1.0f / dot(n, ray_dir);

  const auto u = dot(r, e2) * inv_det;
  const auto v = dot(r, e1) * inv_det;
  const auto w = 1 - u - v;

  if (u >= 0 && v >= 0 && w >= 0) {
    const auto t = dot(n, c) * inv_det;
    if (t >= tmin && t <= info->distance)
      return info->distance = t;
  }

  return tmax;
}

int
main(int argc, char** argv)
{
  const std::string scene_id = (argc > 1) ? argv[1] : "cornell_box";

  auto cmp = [&](const scene_descriptor& s) -> bool { return s.id == scene_id; };

  auto it = std::find_if(g_scene_descriptors.begin(), g_scene_descriptors.end(), cmp);

  if (it == g_scene_descriptors.end()) {
    std::cerr << "Unknown scene '" << scene_id << "'." << std::endl;
    return EXIT_FAILURE;
  }

  scene s;

  if (!s.open(it->path.c_str())) {
    std::cerr << "Failed to open '" << it->path << "'." << std::endl;
    return EXIT_FAILURE;
  }

  const int w = 1280;
  const int h = 720;

  texture<glm::vec3> img(1280, 720, glm::vec3(0, 0, 0));

  const float aspect = static_cast<float>(img.width()) / img.height();

  for (int i = 0; i < img.width() * img.height(); i++) {

    const int x = i % img.width();
    const int y = i / img.width();

    const float u = (x + 0.5f) / img.width();
    const float v = (y + 0.5f) / img.height();

    constexpr auto up = glm::vec3(0, 1, 0);

    const auto right = glm::cross(up, it->camera_direction);

    const auto dir = glm::normalize(((u * 2 - 1) * aspect) * right + (1 - v * 2) * up + it->camera_direction);

    intersection_info info{ INFINITY, 0, s.get_position_buffer() };

#if 1
    cpl_bvh_intersect(s.get_nodes(),
                      s.get_primitive_indices(),
                      glm::value_ptr(it->camera_origin),
                      glm::value_ptr(dir),
                      0,
                      INFINITY,
                      &info,
                      intersect_tri);
#else
    for (int i = 0; i < s.get_primitive_count(); i++)
      intersect_tri(&info, glm::value_ptr(it->camera_origin), glm::value_ptr(dir), 0, info.distance, i);
#endif

    if (info.distance == INFINITY)
      continue;

    // img[i] = glm::vec3(info.distance, info.distance, info.distance);
    img[i] = glm::vec3(info.cost, info.cost, info.cost);
  }

  glm::vec3 min_color(0, 0, 0);

  glm::vec3 max_color(0, 0, 0);

  for (int i = 0; i < img.width() * img.height(); i++) {
    min_color = glm::min(min_color, img[i]);
    max_color = glm::max(max_color, img[i]);
  }

  const auto color_scale = 1.0f / (max_color - min_color);

  std::vector<unsigned char> rgb(img.width() * img.height() * 3);

  for (int i = 0; i < img.width() * img.height(); i++) {
    const auto c = (img[i] - min_color) * color_scale;

    rgb[(i * 3) + 0] = c.r * 255;
    rgb[(i * 3) + 1] = c.g * 255;
    rgb[(i * 3) + 2] = c.b * 255;
  }

  stbi_write_png("result.png", img.width(), img.height(), 3, rgb.data(), img.width() * 3);

  return EXIT_SUCCESS;
}
