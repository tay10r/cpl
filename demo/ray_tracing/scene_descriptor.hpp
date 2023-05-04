#pragma once

#include <glm/glm.hpp>

#include <string>

struct scene_descriptor final
{
  std::string id;

  std::string name;

  std::string path;

  glm::vec3 camera_origin{ 0, 0, 0 };

  glm::vec3 camera_direction{ 0, 0, -1 };
};
