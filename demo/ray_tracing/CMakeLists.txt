cmake_minimum_required(VERSION 3.14.7)

project(cpl_ray_tracing_demo)

include(FetchContent)

FetchContent_Declare(cornell_box
  URL "https://casual-effects.com/g3d/data10/common/model/CornellBox/CornellBox.zip"
  URL_HASH "SHA256=f27eaf47afec74b5236ce3095ac02d4256e8ac8610bfdb5a53560c67c2e0d8f1")

FetchContent_MakeAvailable(cornell_box)

FetchContent_Declare(teapot
  URL "https://casual-effects.com/g3d/data10/common/model/teapot/teapot.zip"
  URL_HASH "SHA256=a6163f5bd662915c21f5a40d2cf6a7231afe25f09cd34679a9d82948ed891ac4")

FetchContent_MakeAvailable(teapot)

FetchContent_Declare(sponza
  URL "https://casual-effects.com/g3d/data10/common/model/crytek_sponza/sponza.zip"
  URL_HASH "SHA256=da005cbee0be2df2abc8513f3ceb61bcb6f69aac112babcd9c00169a27c2770c")

FetchContent_MakeAvailable(sponza)

FetchContent_Declare(tinyobjloader
  URL "https://github.com/tinyobjloader/tinyobjloader/archive/refs/tags/v2.0.0rc10.zip"
  URL_HASH "SHA256=07d5a575395a61f0fa80ea73290dbe5863ef3b67070d3ecda3f8305c40415cd1")

FetchContent_MakeAvailable(tinyobjloader)

FetchContent_Declare(glm
  URL "https://github.com/g-truc/glm/archive/refs/tags/0.9.9.8.zip"
  URL_HASH "SHA256=4605259c22feadf35388c027f07b345ad3aa3b12631a5a316347f7566c6f1839")

FetchContent_MakeAvailable(glm)

add_executable(main
  main.cpp
  scene_descriptor.hpp
  scene.hpp
  scene.cpp
  stb_image_write.h
  stb_image_write.c
  ../../cpl_bvh.h
  cpl_bvh.c)

target_compile_definitions(main
  PRIVATE
  "CORNELL_BOX_PATH=\"${cornell_box_SOURCE_DIR}/CornellBox-Original.obj\""
  "TEAPOT_PATH=\"${teapot_SOURCE_DIR}/teapot.obj\""
  "SPONZA_PATH=\"${sponza_SOURCE_DIR}/sponza.obj\""
  )

target_include_directories(main PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/../..")

target_link_libraries(main
  PRIVATE
    tinyobjloader
    glm)
