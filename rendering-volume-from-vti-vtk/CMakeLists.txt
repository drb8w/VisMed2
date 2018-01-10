cmake_minimum_required(VERSION 2.8)
PROJECT(RenderingVolumeFromVTI)

find_package(VTK REQUIRED)
include(${VTK_USE_FILE})

add_executable(RenderingVolumeFromVTI MACOSX_BUNDLE RenderingVolumeFromVTI)

target_link_libraries(RenderingVolumeFromVTI ${VTK_LIBRARIES})
