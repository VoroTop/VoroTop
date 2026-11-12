# Builds voro++ (vendored, extern/voro) as a static library and the drivers the Python package calls: vb_engine
# (labelings, gluings, frame transport), vb_closepacked (the same for fcc and hcp, from the table of cell topologies),
# voro_faces (Voronoi faces), and vb_meridians (shortest circuits through the shell, for the branch of each atom).
cmake_minimum_required(VERSION 3.15)
project(vorotop LANGUAGES CXX)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()
find_package(OpenMP REQUIRED)

file(GLOB VORO_SOURCES ${CMAKE_SOURCE_DIR}/extern/voro/src/*.cc)
list(REMOVE_ITEM VORO_SOURCES
     ${CMAKE_SOURCE_DIR}/extern/voro/src/cmd_line.cc ${CMAKE_SOURCE_DIR}/extern/voro/src/voro++.cc
     ${CMAKE_SOURCE_DIR}/extern/voro/src/v_base_wl_2d.cc ${CMAKE_SOURCE_DIR}/extern/voro/src/v_base_wl_3d.cc)
add_library(voro STATIC ${VORO_SOURCES})
target_include_directories(voro PUBLIC ${CMAKE_SOURCE_DIR}/extern/voro/src)
target_link_libraries(voro PUBLIC OpenMP::OpenMP_CXX)
set_target_properties(voro PROPERTIES POSITION_INDEPENDENT_CODE ON)

foreach(prog vb_engine voro_faces vb_meridians vb_closepacked)
  add_executable(${prog} src/cpp/${prog}.cc)
  target_link_libraries(${prog} PRIVATE voro OpenMP::OpenMP_CXX)
  install(TARGETS ${prog} RUNTIME DESTINATION vorotop/burgers/bin)
endforeach()
