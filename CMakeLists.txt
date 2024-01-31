cmake_minimum_required (VERSION 3.7)
cmake_policy(VERSION 3.7...3.20)

project(TestOVM)

set(CMAKE_CXX_STANDARD 17)

file(GLOB_RECURSE SOURCES src/*.cpp src/*.c src/*.hpp src/*.h)
include_directories(src)

add_executable(TestOVM ${SOURCES})

add_subdirectory(third_party/OpenVolumeMesh)
add_subdirectory(third_party/glm)
target_link_libraries(TestOVM PRIVATE OpenVolumeMesh::OpenVolumeMesh glm::glm)
