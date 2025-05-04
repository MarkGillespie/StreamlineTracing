cmake_minimum_required(VERSION 2.8.9)

# Maybe stop from CMAKEing in the wrong place
if (CMAKE_BINARY_DIR STREQUAL CMAKE_SOURCE_DIR)
    message(FATAL_ERROR "Source and build directories cannot be the same. Go use the /build directory.")
endif()

add_executable(
  trace-streamlines
  main.cpp
  streamlines.cpp
  utils.cpp
)
target_include_directories(trace-streamlines PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/../deps/")
# borrow args.hxx directly from polyscope
target_include_directories(trace-streamlines PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/../deps/polyscope/deps/args/")
target_link_libraries(trace-streamlines polyscope geometry-central)

option(PACKAGE_TESTS "Build the tests" ON)

if(PACKAGE_TESTS)
  enable_testing()
  add_subdirectory(
    tests
  )
endif()
