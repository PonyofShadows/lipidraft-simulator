cmake_minimum_required(VERSION 3.26)

project(lrsimulator VERSION 1.0)
set(CMAKE_C_STANDARD 23)
set(CMAKE_C_STANDARD_REQUIRED True)
#set(CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} -g")


include_directories(
  "${PROJECT_SOURCE_DIR}/src"
  "${PROJECT_SOURCE_DIR}/include/lrextra"
)

# Specify the minimum version you require.
find_package (mongoc-1.0 1.7 REQUIRED)


set(SRC_DR "${CMAKE_SOURCE_DIR}")
add_executable(lrsimulator run.c)


# add_compile_options(-ftrapping-math)

target_link_libraries (lrsimulator PUBLIC
  m
  mongo::mongoc_shared
)

