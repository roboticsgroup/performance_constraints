cmake_minimum_required(VERSION 2.8.3)
project(performance_constraints)

# find_package(Armadillo  6.500 REQUIRED )  #uncomment if you don't have autharl_core
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

#set(VISION_DIR "")

# Build in release mode
# set(CMAKE_BUILD_TYPE Release)

find_package(catkin REQUIRED COMPONENTS
  roscpp
  roslib
  autharl_core # Remove if you don't have this 
)

catkin_package(INCLUDE_DIRS include)

include_directories(INCLUDE_DIRS
  include
  ${catkin_INCLUDE_DIRS})

add_library(performance-constraints
  src/performanceConstraints.cpp
  src/Jacobian.cpp)

add_dependencies(performance-constraints
  ${${PROJECT_NAME}_EXPORTED_TARGETS}
  ${catkin_EXPORTED_TARGETS})

target_link_libraries(performance-constraints
  ${catkin_LIBRARIES})

add_executable(demo
  src/demo.cpp
  src/performanceConstraints.cpp
  src/Jacobian.cpp)

add_dependencies(demo
  ${${PROJECT_NAME}_EXPORTED_TARGETS}
  ${catkin_EXPORTED_TARGETS})

target_link_libraries(demo
  ${catkin_LIBRARIES}
  pthread
#   armadillo   #uncomment if you don't have autharl_core
)
