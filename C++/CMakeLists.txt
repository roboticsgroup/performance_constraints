cmake_minimum_required(VERSION 3.5)
project(performance_constraints)

find_package(ament_cmake REQUIRED)
find_package(arl_core2 REQUIRED)

set(performance_constraints_sources
   src/performanceConstraints.cpp
)

add_library(
  ${PROJECT_NAME}
  ${performance_constraints_sources})

target_include_directories(${PROJECT_NAME} PUBLIC
"$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>"
"$<INSTALL_INTERFACE:include>")  

ament_target_dependencies(${PROJECT_NAME} arl_core2)
ament_export_include_directories(include) #${osqp_lib_INCLUDE_DIR})
ament_export_libraries(${PROJECT_NAME})

install(
DIRECTORY include/performance_constraints
DESTINATION include
)

install(
  TARGETS performance_constraints
  EXPORT performance_constraints
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib
  RUNTIME DESTINATION bin
  INCLUDES DESTINATION include
)


ament_package()


#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

#set(VISION_DIR "")

# Build in release mode
#set(CMAKE_BUILD_TYPE Release)

# find_package(catkin REQUIRED COMPONENTS
#   roscpp
#   roslib
#   autharl_core # Remove if you don't have this 
# )

# catkin_package(INCLUDE_DIRS include)

# include_directories(INCLUDE_DIRS
#   include
#   ${catkin_INCLUDE_DIRS})

# add_library(performance-constraints
#   src/performanceConstraints.cpp)

# add_dependencies(performance-constraints
#   ${${PROJECT_NAME}_EXPORTED_TARGETS}
#   ${catkin_EXPORTED_TARGETS})

# target_link_libraries(performance-constraints
#   ${catkin_LIBRARIES})

# add_executable(demo
#   src/demo.cpp
#   src/performanceConstraints.cpp)

# add_dependencies(demo
#   ${${PROJECT_NAME}_EXPORTED_TARGETS}
#   ${catkin_EXPORTED_TARGETS})

# target_link_libraries(demo
#   ${catkin_LIBRARIES}
#   pthread
# )
