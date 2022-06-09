cmake_minimum_required(VERSION 3.5)
project(performance_constraints)

find_package(ament_cmake REQUIRED)
find_package(arl_core2 REQUIRED)

set(performance_constraints_sources
   src/performanceConstraints.cpp
)

include_directories(
  include/performance_constraints
)

add_library(
  ${PROJECT_NAME}
  ${performance_constraints_sources})

# target_include_directories(${PROJECT_NAME} PUBLIC
# "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include/performance_constraints>"
# "$<INSTALL_INTERFACE:include>")  

ament_target_dependencies(${PROJECT_NAME} arl_core2)
ament_export_include_directories(include) #${osqp_lib_INCLUDE_DIR})
ament_export_libraries(${PROJECT_NAME})

add_executable(demo src/demo.cpp ${performance_constraints_sources})
ament_target_dependencies(demo rclcpp arl_core2 )

install(
DIRECTORY include/performance_constraints
DESTINATION include
)

install(TARGETS
  demo
  DESTINATION lib/${PROJECT_NAME})

install(
  TARGETS performance_constraints
  EXPORT performance_constraints
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib
  RUNTIME DESTINATION bin
  INCLUDES DESTINATION include
)


ament_package()
