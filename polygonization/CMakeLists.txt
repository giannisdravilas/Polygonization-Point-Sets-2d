# Created by the script cgal_create_CMakeLists
# This is the CMake script for compiling a set of CGAL applications.

cmake_minimum_required(VERSION 3.1...3.15)

project( to_polygon )


# CGAL and its components
find_package( CGAL QUIET COMPONENTS  )

if ( NOT CGAL_FOUND )

  message(STATUS "This project requires the CGAL library, and will not be compiled.")
  return()  

endif()


# Boost and its components
find_package( Boost REQUIRED )

if ( NOT Boost_FOUND )

  message(STATUS "This project requires the Boost library, and will not be compiled.")

  return()  

endif()

# include for local directory

# include for local package


# Creating entries for target: to_polygon
# ############################

add_executable( to_polygon  convex_hull.cpp helper_functions.cpp incremental.cpp polygonization.cpp )

add_to_cached_list( CGAL_EXECUTABLE_TARGETS to_polygon )

# Link the executable to CGAL and third-party libraries
target_link_libraries(to_polygon PRIVATE CGAL::CGAL )

