cmake_minimum_required( VERSION 3.1 )

project( spinorBECdynamics )

enable_language( CXX )
set( CMAKE_CXX_STANDARD 11 )

if( NOT CMAKE_BUILD_TYPE )
   set( CMAKE_BUILD_TYPE Release )
endif()

set( CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -Wall -g" ) 
set( CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_DEBUG} -O3" )

add_executable( main main.cpp )

