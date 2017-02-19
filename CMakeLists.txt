cmake_minimum_required(VERSION 2.8.12 FATAL_ERROR)

project( sparse-ainv )

set( CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake ${CMAKE_MODULE_PATH} )
find_package(MC64 REQUIRED)

file(GLOB HEADERS
	"./source/*.h"
)

set( SRCS source/ldl_driver.cpp )

file(GLOB GFLAGS_SRCS
	"./extern/gflags/*.h"
    "./extern/gflags/*.cc"
)

if( WIN32 )
	add_definitions( -D_CRT_SECURE_NO_DEPRECATE -D_CRT_NONSTDC_NO_DEPRECATE )
else()
	add_definitions( -DUSE_STRINGS_H )
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -std=gnu++11 -mcmodel=large -lgfortran")
	SET(CMAKE_EXE_LINKER_FLAGS  "${CMAKE_EXE_LINKER_FLAGS} -lgfortran" )
endif()
add_definitions( -DSYM_ILDL_DEBUG )

include_directories( ./source ./extern/gflags ./extern/sparsehash/include ${MC64_INCLUDE_DIRS})

add_executable( ldl_driver ${SRCS} ${GFLAGS_SRCS} ${HEADERS} )
target_link_libraries(ldl_driver ${MC64_LIBRARIES})