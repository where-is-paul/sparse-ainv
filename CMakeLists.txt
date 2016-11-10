cmake_minimum_required(VERSION 2.8.12 FATAL_ERROR)

project( sparse-ainv )

file(GLOB HEADERS
	"./source/*.h"
)

set( SRCS source/ldl_driver.cpp )

file(GLOB GFLAGS_SRCS
	"./extern/gflags/*.h"
    "./extern/gflags/*.cc"
)

add_definitions( -D_CRT_SECURE_NO_DEPRECATE -D_CRT_NONSTDC_NO_DEPRECATE )
include_directories( ./source ./extern/gflags )
add_executable( ldl_driver ${SRCS} ${GFLAGS_SRCS} ${HEADERS} )