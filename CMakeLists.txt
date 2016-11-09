cmake_minimum_required(VERSION 2.8.12 FATAL_ERROR)

project( sparse-ainv )

set( HEADERS 	source/block_diag_matrix.h
				source/block_diag_matrix_save.h
				source/block_diag_matrix_to_string.h
				source/lil_sparse_matrix.h
				source/lilc_matrix.h
				source/lilc_matrix_declarations.h
				source/lilc_matrix_find_level_set.h
				source/lilc_matrix_find_root.h
				source/lilc_matrix_ildl.h
				source/lilc_matrix_ildl_helpers.h
				source/lilc_matrix_ildl_inplace.h
				source/lilc_matrix_load.h
				source/lilc_matrix_pivot.h
				source/lilc_matrix_save.h
				source/lilc_matrix_sym_amd.h
				source/lilc_matrix_sym_equil.h
				source/lilc_matrix_sym_perm.h
				source/lilc_matrix_to_string.h
				source/solver.h
				source/solver_minres.h
				source/solver_sqmr.h
				source/swap_struct.h )
				
set( SRCS		source/ldl_driver.cpp )

file(GLOB GFLAGS_SRCS
	"./extern/gflags/*.h"
    "./extern/gflags/*.cc"
)

add_definitions( -D_CRT_SECURE_NO_DEPRECATE -D_CRT_NONSTDC_NO_DEPRECATE )
include_directories( ./source ./extern/gflags )
add_executable( ldl_driver ${SRCS} ${GFLAGS_SRCS} ${HEADERS} )