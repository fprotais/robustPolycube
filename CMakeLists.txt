cmake_minimum_required(VERSION 3.27)
project(robustpolycube)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)



add_subdirectory(ultimaille)
include_directories(ultimaille ultimaille/ext)

find_package(OpenMP)
if(OPENMP_FOUND AND NOT WIN32)
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
  set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()

if (WIN32) 
	set(BUILD_SHARED_LIBS OFF)
endif()

######## LIBS ########
file(GLOB SOURCES_RB_UTILS lib/utils/*.cpp lib/utils/*.h )
add_library(RB_libutils ${SOURCES_RB_UTILS})
target_link_libraries(RB_libutils ultimaille)
if (DEFINED DEBUG_GRAPHITE_PATH) 
    target_compile_definitions(RB_libutils PRIVATE DEBUG_GRAPHITE_PATH="${DEBUG_GRAPHITE_PATH}")
    message("${DEBUG_GRAPHITE_PATH}")
endif()

file(GLOB SOURCES_RB_FLAGGING lib/flagging/*.cpp lib/flagging/*.h )
add_library(RB_libflagging ${SOURCES_RB_FLAGGING})
target_link_libraries(RB_libflagging ultimaille RB_libutils)

file(GLOB SOURCES_RB_DEFORMATION lib/deformation/*.cpp lib/deformation/*.h lib/deformation/*.cxx)
add_library(RB_libdeformation ${SOURCES_RB_DEFORMATION})
target_compile_definitions(RB_libdeformation PRIVATE -DTETLIBRARY)
target_link_libraries(RB_libdeformation ultimaille OpenNL_psm RB_libutils)

file(GLOB SOURCES_RB_QUANTIZATION lib/quantization/*.cpp lib/quantization/*.h)
add_library(RB_libquantization ${SOURCES_RB_QUANTIZATION})
target_link_libraries(RB_libquantization ultimaille OpenNL_psm RB_libutils)

file(GLOB SOURCES_RB_POSTPROCESSING lib/post_processing/*.cpp lib/post_processing/*.h)
add_library(RB_libpostprocessing ${SOURCES_RB_POSTPROCESSING})
target_link_libraries(RB_libpostprocessing ultimaille OpenNL_psm RB_libutils)

include_directories(lib)

######## EXECUTABLES ########
add_executable(rb_generate_flagging bin/generate_flagging.cpp)
target_link_libraries(rb_generate_flagging RB_libflagging)

add_executable(rb_generate_deformation bin/flagged_deformation.cpp)
target_link_libraries(rb_generate_deformation RB_libdeformation)

add_executable(rb_generate_quantization bin/extract_mesh_from_deformation.cpp)
target_link_libraries(rb_generate_quantization RB_libquantization)

add_executable(rb_perform_postprocessing bin/perform_postprocessing.cpp)
target_link_libraries(rb_perform_postprocessing RB_libpostprocessing)

add_executable(rb_fromscratch bin/fullprocess.cpp)
target_link_libraries(rb_fromscratch RB_libflagging RB_libdeformation RB_libquantization RB_libpostprocessing)
######## - ########

if(OpenMP_CXX_FOUND)
    target_link_libraries(RB_libutils OpenMP::OpenMP_CXX)
	target_link_libraries(RB_libflagging OpenMP::OpenMP_CXX)
	target_link_libraries(RB_libdeformation OpenMP::OpenMP_CXX)
	target_link_libraries(RB_libquantization OpenMP::OpenMP_CXX)
	target_link_libraries(RB_libpostprocessing OpenMP::OpenMP_CXX)

    target_link_libraries(rb_generate_flagging OpenMP::OpenMP_CXX)
    target_link_libraries(rb_generate_deformation OpenMP::OpenMP_CXX)
    target_link_libraries(rb_generate_quantization OpenMP::OpenMP_CXX)
	target_link_libraries(rb_perform_postprocessing OpenMP::OpenMP_CXX)
    target_link_libraries(rb_fromscratch OpenMP::OpenMP_CXX)
endif()

