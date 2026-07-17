cmake_minimum_required(VERSION 3.20)
project(cmpuf LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

set(CPLEX_STUDIO_DIR "/opt/ibm/ILOG/CPLEX_Studio201" CACHE PATH "IBM ILOG CPLEX Studio installation directory")
set(CPLEX_SYSTEM "x86-64_linux" CACHE STRING "CPLEX platform directory")
set(CPLEX_LIBFORMAT "static_pic" CACHE STRING "CPLEX library format")

set(CPLEX_INCLUDE_DIR "${CPLEX_STUDIO_DIR}/cplex/include")
set(CONCERT_INCLUDE_DIR "${CPLEX_STUDIO_DIR}/concert/include")
set(CPLEX_LIBRARY_DIR "${CPLEX_STUDIO_DIR}/cplex/lib/${CPLEX_SYSTEM}/${CPLEX_LIBFORMAT}")
set(CONCERT_LIBRARY_DIR "${CPLEX_STUDIO_DIR}/concert/lib/${CPLEX_SYSTEM}/${CPLEX_LIBFORMAT}")

foreach(required_path "${CPLEX_INCLUDE_DIR}" "${CONCERT_INCLUDE_DIR}" "${CPLEX_LIBRARY_DIR}" "${CONCERT_LIBRARY_DIR}")
  if(NOT EXISTS "${required_path}")
    message(FATAL_ERROR "Required CPLEX path not found: ${required_path}")
  endif()
endforeach()

add_library(cplex_runtime INTERFACE)
target_include_directories(cplex_runtime INTERFACE "${CPLEX_INCLUDE_DIR}" "${CONCERT_INCLUDE_DIR}")
target_link_directories(cplex_runtime INTERFACE "${CPLEX_LIBRARY_DIR}" "${CONCERT_LIBRARY_DIR}")
target_link_libraries(cplex_runtime INTERFACE concert ilocplex cplex m pthread dl)
target_compile_definitions(cplex_runtime INTERFACE NDEBUG IL_STD)
target_compile_options(cplex_runtime INTERFACE -m64 -fPIC -fno-strict-aliasing -fexceptions)

set(CMPUF_COMMON_SOURCES cpp/modele.cpp cpp/tri.cpp cpp/pmedian.cpp)

function(add_cmpuf_executable target driver)
  add_executable(${target} ${driver} ${CMPUF_COMMON_SOURCES} ${ARGN})
  target_include_directories(${target} PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/h")
  target_link_libraries(${target} PRIVATE cplex_runtime)
  set_target_properties(${target} PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/bin")
endfunction()

# Active CMPUF solver.
add_cmpuf_executable(cmpuf cpp/cmpuf.cpp cpp/transport.cpp)

