#
# If the user specifies -DCMAKE_BUILD_TYPE on the command line, take their
# definition
# and dump it in the cache along with proper documentation, otherwise set
# CMAKE_BUILD_TYPE
# to Debug prior to calling PROJECT()
#


# Project Definition
project(CPToyMC)
cmake_minimum_required(VERSION 2.6)

set(CMAKE_CXX_FLAGS_DBG "-O0 -ggdb -pg" CACHE STRING "Debug options." FORCE)
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -DNDEBUG" CACHE STRING "Debug options." FORCE)
set(CMAKE_CXX_FLAGS_PROFILING "-O3 -pg" CACHE STRING "Debug options." FORCE)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules/")

if(${APPLE})
  message( STATUS "Building for Mac OS X, switching on C++11 flags for Mac OS X/clang" )
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -stdlib=libc++ -pedantic -Wall -Wextra")
endif(${APPLE})
if(${CMAKE_SYSTEM_NAME} MATCHES "Linux")
  message( STATUS "Building for Linux, switching on C++11 flags for Linux/gcc" )
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++1y -pedantic -Wall -Wextra")
endif(${CMAKE_SYSTEM_NAME} MATCHES "Linux")

message( STATUS "CXX flags are ${CMAKE_CXX_FLAGS}")

set(BASEPATH "${CMAKE_SOURCE_DIR}/src")
include_directories(${BASEPATH})

set(LIBRARY_OUTPUT_PATH "${CMAKE_BINARY_DIR}/lib")
set(EXECUTABLE_OUTPUT_PATH "${CMAKE_BINARY_DIR}/bin")



find_package(Boost 1.55 COMPONENTS program_options filesystem thread system random REQUIRED)
include_directories(SYSTEM ${Boost_INCLUDE_DIR})
link_directories(${Boost_LIBRARY_DIR})

find_package(ROOT REQUIRED)
include_directories(${ROOT_INCLUDE_DIR})
link_directories(${ROOT_LIBRARY_DIR})

set(ALL_LIBRARIES ${ROOT_LIBRARIES} ${ROOFIT_LIBRARIES} ${ADDITIONAL_LIBRARIES} "-lTreePlayer")

if(DEFINED ENV{DOOCORESYS})
  set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "$ENV{DOOCORESYS}/cmake/Modules/")
  find_package(DooCore REQUIRED)
  include_directories(${DooCore_INCLUDE_DIR})
  link_directories(${DooCore_LIBRARY_DIR})
  set(DOOSOFT_LIBRARIES ${DooCore_LIBRARIES})

  if(DEFINED ENV{DOOFITSYS})
    set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "$ENV{DOOFITSYS}/cmake/Modules/")
    find_package(DooFit REQUIRED)
    include_directories(${DooFit_INCLUDE_DIR})
    link_directories(${DooFit_LIBRARY_DIR})
    set(DOOSOFT_LIBRARIES ${DOOSOFT_LIBRARIES} ${DooFit_LIBRARIES})
  endif()

  if(DEFINED ENV{DOOSELECTIONSYS})
    set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "$ENV{DOOSELECTIONSYS}/cmake/Modules/")
    find_package(DooSelection REQUIRED)
    include_directories(${DooSelection_INCLUDE_DIR})
    link_directories(${DooSelection_LIBRARY_DIR})
    set(DOOSOFT_LIBRARIES ${DOOSOFT_LIBRARIES} ${DooSelection_LIBRARIES})
  endif()
else()
  message( STATUS "Could not find DOOCORESYS environment variable. Omitting all DooSoftware dependent parts.")
endif()



add_subdirectory(main)
add_subdirectory(sandbox)
add_subdirectory(src)

if(DEFINED ENV{DOOCORESYS} AND DEFINED ENV{DOOFITSYS})
  add_subdirectory(tests)
endif()
