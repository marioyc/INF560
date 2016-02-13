#
#  @file CMakeLists.txt
#  @brief basic example CMakeLists.txt file
#  @remarks call from:
#  @author Abal-Kassim Cheik Ahamed, Frédéric Magoulès, Sonia Toubaline
#  @author
#  @date Tue Nov 24 16:16:48 CET 2015
#  @version 1.0
#  @remarks
#  @note
#

CMAKE_MINIMUM_REQUIRED(VERSION 2.6)

# -- project name
SET(PROJECT_NAME td1)

PROJECT(${PROJECT_NAME})

# -- include current directory
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})


# -- find package MPI
FIND_PACKAGE(MPI REQUIRED)

# -- include directories
INCLUDE_DIRECTORIES(${MPI_INCLUDE_PATH})


# -- source files
SET(SRC_NAMES
  dllmrg.cpp
  Vector.cpp
  MatrixDense.cpp
  DataTopology.cpp
  BlasMpi.cpp
  BlasMpi1.cpp
  Jacobi.cpp
# -- project
  Schur.cpp
  DirectSolver.cpp
)

# ------------------------------------------------------------------------------
# -- create library
# ------------------------------------------------------------------------------

ADD_LIBRARY(${PROJECT_NAME} ${SRC_NAMES})
TARGET_LINK_LIBRARIES(${PROJECT_NAME} ${MPI_LIBRARIES})

# ------------------------------------------------------------------------------
# -- create executable: projet
# ------------------------------------------------------------------------------

SET(DEMO_DIR demo/projet)


# ------------------------------------------------------------------------------
# -- add executable: DemoSchur
# ------------------------------------------------------------------------------

SET(DEMO_NAME DemoSchur)
ADD_EXECUTABLE(${DEMO_NAME} ${DEMO_DIR}/${DEMO_NAME}.cpp)
TARGET_LINK_LIBRARIES(${DEMO_NAME} ${PROJECT_NAME} ${MPI_LIBRARIES})

# ------------------------------------------------------------------------------
# -- add executable: DemoDirectSolver
# ------------------------------------------------------------------------------

SET(DEMO_NAME DemoDirectSolver)
ADD_EXECUTABLE(${DEMO_NAME} ${DEMO_DIR}/${DEMO_NAME}.cpp)
TARGET_LINK_LIBRARIES(${DEMO_NAME} ${PROJECT_NAME} ${MPI_LIBRARIES})

# ------------------------------------------------------------------------------
# -- add executable: DemoLUSolver
# ------------------------------------------------------------------------------

SET(DEMO_NAME DemoLUSolver)
ADD_EXECUTABLE(${DEMO_NAME} ${DEMO_DIR}/${DEMO_NAME}.cpp)
TARGET_LINK_LIBRARIES(${DEMO_NAME} ${PROJECT_NAME} ${MPI_LIBRARIES})


# ------------------------------------------------------------------------------
# -- create executable: td
# ------------------------------------------------------------------------------

SET(DEMO_DIR demo/td)

## ------------------------------------------------------------------------------
## -- add executable: DemoSplittingVector
## ------------------------------------------------------------------------------

#SET(DEMO_NAME DemoSplittingVector)
#ADD_EXECUTABLE(${DEMO_NAME} ${DEMO_DIR}/${DEMO_NAME}.cpp)
#TARGET_LINK_LIBRARIES(${DEMO_NAME} ${PROJECT_NAME} ${MPI_LIBRARIES})

## ------------------------------------------------------------------------------
## -- add executable: DemoSplittingMatrix
## ------------------------------------------------------------------------------

#SET(DEMO_NAME DemoSplittingMatrix)
#ADD_EXECUTABLE(${DEMO_NAME} ${DEMO_DIR}/${DEMO_NAME}.cpp)
#TARGET_LINK_LIBRARIES(${DEMO_NAME} ${PROJECT_NAME} ${MPI_LIBRARIES})

## ------------------------------------------------------------------------------
## -- add executable: DemoMVPSequential
## ------------------------------------------------------------------------------

#SET(DEMO_NAME DemoMVPSequential)
#ADD_EXECUTABLE(${DEMO_NAME} ${DEMO_DIR}/${DEMO_NAME}.cpp)
#TARGET_LINK_LIBRARIES(${DEMO_NAME} ${PROJECT_NAME} ${MPI_LIBRARIES})


## ------------------------------------------------------------------------------
## -- add executable: DemoMVPBandRow
## ------------------------------------------------------------------------------

#SET(DEMO_NAME DemoMVPBandRow)
#ADD_EXECUTABLE(${DEMO_NAME} ${DEMO_DIR}/${DEMO_NAME}.cpp)
#TARGET_LINK_LIBRARIES(${DEMO_NAME} ${PROJECT_NAME} ${MPI_LIBRARIES})

## ------------------------------------------------------------------------------
## -- add executable: DemoMVPBandColumn
## ------------------------------------------------------------------------------

#SET(DEMO_NAME DemoMVPBandColumn)
#ADD_EXECUTABLE(${DEMO_NAME} ${DEMO_DIR}/${DEMO_NAME}.cpp)
#TARGET_LINK_LIBRARIES(${DEMO_NAME} ${PROJECT_NAME} ${MPI_LIBRARIES})

## ------------------------------------------------------------------------------
## -- add executable: DemoMVPBlock
## ------------------------------------------------------------------------------

#SET(DEMO_NAME DemoMVPBlock)
#ADD_EXECUTABLE(${DEMO_NAME} ${DEMO_DIR}/${DEMO_NAME}.cpp)
#TARGET_LINK_LIBRARIES(${DEMO_NAME} ${PROJECT_NAME} ${MPI_LIBRARIES})

### ------------------------------------------------------------------------------
## -- add executable: DemoMMPBlock
## ------------------------------------------------------------------------------

#SET(DEMO_NAME DemoMMPBlock)
#ADD_EXECUTABLE(${DEMO_NAME} ${DEMO_DIR}/${DEMO_NAME}.cpp)
#TARGET_LINK_LIBRARIES(${DEMO_NAME} ${PROJECT_NAME} ${MPI_LIBRARIES})

## ------------------------------------------------------------------------------
## -- add executable: DemoReadLocal
## ------------------------------------------------------------------------------

#SET(DEMO_NAME DemoReadLocal)
#ADD_EXECUTABLE(${DEMO_NAME} ${DEMO_DIR}/${DEMO_NAME}.cpp)
#TARGET_LINK_LIBRARIES(${DEMO_NAME} ${PROJECT_NAME} ${MPI_LIBRARIES})

## ------------------------------------------------------------------------------
## -- add executable: DemoJacobiSequential
## ------------------------------------------------------------------------------

SET(DEMO_NAME DemoJacobiSequential)
ADD_EXECUTABLE(${DEMO_NAME} ${DEMO_DIR}/${DEMO_NAME}.cpp)
TARGET_LINK_LIBRARIES(${DEMO_NAME} ${PROJECT_NAME} ${MPI_LIBRARIES})


## ------------------------------------------------------------------------------
## -- add executable:  DemoJacobiBandRowSyncISyncC
## ------------------------------------------------------------------------------

#SET(DEMO_NAME  DemoJacobiBandRowSyncISyncC)
#ADD_EXECUTABLE(${DEMO_NAME} ${DEMO_DIR}/${DEMO_NAME}.cpp)
#TARGET_LINK_LIBRARIES(${DEMO_NAME} ${PROJECT_NAME} ${MPI_LIBRARIES})

## ------------------------------------------------------------------------------
## -- add executable:  DemoJacobiBandRowSyncIAsyncC
## ------------------------------------------------------------------------------

#SET(DEMO_NAME  DemoJacobiBandRowSyncIAsyncC)
#ADD_EXECUTABLE(${DEMO_NAME} ${DEMO_DIR}/${DEMO_NAME}.cpp)
#TARGET_LINK_LIBRARIES(${DEMO_NAME} ${PROJECT_NAME} ${MPI_LIBRARIES})

## ------------------------------------------------------------------------------
## -- add executable:  DemoJacobiBandRowAsyncIAsyncC
## ------------------------------------------------------------------------------

#SET(DEMO_NAME  DemoJacobiBandRowAsyncIAsyncC)
#ADD_EXECUTABLE(${DEMO_NAME} ${DEMO_DIR}/${DEMO_NAME}.cpp)
#TARGET_LINK_LIBRARIES(${DEMO_NAME} ${PROJECT_NAME} ${MPI_LIBRARIES})


## ------------------------------------------------------------------------------
## -- documentation
## ------------------------------------------------------------------------------

# The project version number.
set(VERSION_MAJOR   0   CACHE STRING "Project major version number.")
set(VERSION_MINOR   0   CACHE STRING "Project minor version number.")
set(VERSION_PATCH   1   CACHE STRING "Project patch version number.")
mark_as_advanced(VERSION_MAJOR VERSION_MINOR VERSION_PATCH)

# add a target to generate API documentation with Doxygen
find_package(Doxygen)
option(BUILD_DOCUMENTATION "Create and install the HTML based API documentation (requires Doxygen)" ${DOXYGEN_FOUND})

if(BUILD_DOCUMENTATION)
    if(NOT DOXYGEN_FOUND)
        message(FATAL_ERROR "Doxygen is needed to build the documentation.")
    endif()

    set(doxyfile_in ${CMAKE_CURRENT_SOURCE_DIR}/doc/Doxyfile.in)
    set(doxyfile ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile)

    configure_file(${doxyfile_in} ${doxyfile} @ONLY)

    add_custom_target(doc
        COMMAND ${DOXYGEN_EXECUTABLE} ${doxyfile}
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
        COMMENT "Generating API documentation with Doxygen"
        VERBATIM)

    install(DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/html DESTINATION share/doc)
endif()
