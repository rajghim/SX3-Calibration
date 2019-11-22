cmake_minimum_required(VERSION 2.8)
project(Analysis)

set(CMAKE_CXX_STANDARD 11)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules/")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -O3")

#Find ROOT
find_package(ROOT REQUIRED)
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}")
include_directories(${ROOT_INCLUDE_DIR})
link_directories(${ROOT_LIBRARY_DIR})

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/build")

set(SCRIPTS
SX3pedestals.dat
SX3gains.dat
SX3EnCal.dat
SX3PosCalinfo.dat
SX3PosCal.dat
Sectorpedestals.dat
SX3SectorEnCalChannels.dat
SX3SectorEnCal.dat
    )

foreach(_script ${SCRIPTS})
    configure_file(
        ${PROJECT_SOURCE_DIR}/${_script}
        ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/${_script}
        COPYONLY
    )
endforeach()

add_executable(gainmatch Gainmatch.cpp)

target_link_libraries(gainmatch ${ROOT_LIBRARIES})

add_executable(energycalibrate EnergyCalibration.cpp)

target_link_libraries(energycalibrate ${ROOT_LIBRARIES} Spectrum)

add_executable(positioncalibrate PositionCalibration.cpp)

target_link_libraries(positioncalibrate ${ROOT_LIBRARIES})

add_executable(fronttest Fronttest.cpp)

target_link_libraries(fronttest ${ROOT_LIBRARIES})

add_executable(sectorcalibrate SectorCalibration.cpp)

target_link_libraries(sectorcalibrate ${ROOT_LIBRARIES} Spectrum)

