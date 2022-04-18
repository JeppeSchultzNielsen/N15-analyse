cmake_minimum_required(VERSION 3.21)
project(data)

include_directories(.)
include_directories(match)
include_directories(setup)

set(CMAKE_CXX_STANDARD 14)

find_package(AUSALIB REQUIRED)

add_executable(cross AlphaCrossRunner.cpp)
target_link_libraries(cross PRIVATE AUSA::AUSA)
add_executable(AccCalibRunner AccCalibrationRunner.cpp calibrationEcorrector.cpp cmEgaussFitter.cpp factorFitter.cpp)
target_link_libraries(AccCalibRunner PRIVATE AUSA::AUSA)