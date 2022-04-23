cmake_minimum_required(VERSION 3.21)
project(N15data)

include_directories(.)
include_directories(match)
include_directories(setup)

set(CMAKE_CXX_STANDARD 14)

find_package(AUSALIB REQUIRED)

add_executable(AlphaAccCalib AlphaAccCalibRunner.cpp AlphaAccECorrector.cpp cmEAlphaFitter.cpp AlphaFactorFitter.cpp)
target_link_libraries(AlphaAccCalib PRIVATE AUSA::AUSA)
add_executable(AccCalibRunner AccCalibrationRunner.cpp calibrationEcorrector.cpp cmEgaussFitter.cpp factorFitter.cpp)
target_link_libraries(AccCalibRunner PRIVATE AUSA::AUSA)