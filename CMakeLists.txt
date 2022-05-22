cmake_minimum_required(VERSION 3.21)
project(N15data)

include_directories(.)
include_directories(match)
include_directories(setup)

set(CMAKE_CXX_STANDARD 14)

find_package(AUSALIB REQUIRED)

add_executable(AlphaAccCalib src/AlphaAccCalib/AlphaAccCalibRunner.cpp src/AlphaAccCalib/AlphaAccECorrector.cpp src/AlphaAccCalib/AlphacmEFitter.cpp src/AlphaAccCalib/AlphaFactorFitter.cpp)
target_link_libraries(AlphaAccCalib PRIVATE AUSA::AUSA)
add_executable(AccCalibRunner src/AccCalib/AccCalibrationRunner.cpp src/AccCalib/calibrationEcorrector.cpp src/AlphaAccCalib/cmEgaussFitter.cpp src/AccCalib/factorFitter.cpp)
target_link_libraries(AccCalibRunner PRIVATE AUSA::AUSA)
add_executable(Alpha0Cross src/Alpha0Cross/Alpha0crossRunner.cpp src/AlphaAccCalib/AlphacmEFitter.cpp src/VCHARGEfinder.cpp src/AlphaAccCalib/AlphaAccECorrector.cpp src/Alpha0Cross/Alpha0AngularCross.cpp)
target_link_libraries(Alpha0Cross PRIVATE AUSA::AUSA)
add_executable(coinRun src/Coincidence/coinRunner.cpp src/Coincidence/coinECorrector.cpp src/Coincidence/coinCross.cpp src/AlphaAccCalib/AlphacmEFitter.cpp)
target_link_libraries(coinRun PRIVATE AUSA::AUSA)

add_executable(fixFile src/fileFixer.cpp src/fileFixRunner.cpp)
target_link_libraries(fixFile PRIVATE AUSA::AUSA)

add_executable(scatterRun src/AccCalib/ScatterAngRunner.cpp src/AccCalib/ScatterAngularCross.cpp src/AccCalib/calibrationEcorrector.cpp src/VCHARGEfinder.cpp src/AlphaAccCalib/AlphacmEFitter.cpp)
target_link_libraries(scatterRun PRIVATE AUSA::AUSA)