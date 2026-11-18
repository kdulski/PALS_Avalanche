cmake_minimum_required(VERSION 3.16)
project(PALS_avalanche LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

option(PALS_BUILD_GUI "Build the graphical interface PALS_avalanche_gui" ON)

find_package(ROOT 6.22 REQUIRED COMPONENTS Hist MathCore Core)

# ROOT-independent part (config, models, estimates, distribution analysis, text report)
add_library(pals_core STATIC
  src/Config.cpp src/Convolution.cpp src/Distribution.cpp src/FormulaModel.cpp src/LifetimeModel.cpp
  src/ModelSetup.cpp src/Nnls.cpp src/Spectrum.cpp src/TextReport.cpp)
target_include_directories(pals_core PUBLIC src third_party)
target_compile_options(pals_core PRIVATE -Wall -Wextra)

# fitting with ROOT, shared by the command-line program and the GUI
add_library(pals_root STATIC src/FitSession.cpp src/Fitter.cpp src/RootFormula.cpp src/RootIO.cpp)
target_link_libraries(pals_root PUBLIC pals_core ROOT::Hist ROOT::MathCore ROOT::Core)
target_compile_options(pals_root PRIVATE -Wall -Wextra)

add_executable(PALS_avalanche src/main.cpp)
target_link_libraries(PALS_avalanche PRIVATE pals_root)

if(PALS_BUILD_GUI)
  if(TARGET ROOT::Gui AND TARGET ROOT::Gpad)
    add_executable(PALS_avalanche_gui src/gui/main_gui.cpp src/gui/PalsGui.cpp)
    target_include_directories(PALS_avalanche_gui PRIVATE src/gui)
    target_compile_definitions(PALS_avalanche_gui PRIVATE PALS_RESOURCE_DIR="${CMAKE_CURRENT_SOURCE_DIR}/resources")
    target_link_libraries(PALS_avalanche_gui PRIVATE pals_root ROOT::Gui ROOT::Gpad)
    target_compile_options(PALS_avalanche_gui PRIVATE -Wall -Wextra)
  else()
    message(WARNING "ROOT was built without its GUI libraries: PALS_avalanche_gui is not built")
  endif()
endif()

add_executable(test_model tests/test_model.cpp)
target_link_libraries(test_model PRIVATE pals_core)
