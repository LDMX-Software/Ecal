
# Set the minimum version of CMake that's required
cmake_minimum_required(VERSION 3.12)

# Set the project name
project(Ecal VERSION 2.1.0
             DESCRIPTION "Module for all things ECal."
             LANGUAGES CXX
)

# Search and configure ROOT
find_package(ROOT CONFIG REQUIRED)

option(BUILD_EVENT_ONLY "Build the event library." ON)
if(BUILD_EVENT_ONLY)
    
    register_event_object( module_path "Ecal/Event" namespace "ecal::event" class "ClusterAlgoResult" )
    register_event_object( module_path "Ecal/Event" namespace "ecal::event" class "EcalVetoResult" )
    register_event_object( module_path "Ecal/Event" namespace "ecal::event" class "EcalCluster" type "collection" )
    register_event_object( module_path "Ecal/Event" namespace "ecal::event" class "EcalHit"     type "collection" )
   
    # Generate the files needed to build the event classes.
    setup_library(module Ecal
                  name Event
                  dependencies ROOT::Core 
                               Recon::Event
                  register_target)

    return()

endif()

# Search for the ONNXRuntime binaries.  If not found, download the binaries 
# to an external folder and create the targets.
find_package(ONNXRuntime 1.2.0)

setup_library(module Ecal
              dependencies ROOT::Physics 
                           Framework::Framework 
                           Tools::Tools 
                           DetDescr::DetDescr 
                           SimCore::Event
                           Ecal::Event
                           Recon::Event 
                           ONNXRuntime::Interface
)

setup_test(dependencies Ecal::Ecal
           configs ${PROJECT_SOURCE_DIR}/test/ecal_digi_pipeline.py)

setup_python(package_name ${PYTHON_PACKAGE_NAME}/Ecal)

setup_data(module Ecal)

# Add the printing executable
add_executable(print-ecal-hex-readout ${PROJECT_SOURCE_DIR}/src/Ecal/print_ecal_hex_readout.cxx)

# Link to the Framework library
target_link_libraries(print-ecal-hex-readout PRIVATE DetDescr Ecal)

# Install the fire executable
install(TARGETS print-ecal-hex-readout DESTINATION ${CMAKE_INSTALL_PREFIX}/bin)
