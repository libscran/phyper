cmake_minimum_required(VERSION 3.14)

project(phyper
    VERSION 0.1.0
    DESCRIPTION "Hypergeometric tail calculations"
    LANGUAGES CXX)

include(GNUInstallDirs)
include(CMakePackageConfigHelpers)

# Library
add_library(phyper INTERFACE)
add_library(libscran::phyper ALIAS phyper)

target_include_directories(phyper INTERFACE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/phyper>)
target_compile_features(phyper INTERFACE cxx_std_17)

# Tests
if(CMAKE_PROJECT_NAME STREQUAL PROJECT_NAME)
    option(PHYPER_TESTS "Build phyper's test suite." ON)
else()
    option(PHYPER_TESTS "Build phyper's test suite." OFF)
endif()

if(PHYPER_TESTS)
    include(CTest)
    if(BUILD_TESTING)
        add_subdirectory(tests)
    endif() 
endif()

# Install
install(DIRECTORY include/
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/phyper)

install(TARGETS phyper
    EXPORT phyperTargets)

install(EXPORT phyperTargets
    FILE libscran_phyperTargets.cmake
    NAMESPACE libscran::
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/libscran_phyper)

configure_package_config_file(${CMAKE_CURRENT_SOURCE_DIR}/cmake/Config.cmake.in
    "${CMAKE_CURRENT_BINARY_DIR}/libscran_phyperConfig.cmake"
    INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/libscran_phyper)

write_basic_package_version_file(
    "${CMAKE_CURRENT_BINARY_DIR}/libscran_phyperConfigVersion.cmake"
    COMPATIBILITY SameMajorVersion)

install(FILES "${CMAKE_CURRENT_BINARY_DIR}/libscran_phyperConfig.cmake"
    "${CMAKE_CURRENT_BINARY_DIR}/libscran_phyperConfigVersion.cmake"
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/libscran_phyper)
