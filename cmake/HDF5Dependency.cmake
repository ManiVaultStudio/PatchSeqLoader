#set(hdf5_VERSION 1.12.1)
SET(hdf5_VERSION "1.14.2" CACHE STRING "Version of HDF5 Library")
SET_PROPERTY(CACHE hdf5_VERSION PROPERTY STRINGS 1.12.1 1.14.2)

set(USE_HDF5_ARTIFACTORY_LIBS FALSE CACHE BOOL "Use the prebuilt libraries from artifactory")
if(NOT USE_HDF5_ARTIFACTORY_LIBS)
    set(HDF5_ARTIFACTORY_LIBS_INSTALLED FALSE CACHE BOOL "The prebuild libraries from artifactory are installed")
endif()

include(InstallArtifactoryPackage)
set(LIBRARY_INSTALL_DIR ${PROJECT_BINARY_DIR})
if (USE_HDF5_ARTIFACTORY_LIBS)
    if (NOT HDF5_ARTIFACTORY_LIBS_INSTALLED) 
        message(STATUS "Installing artifactory packages to: ${LIBRARY_INSTALL_DIR}")
        # Both HDILib and flann are available prebuilt in the lkeb-artifactory as combined Debug/Release packages
        # lz4 is also available in the lkb-artifactory in separate Debug and |Release packages
        install_artifactory_package(PACKAGE_NAME hdf5 PACKAGE_VERSION ${hdf5_VERSION} PACKAGE_BUILDER lkeb COMBINED_PACKAGE TRUE) 

        message(STATUS "HDF5 root path ${hdf5_ROOT}")
    endif()
    set(HDF5_ROOT ${LIBRARY_INSTALL_DIR}/hdf5)
    set(ZLIB_ROOT ${LIBRARY_INSTALL_DIR}/hdf5/cmake CACHE PATH "Path to zlib root" FORCE)
    set(HDF5_USE_STATIC_LIBRARIES TRUE)
    find_package(HDF5 COMPONENTS CXX C static REQUIRED NO_MODULE)
    find_package (ZLIB CONFIG REQUIRED)

    message(STATUS "Include for HDF5 at ${HDF5_INCLUDE_DIR} - version ${HDF5_VERSION_STRING}")
    set(ARTIFACTORY_LIBS_INSTALLED TRUE CACHE BOOL "Use the prebuilt libraries from artifactory" FORCE)
    message(status "hdf5 include ${HDF5_INCLUDE_DIR}")
endif()
