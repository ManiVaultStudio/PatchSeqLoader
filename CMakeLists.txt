cmake_minimum_required(VERSION 3.17)

# -----------------------------------------------------------------------------
# PatchSeqDataLoader Plugin
# -----------------------------------------------------------------------------
PROJECT("PatchSeqDataLoader")

# -----------------------------------------------------------------------------
# CMake Options
# -----------------------------------------------------------------------------
set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)
set(CMAKE_AUTOMOC ON)

if(MSVC)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /W3 /DWIN32 /EHsc /MP /permissive- /Zc:__cplusplus")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} /NODEFAULTLIB:LIBCMT")
    set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} /MDd")
    set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} /MD")
endif(MSVC)

# -----------------------------------------------------------------------------
# Dependencies
# -----------------------------------------------------------------------------
find_package(Qt6 COMPONENTS Widgets WebEngineWidgets REQUIRED)

find_package(ManiVault COMPONENTS Core PointData ClusterData TextData CONFIG REQUIRED)

find_package(LEAD CONFIG REQUIRED)

# HDF5 and ZLib dependencies
include(HDF5Dependency)

# -----------------------------------------------------------------------------
# Source files
# -----------------------------------------------------------------------------
set(PLUGIN_SOURCES
    src/PatchSeqDataLoader.h
    src/PatchSeqDataLoader.cpp
    src/PatchSeqDataLoader.json
    src/DataFrame.h
    src/DataFrame.cpp
    src/MatrixData.h
    src/MatrixData.cpp
    src/InputDialog.h
    src/InputDialog.cpp
    src/ColorTaxonomy.h
    src/ColorTaxonomy.cpp
    src/Taxonomy.h
    src/Taxonomy.cpp
    src/FeatureNames.h
    src/PatchSeqFilePaths.h
    src/PatchSeqFilePaths.cpp
    src/MatrixDataLoader.h
    src/MatrixDataLoader.cpp
    src/CSVReader.h
    src/CSVReader.cpp
    src/json.hpp
)

set(MORPHOLOGY_SOURCES
    src/Morphology/SWCLoader.h
    src/Morphology/SWCLoader.cpp
)

set(EPHYS_SOURCES
    src/Electrophysiology/NWBLoader.h
    src/Electrophysiology/NWBLoader.cpp
    src/Electrophysiology/StimulusCodeMap.h
)

set(QRESOURCES
    res/met_loader_resources.qrc
)
QT6_ADD_RESOURCES(RESOURCE_FILES ${QRESOURCES})

set(PLUGIN_MOC_HEADERS
    src/PatchSeqDataLoader.h
)

source_group( Plugin FILES ${PLUGIN_SOURCES})
source_group( Morphology FILES ${MORPHOLOGY_SOURCES})
source_group( Electrophysiology FILES ${EPHYS_SOURCES})

# -----------------------------------------------------------------------------
# CMake Target
# -----------------------------------------------------------------------------
add_library(${PROJECT_NAME} SHARED ${PLUGIN_SOURCES} ${MORPHOLOGY_SOURCES} ${EPHYS_SOURCES} ${RESOURCE_FILES})

qt_wrap_cpp(LOADER_MOC ${PLUGIN_MOC_HEADERS} TARGET ${PROJECT_NAME})
target_sources(${PROJECT_NAME} PRIVATE ${LOADER_MOC})

# -----------------------------------------------------------------------------
# Target include directories
# -----------------------------------------------------------------------------
target_include_directories(${PROJECT_NAME} PRIVATE ${PROJECT_SOURCE_DIR}/src)
target_include_directories(${PROJECT_NAME} PRIVATE "${ManiVault_INCLUDE_DIR}")
target_include_directories(${PROJECT_NAME} PRIVATE ${HDF5_INCLUDE_DIR})
target_include_directories(${PROJECT_NAME} PRIVATE ${LIBRARY_INSTALL_DIR}/zlib/$<CONFIG>/include)

# -----------------------------------------------------------------------------
# Target properties
# -----------------------------------------------------------------------------
target_compile_features(${PROJECT_NAME} PRIVATE cxx_std_20)

# -----------------------------------------------------------------------------
# Target library linking
# -----------------------------------------------------------------------------
set(MV_LINK_PATH "${ManiVault_INSTALL_DIR}/$<CONFIGURATION>/lib")
set(PLUGIN_LINK_PATH "${ManiVault_INSTALL_DIR}/$<CONFIGURATION>/$<IF:$<CXX_COMPILER_ID:MSVC>,lib,Plugins>")
set(MV_LINK_SUFFIX $<IF:$<CXX_COMPILER_ID:MSVC>,${CMAKE_LINK_LIBRARY_SUFFIX},${CMAKE_SHARED_LIBRARY_SUFFIX}>)

set(CELL_MORPHOLOGY_DATA_LINK_LIBRARY "${PLUGIN_LINK_PATH}/${CMAKE_SHARED_LIBRARY_PREFIX}CellMorphologyData${MV_LINK_SUFFIX}")
set(EPHYS_DATA_LINK_LIBRARY "${PLUGIN_LINK_PATH}/${CMAKE_SHARED_LIBRARY_PREFIX}EphysData${MV_LINK_SUFFIX}")

target_link_libraries(${PROJECT_NAME} PRIVATE Qt6::Widgets)
target_link_libraries(${PROJECT_NAME} PRIVATE Qt6::WebEngineWidgets)

# Link to ManiVault and data plugins
target_link_libraries(${PROJECT_NAME} PRIVATE ManiVault::Core)
target_link_libraries(${PROJECT_NAME} PRIVATE ManiVault::PointData)
target_link_libraries(${PROJECT_NAME} PRIVATE ManiVault::ClusterData)
target_link_libraries(${PROJECT_NAME} PRIVATE ManiVault::TextData)

target_link_libraries(${PROJECT_NAME} PRIVATE "${CELL_MORPHOLOGY_DATA_LINK_LIBRARY}")
target_link_libraries(${PROJECT_NAME} PRIVATE "${EPHYS_DATA_LINK_LIBRARY}")
target_link_libraries(${PROJECT_NAME} PRIVATE ${HDF5_CXX_STATIC_LIBRARY} ${ZLIB_LIBRARIES})

target_link_libraries(${PROJECT_NAME} PRIVATE LEAD)

# -----------------------------------------------------------------------------
# Target installation
# -----------------------------------------------------------------------------
install(TARGETS ${PROJECT_NAME}
    RUNTIME DESTINATION Plugins COMPONENT PLUGINS # Windows .dll
    LIBRARY DESTINATION Plugins COMPONENT PLUGINS # Linux/Mac .so
)

add_custom_command(TARGET ${PROJECT_NAME} POST_BUILD
    DEPENDS ${PROJECT_NAME}
    COMMAND "${CMAKE_COMMAND}"
        --install ${CMAKE_CURRENT_BINARY_DIR}
        --config $<CONFIGURATION>
        --prefix ${ManiVault_INSTALL_DIR}/$<CONFIGURATION>
)

# -----------------------------------------------------------------------------
# Miscellaneous
# -----------------------------------------------------------------------------
# Automatically set the debug environment (command + working directory) for MSVC
if(MSVC)
    set_property(TARGET ${PROJECT_NAME} PROPERTY VS_DEBUGGER_WORKING_DIRECTORY $<IF:$<CONFIG:DEBUG>,${ManiVault_INSTALL_DIR}/debug,${ManiVault_INSTALL_DIR}/release>)
    set_property(TARGET ${PROJECT_NAME} PROPERTY VS_DEBUGGER_COMMAND $<IF:$<CONFIG:DEBUG>,${ManiVault_INSTALL_DIR}/debug/ManiVault\ Studio.exe,${ManiVault_INSTALL_DIR}/release/ManiVault\ Studio.exe>)
endif()
