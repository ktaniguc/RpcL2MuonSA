cmake_minimum_required(VERSION 2.8)

list(APPEND CMAKE_PREFIX_PATH  $ENV{ROOTSYS})
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} $ENV{ROOTSYS}/etc/cmake)

find_package(ROOT REQUIRED COMPONENTS Hist MathCore Tree RIO Net)

include(${ROOT_USE_FILE})
include_directories(${ROOT_INCLUDE_DIRS})
add_definitions(${ROOT_CXX_FLAGS})
add_library(src 
    SHARED
    EfficiencyFiller.cxx
    EventSelectionTool.cxx
    RPC.cxx
    RPCHitRegionSetter.cxx
    RPC_FCBM.cxx
    RpcClusterInfoFiller.cxx
    RpcClusterSetter.cxx
    RpcClusteringTool.cxx
    RpcDispPreparator.cxx
    RpcDispRZSetter.cxx
    RpcFCBMHistDrawer.cxx
    RpcFCBMHistFiller.cxx
    RpcFitEtaRSetter.cxx
    setLegend.cxx
)
add_library(AtlasStyle
    SHARED
    ../AtlasStyle/AtlasExample.C
    ../AtlasStyle/AtlasLabels.C
    ../AtlasStyle/AtlasStyle.C
    ../AtlasStyle/AtlasUtils.C
)
target_include_directories(src
    PUBLIC 
    ${PROJECT_SOURCE_DIR}/RpcL2MuonSA
    ${ROOT_INCLUDE_DIRS}
    ../AtlasStyle
)
target_link_libraries(src ${ROOT_LIBRARIES} )
