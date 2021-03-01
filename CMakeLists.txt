# setup the cmake version
cmake_minimum_required(VERSION 2.8)

# define project name
project(RpcL2MuonSA CXX)

# make .out from .cxx
add_subdirectory(src AtlasStyle)

add_executable(RpcL2MuonSA.out src/RPC.cxx)
target_link_libraries(RpcL2MuonSA.out src AtlasStyle)
