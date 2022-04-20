cmake_minimum_required(VERSION 3.10)
project(BFV VERSION 1.0)
link_directories(../hash/lib/)

add_executable(low_poly_psi ../low_poly_psi.cpp)


find_package(SEAL 3.6 REQUIRED)

target_link_libraries(low_poly_psi SEAL::seal blake3)
