cmake_minimum_required(VERSION 3.10)
project(BFV VERSION 1.0)
link_directories(../hash/lib/)

add_executable(low_poly_psi ../low_poly_psi.cpp)
add_executable(low_poly_psi_nx28 ../low_poly_psi_nx28.cpp)


find_package(SEAL 3.6 REQUIRED)

target_link_libraries(low_poly_psi SEAL::seal blake3)
target_link_libraries(low_poly_psi_nx28 SEAL::seal blake3)
