cmake_minimum_required(VERSION 3.10)
project(OPRF VERSION 1.0)

set(CMAKE_C_COMPILER "gcc")
set(CMAKE_CXX_COMPILER "g++")

set (CMAKE_CXX_STANDARD 17)
add_compile_options(-g3 -O2 -march=native -funroll-loops -mrdrnd -fpermissive)
set(CMAKE_CXX_FLAGS "-pg -no-pie -msse4.2 -D__SSE4_2_ -m64 -Wno-unused-result -Wno-unused-function -Wno-strict-aliasing")

# link_directories(../externals/miracl_lib/ /usr/lib/x86_64-linux-gnu/glib-2.0)
link_directories(/home/mingli/Documents/SEAL/native/examples/PolyPSI-OPRF/preprocess/src/externals/miracl_lib/ /usr/lib/x86_64-linux-gnu/glib-2.0)


# file(GLOB all_SRC
#      "*.h"
#      "*.cpp"
# )

aux_source_directory(../src/util/crypto CRYPTO)
aux_source_directory(../src/util/cuckoo_filter CUCKOO_FILTER)
aux_source_directory(../src/util/counting_quotient_filter QUOTIENT_FILTER)
aux_source_directory(../src/util UTIL)
aux_source_directory(../externals/miracl_lib MIRACL_SOURCE)
# aux_source_directory(../src/externals/Miracl/source MIRACL_SOURCE)
aux_source_directory(../src/pk-based_gattaca GATTACA_SOURCE)

# file(GLOB ALL_O ../externals/miracl_lib/*.o)


# add_executable(dh_oprf ../dh_oprf.cpp  ${CRYPTO} ${CUCKOO_FILTER} ${QUOTIENT_FILTER} ${UTIL} ${MIRACL_SOURCE})
add_executable(dh_oprf ../dh_oprf.cpp  ${CRYPTO} ${CUCKOO_FILTER} ${QUOTIENT_FILTER} ${UTIL} ${MIRACL_SOURCE} ${GATTACA_SOURCE} ../src/externals/Miracl/source/big.cpp ../src/externals/Miracl/source/ec2.cpp

)


target_link_libraries(dh_oprf /home/mingli/Documents/SEAL/native/examples/PolyPSI-OPRF/preprocess/src/externals/miracl_lib/miracl.a glib-2.0 gmp pthread ssl crypto)

target_include_directories(dh_oprf
    PRIVATE /usr/lib/x86_64-linux-gnu/glib-2.0/include /usr/include/glib-2.0/ /home/mingli/Documents/SEAL/native/examples/PolyPSI-OPRF/preprocess/externals/Miracl/include
)
# target_include_directories(dh_oprf
#     PRIVATE `pkg-config --cflags glib-2.0` /home/mingli/Documents/SEAL/native/examples/PolyPSI-OPRF/preprocess/externals/Miracl/include
# )