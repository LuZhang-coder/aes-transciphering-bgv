cmake_minimum_required(VERSION 3.13)
project(TEST_AES)

find_package(helib REQUIRED)

add_executable(TEST_AES simpleAES.cpp homAES.cpp Test_AES.cpp)
target_link_libraries(TEST_AES PUBLIC helib)
if(USE_PAPER_VERSION)
  target_compile_definitions(TEST_AES PUBLIC USE_PAPER_VERSION)
endif()
