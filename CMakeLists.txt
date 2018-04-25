cmake_minimum_required(VERSION 3.9)
project(KRY_RSA)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_FLAGS "-lgmp")


add_executable(KRY_RSA main.cpp)