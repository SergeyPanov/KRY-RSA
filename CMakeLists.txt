cmake_minimum_required(VERSION 3.9)
project(KRY_RSA)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_FLAGS "-lgmp")
set(OPENSSL_CFLAGS "-I~/ssl")
set(OPENSSL_LIBS "-L~/ssl/lib -lcrypto -lssl")

add_executable(KRY_RSA main.cpp)