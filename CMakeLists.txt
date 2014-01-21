cmake_minimum_required(VERSION 2.8)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}")

#find_package(GSL REQUIRED)
find_package(Boost REQUIRED)
find_package(Threads REQUIRED)
#find_package(ROOT)

include_directories("include" ${Boost_INCLUDE_DIR} )

set(CMAKE_CXX_FLAGS ${Boost_CFLAGS})

file(GLOB sources src/*.cc)


add_executable(TA2 Tiny2.cc ${sources})
target_link_libraries(TA2 ${Boost_LIBRARIES} ${CMAKE_THREAD_LIBS_INIT})

