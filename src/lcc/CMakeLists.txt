cmake_minimum_required(VERSION 3.17)
project(lcc)

set(CMAKE_CXX_STANDARD 11)

if (NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Debug)
endif()
message(STATUS "Build configuration: " ${CMAKE_BUILD_TYPE})

set(CMAKE_INCLUDE_CURRENT_DIR ON)

find_package(FLEX)
find_package(BISON)

if(WIN32)
	set(FLEX_FLAGS "--nounistd")
else()
	set(FLEX_FLAGS "")
endif()

BISON_TARGET(LPPParser ${PROJECT_SOURCE_DIR}/lpp.y ${CMAKE_CURRENT_BINARY_DIR}/lpp.y.cpp)
FLEX_TARGET(LPPScanner ${PROJECT_SOURCE_DIR}/lpp.l ${CMAKE_CURRENT_BINARY_DIR}/lpp.l.cpp COMPILE_FLAGS ${FLEX_FLAGS})
ADD_FLEX_BISON_DEPENDENCY(LPPScanner LPPParser)

set(SOURCEFILE
	main.cpp
	node.cpp
	option.cpp
	context.cpp
	util.cpp
        number.cpp
	simulation.cpp
	state.cpp
	${BISON_LPPParser_OUTPUTS}
	${FLEX_LPPScanner_OUTPUTS}
)

if(WIN32)
	set(SOURCEFILE ${SOURCEFILE} windows/getopt.c)
	add_definitions(-DYY_NO_UNISTD_H)
endif()

add_executable(lcc
	${SOURCEFILE}
)
