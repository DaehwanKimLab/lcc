cmake_minimum_required(VERSION 3.17)
project(lcc)

set(CMAKE_CXX_STANDARD 11)

include_directories(${PROJECT_SOURCE_DIR})

find_package(FLEX)
find_package(BISON)

BISON_TARGET(LPPParser lpp.y ${CMAKE_CURRENT_BINARY_DIR}/lpp.y.cpp)
FLEX_TARGET(LPPScanner lpp.l ${CMAKE_CURRENT_BINARY_DIR}/lpp.l.cpp) 
ADD_FLEX_BISON_DEPENDENCY(LPPScanner LPPParser)

add_executable(lcc-prep
	main.cpp
	node.cpp
	option.cpp
	context.cpp
	util.cpp
	${BISON_LPPParser_OUTPUTS}
	${FLEX_LPPScanner_OUTPUTS}
)
