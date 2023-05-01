cmake_minimum_required(VERSION 3.1...3.23)
project(VW_Simplification)

set(CMAKE_BUILD_TYPE Release)

find_package(CGAL REQUIRED)

# Add header and source files to project
set(HEADER_FILES
    Core.hpp
    Fast_simplifier.hpp
    Slow_simplifier.hpp
    Test.hpp
    Utils.hpp
)

set(SOURCE_FILES
    Fast_simplifier.cpp
    Slow_simplifier.cpp
    Test.cpp
    Utils.cpp
    main.cpp
)

# Create executable target
add_executable(${PROJECT_NAME} ${HEADER_FILES} ${SOURCE_FILES})

# Link target to CGAL libraries
target_link_libraries(${PROJECT_NAME} CGAL::CGAL)
