cmake_minimum_required(VERSION 3.8)
set(CMAKE_CXX_STANDARD 17)
project(SplineSmoother)

set(CMAKE_BUILD_TYPE "Debug")
##test

# matplotlib-cpp library
find_package(Python3 COMPONENTS NumPy)

# spdlog library
#find_package(spdlog REQUIRED)

# osqp library
find_package(OsqpEigen)
find_package(Eigen3)
find_package(qpOASES QUIET)

# compile main.cpp
include_directories(
    ./include/
    /usr/local/include/osqp
    ${EIGEN3_INCLUDE_DIR}
)

add_executable(${CMAKE_PROJECT_NAME}
                src/main.cpp
                src/SplineSmoother.cpp
)

target_link_libraries(${CMAKE_PROJECT_NAME} 
                        Python3::Python
                        #spdlog::spdlog
                        osqp OsqpEigen::OsqpEigen -lxerces-c qpOASES
                        )
