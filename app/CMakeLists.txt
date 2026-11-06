cmake_minimum_required(VERSION 3.10)
set(CMAKE_CXX_STANDARD 17)

# Local
set(SOURCE_FILES
    main.cpp)
    
add_executable(learnSPH_app ${SOURCE_FILES})
add_dependencies(learnSPH_app learnSPH)
target_link_libraries(learnSPH_app learnSPH)

add_executable(learnSPH_benchmark benchmark.cpp)
add_dependencies(learnSPH_benchmark learnSPH)
target_link_libraries(learnSPH_benchmark learnSPH)

add_executable(learnSPH_experiments experiments.cpp)
add_dependencies(learnSPH_experiments learnSPH)
target_link_libraries(learnSPH_experiments learnSPH)
# Assignment 1(b): standalone density experiment.
add_executable(learnSPH_density density.cpp)
target_link_libraries(learnSPH_density learnSPH)