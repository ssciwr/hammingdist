project(covid-tda LANGUAGES CXX)

cmake_minimum_required(VERSION 3.11)

find_package(OpenMP)

find_package(Eigen3)

find_package(blaze)
add_executable(blaze blaze.cc)

find_package(SeqAn3 HINTS "${CMAKE_SOURCE_DIR}/seqan3/build_system")

add_executable(eigen eigen.cc)
target_link_libraries(eigen Eigen3::Eigen)

add_executable(stdvector stdvector.cc)

add_executable(playground playground.cc)
target_link_libraries(playground PUBLIC OpenMP::OpenMP_CXX)

add_executable(seqan seqan.cc)
target_link_libraries(seqan PUBLIC seqan3::seqan3)

include(FeatureSummary)
feature_summary(WHAT ALL)
