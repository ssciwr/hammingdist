pybind11_add_module(hammingdist hammingdist.cc)
target_link_libraries(hammingdist PUBLIC hamming)
install(TARGETS hammingdist DESTINATION .)
