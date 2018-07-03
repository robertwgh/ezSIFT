cmake_minimum_required(VERSION 2.8)

add_subdirectory(feature_extract)
add_subdirectory(image_match)

file(GLOB_RECURSE TEST_IMAGES "../data/*.pgm")
file(COPY ${TEST_IMAGES} DESTINATION ${CMAKE_RUNTIME_OUTPUT_DIRECTORY})
