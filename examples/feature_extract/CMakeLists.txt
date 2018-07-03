cmake_minimum_required(VERSION 2.8)

file(GLOB_RECURSE EXAMPLE_SRC "feature_extract.cpp")
add_executable(feature_extract ${EXAMPLE_SRC})

target_link_libraries(feature_extract
	ezsift
)
