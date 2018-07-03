cmake_minimum_required(VERSION 2.8)

file(GLOB_RECURSE EXAMPLE_SRC "image_match.cpp")
add_executable(image_match ${EXAMPLE_SRC})

target_link_libraries(image_match
	ezsift
)
