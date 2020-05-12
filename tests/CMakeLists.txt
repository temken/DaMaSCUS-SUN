
#######################
# Option A: If a static library is built. (see src/CMakeLists.txt)
add_executable(tests
	test_module.cpp)
target_link_libraries(tests 
	PRIVATE
		module_library
		gtest_main	#contains the main function
)
# (In this case coverage_config is contained in module_library.)
#######################

#######################
# # Option B: If no static library is built. (see src/CMakeLists.txt)
# add_executable(tests
# 	test_module.cpp
# 	${SRC_DIR}/module.cpp)

# target_link_libraries(tests
# 	PUBLIC coverage_config 
# 	PRIVATE gtest_main	#contains the main function
# )
#######################

target_compile_options(tests PUBLIC -Wall -pedantic)

install(TARGETS tests DESTINATION ${TESTS_DIR})

add_test(NAME all_tests COMMAND tests
	WORKING_DIRECTORY ${TESTS_DIR})
