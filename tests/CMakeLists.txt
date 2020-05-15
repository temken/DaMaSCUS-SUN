
add_executable(tests
	test_module.cpp
	${SRC_DIR}/module.cpp)

target_link_libraries(tests
	PUBLIC coverage_config 
	PRIVATE gtest_main	#contains the main function
)

target_compile_options(tests PUBLIC -Wall -pedantic)

target_include_directories(tests
    PRIVATE
        ${GENERATED_DIR} )

install(TARGETS tests DESTINATION ${TESTS_DIR})

add_test(NAME all_tests COMMAND tests
	WORKING_DIRECTORY ${TESTS_DIR})
