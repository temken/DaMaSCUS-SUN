# 1. Test executable
add_test(NAME executable_test COMMAND mpirun DaMaSCUS-SUN config.cfg
	WORKING_DIRECTORY ${BIN_DIR})

# 2. Test all modules
file(GLOB TESTFILES "*.cpp")

foreach(TESTFILE ${TESTFILES})
	get_filename_component(TESTNAME ${TESTFILE} NAME_WLE)
	get_filename_component(FILENAME ${TESTFILE} NAME)

	add_executable(${TESTNAME} ${FILENAME})

	target_link_libraries(${TESTNAME}
		PRIVATE
		lib_damascus_sun
		gtest_main # contains the main function
	)

	target_include_directories(${TESTNAME} PRIVATE ${GENERATED_DIR})

	target_compile_options(${TESTNAME} PUBLIC -Wall -pedantic)

	install(TARGETS ${TESTNAME} DESTINATION ${TESTS_DIR})

	add_test(NAME ${TESTNAME} COMMAND ${TESTNAME}
		WORKING_DIRECTORY ${TESTS_DIR})
endforeach()