add_executable(tests
	test_Simulation_Trajectory.cpp
	${SRC_DIR}/Simulation_Trajectory.cpp
	test_Simulation_Utilities.cpp
	${SRC_DIR}/Simulation_Utilities.cpp
	test_Solar_Model.cpp
	${SRC_DIR}/Solar_Model.cpp)

target_link_libraries(tests
	PUBLIC
		coverage_config 
	PRIVATE
		gtest_main	#contains the main function
		libobscura )

target_include_directories(tests PRIVATE ${GENERATED_DIR} )

target_compile_options(tests PUBLIC -Wall -pedantic)

install(TARGETS tests DESTINATION ${TESTS_DIR})

add_test(NAME all_tests COMMAND tests
	WORKING_DIRECTORY ${TESTS_DIR})