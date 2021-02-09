# 1. Test executable
add_test(NAME executable_test COMMAND mpirun DaMaSCUS-SUN config.cfg
	WORKING_DIRECTORY ${BIN_DIR})

# 2. Unit tests: Dark_Photon.xpp	
add_executable(test_Dark_Photon
	test_Dark_Photon.cpp)
target_include_directories(test_Dark_Photon PRIVATE ${GENERATED_DIR} )
target_link_libraries(test_Dark_Photon
	PRIVATE
		lib_damascus_sun
		gtest_main)	#contains the main function
target_compile_options(test_Dark_Photon PUBLIC -Wall -pedantic)
install(TARGETS test_Dark_Photon DESTINATION ${TESTS_DIR})
add_test(NAME test_Dark_Photon COMMAND test_Dark_Photon
	WORKING_DIRECTORY ${TESTS_DIR})

# 3. Unit tests: Data_Generation.xpp	
add_executable(test_Data_Generation
	test_Data_Generation.cpp)
target_include_directories(test_Data_Generation PRIVATE ${GENERATED_DIR} )
target_link_libraries(test_Data_Generation
	PRIVATE
		lib_damascus_sun
		gtest_main)	#contains the main function
target_compile_options(test_Data_Generation PUBLIC -Wall -pedantic)
install(TARGETS test_Data_Generation DESTINATION ${TESTS_DIR})
add_test(NAME test_Data_Generation COMMAND test_Data_Generation
	WORKING_DIRECTORY ${TESTS_DIR})

# 4. Unit tests: Parameter_Scan.xpp	
add_executable(test_Parameter_Scan
	test_Parameter_Scan.cpp)
target_include_directories(test_Parameter_Scan PRIVATE ${GENERATED_DIR} )
target_link_libraries(test_Parameter_Scan
	PRIVATE
		lib_damascus_sun
		gtest_main)	#contains the main function
target_compile_options(test_Parameter_Scan PUBLIC -Wall -pedantic)
install(TARGETS test_Parameter_Scan DESTINATION ${TESTS_DIR})
add_test(NAME test_Parameter_Scan COMMAND test_Parameter_Scan
	WORKING_DIRECTORY ${TESTS_DIR})

# 5. Unit tests: Reflection_Spectrum.xpp	
add_executable(test_Reflection_Spectrum
	test_Reflection_Spectrum.cpp)
target_include_directories(test_Reflection_Spectrum PRIVATE ${GENERATED_DIR} )
target_link_libraries(test_Reflection_Spectrum
	PRIVATE
		lib_damascus_sun
		gtest_main)	#contains the main function
target_compile_options(test_Reflection_Spectrum PUBLIC -Wall -pedantic)
install(TARGETS test_Reflection_Spectrum DESTINATION ${TESTS_DIR})
add_test(NAME test_Reflection_Spectrum COMMAND test_Reflection_Spectrum
	WORKING_DIRECTORY ${TESTS_DIR})

# 6. Unit tests: Simulation_Trajectory.xpp	
add_executable(test_Simulation_Trajectory
	test_Simulation_Trajectory.cpp)
target_include_directories(test_Simulation_Trajectory PRIVATE ${GENERATED_DIR} )
target_link_libraries(test_Simulation_Trajectory
	PRIVATE
		lib_damascus_sun
		gtest_main)	#contains the main function
target_compile_options(test_Simulation_Trajectory PUBLIC -Wall -pedantic)
install(TARGETS test_Simulation_Trajectory DESTINATION ${TESTS_DIR})
add_test(NAME test_Simulation_Trajectory COMMAND test_Simulation_Trajectory
	WORKING_DIRECTORY ${TESTS_DIR})

# 7. Unit tests: Simulation_Utilties.xpp	
add_executable(test_Simulation_Utilities
	test_Simulation_Utilities.cpp)
target_include_directories(test_Simulation_Utilities PRIVATE ${GENERATED_DIR} )
target_link_libraries(test_Simulation_Utilities
	PRIVATE
		lib_damascus_sun
		gtest_main)	#contains the main function
target_compile_options(test_Simulation_Utilities PUBLIC -Wall -pedantic)
install(TARGETS test_Simulation_Utilities DESTINATION ${TESTS_DIR})
add_test(NAME test_Simulation_Utilities COMMAND test_Simulation_Utilities
	WORKING_DIRECTORY ${TESTS_DIR})

# 8. Unit tests: Solar_Model.xpp	
add_executable(test_Solar_Model
	test_Solar_Model.cpp)
target_include_directories(test_Solar_Model PRIVATE ${GENERATED_DIR} )
target_link_libraries(test_Solar_Model
	PRIVATE
		lib_damascus_sun
		gtest_main)	#contains the main function
target_compile_options(test_Solar_Model PUBLIC -Wall -pedantic)
install(TARGETS test_Solar_Model DESTINATION ${TESTS_DIR})
add_test(NAME test_Solar_Model COMMAND test_Solar_Model
	WORKING_DIRECTORY ${TESTS_DIR})