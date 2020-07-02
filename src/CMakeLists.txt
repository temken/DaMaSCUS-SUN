# Executable
add_executable(DaMaSCUS-SUN
	main.cpp)

target_compile_options(DaMaSCUS-SUN PUBLIC -Wall -pedantic)

target_link_libraries(DaMaSCUS-SUN
	PUBLIC
		lib_damascus_sun)

target_include_directories(DaMaSCUS-SUN
	PRIVATE
		${GENERATED_DIR} )

install(TARGETS DaMaSCUS-SUN DESTINATION ${BIN_DIR})

# Static library
add_library(lib_damascus_sun STATIC
	Data_Generation.cpp
	Reflection_Spectrum.cpp
	Simulation_Trajectory.cpp
	Simulation_Utilities.cpp
	Solar_Model.cpp )

set_target_properties(lib_damascus_sun PROPERTIES PREFIX "")

target_compile_options(lib_damascus_sun PUBLIC -Wall -pedantic)

target_include_directories(lib_damascus_sun 
    PRIVATE
		${GENERATED_DIR})
		
target_link_libraries(lib_damascus_sun
    PUBLIC
        coverage_config 
		libobscura
		${MPI_CXX_LIBRARIES} )	

install(TARGETS lib_damascus_sun DESTINATION ${LIB_DIR})