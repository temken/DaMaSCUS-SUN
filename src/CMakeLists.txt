# Executable
add_executable(DaMaSCUS-SUN
	Data_Generation.cpp
	main.cpp
	Reflection_Spectrum.cpp
	Simulation_Trajectory.cpp
	Simulation_Utilities.cpp
	Solar_Model.cpp )

target_compile_options(DaMaSCUS-SUN PUBLIC -Wall -pedantic)

target_link_libraries(DaMaSCUS-SUN
	PUBLIC
		coverage_config
		libobscura)

target_include_directories(DaMaSCUS-SUN
	PRIVATE
		${GENERATED_DIR} )

install(TARGETS DaMaSCUS-SUN DESTINATION ${BIN_DIR})