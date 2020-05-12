
# Executable
add_executable(DaMaSCUS-SUN
	main.cpp
	module.cpp)

target_compile_options(DaMaSCUS-SUN PUBLIC -Wall -pedantic)

target_link_libraries(DaMaSCUS-SUN
	PUBLIC 
		coverage_config
		libobscura)

install(TARGETS DaMaSCUS-SUN DESTINATION ${BIN_DIR})
