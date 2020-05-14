
# Executable
add_executable(DaMaSCUS-SUN
	main.cpp
	module.cpp)

target_compile_options(DaMaSCUS-SUN PUBLIC -Wall -pedantic)

target_compile_definitions(DaMaSCUS-SUN
	PRIVATE 
		DAMASCUS_SUN_DIR="${PROJECT_SOURCE_DIR}/")

target_link_libraries(DaMaSCUS-SUN
	PUBLIC 
		coverage_config
		libobscura)

install(TARGETS DaMaSCUS-SUN DESTINATION ${BIN_DIR})
