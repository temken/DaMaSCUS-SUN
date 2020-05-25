# Executable
add_executable(DaMaSCUS-SUN
	main.cpp
	module.cpp
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