# Executables
add_executable(DoToyGen DoToyGenMain.cpp)
# Aliases for module libs

# Linker information
target_link_libraries(DoToyGen CPToyMC ${ALL_LIBRARIES})
