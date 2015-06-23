# Executables
add_executable(DoToyGen DoToyGenMain.cpp)
add_executable(P2VV_Generator P2VV_Generator.cpp)
add_executable(P2VV_Validation P2VV_Validation.cpp)

# Aliases for module libs

# Linker information
target_link_libraries(DoToyGen CPToyMC ${ALL_LIBRARIES})
target_link_libraries(P2VV_Generator CPToyMC ${ALL_LIBRARIES})
target_link_libraries(P2VV_Validation CPToyMC ${ALL_LIBRARIES})