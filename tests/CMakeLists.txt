# Executables
add_executable(TaggingComboTest TaggingComboTest.cpp)

# Aliases for module libs

# Linker information
target_link_libraries(TaggingComboTest CPToyMC ${DOOSOFT_LIBRARIES} ${ALL_LIBRARIES})
