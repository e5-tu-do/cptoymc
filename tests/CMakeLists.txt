# Executables
add_executable(TaggingComboTest TaggingComboTest.cpp)
add_executable(OSSignalTest OSSignalTest.cpp)
add_executable(ComboSignalTest ComboSignalTest.cpp)
add_executable(SignalTest SignalTest.cpp)
add_executable(OSTest OSTest.cpp)
add_executable(ComboTest ComboTest.cpp)
add_executable(SpecialComboTest SpecialComboTest.cpp)
add_executable(SelfGenerationComboTest SelfGenerationComboTest.cpp)

# Aliases for module libs

# Linker information
target_link_libraries(TaggingComboTest CPToyMC ${DOOSOFT_LIBRARIES} ${ALL_LIBRARIES})
target_link_libraries(OSSignalTest CPToyMC ${DOOSOFT_LIBRARIES} ${ALL_LIBRARIES})
target_link_libraries(ComboSignalTest CPToyMC ${DOOSOFT_LIBRARIES} ${ALL_LIBRARIES})
target_link_libraries(SignalTest CPToyMC ${DOOSOFT_LIBRARIES} ${ALL_LIBRARIES})
target_link_libraries(OSTest CPToyMC ${DOOSOFT_LIBRARIES} ${ALL_LIBRARIES})
target_link_libraries(ComboTest CPToyMC ${DOOSOFT_LIBRARIES} ${ALL_LIBRARIES})
target_link_libraries(SpecialComboTest CPToyMC ${DOOSOFT_LIBRARIES} ${ALL_LIBRARIES})
target_link_libraries(SelfGenerationComboTest ${DOOSOFT_LIBRARIES} ${ALL_LIBRARIES})
