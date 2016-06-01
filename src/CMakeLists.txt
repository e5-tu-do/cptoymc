# Executables

# Aliases for module libs
add_library( CPToyMC SHARED
  configuration/configuration.h 	configuration/configuration.cpp
  configuration/ToyConfig.h     	configuration/ToyConfig.cpp
  configuration/CompConfig.h    	configuration/CompConfig.cpp
  configuration/ObsConfig.h             configuration/ObsConfig.cpp
  generator/generator.h         	generator/generator.cpp
  generator/ToyGenerator.h      	generator/ToyGenerator.cpp
  generator/Observables.h       	generator/Observables.cpp
  generator/CompGenerator.h     	generator/CompGenerator.cpp
  generator/CompGeneratorFactory.h 	generator/CompGeneratorFactory.cpp
  generator/BSig_CPV_P2VP_Generator.h   generator/BSig_CPV_P2VP_Generator.cpp
  generator/LLBkg_Generator.h           generator/LLBkg_Generator.cpp
  generator/Bkg_2Final_Generator.h      generator/Bkg_2Final_Generator.cpp
 )

install(TARGETS CPToyMC DESTINATION lib)
install(FILES
configuration/configuration.h
configuration/ToyConfig.h
configuration/CompConfig.h
configuration/ObsConfig.h
DESTINATION include/configuration)
install(FILES
generator/generator.h
generator/ToyGenerator.h
generator/Observables.h
generator/CompGenerator.h
generator/CompGeneratorFactory.h
generator/BSig_CPV_P2VP_Generator.h
generator/LLBkg_Generator.h
generator/Bkg_2Final_Generator.h
DESTINATION include/generator)
