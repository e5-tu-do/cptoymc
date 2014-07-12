# Executables

# Aliases for module libs
add_library( CPToyMC
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
 )

# Linker information

