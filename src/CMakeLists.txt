# Executables

# Aliases for module libs
add_library( CPToyMC
  configuration/configuration.h configuration/configuration.cpp 
  configuration/ToyConfig.h     configuration/ToyConfig.cpp
  configuration/CompConfig.h    configuration/CompConfig.cpp
  generator/generator.h         generator/generator.cpp
  generator/ToyGenerator.h      generator/ToyGenerator.cpp
  generator/Observables.h       generator/Observables.cpp
  generator/CompGenerator.h     generator/CompGenerator.cpp
)

# Linker information

