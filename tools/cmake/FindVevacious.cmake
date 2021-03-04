find_library(
  Vevacious_LIBRARY
  NAMES VevaciousPlusPlus
  PATH_SUFFIXES "lib")
find_path(
  Vevacious_INCLUDE_DIR
  NAMES VevaciousPlusPlus.hpp
  PATH_SUFFIXES "include")
find_path(
  Vevacious_LHPC_INCLUDE_DIR
  NAMES SimpleLhaParser.hpp
  PATH_SUFFIXES "include/LHPC")
message(${Vevacious_ROOT})

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(
  Vevacious
  FOUND_VAR Vevacious_FOUND
  REQUIRED_VARS Vevacious_LIBRARY Vevacious_INCLUDE_DIR Vevacious_LHPC_INCLUDE_DIR)

if(Vevacious_FOUND AND NOT TARGET Vevacious::libVevaciousPlusPlus)
  add_library(Vevacious::libVevaciousPlusPlus SHARED IMPORTED)
  set_target_properties(
    Vevacious::libVevaciousPlusPlus
    PROPERTIES INTERFACE_INCLUDE_DIRECTORIES "${Vevacious_INCLUDE_DIR};${Vevacious_LHPC_INCLUDE_DIR}"
               IMPORTED_LOCATION ${Vevacious_LIBRARY})
endif()

mark_as_advanced(Vevacious_LIBRARY Vevacious_INCLUDE_DIRECTORIES)
