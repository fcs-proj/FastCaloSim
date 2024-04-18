# Function to deactivate checks for a specified target
function(deactivate_target_checks target)
  set_target_properties(${target} PROPERTIES
    COMPILE_FLAGS "-w" # Deactivate compiler warnings
    CXX_CLANG_TIDY ""  # Deactivate clang-tidy
    CXX_CPPCHECK ""    # Deactivate cppcheck
  )
endfunction()

# Function that deactivates checks for a list of targets
# E.g. deactivate_checks(target1 target2 target3)
function(deactivate_checks)
  foreach(target IN LISTS ARGN)
    deactivate_target_checks(${target})
  endforeach()
endfunction()