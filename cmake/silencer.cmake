# Function to deactivate checks for specified targets
function(deactivate_checks)
  foreach(target IN LISTS ARGN)
    if (TARGET "${target}")
      # Resolve the real target if it's an alias
      get_target_property(real_target "${target}" ALIASED_TARGET)
      if (NOT real_target)
        set(real_target "${target}")
      endif()

      # Deactivate compiler warnings
      target_compile_options("${real_target}" PRIVATE "-w")

      # Deactivate clang-tidy and cppcheck
      set_target_properties("${real_target}" PROPERTIES
        CXX_CLANG_TIDY ""
        CXX_CPPCHECK  ""
      )
    else()
      message(WARNING "Target '${target}' not found.")
    endif()
  endforeach()
endfunction()
