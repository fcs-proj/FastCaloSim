function(deactivate_checks)
  foreach(target IN LISTS ARGN)
    if (TARGET "${target}")
      # Resolve the actual target if it's an alias
      get_target_property(real_target "${target}" ALIASED_TARGET)
      if (NOT real_target)
        set(real_target "${target}")
      endif()

      # Check if the target is IMPORTED
      get_target_property(is_imported "${real_target}" IMPORTED)
      if (is_imported)
        continue()
      endif()

      # Get the target type
      get_target_property(target_type "${real_target}" TYPE)

      # Determine compile option scope
      if (target_type STREQUAL "INTERFACE_LIBRARY")
        set(compile_option_scope "INTERFACE")
      else()
        set(compile_option_scope "PRIVATE")
      endif()

      # Deactivate compiler warnings
      target_compile_options("${real_target}" ${compile_option_scope} "-w")

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
