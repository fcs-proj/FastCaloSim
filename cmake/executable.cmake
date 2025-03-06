# Adds an executable target with the given name and sources.
# Can be used to add executables to the project.
# Can be executed with: cmake --build . --target run_<NAME>

function(add_exec NAME)
  add_executable("${NAME}" "executables/${NAME}.cxx")

  target_include_directories("${NAME}" PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/include
  )

  target_link_libraries("${NAME}" PRIVATE
    FastCaloSim::FastCaloSim
    FastCaloSimParam::FastCaloSimParam
    fmt::fmt-header-only
  )

  target_compile_features("${NAME}" PRIVATE cxx_std_17)

  add_custom_target("run_${NAME}"
    COMMAND "${NAME}"
    VERBATIM
  )

  add_dependencies("run_${NAME}" "${NAME}")

  deactivate_checks(${NAME})
endfunction()
