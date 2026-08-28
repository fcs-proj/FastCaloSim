cmake_minimum_required(VERSION 3.19)

macro(default name)
  if(NOT DEFINED "${name}")
    set("${name}" "${ARGN}")
  endif()
endmacro()

default(FORMAT_COMMAND clang-format)
# clang-format's output changes between releases, so the whole tree must be
# formatted with one version. The dev containers pin this exact version via
# docker/*/python_packages.txt; keep the two in sync when bumping.
default(EXPECTED_FORMAT_VERSION 20.1.8)
default(
    PATTERNS
    source/*.cxx source/*.h
    include/*.h
    test/*.cxx test/*.h
    example/*.cxx example/*.h
)
default(FIX NO)

set(flag --output-replacements-xml)
set(args OUTPUT_VARIABLE output)
if(FIX)
  set(flag -i)
  set(args "")
endif()

execute_process(
    COMMAND "${FORMAT_COMMAND}" --version
    RESULT_VARIABLE version_result
    OUTPUT_VARIABLE version_output
    OUTPUT_STRIP_TRAILING_WHITESPACE
)
if(NOT version_result EQUAL "0")
  message(FATAL_ERROR "'${FORMAT_COMMAND}': could not be run to query its version")
endif()
if(NOT version_output MATCHES "([0-9]+\\.[0-9]+\\.[0-9]+)")
  message(WARNING "Could not parse a version from '${version_output}'")
elseif(NOT CMAKE_MATCH_1 VERSION_EQUAL EXPECTED_FORMAT_VERSION)
  message(WARNING
      "clang-format ${CMAKE_MATCH_1} does not match the expected version "
      "${EXPECTED_FORMAT_VERSION}, so results will differ from CI. The dev "
      "containers ship the expected version; outside them install it with "
      "'pip install clang-format==${EXPECTED_FORMAT_VERSION}'."
  )
endif()

file(GLOB_RECURSE files ${PATTERNS})
set(badly_formatted "")
set(output "")
string(LENGTH "${CMAKE_SOURCE_DIR}/" path_prefix_length)

foreach(file IN LISTS files)
  execute_process(
      COMMAND "${FORMAT_COMMAND}" --style=file "${flag}" "${file}"
      WORKING_DIRECTORY "${CMAKE_SOURCE_DIR}"
      RESULT_VARIABLE result
      ${args}
  )
  if(NOT result EQUAL "0")
    message(FATAL_ERROR "'${file}': formatter returned with ${result}")
  endif()
  if(NOT FIX AND output MATCHES "\n<replacement offset")
    string(SUBSTRING "${file}" "${path_prefix_length}" -1 relative_file)
    list(APPEND badly_formatted "${relative_file}")
  endif()
  set(output "")
endforeach()

if(NOT badly_formatted STREQUAL "")
  list(JOIN badly_formatted "\n" bad_list)
  message("The following files are badly formatted:\n\n${bad_list}\n")
  message(FATAL_ERROR "Run again with FIX=YES to fix these files.")
endif()
