# Function to create a single test executable out of multiple source files
function(create_test_executable suite_name test_sources)
  # Create the executable from the provided source file
  add_executable(${suite_name} ${test_sources})
  # Link the necessary libraries
  target_link_libraries(${suite_name} PRIVATE FastCaloSim::FastCaloSim GTest::gtest_main fmt::fmt)
  # Add include directory
  target_include_directories(${suite_name} PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include ${CMAKE_CURRENT_SOURCE_DIR})
  # Enable GTest discovery
  gtest_discover_tests(${suite_name} DISCOVERY_TIMEOUT 600)
  # Deactivate compiler warnings and static analyzers
  # Note: this is necessary until we have fixed all FastCaloSim related warnings
  deactivate_checks(${suite_name})
endfunction()
